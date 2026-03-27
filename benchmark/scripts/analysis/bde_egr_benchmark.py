#!/usr/bin/env python3
"""
bde_egr_benchmark.py
====================
Benchmark bond-dissociation-energy based effective graph resistance (BDE-EGR)
features on ToxCast invitrodb and Arp & Hale PTBM endpoints.

Two BDE-EGR representations are tested:
  1. PFG_BDE_EGR+mol  – PFASGroups component-level BDE-EGR (1 value per group ×
                         molecule, matched to PFG_EGR+mol which uses topological EGR).
  2. mol_kf           – Whole-molecule BDE Kirchhoff index (2 scalars: raw value
                         and value normalised by n*(n-1)/2).

Experiments
-----------
  Exp A  – PFAS chemicals only (~808), invitrodb 15 binary ToxCast endpoints.
            CV: RepeatedStratifiedKFold(3 splits × 5 repeats = 15 outer folds)
            Bayesian baselines: PFG_EGR+mol, Morgan
  Exp B  – All invitrodb chemicals (~9 014), same 15 binary endpoints.
            CV: StratifiedKFold(5 folds)
            Bayesian baselines: Morgan+PFG_EGR+mol, Morgan
  PTBM   – Arp & Hale 2022 all chemicals, 9 binary hazard-label endpoints.
            CV: StratifiedKFold(5 folds)
            Bayesian baseline: Morgan

Bayesian comparison
-------------------
  Correlated t-test (Benavoli et al. 2017) with ROPE = ±0.01.
  Metrics: ROC-AUC, Avg. Precision, MCC, Balanced Accuracy.
  Produces bar charts (P(BDE > baseline)) and posterior KDE figures.

Outputs
-------
  data/bde_egr_comparison_results.csv          per-fold metrics
  data/bde_egr_bayesian_<pair>.csv             per-endpoint Bayesian posteriors
  imgs/bde_egr_bayesian_<pair>.png/.pdf        posterior bar chart
  imgs/bde_egr_posteriors_<pair>.png/.pdf      posterior KDE densities

Usage (from benchmark/)
-----------------------
  conda activate chem
  python scripts/analysis/bde_egr_benchmark.py [--db invitrodb_v4_3] [--jobs N]
"""
from __future__ import annotations

import argparse
import hashlib
import math
import os
import sys
import time
import warnings
from pathlib import Path

warnings.filterwarnings("ignore")
os.environ["PYTHONWARNINGS"] = "ignore"

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.font_manager as fm

from rdkit import Chem
from rdkit.Chem import AllChem
from scipy.stats import gaussian_kde
from sklearn.ensemble import HistGradientBoostingClassifier, RandomForestClassifier
from sklearn.metrics import (
    average_precision_score,
    balanced_accuracy_score,
    matthews_corrcoef,
    roc_auc_score,
)
from sklearn.model_selection import (
    GridSearchCV,
    RepeatedStratifiedKFold,
    StratifiedKFold,
)

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
SCRIPT_DIR = Path(__file__).resolve().parent
DATA_DIR   = SCRIPT_DIR.parents[1] / "data"
TEST_DATA  = SCRIPT_DIR.parents[1] / "test_data"
IMGS_DIR   = SCRIPT_DIR.parents[1] / "imgs"
IMGS_DIR.mkdir(exist_ok=True)

STRUCTURES_XLSX = Path(
    r"c:\Users\luc\kDrive\Documents\WORK\PhD\vPM-vPB gap\data\arphale_only_structures.xlsx"
)
SIMPLIFIED_XLSX = TEST_DATA / "arpHale_simplified.xlsx"

# Add workspace root to sys.path so PFASGroups is importable
_WORKSPACE_ROOT = SCRIPT_DIR.resolve().parents[2]
if str(_WORKSPACE_ROOT) not in sys.path:
    sys.path.insert(0, str(_WORKSPACE_ROOT))

# Bayesian tests library
_LIB_DIR = str(SCRIPT_DIR.parent / "lib")
if _LIB_DIR not in sys.path:
    sys.path.insert(0, _LIB_DIR)

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------
RANDOM_STATE   = 42
MORGAN_NBITS   = 512
MORGAN_RADIUS  = 2
INNER_SPLITS   = 3
OUTER_SPLITS   = 5
OUTER_REPEATS  = 5   # Exp A
OUTER_SPLITS_A = 3   # folds for Exp A outer loop (small n)
N_REPEATS_A    = 5   # repeats used in correlated_ttest for Exp A
N_REPEATS_B    = 1   # no repeats for Exp B / PTBM

ROPE    = 0.01
METRICS = ["roc_auc", "avg_prec", "mcc", "bal_acc"]
METRIC_LABELS = {
    "roc_auc":  "ROC-AUC",
    "avg_prec": "Avg. Precision",
    "mcc":      "MCC",
    "bal_acc":  "Bal. Accuracy",
}
METRIC_COLORS = {
    "roc_auc":  "#EB935C",
    "avg_prec": "#759ED1",
    "mcc":      "#BE6A9D",
    "bal_acc":  "#8B61A8",
}
_SIG_COLORS = {"***": "#B8860B", "**": "#777777", "*": "#8B4513"}

LABEL_COLS = [
    "AR_antagonist", "AR_agonist",
    "ERa_antagonist", "ERa_agonist",
    "AhR_agonist", "Aromatase_antagonist",
    "TR_antagonist", "DT40_genotoxicity",
    "MMP_ratio",
    "CYP2D6_antagonist", "CYP2C19_antagonist", "CYP2C9_antagonist",
    "CYP3A4_antagonist", "p53_ratio", "Caspase3_HEPG2",
]
MIN_POS_PFAS = 10
MIN_POS_FULL = 30

MOL_METRICS = [
    "n_components", "total_size", "mean_size", "max_size",
    "mean_branching", "max_branching", "mean_eccentricity",
    "max_diameter", "mean_component_fraction", "max_component_fraction",
]

PTBM_LABEL_COLS = [
    "is_P", "is_M", "is_vPvM", "is_PMT",
    "is_T", "is_T_broad", "is_B",
    "is_P_simp", "is_M_simp",
]

# ---------------------------------------------------------------------------
# Font
# ---------------------------------------------------------------------------
_ubuntu_ok = any("Ubuntu" in f.name for f in fm.fontManager.ttflist)
_FONT = "Ubuntu" if _ubuntu_ok else "sans-serif"
plt.rcParams.update({"font.family": _FONT, "font.size": 9})


# ---------------------------------------------------------------------------
# Whole-molecule BDE Kirchhoff index
# ---------------------------------------------------------------------------

def _mol_kirchhoff_bde(mol) -> tuple[float, float]:
    """
    Compute the BDE-weighted Kirchhoff index of a molecule.

    Edge conductance = BDE(z_i, z_j, bond_order) / ref_bde  (same scheme as
    PFASGroups' effective_graph_resistance_BDE but applied to the full molecular
    graph including H-suppressed non-C atoms and heteroatoms).

    Returns (kf_raw, kf_norm) where kf_norm = kf_raw / (n*(n-1)/2).
    If the molecule is None or too small, returns (nan, nan).
    """
    from PFASGroups.ComponentsSolverModel import _get_bde_scheme  # noqa: PLC0415

    if mol is None:
        return float("nan"), float("nan")

    mol = Chem.RWMol(mol)
    n = mol.GetNumAtoms()
    if n < 2:
        return 0.0, 0.0

    bde = _get_bde_scheme()

    # Build weighted Laplacian
    L = np.zeros((n, n), dtype=np.float64)
    for bond in mol.GetBonds():
        i = bond.GetBeginAtomIdx()
        j = bond.GetEndAtomIdx()
        z_i = mol.GetAtomWithIdx(i).GetAtomicNum()
        z_j = mol.GetAtomWithIdx(j).GetAtomicNum()
        bo  = bond.GetBondTypeAsDouble()
        w   = bde.conductance(z_i, z_j, bo)
        L[i, j] -= w
        L[j, i] -= w
        L[i, i] += w
        L[j, j] += w

    # Pseudoinverse via eigendecomposition (numerically stable for small n)
    eigvals = np.linalg.eigvalsh(L)
    # Kirchhoff index = n * sum(1/λ_k for λ_k > threshold)
    tol = 1e-8 * max(eigvals)
    pos_eigvals = eigvals[eigvals > tol]
    if len(pos_eigvals) == 0:
        return 0.0, 0.0
    kf_raw = n * float(np.sum(1.0 / pos_eigvals))
    kf_norm = kf_raw / (n * (n - 1) / 2.0) if n > 1 else 0.0
    return kf_raw, kf_norm


def build_mol_kf_cache(smiles_list: list[str], cache_path: Path) -> np.ndarray:
    """
    Compute (or load from cache) the [kf_raw, kf_norm] feature matrix.

    Returns shape (n, 2) float64 array.
    """
    if cache_path.exists():
        arr = np.load(cache_path)
        if arr.shape[0] == len(smiles_list):
            print(f"  [cache hit] mol_kf: shape={arr.shape}")
            return arr

    print(f"  Computing whole-molecule BDE Kirchhoff index for {len(smiles_list)} SMILES …")
    t0 = time.time()
    rows = []
    for smi in smiles_list:
        mol = Chem.MolFromSmiles(smi)
        rows.append(_mol_kirchhoff_bde(mol))
    arr = np.array(rows, dtype=np.float64)
    np.save(cache_path, arr)
    print(f"  [cache saved] mol_kf {arr.shape}  ({time.time()-t0:.1f}s)")
    return arr


# ---------------------------------------------------------------------------
# PFASGroups BDE-EGR embedding
# ---------------------------------------------------------------------------

def build_pfg_bde_matrices(
    smiles_list: list[str],
    cache_dir: Path = DATA_DIR,
    cache_prefix: str = "pfg_BDE_EGR",
) -> dict[str, np.ndarray]:
    """
    Compute PFASGroups BDE-EGR component embedding with and without mol metrics.

    Configs produced:
      "BDE_EGR"       – 115 cols (effective_graph_resistance_BDE only)
      "BDE_EGR+mol"   – 125 cols (+ 10 molecule-wide graph metrics)

    Returns dict name → (n, ncols) float64 array.  Results are cached.
    """
    from PFASGroups import parse_smiles as _parse   # noqa: PLC0415
    from PFASGroups.getter import get_compiled_HalogenGroups as _get_groups  # noqa: PLC0415

    configs = {
        "BDE_EGR":     dict(component_metrics=["effective_graph_resistance_BDE"]),
        "BDE_EGR+mol": dict(component_metrics=["effective_graph_resistance_BDE"],
                            molecule_metrics=MOL_METRICS),
    }

    n = len(smiles_list)
    cache_meta  = cache_dir / f"{cache_prefix}_cache_meta.txt"

    if n == 0:
        raise ValueError("smiles_list is empty — nothing to compute.")

    # Check per-config caches (shape-only, no SMILES hash needed for re-use)
    if cache_dir:
        matrices: dict[str, np.ndarray] = {}
        all_ok = True
        for name in configs:
            p = cache_dir / f"{cache_prefix}_{name}_cache.npy"
            if p.exists():
                arr = np.load(p)
                if arr.shape[0] == n:
                    matrices[name] = arr
                    print(f"  [cache hit] PFG {name}: shape={arr.shape}")
                    continue
            all_ok = False
            break
        if all_ok and len(matrices) == len(configs):
            return matrices

    print(f"  [cache miss] Computing PFASGroups BDE-EGR for {n} SMILES …")

    # Deduplicate by InChIKey
    inchikeys: list[str | None] = []
    for smi in smiles_list:
        mol = Chem.MolFromSmiles(smi)
        inchikeys.append(Chem.inchi.MolToInchiKey(mol) if mol else None)

    unique_order: dict[str, int] = {}
    unique_smiles: list[str] = []
    for smi, ik in zip(smiles_list, inchikeys):
        if ik and ik not in unique_order:
            unique_order[ik] = len(unique_smiles)
            unique_smiles.append(smi)

    n_unique = len(unique_smiles)
    n_invalid = sum(1 for ik in inchikeys if ik is None)
    print(f"  Unique: {n_unique}  |  invalid: {n_invalid}")

    t0 = time.time()
    result = _parse(unique_smiles, halogens="F", progress=True)
    pfas_groups = _get_groups(halogens="F")

    result_ik_map: dict[str, int] = {}
    for idx, emb in enumerate(result):
        ik = emb.get("inchikey", "")
        if ik:
            result_ik_map[ik] = idx

    input_to_result = np.full(n, -1, dtype=int)
    for i, ik in enumerate(inchikeys):
        if ik and ik in result_ik_map:
            input_to_result[i] = result_ik_map[ik]
    valid_mask    = input_to_result >= 0
    valid_indices = input_to_result[valid_mask]

    matrices = {}
    for name, kwargs in configs.items():
        print(f"  Computing PFG {name} …")
        arr_result = np.asarray(result.to_array(
            **kwargs, progress=True, pfas_groups=pfas_groups
        ))
        arr = np.zeros((n, arr_result.shape[1]), dtype=np.float64)
        arr[valid_mask] = arr_result[valid_indices]
        matrices[name] = arr
        nonzero = int((np.abs(arr).sum(axis=1) > 0).sum())
        print(f"    shape={arr.shape}  nonzero rows={nonzero}  ({time.time()-t0:.1f}s)")

    # Save caches
    cache_dir.mkdir(parents=True, exist_ok=True)
    for name, arr in matrices.items():
        np.save(cache_dir / f"{cache_prefix}_{name}_cache.npy", arr)
    print(f"  [cache saved] {len(matrices)} BDE-EGR matrices")
    return matrices


# ---------------------------------------------------------------------------
# Morgan fingerprints
# ---------------------------------------------------------------------------

def build_morgan(smiles_list: list[str]) -> np.ndarray:
    X = np.zeros((len(smiles_list), MORGAN_NBITS), dtype=np.uint8)
    failed = 0
    for i, smi in enumerate(smiles_list):
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            failed += 1
            continue
        fp = AllChem.GetMorganFingerprintAsBitVect(mol, MORGAN_RADIUS, MORGAN_NBITS)
        X[i] = np.array(fp)
    if failed:
        print(f"  Warning: {failed} SMILES failed Morgan generation")
    print(f"  Morgan{MORGAN_NBITS}: shape={X.shape}  nonzero rows={int((X.sum(1) > 0).sum())}")
    return X


# ---------------------------------------------------------------------------
# Nested CV
# ---------------------------------------------------------------------------

def _rf_grid(seed: int) -> GridSearchCV:
    return GridSearchCV(
        RandomForestClassifier(class_weight="balanced", random_state=seed, n_jobs=1),
        {"n_estimators": [200, 400], "max_features": ["sqrt", 0.2]},
        cv=StratifiedKFold(INNER_SPLITS, shuffle=True, random_state=seed),
        scoring="roc_auc", n_jobs=-1, refit=True,
    )


def _gb_grid(seed: int) -> GridSearchCV:
    return GridSearchCV(
        HistGradientBoostingClassifier(class_weight="balanced", random_state=seed),
        {"max_depth": [3, 5], "learning_rate": [0.05, 0.1]},
        cv=StratifiedKFold(INNER_SPLITS, shuffle=True, random_state=seed),
        scoring="roc_auc", n_jobs=-1, refit=True,
    )


def nested_cv_metrics(
    X: np.ndarray,
    y: np.ndarray,
    outer_cv,
    label: str,
    fset: str,
    experiment: str,
) -> list[dict]:
    """Run nested CV (GradientBoosting only) and return per-fold dicts with all 4 metrics."""
    rows: list[dict] = []
    splits = list(outer_cv.split(X, y))
    has_repeats = hasattr(outer_cv, "n_repeats")
    n_splits_per = outer_cv.cvargs["n_splits"] if has_repeats else outer_cv.n_splits

    for split_i, (tr, te) in enumerate(splits):
        repeat_idx = split_i // n_splits_per if has_repeats else 0
        fold_idx   = split_i % n_splits_per

        X_tr, X_te = X[tr], X[te]
        y_tr, y_te = y[tr], y[te]

        if len(np.unique(y_tr)) < 2 or len(np.unique(y_te)) < 2:
            continue

        grid = _gb_grid(RANDOM_STATE)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            grid.fit(X_tr, y_tr)
            prob = grid.predict_proba(X_te)[:, 1]
            pred = (prob >= 0.5).astype(int)

        try:
            auc = roc_auc_score(y_te, prob)
        except ValueError:
            auc = float("nan")
        try:
            ap = average_precision_score(y_te, prob)
        except ValueError:
            ap = float("nan")

        rows.append({
            "experiment":  experiment,
            "feature_set": fset,
            "endpoint":    label,
            "model":       "GradientBoosting",
            "repeat":      repeat_idx,
            "fold":        fold_idx,
            "n_total":     len(y),
            "n_pos":       int(y.sum()),
            "roc_auc":     round(float(auc), 4),
            "avg_prec":    round(float(ap),  4),
            "mcc":         round(float(matthews_corrcoef(y_te, pred)), 4),
            "bal_acc":     round(float(balanced_accuracy_score(y_te, pred)), 4),
            "best_params": str(grid.best_params_),
        })

    return rows


# ---------------------------------------------------------------------------
# Bayesian comparison helpers
# ---------------------------------------------------------------------------

def run_bayesian_comparison(
    results_df: pd.DataFrame,
    experiment: str,
    baseline_fset: str,
    compare_fsets: list[str],
    n_repeats: int,
    pair_label: str,
    model: str = "GradientBoosting",
) -> pd.DataFrame:
    """
    Compute Bayesian correlated t-test for each (endpoint, metric, feature_set)
    pair vs baseline. Returns a DataFrame of posterior probabilities.
    """
    from bayesiantests import correlated_ttest  # noqa: PLC0415

    df = results_df[
        (results_df["experiment"] == experiment) &
        (results_df["model"] == model)
    ].copy()

    endpoints = sorted(df["endpoint"].unique())
    rows = []

    for ep in endpoints:
        ep_df = df[df["endpoint"] == ep]
        base_df = ep_df[ep_df["feature_set"] == baseline_fset][
            ["repeat", "fold"] + METRICS
        ]
        if base_df.empty:
            continue

        for fs in compare_fsets:
            cmp_df = ep_df[ep_df["feature_set"] == fs][
                ["repeat", "fold"] + METRICS
            ]
            if cmp_df.empty:
                continue

            for metric in METRICS:
                merged = base_df[["repeat", "fold", metric]].merge(
                    cmp_df[["repeat", "fold", metric]],
                    on=["repeat", "fold"],
                    suffixes=("_base", "_cmp"),
                )
                if merged.empty:
                    continue

                diffs = (merged[f"{metric}_cmp"].values
                         - merged[f"{metric}_base"].values)
                p_base, p_eq, p_cmp = correlated_ttest(diffs, ROPE, runs=n_repeats)

                rows.append({
                    "pair":         pair_label,
                    "experiment":   experiment,
                    "endpoint":     ep,
                    "metric":       metric,
                    "feature_set":  fs,
                    "baseline":     baseline_fset,
                    "mean_diff":    round(float(np.mean(diffs)), 4),
                    "p_better":     round(float(p_cmp),  4),
                    "p_rope":       round(float(p_eq),   4),
                    "p_worse":      round(float(p_base), 4),
                    "n_pairs":      len(diffs),
                })

    return pd.DataFrame(rows)


def _sig_stars(p: float) -> tuple[str, str]:
    """Return (stars, color) for a posterior probability."""
    if p >= 0.999:
        return "***", _SIG_COLORS["***"]
    if p >= 0.99:
        return "**",  _SIG_COLORS["**"]
    if p >= 0.95:
        return "*",   _SIG_COLORS["*"]
    return "", "#444"


def plot_bayesian_bar(
    bayes_df: pd.DataFrame,
    pair_label: str,
    baseline_label: str,
    title_suffix: str = "",
) -> None:
    """
    Horizontal bar chart: P(feature_set > baseline) per endpoint × metric.
    One column per metric, one row per endpoint.
    """
    if bayes_df.empty:
        return

    endpoints = sorted(bayes_df["endpoint"].unique())
    fsets = sorted(bayes_df["feature_set"].unique())
    n_ep = len(endpoints)
    n_fs = len(fsets)

    fig, axes = plt.subplots(
        1, len(METRICS),
        figsize=(len(METRICS) * 3.8, max(4.0, n_ep * 0.5 * n_fs + 1.5)),
        constrained_layout=True,
    )
    if len(METRICS) == 1:
        axes = [axes]

    fig.suptitle(
        f"Bayesian Correlated t-test: P(BDE-EGR > baseline = {baseline_label})\n"
        f"{title_suffix}  (GradientBoosting, ROPE = ±{ROPE})",
        fontsize=13, fontweight="bold",
    )

    cmap = plt.get_cmap("tab10")
    fs_colors = {fs: cmap(i / max(n_fs, 1)) for i, fs in enumerate(fsets)}

    for col_idx, (ax, metric) in enumerate(zip(axes, METRICS)):
        sub = bayes_df[bayes_df["metric"] == metric]

        y_labels = []
        y_pos = []
        y_idx = 0
        for ep in endpoints:
            for fs in fsets:
                row = sub[(sub["endpoint"] == ep) & (sub["feature_set"] == fs)]
                if row.empty:
                    continue
                p     = float(row.iloc[0]["p_better"])
                d     = float(row.iloc[0]["mean_diff"])
                stars, sig_color = _sig_stars(p)

                if y_idx % 2 == 0:
                    ax.axhspan(y_idx - 0.5, y_idx + 0.5,
                               color="#f0f0f0", zorder=0, linewidth=0)

                ax.barh(y_idx, p, color=fs_colors[fs],
                        edgecolor="white", linewidth=0.4, zorder=2)
                sign = "+" if d >= 0 else ""
                ax.text(min(p + 0.02, 0.98), y_idx,
                        f"{sign}{d:.3f}{stars}",
                        va="center", ha="left", fontsize=10,
                        color=sig_color, fontweight="bold" if stars else "normal",
                        zorder=3)
                y_labels.append(f"{ep}\n({fs})" if n_fs > 1 else ep)
                y_pos.append(y_idx)
                y_idx += 1

        ax.axvline(0.5,  color="#aaa", linewidth=0.8, linestyle="--", zorder=1)
        ax.axvline(0.95, color="#888", linewidth=0.6, linestyle=":",  zorder=1)
        ax.set_xlim(0, 1.15)
        ax.set_ylim(-0.5, max(y_pos, default=0) + 0.5)
        ax.set_yticks(y_pos)
        ax.set_xlabel(f"P(BDE > {baseline_label})", fontsize=11)
        ax.set_title(METRIC_LABELS[metric], fontsize=12, fontweight="bold")
        ax.grid(axis="x", linewidth=0.3, alpha=0.5)
        ax.spines[["top", "right"]].set_visible(False)

        if col_idx == 0:
            ax.set_yticklabels(y_labels, fontsize=9)
        else:
            ax.set_yticklabels([])
            ax.tick_params(axis="y", length=0)
            ax.spines["left"].set_visible(False)

    safe = pair_label.replace(" ", "_").replace("+", "plus")
    for ext in ("png", "pdf"):
        out = IMGS_DIR / f"bde_egr_bayesian_{safe}.{ext}"
        fig.savefig(out, dpi=150, bbox_inches="tight")
        print(f"[saved] {out.relative_to(IMGS_DIR.parent)}")
    plt.close(fig)


def plot_bayesian_posteriors(
    results_df: pd.DataFrame,
    bayes_df: pd.DataFrame,
    experiment: str,
    baseline_fset: str,
    compare_fsets: list[str],
    n_repeats: int,
    pair_label: str,
    model: str = "GradientBoosting",
) -> None:
    """
    Posterior KDE figure: rows = endpoints, columns = metrics.
    Shows distribution of mean differences (BDE-EGR − baseline).
    """
    from bayesiantests import correlated_ttest_MC  # noqa: PLC0415

    df = results_df[
        (results_df["experiment"] == experiment) &
        (results_df["model"] == model)
    ].copy()

    # For multi-feature-set comparisons, pick the best fs per (ep, metric)
    endpoints = sorted(bayes_df["endpoint"].unique())
    n_eps = len(endpoints)
    n_met = len(METRICS)

    if n_eps == 0:
        return

    fig, axes = plt.subplots(
        n_eps, n_met,
        figsize=(n_met * 3.2, n_eps * 1.35),
        constrained_layout=True,
        squeeze=False,
    )
    fig.suptitle(
        rf"Posterior distribution of mean difference (BDE-EGR $-$ {baseline_fset})"
        + "\nOrange lines: ROPE boundaries (±0.01). "
        "Red=baseline better, grey=equivalent, colour=BDE-EGR better.",
        fontsize=12, fontweight="bold",
    )

    for col, metric in enumerate(METRICS):
        axes[0, col].set_title(METRIC_LABELS[metric], fontsize=12, fontweight="bold")

    for row, ep in enumerate(endpoints):
        axes[row, 0].set_ylabel(
            ep.replace("_", " "), fontsize=10, rotation=0, ha="right", labelpad=4,
        )
        ep_df = df[df["endpoint"] == ep]

        for col, metric in enumerate(METRICS):
            ax = axes[row, col]

            # Select best feature set (highest mean_diff) for this (ep, metric)
            sub = bayes_df[
                (bayes_df["endpoint"] == ep) & (bayes_df["metric"] == metric)
            ]
            if sub.empty:
                ax.axis("off")
                continue

            best_row = sub.loc[sub["mean_diff"].idxmax()]
            best_fs  = best_row["feature_set"]

            base_df = ep_df[ep_df["feature_set"] == baseline_fset][
                ["repeat", "fold", metric]
            ]
            cmp_df  = ep_df[ep_df["feature_set"] == best_fs][
                ["repeat", "fold", metric]
            ]
            merged = base_df.merge(cmp_df, on=["repeat", "fold"],
                                   suffixes=("_base", "_cmp"))
            if merged.empty:
                ax.axis("off")
                continue

            diffs   = merged[f"{metric}_cmp"].values - merged[f"{metric}_base"].values
            samples = correlated_ttest_MC(diffs, ROPE, runs=n_repeats, nsamples=40000)

            p1, p99 = np.percentile(samples, [0.5, 99.5])
            xs = np.linspace(min(p1, -5 * ROPE), max(p99, 5 * ROPE), 400)
            ys = gaussian_kde(samples)(xs)

            ax.fill_between(xs, ys, where=(xs < -ROPE),
                            color="#E07070", alpha=1, linewidth=0)
            ax.fill_between(xs, ys, where=((xs >= -ROPE) & (xs <= ROPE)),
                            color="#C0C0C0", alpha=1, linewidth=0)
            ax.fill_between(xs, ys, where=(xs > ROPE),
                            color=METRIC_COLORS[metric], alpha=1, linewidth=0)

            ax.plot(xs, ys, color="#222", linewidth=0.9)
            ax.axvline(-ROPE, color="orange", linewidth=0.8)
            ax.axvline( ROPE, color="orange", linewidth=0.8)
            ax.axvline(0,     color="#888",   linewidth=0.5, linestyle="--")

            if row % 2 == 0:
                ax.axhspan(0, max(ys) + 0.5, color="#f0f0f0", zorder=0, linewidth=0)

            mean_d = float(best_row["mean_diff"])
            p_bde  = float(best_row["p_better"])
            stars, sig_color = _sig_stars(p_bde)
            sign  = "+" if mean_d >= 0 else ""
            label = f"{sign}{mean_d:.3f}{stars}\n({p_bde:.0%})"
            ax.text(0.98, 0.97, label, ha="right", va="top", fontsize=10,
                    fontweight="bold" if stars else "normal", color=sig_color,
                    transform=ax.transAxes)
            ax.axvline(mean_d, color=sig_color, linewidth=0.8, alpha=1, zorder=4)

            ax.set_yticks([])
            ax.tick_params(axis="x", labelsize=9)
            ax.spines[["top", "right", "left"]].set_visible(False)

    safe = pair_label.replace(" ", "_").replace("+", "plus")
    for ext in ("png", "pdf"):
        out = IMGS_DIR / f"bde_egr_posteriors_{safe}.{ext}"
        fig.savefig(out, dpi=150, bbox_inches="tight")
        print(f"[saved] {out.relative_to(IMGS_DIR.parent)}")
    plt.close(fig)


# ---------------------------------------------------------------------------
# Data loading
# ---------------------------------------------------------------------------

def load_toxcast(db: str) -> tuple[list[str], np.ndarray, list[str]]:
    """
    Load invitrodb dataset.  Returns (smiles_list, Y, label_cols).
    Y shape: (n_chemicals, n_endpoints), dtype int8.
    """
    parquet = DATA_DIR / f"{db}_dataset.parquet"
    if not parquet.exists():
        raise FileNotFoundError(f"Parquet not found: {parquet}")

    df = pd.read_parquet(parquet)
    # Keep only rows with valid SMILES
    df = df[df["smiles"].notna() & (df["smiles"] != "")].reset_index(drop=True)
    smiles = df["smiles"].tolist()

    # Subset to available label columns
    labels = [c for c in LABEL_COLS if c in df.columns]
    Y = df[labels].values.astype(np.int8)
    print(f"  Loaded {len(smiles)} chemicals, {len(labels)} endpoints from {parquet.name}")
    return smiles, Y, labels


def load_ptbm() -> tuple[list[str], np.ndarray, list[str]]:
    """
    Load Arp & Hale PTBM dataset.  Returns (smiles_list, Y, label_cols).
    Matches the exact loading logic in mqg_ptbm_predictivity.py.
    """
    _VPMV_CONCERNED = {"PM", "PMT", "vPvM", "vPvM & PMT"}
    _VPMV_CLEAR     = {"Not P", "Not M", "Not P, Not M"}
    BCF_B_LOG = 3.30
    BCF_UNQ   = 2.00

    if not STRUCTURES_XLSX.exists():
        raise FileNotFoundError(f"PTBM structures file not found: {STRUCTURES_XLSX}")

    df = pd.read_excel(STRUCTURES_XLSX)
    print(f"  {len(df)} chemicals loaded from {STRUCTURES_XLSX.name}")

    # SMILES column is "Smiles" in the Arp & Hale structures file
    smiles = df["Smiles"].fillna("").tolist()

    # ------------------------------------------------------------------
    # P/M/vPvM/PMT labels from vPvM_EC column
    # ------------------------------------------------------------------
    col = df["vPvM_EC"].fillna("")

    def _p(s: str) -> float:
        if s in _VPMV_CONCERNED:
            return 1.0
        if s in {"Not P", "Not P, Not M"}:
            return 0.0
        return np.nan

    def _m(s: str) -> float:
        if s in _VPMV_CONCERNED:
            return 1.0
        if s in {"Not M", "Not P, Not M"}:
            return 0.0
        return np.nan

    def _vpmv(s: str) -> float:
        if s in {"vPvM", "vPvM & PMT"}:
            return 1.0
        if s in _VPMV_CLEAR:
            return 0.0
        return np.nan

    def _pmt(s: str) -> float:
        if s in {"PMT", "vPvM & PMT"}:
            return 1.0
        if s in _VPMV_CLEAR | {"PM", "vPvM"}:
            return 0.0
        return np.nan

    df["is_P"]    = col.map(_p)
    df["is_M"]    = col.map(_m)
    df["is_vPvM"] = col.map(_vpmv)
    df["is_PMT"]  = col.map(_pmt)

    # ------------------------------------------------------------------
    # T/B/P_simp/M_simp labels via join with arpHale_simplified.xlsx
    # ------------------------------------------------------------------
    if SIMPLIFIED_XLSX.exists():
        ah = pd.read_excel(SIMPLIFIED_XLSX)
        ah["_ec"] = ah["AH_CASRN"].fillna("").str.strip()
        ah_idx    = ah.drop_duplicates("_ec").set_index("_ec")
        ah_lookup = ah_idx[["AH_P", "AH_T", "AH_M_EC", "AH_M_UBA", "AH_B"]].to_dict("index")

        def _first_match(textid) -> "str | None":
            if pd.isna(textid):
                return None
            for part in str(textid).split(";"):
                ec = part.strip()
                if ec in ah_lookup:
                    return ec
            return None

        df["_match_ec"] = df["textid"].map(_first_match)
        n_matched = df["_match_ec"].notna().sum()
        print(f"  EC-number join: {n_matched}/{len(df)} rows matched "
              f"({100*n_matched/len(df):.1f}%)")

        for src_col in ("AH_P", "AH_T", "AH_M_EC", "AH_M_UBA", "AH_B"):
            col_map = {ec: ah_lookup[ec][src_col] for ec in ah_lookup}
            df[src_col] = df["_match_ec"].map(col_map)

        df["is_T"]       = df["AH_T"].map({"T": 1.0, "Not T": 0.0})
        df["is_T_broad"] = df["AH_T"].map({"T": 1.0, "Pot. T": 1.0, "Not T": 0.0})

        def _b(v) -> float:
            if pd.isna(v):
                return np.nan
            return 1.0 if v > BCF_B_LOG else (0.0 if v < BCF_UNQ else np.nan)
        df["is_B"] = df["AH_B"].map(_b)

        p_map = {"P": 1.0, "vP": 1.0, "Potential P/vP": 1.0,
                 "Potential P/vP++": 1.0, "Not P": 0.0}
        df["is_P_simp"] = df["AH_P"].map(p_map)
        m_map = {"M": 1.0, "vM": 1.0, "Pot. M/vM": 1.0, "Not M": 0.0}
        df["is_M_simp"] = df["AH_M_EC"].map(m_map)

        df.drop(columns=["_match_ec"], inplace=True)
    else:
        print(f"  [skip] simplified file not found: {SIMPLIFIED_XLSX}")
        for lbl in ["is_T", "is_T_broad", "is_B", "is_P_simp", "is_M_simp"]:
            df[lbl] = np.nan

    labels = [c for c in PTBM_LABEL_COLS if c in df.columns]
    Y = df[labels].values
    print(f"  Loaded {len(smiles)} chemicals, {len(labels)} PTBM endpoints")
    for lbl in labels:
        n1 = int(df[lbl].eq(1).sum())
        n0 = int(df[lbl].eq(0).sum())
        print(f"    {lbl:<14} pos={n1:4d}  neg={n0:4d}")
    return smiles, Y, labels


# ---------------------------------------------------------------------------
# Main pipeline
# ---------------------------------------------------------------------------

def run_experiment(
    exp_name: str,
    smiles: list[str],
    Y: np.ndarray,
    label_cols: list[str],
    feature_sets: dict[str, np.ndarray],
    outer_cv,
    min_pos: int,
    all_rows: list[dict],
    n_done: list[int],
    n_total_tasks: int,
) -> None:
    """
    Run nested CV for all feature sets and endpoints in one experiment.
    Appends result dicts to all_rows in-place.
    """
    for ep_idx, ep in enumerate(label_cols):
        y_raw = Y[:, ep_idx]
        valid = ~np.isnan(y_raw.astype(float))
        y = y_raw[valid].astype(int)
        n_pos = int(y.sum())
        n_neg = int((y == 0).sum())

        if n_pos < min_pos or n_neg < min_pos:
            print(f"    [{exp_name}] Skip {ep}: {n_pos} pos / {n_neg} neg (< {min_pos})")
            continue

        for fset_name, X_full in feature_sets.items():
            X = X_full[valid]
            t0 = time.time()
            rows = nested_cv_metrics(X, y, outer_cv, ep, fset_name, exp_name)
            all_rows.extend(rows)
            n_done[0] += 1

            n_folds = len(rows)  # 1 model only now
            mean_auc = (
                float(np.nanmean([r["roc_auc"] for r in rows
                                  if r["model"] == "GradientBoosting"]))
                if rows else float("nan")
            )
            print(
                f"    [{exp_name}] {ep} / {fset_name}: "
                f"GB AUC={mean_auc:.3f}  "
                f"folds={n_folds}  ({time.time()-t0:.1f}s)  "
                f"[{n_done[0]}/{n_total_tasks}]"
            )


def main(db: str = "invitrodb_v4_3") -> None:
    t_start = time.time()

    # ------------------------------------------------------------------
    # Load ToxCast data
    # ------------------------------------------------------------------
    print("\n=== Loading ToxCast invitrodb ===")
    smiles_tc, Y_tc, labels_tc = load_toxcast(db)
    n_tc = len(smiles_tc)

    # PFAS mask: chemicals that have ≥1 fluorine in SMILES (fast proxy)
    pfas_mask = np.array(["F" in s for s in smiles_tc], dtype=bool)
    print(f"  PFAS subset: {pfas_mask.sum()} / {n_tc} chemicals")

    smiles_pfas = [s for s, m in zip(smiles_tc, pfas_mask) if m]
    Y_pfas      = Y_tc[pfas_mask]

    # ------------------------------------------------------------------
    # Load PTBM data
    # ------------------------------------------------------------------
    print("\n=== Loading Arp & Hale PTBM dataset ===")
    smiles_ptbm, Y_ptbm, labels_ptbm = load_ptbm()

    # ------------------------------------------------------------------
    # Build / load Morgan fingerprints
    # ------------------------------------------------------------------
    print("\n=== Morgan fingerprints ===")
    morgan_tc   = build_morgan(smiles_tc)
    morgan_pfas = morgan_tc[pfas_mask]
    morgan_ptbm = build_morgan(smiles_ptbm)

    # ------------------------------------------------------------------
    # Existing EGR+mol cache (topological, PFAS-mode F)
    # ------------------------------------------------------------------
    print("\n=== Loading existing EGR+mol cache ===")
    egr_cache_path = DATA_DIR / "pfg_EGR+mol_cache.npy"
    if egr_cache_path.exists():
        egr_mol_tc = np.load(egr_cache_path)
        if egr_mol_tc.shape[0] == n_tc:
            egr_mol_pfas = egr_mol_tc[pfas_mask]
            print(f"  [cache hit] EGR+mol: shape={egr_mol_tc.shape}")
        else:
            print(f"  [cache shape mismatch] EGR+mol cached {egr_mol_tc.shape[0]} "
                  f"vs {n_tc} expected — will not use")
            egr_mol_tc = None
            egr_mol_pfas = None
    else:
        print("  [not found] EGR+mol cache — skipping topological EGR baselines")
        egr_mol_tc   = None
        egr_mol_pfas = None

    # Also load Morgan+PFG_EGR+mol for Exp B baseline
    morgan_egr_cache_path = DATA_DIR / "pfg_EGR+mol_cache.npy"  # same file, used as-is

    # ------------------------------------------------------------------
    # BDE-EGR component embedding
    # ------------------------------------------------------------------
    print("\n=== PFASGroups BDE-EGR (component-level) ===")
    bde_mats_tc   = build_pfg_bde_matrices(smiles_tc,   DATA_DIR, "pfg_BDE_EGR")
    bde_mats_pfas = {k: v[pfas_mask] for k, v in bde_mats_tc.items()}

    print(f"\n=== PFASGroups BDE-EGR (component-level) for PTBM ===")
    # Filter out empty SMILES before computing; map back after
    ptbm_valid_mask = np.array([bool(s and s.strip()) for s in smiles_ptbm])
    smiles_ptbm_valid = [s for s, m in zip(smiles_ptbm, ptbm_valid_mask) if m]
    if len(smiles_ptbm_valid) == 0:
        raise ValueError("No valid PTBM SMILES found — check STRUCTURES_XLSX path and Smiles column.")
    print(f"  PTBM: {ptbm_valid_mask.sum()} / {len(smiles_ptbm)} valid SMILES")
    bde_mats_ptbm_valid = build_pfg_bde_matrices(smiles_ptbm_valid, DATA_DIR, "pfg_BDE_EGR_ptbm")
    # Re-expand to full PTBM size (zeros for empty SMILES rows)
    bde_mats_ptbm: dict[str, np.ndarray] = {}
    for k, arr_v in bde_mats_ptbm_valid.items():
        arr_full = np.zeros((len(smiles_ptbm), arr_v.shape[1]), dtype=np.float64)
        arr_full[ptbm_valid_mask] = arr_v
        bde_mats_ptbm[k] = arr_full

    # ------------------------------------------------------------------
    # Whole-molecule BDE Kirchhoff index
    # ------------------------------------------------------------------
    print("\n=== Whole-molecule BDE Kirchhoff index ===")
    mol_kf_tc   = build_mol_kf_cache(smiles_tc,   DATA_DIR / "mol_kf_tc_cache.npy")
    mol_kf_ptbm = build_mol_kf_cache(smiles_ptbm, DATA_DIR / "mol_kf_ptbm_cache.npy")
    mol_kf_pfas = mol_kf_tc[pfas_mask]

    # ------------------------------------------------------------------
    # Assemble feature set dicts per experiment
    # ------------------------------------------------------------------
    print("\n=== Assembling feature sets ===")

    def _h(X: np.ndarray) -> np.ndarray:
        """Replace NaN with 0 and return float32 (safe for all models)."""
        X = np.where(np.isnan(X), 0.0, X)
        return X.astype(np.float32)

    # ── Exp A: PFAS subset ──────────────────────────────────────────────
    fsets_A: dict[str, np.ndarray] = {
        "Morgan": _h(morgan_pfas),
    }
    if egr_mol_pfas is not None:
        fsets_A["PFG_EGR+mol"]         = _h(egr_mol_pfas)
        fsets_A["PFG_EGR+mol+mol_kf"]  = _h(np.hstack([egr_mol_pfas, mol_kf_pfas]))
    fsets_A["PFG_BDE_EGR+mol"]          = _h(bde_mats_pfas["BDE_EGR+mol"])
    fsets_A["PFG_BDE_EGR+mol+mol_kf"]   = _h(np.hstack([bde_mats_pfas["BDE_EGR+mol"], mol_kf_pfas]))
    fsets_A["Morgan+mol_kf"]             = _h(np.hstack([morgan_pfas, mol_kf_pfas]))

    # ── Exp B: all ToxCast ──────────────────────────────────────────────
    fsets_B: dict[str, np.ndarray] = {
        "Morgan": _h(morgan_tc),
    }
    if egr_mol_tc is not None:
        fsets_B["Morgan+PFG_EGR+mol"]     = _h(np.hstack([morgan_tc, egr_mol_tc]))
    fsets_B["Morgan+PFG_BDE_EGR+mol"]     = _h(np.hstack([morgan_tc, bde_mats_tc["BDE_EGR+mol"]]))
    fsets_B["Morgan+mol_kf"]              = _h(np.hstack([morgan_tc, mol_kf_tc]))

    # ── PTBM: all chemicals ─────────────────────────────────────────────
    fsets_PTBM: dict[str, np.ndarray] = {
        "Morgan":                     _h(morgan_ptbm),
        "Morgan+PFG_BDE_EGR+mol":     _h(np.hstack([morgan_ptbm, bde_mats_ptbm["BDE_EGR+mol"]])),
        "Morgan+mol_kf":              _h(np.hstack([morgan_ptbm, mol_kf_ptbm])),
    }

    for name, X in {**fsets_A, **fsets_B, **fsets_PTBM}.items():
        print(f"  {name}: {X.shape}")

    # ------------------------------------------------------------------
    # CV setup
    # ------------------------------------------------------------------
    outer_cv_A = RepeatedStratifiedKFold(
        n_splits=OUTER_SPLITS_A, n_repeats=OUTER_REPEATS, random_state=RANDOM_STATE
    )
    outer_cv_B    = StratifiedKFold(n_splits=OUTER_SPLITS, shuffle=True, random_state=RANDOM_STATE)
    outer_cv_PTBM = StratifiedKFold(n_splits=OUTER_SPLITS, shuffle=True, random_state=RANDOM_STATE)

    # ------------------------------------------------------------------
    # Run all experiments
    # ------------------------------------------------------------------
    all_rows: list[dict] = []
    n_done = [0]

    # Count total tasks for progress
    def count_tasks(smiles, Y, labels, fsets, min_pos):
        total = 0
        for ep_idx, ep in enumerate(labels):
            y_raw = Y[:, ep_idx]
            valid = ~np.isnan(y_raw.astype(float))
            y = y_raw[valid].astype(int)
            n_pos = int(y.sum())
            n_neg = int((y == 0).sum())
            if n_pos >= min_pos and n_neg >= min_pos:
                total += len(fsets)
        return total

    n_A    = count_tasks(smiles_pfas, Y_pfas,  labels_tc,   fsets_A,    MIN_POS_PFAS)
    n_B    = count_tasks(smiles_tc,   Y_tc,    labels_tc,   fsets_B,    MIN_POS_FULL)
    n_PTBM = count_tasks(smiles_ptbm, Y_ptbm,  labels_ptbm, fsets_PTBM, MIN_POS_FULL)
    n_total = n_A + n_B + n_PTBM

    print(f"\n=== Running nested CV — {n_total} tasks total ===")

    print("\n--- Exp A: PFAS subset (invitrodb) ---")
    run_experiment("Exp A", smiles_pfas, Y_pfas, labels_tc,   fsets_A,    outer_cv_A,
                   MIN_POS_PFAS, all_rows, n_done, n_total)

    print("\n--- Exp B: All chemicals (invitrodb) ---")
    run_experiment("Exp B", smiles_tc,   Y_tc,   labels_tc,   fsets_B,    outer_cv_B,
                   MIN_POS_FULL, all_rows, n_done, n_total)

    print("\n--- PTBM: Arp & Hale all chemicals ---")
    run_experiment("PTBM", smiles_ptbm,  Y_ptbm, labels_ptbm, fsets_PTBM, outer_cv_PTBM,
                   MIN_POS_FULL, all_rows, n_done, n_total)

    # ------------------------------------------------------------------
    # Save per-fold results
    # ------------------------------------------------------------------
    results_df = pd.DataFrame(all_rows)
    out_csv = DATA_DIR / "bde_egr_comparison_results.csv"
    results_df.to_csv(out_csv, index=False)
    print(f"\n[saved] {out_csv} ({len(results_df)} rows)")

    # ------------------------------------------------------------------
    # Bayesian comparisons
    # ------------------------------------------------------------------
    print("\n=== Bayesian Correlated t-tests ===")
    all_bayes: list[pd.DataFrame] = []

    comparison_pairs = [
        # (experiment, baseline_fset, compare_fsets, n_repeats, pair_label, title_suffix)
        (
            "Exp A", "PFG_EGR+mol",
            [f for f in ["PFG_BDE_EGR+mol", "PFG_BDE_EGR+mol+mol_kf", "PFG_EGR+mol+mol_kf"]
             if f in results_df["feature_set"].unique()],
            N_REPEATS_A,
            "ExpA_vs_EGR_mol",
            "Exp A (PFAS); BDE-EGR vs topological EGR",
        ),
        (
            "Exp A", "Morgan",
            [f for f in ["Morgan+mol_kf", "PFG_BDE_EGR+mol"]
             if f in results_df["feature_set"].unique()],
            N_REPEATS_A,
            "ExpA_vs_Morgan",
            "Exp A (PFAS); BDE-EGR features vs Morgan baseline",
        ),
        (
            "Exp B",
            "Morgan+PFG_EGR+mol" if "Morgan+PFG_EGR+mol" in results_df["feature_set"].unique()
            else "Morgan",
            [f for f in ["Morgan+PFG_BDE_EGR+mol", "Morgan+mol_kf"]
             if f in results_df["feature_set"].unique()],
            N_REPEATS_B,
            "ExpB_vs_EGR_mol",
            "Exp B (all invitrodb); BDE-EGR vs topological EGR",
        ),
        (
            "Exp B", "Morgan",
            [f for f in ["Morgan+PFG_BDE_EGR+mol", "Morgan+mol_kf"]
             if f in results_df["feature_set"].unique()],
            N_REPEATS_B,
            "ExpB_vs_Morgan",
            "Exp B (all invitrodb); BDE-EGR features vs Morgan baseline",
        ),
        (
            "PTBM", "Morgan",
            [f for f in ["Morgan+PFG_BDE_EGR+mol", "Morgan+mol_kf"]
             if f in results_df["feature_set"].unique()],
            N_REPEATS_B,
            "PTBM_vs_Morgan",
            "PTBM (Arp & Hale); BDE-EGR features vs Morgan baseline",
        ),
    ]

    for (exp, baseline, compare_fsets, n_rep, pair_label, title_suf) in comparison_pairs:
        if not compare_fsets:
            print(f"  [skip] {pair_label}: no compare feature sets available")
            continue

        bayes_df = run_bayesian_comparison(
            results_df, exp, baseline, compare_fsets, n_rep, pair_label
        )
        if bayes_df.empty:
            print(f"  [skip] {pair_label}: empty Bayesian result")
            continue

        all_bayes.append(bayes_df)

        out = DATA_DIR / f"bde_egr_bayesian_{pair_label}.csv"
        bayes_df.to_csv(out, index=False)
        print(f"[saved] {out.relative_to(DATA_DIR.parent)}")

        plot_bayesian_bar(bayes_df, pair_label, baseline,
                         title_suffix=title_suf)
        plot_bayesian_posteriors(
            results_df, bayes_df, exp, baseline, compare_fsets,
            n_rep, pair_label
        )

    # ------------------------------------------------------------------
    # Summary table
    # ------------------------------------------------------------------
    print("\n=== Summary: mean AUC per experiment × feature set ===")
    gb_df = results_df[results_df["model"] == "GradientBoosting"].copy()
    summary = (
        gb_df.groupby(["experiment", "feature_set", "endpoint"])["roc_auc"]
        .mean()
        .unstack("feature_set")
        .round(3)
    )
    pd.set_option("display.width", 200)
    pd.set_option("display.max_columns", 20)
    print(summary.to_string())

    elapsed = time.time() - t_start
    print(f"\nTotal runtime: {elapsed/60:.1f} min")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------
if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Benchmark BDE-EGR features vs topological EGR on ToxCast + PTBM."
    )
    parser.add_argument(
        "--db", "-d", default="invitrodb_v4_3",
        help="Database name prefix for the parquet file (default: invitrodb_v4_3).",
    )
    args = parser.parse_args()
    main(db=args.db)
