#!/usr/bin/env python3
"""
Compare fingerprint representations on ToxCast endpoint hit-calls.

Experiment A  ─ CF-containing chemicals only (~808 compounds with ≥1
                non-zero TxP_PFAS bit)
  Features compared  (PFAS-specific fingerprint evaluation):
    • TxP_PFAS      – 129 PFAS-structural bits (Richard et al. 2023, v1)
    • PFG_binary    – 115 binary group-presence bits (this work)
    • PFG_binary+mol– PFG_binary + 10 molecule-wide graph metrics (125 cols)
    • PFG_EGR       – effective graph resistance only (115 cols, no binary)
    • PFG_EGR+mol   – PFG_EGR + 10 molecule-wide metrics (125 cols)

Experiment B  ─ Full ToxCast library (n ≈ 9 014)
  Features compared  (general + PFAS augmentation):
    • ToxPrint+TxP_PFAS       – 729 + 129 = 858 bits
    • Morgan+PFG_binary       – ECFP4 512 + PFG binary 115 (627 cols)
    • Morgan+PFG_binary+mol   – ECFP4 512 + PFG binary+mol 125 (637 cols)
    • Morgan+PFG_EGR          – ECFP4 512 + PFG EGR 115 (627 cols)
    • Morgan+PFG_EGR+mol      – ECFP4 512 + PFG EGR+mol 125 (637 cols)
    • Morgan                  – ECFP4 512-bit baseline

Cross-validation strategy  (nested CV)
  Outer loop: StratifiedKFold(5) — unbiased performance estimate
  Inner loop: StratifiedKFold(3) via GridSearchCV — hyperparameter tuning
  For Exp A (small n): outer = RepeatedStratifiedKFold(3 splits, 5 repeats)

Models
  • RandomForestClassifier   (tuned: n_estimators, max_features)
  • HistGradientBoostingClassifier (tuned: max_depth, learning_rate)

Usage
─────
    conda activate chem
    cd benchmark
    python scripts/compare_fingerprints_toxcast.py
"""

from __future__ import annotations

import os
os.environ["PYTHONWARNINGS"] = "ignore"

import warnings
warnings.filterwarnings("ignore")
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from rdkit import Chem
from rdkit.Chem import AllChem
from sklearn.ensemble import (
    HistGradientBoostingClassifier,
    RandomForestClassifier,
)
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
SCRIPT_DIR   = Path(__file__).resolve().parent
DATA_DIR     = SCRIPT_DIR.parent / "data"
TEST_DATA    = SCRIPT_DIR.parent / "test_data"
TXPP_TSV     = TEST_DATA / "TxP_PFAS_v1.tsv"
TOXPRINT_TSV = TEST_DATA / "toxprint_V2.tsv"
DATASET_PKL  = DATA_DIR / "toxcast_dataset.parquet"

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------
RANDOM_STATE  = 42
MORGAN_NBITS  = 512
MORGAN_RADIUS = 2

# Outer / inner CV sizes
OUTER_SPLITS   = 5   # folds for unbiased performance estimate
OUTER_REPEATS  = 5   # repeats for Exp A (small n)
OUTER_SPLITS_A = 3   # folds for Exp A outer loop (small n)
INNER_SPLITS   = 3   # folds for hyperparameter search

META_COLS  = ["chid", "casn", "chnm", "dsstox_substance_id", "smiles"]
LABEL_COLS = [
    "AR_antagonist", "AR_agonist",
    "ERa_antagonist", "ERa_agonist",
    "AhR_agonist", "Aromatase_antagonist",
    "TR_antagonist", "DT40_genotoxicity",
    "MMP_ratio",
    "CYP2D6_antagonist", "CYP2C19_antagonist", "CYP2C9_antagonist",
    "CYP3A4_antagonist", "p53_ratio", "Caspase3_HEPG2",
]

# Minimum positives (and negatives) needed in a task to run CV:
MIN_POS_PFAS = 10   # Exp A (small PFAS subset, ~200 chems)
MIN_POS_FULL = 30   # Exp B (full library, ~9 000 chems)

# Molecule-wide metrics appended after component columns
MOL_METRICS = [
    "n_components", "total_size", "mean_size", "max_size",
    "mean_branching", "max_branching", "mean_eccentricity",
    "max_diameter", "mean_component_fraction", "max_component_fraction",
]

# PFASGroups embedding configurations: name → to_array kwargs
PFG_CONFIGS = {
    "binary":     dict(component_metrics=["binary"]),
    "binary+mol": dict(component_metrics=["binary"], molecule_metrics=MOL_METRICS),
    "EGR":        dict(component_metrics=["effective_graph_resistance"]),
    "EGR+mol":    dict(component_metrics=["effective_graph_resistance"],
                       molecule_metrics=MOL_METRICS),
}

# Visual palette
MODEL_COLORS = {
    "RandomForest":            "#2196F3",  # blue
    "GradientBoosting":        "#FF5722",  # deep-orange
}
FSET_COLORS_A = {
    "TxP_PFAS":       "#E91E63",  # pink
    "PFG_binary":     "#CE93D8",  # light purple
    "PFG_binary+mol": "#9C27B0",  # purple
    "PFG_EGR":        "#80CBC4",  # light teal
    "PFG_EGR+mol":    "#009688",  # teal
}
FSET_COLORS_B = {
    "ToxPrint+TxP_PFAS":  "#607D8B",  # blue-grey
    "Morgan+PFG_binary":     "#CE93D8",  # light purple
    "Morgan+PFG_binary+mol": "#9C27B0",  # purple
    "Morgan+PFG_EGR":        "#80CBC4",  # light teal
    "Morgan+PFG_EGR+mol":    "#009688",  # teal
    "Morgan":                "#2196F3",  # blue
}

# ---------------------------------------------------------------------------
# Nested CV: hyperparameter grids
# ---------------------------------------------------------------------------
# Inner GridSearchCV scoring
_INNER_SCORING = "roc_auc"

def _rf_grid(random_state: int) -> GridSearchCV:
    base = RandomForestClassifier(
        class_weight="balanced",
        random_state=random_state,
        n_jobs=-1,
    )
    param_grid = {
        "n_estimators": [200, 400],
        "max_features": ["sqrt", 0.2],
    }
    return GridSearchCV(
        base, param_grid,
        cv=StratifiedKFold(INNER_SPLITS, shuffle=True, random_state=random_state),
        scoring=_INNER_SCORING,
        n_jobs=-1,
        refit=True,
    )


def _gb_grid(random_state: int) -> GridSearchCV:
    # HistGradientBoosting handles missing values natively (NaN-safe)
    # and is much faster than GradientBoostingClassifier for larger datasets
    base = HistGradientBoostingClassifier(
        class_weight="balanced",
        random_state=random_state,
    )
    param_grid = {
        "max_depth":     [3, 5],
        "learning_rate": [0.05, 0.1],
    }
    return GridSearchCV(
        base, param_grid,
        cv=StratifiedKFold(INNER_SPLITS, shuffle=True, random_state=random_state),
        scoring=_INNER_SCORING,
        n_jobs=-1,
        refit=True,
    )


def make_model_grids(random_state: int = RANDOM_STATE) -> dict[str, GridSearchCV]:
    """Return fresh GridSearchCV wrappers for each model (one call per fold)."""
    return {
        "RandomForest":     _rf_grid(random_state),
        "GradientBoosting": _gb_grid(random_state),
    }


# ---------------------------------------------------------------------------
# Feature builders
# ---------------------------------------------------------------------------

def load_tsv_fingerprints(path: Path) -> np.ndarray:
    """
    Load a pre-computed fingerprint TSV (rows = Structure #1…#9014, first col
    is the row-index label).  Returns (9014, n_bits) uint8 array in row order.
    """
    df = pd.read_csv(path, sep="\t", index_col=0)
    X  = df.values.astype(np.uint8)
    print(f"  Loaded {path.name}: shape={X.shape}  nonzero rows={int((X.sum(1)>0).sum())}")
    return X


def build_pfg_matrices(
    smiles_list: list[str],
    configs: dict[str, dict] | None = None,
    halogens: str = "F",
    cache_dir: Path | None = DATA_DIR,
) -> dict[str, np.ndarray]:
    """
    Compute PFASGroups embedding matrices for several metric configs.

    Handles InChIKey duplicates: unique molecules are parsed once, then
    mapped back so the returned matrices are row-aligned with *smiles_list*.

    Returns dict mapping config name → (n, n_cols) float64 array.
    """
    import hashlib
    import sys as _sys

    if configs is None:
        configs = PFG_CONFIGS

    _sys.path.insert(0, str(SCRIPT_DIR.resolve().parents[1]))
    from PFASGroups import parse_smiles as _parse  # noqa: PLC0415

    n = len(smiles_list)

    # ── Check per-config caches ──────────────────────────────────────
    cache_hash = hashlib.sha256("\n".join(smiles_list).encode()).hexdigest()[:16]
    cache_meta = cache_dir / "pfg_cache_meta.txt" if cache_dir else None

    if cache_dir and cache_meta:
        cached_hash = cache_meta.read_text().strip() if cache_meta.exists() else ""
        if cached_hash == cache_hash:
            all_ok = True
            matrices: dict[str, np.ndarray] = {}
            for name in configs:
                p = cache_dir / f"pfg_{name}_cache.npy"
                if p.exists():
                    arr = np.load(p)
                    if arr.shape[0] == n:
                        matrices[name] = arr
                        continue
                all_ok = False
                break
            if all_ok and len(matrices) == len(configs):
                for name, arr in matrices.items():
                    print(f"  [cache hit] PFG {name}: shape={arr.shape}")
                return matrices
        print("  [cache miss/stale] Recomputing PFASGroups embeddings …")

    # ── Deduplicate by InChIKey ──────────────────────────────────────
    print(f"  Computing InChIKeys for {n} SMILES …")
    inchikeys: list[str | None] = []
    for smi in smiles_list:
        mol = Chem.MolFromSmiles(smi)
        inchikeys.append(Chem.inchi.MolToInchiKey(mol) if mol else None)

    unique_order: dict[str, int] = {}   # inchikey → index in unique list
    unique_smiles: list[str] = []
    for smi, ik in zip(smiles_list, inchikeys):
        if ik and ik not in unique_order:
            unique_order[ik] = len(unique_smiles)
            unique_smiles.append(smi)

    n_unique = len(unique_smiles)
    n_invalid = sum(1 for ik in inchikeys if ik is None)
    print(f"  Unique: {n_unique}  |  duplicates: {n - n_unique - n_invalid}  |  invalid: {n_invalid}")

    # ── Parse unique SMILES once ─────────────────────────────────────
    print(f"  Parsing {n_unique} unique SMILES with PFASGroups …")
    result = _parse(unique_smiles, halogens=halogens, progress=True)

    # Use same halogens filter for group list so to_array produces 115 cols
    # (includes Cl/Br/I group slots for F-mode, which encode mixed-halogen PFAS)
    from PFASGroups.getter import get_compiled_HalogenGroups as _get_groups_fn
    _pfas_groups = _get_groups_fn(halogens=halogens)

    # Build inchikey → result-row-index mapping
    result_ik_map: dict[str, int] = {}
    for idx, emb in enumerate(result):
        ik = emb.get("inchikey", "")
        if ik:
            result_ik_map[ik] = idx
    print(f"  Parsed result: {len(result)} entries, {len(result_ik_map)} with valid InChIKey")

    # Row-index lookup: input row → result row (-1 = no match)
    input_to_result = np.full(n, -1, dtype=int)
    for i, ik in enumerate(inchikeys):
        if ik and ik in result_ik_map:
            input_to_result[i] = result_ik_map[ik]
    valid_mask = input_to_result >= 0
    valid_indices = input_to_result[valid_mask]

    # ── Build matrices for each config ───────────────────────────────
    matrices = {}
    for name, kwargs in configs.items():
        print(f"  Computing PFG {name} embedding …")
        arr_result = np.asarray(result.to_array(**kwargs, progress=True, pfas_groups=_pfas_groups))
        arr = np.zeros((n, arr_result.shape[1]), dtype=np.float64)
        arr[valid_mask] = arr_result[valid_indices]
        matrices[name] = arr
        nonzero = int((np.abs(arr).sum(axis=1) > 0).sum())
        print(f"    shape={arr.shape}  nonzero rows={nonzero}")

    # ── Save caches ──────────────────────────────────────────────────
    if cache_dir:
        cache_dir.mkdir(parents=True, exist_ok=True)
        for name, arr in matrices.items():
            np.save(cache_dir / f"pfg_{name}_cache.npy", arr)
        cache_meta.write_text(cache_hash)  # type: ignore[union-attr]
        print(f"  [cache saved] {len(matrices)} matrices")

    return matrices


def build_morgan(smiles_list: list[str]) -> np.ndarray:
    """Return (n, MORGAN_NBITS) uint8 Morgan fingerprint matrix."""
    X      = np.zeros((len(smiles_list), MORGAN_NBITS), dtype=np.uint8)
    failed = 0
    for i, smi in enumerate(smiles_list):
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            failed += 1
            continue
        fp = AllChem.GetMorganFingerprintAsBitVect(mol, MORGAN_RADIUS, MORGAN_NBITS)
        X[i] = np.array(fp)
    if failed:
        print(f"  Warning: {failed} SMILES failed Morgan generation (left as zeros)")
    print(f"  Morgan{MORGAN_NBITS}: shape={X.shape}  nonzero rows={int((X.sum(1)>0).sum())}")
    return X


# ---------------------------------------------------------------------------
# Nested cross-validation helper
# ---------------------------------------------------------------------------

def nested_cv_metrics(
    X: np.ndarray,
    y: np.ndarray,
    outer_cv,
    label: str,
    fset: str,
    experiment: str,
) -> list[dict]:
    """
    Nested CV: outer loop evaluates; inner loop (via GridSearchCV) tunes.

    The outer loop is a RepeatedStratifiedKFold (Exp A) or StratifiedKFold
    (Exp B).  Each split gets fresh GridSearchCV objects so the inner search
    is independent for every outer fold.

    Returns a list of per-fold result dicts (one per outer-fold × model).
    """
    rows: list[dict] = []

    # Group splits by repeat so we can average per-repeat for repeated KFold
    splits = list(outer_cv.split(X, y))
    has_repeats = hasattr(outer_cv, "n_repeats")
    # RepeatedStratifiedKFold stores n_splits in cvargs; StratifiedKFold exposes it directly
    n_splits_per_repeat = outer_cv.cvargs["n_splits"] if has_repeats else outer_cv.n_splits

    for split_i, (tr, te) in enumerate(splits):
        repeat_idx = split_i // n_splits_per_repeat if has_repeats else 0
        fold_idx   = split_i % n_splits_per_repeat

        X_tr, X_te = X[tr], X[te]
        y_tr, y_te = y[tr], y[te]

        # Skip degenerate folds (only one class in train or test)
        if len(np.unique(y_tr)) < 2 or len(np.unique(y_te)) < 2:
            continue

        for mname, grid_search in make_model_grids().items():
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                grid_search.fit(X_tr, y_tr)
                prob = grid_search.predict_proba(X_te)[:, 1]
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
                "model":       mname,
                "repeat":      repeat_idx,
                "fold":        fold_idx,
                "n_total":     len(y),
                "n_pos":       int(y.sum()),
                "pos_rate":    round(float(y.mean()), 4),
                "roc_auc":     round(float(auc), 4),
                "avg_prec":    round(float(ap),  4),
                "mcc":         round(float(matthews_corrcoef(y_te, pred)), 4),
                "bal_acc":     round(float(balanced_accuracy_score(y_te, pred)), 4),
                "best_params": str(grid_search.best_params_),
            })

    return rows


# ---------------------------------------------------------------------------
# Plotting helpers
# ---------------------------------------------------------------------------

def _mean_std(df: pd.DataFrame, metric: str, endpoints: list[str], fset: str):
    sub = df[df["feature_set"] == fset]
    means, stds = [], []
    for ep in endpoints:
        vals = sub.loc[sub["endpoint"] == ep, metric].dropna().values
        means.append(float(vals.mean()) if len(vals) else float("nan"))
        stds.append(float(vals.std())  if len(vals) > 1 else 0.0)
    return means, stds


def grouped_bar_by_fset(
    results: pd.DataFrame,
    experiment: str,
    metric: str,
    metric_label: str,
    fset_colors: dict[str, str],
    path: Path,
    baseline: float = 0.5,
) -> None:
    """
    Grouped bar: x = endpoint, groups = feature_set.
    Bar height = mean over folds × models; error bar = std over folds × models.
    """
    df = results[results["experiment"] == experiment].copy()
    if df.empty:
        print(f"  [skip plot] no data for {experiment}")
        return

    endpoints = sorted(df["endpoint"].unique())
    fsets_in  = [f for f in fset_colors if f in df["feature_set"].unique()]
    n_ep, n_fs = len(endpoints), len(fsets_in)
    width = 0.8 / max(n_fs, 1)
    x     = np.arange(n_ep)

    fig, ax = plt.subplots(figsize=(max(8, n_ep * 1.1 + 2), 5))
    ax.axhline(baseline, color="grey", linestyle="--", linewidth=0.8, zorder=0)

    for k, fset in enumerate(fsets_in):
        means, stds = _mean_std(df, metric, endpoints, fset)
        offset = (k - (n_fs - 1) / 2) * width
        ax.bar(
            x + offset, means, width * 0.9,
            color=fset_colors[fset], label=fset,
            yerr=stds, capsize=3, error_kw={"elinewidth": 0.8},
        )

    ax.set_xticks(x)
    ax.set_xticklabels(endpoints, rotation=45, ha="right", fontsize=8)
    ax.set_ylabel(metric_label)
    ax.set_ylim(max(0, baseline - 0.05), 1.02)
    ax.set_title(
        f"ToxCast – {metric_label}  [{experiment}]"
        f"\n(mean ± SD across folds × models, nested CV)"
    )
    ax.legend(fontsize=9)
    plt.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)
    print(f"  [plot] {path.name}")


def model_heatmap(
    results: pd.DataFrame,
    metric: str,
    metric_label: str,
    path: Path,
) -> None:
    """
    Heatmap: rows = model × feature_set, columns = endpoint.
    One panel per experiment.
    """
    exps = [e for e in ["Exp A", "Exp B"] if e in results["experiment"].values]
    if not exps:
        return

    fig, axes = plt.subplots(1, len(exps), figsize=(7 * len(exps), 5), squeeze=False)
    for col, exp in enumerate(exps):
        ax = axes[0][col]
        sub = results[results["experiment"] == exp]
        pivot = sub.pivot_table(
            index=["model", "feature_set"],
            columns="endpoint",
            values=metric,
            aggfunc="mean",
        )
        if pivot.empty:
            ax.set_visible(False)
            continue
        im = ax.imshow(pivot.values, aspect="auto", cmap="RdYlGn",
                       vmin=0.5 if "auc" in metric else 0.0, vmax=1.0)
        ax.set_xticks(range(len(pivot.columns)))
        ax.set_xticklabels(pivot.columns, rotation=45, ha="right", fontsize=7)
        ax.set_yticks(range(len(pivot.index)))
        ax.set_yticklabels(
            [f"{m} / {fs}" for m, fs in pivot.index], fontsize=7
        )
        for i in range(pivot.shape[0]):
            for j in range(pivot.shape[1]):
                v = pivot.values[i, j]
                if not np.isnan(v):
                    ax.text(j, i, f"{v:.2f}", ha="center", va="center",
                            fontsize=6, color="black")
        plt.colorbar(im, ax=ax, fraction=0.03)
        ax.set_title(f"{exp} – {metric_label}", fontsize=10)

    fig.suptitle(
        f"Model × Feature-Set Heatmap — {metric_label}\n(nested CV mean over folds)",
        fontsize=11,
    )
    plt.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)
    print(f"  [plot] {path.name}")


def scatter_two_fsets(
    results: pd.DataFrame,
    experiment: str,
    fset_x: str,
    fset_y: str,
    metric: str,
    metric_label: str,
    path: Path,
) -> None:
    """
    Scatter plot: fset_y vs fset_x per endpoint, one point per model.
    Diagonal = equal performance line.
    """
    df    = results[results["experiment"] == experiment].copy()
    means = df.groupby(["feature_set", "endpoint", "model"])[metric].mean().reset_index()
    pivot = means.pivot_table(
        index=["endpoint", "model"],
        columns="feature_set",
        values=metric,
    ).reset_index()

    if fset_x not in pivot.columns or fset_y not in pivot.columns:
        print(f"  [skip scatter] need both '{fset_x}' and '{fset_y}' in results")
        return

    fig, ax = plt.subplots(figsize=(6, 6))
    models = pivot["model"].unique()
    for mname in models:
        sub = pivot[pivot["model"] == mname].dropna(subset=[fset_x, fset_y])
        ax.scatter(
            sub[fset_x], sub[fset_y],
            color=MODEL_COLORS.get(mname, "#888888"),
            label=mname, alpha=0.85, s=60, zorder=3,
        )
        for _, row in sub.iterrows():
            ax.annotate(
                row["endpoint"], (row[fset_x], row[fset_y]),
                fontsize=6, textcoords="offset points", xytext=(3, 3), color="grey",
            )

    lo = min(ax.get_xlim()[0], ax.get_ylim()[0])
    hi = max(ax.get_xlim()[1], ax.get_ylim()[1])
    ax.plot([lo, hi], [lo, hi], "k--", linewidth=0.8, zorder=0, label="equal")
    ax.set_xlabel(f"{fset_x}  {metric_label}")
    ax.set_ylabel(f"{fset_y}  {metric_label}")
    ax.set_title(
        f"{fset_y} vs {fset_x} \u2014 {metric_label}"
        f"\n({experiment}, each point = one endpoint \u00d7 model)"
    )
    ax.legend(fontsize=9)
    plt.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)
    print(f"  [plot] {path.name}")


def model_comparison_bar(
    results: pd.DataFrame,
    experiment: str,
    metric: str,
    metric_label: str,
    path: Path,
) -> None:
    """
    Bar chart comparing models side-by-side, averaged over endpoints and feature-sets.
    """
    df = results[results["experiment"] == experiment].copy()
    if df.empty:
        return
    agg = df.groupby(["model", "endpoint"])[metric].mean().reset_index()
    models = list(MODEL_COLORS.keys())
    models_in = [m for m in models if m in agg["model"].unique()]
    endpoints = sorted(agg["endpoint"].unique())

    n_ep = len(endpoints)
    n_m  = len(models_in)
    width = 0.8 / max(n_m, 1)
    x = np.arange(n_ep)

    fig, ax = plt.subplots(figsize=(max(8, n_ep * 1.1 + 2), 5))
    ax.axhline(0.5 if "auc" in metric else 0.0, color="grey",
               linestyle="--", linewidth=0.8, zorder=0)

    for k, mname in enumerate(models_in):
        sub = agg[agg["model"] == mname]
        means, stds = [], []
        for ep in endpoints:
            vals = df.loc[(df["model"] == mname) & (df["endpoint"] == ep), metric].dropna().values
            means.append(float(vals.mean()) if len(vals) else float("nan"))
            stds.append(float(vals.std())  if len(vals) > 1 else 0.0)
        offset = (k - (n_m - 1) / 2) * width
        ax.bar(x + offset, means, width * 0.9,
               color=MODEL_COLORS[mname], label=mname,
               yerr=stds, capsize=3, error_kw={"elinewidth": 0.8})

    ax.set_xticks(x)
    ax.set_xticklabels(endpoints, rotation=45, ha="right", fontsize=8)
    ax.set_ylabel(metric_label)
    ax.set_ylim(0, 1.05)
    ax.set_title(f"Model comparison — {metric_label}  [{experiment}]\n"
                 "(mean ± SD over folds × feature-sets, nested CV)")
    ax.legend(fontsize=9)
    plt.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)
    print(f"  [plot] {path.name}")


def summary_2x2(
    results: pd.DataFrame,
    fset_colors_a: dict,
    fset_colors_b: dict,
    path: Path,
) -> None:
    """Four-panel overview: Exp A/B × ROC-AUC/AP (mean across models)."""
    panels = [
        ("Exp A", "roc_auc",  "ROC-AUC",         fset_colors_a, 0.5),
        ("Exp A", "avg_prec", "Avg Precision",     fset_colors_a, 0.0),
        ("Exp B", "roc_auc",  "ROC-AUC",         fset_colors_b, 0.5),
        ("Exp B", "avg_prec", "Avg Precision",     fset_colors_b, 0.0),
    ]
    fig, axes = plt.subplots(2, 2, figsize=(18, 10))
    for ax, (exp, metric, mlabel, fcolors, bsl) in zip(axes.flat, panels):
        df = results[results["experiment"] == exp]
        if df.empty:
            ax.set_visible(False)
            continue
        endpoints = sorted(df["endpoint"].unique())
        fsets_in  = [f for f in fcolors if f in df["feature_set"].unique()]
        n_ep = len(endpoints)
        n_fs = len(fsets_in)
        width = 0.8 / max(n_fs, 1)
        x = np.arange(n_ep)
        ax.axhline(bsl, color="grey", linestyle="--", linewidth=0.8)
        for k, fset in enumerate(fsets_in):
            means, stds = _mean_std(df, metric, endpoints, fset)
            offset = (k - (n_fs - 1) / 2) * width
            ax.bar(x + offset, means, width * 0.9, color=fcolors[fset], label=fset,
                   yerr=stds, capsize=2, error_kw={"elinewidth": 0.7})
        ax.set_xticks(x)
        ax.set_xticklabels(endpoints, rotation=45, ha="right", fontsize=7)
        ax.set_ylim(max(0, bsl - 0.05), 1.02)
        ax.set_ylabel(mlabel, fontsize=9)
        ax.set_title(f"{exp} – {mlabel}", fontsize=10)
        ax.legend(fontsize=7, ncol=2)

    fig.suptitle(
        "ToxCast Fingerprint Comparison  (nested CV, mean ± SD over folds × models)",
        fontsize=12,
    )
    plt.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)
    print(f"  [plot] {path.name}")


def delta_bar_plot(
    results: pd.DataFrame,
    experiment: str,
    fset_base: str,
    fset_aug: str,
    metric: str,
    metric_label: str,
    path: Path,
) -> None:
    """
    Bar chart showing delta (augmented - baseline) per endpoint per model.
    Positive = augmented fingerprint is better.
    """
    df = results[results["experiment"] == experiment].copy()
    if df.empty:
        return
    means_ep = (
        df.groupby(["feature_set", "endpoint", "model"])[metric]
        .mean()
        .reset_index()
    )
    pivot = means_ep.pivot_table(
        index=["endpoint", "model"], columns="feature_set", values=metric
    ).reset_index()
    if fset_base not in pivot.columns or fset_aug not in pivot.columns:
        print(f"  [skip delta] need both '{fset_base}' and '{fset_aug}'")
        return
    pivot["delta"] = pivot[fset_aug] - pivot[fset_base]
    endpoints = sorted(pivot["endpoint"].unique())
    models = sorted(pivot["model"].unique())
    n_ep = len(endpoints)
    n_m = len(models)
    width = 0.8 / max(n_m, 1)
    x = np.arange(n_ep)

    fig, ax = plt.subplots(figsize=(max(8, n_ep * 1.1 + 2), 5))
    ax.axhline(0, color="grey", linestyle="--", linewidth=0.8, zorder=0)
    for k, mname in enumerate(models):
        sub = pivot[pivot["model"] == mname]
        deltas = [float(sub.loc[sub["endpoint"] == ep, "delta"].values[0])
                  if ep in sub["endpoint"].values else 0.0
                  for ep in endpoints]
        offset = (k - (n_m - 1) / 2) * width
        colors = [("#4CAF50" if d >= 0 else "#F44336") for d in deltas]
        bars = ax.bar(x + offset, deltas, width * 0.9, color=colors,
                      label=mname, edgecolor="white", linewidth=0.3)

    ax.set_xticks(x)
    ax.set_xticklabels(endpoints, rotation=45, ha="right", fontsize=8)
    ax.set_ylabel(f"\u0394 {metric_label}  ({fset_aug} \u2212 {fset_base})")
    ax.set_title(
        f"Additive value of PFASGroups — \u0394 {metric_label}\n"
        f"({fset_aug} vs {fset_base}, [{experiment}], green=improvement)"
    )
    ax.legend(fontsize=9, title="Model")
    plt.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)
    print(f"  [plot] {path.name}")


def radar_fset_comparison(
    results: pd.DataFrame,
    experiment: str,
    fsets: list[str],
    fset_colors: dict[str, str],
    metric: str,
    metric_label: str,
    path: Path,
) -> None:
    """Radar / spider chart: one axis per endpoint, one line per feature-set."""
    df = results[results["experiment"] == experiment].copy()
    if df.empty:
        return
    means = df.groupby(["feature_set", "endpoint"])[metric].mean().reset_index()
    pivot = means.pivot_table(index="endpoint", columns="feature_set", values=metric)
    fsets_in = [f for f in fsets if f in pivot.columns]
    if len(fsets_in) < 2:
        return
    labels = list(pivot.index)
    n = len(labels)
    angles = np.linspace(0, 2 * np.pi, n, endpoint=False).tolist()
    angles.append(angles[0])

    fig, ax = plt.subplots(figsize=(8, 8), subplot_kw=dict(polar=True))
    for fs in fsets_in:
        vals = pivot[fs].values.tolist()
        vals.append(vals[0])
        ax.plot(angles, vals, "o-", linewidth=1.5, markersize=4,
                label=fs, color=fset_colors.get(fs, "#888888"))
        ax.fill(angles, vals, alpha=0.08, color=fset_colors.get(fs, "#888888"))
    ax.set_xticks(angles[:-1])
    ax.set_xticklabels(labels, fontsize=7)
    ax.set_ylim(0, 1.05)
    ax.set_title(f"{experiment} — {metric_label} per endpoint\n(mean over models & folds)",
                 fontsize=10, pad=20)
    ax.legend(loc="upper right", bbox_to_anchor=(1.3, 1.1), fontsize=8)
    plt.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)
    print(f"  [plot] {path.name}")


# ---------------------------------------------------------------------------
# Summary tables
# ---------------------------------------------------------------------------

def print_and_save_tables(results: pd.DataFrame, out_dir: Path) -> None:
    """Print per-experiment pivot tables and save a summary CSV."""
    summary_rows = []

    for exp in ["Exp A", "Exp B"]:
        sub = results[results["experiment"] == exp]
        if sub.empty:
            continue

        # ── ROC-AUC pivot (feature_set × model rows, endpoint columns) ──
        pivot_auc = sub.pivot_table(
            index=["feature_set", "model"],
            columns="endpoint",
            values="roc_auc",
            aggfunc="mean",
        )
        print(f"\n{'='*70}")
        print(f"{exp} — mean ROC-AUC  (nested CV, averaged over folds)")
        print(f"{'='*70}")
        print(pivot_auc.to_string(float_format="%.3f"))

        # ── Avg Precision pivot ──
        pivot_ap = sub.pivot_table(
            index=["feature_set", "model"],
            columns="endpoint",
            values="avg_prec",
            aggfunc="mean",
        )
        print(f"\n{exp} — mean Average Precision")
        print(pivot_ap.to_string(float_format="%.3f"))

        # ── Macro-averaged summary (mean ± SD across endpoints) ──
        print(f"\n{exp} — macro-average across endpoints")
        macro = (
            sub.groupby(["feature_set", "model"])[["roc_auc", "avg_prec", "mcc", "bal_acc"]]
            .agg(["mean", "std"])
        )
        print(macro.to_string(float_format="%.3f"))

        # Accumulate for CSV
        for (fs, mname), row in macro.iterrows():
            summary_rows.append({
                "experiment":        exp,
                "feature_set":       fs,
                "model":             mname,
                "roc_auc_mean":      round(row[("roc_auc", "mean")], 4),
                "roc_auc_std":       round(row[("roc_auc", "std")],  4),
                "avg_prec_mean":     round(row[("avg_prec", "mean")], 4),
                "avg_prec_std":      round(row[("avg_prec", "std")],  4),
                "mcc_mean":          round(row[("mcc", "mean")], 4),
                "bal_acc_mean":      round(row[("bal_acc", "mean")], 4),
            })

    # ── PFG variants vs TxP_PFAS comparison (Exp A) ──
    exp_a = results[results["experiment"] == "Exp A"]
    pfg_fsets = [f for f in exp_a["feature_set"].unique() if f.startswith("PFG_")]
    if "TxP_PFAS" in exp_a["feature_set"].values and pfg_fsets:
        print("\n" + "="*60)
        print("PFG variants vs TxP_PFAS — macro-average (Exp A)")
        print("="*60)
        fsets_show = ["TxP_PFAS"] + sorted(pfg_fsets)
        comp = (
            exp_a[exp_a["feature_set"].isin(fsets_show)]
            .groupby(["feature_set", "model"])[["roc_auc", "avg_prec"]]
            .mean()
        )
        print(comp.to_string(float_format="%.3f"))

        ep_means = (
            exp_a[exp_a["feature_set"].isin(fsets_show)]
            .groupby(["feature_set", "endpoint"])["roc_auc"]
            .mean()
            .unstack("feature_set")
        )
        for pf in pfg_fsets:
            if pf in ep_means and "TxP_PFAS" in ep_means:
                ep_means[f"delta_{pf}"] = ep_means[pf] - ep_means["TxP_PFAS"]
        print("\nPer-endpoint delta (PFG \u2212 TxP_PFAS) ROC-AUC:")
        print(ep_means.to_string(float_format="%.3f"))

    # ── Morgan+PFG vs ToxPrint+TxP_PFAS comparison (Exp B) ──
    exp_b = results[results["experiment"] == "Exp B"]
    if not exp_b.empty:
        comp_b = exp_b.groupby(["feature_set", "model"])[["roc_auc", "avg_prec"]].mean()
        print("\n" + "="*60)
        print("All feature sets — macro-average (Exp B, full library)")
        print("="*60)
        print(comp_b.to_string(float_format="%.3f"))

        # Morgan+PFG_EGR+mol vs Morgan delta
        if "Morgan" in exp_b["feature_set"].values and "Morgan+PFG_EGR+mol" in exp_b["feature_set"].values:
            ep_m = (
                exp_b[exp_b["feature_set"].isin(["Morgan", "Morgan+PFG_EGR+mol"])]
                .groupby(["feature_set", "endpoint"])["roc_auc"]
                .mean()
                .unstack("feature_set")
            )
            if "Morgan" in ep_m and "Morgan+PFG_EGR+mol" in ep_m:
                ep_m["delta"] = ep_m["Morgan+PFG_EGR+mol"] - ep_m["Morgan"]
                print("\nPer-endpoint delta (Morgan+PFG_EGR+mol \u2212 Morgan) ROC-AUC:")
                print(ep_m.to_string(float_format="%.3f"))
                n_improved = (ep_m["delta"] > 0).sum()
                n_degraded = (ep_m["delta"] < 0).sum()
                print(f"\nEndpoints improved: {n_improved}  |  degraded: {n_degraded}")
                print(f"Mean delta: {ep_m['delta'].mean():.4f}")

    # Save summary CSV
    summary_df = pd.DataFrame(summary_rows)
    csv_path = out_dir / "toxcast_comparison_summary.csv"
    summary_df.to_csv(csv_path, index=False)
    print(f"\n[saved] {csv_path.name}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    # ------------------------------------------------------------------ load
    print("Loading dataset …")
    tox = pd.read_parquet(DATASET_PKL)

    pfg_bin_cols = [c for c in tox.columns if c not in META_COLS + LABEL_COLS]
    print(f"  ToxCast:  {len(tox):,} chemicals  |  {len(pfg_bin_cols)} PFASGroups binary columns")

    smiles_all = tox["smiles"].tolist()

    # ── Pre-compute fingerprint matrices ────────────────────────────────
    print("\nLoading pre-computed TSV fingerprints …")
    X_txpp_all     = load_tsv_fingerprints(TXPP_TSV)        # (9014, 129)
    X_toxprint_all = load_tsv_fingerprints(TOXPRINT_TSV)    # (9014, 729)

    print("\nComputing / loading PFASGroups embeddings …")
    pfg = build_pfg_matrices(smiles_all)
    # pfg keys: 'binary' (n,115), 'binary+mol' (n,125), 'EGR' (n,115), 'EGR+mol' (n,125)

    print("\nGenerating Morgan fingerprints …")
    X_morgan_all = build_morgan(smiles_all)                  # (9014, 512)

    X_toxprint_txpp_all = np.hstack([X_toxprint_all, X_txpp_all]).astype(np.float32)  # (9014, 858)

    # ─────────────────────────────────────────────────────────────────────
    # EXPERIMENT A  ─  CF-containing chemicals (non-zero TxP_PFAS rows)
    # ─────────────────────────────────────────────────────────────────────
    print("\n─── Experiment A: CF-containing chemicals ───")
    cf_mask = (X_txpp_all.sum(axis=1) > 0)   # ~808 rows
    print(f"  CF-containing (TxP_PFAS non-zero): {int(cf_mask.sum())} / {len(cf_mask)}")

    X_txpp_cf      = X_txpp_all[cf_mask].astype(np.float32)
    pfg_cf = {name: arr[cf_mask].astype(np.float32) for name, arr in pfg.items()}
    tox_cf = tox[cf_mask].reset_index(drop=True)

    cv_small = RepeatedStratifiedKFold(
        n_splits=OUTER_SPLITS_A, n_repeats=OUTER_REPEATS, random_state=RANDOM_STATE
    )

    # Build Exp A feature-set list
    fsets_a: list[tuple[str, np.ndarray]] = [
        ("TxP_PFAS",       X_txpp_cf),
        ("PFG_binary",     pfg_cf["binary"]),
        ("PFG_binary+mol", pfg_cf["binary+mol"]),
        ("PFG_EGR",        pfg_cf["EGR"]),
        ("PFG_EGR+mol",    pfg_cf["EGR+mol"]),
    ]

    all_rows: list[dict] = []
    for ep in LABEL_COLS:
        if ep not in tox_cf.columns:
            continue
        mask  = tox_cf[ep].notna()
        y     = tox_cf.loc[mask, ep].values.astype(int)
        n_pos = int(y.sum())
        n_neg = len(y) - n_pos
        if n_pos < MIN_POS_PFAS or n_neg < MIN_POS_PFAS:
            print(f"  [skip A] {ep:<30} n={len(y)} pos={n_pos} neg={n_neg}")
            continue
        print(f"  [A] {ep:<30} n={len(y):3d}  pos={n_pos:3d}  ({100*n_pos/len(y):.0f}%)")
        idx = mask.values

        for fset, X_full in fsets_a:
            all_rows.extend(nested_cv_metrics(X_full[idx], y, cv_small, ep, fset, "Exp A"))

    # ─────────────────────────────────────────────────────────────────────
    # EXPERIMENT B  ─  Full ToxCast library
    # ─────────────────────────────────────────────────────────────────────
    print("\n─── Experiment B: Full ToxCast library ───")
    cv_full = StratifiedKFold(n_splits=OUTER_SPLITS, shuffle=True, random_state=RANDOM_STATE)

    # Build Exp B feature-set list
    X_morgan_f32 = X_morgan_all.astype(np.float32)
    fsets_b: list[tuple[str, np.ndarray]] = [
        ("ToxPrint+TxP_PFAS",  X_toxprint_txpp_all),
        ("Morgan+PFG_binary",     np.hstack([X_morgan_f32, pfg["binary"].astype(np.float32)])),
        ("Morgan+PFG_binary+mol", np.hstack([X_morgan_f32, pfg["binary+mol"].astype(np.float32)])),
        ("Morgan+PFG_EGR",        np.hstack([X_morgan_f32, pfg["EGR"].astype(np.float32)])),
        ("Morgan+PFG_EGR+mol",    np.hstack([X_morgan_f32, pfg["EGR+mol"].astype(np.float32)])),
        ("Morgan",              X_morgan_f32),
    ]

    for ep in LABEL_COLS:
        if ep not in tox.columns:
            continue
        mask  = tox[ep].notna()
        y     = tox.loc[mask, ep].values.astype(int)
        n_pos = int(y.sum())
        n_neg = len(y) - n_pos
        if n_pos < MIN_POS_FULL or n_neg < MIN_POS_FULL:
            print(f"  [skip B] {ep:<30} n={len(y)} pos={n_pos} neg={n_neg}")
            continue
        print(f"  [B] {ep:<30} n={len(y):5d}  pos={n_pos:5d}  ({100*n_pos/len(y):.0f}%)")
        idx = mask.values

        for fset, X_full in fsets_b:
            all_rows.extend(nested_cv_metrics(X_full[idx], y, cv_full, ep, fset, "Exp B"))

    # ─────────────────────────────────────────────────────────────────────
    # Save raw results
    # ─────────────────────────────────────────────────────────────────────
    results = pd.DataFrame(all_rows)
    out_csv = DATA_DIR / "toxcast_comparison_results.csv"
    results.to_csv(out_csv, index=False)
    print(f"\n[saved] {out_csv.name}")

    # ─────────────────────────────────────────────────────────────────────
    # Tables & analysis
    # ─────────────────────────────────────────────────────────────────────
    print_and_save_tables(results, DATA_DIR)

    # ─────────────────────────────────────────────────────────────────────
    # Plots
    # ─────────────────────────────────────────────────────────────────────
    print("\nGenerating plots …")

    # 1. Fingerprint comparison bar charts (mean over models)
    grouped_bar_by_fset(results, "Exp A", "roc_auc",  "ROC-AUC",
                        FSET_COLORS_A, DATA_DIR / "toxcast_comparison_expA_auc.png", baseline=0.5)
    grouped_bar_by_fset(results, "Exp A", "avg_prec", "Average Precision",
                        FSET_COLORS_A, DATA_DIR / "toxcast_comparison_expA_ap.png",  baseline=0.0)
    grouped_bar_by_fset(results, "Exp B", "roc_auc",  "ROC-AUC",
                        FSET_COLORS_B, DATA_DIR / "toxcast_comparison_expB_auc.png", baseline=0.5)
    grouped_bar_by_fset(results, "Exp B", "avg_prec", "Average Precision",
                        FSET_COLORS_B, DATA_DIR / "toxcast_comparison_expB_ap.png",  baseline=0.0)

    # 2. Heatmaps (model × feature_set rows, endpoint columns)
    model_heatmap(results, "roc_auc",  "ROC-AUC",
                  DATA_DIR / "toxcast_comparison_heatmap_auc.png")
    model_heatmap(results, "avg_prec", "Average Precision",
                  DATA_DIR / "toxcast_comparison_heatmap_ap.png")

    # 3. PFG variants vs TxP_PFAS scatter (Exp A)
    for pfg_name in ["PFG_binary", "PFG_EGR", "PFG_EGR+mol"]:
        safe = pfg_name.replace("+", "_plus_")
        scatter_two_fsets(results, "Exp A", "TxP_PFAS", pfg_name, "roc_auc", "ROC-AUC",
                          DATA_DIR / f"toxcast_comparison_{safe}_vs_txp_auc.png")
        scatter_two_fsets(results, "Exp A", "TxP_PFAS", pfg_name, "avg_prec", "Average Precision",
                          DATA_DIR / f"toxcast_comparison_{safe}_vs_txp_ap.png")

    # 4. Delta bar plots — PFG embedding richness comparison (Exp A)
    delta_bar_plot(results, "Exp A", "PFG_binary", "PFG_EGR",
                   "roc_auc", "ROC-AUC",
                   DATA_DIR / "toxcast_comparison_delta_EGR_vs_binary_A.png")
    delta_bar_plot(results, "Exp A", "PFG_EGR", "PFG_EGR+mol",
                   "roc_auc", "ROC-AUC",
                   DATA_DIR / "toxcast_comparison_delta_mol_boost_A.png")

    # 5. Morgan+PFG vs Morgan delta (Exp B)
    delta_bar_plot(results, "Exp B", "Morgan", "Morgan+PFG_binary",
                   "roc_auc", "ROC-AUC",
                   DATA_DIR / "toxcast_comparison_delta_morgan_pfg_binary_B.png")
    delta_bar_plot(results, "Exp B", "Morgan", "Morgan+PFG_EGR+mol",
                   "roc_auc", "ROC-AUC",
                   DATA_DIR / "toxcast_comparison_delta_morgan_pfg_EGR_B.png")

    # 6. Morgan+PFG vs ToxPrint+TxP_PFAS scatter (Exp B)
    scatter_two_fsets(results, "Exp B", "ToxPrint+TxP_PFAS", "Morgan+PFG_EGR+mol",
                      "roc_auc", "ROC-AUC",
                      DATA_DIR / "toxcast_comparison_morgan_pfg_vs_toxprint_auc.png")
    scatter_two_fsets(results, "Exp B", "ToxPrint+TxP_PFAS", "Morgan+PFG_EGR+mol",
                      "avg_prec", "Average Precision",
                      DATA_DIR / "toxcast_comparison_morgan_pfg_vs_toxprint_ap.png")

    # 7. Radar charts — per-endpoint comparison
    radar_fset_comparison(
        results, "Exp A",
        list(FSET_COLORS_A.keys()),
        FSET_COLORS_A, "roc_auc", "ROC-AUC",
        DATA_DIR / "toxcast_comparison_radar_A_auc.png",
    )
    radar_fset_comparison(
        results, "Exp B",
        list(FSET_COLORS_B.keys()),
        FSET_COLORS_B, "roc_auc", "ROC-AUC",
        DATA_DIR / "toxcast_comparison_radar_B_auc.png",
    )

    # 8. Combined 2×2 overview
    summary_2x2(results, FSET_COLORS_A, FSET_COLORS_B,
                DATA_DIR / "toxcast_comparison_summary.png")

    print("\nDone. Outputs:")
    for f in sorted(DATA_DIR.glob("toxcast_comparison_*")):
        print(f"  {f.name}")


# ---------------------------------------------------------------------------
if __name__ == "__main__":
    main()
