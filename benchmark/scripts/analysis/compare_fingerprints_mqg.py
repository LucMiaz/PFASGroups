#!/usr/bin/env python3
"""
MQG (Molecular Quantum Graph) descriptors combined with PFASGroups EGR.

Benchmark using the same ToxCast dataset and cross-validation as
compare_fingerprints_toxcast.py, adding molecule-level quantum graph
descriptors to the PFASGroups embedding.

MQG features (MQG_FP_SIZE + 10 + 24 columns total):
    MQG_FP_SIZE cols — raw BDE-weighted quantum eigenvalue spectrum, zero-padded
                       to a fixed length (normalize_length=True so eigenvalues
                       are size-free: k_n ~ n*pi for all molecules).
                       Computed via mol_to_graph_spectrum() for every molecule.
    10 cols          — size-free spectral ratios derived from the eigenvalues
                       (rel_width, rel_iqr, rel_top, rel_bot, skew_ev,
                        cv_ev, rel_sp_radius, gap_ratio, rel_gap2_gap1,
                        gap_iqr_ratio)
    24 cols          — BDE-based bond-length histogram bins

Feature sets compared (Experiment A — PFAS-containing chemicals):
    PFG_EGR_4X+mol          — all-halogen EGR + molecule metrics (baseline)
    PFG_EGR_4X+mol+MQG_EV   — + raw quantum eigenvalue vector
    PFG_EGR_4X+mol+MQG_RAT  — + 10 spectral ratios
    PFG_EGR_4X+mol+MQG_BL   — + 24 bond-length histogram bins
    PFG_EGR_4X+mol+MQG_full — + all 3 MQG blocks

Feature sets compared (Experiment B — full ToxCast library):
    Morgan+PFG_EGR_4X+mol             (baseline)
    Morgan+PFG_EGR_4X+mol+MQG_EV
    Morgan+PFG_EGR_4X+mol+MQG_RAT
    Morgan+PFG_EGR_4X+mol+MQG_BL
    Morgan+PFG_EGR_4X+mol+MQG_full

Experiment A: CF-containing chemicals (~808 compounds)
Experiment B: Full ToxCast library (~9014 compounds)

Usage:
    conda activate chem
    cd benchmark
    python scripts/compare_fingerprints_mqg.py [--workers N] [--fp-size N]

First run computes the full MQG spectrum for every ToxCast SMILES via
mol_to_graph_spectrum() and caches the result.  This can take several
hours depending on CPU count.  Subsequent runs load from cache (seconds).
The timeout per molecule is 300 s; molecules that exceed it receive a
zero vector (these are large, unusual structures unlikely to appear in
the PFAS-relevant endpoints).
"""

from __future__ import annotations

import argparse
import multiprocessing as mp
import os
import sys
import time
import warnings
from pathlib import Path

os.environ["PYTHONWARNINGS"] = "ignore"
warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from rdkit import Chem
from sklearn.ensemble import HistGradientBoostingClassifier, RandomForestClassifier
from sklearn.metrics import (
    average_precision_score, balanced_accuracy_score,
    matthews_corrcoef, roc_auc_score,
)
from sklearn.model_selection import (
    GridSearchCV, RepeatedStratifiedKFold, StratifiedKFold,
)

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
SCRIPT_DIR   = Path(__file__).resolve().parent
DATA_DIR     = SCRIPT_DIR.parents[1] / "data"
TEST_DATA    = SCRIPT_DIR.parents[1] / "test_data"
DATASET_PKL  = DATA_DIR / "toxcast_dataset.parquet"
MQG_CACHE    = DATA_DIR / "mqg_features_cache.npy"

# Add PFASGroups to path so build_pfg_matrices can be reused
sys.path.insert(0, str(SCRIPT_DIR.resolve().parents[2]))
# Also add MQG
MQG_PATH = "/home/luc/git/molecular_quantum_graph"
if MQG_PATH not in sys.path:
    sys.path.insert(0, MQG_PATH)

# Reuse helpers from the main benchmark script
from compare_fingerprints_toxcast import (  # noqa: E402
    MOL_METRICS, PFG_CONFIGS_4X, RANDOM_STATE,
    OUTER_SPLITS, OUTER_SPLITS_A, OUTER_REPEATS, INNER_SPLITS,
    MIN_POS_PFAS, MIN_POS_FULL,
    META_COLS, LABEL_COLS, MORGAN_NBITS, MORGAN_RADIUS,
    build_pfg_matrices, build_morgan, nested_cv_metrics, make_model_grids,
    grouped_bar_by_fset, summary_2x2, print_and_save_tables,
    load_tsv_fingerprints, radar_fset_comparison, delta_bar_plot,
)

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------
BDE_SINGLE = {
    (1, 1): 104.2, (1, 6): 99.5,  (1, 7): 93.4,  (1, 8): 110.6, (1, 9): 135.7,
    (6, 6): 83.1,  (6, 7): 69.1,  (6, 8): 84.9,  (6, 9): 108.0, (6, 14): 76.0,
    (6, 15): 65.0, (6, 16): 60.0, (6, 17): 79.0, (6, 35): 67.5, (6, 53): 52.0,
    (7, 7): 38.4,  (7, 8): 53.4,  (7, 9): 64.5,  (8, 8): 33.2,  (8, 9): 45.0,
    (9, 9): 36.9,  (16, 16): 54.0,
}
BDE_CC    = 83.1
BL_NBINS  = 24
BL_EDGES  = np.linspace(0.40, 2.80, BL_NBINS + 1)

TXPP_TSV = TEST_DATA / "TxP_PFAS_v1.tsv"

# ── MQG fingerprint layout ────────────────────────────────────────────────
MQG_FP_SIZE   = 16    # eigenvalue columns (zero-padded / truncated)
                      # most ToxCast molecules have ≤10 eigenvalues;
                      # 16 ensures no truncation for typical drug-like molecules
MQG_TIMEOUT_S = 300   # per-molecule timeout in seconds (generous; slow for
                      # large fused ring systems, but that is acceptable)

# Feature-block slices in the combined cache array (width = FP_SIZE + 10 + 24)
_EV_SLICE  = slice(0,              MQG_FP_SIZE)            # raw eigenvalues
_RAT_SLICE = slice(MQG_FP_SIZE,   MQG_FP_SIZE + 10)       # spectral ratios
_BL_SLICE  = slice(MQG_FP_SIZE + 10, MQG_FP_SIZE + 10 + BL_NBINS)  # BL hist
MQG_TOTAL  = MQG_FP_SIZE + 10 + BL_NBINS                  # 50 total

# Column name lists (for documentation)
EV_NAMES  = [f"mqg_ev_{i:02d}" for i in range(MQG_FP_SIZE)]
RAT_NAMES = [
    "qg_rel_width", "qg_rel_iqr", "qg_rel_top", "qg_rel_bot",
    "qg_skew_ev",   "qg_cv_ev",   "qg_rel_sp_radius",
    "qg_gap_ratio", "qg_rel_gap2_gap1", "qg_gap_iqr_ratio",
]
BL_NAMES  = [f"bl_{i:02d}" for i in range(BL_NBINS)]
MQG_NAMES = EV_NAMES + RAT_NAMES + BL_NAMES  # 50 total

# Visual palette
FSET_COLORS_A = {
    "PFG_EGR_4X+mol":          "#00838F",  # dark cyan  (baseline)
    "PFG_EGR_4X+mol+MQG_EV":   "#1565C0",  # dark blue  (+ raw eigenvalues)
    "PFG_EGR_4X+mol+MQG_RAT":  "#F57C00",  # orange     (+ spectral ratios)
    "PFG_EGR_4X+mol+MQG_BL":   "#558B2F",  # dark green (+ BL histogram)
    "PFG_EGR_4X+mol+MQG_full": "#E64A19",  # deep orange (all 3 MQG blocks)
}
FSET_COLORS_B = {
    "Morgan+PFG_EGR_4X+mol":          "#009688",  # teal
    "Morgan+PFG_EGR_4X+mol+MQG_EV":   "#1565C0",  # dark blue
    "Morgan+PFG_EGR_4X+mol+MQG_RAT":  "#F57C00",  # orange
    "Morgan+PFG_EGR_4X+mol+MQG_BL":   "#558B2F",  # dark green
    "Morgan+PFG_EGR_4X+mol+MQG_full": "#E64A19",  # deep orange
}


# ---------------------------------------------------------------------------
# Fast MQG feature extraction
# ---------------------------------------------------------------------------

def _bde_lookup(z1: int, z2: int) -> float:
    return BDE_SINGLE.get((min(z1, z2), max(z1, z2)), 75.0)


def bond_length_histogram(smi: str) -> np.ndarray:
    """24-bin bond-length histogram (BDE-based lengths), fast O(n_bonds)."""
    empty = np.zeros(BL_NBINS, dtype=np.float32)
    mol = Chem.MolFromSmiles(smi) if isinstance(smi, str) else None
    if mol is None or mol.GetNumBonds() == 0:
        return empty
    lengths = []
    for bond in mol.GetBonds():
        z1 = bond.GetBeginAtom().GetAtomicNum()
        z2 = bond.GetEndAtom().GetAtomicNum()
        order = bond.GetBondTypeAsDouble()
        bde = _bde_lookup(z1, z2) * max(1.0 + 0.3 * (order - 1.0), 1e-6)
        lengths.append(BDE_CC / bde)
    if not lengths:
        return empty
    hist, _ = np.histogram(lengths, bins=BL_EDGES)
    return (hist / len(lengths)).astype(np.float32)


def _spectrum_to_ratios(eigenvalues: list[float]) -> np.ndarray:
    """Compute 10 size-free spectral ratios from the eigenvalue list."""
    out = np.zeros(10, dtype=np.float32)
    pos = sorted([x for x in eigenvalues if x > 1e-10])
    if len(pos) < 2:
        return out
    ev = np.asarray(pos)
    mean_c = float(ev.mean()) or 1e-12
    std_c  = float(ev.std())  or 1e-12
    q10, q25, q50, q75, q90 = np.percentile(ev, [10, 25, 50, 75, 90])
    gaps = np.diff(np.concatenate([[0.0], ev]))
    first_gap  = float(gaps[0]) if len(gaps) > 0 else 0.0
    second_gap = float(gaps[1]) if len(gaps) > 1 else 0.0
    g25 = float(np.percentile(gaps, 25)) if len(gaps) > 1 else 0.0
    g75 = float(np.percentile(gaps, 75)) if len(gaps) > 1 else 0.0
    sp_radius = float(ev.max())
    out[0] = sp_radius / (float(ev.min()) or 1e-12)   # rel_width
    out[1] = q75 / (q25 or 1e-12)                      # rel_iqr
    out[2] = q90 / (q50 or 1e-12)                      # rel_top
    out[3] = q50 / (q10 or 1e-12)                      # rel_bot
    out[4] = (mean_c - q50) / std_c                    # skew_ev
    out[5] = std_c / mean_c                            # cv_ev
    out[6] = sp_radius / mean_c                        # rel_sp_radius
    out[7] = first_gap / mean_c                        # gap_ratio
    out[8] = second_gap / (first_gap or 1e-12)         # rel_gap2_gap1
    out[9] = g75 / (g25 or 1e-12)                      # gap_iqr_ratio
    return out


def _spectrum_to_ev_vec(eigenvalues: list[float], fp_size: int = MQG_FP_SIZE) -> np.ndarray:
    """
    Zero-pad (or truncate) the eigenvalue list to a fixed-length vector.

    Because mol_to_graph_spectrum is called with normalize_length=True, the
    eigenvalues satisfy k_n ~ n * pi regardless of molecule size — they are
    directly size-free and can be compared across molecules.
    Only positive eigenvalues are retained (the spectrum may contain a few
    near-zero values from the boundary condition solving).
    """
    pos = sorted([x for x in eigenvalues if x > 1e-10])
    arr = np.zeros(fp_size, dtype=np.float32)
    if pos:
        n = min(len(pos), fp_size)
        arr[:n] = np.asarray(pos[:n], dtype=np.float32)
    return arr


def _compute_one_mqg(smi: str, timeout: int = MQG_TIMEOUT_S) -> np.ndarray:
    """
    Compute the full MQG feature vector for one SMILES using
    mol_to_graph_spectrum() (BDE-weighted quantum graph, normalize_length=True).

    Returns (MQG_FP_SIZE + 10,) float32 array:
        [:MQG_FP_SIZE] — zero-padded eigenvalue vector (size-free)
        [MQG_FP_SIZE:] — 10 size-free spectral ratios

    Returns zeros on failure or timeout.
    """
    zeros = np.zeros(MQG_FP_SIZE + 10, dtype=np.float32)
    try:
        mol = Chem.MolFromSmiles(smi) if isinstance(smi, str) else None
        if mol is None or mol.GetNumBonds() == 0:
            return zeros
        from mqg.core import mol_to_graph_spectrum
        from mqg.schemes import create_weighting_scheme
        import signal

        def _handler(sig, frm):
            raise TimeoutError()
        signal.signal(signal.SIGALRM, _handler)
        signal.alarm(timeout)
        try:
            scheme = create_weighting_scheme("bde", mol=mol)
            # Full spectrum via numerical method with BDE-weighted edge lengths.
            # normalize_length=True fixes total graph length to 1 so eigenvalues
            # are directly size-free (k_n ~ n*pi for all molecules).
            ev = mol_to_graph_spectrum(mol, scheme, normalize_length=True)
            ev_vec = _spectrum_to_ev_vec(ev)
            ratios = _spectrum_to_ratios(ev)
            return np.concatenate([ev_vec, ratios])
        finally:
            signal.alarm(0)
    except Exception:
        return zeros


# Worker function for multiprocessing (must be module-level for pickle)
def _worker(args):
    idx, smi = args
    return idx, _compute_one_mqg(smi)


def build_mqg_features(
    smiles_list: list[str],
    cache_path: Path = MQG_CACHE,
    workers: int = 4,
    force: bool = False,
) -> np.ndarray:
    """
    Compute or load the full MQG feature matrix for every SMILES.

    MQG eigenvalues are computed via mol_to_graph_spectrum() with
    BDE-weighted edge lengths and normalize_length=True.  This is the
    exact same computation as used in the published MQG articles; it is
    slow for large molecules but produces physically meaningful features.
    The per-molecule timeout (MQG_TIMEOUT_S = 300 s) only affects the
    rare, very large structures; typical drug-like molecules finish in
    <10 s.

    Returns (n, MQG_TOTAL) float32 array where MQG_TOTAL = MQG_FP_SIZE + 10 + 24:
        cols [_EV_SLICE]  — MQG_FP_SIZE raw eigenvalues (zero-padded)
        cols [_RAT_SLICE] — 10 size-free spectral ratios
        cols [_BL_SLICE]  — 24 bond-length histogram bins (always computed
                            analytically; no timeout risk)
    """
    n = len(smiles_list)

    # ── Per-block cache files ────────────────────────────────────────────
    spec_cache = cache_path.with_suffix(".spectrum.npy")   # ev + ratios combined
    bl_cache   = cache_path.with_suffix(".bl.npy")

    spec_arr = None
    bl_arr   = None

    if not force:
        if spec_cache.exists():
            arr = np.load(spec_cache)
            if arr.shape == (n, MQG_FP_SIZE + 10):
                spec_arr = arr
                print(f"  [cache hit] MQG spectrum features: shape={arr.shape}")
        if bl_cache.exists():
            arr = np.load(bl_cache)
            if arr.shape[0] == n:
                bl_arr = arr
                print(f"  [cache hit] MQG bond-length histogram: shape={arr.shape}")

    # ── Bond-length histogram (fast) ─────────────────────────────────────
    if bl_arr is None:
        print(f"  Computing bond-length histograms for {n} SMILES …")
        bl_arr = np.vstack([bond_length_histogram(s) for s in smiles_list])
        np.save(bl_cache, bl_arr)
        print(f"  [cache saved] bl_arr shape={bl_arr.shape}")

    # ── Full quantum spectrum (slow, parallel, with per-molecule timeout) ─
    if spec_arr is None:
        print(f"  Computing MQG full quantum spectra for {n} SMILES …")
        print(f"  (mol_to_graph_spectrum, BDE-weighted, normalize_length=True)")
        print(f"  {workers} parallel workers, {MQG_TIMEOUT_S}s timeout per molecule.")
        print(f"  ETA rough estimate: several hours on first run. Safe to interrupt.")

        spec_arr  = np.zeros((n, MQG_FP_SIZE + 10), dtype=np.float32)
        args_list = [(i, smi) for i, smi in enumerate(smiles_list)]

        try:
            t0      = time.time()
            done    = 0
            skipped = 0
            with mp.Pool(processes=workers) as pool:
                for idx, vec in pool.imap_unordered(_worker, args_list, chunksize=4):
                    spec_arr[idx] = vec
                    done += 1
                    if np.all(vec == 0):
                        skipped += 1
                    if done % 100 == 0:
                        elapsed = time.time() - t0
                        rate    = done / max(elapsed, 1e-6)
                        eta_s   = (n - done) / rate
                        print(
                            f"  [{done}/{n}]  timeout/fail={skipped}  "
                            f"rate={rate:.1f}/s  ETA={eta_s/60:.0f}min"
                        )
        except KeyboardInterrupt:
            print("\n  [interrupted] Saving partial cache …")

        np.save(spec_cache, spec_arr)
        nonzero = int((spec_arr.sum(axis=1) != 0).sum())
        print(f"  [cache saved] spectrum shape={spec_arr.shape}  nonzero={nonzero}")

    # Combine: [ev_vec | ratios | bl_hist]
    return np.hstack([spec_arr, bl_arr]).astype(np.float32)  # (n, MQG_TOTAL)


# ---------------------------------------------------------------------------
# Helpers (re-used from compare_fingerprints_toxcast)
# ---------------------------------------------------------------------------

def _mean_std(df, metric, endpoints, fset):
    sub = df[df["feature_set"] == fset]
    means, stds = [], []
    for ep in endpoints:
        vals = sub.loc[sub["endpoint"] == ep, metric].dropna().values
        means.append(float(vals.mean()) if len(vals) else float("nan"))
        stds.append(float(vals.std())  if len(vals) > 1 else 0.0)
    return means, stds


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main(workers: int = 4) -> None:
    # ------------------------------------------------------------------ load
    print("Loading dataset …")
    tox   = pd.read_parquet(DATASET_PKL)
    smiles_all = tox["smiles"].tolist()
    n          = len(smiles_all)
    print(f"  ToxCast: {n:,} chemicals")

    # ── Pre-compute fingerprint matrices ────────────────────────────────
    print("\nLoading TxP_PFAS fingerprints …")
    X_txpp_all = load_tsv_fingerprints(TXPP_TSV)        # (n, 129)

    print("\nComputing / loading all-halogen PFASGroups embeddings …")
    pfg_4x = build_pfg_matrices(smiles_all, configs=PFG_CONFIGS_4X, halogens=None)

    print("\nGenerating Morgan fingerprints …")
    X_morgan_all = build_morgan(smiles_all)             # (n, 512)

    print("\nComputing / loading MQG features …")
    X_mqg = build_mqg_features(smiles_all, workers=workers)  # (n, MQG_TOTAL)
    X_ev   = X_mqg[:, _EV_SLICE].astype(np.float32)   # raw eigenvalue vector
    X_rat  = X_mqg[:, _RAT_SLICE].astype(np.float32)  # spectral ratios
    X_bl   = X_mqg[:, _BL_SLICE].astype(np.float32)   # bond-length histogram
    print(f"  MQG total shape={X_mqg.shape}  "
          f"(EV={X_ev.shape[1]}, RAT={X_rat.shape[1]}, BL={X_bl.shape[1]})"
          f"  EV nonzero rows={int((X_ev.sum(1)!=0).sum())}")

    # ─────────────────────────────────────────────────────────────────────
    # EXPERIMENT A  ─  CF-containing chemicals
    # ─────────────────────────────────────────────────────────────────────
    print("\n─── Experiment A: CF-containing chemicals ───")
    cf_mask = X_txpp_all.sum(axis=1) > 0
    print(f"  CF-containing: {int(cf_mask.sum())} / {n}")

    pfg4x_cf = {k: v[cf_mask].astype(np.float32) for k, v in pfg_4x.items()}
    X_ev_cf  = X_ev[cf_mask]
    X_rat_cf = X_rat[cf_mask]
    X_bl_cf  = X_bl[cf_mask]
    X_mqg_cf = X_mqg[cf_mask]
    tox_cf   = tox[cf_mask].reset_index(drop=True)

    base_cf = pfg4x_cf["EGR_4X+mol"]

    fsets_a: list[tuple[str, np.ndarray]] = [
        ("PFG_EGR_4X+mol",
         base_cf),
        ("PFG_EGR_4X+mol+MQG_EV",
         np.hstack([base_cf, X_ev_cf])),
        ("PFG_EGR_4X+mol+MQG_RAT",
         np.hstack([base_cf, X_rat_cf])),
        ("PFG_EGR_4X+mol+MQG_BL",
         np.hstack([base_cf, X_bl_cf])),
        ("PFG_EGR_4X+mol+MQG_full",
         np.hstack([base_cf, X_mqg_cf])),
    ]

    cv_small = RepeatedStratifiedKFold(
        n_splits=OUTER_SPLITS_A, n_repeats=OUTER_REPEATS, random_state=RANDOM_STATE
    )

    all_rows: list[dict] = []
    for ep in LABEL_COLS:
        if ep not in tox_cf.columns:
            continue
        mask  = tox_cf[ep].notna()
        y     = tox_cf.loc[mask, ep].values.astype(int)
        n_pos = int(y.sum())
        n_neg = len(y) - n_pos
        if n_pos < MIN_POS_PFAS or n_neg < MIN_POS_PFAS:
            continue
        print(f"  [A] {ep:<30} n={len(y):3d}  pos={n_pos:3d}")
        idx = mask.values
        for fset, X_full in fsets_a:
            all_rows.extend(nested_cv_metrics(X_full[idx], y, cv_small, ep, fset, "Exp A"))

    # ─────────────────────────────────────────────────────────────────────
    # EXPERIMENT B  ─  Full ToxCast library
    # ─────────────────────────────────────────────────────────────────────
    print("\n─── Experiment B: Full ToxCast library ───")
    cv_full = StratifiedKFold(n_splits=OUTER_SPLITS, shuffle=True, random_state=RANDOM_STATE)

    base_all = pfg_4x["EGR_4X+mol"].astype(np.float32)
    X_m_all  = X_morgan_all.astype(np.float32)

    fsets_b: list[tuple[str, np.ndarray]] = [
        ("Morgan+PFG_EGR_4X+mol",
         np.hstack([X_m_all, base_all])),
        ("Morgan+PFG_EGR_4X+mol+MQG_EV",
         np.hstack([X_m_all, base_all, X_ev])),
        ("Morgan+PFG_EGR_4X+mol+MQG_RAT",
         np.hstack([X_m_all, base_all, X_rat])),
        ("Morgan+PFG_EGR_4X+mol+MQG_BL",
         np.hstack([X_m_all, base_all, X_bl])),
        ("Morgan+PFG_EGR_4X+mol+MQG_full",
         np.hstack([X_m_all, base_all, X_mqg])),
    ]

    for ep in LABEL_COLS:
        if ep not in tox.columns:
            continue
        mask  = tox[ep].notna()
        y     = tox.loc[mask, ep].values.astype(int)
        n_pos = int(y.sum())
        n_neg = len(y) - n_pos
        if n_pos < MIN_POS_FULL or n_neg < MIN_POS_FULL:
            continue
        print(f"  [B] {ep:<30} n={len(y):5d}  pos={n_pos:5d}")
        idx = mask.values
        for fset, X_full in fsets_b:
            all_rows.extend(nested_cv_metrics(X_full[idx], y, cv_full, ep, fset, "Exp B"))

    # ─────────────────────────────────────────────────────────────────────
    # Save results
    # ─────────────────────────────────────────────────────────────────────
    results  = pd.DataFrame(all_rows)
    out_csv  = DATA_DIR / "mqg_comparison_results.csv"
    results.to_csv(out_csv, index=False)
    print(f"\n[saved] {out_csv.name}")

    print_and_save_tables(results, DATA_DIR)

    # ─────────────────────────────────────────────────────────────────────
    # Plots
    # ─────────────────────────────────────────────────────────────────────
    print("\nGenerating plots …")

    grouped_bar_by_fset(results, "Exp A", "roc_auc", "ROC-AUC",
                        FSET_COLORS_A,
                        DATA_DIR / "mqg_comparison_expA_auc.png", baseline=0.5)
    grouped_bar_by_fset(results, "Exp A", "avg_prec", "Average Precision",
                        FSET_COLORS_A,
                        DATA_DIR / "mqg_comparison_expA_ap.png", baseline=0.0)
    grouped_bar_by_fset(results, "Exp B", "roc_auc", "ROC-AUC",
                        FSET_COLORS_B,
                        DATA_DIR / "mqg_comparison_expB_auc.png", baseline=0.5)
    grouped_bar_by_fset(results, "Exp B", "avg_prec", "Average Precision",
                        FSET_COLORS_B,
                        DATA_DIR / "mqg_comparison_expB_ap.png", baseline=0.0)

    delta_bar_plot(results, "Exp A", "PFG_EGR_4X+mol", "PFG_EGR_4X+mol+MQG_EV",
                   "roc_auc", "ROC-AUC",
                   DATA_DIR / "mqg_delta_EV_A.png")
    delta_bar_plot(results, "Exp A", "PFG_EGR_4X+mol", "PFG_EGR_4X+mol+MQG_RAT",
                   "roc_auc", "ROC-AUC",
                   DATA_DIR / "mqg_delta_RAT_A.png")
    delta_bar_plot(results, "Exp A", "PFG_EGR_4X+mol", "PFG_EGR_4X+mol+MQG_full",
                   "roc_auc", "ROC-AUC",
                   DATA_DIR / "mqg_delta_full_A.png")
    delta_bar_plot(results, "Exp B", "Morgan+PFG_EGR_4X+mol",
                   "Morgan+PFG_EGR_4X+mol+MQG_EV",
                   "roc_auc", "ROC-AUC",
                   DATA_DIR / "mqg_delta_EV_B.png")
    delta_bar_plot(results, "Exp B", "Morgan+PFG_EGR_4X+mol",
                   "Morgan+PFG_EGR_4X+mol+MQG_RAT",
                   "roc_auc", "ROC-AUC",
                   DATA_DIR / "mqg_delta_RAT_B.png")
    delta_bar_plot(results, "Exp B", "Morgan+PFG_EGR_4X+mol",
                   "Morgan+PFG_EGR_4X+mol+MQG_full",
                   "roc_auc", "ROC-AUC",
                   DATA_DIR / "mqg_delta_full_B.png")

    radar_fset_comparison(results, "Exp A", list(FSET_COLORS_A.keys()),
                          FSET_COLORS_A, "roc_auc", "ROC-AUC",
                          DATA_DIR / "mqg_radar_A_auc.png")
    radar_fset_comparison(results, "Exp B", list(FSET_COLORS_B.keys()),
                          FSET_COLORS_B, "roc_auc", "ROC-AUC",
                          DATA_DIR / "mqg_radar_B_auc.png")

    print("\nDone. Outputs:")
    for f in sorted(DATA_DIR.glob("mqg_*")):
        print(f"  {f.name}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="MQG + PFASGroups combined benchmark on ToxCast"
    )
    parser.add_argument(
        "--workers", type=int, default=max(1, mp.cpu_count() - 2),
        help="Number of parallel worker processes for MQG computation "
             "(default: n_cpu − 2)"
    )
    args = parser.parse_args()
    main(workers=args.workers)
