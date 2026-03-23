#!/usr/bin/env python3
"""
Patch script: rerun Experiment B for the 4X feature sets only
(Morgan+PFG_EGR_4X+mol and ToxPrint+PFG_EGR_4X+mol) using the corrected
build_pfg_matrices_4x() that produces 117×4=468 component columns.

The existing toxcast_comparison_results.csv is loaded, the stale 4X rows for
Exp B are dropped, new rows are computed and merged back, and both
toxcast_comparison_results.csv and toxcast_comparison_summary.csv are
re-saved.

Usage
─────
    conda activate chem
    cd benchmark
    python scripts/patch_4x_expB.py
"""

from __future__ import annotations

import os
os.environ["PYTHONWARNINGS"] = "ignore"

import warnings
warnings.filterwarnings("ignore")

from pathlib import Path
import numpy as np
import pandas as pd

# ---------------------------------------------------------------------------
# Re-use everything from compare_fingerprints_toxcast
# ---------------------------------------------------------------------------
import sys
sys.path.insert(0, str(Path(__file__).resolve().parent))
from compare_fingerprints_toxcast import (
    DATA_DIR, TEST_DATA, DATASET_PKL,
    TXPP_TSV, TOXPRINT_TSV,
    LABEL_COLS, MIN_POS_FULL, OUTER_SPLITS, RANDOM_STATE,
    load_tsv_fingerprints, build_pfg_matrices, build_pfg_matrices_4x,
    build_morgan, nested_cv_metrics, print_and_save_tables,
)
from sklearn.model_selection import StratifiedKFold

# Feature sets in Exp B that use the 4X embeddings
FSETS_4X_B = {"Morgan+PFG_EGR_4X+mol", "ToxPrint+PFG_EGR_4X+mol"}

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    print("Loading dataset …")
    tox = pd.read_parquet(DATASET_PKL)
    smiles_all = tox["smiles"].tolist()

    # ── Build required fingerprint matrices ─────────────────────────
    print("\nLoading pre-computed TSV fingerprints …")
    X_txpp_all     = load_tsv_fingerprints(TXPP_TSV)
    X_toxprint_all = load_tsv_fingerprints(TOXPRINT_TSV)

    print("\nComputing / loading EGR+mol embedding (for ToxPrint combo) …")
    pfg = build_pfg_matrices(smiles_all)

    print("\nComputing corrected 4X embeddings …")
    pfg_4x = build_pfg_matrices_4x(smiles_all)
    print(f"  EGR_4X      shape: {pfg_4x['EGR_4X'].shape}")
    print(f"  EGR_4X+mol  shape: {pfg_4x['EGR_4X+mol'].shape}")

    print("\nGenerating Morgan fingerprints …")
    X_morgan_all = build_morgan(smiles_all)

    X_morgan_f32         = X_morgan_all.astype(np.float32)
    X_toxprint_pfgegr4x  = np.hstack([X_toxprint_all,
                                       pfg_4x["EGR_4X+mol"].astype(np.float32)])

    # ── Exp B CV setup ───────────────────────────────────────────────
    cv_full = StratifiedKFold(n_splits=OUTER_SPLITS, shuffle=True,
                              random_state=RANDOM_STATE)

    fsets_4x_b: list[tuple[str, np.ndarray]] = [
        ("Morgan+PFG_EGR_4X+mol",   np.hstack([X_morgan_f32,
                                                pfg_4x["EGR_4X+mol"].astype(np.float32)])),
        ("ToxPrint+PFG_EGR_4X+mol", X_toxprint_pfgegr4x),
    ]

    print("\n─── Experiment B (4X feature sets only) ───")
    new_rows: list[dict] = []
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

        for fset, X_full in fsets_4x_b:
            new_rows.extend(
                nested_cv_metrics(X_full[idx], y, cv_full, ep, fset, "Exp B")
            )

    # ── Merge with existing results ──────────────────────────────────
    results_csv = DATA_DIR / "toxcast_comparison_results.csv"
    existing = pd.read_csv(results_csv)
    print(f"\nExisting rows: {len(existing)}")

    # Drop stale 4X rows for Exp B
    mask_stale = (existing["experiment"] == "Exp B") & (existing["feature_set"].isin(FSETS_4X_B))
    dropped = mask_stale.sum()
    existing = existing[~mask_stale].copy()
    print(f"Dropped stale rows: {dropped}")

    new_df  = pd.DataFrame(new_rows)
    results = pd.concat([existing, new_df], ignore_index=True)
    print(f"New total rows: {len(results)}")

    results.to_csv(results_csv, index=False)
    print(f"[saved] {results_csv.name}")

    # Re-generate summary CSV
    print_and_save_tables(results, DATA_DIR)
    print("\nDone.")


if __name__ == "__main__":
    main()
