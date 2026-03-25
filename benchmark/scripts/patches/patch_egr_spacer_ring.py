#!/usr/bin/env python3
"""
Patch script: add benchmark results for PFG_EGR+spacer+ring+mol.

Feature set added
──────────────────
  Experiment A:
    PFG_EGR+spacer+ring+mol   (EGR + n_spacer + ring_size + 10 mol metrics)

  Experiment B (ToxPrint augmentation):
    ToxPrint+PFG_EGR+spacer+ring+mol

Usage
─────
    conda activate chem
    cd benchmark
    python scripts/patch_egr_spacer_ring.py
"""

from __future__ import annotations

import os
os.environ["PYTHONWARNINGS"] = "ignore"

import warnings
warnings.filterwarnings("ignore")

from pathlib import Path
import numpy as np
import pandas as pd
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent))
from compare_fingerprints_toxcast import (
    DATA_DIR, TEST_DATA, DATASET_PKL,
    TXPP_TSV, TOXPRINT_TSV, MOL_METRICS,
    LABEL_COLS, MIN_POS_PFAS, MIN_POS_FULL,
    OUTER_SPLITS, OUTER_SPLITS_A, OUTER_REPEATS, RANDOM_STATE,
    load_tsv_fingerprints, build_pfg_matrices,
    nested_cv_metrics, print_and_save_tables,
)
from sklearn.model_selection import RepeatedStratifiedKFold, StratifiedKFold

NEW_CONFIG_KEY = "EGR+spacer+ring+mol"
NEW_CONFIG = {
    NEW_CONFIG_KEY: dict(
        component_metrics=["effective_graph_resistance", "n_spacer", "ring_size"],
        molecule_metrics=MOL_METRICS,
    )
}

FSET_A = f"PFG_{NEW_CONFIG_KEY}"
FSET_B = f"ToxPrint+PFG_{NEW_CONFIG_KEY}"


def main() -> None:
    results_csv = DATA_DIR / "toxcast_comparison_results.csv"

    print("Loading dataset …")
    tox = pd.read_parquet(DATASET_PKL)
    smiles_all = tox["smiles"].tolist()

    print("\nLoading pre-computed TSV fingerprints …")
    X_txpp_all     = load_tsv_fingerprints(TXPP_TSV)       # (9014, 129)
    X_toxprint_all = load_tsv_fingerprints(TOXPRINT_TSV)   # (9014, 729)

    print(f"\nComputing PFASGroups embedding: {NEW_CONFIG_KEY} …")
    pfg = build_pfg_matrices(smiles_all, configs=NEW_CONFIG)
    X_pfg = pfg[NEW_CONFIG_KEY]
    print(f"  shape={X_pfg.shape}")

    # ── Experiment A — CF-containing subset ─────────────────────────────
    cf_mask = (X_txpp_all.sum(axis=1) > 0)
    print(f"\nCF-containing chemicals: {int(cf_mask.sum())} / {len(cf_mask)}")
    tox_cf  = tox[cf_mask].reset_index(drop=True)
    X_pfg_cf = X_pfg[cf_mask].astype(np.float32)

    cv_small = RepeatedStratifiedKFold(
        n_splits=OUTER_SPLITS_A, n_repeats=OUTER_REPEATS, random_state=RANDOM_STATE
    )

    new_rows: list[dict] = []

    print(f"\n─── Experiment A — {FSET_A} (shape {X_pfg_cf.shape}) ───")
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
        new_rows.extend(nested_cv_metrics(X_pfg_cf[idx], y, cv_small, ep, FSET_A, "Exp A"))

    # ── Experiment B — full ToxCast library ─────────────────────────────
    X_toxprint_pfg = np.hstack(
        [X_toxprint_all, X_pfg.astype(np.float32)]
    )
    cv_full = StratifiedKFold(n_splits=OUTER_SPLITS, shuffle=True, random_state=RANDOM_STATE)

    print(f"\n─── Experiment B — {FSET_B} (shape {X_toxprint_pfg.shape}) ───")
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
        new_rows.extend(
            nested_cv_metrics(X_toxprint_pfg[idx], y, cv_full, ep, FSET_B, "Exp B")
        )

    # ── Merge with existing results ──────────────────────────────────────
    existing = pd.read_csv(results_csv)
    print(f"\nExisting rows: {len(existing)}")

    stale_mask = (existing["feature_set"].isin([FSET_A, FSET_B]))
    existing = existing[~stale_mask].copy()
    print(f"Dropped stale rows: {int(stale_mask.sum())}")

    new_df   = pd.DataFrame(new_rows)
    results  = pd.concat([existing, new_df], ignore_index=True)
    print(f"New total rows: {len(results)}  (added {len(new_df)})")

    results.to_csv(results_csv, index=False)
    print(f"[saved] {results_csv.name}")

    print_and_save_tables(results, DATA_DIR)
    print("\nDone.")


if __name__ == "__main__":
    main()
