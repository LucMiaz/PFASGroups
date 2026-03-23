#!/usr/bin/env python3
"""
Patch script: add Experiment A results for TxP_PFAS+g51g52_total.

This feature set augments the 129-bit TxP_PFAS fingerprint with 2 extra columns:
  - Group 51 (perhalogenated alkyl)  total_component
  - Group 52 (polyhalogenated alkyl) total_component
producing a 131-column hybrid.

Usage
─────
    conda activate chem
    cd benchmark
    python scripts/patch_txpp_g51g52_expA.py
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
    TXPP_TSV,
    LABEL_COLS, MIN_POS_PFAS, OUTER_SPLITS_A, OUTER_REPEATS, RANDOM_STATE,
    load_tsv_fingerprints, build_pfg_matrices,
    nested_cv_metrics, print_and_save_tables,
)
from sklearn.model_selection import RepeatedStratifiedKFold

NEW_FSET = "TxP_PFAS+g51g52_total"


def main() -> None:
    results_csv = DATA_DIR / "toxcast_comparison_results.csv"

    print("Loading dataset …")
    tox = pd.read_parquet(DATASET_PKL)
    smiles_all = tox["smiles"].tolist()

    print("\nLoading TxP_PFAS fingerprints …")
    X_txpp_all = load_tsv_fingerprints(TXPP_TSV)  # (n, 129)

    # ── Build the 2-column total_component embedding for groups 51 & 52 ─
    print("\nComputing g51g52_total embedding (2 cols) …")
    pfg = build_pfg_matrices(
        smiles_all,
        configs={"g51g52_total": dict(component_metrics=["total_component"],
                                      selected_group_ids=[51, 52])},
    )
    X_g5152 = pfg["g51g52_total"]  # (n, 2)
    print(f"  g51g52_total shape: {X_g5152.shape}")

    # ── Exp A: CF-containing subset ──────────────────────────────────
    cf_mask = (X_txpp_all.sum(axis=1) > 0)
    print(f"\nCF-containing chemicals: {int(cf_mask.sum())} / {len(cf_mask)}")

    X_hybrid_cf = np.hstack([
        X_txpp_all[cf_mask].astype(np.float32),
        X_g5152[cf_mask].astype(np.float32),
    ])  # (n_cf, 131)
    print(f"  Hybrid TxP+g51g52_total shape: {X_hybrid_cf.shape}")

    tox_cf = tox[cf_mask].reset_index(drop=True)
    cv_small = RepeatedStratifiedKFold(
        n_splits=OUTER_SPLITS_A, n_repeats=OUTER_REPEATS, random_state=RANDOM_STATE
    )

    print(f"\n─── Experiment A — {NEW_FSET} ───")
    new_rows: list[dict] = []
    for ep in LABEL_COLS:
        if ep not in tox_cf.columns:
            continue
        mask  = tox_cf[ep].notna()
        y     = tox_cf.loc[mask, ep].values.astype(int)
        n_pos = int(y.sum())
        n_neg = len(y) - n_pos
        if n_pos < MIN_POS_PFAS or n_neg < MIN_POS_PFAS:
            print(f"  [skip] {ep:<30} n={len(y)} pos={n_pos} neg={n_neg}")
            continue
        print(f"  [A] {ep:<30} n={len(y):3d}  pos={n_pos:3d}  ({100*n_pos/len(y):.0f}%)")
        idx = mask.values
        new_rows.extend(
            nested_cv_metrics(X_hybrid_cf[idx], y, cv_small, ep, NEW_FSET, "Exp A")
        )

    # ── Merge with existing results ──────────────────────────────────
    existing = pd.read_csv(results_csv)
    print(f"\nExisting rows: {len(existing)}")

    # Drop any stale rows for this feature set
    stale = (existing["experiment"] == "Exp A") & (existing["feature_set"] == NEW_FSET)
    existing = existing[~stale].copy()
    print(f"Dropped stale rows: {int(stale.sum())}")

    new_df  = pd.DataFrame(new_rows)
    results = pd.concat([existing, new_df], ignore_index=True)
    print(f"New total rows: {len(results)}  (added {len(new_df)})")

    results.to_csv(results_csv, index=False)
    print(f"[saved] {results_csv.name}")

    print_and_save_tables(results, DATA_DIR)
    print("\nDone.")


if __name__ == "__main__":
    main()
