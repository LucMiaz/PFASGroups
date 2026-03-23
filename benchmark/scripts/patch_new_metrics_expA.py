#!/usr/bin/env python3
"""
Patch script: add Experiment A results for total_component and
min_dist_to_barycenter PFASGroups embeddings (with and without mol metrics).

Feature sets added
──────────────────
  PFG_total_component            117 cols  (total component size per group)
  PFG_total_component+mol        127 cols  (+ 10 molecule-wide metrics)
  PFG_min_dist_to_barycenter     117 cols  (min dist to EGR-space centroid)
  PFG_min_dist_to_barycenter+mol 127 cols  (+ 10 molecule-wide metrics)

Usage
─────
    conda activate chem
    cd benchmark
    python scripts/patch_new_metrics_expA.py
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
    TXPP_TSV, MOL_METRICS,
    LABEL_COLS, MIN_POS_PFAS, OUTER_SPLITS_A, OUTER_REPEATS, RANDOM_STATE,
    load_tsv_fingerprints, build_pfg_matrices,
    nested_cv_metrics, print_and_save_tables,
)
from sklearn.model_selection import RepeatedStratifiedKFold

NEW_CONFIGS = {
    "total_component":            dict(component_metrics=["total_component"]),
    "total_component+mol":        dict(component_metrics=["total_component"],
                                       molecule_metrics=MOL_METRICS),
    "min_dist_to_barycenter":     dict(component_metrics=["min_dist_to_barycenter"]),
    "min_dist_to_barycenter+mol": dict(component_metrics=["min_dist_to_barycenter"],
                                       molecule_metrics=MOL_METRICS),
}

# Feature-set name = "PFG_" + config key (consistent with compare_fingerprints_toxcast.py)
FSET_NAMES = {k: f"PFG_{k}" for k in NEW_CONFIGS}


def main() -> None:
    results_csv = DATA_DIR / "toxcast_comparison_results.csv"

    print("Loading dataset …")
    tox = pd.read_parquet(DATASET_PKL)
    smiles_all = tox["smiles"].tolist()

    print("\nLoading TxP_PFAS fingerprints (for CF mask) …")
    X_txpp_all = load_tsv_fingerprints(TXPP_TSV)

    # ── Compute new embeddings ──────────────────────────────────────
    print("\nComputing new PFASGroups embeddings …")
    pfg = build_pfg_matrices(smiles_all, configs=NEW_CONFIGS)
    for name, arr in pfg.items():
        print(f"  {name}: shape={arr.shape}")

    # ── Exp A: CF-containing subset ──────────────────────────────────
    cf_mask = (X_txpp_all.sum(axis=1) > 0)
    print(f"\nCF-containing chemicals: {int(cf_mask.sum())} / {len(cf_mask)}")
    tox_cf = tox[cf_mask].reset_index(drop=True)

    pfg_cf = {name: arr[cf_mask].astype(np.float32) for name, arr in pfg.items()}

    cv_small = RepeatedStratifiedKFold(
        n_splits=OUTER_SPLITS_A, n_repeats=OUTER_REPEATS, random_state=RANDOM_STATE
    )

    # ── Run nested CV for each new feature set ───────────────────────
    new_rows: list[dict] = []
    for cfg_key, fset_name in FSET_NAMES.items():
        X_cf = pfg_cf[cfg_key]
        print(f"\n─── Experiment A — {fset_name} (shape {X_cf.shape}) ───")
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
                nested_cv_metrics(X_cf[idx], y, cv_small, ep, fset_name, "Exp A")
            )

    # ── Merge with existing results ──────────────────────────────────
    existing = pd.read_csv(results_csv)
    print(f"\nExisting rows: {len(existing)}")

    # Drop any stale rows for the new feature sets
    stale_mask = (existing["experiment"] == "Exp A") & (
        existing["feature_set"].isin(list(FSET_NAMES.values()))
    )
    existing = existing[~stale_mask].copy()
    print(f"Dropped stale rows: {int(stale_mask.sum())}")

    new_df  = pd.DataFrame(new_rows)
    results = pd.concat([existing, new_df], ignore_index=True)
    print(f"New total rows: {len(results)}  (added {len(new_df)})")

    results.to_csv(results_csv, index=False)
    print(f"[saved] {results_csv.name}")

    print_and_save_tables(results, DATA_DIR)
    print("\nDone.")


if __name__ == "__main__":
    main()
