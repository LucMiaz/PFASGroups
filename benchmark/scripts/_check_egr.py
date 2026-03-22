"""Diagnostic: check how much information the EGR block actually adds over binary."""
import numpy as np
import pandas as pd
from pathlib import Path

ROOT     = Path(__file__).resolve().parents[1]
DATA_DIR = ROOT / "data"
TEST_DATA = ROOT / "test_data"

binary = np.load(DATA_DIR / "pfg_binary_cache.npy")
egr    = np.load(DATA_DIR / "pfg_EGR_cache.npy")

print(f"binary shape : {binary.shape}")
print(f"EGR shape    : {egr.shape}")

n_groups = binary.shape[1]          # should be 112
egr_block = egr[:, n_groups:]       # second half = EGR values

print(f"\n=== Full dataset (n={egr.shape[0]}) ===")
nonzero_rows = int((egr_block != 0).any(axis=1).sum())
nan_rows     = int(np.isnan(egr_block).any(axis=1).sum())
all_zero_cols = int((egr_block == 0).all(axis=0).sum())
print(f"  Rows with >=1 nonzero EGR  : {nonzero_rows} / {egr.shape[0]}")
print(f"  Rows with >=1 NaN EGR      : {nan_rows} / {egr.shape[0]}")
print(f"  Columns that are all-zero  : {all_zero_cols} / {n_groups}")
valid = egr_block[egr_block != 0]
if valid.size:
    print(f"  Nonzero values: min={valid.min():.3f}  mean={valid.mean():.3f}  max={valid.max():.3f}")
zero_frac = (egr_block == 0).mean()
print(f"  Zero fraction: {zero_frac*100:.1f}%")
print(f"  First 112 cols == binary: {np.array_equal(egr[:, :n_groups], binary)}")

# CF-containing subset (Experiment A)
txpp_path = TEST_DATA / "TxP_PFAS_v1.tsv"
if txpp_path.exists():
    txpp = pd.read_csv(txpp_path, sep="\t", index_col=0)
    cf_mask = txpp.values.astype(np.uint8).sum(axis=1) > 0
    print(f"\n=== CF-containing subset / Exp A (n={cf_mask.sum()}) ===")
    egr_cf = egr_block[cf_mask]
    nonzero_cf = int((egr_cf != 0).any(axis=1).sum())
    print(f"  Rows with >=1 nonzero EGR  : {nonzero_cf} / {cf_mask.sum()}  ({100*nonzero_cf/cf_mask.sum():.1f}%)")
    vals_cf = egr_cf[egr_cf != 0]
    if vals_cf.size:
        print(f"  Nonzero values: min={vals_cf.min():.2f}  median={np.median(vals_cf):.2f}  mean={vals_cf.mean():.2f}  max={vals_cf.max():.2f}")
    # How many EGR cols have >2 distinct values (i.e., carry information beyond binary)?
    multi_val_cols = sum(np.unique(egr_cf[:, i]).size > 2 for i in range(n_groups))
    print(f"  EGR cols with >2 distinct values: {multi_val_cols} / {n_groups}")
    zero_frac_cf = (egr_cf == 0).mean()
    print(f"  Zero fraction: {zero_frac_cf*100:.1f}%")
else:
    print("\n  (TxP_PFAS_v1.tsv not found — skipping CF subset)")


binary = np.load(DATA_DIR / "pfg_binary_cache.npy")
egr    = np.load(DATA_DIR / "pfg_EGR_cache.npy")

print(f"binary shape : {binary.shape}")
print(f"EGR shape    : {egr.shape}")

n_groups = binary.shape[1]          # should be 112
egr_block = egr[:, n_groups:]       # second half = EGR values

print(f"\n--- EGR block (cols {n_groups}:{egr.shape[1]}) ---")
nonzero_rows = int((egr_block != 0).any(axis=1).sum())
nan_rows     = int(np.isnan(egr_block).any(axis=1).sum())
all_zero_cols = int((egr_block == 0).all(axis=0).sum())
print(f"  Rows with ≥1 nonzero EGR  : {nonzero_rows} / {egr.shape[0]}")
print(f"  Rows with ≥1 NaN EGR      : {nan_rows} / {egr.shape[0]}")
print(f"  Columns that are all-zero : {all_zero_cols} / {n_groups}")
valid = egr_block[egr_block != 0]
if valid.size:
    print(f"  Nonzero values — min={valid.min():.3f}  mean={valid.mean():.3f}  max={valid.max():.3f}")

print(f"\n--- First 112 cols: EGR == binary? ---")
print(f"  Identical: {np.array_equal(egr[:, :n_groups], binary)}")
diff = np.abs(egr[:, :n_groups] - binary).max()
print(f"  Max abs diff: {diff}")

print(f"\n--- Fraction of EGR block that is zero ---")
zero_frac = (egr_block == 0).mean()
print(f"  {zero_frac*100:.1f}% zero entries in EGR block")

# Per-column stats
col_means = np.nanmean(egr_block, axis=0)
col_stds  = np.nanstd(egr_block, axis=0)
print(f"\n--- Per-column stats for EGR block ---")
print(f"  Cols with mean > 0.01 : {(col_means > 0.01).sum()}")
print(f"  Cols with std  > 0.01 : {(col_stds  > 0.01).sum()}")
print(f"\n  Top-10 columns by mean EGR value:")
top10 = np.argsort(col_means)[-10:][::-1]
for i in top10:
    print(f"    col {i+n_groups:3d} (group idx {i:3d}): mean={col_means[i]:.4f}  std={col_stds[i]:.4f}")
