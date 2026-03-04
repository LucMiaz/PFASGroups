"""Quick test of fp_benchmark_deep.py building on a small subset."""
import sys, traceback
sys.path.insert(0, r'c:\Users\luc\git\PFASGroups')

try:
    print("Step 1: imports")
    import numpy as np
    import pandas as pd
    import matplotlib
    matplotlib.use('Agg')
    print("  numpy, pandas, matplotlib OK")

    print("Step 2: scipy")
    from scipy.stats import kruskal, mannwhitneyu
    print("  scipy OK")

    print("Step 3: HalogenGroups")
    from HalogenGroups import generate_fingerprint
    fp, info = generate_fingerprint('FC(F)(F)C(F)(F)C(=O)O',
                                    representation='vector',
                                    count_mode='binary',
                                    halogens='F', saturation='per')
    print(f"  FP shape: {fp.shape}, active: {int(fp.sum())}")

    print("Step 4: build_dataset (small)")
    from benchmark.scripts.fp_benchmark_deep import build_dataset, compute_fp_variant
    df = build_dataset(verbose=True)
    print(f"  Dataset: {len(df)} compounds")

    print("Step 5: compute binary_F_per on first 10 Dataset A compounds")
    smi_10 = df[df['dataset'] == 'A_F']['smiles'].head(10).tolist()
    X = compute_fp_variant(smi_10, 'F', 'per', 'binary', verbose=False)
    print(f"  FP matrix shape: {X.shape}")

    print("\nAll tests PASSED")

except Exception:
    traceback.print_exc()
