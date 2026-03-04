"""
Generate the two missing deep-analysis figures:
- fp_deep_multihalogen.pdf/.png
- fp_deep_pca_compare.pdf/.png

These were not produced in the main run because a Unicode crash (eta2 in summary print)
occurred after generating the first 3 figures.
"""
import sys
sys.path.insert(0, r'c:\Users\luc\git\PFASGroups')
from pathlib import Path
import pandas as pd
import numpy as np

OUTDIR = Path(r'C:\Users\luc\git\PFASGroups\benchmark\data\figures_fp_benchmarks')
ARTICLE_IMGS = Path(r'C:\Users\luc\git\overleaf_PFASgroups_article\imgs')

from benchmark.scripts.fp_benchmark_deep import (
    build_dataset, compute_fp_variant, compute_descriptive,
    plot_multihalogen, plot_pca_compare, FP_DISPLAY_NAMES
)

print("Building dataset...", flush=True)
df = build_dataset(verbose=True)
df_A = df[df['dataset'] == 'A_F'].reset_index(drop=True)
idx_A = np.where(df['dataset'].values == 'A_F')[0]
smiles_all = df['smiles'].tolist()

print("Computing fingerprints...", flush=True)
fps = {}
print("  binary_F_per", flush=True)
fps['binary_F_per'] = compute_fp_variant(smiles_all, 'F', 'per', 'binary', False)
print("  count_F_per", flush=True)
fps['count_F_per'] = compute_fp_variant(smiles_all, 'F', 'per', 'count', False)
print("  maxcomp_F_per", flush=True)
fps['maxcomp_F_per'] = compute_fp_variant(smiles_all, 'F', 'per', 'max_component', False)
print("  binary_F_poly", flush=True)
fps['binary_F_poly'] = compute_fp_variant(smiles_all, 'F', 'poly', 'binary', False)
print("  binary_F_none", flush=True)
fps['binary_F_none'] = compute_fp_variant(smiles_all, 'F', None, 'binary', False)
print("  binary_FCl_per", flush=True)
fps['binary_FCl_per'] = compute_fp_variant(smiles_all, ['F', 'Cl'], 'per', 'binary', False)
print("  descriptive", flush=True)
fps['descriptive'] = compute_descriptive(smiles_all, fps['binary_F_per'],
                                          fps['count_F_per'], fps['maxcomp_F_per'])

print("Generating fp_deep_multihalogen...", flush=True)
plot_multihalogen(df, fps, OUTDIR, verbose=True)

print("Generating fp_deep_pca_compare...", flush=True)
plot_pca_compare(df_A, fps, idx_A, OUTDIR, verbose=True)

# Copy to article imgs
import shutil
for stem in ('fp_deep_multihalogen', 'fp_deep_pca_compare',
             'fp_deep_split_violin', 'fp_deep_heatmap_summary', 'fp_deep_significance'):
    for ext in ('.pdf', '.png'):
        src = OUTDIR / (stem + ext)
        dst = ARTICLE_IMGS / (stem + ext)
        if src.exists():
            shutil.copy2(src, dst)
            print(f"  Copied {stem}{ext}", flush=True)

print("Done!", flush=True)
