# OECD PFAS Dataset - ResultsFingerprint Showcase
# v2.2.4

# Import libraries
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from rdkit import Chem
import os

# Import HalogenGroups
from HalogenGroups import parse_smiles

# Load OECD dataset
print("Loading OECD PFAS dataset...")
df = pd.read_csv('tests/test_data/OECDPFAS_list_22012019.csv', encoding='latin-1')
smiles_list = df[df['SMILES'].notna() & (df['SMILES'] != '-')]['SMILES'].tolist()[:200]
print(f"Loaded {len(smiles_list)} SMILES")

# Filter valid SMILES
print("Filtering valid SMILES...")
valid_smiles = []
for smi in smiles_list:
    try:
        mol = Chem.MolFromSmiles(smi)
        if mol is not None:
            valid_smiles.append(smi)
    except:
        pass
print(f"Valid SMILES: {len(valid_smiles)}")

# Parse with HalogenGroups
print("Parsing with HalogenGroups...")
results = parse_smiles(valid_smiles)
print(f"Parsed {len(results)} molecules")

# Convert to fingerprints
print("\nConverting to ResultsFingerprint...")
fp = results.to_fingerprint(group_selection='all', count_mode='binary')
print(fp.summary())

# PCA Analysis
print("\nPerforming PCA...")
pca_result = fp.perform_pca(n_components=10, plot=True, output_file='oecd_pca.png')
cum_var = np.cumsum(pca_result['explained_variance'])
print(f"Variance explained by PC1-2: {cum_var[1]:.1%}")
print(f"Variance explained by PC1-5: {cum_var[4]:.1%}")
print(f"Variance explained by PC1-10: {cum_var[9]:.1%}")

# t-SNE Analysis
print("\nPerforming t-SNE...")
tsne_result = fp.perform_tsne(perplexity=30, plot=True, output_file='oecd_tsne.png')
print(f"t-SNE transformed shape: {tsne_result['transformed'].shape}")

# UMAP Analysis
print("\nPerforming UMAP...")
try:
    umap_result = fp.perform_umap(n_neighbors=15, plot=True, output_file='oecd_umap.png')
    print(f"UMAP transformed shape: {umap_result['transformed'].shape}")
except ImportError:
    print("UMAP not available - install with: pip install umap-learn")

# KL Divergence comparison
print("\nKL Divergence comparison...")
mid = len(valid_smiles) // 2
results_1 = parse_smiles(valid_smiles[:mid])
results_2 = parse_smiles(valid_smiles[mid:])
fp_1 = results_1.to_fingerprint(group_selection='all')
fp_2 = results_2.to_fingerprint(group_selection='all')
kl_div = fp_1.compare_kld(fp_2, method='minmax')
similarity = 1 - kl_div
print(f"KL divergence: {kl_div:.4f}")
print(f"Similarity: {similarity:.1%}")

# SQL Persistence
print("\nSQL Persistence...")
results.to_sql(filename='oecd_results.db', if_exists='replace')
fp.to_sql(filename='oecd_fingerprints.db', if_exists='replace')
print("✓ Results saved to oecd_results.db")
print("✓ Fingerprints saved to oecd_fingerprints.db")

# Load back
from PFASGroups.results_model import ResultsModel, ResultsFingerprint
loaded_results = ResultsModel.from_sql(filename='oecd_results.db')
loaded_fp = ResultsFingerprint.from_sql(filename='oecd_fingerprints.db')
print(f"✓ Loaded {len(loaded_results)} results and {loaded_fp.fingerprint_matrix.shape[0]} fingerprints")

print("\n" + "="*60)
print("SHOWCASE COMPLETE!")
print("="*60)
print("\nGenerated files:")
for f in ['oecd_pca.png', 'oecd_tsne.png', 'oecd_umap.png', 'oecd_results.db', 'oecd_fingerprints.db']:
    if os.path.exists(f):
        print(f"  ✓ {f}")
