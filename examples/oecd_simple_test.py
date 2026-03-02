"""Simple test of ResultsFingerprint on OECD dataset"""
import pandas as pd
import sys
from rdkit import Chem

print("Loading dataset...")
df = pd.read_csv('tests/test_data/OECDPFAS_list_22012019.csv', encoding='latin-1')
smiles_list = df[df['SMILES'].notna() & (df['SMILES'] != '-')]['SMILES'].tolist()[:500]

print(f"Loaded {len(smiles_list)} SMILES")

# Filter valid
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

# Parse
print("Parsing with HalogenGroups...")
from HalogenGroups import parse_smiles
results = parse_smiles(valid_smiles[:100])  # Start with just 100
print(f"Parsed {len(results)} molecules")

# Convert to fingerprints
print("Converting to fingerprints...")
fp = results.to_fingerprint(group_selection='all', count_mode='binary')
print(fp.summary())

# PCA
print("\nPerforming PCA...")
pca = fp.perform_pca(n_components=5, plot=False)
print(f"PCA shape: {pca['transformed'].shape}")
print(f"Explained variance: {pca['explained_variance'][:3]}")

print("\n✓ Test successful!")
