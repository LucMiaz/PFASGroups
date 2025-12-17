"""
Script to draw misclassified molecules from the specificity testing.
Uses RDKit to generate molecular structure images for inclusion in the LaTeX document.
"""

import pandas as pd
from rdkit import Chem
from rdkit.Chem import Draw
import os
import ast

# Get the directory of this script
script_dir = os.path.dirname(os.path.abspath(__file__))

# Load misclassified molecules
df = pd.read_csv(os.path.join(script_dir, 'misclassified_molecules.csv'))

# Parse group lists from string representation
def parse_group_list(group_str):
    """Parse group list from string or list representation"""
    if pd.isna(group_str):
        return []
    if isinstance(group_str, list):
        return group_str
    try:
        # Try to evaluate as Python literal
        return ast.literal_eval(str(group_str))
    except:
        # Fallback: remove brackets and parse
        cleaned = str(group_str).strip('[]')
        if not cleaned:
            return []
        return [int(x.strip()) for x in cleaned.split(',')]

# Apply parsing to group columns
df['expected_groups'] = df['expected_groups'].apply(parse_group_list)
df['detected_groups'] = df['detected_groups'].apply(parse_group_list)

print(f"Total misclassified molecules: {len(df)}")
print(f"\nBreakdown by failure type:")
print(df['failure_type'].value_counts())
print(f"\nBreakdown by origin (top 20):")
print(df['origin'].value_counts().head(20))

# Function to draw molecules by category
def draw_molecules_by_category(df, origin_filter, output_file, n_samples=6, mols_per_row=3):
    """
    Draw a grid of molecules matching the origin filter.
    
    Args:
        df: DataFrame with misclassified molecules
        origin_filter: String to filter origin column (e.g., 'azole-per')
        output_file: Output PNG file path
        n_samples: Number of molecules to draw
        mols_per_row: Number of molecules per row in the grid
    """
    # Filter molecules
    filtered = df[df['origin'] == origin_filter].head(n_samples)
    
    if len(filtered) == 0:
        print(f"No molecules found for filter: {origin_filter}")
        return
    
    # Convert SMILES to RDKit molecules
    mols = []
    legends = []
    for idx, row in filtered.iterrows():
        mol = Chem.MolFromSmiles(row['smiles'])
        if mol is not None:
            mols.append(mol)
            expected_str = str(row['expected_groups']) if isinstance(row['expected_groups'], list) else row['expected_groups']
            detected_str = str(row['detected_groups']) if isinstance(row['detected_groups'], list) else row['detected_groups']
            legends.append(f"{row['origin']}\nExpected: {expected_str}\nDetected: {detected_str}")
    
    # Draw molecules
    if mols:
        img = Draw.MolsToGridImage(
            mols,
            molsPerRow=mols_per_row,
            subImgSize=(400, 400),
            legends=legends,
            returnPNG=False
        )
        img.save(output_file)
        print(f"Saved {len(mols)} molecules to {output_file}")
    else:
        print(f"No valid molecules to draw for {origin_filter}")

# Generate figures for each category
categories = [
    ('azole-per', 'azole_perfluorinated.png'),
    ('azole-poly', 'azole_polyfluorinated.png'),
    ('azine-per', 'azine_perfluorinated.png'),
    ('azine-poly', 'azine_polyfluorinated.png')
]

print("\nGenerating molecular structure images...")
for origin, filename in categories:
    output_path = os.path.join(script_dir, filename)
    draw_molecules_by_category(df, origin, output_path, n_samples=9, mols_per_row=3)

# Create a combined figure with representatives from each category
print("\nGenerating combined representative figure...")
representatives = []
rep_mols = []
rep_legends = []

for origin in ['azole-per', 'azole-poly', 'azine-per', 'azine-poly']:
    sample = df[df['origin'] == origin].head(2)
    for idx, row in sample.iterrows():
        mol = Chem.MolFromSmiles(row['smiles'])
        if mol is not None:
            rep_mols.append(mol)
            rep_legends.append(f"{row['origin']}")

if rep_mols:
    img = Draw.MolsToGridImage(
        rep_mols,
        molsPerRow=4,
        subImgSize=(350, 350),
        legends=rep_legends,
        returnPNG=False
    )
    img.save(os.path.join(script_dir, 'misclassified_representatives.png'))
    print(f"Saved combined representative figure with {len(rep_mols)} molecules")

print("\nDone! Images saved to:")
print(script_dir)
print("\nTo include in LaTeX, copy the PNG files to your LaTeX figures directory and update:")
print("  - Figure \\ref{fig:misclassified_azole}")
print("  - Figure \\ref{fig:misclassified_azine}")
