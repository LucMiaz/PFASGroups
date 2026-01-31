import sys
sys.path.insert(0, r'C:\Users\luc\git\PFASGroups')

from rdkit import Chem
from PFASgroups import parse_smiles
import json

# Test molecule that should NOT match group 89
smiles = "C[N+](C)(C)CCCNS(=O)(=O)CCC(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)F"
mol = Chem.MolFromSmiles(smiles)

# Check linker_smarts pattern
linker_smarts = Chem.MolFromSmarts("[CH2X4]")
matches = mol.GetSubstructMatches(linker_smarts)
print(f"Linker [CH2X4] matches: {len(matches)} matches")
print(f"Match indices: {matches}")

# Check the path between functional group and fluorinated chain
# Functional group pattern for group 89
fg_smarts = Chem.MolFromSmarts("[CX4$([#6]N([CH3])([CH3])C)]")
fg_matches = mol.GetSubstructMatches(fg_smarts)
print(f"\nFunctional group matches: {fg_matches}")

# Check perfluoroalkyl chain
pf_smarts = Chem.MolFromSmarts("[C;X4](F)(F)!@!=!#[C;X4](F)(F)")
pf_matches = mol.GetSubstructMatches(pf_smarts)
print(f"Perfluoroalkyl chain matches: {len(pf_matches)} matches")

# Parse with PFASgroups
results = parse_smiles(smiles, output_format='list')
print(f"\n=== PFASgroups Results ===")
for result in results:
    for match in result.get('matches', []):
        if match.get('id') == 89:
            print(f"Group {match['id']}: {match['group_name']}")
            print(f"  Match count: {match['match_count']}")
            print(f"  Components: {match.get('num_components', 0)}")

# Let's also check what's between the functional group and the chain
print("\n=== Molecule structure analysis ===")
for atom in mol.GetAtoms():
    symbol = atom.GetSymbol()
    idx = atom.GetIdx()
    print(f"Atom {idx}: {symbol}")
