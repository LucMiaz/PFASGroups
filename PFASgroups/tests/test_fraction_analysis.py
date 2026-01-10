"""Compare component fractions with and without SMARTS precomputation"""
import sys
import os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))

from PFASgroups.core import parse_smiles
from rdkit import Chem

# Test with PFOA (Perfluorooctanoic acid)
pfoa_smiles = "C(C(C(C(C(C(C(C(=O)O)(F)F)(F)F)(F)F)(F)F)(F)F)(F)F)(F)(F)F"

mol = Chem.MolFromSmiles(pfoa_smiles)
mol = Chem.AddHs(mol)
total_atoms = mol.GetNumAtoms()

print(f"PFOA Analysis:")
print(f"SMILES: {pfoa_smiles}")
print(f"Total atoms (with H): {total_atoms}")
print()

# Count atoms manually
c_count = sum(1 for atom in mol.GetAtoms() if atom.GetSymbol() == 'C')
f_count = sum(1 for atom in mol.GetAtoms() if atom.GetSymbol() == 'F')
o_count = sum(1 for atom in mol.GetAtoms() if atom.GetSymbol() == 'O')
h_count = sum(1 for atom in mol.GetAtoms() if atom.GetSymbol() == 'H')

print(f"Atom counts: C={c_count}, F={f_count}, O={o_count}, H={h_count}")
print(f"Total: {c_count + f_count + o_count + h_count}")
print()

# Parse with PFASgroups
results = parse_smiles([pfoa_smiles], bycomponent=True)

if results:
    result = results[0]
    if 'matches' in result:
        for match in result['matches']:
            if 'Perfluoroalkyl carboxylic acid' in match.get('group_name', ''):
                print(f"Group: {match.get('group_name', 'Unknown')}")
                print(f"Total components fraction: {match.get('total_components_fraction', 0):.3f}")
                print()
                
                if 'components' in match:
                    for i, comp in enumerate(match['components']):
                        carbon_backbone_size = comp['size']
                        component_fraction = comp['component_fraction']
                        implied_atoms = int(component_fraction * total_atoms)
                        
                        print(f"Component {i+1}:")
                        print(f"  Carbon backbone size: {carbon_backbone_size}")
                        print(f"  Component fraction: {component_fraction:.3f}")
                        print(f"  Implied atom count: {implied_atoms} / {total_atoms}")
                        
                        # Calculate expected coverage
                        # 7 carbons in backbone (C8 - 1 for COOH)
                        # Each C has 2 F attached = 14 F
                        # Plus 1 COOH carbon with 1 O attached
                        # Plus 1 OH
                        expected = carbon_backbone_size + (carbon_backbone_size * 2) + 1  # 7 C + 14 F + 1 C(COOH)
                        print(f"  Expected coverage (manual): ~{expected} atoms")
                        print(f"    {carbon_backbone_size} carbons + {carbon_backbone_size * 2} fluorines + 1 COOH carbon")
                        print()
                break
