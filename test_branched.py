#!/usr/bin/env python
"""Test if Python PFASGroups detects groups for branched PFAS molecules"""

from rdkit import Chem
from PFASgroups.parser import parse_mol

# Test branched PFAS molecule
smiles = "OC(F)(F)C(C(F)(F)F)(C(F)(F)F)C(F)(F)F"
mol = Chem.MolFromSmiles(smiles)
result = parse_mol(mol)
groups = sorted([m['id'] for m in result['matches'] if m['type'] == 'PFASgroup'])
print(f"SMILES: {smiles}")
print(f"Groups: {groups}")
print(f"Expected: [14, 29, 49, 50]")
print(f"Match: {groups == [14, 29, 49, 50]}")
