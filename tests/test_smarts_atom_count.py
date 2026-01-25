"""Test SMARTS atom count precomputation"""
import sys
import os
import json
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))

from PFASgroups.PFASGroupModel import PFASGroup
from rdkit import Chem

# Load PFAS groups from JSON
data_dir = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'PFASgroups','data')
with open(os.path.join(data_dir, 'PFAS_groups_smarts.json'), 'r') as f:
    groups_data = json.load(f)

print("Testing SMARTS atom count precomputation:\n")
print("="*70)

# Test first 10 groups
for group_data in groups_data[:10]:
    group = PFASGroup(**group_data)
    print(f"\nGroup {group.id}: {group.name}")
    if not group.smarts or len(group.smarts) == 0:
        print("  No SMARTS patterns defined.")
        continue
    for i in range(len(group.smarts)):
        smarts_mol = group.smarts[i]
        num_atoms = smarts_mol.GetNumAtoms()
        min_count = group.smarts_count[i]
        smarts_pattern = group.smarts_str[i]
        print(f"  SMARTS {i+1}: {smarts_pattern}")
        print(f"    Total atoms: {num_atoms}")
        print(f"    Minimum required matches: {min_count}")
    
print("\n" + "="*70)
