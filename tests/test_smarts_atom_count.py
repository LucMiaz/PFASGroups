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
    
    if group.smarts1 is not None:
        num_atoms1 = group.smarts1.GetNumAtoms()
        print(f"  smarts1: {group_data['smarts1']}")
        print(f"    Total atoms: {num_atoms1}")
        print(f"    Extra atoms (total - 1): {group.smarts1_extra_atoms}")
    else:
        print(f"  smarts1: None")
        print(f"    Extra atoms: {group.smarts1_extra_atoms}")
    
    if group.smarts2 is not None:
        num_atoms2 = group.smarts2.GetNumAtoms()
        print(f"  smarts2: {group_data['smarts2']}")
        print(f"    Total atoms: {num_atoms2}")
        print(f"    Extra atoms (total - 1): {group.smarts2_extra_atoms}")
    elif group_data.get('smarts2'):
        print(f"  smarts2: {group_data['smarts2']} (not compiled)")
    else:
        print(f"  smarts2: None")
        print(f"    Extra atoms: {group.smarts2_extra_atoms}")

print("\n" + "="*70)
