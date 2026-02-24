#!/usr/bin/env python3
"""Test detection of problematic groups"""

from rdkit import Chem
from HalogenGroups.parser import parse_mol
import json

# Load group definitions
with open('/home/luc/git/HalogenGroups/HalogenGroups/data/PFAS_groups_smarts.json', 'r') as f:
    groups_data = json.load(f)
    groups_dict = {g['id']: g for g in groups_data}

# Test problematic groups
test_cases = [
    {
        'group_id': 52,
        'name': 'alkene',
        'smiles': 'C(F)(F)C(F)=C(F)C(F)(F)C(F)(F)F',
        'description': 'Perfluorinated alkene'
    },
    {
        'group_id': 61,
        'name': 'Silane',
        'smiles': 'C(F)(F)C(F)(F)C([SiH3])(F)F',
        'description': 'Perfluorinated with SiH3'
    },
    {
        'group_id': 63,
        'name': 'Acrylate',
        'smiles': 'C(F)(F)C(F)(F)COC(=O)C=C',
        'description': 'Perfluorinated acrylate'
    },
    {
        'group_id': 64,
        'name': 'Methacrylate',
        'smiles': 'C(F)(F)C(F)(F)COC(=O)C(=C)C',
        'description': 'Perfluorinated methacrylate'
    },
    {
        'group_id': 67,
        'name': 'Thia keto propanoic acid',
        'smiles': 'C(F)(F)C(F)(F)C(=O)SC(=O)C(O)C(=O)O',
        'description': 'Perfluorinated thia keto propanoic acid'
    },
    {
        'group_id': 98,
        'name': 'Glycine',
        'smiles': 'C(F)(F)C(F)(F)CNCC(=O)O',
        'description': 'Perfluorinated glycine'
    },
    {
        'group_id': 100,
        'name': 'Sulfonamidoethanol',
        'smiles': 'C(F)(F)C(F)(F)CS(=O)(=O)NCCO',
        'description': 'Perfluorinated sulfonamidoethanol'
    },
    {
        'group_id': 104,
        'name': 'Sulfonyl propanoic acid',
        'smiles': 'C(F)(F)C(F)(F)CS(=O)(=O)CCC(=O)O',
        'description': 'Perfluorinated sulfonyl propanoic acid'
    },
    {
        'group_id': 106,
        'name': 'Sulfinyl amido sulfonic acid',
        'smiles': 'C(F)(F)C(F)(F)CS(=O)CCC(=O)NCC(C)(C)CS(=O)(=O)O',
        'description': 'Perfluorinated sulfinyl amido sulfonic acid'
    }
]

print("Testing problematic group detection")
print("=" * 80)

for test in test_cases:
    group_id = test['group_id']
    group_info = groups_dict[group_id]
    
    print(f"\nGroup {group_id}: {test['name']}")
    print(f"Description: {test['description']}")
    print(f"Test SMILES: {test['smiles']}")
    print(f"Group SMARTS: {group_info['smarts']}")
    print(f"max_dist_from_comp: {group_info.get('max_dist_from_comp', 'N/A')}")
    print(f"componentSmarts: {group_info.get('componentSmarts')}")
    print(f"Constraints: {group_info.get('constraints', {})}")
    
    # Try to parse
    mol = Chem.MolFromSmiles(test['smiles'])
    if mol is None:
        print("ERROR: Invalid SMILES")
        continue
    
    # Test with parser
    result = parse_mol(mol, include_PFAS_definitions=True)
    detected_groups = [m['group_id'] for m in result['matches'] if 'group_id' in m]
    
    print(f"Detected groups: {detected_groups}")
    print(f"Target group detected: {'YES ✓' if group_id in detected_groups else 'NO ✗'}")
    
    # Test SMARTS pattern directly
    for smarts_pattern, count in group_info['smarts'].items():
        try:
            pattern = Chem.MolFromSmarts(smarts_pattern)
            if pattern:
                matches = mol.GetSubstructMatches(pattern)
                print(f"  SMARTS '{smarts_pattern[:50]}...' matches: {len(matches)}")
            else:
                print(f"  SMARTS '{smarts_pattern[:50]}...' INVALID PATTERN")
        except Exception as e:
            print(f"  SMARTS '{smarts_pattern[:50]}...' ERROR: {e}")
    
    print("-" * 80)

print("\nAnalysis complete!")
