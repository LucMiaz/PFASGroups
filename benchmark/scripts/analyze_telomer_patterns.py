#!/usr/bin/env python3
"""Analyze fluorotelomer end groups to generate better SMARTS patterns."""

import json
from collections import Counter

# Load validation results
with open('data/telomer_validation_results.json') as f:
    data = json.load(f)

# Analyze detected telomer groups (excluding generic 86/87)
end_groups = Counter()
for result in data['results']:
    if result.get('detected'):
        for group in result.get('telomer_groups', []):
            gid = group['id']
            if gid != 86 and gid != 87:
                end_groups[f"{gid}: {group['name']}"] += 1

print('Specific telomer end groups detected in PubChem dataset:')
print('='*70)
for group, count in end_groups.most_common():
    print(f'  {group}: {count} detections')
print()

# Load PFAS groups to analyze patterns
pfas_json = r'c:\Users\luc\git\HalogenGroups\HalogenGroups\data\PFAS_groups_smarts.json'
with open(pfas_json) as f:
    groups = json.load(f)

print('Existing saturated fluorotelomer SMARTS patterns (62-74):')
print('='*70)
for g in groups:
    if 62 <= g['id'] <= 74:
        print(f"ID {g['id']:2d} ({g['name'][:30]:30s}): {g['smarts']}")
print()

print('Existing unsaturated fluorotelomer SMARTS patterns (75-85):')
print('='*70)
for g in groups:
    if 75 <= g['id'] <= 85:
        print(f"ID {g['id']:2d} ({g['name'][:30]:30s}): {g['smarts']}")
print()

# Extract just the functional group part (after CH2$)
print('Functional group patterns (extracted):')
print('='*70)
patterns = []
for g in groups:
    if 62 <= g['id'] <= 85 and g['id'] not in [86, 87]:
        smarts = g['smarts']
        # Extract the functional group part
        if 'CH2$' in smarts or 'CH2\\$' in smarts:
            # For saturated telomers
            start = smarts.find('(') if '(' in smarts else 0
            end = smarts.find('))') if '))' in smarts else len(smarts)
            if start > 0:
                fg = smarts[start+1:end]
                if fg and fg not in ['CH2', '[CH2]']:
                    patterns.append(f"{g['id']:2d}: {fg}")
                    
print('Saturated telomer end groups to combine:')
for p in patterns:
    print(f'  {p}')
print()

# Generate combined SMARTS
print('Suggested SMARTS patterns:')
print('='*70)
print()
print('For Group 86 (generic saturated telomers):')
print('Current: [CH2]')
print('Improved: Combine common end groups:')
print('  [CH2$([CH2][OH1,Oh1,O-,I,SiH3,Si,#16,#15,#8])]')
print('  Or more specific:')
saturated_fgs = [
    '[OH1,Oh1,O-]',  # Alcohols, ethoxylates
    '[SiH3]',  # Silanes
    '[Si](Cl)(Cl)Cl',  # Trichlorosilanes
    'I',  # Iodides
    'C(=O)[H]',  # Aldehydes
    'C(=O)[OH1,Oh1,O-]',  # Carboxylic acids
    '[#16](=O)(=O)[OH1,Oh1,O-]',  # Sulfonic acids
    'OP(=O)(=O)[OH1,Oh1,O-]',  # Monophosphates
    'OP(=O)(=O)O[CH2][CH2]',  # Diphosphates
    'OC(O)[CH1]=[CH2]',  # Acrylates
    'OC(O)C(C)=[CH2]',  # Methacrylates
]
combined = '|'.join(saturated_fgs)
print(f'  [CH2$([CH2]({combined}))]')
print()
print('For Group 87 (generic unsaturated telomers):')
print('Current: [CH2$(C(F)=C)]')
print('Could remain as is (captures the unsaturation) or add end groups')
