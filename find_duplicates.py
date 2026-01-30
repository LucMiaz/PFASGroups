#!/usr/bin/env python
"""Find duplicate IDs in PFAS_groups_smarts.json."""
import json

with open('PFASgroups/data/PFAS_groups_smarts.json', 'r') as f:
    data = json.load(f)

print('Groups with duplicate IDs:')
print('='*70)
for g in data:
    id_val = g['id']
    if id_val in [71, 74, 75]:
        print(f"ID {id_val:2d}: {g['name']:45s} ({g.get('alias','')})")

max_id = max([g['id'] for g in data])
print(f'\nMax ID currently: {max_id}')
print(f'Next available IDs: {max_id+1} onwards')

# Show what original IDs 71, 74, 75 should be
print('\n' + '='*70)
print('Original telomer groups (should keep their IDs):')
for g in data:
    if g['id'] in [71, 74, 75] and 'unsaturated' not in g['name'].lower():
        print(f"ID {g['id']:2d}: {g['name']:45s} ({g.get('alias','')})")
