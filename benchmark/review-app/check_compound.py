import json

# Load benchmark data
with open('../data/pfas_oecd_benchmark_20260204_085552.json', 'r') as f:
    data = json.load(f)

# Find the trichlorosilane compound
target_smiles = 'FC(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)CC[Si](Cl)(Cl)Cl'
mol = None
for m in data:
    if m['molecule_data']['smiles'] == target_smiles:
        mol = m
        break

if mol:
    print('Found molecule:')
    print(f'SMILES: {mol["molecule_data"]["smiles"]}')
    print(f'\nDetected groups: {mol["HalogenGroups_result"]["detected_groups"]}')
    
    print('\nGroup matches:')
    for match in mol['HalogenGroups_result']['matches']:
        if match['type'] == 'group':
            print(f'  ID {match["id"]}: {match["name"]}')
    
    # Check against current mapping
    print('\n' + '='*60)
    print('Checking against current pfas_groups_map.json:')
    with open('pfas_groups_map.json', 'r') as f:
        groups_map = json.load(f)
    
    groups_dict = {g['id']: g for g in groups_map}
    
    for gid in mol["HalogenGroups_result"]["detected_groups"]:
        if gid in groups_dict:
            print(f'  ID {gid}: {groups_dict[gid]["name"]} (alias: {groups_dict[gid]["alias"]})')
        else:
            print(f'  ID {gid}: NOT FOUND IN MAP')
else:
    print('Molecule not found')
