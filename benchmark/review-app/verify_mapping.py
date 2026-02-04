import json

# Read source data
with open('../../PFASgroups/data/PFAS_groups_smarts.json', 'r', encoding='utf-8') as f:
    source_data = json.load(f)

print(f'Total groups in source: {len(source_data)}')
print('\nFirst 5 groups:')
for g in source_data[:5]:
    print(f'  ID {g["id"]}: {g["name"]} (alias: {g.get("alias", "N/A")})')

print('\nGroups 73-76:')
for g in source_data:
    if g['id'] in [73, 74, 75, 76]:
        print(f'  ID {g["id"]}: {g["name"]} (alias: {g.get("alias", "N/A")})')

print('\nGroups 98-102:')
for g in source_data:
    if g['id'] in [98, 99, 100, 101, 102]:
        print(f'  ID {g["id"]}: {g["name"]} (alias: {g.get("alias", "N/A")})')

print('\nLast 5 groups:')
for g in source_data[-5:]:
    print(f'  ID {g["id"]}: {g["name"]} (alias: {g.get("alias", "N/A")})')

# Check current mapping
print('\n' + '='*60)
print('Current pfas_groups_map.json:')
with open('pfas_groups_map.json', 'r', encoding='utf-8') as f:
    current_map = json.load(f)

print(f'Total groups in current map: {len(current_map)}')

print('\nGroups 73-76 in current map:')
for g in current_map:
    if g['id'] in [73, 74, 75, 76]:
        print(f'  ID {g["id"]}: {g["name"]} (alias: {g.get("alias", "N/A")})')

print('\nGroups 98-102 in current map:')
for g in current_map:
    if g['id'] in [98, 99, 100, 101, 102]:
        print(f'  ID {g["id"]}: {g["name"]} (alias: {g.get("alias", "N/A")})')

# Check for mismatches
print('\n' + '='*60)
print('Checking for mismatches:')
source_dict = {g['id']: g for g in source_data}
current_dict = {g['id']: g for g in current_map}

mismatches = []
for gid in source_dict:
    if gid in current_dict:
        if source_dict[gid]['name'] != current_dict[gid]['name']:
            mismatches.append(gid)
            print(f'ID {gid}: Source="{source_dict[gid]["name"]}" vs Current="{current_dict[gid]["name"]}"')

if not mismatches:
    print('No mismatches found!')
else:
    print(f'\nTotal mismatches: {len(mismatches)}')
