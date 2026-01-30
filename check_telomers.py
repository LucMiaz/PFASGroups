import json

with open('PFASgroups/data/PFAS_groups_smarts.json', 'r') as f:
    data = json.load(f)

print(f'Total groups: {len(data)}')
max_id = max(g['id'] for g in data)
print(f'Max ID: {max_id}')

telomers = [g for g in data if g.get('linker_smarts') is not None]
print(f'Telomer groups: {len(telomers)}')
print(f'Telomer IDs: {sorted([g["id"] for g in telomers])}')

non_telomers_60_65 = [g for g in data if 60 <= g['id'] <= 65]
print(f'\nGroups 60-65 (should NOT be telomers):')
for g in non_telomers_60_65:
    has_linker = g.get('linker_smarts') is not None
    print(f"  ID {g['id']}: {g['name']:40s} linker_smarts={'YES' if has_linker else 'NO '}")
