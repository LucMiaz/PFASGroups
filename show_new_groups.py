import json

with open(r'PFASgroups\data\PFAS_groups_smarts.json', 'r') as f:
    data = json.load(f)

print(f"Total groups: {len(data)}")
print(f"\nNew groups added (IDs 75-87):\n")

for g in data:
    if g['id'] >= 75:
        print(f"ID {g['id']:2d}: {g['name']:<50s} ({g['alias']})")
        print(f"       smartsPath: {g['smartsPath']}, linker_smarts: {g['linker_smarts']}, max_dist: {g['max_dist_from_CF']}")
        print(f"       SMARTS: {list(g['smarts'].keys())[0]}")
        print()
