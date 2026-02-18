import json

# Load PFAS groups from the main data file
with open('c:/Users/luc/git/HalogenGroups/HalogenGroups/data/PFAS_groups_smarts.json', 'r') as f:
    data = json.load(f)

# Extract only id, name, and alias for the lookup map
groups = [{'id': g['id'], 'name': g['name'], 'alias': g.get('alias', '')} for g in data]

# Save to review-app directory
with open('c:/Users/luc/git/HalogenGroups/benchmark/review-app/pfas_groups_map.json', 'w') as f:
    json.dump(groups, f, indent=2)

print(f'Created pfas_groups_map.json with {len(groups)} groups')
print('First few groups:')
for g in groups[:5]:
    print(f"  ID {g['id']}: {g['name']}")
