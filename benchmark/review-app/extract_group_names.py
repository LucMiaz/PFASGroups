import json

# Load the PFAS groups data
with open('../../PFASgroups/data/PFAS_groups_smarts.json', 'r') as f:
    data = json.load(f)

# Create the mapping
group_names = {}
for group in data:
    group_id = str(group['id'])
    name = group['name']
    alias = group.get('alias')
    
    if alias:
        full_name = f"{name} ({alias})"
    else:
        full_name = name
    
    group_names[group_id] = full_name

# Print as JavaScript object
print("export const PFAS_GROUP_NAMES = {")
for group_id, name in group_names.items():
    # Escape quotes in name
    escaped_name = name.replace('"', '\\"')
    print(f'  {group_id}: "{escaped_name}",')
print("};")
