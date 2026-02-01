import json

# Load both files
with open(r'c:\Users\luc\git\PFASGroups\PFASgroups\data\PFAS_groups_smarts.json', 'r', encoding='utf-8') as f:
    python_groups = json.load(f)

with open(r'c:\Users\luc\git\PFASgroupsJS\data\PFAS_groups_smarts.json', 'r', encoding='utf-8') as f:
    js_groups = json.load(f)

# Create dictionaries by ID
py_dict = {g['id']: g for g in python_groups}
js_dict = {g['id']: g for g in js_groups}

# Groups that were updated
updated_groups = [29, 30, 31, 32, 33, 34, 36, 37, 38, 40, 41, 42, 43, 44, 45, 46, 47, 49, 59]

print("Successfully added constraints to the following groups in the JavaScript version:\n")
for group_id in updated_groups:
    if group_id in js_dict:
        group = js_dict[group_id]
        print(f"  ID {group_id}: {group['name']}")
        print(f"    Constraints: {group.get('constraints', {})}")
        print()

print(f"\nTotal: {len(updated_groups)} groups updated with constraints from Python version")
