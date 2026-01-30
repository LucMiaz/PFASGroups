import json

# Load both files
with open(r'c:\Users\luc\git\PFASGroups\PFASgroups\data\PFAS_groups_smarts.json', 'r', encoding='utf-8') as f:
    python_groups = json.load(f)

with open(r'c:\Users\luc\git\PFASgroupsJS\data\PFAS_groups_smarts.json', 'r', encoding='utf-8') as f:
    js_groups = json.load(f)

# Create dictionaries by ID
python_dict = {g['id']: g for g in python_groups}
js_dict = {g['id']: g for g in js_groups}

# Find groups with constraints in Python but not in JS (or empty)
missing_constraints = []
for group_id in python_dict:
    if group_id in js_dict:
        py_group = python_dict[group_id]
        js_group = js_dict[group_id]
        
        py_constraints = py_group.get('constraints', {})
        js_constraints = js_group.get('constraints', {})
        
        # Check if Python has non-empty constraints but JS doesn't or is empty
        if py_constraints and (not js_constraints or js_constraints == {}):
            missing_constraints.append({
                'id': group_id,
                'name': py_group['name'],
                'constraints': py_constraints
            })

print(f'Found {len(missing_constraints)} groups with missing constraints in JS version:\n')
for item in missing_constraints:
    print(f"ID {item['id']}: {item['name']}")
    print(f"  Constraints: {json.dumps(item['constraints'], indent=4)}")
    print()
