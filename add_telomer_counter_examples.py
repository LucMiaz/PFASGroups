#!/usr/bin/env python3
"""
Add counter-examples to telomer groups by copying examples.
The CH2 chains will need to be manually removed later.
"""
import json
from pathlib import Path

# Load the JSON file
file_path = Path('PFASgroups/data/PFAS_groups_smarts.json')
with open(file_path, 'r', encoding='utf-8') as f:
    groups = json.load(f)

# Process telomer groups
modified_count = 0
for group in groups:
    # Check if it's a telomer group (has 'telomer' in name or test category)
    is_telomer = False
    if 'telomer' in group.get('name', '').lower():
        is_telomer = True
    elif group.get('test', {}).get('category') == 'telomer':
        is_telomer = True
    elif group.get('test', {}).get('generate', {}).get('is_telomer'):
        is_telomer = True
    
    if is_telomer and 'test' in group and 'examples' in group['test']:
        examples = group['test']['examples']
        if examples:  # Only if there are examples
            # Check if counter-examples already exist
            if 'counter-examples' not in group['test'] or not group['test']['counter-examples']:
                # Copy examples to counter-examples
                group['test']['counter-examples'] = examples.copy()
                modified_count += 1
                print(f"Group {group['id']}: {group['name']} - Added {len(examples)} counter-examples (copied from examples)")

print(f"\nModified {modified_count} telomer groups")
print("\nNote: Counter-examples are copies of examples.")
print("You'll need to manually edit them to remove CH2 chains to make them non-telomers.")

# Save the modified JSON
with open(file_path, 'w', encoding='utf-8') as f:
    json.dump(groups, f, indent=4)
    
print(f"\nSaved modified file to {file_path}")
