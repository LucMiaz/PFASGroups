"""Add smarts1_size and smarts2_size fields to PFAS_groups_smarts.json"""
import json
import os

# Load the JSON file
json_path = os.path.join(os.path.dirname(__file__), '..', 'data', 'PFAS_groups_smarts.json')
with open(json_path, 'r') as f:
    groups = json.load(f)

# Add the new fields to each group
for group in groups:
    # Add smarts1_size if not already present
    if 'smarts1_size' not in group:
        group['smarts1_size'] = None
    
    # Add smarts2_size if not already present
    if 'smarts2_size' not in group:
        group['smarts2_size'] = None

# Save the updated JSON with proper formatting
with open(json_path, 'w') as f:
    json.dump(groups, f, indent=4)

print(f"Updated {len(groups)} groups with smarts1_size and smarts2_size fields")
