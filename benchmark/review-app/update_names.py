import json
import os

# Read the PFAS groups data
json_path = os.path.join('..', '..', 'HalogenGroups', 'data', 'PFAS_groups_smarts.json')
with open(json_path, 'r', encoding='utf-8') as f:
    data = json.load(f)

# Generate the JavaScript mapping
js_content = """// PFAS Group Names mapping (ID to Name)
// Auto-generated from PFAS_groups_smarts.json - DO NOT EDIT MANUALLY
export const PFAS_GROUP_NAMES = {
"""

for item in data:
    name = item['name'].replace('"', '\\"')
    js_content += f'  {item["id"]}: "{name}",\n'

js_content += "};\n"

# Write to the HalogenGroupNames.js file
output_path = os.path.join('client', 'src', 'data', 'HalogenGroupNames.js')
with open(output_path, 'w', encoding='utf-8') as f:
    f.write(js_content)

print(f"Updated {output_path}")

# Also generate the pfas_groups_map.json for the server
server_map = []
for item in data:
    server_map.append({
        'id': item['id'],
        'name': item['name'],
        'alias': item.get('alias', item['name'])
    })

output_path_json = 'pfas_groups_map.json'
with open(output_path_json, 'w', encoding='utf-8') as f:
    json.dump(server_map, f, indent=2, ensure_ascii=False)

print(f"Updated {output_path_json}")
print(f"Total groups: {len(data)}")
