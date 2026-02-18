import json
with open('../../HalogenGroups/data/PFAS_groups_smarts.json','r') as f:
    data = json.load(f)
print(f"Loaded {len(data)} PFAS groups")

OECD = []
generic = []
telomers = []
for group in data:
    if group.get('id') <29:
        OECD.append(group)
    elif 'telomer' in group.get('name','').lower():
        telomers.append(group)
    else:
        generic.append(group)
print(f"OECD groups: {len(OECD)}")
print(f"Generic groups: {len(generic)}")
print(f"Telomer groups: {len(telomers)}")

new_data = []
new_data_JS = []
id =1
old_to_new = {}
for og in [OECD, generic]:
    for group in og:
        old_to_new[group['id']] = id
        group['id'] = id
        new_data.append({k:v for k,v in group.items() if k not in ['main_group','base_functional_groups']})
        new_data_JS.append({k:v for k,v in group.items() if k not in ['test','linker_smarts','main_group','base_functional_groups']})
        id +=1
for group in telomers:
    old_to_new[group['id']] = id
    group['id'] = id
    new_data.append({k:v for k,v in group.items() if k not in ['main_group','base_functional_groups']})
    id +=1

with open('../../HalogenGroups/data/PFAS_groups_smarts_old.json','w') as f:
    json.dump(data,f, indent=4)
with open('../../HalogenGroups/data/PFAS_groups_smarts.json','w') as f:
    json.dump(new_data,f, indent=4)
with open('../../../HalogenGroupsJS/data/PFAS_groups_smarts_reordered.json','w') as f:
    json.dump(new_data_JS,f, indent=4)
print(f"Reordered PFAS groups saved. Total groups: {len(new_data)}")
with open('../data/PFAS_group_id_mapping.json','w') as f:
    json.dump(old_to_new,f, indent=4)