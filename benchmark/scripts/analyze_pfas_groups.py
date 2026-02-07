"""
Analyze PFAS groups structure and categorization.
"""

import json

# Load PFAS groups
with open('PFASgroups/data/PFAS_groups_smarts.json', 'r') as f:
    pfas_groups = json.load(f)

print("="*80)
print("PFAS GROUPS ANALYSIS")
print("="*80)

# Count by category
total_groups = len(pfas_groups)
compute_groups = [g for g in pfas_groups if g.get('compute', True)]
aggregate_groups = [g for g in pfas_groups if not g.get('compute', True)]

print(f"\nTotal groups defined: {total_groups}")
print(f"  Computed groups: {len(compute_groups)}")
print(f"  Aggregate groups: {len(aggregate_groups)}")

# Categorize by main_group
main_groups = {}
for g in compute_groups:
    mg = g.get('main_group', 'uncategorized')
    if mg not in main_groups:
        main_groups[mg] = []
    main_groups[mg].append(g)

print(f"\nGroups by main category:")
for mg in sorted(main_groups.keys()):
    print(f"  {mg}: {len(main_groups[mg])} groups")

# Count telomer groups (those with linker_smarts)
telomer_groups = [g for g in compute_groups if g.get('linker_smarts') is not None]
non_telomer_groups = [g for g in compute_groups if g.get('linker_smarts') is None]

print(f"\nGroups with linker validation (telomers): {len(telomer_groups)}")
print(f"Groups without linker validation: {len(non_telomer_groups)}")

# Count OECD vs generic
# OECD groups typically have strict constraints and specific pathways
oecd_like = []
generic_like = []
for g in compute_groups:
    # OECD groups typically have specific aliases and strict constraints
    alias = g.get('alias', '')
    constraints = g.get('constraints', {})
    componentSmarts = g.get('componentSmarts')
    
    # Heuristic: OECD groups have acronyms in alias and strict pathways
    if componentSmarts and 'fluoroalkyl' in componentSmarts.lower() and constraints:
        if any(k in alias for k in ['PF', 'Poly', 'Per']) or 'acid' in g.get('name', '').lower():
            oecd_like.append(g)
        else:
            generic_like.append(g)
    else:
        generic_like.append(g)

print(f"\nOECD-like groups (strict definitions): {len(oecd_like)}")
print(f"Generic groups (broader definitions): {len(generic_like)}")

# Analyze telomer groups specifically
print(f"\nTelomer groups breakdown:")
telomer_by_type = {}
for g in telomer_groups:
    name = g.get('name', '')
    if 'Fluorotelomer' in name:
        # Extract the functional group part
        func = name.replace('Fluorotelomer ', '')
        telomer_by_type[func] = telomer_by_type.get(func, 0) + 1

for func in sorted(telomer_by_type.keys()):
    print(f"  Fluorotelomer {func}: {telomer_by_type[func]}")

# Show aggregate groups
print(f"\nAggregate groups (not computed, pattern matchers):")
for g in aggregate_groups:
    print(f"  ID {g['id']}: {g['name']} - matches pattern: '{g.get('re_search', '')}'")

# Group ID ranges
compute_ids = [g['id'] for g in compute_groups]
telomer_ids = [g['id'] for g in telomer_groups]
print(f"\nGroup ID ranges:")
print(f"  All computed groups: {min(compute_ids)} - {max(compute_ids)}")
print(f"  Telomer groups: {min(telomer_ids)} - {max(telomer_ids)}")

# Save categorization
categorization = {
    'total_groups': total_groups,
    'computed_groups': len(compute_groups),
    'aggregate_groups': len(aggregate_groups),
    'telomer_groups': len(telomer_groups),
    'non_telomer_groups': len(non_telomer_groups),
    'oecd_like_groups': len(oecd_like),
    'generic_groups': len(generic_like),
    'main_categories': {mg: len(groups) for mg, groups in main_groups.items()},
    'telomer_by_type': telomer_by_type,
    'aggregate_group_list': [{'id': g['id'], 'name': g['name'], 'pattern': g.get('re_search', '')} for g in aggregate_groups]
}

with open('benchmark/reports/pfas_groups_categorization.json', 'w') as f:
    json.dump(categorization, f, indent=2)

print(f"\nCategorization saved to: benchmark/reports/pfas_groups_categorization.json")
print("="*80)
