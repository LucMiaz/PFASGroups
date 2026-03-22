"""Probe group counts and excludeHalogens."""
import sys
sys.path.insert(0, r'C:\Users\luc\git\PFASGroups')
import json
from PFASGroups.core import HALOGEN_GROUPS_FILE

with open(HALOGEN_GROUPS_FILE) as f:
    data = json.load(f)

compute_groups = [g for g in data if g.get('compute', True)]
print(f'Total compute=True groups: {len(compute_groups)}')

excl = [g for g in compute_groups if g.get('excludeHalogens')]
print(f'Groups with excludeHalogens set: {len(excl)}')
for g in excl:
    print(f"  id={g['id']:3d}  excludeHalogens={g['excludeHalogens']}  name={g['name']!r}")

# Reproduce the load_HalogenGroups filter with halogens=['F','Cl','Br','I'] (default)
default_halogens = ['F','Cl','Br','I']
included_default = [g for g in compute_groups if
    g.get('excludeHalogens') is None or
    set(g['excludeHalogens']).isdisjoint(set(default_halogens))]
print(f'\nWith default halogens {default_halogens}: {len(included_default)} groups included')

# Reproduce the filter with halogens='F'
f_only = 'F'
included_f = [g for g in compute_groups if
    g.get('excludeHalogens') is None or
    set(g['excludeHalogens']).isdisjoint(set(f_only) if isinstance(f_only, str) else set(f_only))]
print(f'With halogens="F": {len(included_f)} groups included')

# Show the 4 missing ones
missing = [g for g in included_f if g not in included_default]
print(f'\nGroups present with halogens="F" but missing from default: {len(missing)}')
for g in missing:
    print(f"  id={g['id']:3d}  excludeHalogens={g['excludeHalogens']}  name={g['name']!r}")
