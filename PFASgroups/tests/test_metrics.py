from rdkit import Chem
from PFASgroups.core import parse_groups_in_mol

# Test with PFOA (linear chain)
print("Testing PFOA (linear chain):")
mol1 = Chem.MolFromSmiles('C(C(C(C(C(C(C(C(F)(F)F)(F)F)(F)F)(F)F)(F)F)(F)F)(F)F)(=O)O')
matches1 = parse_groups_in_mol(mol1)

for pfas_group, count, comp_lengths, components in matches1:
    if components and len(components) > 0:
        comp = components[0]
        print(f"  {pfas_group.name}:")
        print(f"    Length: {comp['length']}")
        print(f"    Eccentricity: {comp.get('eccentricity', 'N/A'):.3f}")
        print(f"    SMARTS Centrality: {comp.get('smarts_centrality', 'N/A'):.3f}")
        if len(components) > 1:
            print(f"    (showing first of {len(components)} components)")
        break

# Test with branched PFAS
print("\nTesting Branched PFAS:")
mol2 = Chem.MolFromSmiles('FC(F)(C(F)(F)F)C(F)(C(F)(F)F)C(=O)O')
matches2 = parse_groups_in_mol(mol2)

for pfas_group, count, comp_lengths, components in matches2:
    if components and len(components) > 0:
        comp = components[0]
        print(f"  {pfas_group.name}:")
        print(f"    Length: {comp['length']}")
        print(f"    Eccentricity: {comp.get('eccentricity', 'N/A'):.3f}")
        print(f"    SMARTS Centrality: {comp.get('smarts_centrality', 'N/A'):.3f}")
        if len(components) > 1:
            print(f"    (showing first of {len(components)} components)")
        break

# Test with PFOS (sulfonic acid at end)
print("\nTesting PFOS (terminal functional group):")
mol3 = Chem.MolFromSmiles('FC(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)S(=O)(=O)O')
matches3 = parse_groups_in_mol(mol3)

for pfas_group, count, comp_lengths, components in matches3:
    if components and len(components) > 0:
        comp = components[0]
        print(f"  {pfas_group.name}:")
        print(f"    Length: {comp['length']}")
        print(f"    Eccentricity: {comp.get('eccentricity', 'N/A'):.3f}")
        print(f"    SMARTS Centrality: {comp.get('smarts_centrality', 'N/A'):.3f}")
        if len(components) > 1:
            print(f"    (showing first of {len(components)} components)")
        break

print("\nMetrics explanation:")
print("- Eccentricity: 1.0 = linear chain, 0.0 = highly branched")
print("- SMARTS Centrality: 1.0 = functional group at center, 0.0 = at periphery")
