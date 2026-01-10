from rdkit import Chem
from PFASgroups.core import parse_groups_in_mol

# Test with functional group at the end (peripheral)
print("Test 1: Functional group at END (peripheral):")
mol1 = Chem.MolFromSmiles('FC(F)(F)C(F)(F)C(F)(F)C(F)(F)C(=O)O')
matches1 = parse_groups_in_mol(mol1)
for pfas_group, count, comp_lengths, components in matches1[:1]:
    if components:
        comp = components[0]
        print(f"  {pfas_group.name}")
        print(f"    Length: {comp['length']}, Ecc: {comp.get('eccentricity', 0):.3f}, Centr: {comp.get('smarts_centrality', 0):.3f}")

# Test with functional group in the middle (central)
print("\nTest 2: Functional group in MIDDLE (central):")
mol2 = Chem.MolFromSmiles('FC(F)(F)C(F)(F)C(=O)C(F)(F)C(F)(F)F')
matches2 = parse_groups_in_mol(mol2)
for pfas_group, count, comp_lengths, components in matches2[:1]:
    if components:
        comp = components[0]
        print(f"  {pfas_group.name}")
        print(f"    Length: {comp['length']}, Ecc: {comp.get('eccentricity', 0):.3f}, Centr: {comp.get('smarts_centrality', 0):.3f}")

# Test with highly branched structure
print("\nTest 3: HIGHLY BRANCHED structure:")
mol3 = Chem.MolFromSmiles('FC(F)(C(F)(F)F)C(F)(C(F)(F)C(F)(F)F)C(F)(C(F)(F)F)C(=O)O')
matches3 = parse_groups_in_mol(mol3)
for pfas_group, count, comp_lengths, components in matches3[:1]:
    if components:
        comp = components[0]
        print(f"  {pfas_group.name}")
        print(f"    Length: {comp['length']}, Ecc: {comp.get('eccentricity', 0):.3f}, Centr: {comp.get('smarts_centrality', 0):.3f}")

# Test with linear structure
print("\nTest 4: LINEAR structure:")
mol4 = Chem.MolFromSmiles('FC(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(=O)O')
matches4 = parse_groups_in_mol(mol4)
for pfas_group, count, comp_lengths, components in matches4[:1]:
    if components:
        comp = components[0]
        print(f"  {pfas_group.name}")
        print(f"    Length: {comp['length']}, Ecc: {comp.get('eccentricity', 0):.3f}, Centr: {comp.get('smarts_centrality', 0):.3f}")

print("\n" + "="*60)
print("METRIC INTERPRETATION:")
print("="*60)
print("Eccentricity:")
print("  1.0 = Perfectly linear (no branch points)")
print("  0.5 = Moderately branched")
print("  0.0 = Highly branched (all nodes are branch points)")
print()
print("SMARTS Centrality:")
print("  1.0 = Functional group at exact center")
print("  0.5 = Functional group halfway between center and edge")
print("  0.0 = Functional group at periphery (end of chain)")
