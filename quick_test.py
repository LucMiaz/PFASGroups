from PFASgroups import parse_smiles

mol1 = "FC(F)(F)C(F)(F)C(F)(F)CCI"
mol2 = "FC(F)(F)C(F)(F)CCI"

print("Mol1:", mol1)
r1 = parse_smiles(mol1)
ids1 = [g['pfasgroup'].id for g in r1]
print(f"  Groups: {ids1}")
print(f"  Group 71: {71 in ids1}")

print("\nMol2:", mol2)
r2 = parse_smiles(mol2)
ids2 = [g['pfasgroup'].id for g in r2]
print(f"  Groups: {ids2}")
print(f"  Group 71: {71 in ids2}")
