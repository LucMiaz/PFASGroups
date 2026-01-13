"""Debug SMARTS matching to understand what atom is matched"""
from rdkit import Chem

# Test molecule: FC(F)(F)C(F)(F)C(F)(F)C(=O)O
smiles = "FC(F)(F)C(F)(F)C(F)(F)C(=O)O"
mol = Chem.MolFromSmiles(smiles)

# SMARTS for perfluoroalkyl carboxylic acid
smarts_str = "[#6$([#6][#6](=O)([OH1,Oh1,O-]))]"
smarts = Chem.MolFromSmarts(smarts_str)

print(f"Molecule: {smiles}")
print(f"SMARTS: {smarts_str}")
print(f"\nAtoms in molecule:")
for i, atom in enumerate(mol.GetAtoms()):
    print(f"  Atom {i}: {atom.GetSymbol()} (neighbors: {[n.GetIdx() for n in atom.GetNeighbors()]})")

matches = mol.GetSubstructMatches(smarts)
print(f"\nSMARTS matches: {matches}")

if matches:
    for match_idx in matches[0]:
        atom = mol.GetAtomWithIdx(match_idx)
        print(f"\nMatched atom {match_idx}: {atom.GetSymbol()}")
        print(f"  Neighbors: {[(n.GetIdx(), n.GetSymbol()) for n in atom.GetNeighbors()]}")

# Count carbons
total_carbons = sum(1 for atom in mol.GetAtoms() if atom.GetSymbol() == 'C')
print(f"\nTotal carbons in molecule: {total_carbons}")

# Identify the carboxylic acid group atoms
print("\nCarboxylic acid group should be:")
print("  C=O carbon at index:", mol.GetSubstructMatch(Chem.MolFromSmarts("C(=O)O")))
