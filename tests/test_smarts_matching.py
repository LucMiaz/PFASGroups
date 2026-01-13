"""Understand what SMARTS patterns match"""
from rdkit import Chem

# Test different SMARTS patterns
smarts_patterns = {
    "Carboxylic acid (current)": "[#6$([#6][#6](=O)([OH1,Oh1,O-]))]",
    "Carboxylic acid (full)": "C(=O)O",
    "Sulfonic acid (current)": "[#6$([#6][#16](=O)(=O)[OH1,Oh1,O-])]",
    "Sulfonic acid (full)": "S(=O)(=O)O",
}

test_mol = Chem.MolFromSmiles("CC(=O)O")  # Acetic acid
test_mol = Chem.AddHs(test_mol)

print("Test molecule: Acetic acid (CC(=O)O)")
print(f"Total atoms: {test_mol.GetNumAtoms()}")
print()

for name, smarts_str in smarts_patterns.items():
    smarts = Chem.MolFromSmarts(smarts_str)
    matches = test_mol.GetSubstructMatches(smarts)
    
    print(f"{name}:")
    print(f"  SMARTS: {smarts_str}")
    print(f"  Atoms in SMARTS pattern: {smarts.GetNumAtoms()}")
    print(f"  Number of matches: {len(matches)}")
    
    if matches:
        for i, match in enumerate(matches):
            print(f"  Match {i+1}: atom indices {match}")
            for atom_idx in match:
                atom = test_mol.GetAtomWithIdx(atom_idx)
                print(f"    Atom {atom_idx}: {atom.GetSymbol()}")
    print()
