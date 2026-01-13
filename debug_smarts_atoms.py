"""Count carbon atoms in SMARTS patterns"""
from rdkit import Chem

# Original SMARTS from group 1
smarts_str = "[#6$([#6][#6](=O)([OH1,Oh1,O-]))]"
smarts = Chem.MolFromSmarts(smarts_str)

print(f"SMARTS: {smarts_str}")
print(f"Number of atoms in SMARTS: {smarts.GetNumAtoms()}")
print("\nAtoms in SMARTS:")
for i, atom in enumerate(smarts.GetAtoms()):
    print(f"  Atom {i}: {atom.GetSmarts() if hasattr(atom, 'GetSmarts') else atom.GetSymbol()}")

# Count carbons in SMARTS
print(f"\nCarbons in SMARTS pattern: counting by atomic num = {sum(1 for atom in smarts.GetAtoms() if atom.GetAtomicNum() == 6)}")

# Let me test with a simpler SMARTS
simple_cooh = "C(=O)O"
simple_smarts = Chem.MolFromSmarts(simple_cooh)
print(f"\nSimple COOH SMARTS: {simple_cooh}")
print(f"Atoms: {simple_smarts.GetNumAtoms()}")
print("Atom types:")
for i, atom in enumerate(simple_smarts.GetAtoms()):
    print(f"  Atom {i}: AtomicNum={atom.GetAtomicNum()}, Symbol={atom.GetSymbol()}")
