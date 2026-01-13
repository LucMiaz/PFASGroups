"""Test Polyfluorobr SMARTS pattern"""
from rdkit import Chem

# Test molecule 1: NHS ester with perfluoro chain
smiles1 = "FC(F)(F)C(F)(F)C(F)(F)C(F)(F)CCC(=O)ON1C(=O)CCC1=O"
mol1 = Chem.MolFromSmiles(smiles1)

# Polyfluorobr chain SMARTS
polyfluorobr_chain = "[#6$([#6X4](~[#6,#9])~[#6,#7,#8,#16,F,#1,#17,#35,#53]),#6$([#6X3](~[#9,#6])~[#6,#7,#8,#16,F,#1,#17,#35,#53]),#6$([#6X2]~[#9,#6]),$([#8X2,#7]([$([#6][#9])])[$([#6][#9])]),$([#8X2,#7,#15,#16]([$([#6][#9])])[$([#6][#9])]),$([#6](=O)([$(*[#9])])[$(*[#9])])]"

# Compile and find matches
smarts_mol = Chem.MolFromSmarts(polyfluorobr_chain)
matches = mol1.GetSubstructMatches(smarts_mol)

print("Polyfluorobr chain matches:")
print(f"Found {len(matches)} matches")

all_matched_atoms = set()
for match in matches:
    all_matched_atoms.update(match)

print(f"\nAll matched atoms: {sorted(all_matched_atoms)}")

# Show which atoms these are
print("\nMatched atom details:")
for idx in sorted(all_matched_atoms):
    atom = mol1.GetAtomWithIdx(idx)
    neighbors = [f"{mol1.GetAtomWithIdx(n).GetSymbol()}{n}" for n in [x.GetIdx() for x in atom.GetNeighbors()]]
    print(f"  Atom {idx}: {atom.GetSymbol()} (neighbors: {neighbors})")

# Check if NHS ring carbons (19, 21, 22, 23) are matched
nhs_carbons = [19, 21, 22, 23]
print(f"\nNHS ring carbons {nhs_carbons} in matches: {all(c in all_matched_atoms for c in nhs_carbons)}")

# Test the problematic sub-pattern
print("\n\nTesting sub-patterns:")

# Pattern for nitrogen with two carbons having fluorine
# This SHOULD NOT match the NHS nitrogen
n_pattern = "$([#8X2,#7]([$([#6][#9])])[$([#6][#9])])"
n_smarts = Chem.MolFromSmarts(n_pattern)
n_matches = mol1.GetSubstructMatches(n_smarts)
print(f"\nNitrogen pattern matches: {n_matches}")
if n_matches:
    for match in n_matches:
        for idx in match:
            atom = mol1.GetAtomWithIdx(idx)
            print(f"  Matched atom {idx}: {atom.GetSymbol()}")

# Check if the carbons in NHS ring are directly connected to fluorine
print("\n\nChecking if NHS ring carbons have fluorine neighbors:")
for c_idx in nhs_carbons:
    carbon = mol1.GetAtomWithIdx(c_idx)
    has_fluorine = any(mol1.GetAtomWithIdx(n).GetSymbol() == 'F' for n in [x.GetIdx() for x in carbon.GetNeighbors()])
    print(f"  Carbon {c_idx}: has_fluorine = {has_fluorine}")

# Test a simpler pattern: carbon connected to fluorine
c_with_f = Chem.MolFromSmarts("[#6][#9]")
cf_matches = mol1.GetSubstructMatches(c_with_f)
cf_carbons = set(m[0] for m in cf_matches)
print(f"\n\nCarbons directly connected to fluorine: {sorted(cf_carbons)}")
print(f"NHS carbons in this list: {[c for c in nhs_carbons if c in cf_carbons]}")
