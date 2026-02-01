#!/usr/bin/env python
"""Test quaternary carbon pattern matching in Python RDKit"""

from rdkit import Chem
import json

# Test molecule: highly branched fluoroalcohol
smiles = "OC(F)(F)C(C(F)(F)F)(C(F)(F)F)C(F)(F)F"
print(f"Test molecule: {smiles}")
print(f"Expected: Python detects groups [14, 29, 49, 50]")
print(f"Currently JS detects: [49, 50]\n")

mol = Chem.MolFromSmiles(smiles)
mol = Chem.AddHs(mol)

# Show atom indices with symbols
print("Molecule structure (with explicit H):")
for i, atom in enumerate(mol.GetAtoms()):
    symbol = atom.GetSymbol()
    degree = atom.GetDegree()
    totalH = atom.GetTotalNumHs()
    print(f"  Atom {i}: {symbol}, degree={degree}, totalH={totalH}")

# Load fpaths
with open('c:/Users/luc/git/PFASGroups/PFASgroups/data/fpaths.json', 'r') as f:
    fpaths = json.load(f)

print('\n' + '='*80)

# Test patterns
patterns = {
    "Quaternary carbon (simple)": "[#6H0X4]",
    "Quaternary with 4 C neighbors": "[#6$([#6H0X4]([#6])([#6])([#6])[#6])]",
    "Quaternary with 4 quaternary neighbors": "[#6$([#6H0X4]([#6H0X4])([#6H0X4])([#6H0X4])[#6H0X4])]",
    "Quaternary with 4 quaternary neighbors (acyclic)": "[C$([#6H0X4]([#6H0X4])([#6H0X4])[#6H0X4])!R]",
    "Perfluoroalkyl chain (full, from fpaths)": fpaths['Perfluoroalkyl']['chain']
}

for name, pattern in patterns.items():
    print(f"\n{name}:")
    print(f"  Pattern: {pattern}")
    
    try:
        patt = Chem.MolFromSmarts(pattern)
        if patt:
            matches = mol.GetSubstructMatches(patt)
            if matches:
                # Collect all matched atoms
                all_atoms = set()
                for match in matches:
                    for atom_idx in match:
                        all_atoms.add(atom_idx)
                
                sorted_atoms = sorted(all_atoms)
                print(f"  ✓ Found {len(matches)} match(es)")
                print(f"    Matched atoms: {sorted_atoms}")
                
                # Show what these atoms are
                atom_details = [f"{idx}:{mol.GetAtomWithIdx(idx).GetSymbol()}" for idx in sorted_atoms]
                print(f"    Atom types: {', '.join(atom_details)}")
            else:
                print(f"  ✗ No matches")
        else:
            print(f"  ✗ Failed to compile pattern")
    except Exception as e:
        print(f"  ✗ Error: {e}")

print('\n' + '='*80)
print('\nANALYSIS:')
print('If atom 4 is matched by the quaternary carbon pattern, what atom is it?')
print('Is it the central quaternary carbon that connects CF2-OH to the CF3 branches?')
