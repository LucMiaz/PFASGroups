"""Test generation of group 11 molecules"""
from rdkit import Chem
from rdkit.Chem.rdMolDescriptors import CalcMolFormula
import sys
import os

parent_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, parent_dir)

from PFASgroups.core import parse_PFAS_groups
from PFASgroups.generate_mol import generate_random_mol

# Group 11 definition from test_examples.py:
# (11, "Perfluoroalkyl sulfinic acids", [{"group_smiles":"S(=O)O[H]", 'n':1, 'mode':'attach'}], 'Perfluoroalkyl')

template = [{"group_smiles":"S(=O)O[H]", 'n':1, 'mode':'attach'}]
pathtype = 'Perfluoroalkyl'

print("Generating group 11 test molecules...")
print("="*60)

for i in range(5):
    try:
        mol = generate_random_mol(8, template, perfluorinated=(pathtype == 'Perfluoroalkyl'))
        if mol is None:
            continue
            
        smiles = Chem.MolToSmiles(mol)
        formula = CalcMolFormula(mol)
        
        # Check for sulfur
        atoms = [atom.GetSymbol() for atom in mol.GetAtoms()]
        has_sulfur = 'S' in atoms
        
        # Test classification
        matches_default = parse_PFAS_groups(mol, formula, bycomponent=False)
        matches_bycomp = parse_PFAS_groups(mol, formula, bycomponent=True)
        
        detected_default = [m[0].id for m in matches_default]
        detected_bycomp = [m[0].id for m in matches_bycomp]
        
        group11_in_default = 11 in detected_default
        group11_in_bycomp = 11 in detected_bycomp
        
        print(f"\nMolecule {i+1}:")
        print(f"  SMILES: {smiles}")
        print(f"  Formula: {formula}")
        print(f"  Has Sulfur: {has_sulfur}")
        print(f"  Detected (default): {detected_default}")
        print(f"  Detected (bycomponent): {detected_bycomp}")
        print(f"  Group 11 detected (default): {group11_in_default}")
        print(f"  Group 11 detected (bycomponent): {group11_in_bycomp}")
        
        if not group11_in_default or not group11_in_bycomp:
            print(f"  ⚠️  WARNING: Group 11 not detected!")
        else:
            print(f"  ✓ Group 11 correctly detected")
            
    except Exception as e:
        print(f"Error generating molecule {i+1}: {e}")

print("\n" + "="*60)
