#!/usr/bin/env python3
"""Test new SMILES patterns for groups 67, 104, 106"""

import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from rdkit import Chem
from rdkit.Chem import AllChem

# Simulate molecule generation with new SMILES
test_groups = {
    67: {'smiles': 'C(=O)SC(=O)C(O)C(=O)O', 'mode': 'attach', 'name': 'Thia keto propanoic acid'},
    104: {'smiles': 'S(=O)(=O)CCC(=O)O', 'mode': 'attach', 'name': 'Sulfonyl propanoic acid'},
    106: {'smiles': 'S(=O)CCC(=O)NCC(C)(C)CS(=O)(=O)O', 'mode': 'attach', 'name': 'Sulfinyl amido sulfonic acid'},
}

print("🧪 TESTING NEW SMILES PATTERNS\n")
print("=" * 80)

for group_id, group_info in test_groups.items():
    print(f"\n📊 Group {group_id}: {group_info['name']}")
    print(f"   Functional group SMILES: {group_info['smiles']}")
    
    # Simulate molecule generation by attaching to perfluoro chain
    # For max_dist_from_CF:0, the functional group attaches directly to a CF carbon
    
    # Create a simple perfluorinated molecule with the functional group attached
    # Mode 'attach' means the functional group is attached to the end of the perfluoro chain
    
    # Example: C(F)(F)C(F)(F)C(F)(F)F with functional group attached
    # The functional group connects at the terminal carbon
    
    # Construct test SMILES by concatenating perfluoro chain + functional group
    perfluoro_chain = "FC(F)(F)C(F)(F)C(F)(F)"  # Terminal perfluoro unit
    full_smiles = perfluoro_chain + group_info['smiles']
    
    print(f"   Generated test molecule: {full_smiles}")
    
    mol = Chem.MolFromSmiles(full_smiles)
    if mol:
        canonical = Chem.MolToSmiles(mol)
        print(f"   Canonical SMILES: {canonical}")
        
        # Check if molecule is valid
        try:
            AllChem.Compute2DCoords(mol)
            print(f"   ✅ Valid molecule structure")
            
            # Count atoms
            num_atoms = mol.GetNumAtoms()
            num_f = sum(1 for atom in mol.GetAtoms() if atom.GetSymbol() == 'F')
            num_c = sum(1 for atom in mol.GetAtoms() if atom.GetSymbol() == 'C')
            
            print(f"   Atoms: {num_atoms} total, {num_c} C, {num_f} F")
            
        except Exception as e:
            print(f"   ⚠️ Structure generation issue: {e}")
    else:
        print(f"   ❌ Invalid SMILES generated")

print("\n" + "=" * 80)
print("💡 ANALYSIS:\n")
print("The new patterns attach functional groups directly to perfluoro carbons.")
print("For groups with max_dist_from_CF: 0, this should satisfy the constraint.")
print("\nNext step: Re-run benchmarks to validate these fixes.")
