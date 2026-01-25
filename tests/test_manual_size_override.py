"""
Test demonstrating how to use manual SMARTS size overrides.

This test shows how the smarts1_size and smarts2_size fields in the JSON
can be used to manually specify the number of atoms, overriding the automatic
pattern counting.
"""
import sys
import os
import json
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))

from rdkit import Chem
from PFASgroups.PFASGroupModel import PFASGroup
from PFASgroups.parser import parse_groups_in_mol

def test_with_manual_override():
    """Test component fraction with manual SMARTS size specification."""
    
    # Load groups from JSON
    data_dir = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'PFASgroups','data')
    with open(os.path.join(data_dir, 'PFAS_groups_smarts.json'), 'r') as f:
        groups_data = json.load(f)
    
    pfas_groups = [PFASGroup(**group_data) for group_data in groups_data]
    
    # Test molecule: PFOA (C8HF15O2)
    smiles = "C(=O)(O)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)F"
    mol = Chem.MolFromSmiles(smiles)
    
    print("PFOA Test - Current Behavior (automatic counting):")
    print(f"SMILES: {smiles}")
    print(f"Total atoms: {mol.GetNumAtoms()}")
    print()
    
    # Parse with current automatic counting
    results = parse_groups_in_mol(mol, pfas_groups=pfas_groups)
    
    # Check if results is a dict (group_id -> components) or list
    if isinstance(results, dict):
        for group_id, components in results.items():
            group = pfas_groups[group_id - 1]
            print(f"Group {group_id}: {group.group_name}")
            print(f"  smarts1_extra_atoms (automatic): {group.smarts1_extra_atoms}")
            if group.smarts2:
                print(f"  smarts2_extra_atoms (automatic): {group.smarts2_extra_atoms}")
            
            for comp in components:
                if 'component_fraction' in comp:
                    print(f"  Component fraction: {comp['component_fraction']:.3f}")
                    break
            print()
    else:
        # Results is a list of all components
        for comp in results:
            if 'group_id' in comp and 'component_fraction' in comp:
                group_id = comp['group_id']
                group = pfas_groups[group_id - 1]
                print(f"Group {group_id}: {group.group_name}")
                print(f"  smarts1_extra_atoms (automatic): {group.smarts1_extra_atoms}")
                if group.smarts2:
                    print(f"  smarts2_extra_atoms (automatic): {group.smarts2_extra_atoms}")
                print(f"  Component fraction: {comp['component_fraction']:.3f}")
                print()
                break  # Only show first matching group
    
    print("\n" + "="*70)
    print("How to use manual override:")
    print("="*70)
    print("""
To manually specify the SMARTS size for a group, edit PFAS_groups_smarts.json:

Example for Group 1 (Perfluoroalkyl carboxylic acids):
{
    "id": 1,
    "name": "Perfluoroalkyl carboxylic acids",
    "smarts1": "[#6$([#6][#6](=O)([OH1,Oh1,O-]))]",
    "smarts2": null,
    "constraints": {...},
    "smarts1_size": 3,  // <-- Set this to total atoms (C + 2O = 3)
    "smarts2_size": null
}

The code will then use: extra_atoms = smarts1_size - 1 = 3 - 1 = 2

Notes:
- Leave as null to use automatic pattern counting
- Specify total atoms (including the matched atom)
- Code automatically subtracts 1 for the matched atom
- Useful for edge cases where automatic counting is inaccurate
""")

if __name__ == "__main__":
    test_with_manual_override()
