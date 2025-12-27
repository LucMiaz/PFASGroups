"""
Simplified test to check JavaScript implementation directly using file import.
"""
from rdkit import Chem
from PFASgroups.core import parse_groups_in_mol
import json

# Test molecules
TEST_MOLECULES = [
    {
        'name': 'PFOA',
        'smiles': 'C(=O)(O)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)F',
    },
    {
        'name': 'PFOS',
        'smiles': 'C(C(C(C(C(C(C(C(F)(F)S(=O)(=O)O)(F)F)(F)F)(F)F)(F)F)(F)F)(F)F)(F)F',
    },
    {
        'name': '6:2 FTOH',
        'smiles': 'C(CCO)(C(C(C(C(C(F)(F)F)(F)F)(F)F)(F)F)(F)F)(F)F',
    },
]

def test_python_both_flavors():
    """Test that Python implementation works with both flavors."""
    print("="*80)
    print("Testing Python Implementation - Both Flavors")
    print("="*80)
    
    all_tests_pass = True
    
    for test_case in TEST_MOLECULES:
        print(f"\n{'='*80}")
        print(f"Molecule: {test_case['name']}")
        print(f"SMILES: {test_case['smiles']}")
        print(f"{'='*80}")
        
        mol = Chem.MolFromSmiles(test_case['smiles'])
        if mol is None:
            print(f"❌ Failed to parse SMILES")
            all_tests_pass = False
            continue
        
        formula = Chem.rdMolDescriptors.CalcMolFormula(mol)
        
        # Test bycomponent=False
        print(f"\n--- Testing bycomponent=False ---")
        result_default = parse_groups_in_mol(mol, formula, bycomponent=False)
        groups_default = result_default.get('detected_groups', [])
        print(f"Detected groups: {[g['id'] for g in groups_default]}")
        for group in groups_default:
            print(f"  Group {group['id']}: {group['name']}")
            print(f"    match_count: {group.get('match_count', 0)}, chain_lengths: {group.get('chain_lengths', [])}")
        
        # Test bycomponent=True
        print(f"\n--- Testing bycomponent=True ---")
        result_bycomp = parse_groups_in_mol(mol, formula, bycomponent=True)
        groups_bycomp = result_bycomp.get('detected_groups', [])
        print(f"Detected groups: {[g['id'] for g in groups_bycomp]}")
        for group in groups_bycomp:
            print(f"  Group {group['id']}: {group['name']}")
            print(f"    match_count: {group.get('match_count', 0)}, chain_lengths: {group.get('chain_lengths', [])}")
        
        # Check adequation
        ids_default = set(g['id'] for g in groups_default)
        ids_bycomp = set(g['id'] for g in groups_bycomp)
        
        if ids_default == ids_bycomp:
            print(f"\n✅ Perfect adequation - both flavors detect the same groups")
        else:
            print(f"\n⚠️  Adequation mismatch:")
            print(f"   Only in default: {ids_default - ids_bycomp}")
            print(f"   Only in bycomponent: {ids_bycomp - ids_default}")
            all_tests_pass = False
    
    print(f"\n\n{'='*80}")
    if all_tests_pass:
        print("✅ ALL PYTHON TESTS PASSED")
    else:
        print("❌ SOME PYTHON TESTS FAILED")
    print(f"{'='*80}\n")
    
    return all_tests_pass

if __name__ == '__main__':
    success = test_python_both_flavors()
    exit(0 if success else 1)
