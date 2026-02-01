"""Test script to verify the multi-SMARTS refactoring works correctly."""

from PFASgroups import parse_smiles, get_PFASGroups
from rdkit import Chem

def test_simple_pfca():
    """Test a simple PFCA (perfluoroalkyl carboxylic acid)"""
    # PFOA - Perfluorooctanoic acid
    smiles = "C(=O)(O)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)F"
    results = parse_smiles(smiles)
    
    print("Test 1: PFOA (simple PFCA)")
    print(f"  SMILES: {smiles}")
    
    if len(results) > 0 and 'matches' in results[0]:
        matches = results[0]['matches']
        print(f"  Number of matches: {len(matches)}")
        if len(matches) > 0:
            print(f"  First match: {matches[0]['group_name']}")
            print("  ✓ Test passed")
            print()
            return True
        else:
            print("  ✗ Test failed - no matches found")
            print()
            return False
    else:
        print("  ✗ Test failed - unexpected result format")
        print()
        return False

def test_diacid():
    """Test a diacid that requires 2 copies of the carboxylic acid SMARTS"""
    # Perfluorooctanedicarboxylic acid
    smiles = "C(=O)(O)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(=O)(O)"
    results = parse_smiles(smiles)
    
    print("Test 2: Perfluorooctanedicarboxylic acid (diacid)")
    print(f"  SMILES: {smiles}")
    
    if len(results) > 0 and 'matches' in results[0]:
        matches = results[0]['matches']
        print(f"  Number of matches: {len(matches)}")
        
        # Look for diacid match
        found_diacid = False
        for match in matches:
            if match.get('type') == 'PFASgroup':
                group_name = match['group_name']
                print(f"  Match: {group_name}")
                if 'dicarboxylic' in group_name.lower():
                    found_diacid = True
                    print("  ✓ Test passed - diacid detected")
        
        if not found_diacid:
            print("  ✗ Test failed - diacid not detected")
        print()
        return found_diacid
    else:
        print("  ✗ Test failed - no matches found")
        print()
        return False

def test_pfas_group_structure():
    """Test that PFASGroup objects have the correct new structure"""
    groups = get_PFASGroups()
    
    print("Test 3: PFASGroup structure")
    print(f"  Total groups: {len(groups)}")
    
    # Check a few groups for the new structure
    test_passed = True
    for i, group in enumerate(groups[:3]):
        print(f"  Group {i+1}: {group.name}")
        
        # Check for new attributes
        has_smarts_str = hasattr(group, 'smarts_str')
        has_smarts_count = hasattr(group, 'smarts_count')
        has_smarts_extra_atoms = hasattr(group, 'smarts_extra_atoms')
        
        # Check that old attributes don't exist
        has_old_smarts1 = hasattr(group, 'smarts1')
        has_old_smarts2 = hasattr(group, 'smarts2')
        
        if has_smarts_str and has_smarts_count and has_smarts_extra_atoms:
            print(f"    ✓ Has new attributes (smarts_str, smarts_count, smarts_extra_atoms)")
        else:
            print(f"    ✗ Missing new attributes")
            test_passed = False
            
        if has_old_smarts1 or has_old_smarts2:
            print(f"    ✗ Still has old attributes (smarts1/smarts2)")
            test_passed = False
        else:
            print(f"    ✓ Old attributes removed")
            
        if group.smarts_str:
            print(f"    Number of SMARTS: {len(group.smarts_str)}")
            print(f"    SMARTS counts: {group.smarts_count}")
        print()
    
    if test_passed:
        print("  ✓ All structure tests passed")
    else:
        print("  ✗ Some structure tests failed")
    print()
    return test_passed

def main():
    print("=" * 60)
    print("Testing Multi-SMARTS Refactoring")
    print("=" * 60)
    print()
    
    results = []
    results.append(test_pfas_group_structure())
    results.append(test_simple_pfca())
    results.append(test_diacid())
    
    print("=" * 60)
    print("Summary")
    print("=" * 60)
    print(f"Tests passed: {sum(results)}/{len(results)}")
    
    if all(results):
        print("✓ All tests passed!")
    else:
        print("✗ Some tests failed")
    
    return all(results)

if __name__ == "__main__":
    main()
