"""Test script for linker_smarts functionality - Testing Fluorotelomer Alcohols (Group 15).

Fluorotelomer alcohols have the structure: CF3(CF2)n-(CH2)m-OH
where the CH2 chain connects the perfluorinated part to the alcohol group.

Group 15 has linker_smarts set to "[#6H2X4]" (only CH2 groups allowed as linkers)
and max_dist_from_CF=12 to allow for various chain lengths.
"""

from rdkit import Chem
from PFASgroups import parse_smiles

# Test cases for fluorotelomer alcohols with different CH2 chain lengths
test_cases = [
    # (SMILES, Description, CH2_count, Expected_Group_15)
    ("FC(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)CCO", "6:2 Fluorotelomer alcohol (CF3-CF2-CF2-CH2-CH2-OH)", 2, True),
    ("C(F)(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)CCCO", "6:3 Fluorotelomer alcohol (CF3-CF2-CF2-CH2-CH2-CH2-OH)", 3, True),
    ("C(F)(F)(F)C(F)(F)C(F)(F)C(F)(F)CO", "4:1 Fluorotelomer alcohol (CF3-CF2-CH2-OH)", 1, True),
    ("C(F)(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)CCCCO", "6:4 Fluorotelomer alcohol (CF3-CF2-CF2-CH2-CH2-CH2-CH2-OH)", 4, True),
    ("FC(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)CCCO", "10:3 Fluorotelomer alcohol", 3, True),
    
    # Negative test: oxygen linker instead of CH2 (should NOT match group 15)
    ("C(F)(F)(F)CO", "CF3-CF2-CH2-OH (valid)", 1, True),
    ("C(F)(F)(F)OCCO", "CF3-CF2-O-CH2-CH2-OH (oxygen linker - should fail)", -1, False),
    ("C(F)(F)(F)C(F)(F)CC(F)CO", "CF3-CF2-CH2-CFH-CH2-OH (unsaturated fluorinated linker)", -1, False),
    
    # Direct connection (no CH2 linker between perfluorinated and alcohol)
    ("C(F)(F)CO", "CF2-CH2-OH (too short perfluoro chain)", 1, False),
]

print("Testing Fluorotelomer Alcohols (Group 15) with linker_smarts validation")
print("Group 15 linker_smarts: '[#6H2X4]' (only CH2 groups allowed)")
print("=" * 80)

for smiles, description, ch2_count, should_match in test_cases:
    print(f"\n{description}")
    print(f"  SMILES: {smiles}")
    
    mol = Chem.MolFromSmiles(smiles)
    if not mol:
        print("  ERROR: Invalid SMILES")
        continue
    
    # Parse the molecule
    result = parse_smiles(smiles)
    
    if len(result) == 0 or 'matches' not in result[0] or len(result[0]['matches']) == 0:
        print("  Result: No PFAS groups found")
        group_15_found = False
    else:
        # Check if group 15 (fluorotelomer alcohols) was found
        matches = result[0]['matches']
        pfas_groups = [m for m in matches if m['type'] == 'PFASgroup']
        group_15_found = any(m['id'] == 15 for m in pfas_groups)
        
        print(f"  Result: {len(pfas_groups)} PFAS group(s) found")
        for match in pfas_groups:
            if match['id'] == 15:
                print(f"    ✓ Group 15 (fluorotelomer alcohols) - MATCHED")
                print(f"      Match count: {match['match_count']}")
                print(f"      Component sizes: {match['components_sizes']}")
    
    # Validate expected result
    status = "✓ PASS" if group_15_found == should_match else "✗ FAIL"
    expected = "MATCH" if should_match else "NO MATCH"
    actual = "MATCH" if group_15_found else "NO MATCH"
    print(f"  {status}: Expected {expected}, Got {actual}")

print("\n" + "=" * 80)
print("\nSummary:")
print("  • Fluorotelomer alcohols require CH2 linkers only (linker_smarts: '[#6H2X4]')")
print("  • Molecules with O or other atom types in the linker should NOT match Group 15")
print("  • Different CH2 chain lengths (1-4+) should all match if structure is valid")

