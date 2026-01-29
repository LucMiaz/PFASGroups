import sys
sys.path.insert(0, '.')

from PFASgroups import parse_smiles, get_PFASGroups

print("Testing new fluorotelomer groups...\n")

# Get all groups
groups = get_PFASGroups()
print(f"✓ Loaded {len(groups)} PFAS groups\n")

# Test molecules
test_cases = [
    # Saturated fluorotelomer alcohol (should match ID 74 and generic ID 86)
    ("FC(F)(F)C(F)(F)CCOH", "Saturated fluorotelomer alcohol"),
    
    # Unsaturated fluorotelomer carboxylic acid (should match ID 67 and generic ID 87)
    ("FC(F)(F)C(F)(F)C(F)=CC(=O)O", "Unsaturated fluorotelomer carboxylic acid"),
    
    # Unsaturated fluorotelomer alcohol (should match ID 85 and generic ID 87)
    ("FC(F)(F)C(F)(F)C(F)=CCOH", "Unsaturated fluorotelomer alcohol"),
]

for smiles, description in test_cases:
    print(f"Testing: {description}")
    print(f"SMILES: {smiles}")
    
    try:
        results = parse_smiles(smiles)
        
        if results and len(results) > 0:
            mol_result = results[0]
            
            if isinstance(mol_result, dict) and 'matches' in mol_result:
                matches = mol_result['matches']
                telomer_matches = [m for m in matches if 'telomer' in m.get('group_name', '').lower()]
                
                if telomer_matches:
                    print(f"  ✓ Found {len(telomer_matches)} telomer group(s):")
                    for match in telomer_matches:
                        print(f"    - ID {match['id']}: {match['group_name']}")
                else:
                    print(f"  ⚠ No telomer groups detected")
            else:
                print(f"  Results: {mol_result}")
        else:
            print(f"  ⚠ No results")
            
    except Exception as e:
        print(f"  ✗ Error: {e}")
    
    print()

print("✅ Testing complete!")
