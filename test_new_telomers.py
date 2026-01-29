import json
import sys

# Test JSON validity
try:
    with open(r'PFASgroups\data\PFAS_groups_smarts.json', 'r') as f:
        data = json.load(f)
    
    print(f"✓ Valid JSON with {len(data)} groups")
    
    # Find new groups
    new_groups = [g for g in data if g['id'] >= 75]
    print(f"\n✓ Added {len(new_groups)} new groups:")
    
    for g in new_groups:
        print(f"  ID {g['id']}: {g['name']} ({g['alias']})")
    
    # Test that we can import and use them
    sys.path.insert(0, '.')
    from PFASgroups import get_PFASGroups
    
    groups = get_PFASGroups()
    print(f"\n✓ Successfully loaded {len(groups)} PFAS groups into PFASgroups")
    
    # Find telomer groups
    telomer_groups = [g for g in groups if 'telomer' in g.name.lower()]
    unsaturated_telomers = [g for g in telomer_groups if 'unsaturated' in g.name.lower()]
    
    print(f"\n✓ Total fluorotelomer groups: {len(telomer_groups)}")
    print(f"✓ Unsaturated fluorotelomer groups: {len(unsaturated_telomers)}")
    
    # Show generic groups
    generic_groups = [g for g in groups if g.id in [86, 87]]
    print(f"\n✓ Generic telomer groups:")
    for g in generic_groups:
        print(f"  ID {g.id}: {g.name}")
    
    print("\n✅ All tests passed!")
    
except json.JSONDecodeError as e:
    print(f"✗ JSON Error: {e}")
    sys.exit(1)
except Exception as e:
    print(f"✗ Error: {e}")
    import traceback
    traceback.print_exc()
    sys.exit(1)
