#!/usr/bin/env python3
"""Test if HalogenGroups detects the new molecules correctly"""

import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from rdkit import Chem
from HalogenGroups.parser import parse_mol

# New test molecules with direct attachment
test_molecules = {
    67: 'O=C(O)C(O)C(=O)SC(=O)C(F)(F)C(F)(F)C(F)(F)F',
    104: 'O=C(O)CCS(=O)(=O)C(F)(F)C(F)(F)C(F)(F)F',
    106: 'S(=O)(CCC(=O)NC(C)(C)CS(=O)(=O)O)C(F)(F)C(F)(F)C(F)(F)F',
}

print("🧪 TESTING HalogenGroupS DETECTION WITH NEW MOLECULES\n")
print("=" * 80)

for group_id, smiles in test_molecules.items():
    print(f"\n📊 Group {group_id}:")
    print(f"   SMILES: {smiles}")
    print(f"   (Functional group attached DIRECTLY to perfluoro carbon)")
    
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        print("   ❌ Invalid SMILES")
        continue
    
    # Test with HalogenGroups
    result = parse_mol(mol, include_PFAS_definitions=True)
    
    if isinstance(result, dict) and 'matches' in result:
        detected_groups = [m['id'] for m in result['matches'] if m.get('type') == 'HalogenGroup']
        print(f"   Detected groups: {detected_groups}")
        
        if group_id in detected_groups:
            print(f"   ✅ Group {group_id} WAS DETECTED!")
            # Find the match details
            for match in result['matches']:
                if match.get('type') == 'HalogenGroup' and match['id'] == group_id:
                    print(f"      Match count: {match.get('match_count', 1)}")
        else:
            print(f"   ❌ Group {group_id} was NOT detected")
            
            # Check if fluorotelomer version was detected instead
            telomer_versions = {67: 97, 104: 103, 106: 105}
            telomer_id = telomer_versions.get(group_id)
            if telomer_id and telomer_id in detected_groups:
                print(f"   ⚠️  Detected as fluorotelomer version (Group {telomer_id}) instead")
            
            # Check what groups were actually found
            if detected_groups:
                print(f"      Groups found instead: {detected_groups}")
    else:
        print(f"   ❌ Unexpected result format")

print("\n" + "=" * 80)
print("💡 EXPECTED RESULT:\n")
print("All three groups (67, 104, 106) should now be detected because:")
print("1. Functional groups are attached directly to perfluoro carbons")
print("2. Satisfies max_dist_from_CF: 0 constraint")
print("3. SMARTS patterns match the generated structures")
