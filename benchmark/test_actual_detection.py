#!/usr/bin/env python3
"""Test actual HalogenGroups detection on problem molecules"""

import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from rdkit import Chem
from HalogenGroups.parser import parse_mol

test_molecules = {
    67: 'O=C(O)C(O)C(=O)SC(=O)C(F)(C(F)(F)F)C(F)(F)F',
    104: 'O=C(O)CCS(=O)(=O)CC(C(F)(F)F)(C(F)(F)F)C(F)(F)C(F)(F)C(F)(F)F',
    106: 'CC(C)(CNC(=O)CCS(=O)CC(F)(F)C(F)(F)C(F)(F)F)CS(=O)(=O)O',
}

print("🧪 TESTING HalogenGroupS DETECTION\n")
print("=" * 80)

for group_id, smiles in test_molecules.items():
    print(f"\n📊 Group {group_id}:")
    print(f"   SMILES: {smiles}")
    
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
            print(f"   ✅ Group {group_id} WAS detected!")
            # Find the match details
            for match in result['matches']:
                if match.get('type') == 'HalogenGroup' and match['id'] == group_id:
                    print(f"      Match details: {match}")
        else:
            print(f"   ❌ Group {group_id} was NOT detected")
            print(f"   🔍 Checking why...")
            
            # Look for any error messages or filtering
            if 'errors' in result:
                print(f"      Errors: {result['errors']}")
            
            # Check what groups were actually found
            print(f"      Groups found instead: {detected_groups}")
    else:
        print(f"   ❌ Unexpected result format: {result}")

print("\n" + "=" * 80)
