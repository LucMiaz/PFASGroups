#!/usr/bin/env python3
"""
Test OECD robustness analysis with a small subset
"""

import sys
import pandas as pd
from collections import defaultdict

# Add paths
sys.path.append('/home/luc/git/PFAS-atlas')
sys.path.append('/home/luc/git/PFASGroups')

try:
    from classification_helper import classify_pfas_molecule
    print("✅ Successfully imported PFAS-Atlas classify_pfas_molecule")
except ImportError as e:
    print(f"❌ Failed to import PFAS-Atlas: {e}")
    sys.exit(1)

try:
    from PFASgroups import parse_pfas
    print("✅ Successfully imported PFASGroups parse_pfas")
except ImportError as e:
    print(f"❌ Failed to import PFASGroups: {e}")
    sys.exit(1)

def test_oecd_analysis():
    """Test OECD analysis with a small subset."""
    
    # Load a small subset of OECD data
    oecd_csv_path = '/home/luc/git/PFAS-atlas/input_data/OECD_4000/step3_OECD_Class_0812.csv'
    
    try:
        df = pd.read_csv(oecd_csv_path)
        print(f"📊 Loaded {len(df)} OECD molecules")
        
        # Test with first 5 molecules only
        test_molecules = []
        for i, row in df.head(5).iterrows():
            test_molecules.append({
                'rdkit_smiles': row['RDKIT_SMILES'].strip(),
                'original_smiles': row['SMILES'].strip() if pd.notna(row['SMILES']) else row['RDKIT_SMILES'].strip(), 
                'first_class': row['First_Class'],
                'second_class': row['Second_Class'],
                'type': row['type']
            })
        
        print(f"🧪 Testing with {len(test_molecules)} molecules")
        
        atlas_agreements = 0
        pfasgroups_detections = 0
        oecd_groups_correspondence = defaultdict(lambda: defaultdict(int))
        
        for i, mol in enumerate(test_molecules):
            smiles = mol['rdkit_smiles']
            expected_first = mol['first_class']
            expected_second = mol['second_class']
            expected_type = mol['type']
            
            print(f"\n🔬 Testing molecule {i+1}: {smiles[:50]}...")
            print(f"   Expected: {expected_first} / {expected_second} (Type: {expected_type})")
            
            # Test PFAS-Atlas
            try:
                atlas_result = classify_pfas_molecule(smiles)
                atlas_first = atlas_result[0] if len(atlas_result) > 0 else 'Not detected'
                atlas_second = atlas_result[1] if len(atlas_result) > 1 else 'Not detected'
                
                print(f"   Atlas:    {atlas_first} / {atlas_second}")
                
                if atlas_first == expected_first and atlas_second == expected_second:
                    atlas_agreements += 1
                    print("   ✅ Atlas match!")
                else:
                    print("   ❌ Atlas mismatch")
                    
            except Exception as e:
                print(f"   ❌ Atlas error: {e}")
            
            # Test PFASGroups
            try:
                print(f"   Testing PFASGroups with: {smiles}")
                pfas_result = parse_pfas(smiles)
                print(f"   PFASGroups raw result: {pfas_result}")
                
                if pfas_result and len(pfas_result) > 0:
                    pfasgroups_detections += 1
                    # Extract group IDs from the result
                    detected_groups = []
                    for group_dict in pfas_result:
                        if isinstance(group_dict, dict) and 'id' in group_dict:
                            detected_groups.append(group_dict['id'])
                    
                    oecd_detected_groups = [g for g in detected_groups if g < 29]
                    
                    print(f"   PFASGroups: ✅ Detected (Groups: {oecd_detected_groups})")
                    
                    if oecd_detected_groups:
                        for group in oecd_detected_groups:
                            oecd_groups_correspondence[expected_type][group] += 1
                    else:
                        oecd_groups_correspondence[expected_type]['no_detection'] += 1
                else:
                    print(f"   PFASGroups: ❌ Not detected")
                    oecd_groups_correspondence[expected_type]['no_detection'] += 1
                    
            except Exception as e:
                print(f"   ❌ PFASGroups error: {e}")
        
        print(f"\n📊 Test Results:")
        print(f"   Atlas accuracy: {atlas_agreements}/{len(test_molecules)} ({100*atlas_agreements/len(test_molecules):.1f}%)")
        print(f"   PFASGroups detections: {pfasgroups_detections}/{len(test_molecules)} ({100*pfasgroups_detections/len(test_molecules):.1f}%)")
        
        print(f"\n🔗 OECD Groups Correspondence:")
        for oecd_type, detections in oecd_groups_correspondence.items():
            print(f"   {oecd_type}: {dict(detections)}")
        
        return True
        
    except Exception as e:
        print(f"❌ Error in test: {e}")
        return False

if __name__ == "__main__":
    success = test_oecd_analysis()
    if success:
        print("\n🎉 Test completed successfully!")
    else:
        print("\n❌ Test failed!")
        sys.exit(1)