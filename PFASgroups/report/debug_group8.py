"""
Debug script to analyze group 8 detection issue without requiring RDKit
"""

import csv
import json

def analyze_group8_detection():
    """Analyze group 8 detection patterns in the test results"""
    
    # Read the test results
    results_file = '../tests/specificity_test_results.csv'
    
    try:
        with open(results_file, 'r') as f:
            reader = csv.DictReader(f)
            
            print("🔍 ANALYZING GROUP 8 DETECTION PATTERNS")
            print("=" * 60)
            
            group8_molecules = []
            target_molecule_found = False
            target_smiles = "O=S(=O)(O)C(F)(F)C(OC(F)(F)F)(C(F)(F)C(F)(F)F)C(C(F)(F)F)(C(F)(F)F)C(F)(F)C(F)(F)F"
            
            for row in reader:
                smiles = row['smiles']
                expected_groups = eval(row['group_ids'])  # Parse the list
                detected_groups = eval(row['detected_groups'])  # Parse the list
                
                # Check for target molecule
                if smiles == target_smiles:
                    target_molecule_found = True
                    print(f"🎯 TARGET MOLECULE FOUND:")
                    print(f"   SMILES: {smiles}")
                    print(f"   Expected: {expected_groups}")
                    print(f"   Detected: {detected_groups}")
                    print(f"   Has Group 8: {8 in detected_groups}")
                    print()
                
                # Collect molecules with group 8
                if 8 in detected_groups:
                    group8_molecules.append({
                        'smiles': smiles,
                        'expected': expected_groups,
                        'detected': detected_groups,
                        'sulfonic_count': smiles.count('S(=O)(=O)O')
                    })
            
            print(f"📊 GROUP 8 DETECTION SUMMARY:")
            print(f"   Total molecules with Group 8: {len(group8_molecules)}")
            print(f"   Target molecule has Group 8: {not target_molecule_found or 8 in detected_groups if target_molecule_found else 'Not found'}")
            print()
            
            if group8_molecules:
                print("🧬 MOLECULES WITH GROUP 8 DETECTED:")
                print("-" * 40)
                for i, mol in enumerate(group8_molecules[:5]):  # Show first 5
                    print(f"{i+1}. Sulfonic groups: {mol['sulfonic_count']}")
                    print(f"   SMILES: {mol['smiles'][:80]}...")
                    print(f"   Expected: {mol['expected']}")
                    print(f"   Detected: {mol['detected']}")
                    print()
                    
                if len(group8_molecules) > 5:
                    print(f"   ... and {len(group8_molecules) - 5} more")
                    
                # Analyze sulfonic acid count patterns
                sulfonic_counts = [mol['sulfonic_count'] for mol in group8_molecules]
                print(f"📈 SULFONIC ACID GROUP COUNTS IN GROUP 8 MOLECULES:")
                for count in set(sulfonic_counts):
                    num_mols = sulfonic_counts.count(count)
                    print(f"   {count} sulfonic groups: {num_mols} molecules")
            
            if not target_molecule_found:
                print("❌ Target molecule not found in results!")
                
    except FileNotFoundError:
        print(f"❌ Could not find results file: {results_file}")
    except Exception as e:
        print(f"❌ Error analyzing results: {e}")

if __name__ == "__main__":
    analyze_group8_detection()