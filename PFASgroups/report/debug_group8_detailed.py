"""
Advanced debugging script to understand group 8 detection issue
"""

import csv

def analyze_problematic_molecules():
    """Analyze the specific molecules that shouldn't have group 8"""
    
    results_file = '../tests/specificity_test_results.csv'
    
    try:
        with open(results_file, 'r') as f:
            reader = csv.DictReader(f)
            
            print("🔍 DETAILED GROUP 8 PROBLEM ANALYSIS")
            print("=" * 60)
            
            problematic_cases = []
            
            for row in reader:
                smiles = row['smiles']
                expected_groups = eval(row['group_ids'])
                detected_groups = eval(row['detected_groups'])
                
                # Check if group 8 is detected but not expected
                if 8 in detected_groups and 8 not in expected_groups:
                    sulfonic_count = smiles.count('S(=O)(=O)O')
                    problematic_cases.append({
                        'smiles': smiles,
                        'expected': expected_groups,
                        'detected': detected_groups,
                        'sulfonic_count': sulfonic_count,
                        'has_ether': 'O' in smiles and ('OC(' in smiles or ')O' in smiles)
                    })
            
            print(f"📊 FOUND {len(problematic_cases)} PROBLEMATIC CASES")
            print("   (Group 8 detected but not expected)")
            print()
            
            if problematic_cases:
                print("🧬 DETAILED ANALYSIS OF FIRST 10 CASES:")
                print("-" * 50)
                
                for i, case in enumerate(problematic_cases[:10]):
                    print(f"\n{i+1}. CASE ANALYSIS:")
                    print(f"   SMILES: {case['smiles']}")
                    print(f"   Expected groups: {case['expected']}")
                    print(f"   Detected groups: {case['detected']}")
                    print(f"   Sulfonic acid count: {case['sulfonic_count']}")
                    print(f"   Has ether groups: {case['has_ether']}")
                    
                    # Analyze what other sulfonic acid groups are detected
                    sulfonic_groups = [g for g in case['detected'] if g in [6, 7, 8, 9, 10, 11]]
                    print(f"   All sulfonic-related groups: {sulfonic_groups}")
                    
                    # Count actual S=O patterns
                    so2_count = case['smiles'].count('S(=O)(=O)')
                    print(f"   S(=O)(=O) pattern count: {so2_count}")
                    
                print(f"\n... and {len(problematic_cases) - 10} more cases" if len(problematic_cases) > 10 else "")
                
                # Summary statistics
                sulfonic_counts = [case['sulfonic_count'] for case in problematic_cases]
                so2_counts = [case['smiles'].count('S(=O)(=O)') for case in problematic_cases]
                
                print(f"\n📈 SUMMARY STATISTICS:")
                print(f"   Cases with 0 sulfonic groups: {sulfonic_counts.count(0)}")
                print(f"   Cases with 1 sulfonic group: {sulfonic_counts.count(1)}")
                print(f"   Cases with 2+ sulfonic groups: {sum(1 for x in sulfonic_counts if x >= 2)}")
                print()
                print(f"   Cases with 0 S(=O)(=O) patterns: {so2_counts.count(0)}")
                print(f"   Cases with 1 S(=O)(=O) pattern: {so2_counts.count(1)}")
                print(f"   Cases with 2+ S(=O)(=O) patterns: {sum(1 for x in so2_counts if x >= 2)}")
                
                # Show examples of each type
                print(f"\n🔬 EXAMPLES BY TYPE:")
                
                examples_shown = set()
                for label, condition in [
                    ("0 sulfonic groups", lambda x: x['sulfonic_count'] == 0),
                    ("1 sulfonic group", lambda x: x['sulfonic_count'] == 1),
                    ("2+ sulfonic groups", lambda x: x['sulfonic_count'] >= 2)
                ]:
                    matching = [c for c in problematic_cases if condition(c)]
                    if matching and label not in examples_shown:
                        example = matching[0]
                        print(f"   {label}: {example['smiles'][:60]}...")
                        examples_shown.add(label)
            
            else:
                print("✅ No problematic cases found!")
                
        print(f"\n🎯 CONCLUSION:")
        print(f"   The issue appears to be that Group 8 (disulfonic acids)")
        print(f"   is being detected in molecules that don't have 2 distinct")
        print(f"   sulfonic acid groups. This suggests a bug in the pathway")
        print(f"   validation logic in the core algorithm.")
                
    except FileNotFoundError:
        print(f"❌ Could not find results file: {results_file}")
    except Exception as e:
        print(f"❌ Error analyzing results: {e}")

if __name__ == "__main__":
    analyze_problematic_molecules()