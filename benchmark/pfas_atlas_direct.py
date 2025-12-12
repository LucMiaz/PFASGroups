#!/usr/bin/env python3
"""
Direct PFAS-atlas benchmark using working import.

Since we know RDKit imports work but subprocess JSON parsing fails,
this version uses direct imports within the same process.
"""

import sys
import os
import pandas as pd
import json
import ast

# Add PFAS-atlas to path
sys.path.insert(0, '/home/luc/git/PFAS-atlas')

def test_pfas_atlas_import():
    """Test if we can import and use PFAS-atlas directly."""
    try:
        from classification_helper import classify_pfas_molecule
        # Test with a simple molecule
        test_result = classify_pfas_molecule("FC(F)(F)C(F)(F)F")
        print(f"✓ PFAS-atlas import successful, test result: {test_result}")
        return classify_pfas_molecule, True
    except Exception as e:
        print(f"✗ PFAS-atlas import failed: {e}")
        return None, False

def classify_molecules_direct(classify_func, smiles_list):
    """Classify molecules directly using imported function."""
    if classify_func is None:
        return [("Unavailable", "Unavailable")] * len(smiles_list)
    
    results = []
    successful = 0
    failed = 0
    
    for i, smiles in enumerate(smiles_list):
        try:
            if pd.isna(smiles) or smiles == '':
                results.append(("Invalid", "Invalid"))
                continue
                
            result = classify_func(smiles)
            if result and len(result) >= 2:
                first_class = str(result[0]) if result[0] else "Unclassified"
                second_class = str(result[1]) if result[1] else "Unclassified"
                results.append((first_class, second_class))
                successful += 1
            else:
                results.append(("Unclassified", "Unclassified"))
                
        except Exception as e:
            results.append(("Error", "Error"))
            failed += 1
            if failed <= 5:  # Only show first 5 errors
                print(f"Warning: Classification failed for molecule {i}: {e}")
        
        # Progress indicator
        if (i + 1) % 200 == 0:
            print(f"  Classified {i + 1}/{len(smiles_list)} molecules (Success: {successful}, Failed: {failed})...")
    
    print(f"✓ Classification complete: {successful} successful, {failed} failed")
    return results

def main():
    """Run the benchmark with direct classification."""
    print("🧪 PFAS Classification Benchmark (Direct Import)")
    print("=" * 50)
    
    # Test PFAS-atlas import
    classify_func, atlas_available = test_pfas_atlas_import()
    
    # Load data
    csv_path = "/home/luc/git/PFASGroups/PFASgroups/tests/results/specificity_test_results.csv"
    df = pd.read_csv(csv_path)
    print(f"✓ Loaded {len(df)} molecules")
    
    # Process PFASGroups results
    df['detected_groups_parsed'] = df['detected_groups'].apply(
        lambda x: ast.literal_eval(x) if pd.notna(x) and x.strip() else []
    )
    df['pfasgroups_classified'] = df['detected_groups_parsed'].apply(lambda x: len(x) > 0)
    
    if atlas_available:
        print("🤖 Running direct PFAS-atlas classification...")
        smiles_list = df['smiles'].tolist()
        atlas_results = classify_molecules_direct(classify_func, smiles_list)
        
        df['atlas_first_class'] = [r[0] for r in atlas_results]
        df['atlas_second_class'] = [r[1] for r in atlas_results]
        df['atlas_classified'] = df['atlas_first_class'].apply(
            lambda x: x not in ['Invalid', 'Error', 'Unclassified', 'Unavailable']
        )
    else:
        df['atlas_first_class'] = 'Unavailable'
        df['atlas_second_class'] = 'Unavailable'  
        df['atlas_classified'] = False
    
    # Compute statistics
    total = len(df)
    pfasgroups_count = int(df['pfasgroups_classified'].sum())
    
    stats = {
        'total_molecules': total,
        'pfasgroups_classified': pfasgroups_count,
        'pfasgroups_coverage': float(pfasgroups_count / total * 100)
    }
    
    if atlas_available:
        atlas_count = int(df['atlas_classified'].sum())
        both_count = int(((df['pfasgroups_classified']) & (df['atlas_classified'])).sum())
        
        stats.update({
            'atlas_classified': atlas_count,
            'atlas_coverage': float(atlas_count / total * 100),
            'both_classified': both_count,
            'agreement_rate': float(both_count / max(pfasgroups_count, atlas_count, 1) * 100)
        })
    
    # Generate summary report
    print("\n📊 Results Summary:")
    print(f"Total molecules: {stats['total_molecules']}")
    print(f"PFASGroups classified: {stats['pfasgroups_classified']} ({stats['pfasgroups_coverage']:.1f}%)")
    
    if atlas_available:
        print(f"PFAS-atlas classified: {stats['atlas_classified']} ({stats['atlas_coverage']:.1f}%)")
        print(f"Both systems classified: {stats['both_classified']}")
        print(f"Agreement rate: {stats['agreement_rate']:.1f}%")
        
        # Top classifications
        print("\\nTop PFAS-atlas classifications:")
        atlas_counts = df['atlas_first_class'].value_counts().head(10)
        for class_name, count in atlas_counts.items():
            if class_name not in ['Invalid', 'Error', 'Unclassified', 'Unavailable']:
                print(f"  {class_name}: {count} molecules")
    
    # Save results
    output_dir = "/home/luc/git/PFASGroups/benchmark"
    os.makedirs(output_dir, exist_ok=True)
    
    results_path = os.path.join(output_dir, "direct_benchmark_results.csv")
    df.to_csv(results_path, index=False)
    print(f"\\n✓ Results saved to {results_path}")
    
    summary_path = os.path.join(output_dir, "direct_benchmark_summary.json")
    with open(summary_path, 'w') as f:
        json.dump(stats, f, indent=2)
    print(f"✓ Summary saved to {summary_path}")
    
    print("\\n✅ Direct benchmark complete!")

if __name__ == "__main__":
    main()