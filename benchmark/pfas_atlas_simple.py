#!/usr/bin/env python3
"""
Simple PFAS-atlas benchmark using direct import approach.

This script attempts to import and use PFAS-atlas directly,
with fallback to PFASGroups-only analysis.
"""

import sys
import os
import pandas as pd
import json
import ast

# Add PFAS-atlas to Python path
sys.path.insert(0, '/home/luc/git/PFAS-atlas')

def try_import_pfas_atlas():
    """Try to import PFAS-atlas classification function."""
    try:
        from classification_helper import classify_pfas_molecule
        return classify_pfas_molecule, True
    except Exception as e:
        print(f"Warning: Could not import PFAS-atlas: {e}")
        return None, False

def classify_with_pfas_atlas_safe(classify_func, smiles_list):
    """Safely classify molecules with PFAS-atlas, handling errors."""
    if classify_func is None:
        return [("Unavailable", "Unavailable")] * len(smiles_list)
    
    results = []
    failed_count = 0
    
    for i, smiles in enumerate(smiles_list):
        try:
            if pd.isna(smiles) or smiles == '':
                results.append(("Invalid", "Invalid"))
                continue
            
            result = classify_func(smiles)
            if result and len(result) >= 2:
                first_class = result[0] if result[0] else "Unclassified"
                second_class = result[1] if result[1] else "Unclassified"
                results.append((first_class, second_class))
            else:
                results.append(("Unclassified", "Unclassified"))
                
        except Exception as e:
            results.append(("Error", "Error"))
            failed_count += 1
            if failed_count < 5:  # Only show first few errors
                print(f"Warning: Classification failed for '{smiles}': {e}")
        
        # Progress indicator
        if (i + 1) % 200 == 0:
            print(f"  Classified {i + 1}/{len(smiles_list)} molecules...")
    
    return results

def main():
    """Run the comparison benchmark."""
    print("🧪 PFAS Classification Comparison")
    print("=" * 40)
    
    # Try to import PFAS-atlas
    classify_func, atlas_available = try_import_pfas_atlas()
    
    if atlas_available:
        print("✓ PFAS-atlas imported successfully")
    else:
        print("⚠️  PFAS-atlas not available - running PFASGroups analysis only")
    
    # Load data
    csv_path = "/home/luc/git/PFASGroups/PFASgroups/tests/results/specificity_test_results.csv"
    df = pd.read_csv(csv_path)
    print(f"✓ Loaded {len(df)} molecules")
    
    # Process PFASGroups results
    df['detected_groups_parsed'] = df['detected_groups'].apply(
        lambda x: ast.literal_eval(x) if pd.notna(x) and x.strip() else []
    )
    df['pfasgroups_classified'] = df['detected_groups_parsed'].apply(lambda x: len(x) > 0)
    
    # Run PFAS-atlas if available
    if atlas_available:
        print("\n🤖 Running PFAS-atlas classification...")
        smiles_list = df['smiles'].tolist()
        atlas_results = classify_with_pfas_atlas_safe(classify_func, smiles_list)
        
        df['atlas_first_class'] = [r[0] for r in atlas_results]
        df['atlas_second_class'] = [r[1] for r in atlas_results]
        df['atlas_classified'] = df['atlas_first_class'].apply(
            lambda x: x not in ['Invalid', 'Error', 'Unclassified', 'Unavailable']
        )
    else:
        df['atlas_first_class'] = 'Unavailable'
        df['atlas_second_class'] = 'Unavailable'
        df['atlas_classified'] = False
    
    # Analysis
    print("\n📊 Computing statistics...")
    
    total = len(df)
    pfasgroups_count = df['pfasgroups_classified'].sum()
    
    stats = {
        'total_molecules': total,
        'pfasgroups_classified': int(pfasgroups_count),
        'pfasgroups_coverage': float(pfasgroups_count / total * 100)
    }
    
    if atlas_available:
        atlas_count = df['atlas_classified'].sum()
        both_count = ((df['pfasgroups_classified']) & (df['atlas_classified'])).sum()
        
        stats.update({
            'atlas_classified': int(atlas_count),
            'atlas_coverage': float(atlas_count / total * 100),
            'both_classified': int(both_count),
            'agreement_rate': float(both_count / max(pfasgroups_count, atlas_count, 1) * 100)
        })
    
    # Generate report
    output_dir = "/home/luc/git/PFASGroups/benchmark"
    os.makedirs(output_dir, exist_ok=True)
    
    report = f"""# PFAS Classification Benchmark Results

## Summary
- **Total molecules**: {stats['total_molecules']}
- **PFASGroups classified**: {stats['pfasgroups_classified']} ({stats['pfasgroups_coverage']:.1f}%)
"""
    
    if atlas_available:
        report += f"""- **PFAS-atlas classified**: {stats['atlas_classified']} ({stats['atlas_coverage']:.1f}%)
- **Both systems classified**: {stats['both_classified']}
- **Agreement rate**: {stats['agreement_rate']:.1f}%

## Top PFAS-atlas Classifications
"""
        atlas_counts = df['atlas_first_class'].value_counts().head(10)
        for class_name, count in atlas_counts.items():
            if class_name not in ['Invalid', 'Error', 'Unclassified', 'Unavailable']:
                report += f"- {class_name}: {count} molecules\n"
    
    else:
        report += "\n- **PFAS-atlas**: Not available\n"
    
    # Top PFASGroups
    group_names = {48: "Perfluoroalkane-per", 23: "alkane-per", 21: "halide", 
                   41: "alkene-per", 42: "iodide-per", 35: "acyl halide-per"}
    
    all_groups = []
    for groups in df['detected_groups_parsed']:
        all_groups.extend(groups)
    
    group_counts = {}
    for group in all_groups:
        group_counts[group] = group_counts.get(group, 0) + 1
    
    report += f"""
## Top PFASGroups Classifications
"""
    for group_id, count in sorted(group_counts.items(), key=lambda x: x[1], reverse=True)[:10]:
        group_name = group_names.get(group_id, f"Group_{group_id}")
        report += f"- {group_name} (ID: {group_id}): {count} detections\n"
    
    # Save results
    report_path = os.path.join(output_dir, "benchmark_final_report.md")
    with open(report_path, 'w') as f:
        f.write(report)
    
    results_path = os.path.join(output_dir, "benchmark_final_results.csv")
    df.to_csv(results_path, index=False)
    
    summary_path = os.path.join(output_dir, "benchmark_final_summary.json")
    with open(summary_path, 'w') as f:
        json.dump(stats, f, indent=2)
    
    # Print summary
    print("\n📋 Final Results:")
    print(f"Total molecules: {stats['total_molecules']}")
    print(f"PFASGroups coverage: {stats['pfasgroups_coverage']:.1f}%")
    if atlas_available:
        print(f"PFAS-atlas coverage: {stats['atlas_coverage']:.1f}%")
        print(f"Agreement rate: {stats['agreement_rate']:.1f}%")
    
    print(f"\n✓ Results saved to {output_dir}")
    print("✅ Benchmark complete!")

if __name__ == "__main__":
    main()