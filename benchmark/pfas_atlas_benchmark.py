#!/usr/bin/env python3
"""
Benchmark test comparing PFASGroups classification with PFAS-atlas classification.

This script reads molecules from the specificity test results and compares:
- PFASGroups classification results (from CSV file)
- PFAS-atlas classification results (using their model)

The comparison will show agreement/disagreement between the two classification systems.
"""

import sys
import os
import pandas as pd
import numpy as np
from typing import List, Dict, Tuple, Any
import json
import ast
from pathlib import Path

# Add PFAS-atlas to path to import their classification module
sys.path.insert(0, '/home/luc/git/PFAS-atlas')
sys.path.insert(0, '/home/luc/git/PFAS-atlas/classification_helper')

try:
    from classification_helper import classify_pfas_molecule
    from step3_classify import get_classify
    PFAS_ATLAS_AVAILABLE = True
    print("✓ PFAS-atlas classification module loaded successfully")
except ImportError as e:
    PFAS_ATLAS_AVAILABLE = False
    print(f"✗ Failed to import PFAS-atlas classification: {e}")
    print("Make sure you're running this in the 'pfasatlas' environment")

def load_specificity_test_data(csv_path: str) -> pd.DataFrame:
    """Load the specificity test results CSV file."""
    try:
        df = pd.read_csv(csv_path)
        print(f"✓ Loaded {len(df)} molecules from {csv_path}")
        return df
    except Exception as e:
        print(f"✗ Failed to load CSV file: {e}")
        return None

def parse_pfasgroups_classification(detected_groups_str: str) -> List[int]:
    """Parse the detected_groups string from PFASGroups results."""
    try:
        if pd.isna(detected_groups_str) or detected_groups_str == '':
            return []
        # Handle string representation of list
        if isinstance(detected_groups_str, str):
            # Remove extra spaces and parse as list
            detected_groups_str = detected_groups_str.strip()
            if detected_groups_str.startswith('[') and detected_groups_str.endswith(']'):
                return ast.literal_eval(detected_groups_str)
            else:
                # Try to split by comma if it's a simple comma-separated string
                return [int(x.strip()) for x in detected_groups_str.split(',') if x.strip()]
        elif isinstance(detected_groups_str, list):
            return detected_groups_str
        else:
            return []
    except Exception as e:
        print(f"Warning: Failed to parse detected_groups '{detected_groups_str}': {e}")
        return []

def classify_with_pfas_atlas(smiles_list: List[str]) -> List[Tuple[str, str]]:
    """Classify molecules using PFAS-atlas classification system."""
    if not PFAS_ATLAS_AVAILABLE:
        print("✗ PFAS-atlas not available, returning empty classifications")
        return [("Unknown", "Unknown")] * len(smiles_list)
    
    classifications = []
    failed_count = 0
    
    for i, smiles in enumerate(smiles_list):
        try:
            if pd.isna(smiles) or smiles == '':
                classifications.append(("Invalid", "Invalid"))
                failed_count += 1
                continue
                
            # Use the PFAS-atlas classification function
            result = classify_pfas_molecule(smiles)
            if result and len(result) >= 2:
                first_class = result[0] if result[0] else "Unclassified"
                second_class = result[1] if result[1] else "Unclassified"
                classifications.append((first_class, second_class))
            else:
                classifications.append(("Unclassified", "Unclassified"))
                failed_count += 1
                
        except Exception as e:
            print(f"Warning: Failed to classify SMILES '{smiles}': {e}")
            classifications.append(("Error", "Error"))
            failed_count += 1
            
        # Progress indicator
        if (i + 1) % 100 == 0:
            print(f"  Classified {i + 1}/{len(smiles_list)} molecules...")
    
    print(f"✓ PFAS-atlas classification complete. {failed_count}/{len(smiles_list)} failed.")
    return classifications

def map_pfasgroups_to_functional_names() -> Dict[int, str]:
    """Map PFASGroups group IDs to functional group names."""
    # This mapping would need to be updated based on the actual PFASGroups definitions
    # For now, using common PFAS functional groups
    return {
        18: "ethene",
        19: "alkene", 
        20: "alkane",
        21: "halide",
        23: "alkane-per",
        25: "iodide-per", 
        27: "ketone",
        28: "ketone-per",
        30: "Perfluoroalkyl ketones-per",
        35: "acyl halide-per",
        41: "alkene-per",
        42: "iodide-per",
        48: "Perfluoroalkane-per",
        49: "Perfluoroalkene-per",
        50: "alkyne-per"
    }

def analyze_classification_agreement(df_results: pd.DataFrame) -> Dict[str, Any]:
    """Analyze agreement between PFASGroups and PFAS-atlas classifications."""
    
    # Count molecules by classification type
    pfasgroups_classified = df_results['pfasgroups_has_classification'].sum()
    pfas_atlas_classified = df_results['pfas_atlas_classified'].sum()
    total_molecules = len(df_results)
    
    # Agreement analysis
    both_classified = ((df_results['pfasgroups_has_classification']) & 
                      (df_results['pfas_atlas_classified'])).sum()
    
    neither_classified = ((~df_results['pfasgroups_has_classification']) & 
                         (~df_results['pfas_atlas_classified'])).sum()
    
    pfasgroups_only = ((df_results['pfasgroups_has_classification']) & 
                      (~df_results['pfas_atlas_classified'])).sum()
    
    pfas_atlas_only = ((~df_results['pfasgroups_has_classification']) & 
                      (df_results['pfas_atlas_classified'])).sum()
    
    analysis = {
        'total_molecules': total_molecules,
        'pfasgroups_classified': pfasgroups_classified,
        'pfas_atlas_classified': pfas_atlas_classified,
        'both_classified': both_classified,
        'neither_classified': neither_classified,
        'pfasgroups_only': pfasgroups_only,
        'pfas_atlas_only': pfas_atlas_only,
        'pfasgroups_coverage': pfasgroups_classified / total_molecules * 100,
        'pfas_atlas_coverage': pfas_atlas_classified / total_molecules * 100,
        'agreement_rate': both_classified / max(pfasgroups_classified, pfas_atlas_classified, 1) * 100
    }
    
    return analysis

def generate_benchmark_report(df_results: pd.DataFrame, analysis: Dict[str, Any], 
                            output_path: str = None) -> str:
    """Generate a detailed benchmark report."""
    
    report = f"""
# PFAS Classification Benchmark Report
## PFASGroups vs PFAS-atlas

### Summary Statistics
- Total molecules analyzed: {analysis['total_molecules']}
- PFASGroups classified: {analysis['pfasgroups_classified']} ({analysis['pfasgroups_coverage']:.1f}%)
- PFAS-atlas classified: {analysis['pfas_atlas_classified']} ({analysis['pfas_atlas_coverage']:.1f}%)

### Classification Agreement
- Both systems classified: {analysis['both_classified']}
- Neither system classified: {analysis['neither_classified']}
- Only PFASGroups classified: {analysis['pfasgroups_only']}
- Only PFAS-atlas classified: {analysis['pfas_atlas_only']}
- Agreement rate: {analysis['agreement_rate']:.1f}%

### Top PFASGroups Classifications
"""
    
    # Add top PFASGroups classifications
    pfasgroups_counts = {}
    for _, row in df_results.iterrows():
        if row['pfasgroups_has_classification']:
            groups = row['detected_groups_parsed']
            for group in groups:
                pfasgroups_counts[group] = pfasgroups_counts.get(group, 0) + 1
    
    if pfasgroups_counts:
        group_names = map_pfasgroups_to_functional_names()
        report += "\n"
        for group_id, count in sorted(pfasgroups_counts.items(), key=lambda x: x[1], reverse=True)[:10]:
            group_name = group_names.get(group_id, f"Group_{group_id}")
            report += f"- {group_name} (ID: {group_id}): {count} molecules\n"
    
    # Add top PFAS-atlas classifications
    report += "\n### Top PFAS-atlas Classifications\n"
    pfas_atlas_first = df_results['pfas_atlas_first_class'].value_counts().head(10)
    for class_name, count in pfas_atlas_first.items():
        if class_name not in ['Invalid', 'Error', 'Unclassified']:
            report += f"- {class_name}: {count} molecules\n"
    
    if output_path:
        with open(output_path, 'w') as f:
            f.write(report)
        print(f"✓ Report saved to {output_path}")
    
    return report

def main():
    """Main benchmark function."""
    print("🧪 PFAS Classification Benchmark: PFASGroups vs PFAS-atlas")
    print("=" * 60)
    
    # File paths
    csv_path = "/home/luc/git/PFASGroups/PFASgroups/tests/results/specificity_test_results.csv"
    output_dir = "/home/luc/git/PFASGroups/benchmark"
    
    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)
    
    # Load data
    print("\n📂 Loading data...")
    df = load_specificity_test_data(csv_path)
    if df is None:
        return
    
    # Prepare data
    print("\n🔄 Processing PFASGroups classifications...")
    df['detected_groups_parsed'] = df['detected_groups'].apply(parse_pfasgroups_classification)
    df['pfasgroups_has_classification'] = df['detected_groups_parsed'].apply(lambda x: len(x) > 0)
    
    # Classify with PFAS-atlas
    print("\n🤖 Running PFAS-atlas classification...")
    smiles_list = df['smiles'].tolist()
    pfas_atlas_results = classify_with_pfas_atlas(smiles_list)
    
    # Add PFAS-atlas results to dataframe
    df['pfas_atlas_first_class'] = [result[0] for result in pfas_atlas_results]
    df['pfas_atlas_second_class'] = [result[1] for result in pfas_atlas_results]
    df['pfas_atlas_classified'] = df['pfas_atlas_first_class'].apply(
        lambda x: x not in ['Invalid', 'Error', 'Unclassified', 'Unknown']
    )
    
    # Analyze results
    print("\n📊 Analyzing classification agreement...")
    analysis = analyze_classification_agreement(df)
    
    # Generate report
    print("\n📄 Generating benchmark report...")
    report_path = os.path.join(output_dir, "pfas_atlas_benchmark_report.md")
    report = generate_benchmark_report(df, analysis, report_path)
    
    # Save detailed results
    results_path = os.path.join(output_dir, "pfas_atlas_benchmark_results.csv")
    df.to_csv(results_path, index=False)
    print(f"✓ Detailed results saved to {results_path}")
    
    # Save summary (convert numpy int64 to regular int for JSON serialization)
    summary_path = os.path.join(output_dir, "pfas_atlas_benchmark_summary.json")
    # Convert numpy types to native Python types for JSON serialization
    analysis_json = {}
    for key, value in analysis.items():
        if hasattr(value, 'item'):  # numpy scalar
            analysis_json[key] = value.item()
        elif isinstance(value, np.integer):
            analysis_json[key] = int(value)
        elif isinstance(value, np.floating):
            analysis_json[key] = float(value)
        else:
            analysis_json[key] = value
    
    with open(summary_path, 'w') as f:
        json.dump(analysis_json, f, indent=2)
    print(f"✓ Summary statistics saved to {summary_path}")
    
    # Print summary
    print("\n📋 Benchmark Summary:")
    print("=" * 40)
    print(f"Total molecules: {analysis['total_molecules']}")
    print(f"PFASGroups coverage: {analysis['pfasgroups_coverage']:.1f}%")
    print(f"PFAS-atlas coverage: {analysis['pfas_atlas_coverage']:.1f}%")
    print(f"Agreement rate: {analysis['agreement_rate']:.1f}%")
    
    print("\n✅ Benchmark complete!")

if __name__ == "__main__":
    main()