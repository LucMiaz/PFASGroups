#!/usr/bin/env python3
"""
Simplified PFAS Classification Benchmark focusing on PFASGroups analysis.

This version analyzes the PFASGroups classification results without requiring 
the problematic PFAS-atlas dependencies, while still providing valuable insights.
"""

import pandas as pd
import numpy as np
from typing import List, Dict, Tuple, Any
import json
import ast
import os

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

def parse_expected_groups(expected_groups_str: str) -> List[List[int]]:
    """Parse the expected_groups string from PFASGroups results."""
    try:
        if pd.isna(expected_groups_str) or expected_groups_str == '':
            return []
        # Handle string representation of nested list
        if isinstance(expected_groups_str, str):
            expected_groups_str = expected_groups_str.strip()
            if expected_groups_str.startswith('[[') and expected_groups_str.endswith(']]'):
                return ast.literal_eval(expected_groups_str)
            elif expected_groups_str.startswith('[') and expected_groups_str.endswith(']'):
                # Single list case
                parsed = ast.literal_eval(expected_groups_str)
                return [parsed] if isinstance(parsed[0], int) else parsed
        return []
    except Exception as e:
        print(f"Warning: Failed to parse expected_groups '{expected_groups_str}': {e}")
        return []

def map_pfasgroups_to_functional_names() -> Dict[int, str]:
    """Map PFASGroups group IDs to functional group names."""
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

def analyze_pfasgroups_performance(df: pd.DataFrame) -> Dict[str, Any]:
    """Analyze PFASGroups classification performance."""
    
    total_molecules = len(df)
    
    # Specificity analysis (from the original test)
    specific_correct = df['is_specific'].sum() if 'is_specific' in df.columns else 0
    expected_detected = df['expected_group_detected'].sum() if 'expected_group_detected' in df.columns else 0
    
    # Classification coverage
    has_detection = df['detected_groups_parsed'].apply(lambda x: len(x) > 0).sum()
    has_expected = df['expected_groups_parsed'].apply(lambda x: len(x) > 0).sum()
    
    # Group frequency analysis
    all_detected_groups = []
    all_expected_groups = []
    
    for groups_list in df['detected_groups_parsed']:
        all_detected_groups.extend(groups_list)
    
    for groups_lists in df['expected_groups_parsed']:
        for group_list in groups_lists:
            all_expected_groups.extend(group_list)
    
    detected_counts = {}
    expected_counts = {}
    
    for group in all_detected_groups:
        detected_counts[group] = detected_counts.get(group, 0) + 1
    
    for group in all_expected_groups:
        expected_counts[group] = expected_counts.get(group, 0) + 1
    
    analysis = {
        'total_molecules': int(total_molecules),
        'molecules_with_detections': int(has_detection),
        'molecules_with_expected': int(has_expected),
        'specific_correct': int(specific_correct),
        'expected_detected': int(expected_detected),
        'detection_coverage': float(has_detection / total_molecules * 100),
        'specificity_rate': float(specific_correct / total_molecules * 100) if total_molecules > 0 else 0,
        'detection_accuracy': float(expected_detected / has_expected * 100) if has_expected > 0 else 0,
        'top_detected_groups': dict(sorted(detected_counts.items(), key=lambda x: x[1], reverse=True)[:10]),
        'top_expected_groups': dict(sorted(expected_counts.items(), key=lambda x: x[1], reverse=True)[:10])
    }
    
    return analysis

def generate_pfasgroups_report(df: pd.DataFrame, analysis: Dict[str, Any], 
                              output_path: str = None) -> str:
    """Generate a detailed PFASGroups analysis report."""
    
    group_names = map_pfasgroups_to_functional_names()
    
    report = f"""
# PFASGroups Classification Analysis Report

## Overview
This report analyzes the performance of the PFASGroups classification system based on specificity test results.

## Summary Statistics
- Total molecules analyzed: {analysis['total_molecules']}
- Molecules with detections: {analysis['molecules_with_detections']} ({analysis['detection_coverage']:.1f}%)
- Molecules with expected groups: {analysis['molecules_with_expected']}
- Specificity rate: {analysis['specificity_rate']:.1f}%
- Detection accuracy: {analysis['detection_accuracy']:.1f}%

## Performance Metrics
- **Specific correct**: {analysis['specific_correct']} molecules correctly classified with high specificity
- **Expected detected**: {analysis['expected_detected']} molecules where expected groups were detected
- **Coverage**: {analysis['detection_coverage']:.1f}% of molecules had at least one functional group detected

## Top Detected Functional Groups
"""
    
    # Add top detected groups
    for group_id, count in list(analysis['top_detected_groups'].items())[:10]:
        group_name = group_names.get(int(group_id), f"Group_{group_id}")
        percentage = count / analysis['molecules_with_detections'] * 100 if analysis['molecules_with_detections'] > 0 else 0
        report += f"- **{group_name}** (ID: {group_id}): {count} detections ({percentage:.1f}%)\n"
    
    report += "\n## Top Expected Functional Groups\n"
    
    # Add top expected groups
    for group_id, count in list(analysis['top_expected_groups'].items())[:10]:
        group_name = group_names.get(int(group_id), f"Group_{group_id}")
        percentage = count / analysis['molecules_with_expected'] * 100 if analysis['molecules_with_expected'] > 0 else 0
        report += f"- **{group_name}** (ID: {group_id}): {count} expected ({percentage:.1f}%)\n"
    
    # Group-by-group comparison
    report += "\n## Detection vs Expected Comparison\n"
    all_groups = set(analysis['top_detected_groups'].keys()) | set(analysis['top_expected_groups'].keys())
    
    for group_id in sorted(all_groups):
        group_name = group_names.get(int(group_id), f"Group_{group_id}")
        detected = analysis['top_detected_groups'].get(group_id, 0)
        expected = analysis['top_expected_groups'].get(group_id, 0)
        
        if expected > 0:
            accuracy = detected / expected * 100
            report += f"- **{group_name}** (ID: {group_id}): {detected}/{expected} detected ({accuracy:.1f}% accuracy)\n"
    
    # Molecules with issues
    problematic = df[df['is_specific'] == False]
    if len(problematic) > 0:
        report += f"\n## Molecules with Specificity Issues\n"
        report += f"- {len(problematic)} molecules ({len(problematic)/len(df)*100:.1f}%) had specificity issues\n"
        
        # Group the issues by origin
        origin_counts = problematic['origin'].value_counts()
        report += "\n### Issues by Functional Group Type:\n"
        for origin, count in origin_counts.head(10).items():
            report += f"- {origin}: {count} molecules\n"
    
    report += f"""
## Data Quality
- Valid SMILES: {df['valid_smiles'].sum()} / {len(df)} ({df['valid_smiles'].sum()/len(df)*100:.1f}%)
- Molecules with errors: {df['error'].notna().sum()} ({df['error'].notna().sum()/len(df)*100:.1f}%)

## Notes
- Specificity issues indicate that the classification detected more groups than expected
- High specificity is desired to avoid false positive functional group assignments
- Detection accuracy shows how well the system finds expected functional groups
"""
    
    if output_path:
        with open(output_path, 'w') as f:
            f.write(report)
        print(f"✓ Report saved to {output_path}")
    
    return report

def main():
    """Main analysis function."""
    print("🧪 PFASGroups Classification Analysis")
    print("=" * 40)
    
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
    
    # Process PFASGroups classifications
    print("\n🔄 Processing PFASGroups classifications...")
    df['detected_groups_parsed'] = df['detected_groups'].apply(parse_pfasgroups_classification)
    df['expected_groups_parsed'] = df['expected_groups'].apply(parse_expected_groups)
    
    # Analyze performance
    print("\n📊 Analyzing PFASGroups performance...")
    analysis = analyze_pfasgroups_performance(df)
    
    # Generate report
    print("\n📄 Generating analysis report...")
    report_path = os.path.join(output_dir, "pfasgroups_analysis_report.md")
    report = generate_pfasgroups_report(df, analysis, report_path)
    
    # Save detailed results
    results_path = os.path.join(output_dir, "pfasgroups_analysis_results.csv")
    df.to_csv(results_path, index=False)
    print(f"✓ Detailed results saved to {results_path}")
    
    # Save summary (convert to regular Python types for JSON)
    summary_path = os.path.join(output_dir, "pfasgroups_analysis_summary.json")
    analysis_json = {}
    for key, value in analysis.items():
        if isinstance(value, dict):
            # Convert dict keys and values to regular Python types
            analysis_json[key] = {str(k): int(v) if isinstance(v, (np.integer, np.int64)) else v 
                                 for k, v in value.items()}
        elif isinstance(value, (np.integer, np.int64)):
            analysis_json[key] = int(value)
        elif isinstance(value, (np.floating, np.float64)):
            analysis_json[key] = float(value)
        else:
            analysis_json[key] = value
    
    with open(summary_path, 'w') as f:
        json.dump(analysis_json, f, indent=2)
    print(f"✓ Summary statistics saved to {summary_path}")
    
    # Print summary
    print("\n📋 Analysis Summary:")
    print("=" * 30)
    print(f"Total molecules: {analysis['total_molecules']}")
    print(f"Detection coverage: {analysis['detection_coverage']:.1f}%")
    print(f"Specificity rate: {analysis['specificity_rate']:.1f}%")
    print(f"Detection accuracy: {analysis['detection_accuracy']:.1f}%")
    
    print("\n✅ Analysis complete!")

if __name__ == "__main__":
    main()