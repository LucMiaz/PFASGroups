"""
Export misclassified molecules from specificity test results.
Creates a CSV file that can be used with draw_misclassified.py to visualize problematic molecules.
"""

import pandas as pd
import os
import sys

# Add parent directory to path for imports
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

def export_misclassified_molecules(
    input_file='results/specificity_test_results.csv',
    output_file='results/misclassified_molecules.csv',
    detection_failures_only=False,
    specificity_failures_only=False
):
    """
    Export misclassified molecules from specificity test results.
    
    Args:
        input_file: Path to specificity test results CSV
        output_file: Path to output misclassified molecules CSV
        detection_failures_only: Only export molecules where expected groups were not detected
        specificity_failures_only: Only export molecules with extra unexpected groups
    
    Returns:
        DataFrame with misclassified molecules
    """
    
    # Get the directory of this script
    script_dir = os.path.dirname(os.path.abspath(__file__))
    
    # Build full paths
    input_path = os.path.join(script_dir, input_file)
    output_path = os.path.join(script_dir, output_file)
    
    # Load specificity test results
    print(f"Loading specificity test results from: {input_path}")
    df = pd.read_csv(input_path)
    
    print(f"Total test molecules: {len(df)}")
    print(f"Valid SMILES: {df['valid_smiles'].sum()}")
    
    # Filter for misclassified molecules
    valid_df = df[df['valid_smiles'] == True].copy()
    
    # Identify different types of failures
    detection_failures = valid_df[valid_df['expected_group_detected'] == False]
    specificity_failures = valid_df[valid_df['is_specific'] == False]
    
    print(f"\nDetection failures: {len(detection_failures)} (expected groups not detected)")
    print(f"Specificity failures: {len(specificity_failures)} (extra unexpected groups detected)")
    
    # Apply filters if specified
    if detection_failures_only:
        misclassified = detection_failures
        print("\nExporting detection failures only")
    elif specificity_failures_only:
        misclassified = specificity_failures
        print("\nExporting specificity failures only")
    else:
        # Export all misclassified (either type of failure)
        misclassified = valid_df[
            (valid_df['expected_group_detected'] == False) | 
            (valid_df['is_specific'] == False)
        ]
        print("\nExporting all misclassified molecules")
    
    print(f"Total misclassified molecules to export: {len(misclassified)}")
    
    # Parse expected and detected groups from string representation
    def parse_groups(group_str):
        """Parse group list from string representation like '[1, 2, 3]'"""
        if pd.isna(group_str) or group_str == '[]':
            return []
        try:
            # Remove brackets and split by comma
            cleaned = str(group_str).strip('[]')
            if not cleaned:
                return []
            return [int(x.strip()) for x in cleaned.split(',')]
        except:
            return []
    
    # Create export dataframe with useful columns
    export_df = pd.DataFrame({
        'smiles': misclassified['smiles'],
        'inchi': misclassified['inchi'],
        'inchikey': misclassified['inchikey'],
        'formula': misclassified['formula'],
        'origin': misclassified['origin'],
        'pathtype': misclassified['pathtype'],
        'expected_groups': misclassified['expected_groups'].apply(parse_groups),
        'detected_groups': misclassified['detected_groups'].apply(parse_groups),
        'n_expected': misclassified['n_expected_groups'],
        'n_detected': misclassified['n_detected_groups'],
        'expected_detected': misclassified['expected_group_detected'],
        'is_specific': misclassified['is_specific'],
        'failure_type': misclassified.apply(
            lambda row: 'detection' if not row['expected_group_detected'] 
                       else ('specificity' if not row['is_specific'] else 'none'),
            axis=1
        )
    })
    
    # Sort by origin for easier analysis
    export_df = export_df.sort_values(['origin', 'pathtype', 'failure_type'])
    
    # Save to CSV
    export_df.to_csv(output_path, index=False)
    print(f"\nMisclassified molecules exported to: {output_path}")
    
    # Print summary statistics
    print("\n" + "="*60)
    print("MISCLASSIFICATION SUMMARY")
    print("="*60)
    
    print("\nBy failure type:")
    print(export_df['failure_type'].value_counts())
    
    print("\nBy origin (top 20):")
    origin_counts = export_df['origin'].value_counts()
    print(origin_counts.head(20))
    
    print("\nBy pathtype:")
    print(export_df['pathtype'].value_counts())
    
    # Identify most problematic origins
    print("\n" + "-"*60)
    print("MOST PROBLEMATIC GROUPS")
    print("-"*60)
    
    # Detection failures by origin
    detection_by_origin = export_df[export_df['failure_type'] == 'detection']['origin'].value_counts()
    if len(detection_by_origin) > 0:
        print("\nMost detection failures:")
        print(detection_by_origin.head(10))
    
    # Specificity failures by origin
    specificity_by_origin = export_df[export_df['failure_type'] == 'specificity']['origin'].value_counts()
    if len(specificity_by_origin) > 0:
        print("\nMost specificity failures:")
        print(specificity_by_origin.head(10))
    
    print("\n" + "="*60)
    print(f"\nTo visualize these molecules, run:")
    print(f"  cd {script_dir}")
    print(f"  python results/draw_misclassified.py")
    print("="*60)
    
    return export_df


if __name__ == "__main__":
    import argparse
    
    parser = argparse.ArgumentParser(
        description='Export misclassified molecules from specificity test results'
    )
    parser.add_argument(
        '--input',
        default='results/specificity_test_results.csv',
        help='Input specificity test results CSV file'
    )
    parser.add_argument(
        '--output',
        default='results/misclassified_molecules.csv',
        help='Output misclassified molecules CSV file'
    )
    parser.add_argument(
        '--detection-only',
        action='store_true',
        help='Only export detection failures (expected groups not detected)'
    )
    parser.add_argument(
        '--specificity-only',
        action='store_true',
        help='Only export specificity failures (unexpected extra groups)'
    )
    
    args = parser.parse_args()
    
    # Export misclassified molecules
    export_misclassified_molecules(
        input_file=args.input,
        output_file=args.output,
        detection_failures_only=args.detection_only,
        specificity_failures_only=args.specificity_only
    )
