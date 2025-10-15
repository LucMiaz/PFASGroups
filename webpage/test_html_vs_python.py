"""
Test script to compare HTML JavaScript implementation vs Python implementation.

This script:
1. Generates test molecules for all PFAS groups (using test_examples.py)
2. Runs them through the Python parse_PFAS_groups function
3. Exports test molecules to a format the HTML can process
4. Provides instructions for running HTML analysis
5. Compares HTML results with Python results

Usage:
    python test_html_vs_python.py --generate    # Generate test molecules and Python results
    python test_html_vs_python.py --compare results.csv  # Compare HTML output with Python baseline
"""

import os
import sys
import json
import argparse
import pandas as pd
from datetime import datetime
from rdkit import Chem
from rdkit.Chem.rdMolDescriptors import CalcMolFormula

# Add parent directory to path for imports
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

try:
    from PFASgroups.core import parse_PFAS_groups
    from PFASgroups.tests.test_examples import (
        create_specificity_test_molecules,
        OECD_PFAS_GROUPS,
        GENERIC_PFAS_GROUPS,
        IGNORE_GROUPS
    )
    from PFASgroups.generate_mol import (
        generate_random_carbon_chain,
        fluorinate_mol,
        append_functional_group,
        get_attachment
    )
except ImportError as e:
    print(f"Import error: {e}")
    print("Make sure PFASgroups package is installed or in PYTHONPATH")
    sys.exit(1)


def generate_test_molecules_from_groups(n_per_group=5):
    """
    Generate test molecules for OECD and generic PFAS groups.
    
    Args:
        n_per_group: Number of molecules to generate per group
        
    Returns:
        DataFrame with columns: test_id, group_id, group_name, smiles, formula, inchi
    """
    test_molecules = []
    test_id = 1
    
    print("Generating OECD PFAS group test molecules...")
    for group_id, group_name, functional_groups, pathway in OECD_PFAS_GROUPS:
        if group_id in IGNORE_GROUPS:
            continue
            
        print(f"  Group {group_id}: {group_name}")
        for i in range(n_per_group):
            try:
                # Generate random carbon chain
                n_carbons = 4 + (i % 8)  # Vary chain length
                chain = generate_random_carbon_chain(n_carbons, branching_probability=0.2)
                
                # Fluorinate according to pathway
                if pathway == 'Perfluoroalkyl':
                    mol = fluorinate_mol(chain, perfluorinated=True, p=1.0)
                elif pathway == 'Polyfluoroalkyl':
                    mol = fluorinate_mol(chain, perfluorinated=False, p=0.7)
                else:
                    mol = fluorinate_mol(chain, perfluorinated=False, p=0.5)
                
                # Append functional groups
                for fg in functional_groups:
                    group_smiles = fg['group_smiles']
                    n = fg['n'] if isinstance(fg['n'], int) else 1
                    mode = fg.get('mode', 'attach')
                    neighbours = fg.get('neighbours', None)
                    
                    for _ in range(n):
                        attachment_point = get_attachment(mol, neighbours=neighbours)
                        if attachment_point is not None:
                            mol = append_functional_group(
                                mol, 
                                group_smiles, 
                                attachment_point, 
                                mode=mode
                            )
                
                if mol is not None:
                    smiles = Chem.MolToSmiles(mol)
                    inchi = Chem.MolToInchi(mol)
                    formula = CalcMolFormula(mol)
                    
                    test_molecules.append({
                        'test_id': test_id,
                        'group_id': group_id,
                        'group_name': group_name,
                        'group_type': 'OECD',
                        'smiles': smiles,
                        'formula': formula,
                        'inchi': inchi
                    })
                    test_id += 1
            except Exception as e:
                print(f"    Error generating molecule for group {group_id}: {e}")
                continue
    
    print("\nGenerating generic PFAS group test molecules...")
    for group_id, group_name, group_smiles, mode in GENERIC_PFAS_GROUPS:
        if group_id in IGNORE_GROUPS:
            continue
            
        print(f"  Group {group_id}: {group_name}")
        for i in range(n_per_group):
            try:
                # Generate random carbon chain
                n_carbons = 4 + (i % 8)
                chain = generate_random_carbon_chain(n_carbons, branching_probability=0.2)
                
                # Fluorinate with varying degrees
                if i % 3 == 0:
                    mol = fluorinate_mol(chain, perfluorinated=True, p=1.0)
                else:
                    mol = fluorinate_mol(chain, perfluorinated=False, p=0.6)
                
                # Append functional group
                attachment_point = get_attachment(mol, neighbours=None)
                if attachment_point is not None:
                    mol = append_functional_group(
                        mol, 
                        group_smiles, 
                        attachment_point, 
                        mode=mode
                    )
                
                if mol is not None:
                    # Remove explicit hydrogens before generating SMILES to avoid [H] notation
                    # This prevents issues with formula calculation
                    mol = Chem.RemoveHs(mol)
                    smiles = Chem.MolToSmiles(mol)
                    inchi = Chem.MolToInchi(mol)
                    formula = CalcMolFormula(mol)
                    
                    test_molecules.append({
                        'test_id': test_id,
                        'group_id': group_id,
                        'group_name': group_name,
                        'group_type': 'Generic',
                        'smiles': smiles,
                        'formula': formula,
                        'inchi': inchi
                    })
                    test_id += 1
                    
            except Exception as e:
                print(f"    Error generating molecule for group {group_id}: {e}")
                continue
    
    # Also add specificity test molecules
    print("\nAdding specificity test molecules...")
    spec_molecules = create_specificity_test_molecules()
    for _, row in spec_molecules.iterrows():
        for group_id in row['group_ids']:
            test_molecules.append({
                'test_id': test_id,
                'group_id': group_id,
                'group_name': f"Specificity test for group {group_id}",
                'group_type': 'Specificity',
                'smiles': row['smiles'],
                'formula': row['formula'],
                'inchi': row['inchi']
            })
            test_id += 1
    
    df = pd.DataFrame(test_molecules)
    print(f"\nGenerated {len(df)} test molecules across {df['group_id'].nunique()} groups")
    return df


def run_python_analysis(test_molecules_df):
    """
    Run Python parse_PFAS_groups on test molecules.
    
    Args:
        test_molecules_df: DataFrame with test molecules
        
    Returns:
        DataFrame with Python analysis results
    """
    print("\nRunning Python parse_PFAS_groups analysis...")
    results = []
    
    for idx, row in test_molecules_df.iterrows():
        try:
            mol = Chem.MolFromSmiles(row['smiles'])
            if mol is None:
                results.append({
                    'test_id': row['test_id'],
                    'python_success': False,
                    'python_groups': None,
                    'python_n_groups': 0,
                    'python_error': 'Invalid SMILES'
                })
                continue
            
            # Run parse_PFAS_groups with both mol and formula
            # Note: parse_PFAS_groups will call Chem.AddHs internally, so don't do it here
            # IMPORTANT: Recalculate formula from SMILES to avoid issues with explicit [H] atoms
            formula = CalcMolFormula(mol)
            detected_groups = parse_PFAS_groups(mol, formula)
            
            # parse_PFAS_groups returns list of tuples: (PFASGroup_object, n, n_CFchain, chains)
            if detected_groups:
                group_ids = [int(g[0].id) for g in detected_groups if len(g[3])==0 or any([k['SMARTS'] in ['Perfluoroalkyl','Polyfluoroalkyl'] for k in g[3]])]
            else:
                group_ids = []
            
            results.append({
                'test_id': row['test_id'],
                'python_success': True,
                'python_groups': json.dumps(sorted(group_ids)),
                'python_n_groups': len(group_ids),
                'python_error': None
            })
            
        except Exception as e:
            results.append({
                'test_id': row['test_id'],
                'python_success': False,
                'python_groups': None,
                'python_n_groups': 0,
                'python_error': str(e)
            })
    
    results_df = pd.DataFrame(results)
    return results_df


def export_for_html(test_molecules_df, output_csv='test_molecules_for_html.csv'):
    """
    Export test molecules in a format ready for HTML analyzer.
    
    Args:
        test_molecules_df: DataFrame with test molecules
        output_csv: Output CSV filename
    """
    # Create a simple CSV with SMILES and formula for HTML import
    html_input = test_molecules_df[['test_id', 'smiles', 'formula']].copy()
    html_input.columns = ['TestID', 'SMILES', 'Formula']
    
    html_input.to_csv(output_csv, index=False)
    print(f"\nExported {len(html_input)} molecules to {output_csv}")
    print("\nNEXT STEPS:")
    print("1. Open pfas_analyzer.html in your browser")
    print("2. Upload the file: test_molecules_for_html.csv")
    print("3. Check the 'Use molecular formula from file column' checkbox")
    print("4. Enter 'Formula' in the formula column name textbox")
    print("5. Click 'Analyze All Molecules'")
    print("6. After analysis completes, click 'Export to CSV'")
    print("7. Save the exported file as 'html_results.csv'")
    print("8. Run: python test_html_vs_python.py --compare html_results.csv")


def compare_results(python_baseline_file, html_results_file, output_report='comparison_report.json'):
    """
    Compare HTML results with Python baseline.
    
    Args:
        python_baseline_file: CSV file with Python results
        html_results_file: CSV file exported from HTML analyzer
        output_report: JSON file for detailed comparison report
        
    Returns:
        dict: Comparison statistics
    """
    print(f"\nComparing HTML results with Python baseline...")
    
    # Load Python baseline
    python_df = pd.read_csv(python_baseline_file)
    print(f"Loaded {len(python_df)} Python baseline results")
    
    # Load HTML results
    html_df = pd.read_csv(html_results_file)
    print(f"Loaded {len(html_df)} HTML results")
    
    # Check column names in HTML export
    print(f"HTML columns: {list(html_df.columns)}")
    
    # Parse HTML results to extract group IDs
    # HTML export format has columns: Input, Formula, PFAS_Groups_Count, PFAS_Groups
    # We'll match by SMILES (Input column) to Python baseline
    
    # First, load group name to ID mapping (including aliases)
    data_folder = os.path.join(os.path.dirname(os.path.dirname(__file__)), 'PFASgroups', 'data')
    with open(os.path.join(data_folder, 'PFAS_groups_smarts.json'), 'r') as f:
        group_data = json.load(f)
    name_to_id = {}
    id_to_name = {}
    for g in group_data:
        gid = int(g['id'])
        name_to_id[g['name']] = gid
        id_to_name[gid] = g['name']
        # Also map alias to ID
        if 'alias' in g and g['alias']:
            name_to_id[g['alias']] = gid
    
    comparison_results = []
    
    for idx, html_row in html_df.iterrows():
        try:
            # Match by SMILES (Input column in HTML export)
            smiles = html_row.get('Input', '')
            if not smiles:
                print(f"Warning: No SMILES in HTML row {idx}")
                continue
            
            # Get Python baseline for this SMILES
            python_row = python_df[python_df['smiles'] == smiles]
            if len(python_row) == 0:
                print(f"Warning: No Python baseline for SMILES: {smiles[:50]}...")
                continue
            python_row = python_row.iloc[0]
            test_id = python_row['test_id']
            
            # Extract HTML detected groups
            # HTML export has "PFAS_Groups" column with semicolon-separated group names
            html_groups_str = html_row.get('PFAS_Groups', '')
            html_groups = set()
            if isinstance(html_groups_str, str) and html_groups_str:
                # Parse group names from HTML output
                for group_name in html_groups_str.split(';'):
                    group_name = group_name.strip()
                    if group_name in name_to_id:
                        html_groups.add(name_to_id[group_name])
                    else:
                        print(f"Warning: Unknown group name '{group_name}' in HTML results")
            
            # Extract Python detected groups
            python_groups = set()
            if python_row['python_success'] and python_row['python_groups']:
                try:
                    python_groups = set(json.loads(python_row['python_groups']))
                except:
                    pass
            
            # Compare
            match = html_groups == python_groups
            
            # Get group names for better readability
            python_group_names = [f"{gid}:{id_to_name.get(gid, 'Unknown')}" for gid in sorted(python_groups)]
            html_group_names = [f"{gid}:{id_to_name.get(gid, 'Unknown')}" for gid in sorted(html_groups)]
            missing_names = [f"{gid}:{id_to_name.get(gid, 'Unknown')}" for gid in sorted(python_groups - html_groups)]
            extra_names = [f"{gid}:{id_to_name.get(gid, 'Unknown')}" for gid in sorted(html_groups - python_groups)]
            
            comparison_results.append({
                'test_id': test_id,
                'expected_group': python_row['group_id'],
                'expected_group_name': python_row['group_name'],
                'python_groups': sorted(list(python_groups)),
                'python_group_names': '; '.join(python_group_names),
                'html_groups': sorted(list(html_groups)),
                'html_group_names': '; '.join(html_group_names),
                'exact_match': match,
                'python_success': python_row['python_success'],
                'missing_in_html': sorted(list(python_groups - html_groups)),
                'missing_in_html_names': '; '.join(missing_names),
                'extra_in_html': sorted(list(html_groups - python_groups)),
                'extra_in_html_names': '; '.join(extra_names),
                'smiles': python_row['smiles'],
                'formula': python_row['formula']
            })
            
        except Exception as e:
            print(f"Error processing row {idx}: {e}")
            continue
    
    # Create comparison DataFrame
    comparison_df = pd.DataFrame(comparison_results)
    
    # Calculate statistics
    n_total = len(comparison_df)
    n_exact_match = len(comparison_df[comparison_df['exact_match']])
    n_python_success = len(comparison_df[comparison_df['python_success']])
    
    match_rate = n_exact_match / n_total if n_total > 0 else 0
    
    stats = {
        'timestamp': datetime.now().isoformat(),
        'total_comparisons': n_total,
        'exact_matches': n_exact_match,
        'match_rate': round(match_rate, 4),
        'python_success_rate': round(n_python_success / n_total, 4) if n_total > 0 else 0,
        'mismatches': []
    }
    
    # Collect mismatches
    mismatches = comparison_df[comparison_df['exact_match'] == False]
    
    # Save ALL mismatches to a separate detailed file
    if len(mismatches) > 0:
        mismatches_detailed = mismatches[[
            'test_id', 'expected_group', 'expected_group_name',
            'python_groups', 'python_group_names',
            'html_groups', 'html_group_names',
            'missing_in_html', 'missing_in_html_names',
            'extra_in_html', 'extra_in_html_names',
            'smiles', 'formula'
        ]].copy()
        mismatches_detailed.to_csv('mismatches_detailed.csv', index=False)
        print(f"\nDetailed mismatches saved to mismatches_detailed.csv ({len(mismatches)} rows)")
    
    for _, row in mismatches.head(20).iterrows():  # Limit to first 20 for JSON report
        stats['mismatches'].append({
            'test_id': int(row['test_id']),
            'expected_group': int(row['expected_group']),
            'expected_group_name': row['expected_group_name'],
            'python_groups': row['python_groups'],
            'python_group_names': row['python_group_names'],
            'html_groups': row['html_groups'],
            'html_group_names': row['html_group_names'],
            'missing_in_html': row['missing_in_html'],
            'missing_in_html_names': row['missing_in_html_names'],
            'extra_in_html': row['extra_in_html'],
            'extra_in_html_names': row['extra_in_html_names'],
            'smiles': row['smiles']
        })
    
    # Save detailed comparison
    comparison_df.to_csv('detailed_comparison.csv', index=False)
    print(f"\nDetailed comparison saved to detailed_comparison.csv")
    
    # Save summary report
    with open(output_report, 'w') as f:
        json.dump(stats, f, indent=2)
    print(f"Summary report saved to {output_report}")
    
    # Print summary
    print("\n" + "="*60)
    print("HTML vs PYTHON COMPARISON SUMMARY")
    print("="*60)
    print(f"Total comparisons: {n_total}")
    print(f"Exact matches: {n_exact_match} ({match_rate:.1%})")
    print(f"Mismatches: {n_total - n_exact_match}")
    print("="*60)
    
    if len(mismatches) > 0:
        print("\nFirst 10 mismatches:")
        for _, row in mismatches.head(10).iterrows():
            print(f"\nTest ID {row['test_id']}: {row['expected_group_name']}")
            print(f"  Python detected: {row['python_group_names']}")
            print(f"  HTML detected: {row['html_group_names']}")
            if row['missing_in_html']:
                print(f"  Missing in HTML: {row['missing_in_html_names']}")
            if row['extra_in_html']:
                print(f"  Extra in HTML: {row['extra_in_html_names']}")
            print(f"  SMILES: {row['smiles'][:80]}...")
        
        print(f"\n\nSee 'mismatches_detailed.csv' for all {len(mismatches)} mismatches with full details.")
    
    return stats


def main():
    parser = argparse.ArgumentParser(
        description='Test HTML JavaScript implementation vs Python implementation'
    )
    parser.add_argument(
        '--generate',
        action='store_true',
        help='Generate test molecules and run Python analysis'
    )
    parser.add_argument(
        '--compare',
        type=str,
        metavar='HTML_RESULTS',
        help='Compare HTML results CSV with Python baseline'
    )
    parser.add_argument(
        '--n-per-group',
        type=int,
        default=5,
        help='Number of molecules to generate per group (default: 5)'
    )
    
    args = parser.parse_args()
    
    if args.generate:
        # Generate test molecules
        test_molecules = generate_test_molecules_from_groups(n_per_group=args.n_per_group)
        
        # Run Python analysis
        python_results = run_python_analysis(test_molecules)
        
        # Merge results
        baseline = test_molecules.merge(python_results, on='test_id')
        baseline.to_csv('python_baseline.csv', index=False)
        print(f"\nPython baseline saved to python_baseline.csv")
        
        # Export for HTML
        export_for_html(test_molecules, 'test_molecules_for_html.csv')
        
    elif args.compare:
        # Compare HTML results with Python baseline
        if not os.path.exists('python_baseline.csv'):
            print("ERROR: python_baseline.csv not found!")
            print("Run with --generate first to create baseline")
            sys.exit(1)
        
        if not os.path.exists(args.compare):
            print(f"ERROR: HTML results file not found: {args.compare}")
            sys.exit(1)
        
        compare_results('python_baseline.csv', args.compare)
        
    else:
        parser.print_help()


if __name__ == '__main__':
    main()
