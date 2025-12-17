#!/usr/bin/env python3
"""
Generate molecular structure drawings for misclassified molecules in PFASGroups testing.
"""

import pandas as pd
import ast
import json
import sys
import os
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem.Draw import rdMolDraw2D
import base64
from io import BytesIO

# Load group names mapping
try:
    with open('/home/luc/git/PFASGroups/PFASgroups/data/PFAS_groups_smarts.json', 'r') as f:
        groups_data = json.load(f)
    group_names = {int(g['id']): g['name'] for g in groups_data}
except:
    group_names = {}

def draw_molecule_svg(smiles, width=400, height=300):
    """Draw molecule as SVG string."""
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None
        
        drawer = rdMolDraw2D.MolDraw2DSVG(width, height)
        drawer.DrawMolecule(mol)
        drawer.FinishDrawing()
        svg = drawer.GetDrawingText()
        return svg
    except:
        return None

def generate_misclassification_table():
    """Generate comprehensive misclassification table with molecular structures."""
    
    misclassified_data = []
    
    # Process OECD misclassified molecules
    print("Processing OECD misclassified molecules...")
    oecd_df = pd.read_csv('oecd_test_results.csv')
    oecd_df['all_matches'] = oecd_df['all_matches'].apply(ast.literal_eval)
    oecd_failed = oecd_df[oecd_df['detected'] == False]
    
    for i, row in oecd_failed.iterrows():
        expected_names = [group_names.get(row['group_ids'], f"Unknown_{row['group_ids']}")]
        detected_names = [group_names.get(g, f"Unknown_{g}") for g in row['all_matches']]
        
        svg = draw_molecule_svg(row['smiles'])
        
        misclassified_data.append({
            'Test_Suite': 'OECD',
            'Expected_Groups': f"{row['group_ids']} ({expected_names[0]})",
            'Detected_Groups': f"{row['all_matches']} ({', '.join(detected_names)})",
            'SMILES': row['smiles'],
            'Chain_Length': row['chain_length'],
            'Pathtype': row['pathtype'],
            'Error_Type': 'False Negative',
            'SVG': svg if svg else 'Structure drawing failed',
            'Origin': row['group_name']
        })
    
    # Process Generic misclassified molecules
    print("Processing Generic misclassified molecules...")
    generic_df = pd.read_csv('generic_test_results.csv')
    generic_df['all_matches'] = generic_df['all_matches'].apply(ast.literal_eval)
    generic_failed = generic_df[generic_df['detected'] == False]
    
    for i, row in generic_failed.iterrows():
        expected_names = [group_names.get(row['group_ids'], f"Unknown_{row['group_ids']}")]
        detected_names = [group_names.get(g, f"Unknown_{g}") for g in row['all_matches']]
        
        svg = draw_molecule_svg(row['smiles'])
        
        misclassified_data.append({
            'Test_Suite': 'Generic',
            'Expected_Groups': f"{row['group_ids']} ({expected_names[0]})",
            'Detected_Groups': f"{row['all_matches']} ({', '.join(detected_names)})",
            'SMILES': row['smiles'],
            'Chain_Length': row['chain_length'],
            'Pathtype': row['pathtype'],
            'Error_Type': 'False Negative',
            'SVG': svg if svg else 'Structure drawing failed',
            'Origin': row['group_name']
        })
    
    # Process Specificity detection failures (first 20 for brevity)
    print("Processing Specificity detection failures...")
    spec_df = pd.read_csv('specificity_test_results.csv')
    spec_df['group_ids'] = spec_df['group_ids'].apply(ast.literal_eval)
    spec_df['detected_groups'] = spec_df['detected_groups'].apply(ast.literal_eval)
    spec_failed = spec_df[spec_df['expected_group_detected'] == False].head(20)
    
    for i, row in spec_failed.iterrows():
        expected_names = [group_names.get(g, f'Unknown_{g}') for g in row['group_ids']]
        detected_names = [group_names.get(g, f'Unknown_{g}') for g in row['detected_groups']]
        
        svg = draw_molecule_svg(row['smiles']) if pd.notna(row['smiles']) else None
        
        misclassified_data.append({
            'Test_Suite': 'Specificity',
            'Expected_Groups': f"{row['group_ids']} ({', '.join(expected_names)})",
            'Detected_Groups': f"{row['detected_groups']} ({', '.join(detected_names)})",
            'SMILES': row['smiles'],
            'Chain_Length': 'Variable',
            'Pathtype': 'Mixed',
            'Error_Type': 'False Negative',
            'SVG': svg if svg else 'Structure drawing failed',
            'Origin': row['origin']
        })
    
    # Convert to DataFrame and save
    df = pd.DataFrame(misclassified_data)
    
    # Save detailed results
    df.to_csv('misclassified_molecules_detailed.csv', index=False)
    
    # Generate HTML table with structures
    html_content = """
<!DOCTYPE html>
<html>
<head>
    <title>PFASGroups Misclassified Molecules</title>
    <style>
        body { font-family: Arial, sans-serif; margin: 20px; }
        table { border-collapse: collapse; width: 100%; margin-top: 20px; }
        th, td { border: 1px solid #ddd; padding: 12px; text-align: left; }
        th { background-color: #f2f2f2; font-weight: bold; }
        .structure { max-width: 300px; }
        .smiles { font-family: monospace; font-size: 12px; max-width: 400px; word-break: break-all; }
        .error-type { font-weight: bold; }
        .false-negative { color: #d32f2f; }
        .test-suite { font-weight: bold; padding: 4px 8px; border-radius: 4px; color: white; }
        .oecd { background-color: #1976d2; }
        .generic { background-color: #388e3c; }
        .specificity { background-color: #f57c00; }
    </style>
</head>
<body>
    <h1>PFASGroups Misclassified Molecules Analysis</h1>
    <p><strong>Generated:</strong> """ + pd.Timestamp.now().strftime('%Y-%m-%d %H:%M:%S') + """</p>
    <p><strong>Total Misclassified:</strong> """ + str(len(df)) + """ molecules</p>
    
    <h2>Summary by Test Suite</h2>
    <ul>
        <li><strong>OECD Tests:</strong> """ + str(len(df[df['Test_Suite'] == 'OECD'])) + """ misclassified</li>
        <li><strong>Generic Tests:</strong> """ + str(len(df[df['Test_Suite'] == 'Generic'])) + """ misclassified</li>
        <li><strong>Specificity Tests:</strong> """ + str(len(df[df['Test_Suite'] == 'Specificity'])) + """ misclassified (top 20 shown)</li>
    </ul>
    
    <table>
        <thead>
            <tr>
                <th>Test Suite</th>
                <th>Structure</th>
                <th>Expected Groups</th>
                <th>Detected Groups</th>
                <th>SMILES</th>
                <th>Error Type</th>
                <th>Origin</th>
            </tr>
        </thead>
        <tbody>
"""
    
    for i, row in df.iterrows():
        test_suite_class = row['Test_Suite'].lower()
        
        # Handle SVG structure
        structure_content = row['SVG'] if row['SVG'] and row['SVG'] != 'Structure drawing failed' else 'No structure available'
        if structure_content != 'No structure available':
            # Clean up SVG for embedding
            structure_content = structure_content.replace('<?xml version="1.0" encoding="UTF-8"?>', '')
        
        html_content += f"""
            <tr>
                <td><span class="test-suite {test_suite_class}">{row['Test_Suite']}</span></td>
                <td class="structure">{structure_content}</td>
                <td><strong>{row['Expected_Groups']}</strong></td>
                <td>{row['Detected_Groups']}</td>
                <td class="smiles">{row['SMILES']}</td>
                <td class="error-type false-negative">{row['Error_Type']}</td>
                <td>{row['Origin']}</td>
            </tr>
        """
    
    html_content += """
        </tbody>
    </table>
    
    <h2>Notes</h2>
    <ul>
        <li><strong>False Negative:</strong> Expected group not detected by PFASGroups</li>
        <li><strong>OECD Tests:</strong> Perfluoroalkyl disulfonic acids (Group 8) consistently misclassified as polyfluoroalkyl sulfonic acid</li>
        <li><strong>Generic Tests:</strong> One alkane group missed due to complex polyfluoroalkyl structure</li>
        <li><strong>Specificity Tests:</strong> Many cases where alkane groups (Group 48) are not detected due to complex fluorination patterns</li>
    </ul>
    
</body>
</html>
"""
    
    # Save HTML report
    with open('misclassified_molecules_report.html', 'w') as f:
        f.write(html_content)
    
    print(f"Generated misclassification report:")
    print(f"- CSV: misclassified_molecules_detailed.csv ({len(df)} molecules)")
    print(f"- HTML: misclassified_molecules_report.html")
    
    # Generate markdown table for summary
    markdown_table = "| Test Suite | Expected Groups | Detected Groups | SMILES | Error Type | Origin |\n"
    markdown_table += "|------------|----------------|----------------|---------|------------|--------|\n"
    
    for i, row in df.iterrows():
        # Truncate SMILES for readability
        smiles_short = row['SMILES'][:50] + "..." if len(row['SMILES']) > 50 else row['SMILES']
        markdown_table += f"| {row['Test_Suite']} | {row['Expected_Groups']} | {row['Detected_Groups']} | `{smiles_short}` | {row['Error_Type']} | {row['Origin']} |\n"
    
    # Save markdown table
    with open('misclassified_molecules_table.md', 'w') as f:
        f.write(markdown_table)
    
    print(f"- Markdown table: misclassified_molecules_table.md")
    
    return df

if __name__ == "__main__":
    print("Generating misclassification analysis...")
    df = generate_misclassification_table()
    print(f"\nAnalysis complete! {len(df)} misclassified molecules processed.")