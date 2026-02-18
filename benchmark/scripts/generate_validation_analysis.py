#!/usr/bin/env python3
"""
Generate detailed benchmark validation analysis with molecular structures and Sankey diagram.
"""

import pandas as pd
import ast
import json
import sys
import os
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem.Draw import rdMolDraw2D
from rdkit.Chem import rdMolDescriptors
import base64
from io import BytesIO
import plotly.graph_objects as go
from collections import defaultdict

def draw_molecule_svg(smiles, width=300, height=200):
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

def analyze_benchmark_correspondences():
    """Analyze OECD benchmark correspondences."""
    
    # Load OECD benchmark data
    benchmark_file = '/home/luc/git/HalogenGroups/benchmark/data/pfas_oecd_benchmark_20251217_142445.json'
    with open(benchmark_file, 'r') as f:
        oecd_data = json.load(f)
    
    # Load group names
    try:
        with open('/home/luc/git/HalogenGroups/HalogenGroups/data/PFAS_groups_smarts.json', 'r') as f:
            groups_data = json.load(f)
        group_names = {int(g['id']): g['name'] for g in groups_data}
    except:
        group_names = {}
    
    # Analyze correspondences for Sankey diagram
    sankey_data = defaultdict(int)
    disagreement_examples = []
    total_molecules = len(oecd_data)
    success_count = 0
    
    for entry in oecd_data:
        if entry['HalogenGroups_result']['success'] and entry['atlas_result']['success']:
            success_count += 1
            atlas_class = f"{entry['atlas_result']['first_class']} / {entry['atlas_result']['second_class']}"
            pfas_groups = entry['HalogenGroups_result']['detected_groups']
            pfas_names = [group_names.get(g, f'Group_{g}') for g in pfas_groups]
            pfas_class = ', '.join(pfas_names) if pfas_names else 'No groups detected'
            
            sankey_data[(atlas_class, pfas_class)] += 1
            
            # Track significant disagreements for examples
            if len(pfas_groups) == 0:
                disagreement_examples.append({
                    'type': 'No_HalogenGroups_Detected',
                    'smiles': entry['molecule_data']['smiles'],
                    'atlas': atlas_class,
                    'HalogenGroups': pfas_class,
                    'HalogenGroups_ids': pfas_groups
                })
    
    return sankey_data, disagreement_examples, success_count, total_molecules

def create_sankey_diagram(sankey_data):
    """Create a Sankey diagram for Atlas vs HalogenGroups correspondence."""
    
    # Prepare data for Sankey
    atlas_classes = list(set(pair[0] for pair in sankey_data.keys()))
    pfas_classes = list(set(pair[1] for pair in sankey_data.keys()))
    
    # Create node labels and indices
    all_nodes = atlas_classes + pfas_classes
    node_indices = {node: i for i, node in enumerate(all_nodes)}
    
    # Create links
    sources = []
    targets = []
    values = []
    
    for (atlas_class, pfas_class), count in sankey_data.items():
        sources.append(node_indices[atlas_class])
        targets.append(node_indices[pfas_class])
        values.append(count)
    
    # Create Sankey diagram
    fig = go.Figure(data=[go.Sankey(
        node = dict(
            pad = 15,
            thickness = 20,
            line = dict(color = "black", width = 0.5),
            label = all_nodes,
            color = ["lightblue" if i < len(atlas_classes) else "lightgreen" for i in range(len(all_nodes))]
        ),
        link = dict(
            source = sources,
            target = targets,
            value = values
        )
    )])
    
    fig.update_layout(
        title_text="PFAS Classification Correspondence: PFAS-Atlas vs HalogenGroups",
        font_size=10,
        width=1200,
        height=800
    )
    
    return fig

def analyze_disulfonic_acid_failures():
    """Analyze the perfluoroalkyl disulfonic acid classification failures."""
    
    # Load test results
    oecd_df = pd.read_csv('/home/luc/git/HalogenGroups/HalogenGroups/tests/results/oecd_test_results.csv')
    oecd_df['all_matches'] = oecd_df['all_matches'].apply(ast.literal_eval)
    
    # Load group names
    try:
        with open('/home/luc/git/HalogenGroups/HalogenGroups/data/PFAS_groups_smarts.json', 'r') as f:
            groups_data = json.load(f)
        group_names = {int(g['id']): g['name'] for g in groups_data}
    except:
        group_names = {}
    
    # Find Group 8 failures
    group_8_failed = oecd_df[(oecd_df['group_ids'] == 8) & (oecd_df['detected'] == False)]
    group_8_success = oecd_df[(oecd_df['group_ids'] == 8) & (oecd_df['detected'] == True)]
    
    # Analyze failed molecules
    failed_analysis = []
    for i, row in group_8_failed.iterrows():
        smiles = row['smiles']
        
        # Analyze structure
        try:
            mol = Chem.MolFromSmiles(smiles)
            if mol:
                formula = rdMolDescriptors.CalcMolFormula(mol)
                # Count sulfonic acid groups
                sulfonic_pattern = Chem.MolFromSmarts('S(=O)(=O)O')
                sulfonic_matches = mol.GetSubstructMatches(sulfonic_pattern) if sulfonic_pattern else []
                sulfonic_count = len(sulfonic_matches)
                
                # Generate structure drawing
                svg = draw_molecule_svg(smiles)
                
                detected_names = [group_names.get(g, f'Group_{g}') for g in row['all_matches']]
                
                failed_analysis.append({
                    'smiles': smiles,
                    'formula': formula,
                    'sulfonic_count': sulfonic_count,
                    'detected_groups': row['all_matches'],
                    'detected_names': detected_names,
                    'svg': svg,
                    'chain_length': row['chain_length'],
                    'pathtype': row['pathtype']
                })
        except Exception as e:
            failed_analysis.append({
                'smiles': smiles,
                'formula': 'Error',
                'sulfonic_count': 'Error',
                'detected_groups': row['all_matches'],
                'detected_names': 'Error',
                'svg': None,
                'error': str(e)
            })
    
    # Analyze successful molecules for comparison
    success_analysis = []
    for i, row in group_8_success.head(3).iterrows():  # Just first 3 for comparison
        smiles = row['smiles']
        
        try:
            mol = Chem.MolFromSmiles(smiles)
            if mol:
                formula = rdMolDescriptors.CalcMolFormula(mol)
                sulfonic_pattern = Chem.MolFromSmarts('S(=O)(=O)O')
                sulfonic_matches = mol.GetSubstructMatches(sulfonic_pattern) if sulfonic_pattern else []
                sulfonic_count = len(sulfonic_matches)
                
                svg = draw_molecule_svg(smiles)
                detected_names = [group_names.get(g, f'Group_{g}') for g in row['all_matches']]
                
                success_analysis.append({
                    'smiles': smiles,
                    'formula': formula,
                    'sulfonic_count': sulfonic_count,
                    'detected_groups': row['all_matches'],
                    'detected_names': detected_names,
                    'svg': svg,
                    'chain_length': row['chain_length'],
                    'pathtype': row['pathtype']
                })
        except:
            pass
    
    return failed_analysis, success_analysis, len(group_8_failed), len(group_8_success)

def generate_validation_report():
    """Generate comprehensive validation report."""
    
    print("Generating benchmark validation analysis...")
    
    # Analyze benchmark correspondences
    sankey_data, disagreements, success_count, total_count = analyze_benchmark_correspondences()
    
    # Create Sankey diagram
    sankey_fig = create_sankey_diagram(sankey_data)
    sankey_html = sankey_fig.to_html(include_plotlyjs='cdn')
    
    # Analyze disulfonic acid failures
    failed_disulfonic, success_disulfonic, n_failed, n_success = analyze_disulfonic_acid_failures()
    
    # Convert sankey_data tuples to strings for JSON serialization
    sankey_data_serializable = {f"{k[0]} -> {k[1]}": v for k, v in sankey_data.items()}
    
    # Generate results
    results = {
        'benchmark_summary': {
            'total_molecules': total_count,
            'successful_processing': success_count,
            'success_rate': success_count / total_count,
            'unique_correspondences': len(sankey_data)
        },
        'disulfonic_analysis': {
            'failed_count': n_failed,
            'success_count': n_success,
            'failed_molecules': failed_disulfonic,
            'success_molecules': success_disulfonic
        },
        'sankey_data': sankey_data_serializable,
        'sankey_html': sankey_html,
        'disagreement_examples': disagreements[:10]  # First 10 examples
    }
    
    # Save results
    with open('benchmark_validation_analysis.json', 'w') as f:
        json.dump({k: v for k, v in results.items() if k != 'sankey_html'}, f, indent=2, default=str)
    
    # Save Sankey HTML separately
    with open('atlas_HalogenGroups_sankey.html', 'w') as f:
        f.write(sankey_html)
    
    print(f"Generated benchmark validation analysis:")
    print(f"- Total molecules: {total_count}")
    print(f"- Success rate: {success_count/total_count:.1%}")
    print(f"- Disulfonic acid failures: {n_failed}/{n_failed+n_success}")
    print(f"- Saved to: benchmark_validation_analysis.json")
    print(f"- Sankey diagram: atlas_HalogenGroups_sankey.html")
    
    return results

if __name__ == "__main__":
    results = generate_validation_report()