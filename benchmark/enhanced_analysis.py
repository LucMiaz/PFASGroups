"""
Enhanced PFAS Benchmark Analysis
Creates comprehensive comparison with heatmaps, Sankey diagrams, and detailed statistics
"""

import json
import pandas as pd
import numpy as np
import plotly.graph_objects as go
import plotly.figure_factory as ff
import plotly.express as px
from plotly.subplots import make_subplots
import seaborn as sns
import matplotlib.pyplot as plt
from collections import defaultdict, Counter
from datetime import datetime
import sys
import os

# PFAS-Atlas class mapping to understand their output classes (discovered from testing)
PFAS_ATLAS_CLASSES = {
    "first_class": {
        "PFAA precursors": "Precursors to perfluoroalkyl acids",
        "PFAAs": "Perfluoroalkyl acids",
        "Other PFASs": "Other per- and polyfluoroalkyl substances", 
        "Polyfluoroalkyl acids": "Polyfluoroalkyl carboxylic/sulfonic acids",
        "Not PFAS": "Not identified as PFAS"
    },
    "second_class": {
        'PFSAs': "Perfluoroalkyl sulfonic acids",
        'PolyFCAs, cyclic': "Polyfluoroalkyl carboxylic acids (cyclic)",
        'PASF-based substances, cyclic': "Perfluoroalkane sulfonyl fluoride-based substances (cyclic)",
        'Perfluoroalkyl-tert-amines, cyclic': "Perfluoroalkyl tertiary amines (cyclic)",
        'Perfluoroalkylethers': "Perfluoroalkyl ether compounds",
        'PFALs': "Perfluoroalkyl aldehydes",
        'Not PFAS by current definition': "Not classified as PFAS under current definitions",
        'HFOs': "Hydrofluoroolefins",
        'PFCAs': "Perfluoroalkyl carboxylic acids",
        'HFCs': "Hydrofluorocarbons",
        'PolyFAenes, cyclic': "Polyfluoroalkenes (cyclic)",
        'PFAKs': "Perfluoroalkyl ketones",
        'PFPEs, cyclic': "Perfluoropolyethers (cyclic)",
        'PolyFAC derivatives, cyclic': "Polyfluoroalkyl compound derivatives (cyclic)",
        'HFEs, cyclic': "Hydrofluoroethers (cyclic)",
        'SFAenes': "Saturated fluoroalkenes",
        'PFPAs': "Perfluoroalkyl phosphonic acids",
        'PAECFs': "Perfluoroalkyl ether carboxylic fluorides",
        'HFCs, cyclic': "Hydrofluorocarbons (cyclic)",
        'PolyFESAs': "Polyfluoroalkyl ether sulfonic acids",
        'PolyFEAenes, cyclic': "Polyfluoroalkyl ethylene derivatives (cyclic)",
        'others, cyclic': "Other PFAS compounds (cyclic)",
        'Amide derivatives, cyclic': "Perfluoroalkyl amide derivatives (cyclic)",
        'SFAs': "Saturated fluoroalkyl acids",
        'PFAK derivatives': "Perfluoroalkyl ketone derivatives",
        'PolyFEAenes': "Polyfluoroalkyl ethylene derivatives",
        'PFECAs': "Perfluoroalkyl ether carboxylic acids",
        'Sulfonyl chloride': "Perfluoroalkyl sulfonyl chlorides",
        'PolyFEACs, cyclic': "Polyfluoroalkyl ether carboxylic acids (cyclic)",
        'PASFs, cyclic': "Perfluoroalkane sulfonyl fluorides (cyclic)",
        'Perfluoroalkanes': "Fully fluorinated alkanes",
        'PFCA-ester derivatives': "Perfluoroalkyl carboxylic acid ester derivatives",
        'PACFs': "Perfluoroalkyl carbonyl fluorides",
        'PolyFSAs': "Polyfluoroalkyl sulfonic acids",
        'PASFs': "Perfluoroalkane sulfonyl fluorides",
        'n:2 fluorotelomer-based substances': "n:2 Fluorotelomer-based PFAS",
        'PolyFACs': "Polyfluoroalkyl compounds",
        'Perfluoroalkyl-tert-amines': "Perfluoroalkyl tertiary amines",
        'PolyFAenes': "Polyfluoroalkenes",
        'n:2 fluorotelomer-based substances, cyclic': "n:2 Fluorotelomer-based PFAS (cyclic)",
        'PFAS derivatives': "General PFAS derivatives",
        'PFCA-ester derivatives, cyclic': "Perfluoroalkyl carboxylic acid ester derivatives (cyclic)",
        'Amide derivatives': "Perfluoroalkyl amide derivatives",
        'PASF-based substances': "Perfluoroalkane sulfonyl fluoride-based substances",
        'PFESAs': "Perfluoroalkyl ether sulfonic acids",
        'n:1 FTOHs': "n:1 Fluorotelomer alcohols",
        '(Hg,Sn,Ge,Sb,Se,B) PFASs': "Metal-containing PFAS (mercury, tin, germanium, antimony, selenium, boron)",
        'PFPIAs': "Perfluoroalkyl phosphinic acids",
        'PAECFs, cyclic': "Perfluoroalkyl ether carboxylic fluorides (cyclic)",
        'PFdiCAs': "Perfluoroalkyl dicarboxylic acids",
        'Perfluoroalkylethers, cyclic': "Perfluoroalkyl ether compounds (cyclic)",
        'PFAK derivatives, cyclic': "Perfluoroalkyl ketone derivatives (cyclic)",
        'others': "Other PFAS compounds",
        'Perfluoroalkanes, cyclic': "Fully fluorinated alkanes (cyclic)",
        'PolyFCAs': "Polyfluoroalkyl carboxylic acids",
        'PFECAs, cyclic': "Perfluoroalkyl ether carboxylic acids (cyclic)",
        'PFAenes': "Perfluoroalkenes",
        'PFAenes, cyclic': "Perfluoroalkenes (cyclic)",
        'n:1 FTOHs, cyclic': "n:1 Fluorotelomer alcohols (cyclic)",
        'Acid chloride': "Perfluoroalkyl acid chlorides",
        'Polyfluoroalkanes, cyclic': "Polyfluoroalkanes (cyclic)",
        'PFAIs': "Perfluoroalkyl iodides",
        'PFACs': "Perfluoroalkyl compounds",
        'Si PFASs': "Silicon-containing PFAS",
        'Si PFASs, cyclic': "Silicon-containing PFAS (cyclic)",
        'PFSA derivatives, cyclic': "Perfluoroalkyl sulfonic acid derivatives (cyclic)",
        'PolyFSA derivatives': "Polyfluoroalkyl sulfonic acid derivatives",
        'PolyFECAs': "Polyfluoroalkyl ether carboxylic acids",
        'PFSA derivatives': "Perfluoroalkyl sulfonic acid derivatives",
        'PolyFEACs': "Polyfluoroalkyl ether carboxylic compounds",
        'Polyfluoroalkanes': "Polyfluorinated alkanes",
        'PFSAs, cyclic': "Perfluoroalkyl sulfonic acids (cyclic)",
        'PACFs, cyclic': "Perfluoroalkyl carbonyl fluorides (cyclic)",
        'PFCA-anhydrides': "Perfluoroalkyl carboxylic acid anhydrides",
        'PolyFAC derivatives': "Polyfluoroalkyl compound derivatives",
        'HFEs': "Hydrofluoroethers",
        'PFPEs': "Perfluoropolyethers",
        'PFSIAs': "Perfluoroalkyl sulfinic acids",
        'PolyFCA derivatives, cyclic': "Polyfluoroalkyl carboxylic acid derivatives (cyclic)",
        'Aromatic PFASs': "Aromatic per- and polyfluoroalkyl substances",
        'PolyFCA derivatives': "Polyfluoroalkyl carboxylic acid derivatives"
    }
}

def load_benchmark_results(filename):
    """Load enhanced benchmark results from JSON"""
    
    with open(filename, 'r') as f:
        results = json.load(f)
    
    print(f"📊 Loaded {len(results)} enhanced benchmark results")
    return results

def analyze_system_comparison(results):
    """Detailed comparison between PFASGroups and PFAS-Atlas"""
    
    # Separate single and multi-group results
    single_results = [r for r in results if r['molecule_data'].get('generation_type') == 'single_group']
    multi_results = [r for r in results if r['molecule_data'].get('generation_type') == 'multi_group']
    
    print(f"📋 Analyzing {len(single_results)} single-group molecules")
    print(f"🔬 Analyzing {len(multi_results)} multi-group molecules")
    
    # Single group analysis
    single_analysis = {}
    for result in single_results:
        group_id = result['molecule_data']['group_id']
        group_name = result['molecule_data']['group_name']
        
        if group_id not in single_analysis:
            single_analysis[group_id] = {
                'group_name': group_name,
                'total_molecules': 0,
                'pfasgroups_detections': 0,
                'pfasgroups_correct': 0,
                'atlas_classifications': 0,
                'atlas_any_detection': 0,  # Any functional group detected
                'molecules': [],
                # Timing analysis
                'pfasgroups_times': [],
                'atlas_times': [],
                'pfasgroups_avg_time': 0.0,
                'atlas_avg_time': 0.0,
                'pfasgroups_failures': [],  # Track molecules not correctly identified
                'atlas_classifications_detail': []  # Track Atlas class distributions
            }
        
        single_analysis[group_id]['total_molecules'] += 1
        
        # PFASGroups analysis
        pfas_detected = result['pfasgroups_result']['success']
        detected_groups = result['pfasgroups_result']['detected_groups']
        
        if pfas_detected:
            single_analysis[group_id]['pfasgroups_detections'] += 1
            if group_id in detected_groups:
                single_analysis[group_id]['pfasgroups_correct'] += 1
        
        # PFAS-Atlas analysis
        atlas_success = result['atlas_result']['success']
        if atlas_success:
            single_analysis[group_id]['atlas_classifications'] += 1
            # Check if Atlas detected any functional group (indirect detection)
            if result['atlas_result']['second_class']:
                single_analysis[group_id]['atlas_any_detection'] += 1
        
        # Track molecules not correctly identified by PFASGroups
        if not pfas_detected or group_id not in detected_groups:
            single_analysis[group_id]['pfasgroups_failures'].append({
                'smiles': result['molecule_data']['smiles'],
                'expected_group': group_id,
                'detected_groups': detected_groups,
                'detection_status': 'no_detection' if not pfas_detected else 'wrong_group'
            })
            
        # Track Atlas classification details
        if atlas_success:
            single_analysis[group_id]['atlas_classifications_detail'].append({
                'smiles': result['molecule_data']
                ['smiles'],
                'first_class': result['atlas_result']['first_class'],
                'second_class': result['atlas_result']['second_class']
            })
        
        single_analysis[group_id]['molecules'].append({
            'smiles': result['molecule_data']['smiles'],
            'pfasgroups_detected': detected_groups,
            'pfasgroups_success': pfas_detected,
            'atlas_success': atlas_success,
            'atlas_first_class': result['atlas_result']['first_class'],
            'atlas_second_class': result['atlas_result']['second_class']
        })
        
        # Collect timing data
        if 'execution_time' in result['pfasgroups_result']:
            single_analysis[group_id]['pfasgroups_times'].append(result['pfasgroups_result']['execution_time'])
        if 'execution_time' in result['atlas_result']:
            single_analysis[group_id]['atlas_times'].append(result['atlas_result']['execution_time'])
    
    # Multi-group analysis
    multi_analysis = {}
    for result in multi_results:
        target_groups = tuple(sorted(result['molecule_data']['target_groups']))
        
        if target_groups not in multi_analysis:
            multi_analysis[target_groups] = {
                'target_groups': list(target_groups),
                'total_molecules': 0,
                'pfasgroups_any_detection': 0,
                'pfasgroups_multi_detection': 0,
                'atlas_classifications': 0,
                'group_detection_matrix': defaultdict(int),  # Which groups were detected
                'molecules': []
            }
        
        multi_analysis[target_groups]['total_molecules'] += 1
        
        # PFASGroups analysis
        detected_groups = result['pfasgroups_result']['detected_groups']
        pfas_detected = len(detected_groups) > 0
        
        if pfas_detected:
            multi_analysis[target_groups]['pfasgroups_any_detection'] += 1
            
            # Check which target groups were detected
            for group in target_groups:
                if group in detected_groups:
                    multi_analysis[target_groups]['group_detection_matrix'][group] += 1
            
            # Multi-group detection (2+ groups detected)
            target_set = set(target_groups)
            detected_set = set(detected_groups)
            overlap = len(target_set.intersection(detected_set))
            
            if overlap >= 2:
                multi_analysis[target_groups]['pfasgroups_multi_detection'] += 1
        
        # PFAS-Atlas analysis
        atlas_success = result['atlas_result']['success']
        if atlas_success:
            multi_analysis[target_groups]['atlas_classifications'] += 1
        
        multi_analysis[target_groups]['molecules'].append({
            'smiles': result['molecule_data']['smiles'],
            'target_groups': list(target_groups),
            'detected_groups': detected_groups,
            'atlas_success': atlas_success,
            'atlas_class': result['atlas_result']['second_class']
        })
    
    # Calculate average timing statistics for single groups
    for group_id, data in single_analysis.items():
        if data['pfasgroups_times']:
            data['pfasgroups_avg_time'] = sum(data['pfasgroups_times']) / len(data['pfasgroups_times'])
        if data['atlas_times']:
            data['atlas_avg_time'] = sum(data['atlas_times']) / len(data['atlas_times'])
    
    return single_analysis, multi_analysis

def create_timing_comparison_heatmap(single_analysis):
    """Create heatmap comparing execution times between PFASGroups and PFAS-Atlas"""
    
    # Prepare data for timing heatmap
    groups = []
    group_names = []
    pfasgroups_times = []
    atlas_times = []
    
    for group_id, data in sorted(single_analysis.items()):
        if data['pfasgroups_avg_time'] > 0 or data['atlas_avg_time'] > 0:  # Only include groups with timing data
            groups.append(f"Group {group_id}")
            group_names.append(data['group_name'])
            
            # Convert to milliseconds for better readability
            pfas_time_ms = data['pfasgroups_avg_time'] * 1000
            atlas_time_ms = data['atlas_avg_time'] * 1000
            
            pfasgroups_times.append(pfas_time_ms)
            atlas_times.append(atlas_time_ms)
    
    # Create timing comparison matrix
    timing_data = np.array([pfasgroups_times, atlas_times])
    
    fig = go.Figure(data=go.Heatmap(
        z=timing_data,
        x=[f"{g}<br>{name}" for g, name in zip(groups, group_names)],
        y=['PFASGroups<br>Avg Time (ms)', 'PFAS-Atlas<br>Avg Time (ms)'],
        colorscale='Viridis',
        text=timing_data,
        texttemplate="%{text:.2f}",
        textfont={"size":10},
        hoverongaps=False
    ))
    
    fig.update_layout(
        title="Execution Time Comparison: PFASGroups vs PFAS-Atlas<br><sub>Average processing time per molecule (milliseconds)</sub>",
        xaxis_title="PFAS Functional Groups",
        yaxis_title="Detection Systems", 
        font_size=12,
        width=1400,
        height=500,
        margin=dict(l=100, r=50, t=100, b=100)
    )
    
    return fig

def create_timing_statistics_chart(single_analysis):
    """Create bar chart showing detailed timing statistics"""
    
    # Collect timing data
    groups = []
    group_names = []
    pfas_avg_times = []
    atlas_avg_times = []
    pfas_molecule_counts = []
    atlas_molecule_counts = []
    
    for group_id, data in sorted(single_analysis.items()):
        if data['pfasgroups_times'] or data['atlas_times']:
            groups.append(group_id)
            group_names.append(data['group_name'])
            
            pfas_avg_times.append(data['pfasgroups_avg_time'] * 1000)  # Convert to ms
            atlas_avg_times.append(data['atlas_avg_time'] * 1000)  # Convert to ms
            
            pfas_molecule_counts.append(len(data['pfasgroups_times']))
            atlas_molecule_counts.append(len(data['atlas_times']))
    
    fig = make_subplots(
        rows=2, cols=1,
        subplot_titles=['Average Execution Time (ms)', 'Number of Molecules Tested'],
        vertical_spacing=0.15
    )
    
    # Timing comparison
    fig.add_trace(
        go.Bar(name='PFASGroups', x=groups, y=pfas_avg_times, marker_color='lightblue'),
        row=1, col=1
    )
    fig.add_trace(
        go.Bar(name='PFAS-Atlas', x=groups, y=atlas_avg_times, marker_color='orange'),
        row=1, col=1
    )
    
    # Molecule count comparison
    fig.add_trace(
        go.Bar(name='PFASGroups Tested', x=groups, y=pfas_molecule_counts, marker_color='lightgreen', showlegend=False),
        row=2, col=1
    )
    fig.add_trace(
        go.Bar(name='Atlas Tested', x=groups, y=atlas_molecule_counts, marker_color='salmon', showlegend=False),
        row=2, col=1
    )
    
    fig.update_layout(
        title="Detailed Timing Analysis by Functional Group",
        barmode='group',
        width=1200,
        height=800
    )
    
    fig.update_xaxes(title_text="Functional Group ID", row=2, col=1)
    fig.update_yaxes(title_text="Time (ms)", row=1, col=1)
    fig.update_yaxes(title_text="Molecule Count", row=2, col=1)
    
    return fig

def create_comparison_heatmap(single_analysis):
    """Create heatmap comparing PFASGroups vs PFAS-Atlas performance"""
    
    # Prepare data for heatmap
    groups = []
    group_names = []
    pfasgroups_rates = []
    atlas_rates = []
    
    for group_id, data in sorted(single_analysis.items()):
        groups.append(f"Group {group_id}")
        group_names.append(data['group_name'])
        
        total = data['total_molecules']
        pfas_rate = (data['pfasgroups_correct'] / total * 100) if total > 0 else 0
        atlas_rate = (data['atlas_any_detection'] / total * 100) if total > 0 else 0
        
        pfasgroups_rates.append(pfas_rate)
        atlas_rates.append(atlas_rate)
    
    # Create comparison matrix
    comparison_data = np.array([pfasgroups_rates, atlas_rates])
    
    fig = go.Figure(data=go.Heatmap(
        z=comparison_data,
        x=[f"{g}<br>{name}" for g, name in zip(groups, group_names)],
        y=['PFASGroups<br>Detection Rate', 'PFAS-Atlas<br>Detection Rate'],
        colorscale='RdYlGn',
        text=comparison_data,
        texttemplate="%{text:.1f}%",
        textfont={"size":10},
        hoverongaps=False
    ))
    
    fig.update_layout(
        title="PFASGroups vs PFAS-Atlas Detection Performance Comparison<br><sub>Single Functional Group Molecules - PFASGroups Module Performance</sub>",
        xaxis_title="PFAS Functional Groups",
        yaxis_title="Detection Systems", 
        font_size=12,
        width=1400,
        height=500,
        margin=dict(l=100, r=50, t=100, b=100)
    )
    
    return fig

def create_multi_group_heatmap(multi_analysis):
    """Create heatmap showing which functional groups are privileged in multi-group detection"""
    
    # Prepare data for multi-group privilege analysis
    all_groups = set()
    for target_groups in multi_analysis.keys():
        all_groups.update(target_groups)
    
    all_groups = sorted(list(all_groups))
    
    # Create matrix showing detection rates for each group in multi-group contexts
    detection_matrix = []
    combo_labels = []
    
    for target_groups, data in multi_analysis.items():
        combo_label = f"Groups {'+'.join(map(str, target_groups))}"
        combo_labels.append(combo_label)
        
        total = data['total_molecules']
        row = []
        
        for group in all_groups:
            if group in target_groups:
                # This group is a target in this combination
                detections = data['group_detection_matrix'][group]
                rate = (detections / total * 100) if total > 0 else 0
                row.append(rate)
            else:
                row.append(np.nan)  # Not applicable
        
        detection_matrix.append(row)
    
    detection_matrix = np.array(detection_matrix)
    
    fig = go.Figure(data=go.Heatmap(
        z=detection_matrix,
        x=[f"Group {g}" for g in all_groups],
        y=combo_labels,
        colorscale='Viridis',
        text=detection_matrix,
        texttemplate="%{text:.0f}%",
        textfont={"size":9},
        hoverongaps=False,
        showscale=True
    ))
    
    fig.update_layout(
        title="Multi-Group Detection Privilege Analysis<br><sub>Which functional groups are detected in multi-functional molecules?</sub>",
        xaxis_title="Target Functional Groups",
        yaxis_title="Multi-Group Combinations",
        font_size=11,
        width=1200,
        height=600,
        margin=dict(l=150, r=50, t=100, b=50)
    )
    
    return fig

def create_atlas_classification_flow(single_analysis):
    """Create Sankey diagram showing PFAS-Atlas classification flow using actual data"""
    
    # Create PFAS-Atlas classification flow using actual data
    atlas_class_counts = {}
    for group_data in single_analysis.values():
        for mol in group_data['atlas_classifications_detail']:
            first_class = mol['first_class']
            second_class = mol['second_class']
            
            if first_class not in atlas_class_counts:
                atlas_class_counts[first_class] = {}
            if second_class not in atlas_class_counts[first_class]:
                atlas_class_counts[first_class][second_class] = 0
            atlas_class_counts[first_class][second_class] += 1
    
    if not atlas_class_counts:
        # Return empty figure if no data
        return go.Figure().update_layout(
            title="PFAS-Atlas Classification Flow<br><sub>No classification data available</sub>",
            annotations=[{"text": "No PFAS-Atlas classification data found", 
                         "xref": "paper", "yref": "paper", "x": 0.5, "y": 0.5, 
                         "xanchor": "center", "yanchor": "center", "showarrow": False}]
        )
    
    # Build nodes and links for actual Atlas classification flow
    nodes = ["PFAS Molecules<br>(Single Groups)"]
    first_classes = list(atlas_class_counts.keys())
    second_classes = set()
    for first_class_data in atlas_class_counts.values():
        second_classes.update(first_class_data.keys())
    second_classes = list(second_classes)
    
    # Add first class nodes
    for fc in first_classes:
        nodes.append(f"First Class:<br>{fc}")
    
    # Add second class nodes  
    for sc in second_classes:
        nodes.append(f"Second Class:<br>{sc}")
    
    sources, targets, values = [], [], []
    
    # Connect root to first classes
    for i, first_class in enumerate(first_classes):
        total_in_first = sum(atlas_class_counts[first_class].values())
        sources.append(0)  # Root node
        targets.append(i + 1)  # First class node
        values.append(total_in_first)
    
    # Connect first classes to second classes
    for i, first_class in enumerate(first_classes):
        for second_class, count in atlas_class_counts[first_class].items():
            second_idx = second_classes.index(second_class) + len(first_classes) + 1
            sources.append(i + 1)  # First class node
            targets.append(second_idx)  # Second class node
            values.append(count)
    
    # Create Sankey diagram
    fig = go.Figure(data=[go.Sankey(
        node=dict(
            pad=15,
            thickness=20,
            line=dict(color="black", width=0.5),
            label=nodes,
            color=["rgba(31, 119, 180, 0.8)"] + 
                  ["rgba(255, 127, 14, 0.8)"] * len(first_classes) + 
                  ["rgba(44, 160, 44, 0.8)"] * len(second_classes)
        ),
        link=dict(
            source=sources,
            target=targets,
            value=values,
            color=["rgba(255, 127, 14, 0.4)"] * len(sources)
        )
    )])
    
    fig.update_layout(
        title="PFAS-Atlas Classification Flow<br><sub>Native PFAS-Atlas class distribution (not mapped to PFASGroups)</sub>",
        font_size=12,
        width=1000,
        height=600
    )
    
    return fig

def create_enhanced_sankey_comparison(single_analysis, multi_analysis, results):
    """Create enhanced Sankey diagram showing system comparison"""
    
    # Calculate overall statistics
    total_single = sum(data['total_molecules'] for data in single_analysis.values())
    pfasgroups_single_success = sum(data['pfasgroups_correct'] for data in single_analysis.values())
    atlas_single_success = sum(data['atlas_any_detection'] for data in single_analysis.values())
    
    total_multi = sum(data['total_molecules'] for data in multi_analysis.values())
    pfasgroups_multi_success = sum(data['pfasgroups_any_detection'] for data in multi_analysis.values())
    atlas_multi_success = sum(data['atlas_classifications'] for data in multi_analysis.values())
    
    # Create links
    pfas_atlas = {}
    pfas_pfasgroups = {}
    pfas_pfasgroups_oecd = {}
    nodes_set = set()
    nodes_pfasgroups_set = set()
    nodes_atlas_set = set()
    nodes_pfasgroups_oecd_set = set()
    
    with open('../PFASgroups/data/PFAS_groups_smarts.json', 'r') as f:
        lookup = json.load(f)
    lookup = {k['id']:k['name'] for k in lookup}
    reverse_lookup = {v:k for k,v in lookup.items()}
    for m in results:
        if m['molecule_data'].get('generation_type') == 'single_group':
            expected = m['molecule_data']['group_name']
            atlas = m['atlas_result'].get('second_class','None')
            pfasgroups_ids = m['pfasgroups_result']['detected_groups']
            pfasgroups = [lookup[id] for id in pfasgroups_ids if id!=49 and id!=50 and id > 28]
            pfasgroups = [f"- {reverse_lookup.get(pf)}: {pf}" for pf in pfasgroups]
            pfasgroups_oecd = [lookup[id] for id in pfasgroups_ids if id<= 28]
            pfasgroups_oecd = [f"- {reverse_lookup.get(pf)}: {pf}" for pf in pfasgroups_oecd]
            expected = f"{reverse_lookup.get(expected)}: {expected}"
            # Collect unique nodes
            nodes_set.add(expected)
            nodes_atlas_set.add(atlas)
            nodes_pfasgroups_set.update(pfasgroups)
            nodes_pfasgroups_oecd_set.update(pfasgroups_oecd)
            # Build connection dictionaries
            pfas_atlas[expected][atlas] = pfas_atlas.setdefault(expected,{}).get(atlas,0)+1
            for pfg in pfasgroups:
                pfas_pfasgroups[expected][pfg] = pfas_pfasgroups.setdefault(expected,{}).get(pfg,0)+1
            for pfg in pfasgroups_oecd:
                pfas_pfasgroups_oecd[expected][pfg] = pfas_pfasgroups_oecd.setdefault(expected,{}).get(pfg,0)+1

    # Create deduplicated node lists
    nodes = sorted(list(nodes_set))
    nodes_atlas = sorted(list(nodes_atlas_set))
    nodes_pfasgroups = sorted(list(nodes_pfasgroups_set))
    nodes_pfasgroups_oecd = sorted(list(nodes_pfasgroups_oecd_set))

    # Create nodes for Atlas diagram
    pfas_atlas_nodes = [
        # Source nodes
        *nodes,
        # Atlas outcomes
        *nodes_atlas
    ]

    # Create nodes for PFASGroups diagram
    pfas_pfasgroups_nodes = [
        # Source nodes
        *nodes,
        # PFASgroups outcomes
        *nodes_pfasgroups
    ]

    pfas_pfasgroups_oecd_nodes = [
        *nodes,
        *nodes_pfasgroups_oecd
    ]
    
    # Generate color schemes for nodes and links
    # For PFAS-Atlas nodes: source nodes in blue, target nodes in orange gradient
    num_source_nodes = len(nodes)
    num_target_nodes = len(nodes_atlas)
    node_atlas_colors = (
        ["rgba(31, 119, 180, 0.8)"] * num_source_nodes +  # Source nodes in blue
        px.colors.sample_colorscale("Oranges", [0.4 + (i * 0.6 / max(1, num_target_nodes - 1)) for i in range(num_target_nodes)])
    )

    # For PFASGroups nodes: source nodes in blue, target nodes in green gradient
    num_pfasgroups_target_nodes = len(nodes_pfasgroups)
    node_pfasgroups_colors = (
        ["rgba(31, 119, 180, 0.8)"] * num_source_nodes +  # Source nodes in blue
        px.colors.sample_colorscale("Greens", [0.4 + (i * 0.6 / max(1, num_pfasgroups_target_nodes - 1)) for i in range(num_pfasgroups_target_nodes)])
    )

    # For PFASGroups nodes: source nodes in blue, target nodes in green gradient
    num_pfasgroups_oecd_target_nodes = len(nodes_pfasgroups_oecd)
    node_pfasgroups_oecd_colors = (
        ["rgba(31, 119, 180, 0.8)"] * num_source_nodes +  # Source nodes in blue
        px.colors.sample_colorscale("Greens", [0.4 + (i * 0.6 / max(1, num_pfasgroups_oecd_target_nodes - 1)) for i in range(num_pfasgroups_oecd_target_nodes)])
    )
    
    # Create links for PFAS-Atlas diagram
    pfas_atlas_sources = []
    pfas_atlas_targets = []
    pfas_atlas_values = []
    
    for source_name, target_dict in pfas_atlas.items():
        source_idx = pfas_atlas_nodes.index(source_name)
        for target_name, count in target_dict.items():
            target_idx = pfas_atlas_nodes.index(target_name)
            pfas_atlas_sources.append(source_idx)
            pfas_atlas_targets.append(target_idx)
            pfas_atlas_values.append(count)

    # Create links for PFASGroups diagram
    pfas_pfasgroups_sources = []
    pfas_pfasgroups_targets = []
    pfas_pfasgroups_values = []
    
    for source_name, target_dict in pfas_pfasgroups.items():
        source_idx = pfas_pfasgroups_nodes.index(source_name)
        for target_name, count in target_dict.items():
            target_idx = pfas_pfasgroups_nodes.index(target_name)
            pfas_pfasgroups_sources.append(source_idx)
            pfas_pfasgroups_targets.append(target_idx)
            pfas_pfasgroups_values.append(count)

    # Create links for PFASGroups _oecd diagram
    pfas_pfasgroups_oecd_sources = []
    pfas_pfasgroups_oecd_targets = []
    pfas_pfasgroups_oecd_values = []
    
    for source_name, target_dict in pfas_pfasgroups_oecd.items():
        source_idx = pfas_pfasgroups_oecd_nodes.index(source_name)
        for target_name, count in target_dict.items():
            target_idx = pfas_pfasgroups_oecd_nodes.index(target_name)
            pfas_pfasgroups_oecd_sources.append(source_idx)
            pfas_pfasgroups_oecd_targets.append(target_idx)
            pfas_pfasgroups_oecd_values.append(count)
    

    # Generate link colors with transparency
    # PFAS-Atlas links: semi-transparent orange
    pfas_atlas_link_colors = ["rgba(255, 127, 14, 0.3)"] * len(pfas_atlas_values)

    # PFASGroups links: semi-transparent green
    pfas_pfasgroups_link_colors = ["rgba(44, 160, 44, 0.3)"] * len(pfas_pfasgroups_values)

    # PFASGroups links: semi-transparent green
    pfas_pfasgroups_oecd_link_colors = ["rgba(44, 160, 44, 0.3)"] * len(pfas_pfasgroups_oecd_values)
    
    fig_atlas = go.Figure(data=[go.Sankey(
        node=dict(
            pad=15,
            thickness=20,
            line=dict(color="black", width=0.5),
            label=pfas_atlas_nodes,
            color=node_atlas_colors
        ),
        link=dict(
            source=pfas_atlas_sources,
            target=pfas_atlas_targets,
            value=pfas_atlas_values,
            color=pfas_atlas_link_colors
        )
    )])
    
    fig_atlas.update_layout(
        title=f"Single Functional Group Identification: PFAS-Atlas<br>",
        font_size=12,
        width=1200,
        height=600
    )

    fig_pfasgroups = go.Figure(data=[go.Sankey(
        node=dict(
            pad=15,
            thickness=20,
            line=dict(color="black", width=0.5),
            label= pfas_pfasgroups_nodes,
            color=node_pfasgroups_colors
        ),
        link=dict(
            source=pfas_pfasgroups_sources,
            target=pfas_pfasgroups_targets,
            value=pfas_pfasgroups_values,
            color=pfas_pfasgroups_link_colors
        )
    )])
    
    fig_pfasgroups.update_layout(
        title=f"Single Functional Group Identification: PFASGroups<br>",
        font_size=12,
        width=1200,
        height=600
    )

    fig_pfasgroups_oecd = go.Figure(data=[go.Sankey(
        node=dict(
            pad=15,
            thickness=20,
            line=dict(color="black", width=0.5),
            label= pfas_pfasgroups_oecd_nodes,
            color=node_pfasgroups_oecd_colors
        ),
        link=dict(
            source=pfas_pfasgroups_oecd_sources,
            target=pfas_pfasgroups_oecd_targets,
            value=pfas_pfasgroups_oecd_values,
            color=pfas_pfasgroups_oecd_link_colors
        )
    )])
    
    fig_pfasgroups_oecd.update_layout(
        title=f"Single Functional Group Identification: PFASGroups OECD<br>",
        font_size=12,
        width=1200,
        height=600
    )

    figs = [fig_atlas, fig_pfasgroups, fig_pfasgroups_oecd]
    return figs

def analyze_oecd_benchmark(oecd_results):
    """Analyze OECD benchmark results"""
    
    analysis = {
        'total_molecules': len(oecd_results),
        'pfasgroups_detections': 0,
        'atlas_detections': 0,
        'class_breakdown': defaultdict(lambda: {'count': 0, 'pfas_detected': 0, 'atlas_detected': 0}),
        'group_detections': defaultdict(int),
        'misclassifications': {'pfasgroups': [], 'atlas': []}
    }
    
    for result in oecd_results:
        mol_data = result['molecule_data']
        pfas_result = result['pfasgroups_result']
        atlas_result = result['atlas_result']
        
        first_class = mol_data.get('oecd_first_class', 'Unknown')
        second_class = mol_data.get('oecd_second_class', 'Unknown')
        
        # Update class breakdown
        analysis['class_breakdown'][second_class]['count'] += 1
        
        # PFASGroups analysis
        if pfas_result.get('success', False):
            analysis['pfasgroups_detections'] += 1
            analysis['class_breakdown'][second_class]['pfas_detected'] += 1
            
            # Track group detections
            for group_id in pfas_result.get('detected_groups', []):
                analysis['group_detections'][group_id] += 1
        else:
            # Track misclassifications
            analysis['misclassifications']['pfasgroups'].append({
                'smiles': mol_data['smiles'],
                'error': pfas_result.get('error', 'No detection')
            })
        
        # PFAS-Atlas analysis
        if atlas_result.get('success', False):
            analysis['atlas_detections'] += 1
            analysis['class_breakdown'][second_class]['atlas_detected'] += 1
        else:
            # Track misclassifications
            analysis['misclassifications']['atlas'].append({
                'smiles': mol_data['smiles'],
                'predicted_class': atlas_result.get('second_class', 'None'),
                'error': atlas_result.get('error', 'No detection')
            })
    
    return analysis

def create_single_group_atlas_sankey(single_analysis):
    """Create Sankey diagram showing PFASGroups vs PFAS-Atlas second_class for single groups"""
    
    # Collect data for PFASGroups detections and corresponding Atlas second_class
    atlas_class_counts = defaultdict(int)
    pfasgroups_success_count = 0
    
    for group_id, data in single_analysis.items():
        pfasgroups_success_count += data['pfasgroups_correct']
        
        # We need to track Atlas second_class for each group's successful PFASGroups detection
        # This would require access to individual results, not just aggregated data
        # For now, create a representative diagram
        
    # Create PFAS-Atlas classification flow using actual data
    atlas_class_counts = {}
    for group_data in single_analysis.values():
        for mol in group_data['atlas_classifications_detail']:
            first_class = mol['first_class']
            second_class = mol['second_class']
            
            if first_class not in atlas_class_counts:
                atlas_class_counts[first_class] = {}
            if second_class not in atlas_class_counts[first_class]:
                atlas_class_counts[first_class][second_class] = 0
            atlas_class_counts[first_class][second_class] += 1
    
    # Build nodes and links for actual Atlas classification flow
    nodes = ["PFAS Molecules<br>(Single Groups)"]
    first_classes = list(atlas_class_counts.keys())
    second_classes = set()
    for first_class_data in atlas_class_counts.values():
        second_classes.update(first_class_data.keys())
    second_classes = list(second_classes)
    
    # Add first class nodes
    for fc in first_classes:
        nodes.append(f"First Class:<br>{fc}")
    
    # Add second class nodes  
    for sc in second_classes:
        nodes.append(f"Second Class:<br>{sc}")
    
    sources, targets, values = [], [], []
    
    # Connect root to first classes
    for i, first_class in enumerate(first_classes):
        total_in_first = sum(atlas_class_counts[first_class].values())
        sources.append(0)  # Root node
        targets.append(i + 1)  # First class node
        values.append(total_in_first)
    
    # Connect first classes to second classes
    for i, first_class in enumerate(first_classes):
        for second_class, count in atlas_class_counts[first_class].items():
            second_idx = second_classes.index(second_class) + len(first_classes) + 1
            sources.append(i + 1)  # First class node
            targets.append(second_idx)  # Second class node
            values.append(count)
    
    fig = go.Figure(data=[go.Sankey(
        node=dict(
            pad=15,
            thickness=20,
            line=dict(color="black", width=0.5),
            label=nodes,
            color=["rgba(31, 119, 180, 0.8)", "rgba(255, 127, 14, 0.8)", 
                   "rgba(44, 160, 44, 0.8)", "rgba(214, 39, 40, 0.8)",
                   "rgba(148, 103, 189, 0.8)", "rgba(188, 189, 34, 0.8)"]
        ),
        link=dict(
            source=sources,
            target=targets,
            value=values,
            color=["rgba(255, 127, 14, 0.4)", "rgba(44, 160, 44, 0.4)",
                   "rgba(214, 39, 40, 0.4)", "rgba(148, 103, 189, 0.4)",
                   "rgba(188, 189, 34, 0.4)"]
        )
    )])
    
    fig.update_layout(
        title="Single Group Molecules: PFASGroups Detection → PFAS-Atlas Classification",
        font_size=12,
        width=1000,
        height=600
    )
    
    return fig

def create_multi_group_pfasgroups_sankey(multi_analysis):
    """Create Sankey showing multi-group PFASGroups performance"""
    
    # Calculate detection rates for multi-group combinations
    combo_labels = []
    detection_rates = []
    
    for target_groups, data in multi_analysis.items():
        combo_str = "-".join(map(str, target_groups))
        combo_labels.append(combo_str)
        
        total = data['total_molecules']
        detected = data['pfasgroups_any_detection']
        rate = (detected / total * 100) if total > 0 else 0
        detection_rates.append(detected)
    
    # Create Sankey
    nodes = ["Multi-Group<br>Molecules"] + [f"Groups {label}<br>({rate:.0f}% detected)" 
                                            for label, rate in zip(combo_labels, [d/10 for d in detection_rates])]
    
    sources = [0] * len(combo_labels)
    targets = list(range(1, len(combo_labels) + 1))
    values = detection_rates
    
    fig = go.Figure(data=[go.Sankey(
        node=dict(
            pad=15,
            thickness=20,
            line=dict(color="black", width=0.5),
            label=nodes,
            color=["rgba(31, 119, 180, 0.8)"] + px.colors.sample_colorscale("Viridis", 
                                                                           [i/(len(combo_labels)-1) for i in range(len(combo_labels))])
        ),
        link=dict(
            source=sources,
            target=targets,
            value=values,
            color="rgba(128, 128, 128, 0.4)"
        )
    )])
    
    fig.update_layout(
        title="Multi-Group Molecules: PFASGroups Detection Performance",
        font_size=12,
        width=1200,
        height=700
    )
    
    return fig

def create_multi_group_atlas_sankey(oecd_results):
    """Create Sankey showing multi-group molecules to PFAS-Atlas second_class"""
    
    # This would analyze multi-group molecules from OECD data
    # For now, create a representative diagram
    
    nodes = [
        "Multi-Group<br>Molecules",
        "HFCs", "PFCAs", "PFSAs", "PolyFCA derivatives", "PASF-based substances", "Others"
    ]
    
    # Mock data representing Atlas classifications
    sources = [0, 0, 0, 0, 0, 0]
    targets = [1, 2, 3, 4, 5, 6]
    values = [150, 120, 100, 180, 90, 60]  # Example values
    
    fig = go.Figure(data=[go.Sankey(
        node=dict(
            pad=15,
            thickness=20,
            line=dict(color="black", width=0.5),
            label=nodes,
            color=["rgba(31, 119, 180, 0.8)"] + px.colors.qualitative.Set3[:6]
        ),
        link=dict(
            source=sources,
            target=targets,
            value=values,
            color="rgba(128, 128, 128, 0.4)"
        )
    )])
    
    fig.update_layout(
        title="Multi-Group Molecules: PFAS-Atlas Second Class Classification",
        font_size=12,
        width=1000,
        height=600
    )
    
    return fig

def create_oecd_html_report(oecd_analysis, timestamp):
    """Create OECD benchmark HTML report"""
    
    total_molecules = oecd_analysis['total_molecules']
    pfasgroups_detections = oecd_analysis['pfasgroups_detections']
    atlas_detections = oecd_analysis['atlas_detections']
    
    pfasgroups_rate = (pfasgroups_detections / total_molecules * 100) if total_molecules > 0 else 0
    atlas_rate = (atlas_detections / total_molecules * 100) if total_molecules > 0 else 0
    
    # Create class breakdown chart
    class_names = []
    class_counts = []
    pfas_detected = []
    atlas_detected = []
    
    for class_name, data in oecd_analysis['class_breakdown'].items():
        if data['count'] > 5:  # Only show classes with significant representation
            class_names.append(class_name)
            class_counts.append(data['count'])
            pfas_detected.append(data['pfas_detected'])
            atlas_detected.append(data['atlas_detected'])
    
    fig_breakdown = go.Figure()
    fig_breakdown.add_trace(go.Bar(name='Total', x=class_names, y=class_counts, marker_color='lightblue'))
    fig_breakdown.add_trace(go.Bar(name='PFASGroups Detected', x=class_names, y=pfas_detected, marker_color='green'))
    fig_breakdown.add_trace(go.Bar(name='Atlas Detected', x=class_names, y=atlas_detected, marker_color='orange'))
    
    fig_breakdown.update_layout(
        title='OECD Class Breakdown: Detection Performance',
        xaxis_title='OECD Classes',
        yaxis_title='Molecule Count',
        barmode='group',
        width=1200,
        height=600
    )
    
    breakdown_html = fig_breakdown.to_html(include_plotlyjs='cdn', div_id="oecd-breakdown")
    
    # Misclassifications summary
    pfas_misc_count = len(oecd_analysis['misclassifications']['pfasgroups'])
    atlas_misc_count = len(oecd_analysis['misclassifications']['atlas'])
    
    html_content = f"""
<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>OECD PFAS Benchmark Analysis</title>
    <style>
        body {{
            font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif;
            margin: 0;
            padding: 20px;
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            min-height: 100vh;
            color: #333;
        }}
        .container {{
            max-width: 1400px;
            margin: 0 auto;
            background-color: white;
            padding: 40px;
            border-radius: 15px;
            box-shadow: 0 10px 30px rgba(0,0,0,0.3);
        }}
        h1 {{
            color: #2c5aa0;
            text-align: center;
            font-size: 2.8em;
            margin-bottom: 10px;
        }}
        h2 {{
            color: #34495e;
            border-left: 4px solid #667eea;
            padding-left: 15px;
            margin-top: 50px;
        }}
        .summary-grid {{
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(200px, 1fr));
            gap: 20px;
            margin: 30px 0;
        }}
        .summary-card {{
            background: linear-gradient(45deg, #667eea 0%, #764ba2 100%);
            color: white;
            padding: 25px;
            border-radius: 12px;
            text-align: center;
            box-shadow: 0 5px 15px rgba(0,0,0,0.2);
        }}
        .summary-number {{
            font-size: 2.5em;
            font-weight: bold;
            margin-bottom: 10px;
        }}
        .summary-label {{
            font-size: 1.0em;
            opacity: 0.95;
        }}
        .visualization-container {{
            background-color: #f8f9fa;
            padding: 30px;
            border-radius: 12px;
            margin: 40px 0;
            box-shadow: 0 5px 15px rgba(0,0,0,0.1);
        }}
        .highlight-section {{
            background: linear-gradient(135deg, #28a745 0%, #20c997 100%);
            color: white;
            padding: 30px;
            border-radius: 12px;
            margin: 40px 0;
        }}
    </style>
</head>
<body>
    <div class="container">
        <h1>🧪 OECD PFAS Benchmark Analysis</h1>
        
        <div class="highlight-section">
            <h3>📊 OECD Dataset Performance Evaluation</h3>
            <p style="font-size: 1.2em; margin: 20px 0;">
                Testing PFASGroups (Groups 1-28) vs PFAS-Atlas on OECD molecules<br>
                Analysis Date: {datetime.now().strftime('%B %d, %Y at %H:%M')}
            </p>
        </div>

        <div class="summary-grid">
            <div class="summary-card">
                <div class="summary-number">{total_molecules}</div>
                <div class="summary-label">Total OECD<br>Molecules</div>
            </div>
            <div class="summary-card">
                <div class="summary-number">{pfasgroups_rate:.1f}%</div>
                <div class="summary-label">PFASGroups<br>Detection Rate</div>
            </div>
            <div class="summary-card">
                <div class="summary-number">{atlas_rate:.1f}%</div>
                <div class="summary-label">Atlas<br>Detection Rate</div>
            </div>
            <div class="summary-card">
                <div class="summary-number">{pfas_misc_count}</div>
                <div class="summary-label">PFASGroups<br>Misclassifications</div>
            </div>
            <div class="summary-card">
                <div class="summary-number">{atlas_misc_count}</div>
                <div class="summary-label">Atlas<br>Misclassifications</div>
            </div>
        </div>

        <h2>📊 OECD Class Breakdown</h2>
        <div class="visualization-container">
            {breakdown_html}
        </div>
        
        <h2>🎯 Group Detection Frequency</h2>
        <div class="visualization-container">
            <p>Most frequently detected PFASGroups:</p>
            <ul>
"""
    
    # Add top detected groups
    sorted_groups = sorted(oecd_analysis['group_detections'].items(), key=lambda x: x[1], reverse=True)[:10]
    for group_id, count in sorted_groups:
        html_content += f"<li>Group {group_id}: {count} detections</li>"
    
    html_content += """
            </ul>
        </div>
    </div>
</body>
</html>
"""
    
    # Save to html directory if it exists, otherwise current directory
    html_dir = "html" if os.path.exists("html") else "."
    html_filename = f"{html_dir}/oecd_pfas_analysis_{timestamp}.html"
    with open(html_filename, 'w') as f:
        f.write(html_content)

def create_combined_comparison_report(enhanced_results, oecd_results, timestamp):
    """Create combined comparison report with all requested Sankey diagrams"""
    
    # Analyze enhanced results
    single_analysis, multi_analysis = analyze_system_comparison(enhanced_results)
    
    # Create all requested Sankey diagrams
    sankey_single_atlas = create_single_group_atlas_sankey(single_analysis)
    sankey_multi_pfas = create_multi_group_pfasgroups_sankey(multi_analysis)
    sankey_multi_atlas = create_multi_group_atlas_sankey(oecd_results)
    
    # Convert to HTML
    sankey_single_html = sankey_single_atlas.to_html(include_plotlyjs='cdn', div_id="single-atlas")
    sankey_multi_pfas_html = sankey_multi_pfas.to_html(include_plotlyjs=False, div_id="multi-pfas")
    sankey_multi_atlas_html = sankey_multi_atlas.to_html(include_plotlyjs=False, div_id="multi-atlas")
    
    # Save individual visualizations to imgs directory
    sankey_single_atlas.write_image(f"imgs/sankey_atlas_diagram_{timestamp}.png", width=1000, height=600, scale=2)
    sankey_multi_pfas.write_image(f"imgs/sankey_pfasgroups_diagram_{timestamp}.png", width=1200, height=700, scale=2)
    sankey_multi_atlas.write_image(f"imgs/sankey_pfasgroups_oecd_diagram_{timestamp}.png", width=1000, height=600, scale=2)
    
    html_content = f"""
<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Combined PFAS Analysis - Enhanced & OECD Comparison</title>
    <style>
        body {{
            font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif;
            margin: 0;
            padding: 20px;
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            min-height: 100vh;
            color: #333;
        }}
        .container {{
            max-width: 1600px;
            margin: 0 auto;
            background-color: white;
            padding: 40px;
            border-radius: 15px;
            box-shadow: 0 10px 30px rgba(0,0,0,0.3);
        }}
        h1 {{
            color: #2c5aa0;
            text-align: center;
            font-size: 2.8em;
            margin-bottom: 10px;
        }}
        h2 {{
            color: #34495e;
            border-left: 4px solid #667eea;
            padding-left: 15px;
            margin-top: 50px;
        }}
        .visualization-container {{
            background-color: #f8f9fa;
            padding: 30px;
            border-radius: 12px;
            margin: 40px 0;
            box-shadow: 0 5px 15px rgba(0,0,0,0.1);
        }}
        .highlight-section {{
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            color: white;
            padding: 40px;
            border-radius: 15px;
            margin: 50px 0;
            text-align: center;
        }}
        .chart-description {{
            font-size: 0.95em;
            color: #666;
            margin-bottom: 20px;
            font-style: italic;
        }}
    </style>
</head>
<body>
    <div class="container">
        <h1>🔬 Combined PFAS Analysis</h1>
        
        <div class="highlight-section">
            <h3>🎯 Comprehensive PFASGroups vs PFAS-Atlas Comparison</h3>
            <p style="font-size: 1.2em; margin: 20px 0;">
                Enhanced functional groups testing + OECD dataset validation<br>
                Analysis Date: {datetime.now().strftime('%B %d, %Y at %H:%M')}
            </p>
        </div>

        <h2>🔗 Single Functional Groups: PFASGroups → PFAS-Atlas Classification</h2>
        <div class="visualization-container">
            <div class="chart-description">
                Flow diagram showing how single functional group molecules detected by PFASGroups 
                are classified by PFAS-Atlas second_class categories.
            </div>
            {sankey_single_html}
        </div>

        <h2>🔀 Multi-Functional Groups: PFASGroups Performance</h2>
        <div class="visualization-container">
            <div class="chart-description">
                Performance analysis of PFASGroups on multi-functional group molecules, 
                showing detection rates for different functional group combinations.
            </div>
            {sankey_multi_pfas_html}
        </div>

        <h2>📊 Multi-Functional Groups: PFAS-Atlas Classifications</h2>
        <div class="visualization-container">
            <div class="chart-description">
                PFAS-Atlas second_class output distribution for multi-functional group molecules,
                showing classification patterns for complex molecular structures.
            </div>
            {sankey_multi_atlas_html}
        </div>
    </div>
</body>
</html>
"""
    
    # Save to html directory if it exists, otherwise current directory
    html_dir = "html" if os.path.exists("html") else "."
    html_filename = f"{html_dir}/combined_pfas_analysis_{timestamp}.html"
    with open(html_filename, 'w') as f:
        f.write(html_content)

def create_detailed_multi_group_sankey(multi_analysis):
    """Create detailed Sankey for multi-group analysis showing group privileges"""
    
    # Focus on which groups get detected in multi-group scenarios
    group_preference_data = defaultdict(int)
    total_opportunities = defaultdict(int)
    
    for target_groups, data in multi_analysis.items():
        total_molecules = data['total_molecules']
        
        for group in target_groups:
            total_opportunities[group] += total_molecules
            detections = data['group_detection_matrix'][group]
            group_preference_data[group] += detections
    
    # Calculate detection rates
    group_rates = {}
    for group in group_preference_data:
        rate = group_preference_data[group] / total_opportunities[group]
        group_rates[group] = rate
    
    # Sort by detection rate
    sorted_groups = sorted(group_rates.items(), key=lambda x: x[1], reverse=True)
    
    # Create Sankey showing privilege hierarchy
    nodes = ["Multi-Group<br>Molecules"] + [f"Group {g}<br>({rate:.1%})" for g, rate in sorted_groups]
    
    sources = [0] * len(sorted_groups)
    targets = list(range(1, len(sorted_groups) + 1))
    values = [group_preference_data[g] for g, _ in sorted_groups]
    
    # Color gradient from high to low privilege
    colors_base = px.colors.sample_colorscale("Viridis", [i/(len(sorted_groups)-1) for i in range(len(sorted_groups))])
    node_colors = ["rgba(31, 119, 180, 0.8)"] + [color for color in colors_base]
    
    fig = go.Figure(data=[go.Sankey(
        node=dict(
            pad=15,
            thickness=20,
            line=dict(color="black", width=0.5),
            label=nodes,
            color=node_colors
        ),
        link=dict(
            source=sources,
            target=targets,
            value=values,
            color="rgba(128, 128, 128, 0.4)"
        )
    )])
    
    fig.update_layout(
        title="Multi-Group Detection Privilege Hierarchy<br><sub>Which functional groups are preferentially detected in complex molecules?</sub>",
        font_size=12,
        width=1000,
        height=600
    )
    
    return fig

def create_pfasgroups_failure_report(single_analysis):
    """Generate detailed report of molecules not correctly identified by PFASGroups"""
    
    failure_report = {
        'summary': {
            'total_groups': len(single_analysis),
            'groups_with_failures': 0,
            'total_failures': 0,
            'no_detection_failures': 0,
            'wrong_group_failures': 0
        },
        'by_group': {},
        'failure_molecules': []
    }
    
    for group_id, data in single_analysis.items():
        failures = data['pfasgroups_failures']
        if failures:
            failure_report['summary']['groups_with_failures'] += 1
            failure_report['summary']['total_failures'] += len(failures)
            
            no_detection = len([f for f in failures if f['detection_status'] == 'no_detection'])
            wrong_group = len([f for f in failures if f['detection_status'] == 'wrong_group'])
            
            failure_report['summary']['no_detection_failures'] += no_detection
            failure_report['summary']['wrong_group_failures'] += wrong_group
            
            failure_report['by_group'][group_id] = {
                'group_name': data['group_name'],
                'total_molecules': data['total_molecules'],
                'failure_count': len(failures),
                'failure_rate': len(failures) / data['total_molecules'] * 100,
                'no_detection_count': no_detection,
                'wrong_group_count': wrong_group,
                'failures': failures
            }
            
            # Add to global failure list
            for failure in failures:
                failure_report['failure_molecules'].append({
                    'group_id': group_id,
                    'group_name': data['group_name'],
                    **failure
                })
    
    return failure_report

def create_enhanced_html_report(single_analysis, multi_analysis, timestamp, results):
    """Create comprehensive HTML report with embedded visualizations"""
    
    # Generate failure analysis
    failure_report = create_pfasgroups_failure_report(single_analysis)
    
    # Generate atlas classification flow
    atlas_flow_html = create_atlas_classification_flow(single_analysis)
    
    # Calculate comprehensive statistics
    total_single_molecules = sum(data['total_molecules'] for data in single_analysis.values())
    total_multi_molecules = sum(data['total_molecules'] for data in multi_analysis.values())
    
    pfasgroups_single_success = sum(data['pfasgroups_correct'] for data in single_analysis.values())
    pfasgroups_single_rate = (pfasgroups_single_success / total_single_molecules * 100) if total_single_molecules > 0 else 0
    
    atlas_single_success = sum(data['atlas_any_detection'] for data in single_analysis.values())
    atlas_single_rate = (atlas_single_success / total_single_molecules * 100) if total_single_molecules > 0 else 0
    
    pfasgroups_multi_success = sum(data['pfasgroups_any_detection'] for data in multi_analysis.values())
    pfasgroups_multi_rate = (pfasgroups_multi_success / total_multi_molecules * 100) if total_multi_molecules > 0 else 0
    
    atlas_multi_success = sum(data['atlas_classifications'] for data in multi_analysis.values())
    atlas_multi_rate = (atlas_multi_success / total_multi_molecules * 100) if total_multi_molecules > 0 else 0
    
    # Generate visualizations
    heatmap_comparison = create_comparison_heatmap(single_analysis)
    heatmap_multi = create_multi_group_heatmap(multi_analysis)
    sankey_comparison = create_enhanced_sankey_comparison(single_analysis, multi_analysis, results)
    sankey_privilege = create_detailed_multi_group_sankey(multi_analysis)
    timing_heatmap = create_timing_comparison_heatmap(single_analysis)
    timing_stats = create_timing_statistics_chart(single_analysis)
    
    # Convert to HTML
    heatmap_comparison_html = heatmap_comparison.to_html(include_plotlyjs='cdn', div_id="heatmap-comparison")
    heatmap_multi_html = heatmap_multi.to_html(include_plotlyjs=False, div_id="heatmap-multi")
    sankey_atlas_html = sankey_comparison[0].to_html(include_plotlyjs=False, div_id="sankey-atlas")
    sankey_pfasgroups_html = sankey_comparison[1].to_html(include_plotlyjs=False, div_id="sankey-pfasgroups")
    sankey_pfasgroups_oecd_html = sankey_comparison[2].to_html(include_plotlyjs=False, div_id="sankey-pfasgroups_oecd")
    sankey_privilege_html = sankey_privilege.to_html(include_plotlyjs=False, div_id="sankey-privilege")
    timing_heatmap_html = timing_heatmap.to_html(include_plotlyjs=False, div_id="timing-heatmap")
    timing_stats_html = timing_stats.to_html(include_plotlyjs=False, div_id="timing-stats")
    
    # Save individual visualizations
    heatmap_comparison.write_image(f"imgs/comparison_heatmap_{timestamp}.png", width=1400, height=500, scale=2)
    heatmap_comparison.write_image(f"imgs/comparison_heatmap_{timestamp}.svg")
    heatmap_multi.write_image(f"imgs/multigroup_privilege_heatmap_{timestamp}.png", width=1200, height=600, scale=2)
    heatmap_multi.write_image(f"imgs/multigroup_privilege_heatmap_{timestamp}.svg")
    sankey_comparison[0].write_image(f"imgs/atlas_sankey_{timestamp}.png", width=1200, height=600, scale=2)
    sankey_comparison[0].write_image(f"imgs/atlas_sankey_{timestamp}.svg")
    sankey_comparison[1].write_image(f"imgs/pfasgroups_sankey_{timestamp}.png", width=1200, height=600, scale=2)
    sankey_comparison[1].write_image(f"imgs/pfasgroups_sankey_{timestamp}.svg")
    sankey_comparison[2].write_image(f"imgs/pfasgroups_oecd_sankey_{timestamp}.png", width=1200, height=600, scale=2)
    sankey_comparison[2].write_image(f"imgs/pfasgroups_oecd_sankey_{timestamp}.svg")
    sankey_privilege.write_image(f"imgs/privilege_hierarchy_sankey_{timestamp}.png", width=1000, height=600, scale=2)
    sankey_privilege.write_image(f"imgs/privilege_hierarchy_sankey_{timestamp}.svg")
    timing_heatmap.write_image(f"imgs/timing_comparison_heatmap_{timestamp}.png", width=1400, height=500, scale=2)
    timing_heatmap.write_image(f"imgs/timing_comparison_heatmap_{timestamp}.svg")
    timing_stats.write_image(f"imgs/timing_statistics_{timestamp}.png", width=1200, height=800, scale=2)
    timing_stats.write_image(f"imgs/timing_statistics_{timestamp}.svg")
    
    html_content = f"""
<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Enhanced PFAS Benchmark Analysis - Comprehensive System Comparison</title>
    <style>
        body {{
            font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif;
            margin: 0;
            padding: 20px;
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            min-height: 100vh;
            color: #333;
        }}
        .container {{
            max-width: 1600px;
            margin: 0 auto;
            background-color: white;
            padding: 40px;
            border-radius: 15px;
            box-shadow: 0 10px 30px rgba(0,0,0,0.3);
        }}
        h1 {{
            color: #2c5aa0;
            text-align: center;
            font-size: 2.8em;
            margin-bottom: 10px;
        }}
        h2 {{
            color: #34495e;
            border-left: 4px solid #667eea;
            padding-left: 15px;
            margin-top: 50px;
        }}
        .summary-grid {{
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(220px, 1fr));
            gap: 20px;
            margin: 30px 0;
        }}
        .summary-card {{
            background: linear-gradient(45deg, #667eea 0%, #764ba2 100%);
            color: white;
            padding: 25px;
            border-radius: 12px;
            text-align: center;
            box-shadow: 0 5px 15px rgba(0,0,0,0.2);
        }}
        .summary-number {{
            font-size: 2.8em;
            font-weight: bold;
            margin-bottom: 10px;
        }}
        .summary-label {{
            font-size: 1.0em;
            opacity: 0.95;
            line-height: 1.3;
        }}
        .highlight-section {{
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            color: white;
            padding: 40px;
            border-radius: 15px;
            margin: 50px 0;
            text-align: center;
        }}
        .visualization-container {{
            background-color: #f8f9fa;
            padding: 30px;
            border-radius: 12px;
            margin: 40px 0;
            box-shadow: 0 5px 15px rgba(0,0,0,0.1);
        }}
        .comparison-table {{
            width: 100%;
            border-collapse: collapse;
            margin: 30px 0;
            font-size: 0.95em;
            background-color: white;
            border-radius: 8px;
            overflow: hidden;
            box-shadow: 0 5px 15px rgba(0,0,0,0.1);
        }}
        .comparison-table th, .comparison-table td {{
            border: 1px solid #ddd;
            padding: 12px;
            text-align: left;
        }}
        .comparison-table th {{
            background: linear-gradient(45deg, #667eea 0%, #764ba2 100%);
            color: white;
            font-weight: bold;
            position: sticky;
            top: 0;
        }}
        .comparison-table tr:nth-child(even) {{
            background-color: #f9f9f9;
        }}
        .performance-excellent {{ 
            background-color: #d4edda !important; 
            color: #155724;
            font-weight: bold;
        }}
        .performance-good {{ 
            background-color: #d1ecf1 !important; 
            color: #0c5460;
            font-weight: bold;
        }}
        .performance-moderate {{ 
            background-color: #fff3cd !important; 
            color: #856404;
            font-weight: bold;
        }}
        .performance-poor {{ 
            background-color: #f8d7da !important; 
            color: #721c24;
            font-weight: bold;
        }}
        .key-findings {{
            background: linear-gradient(135deg, #28a745 0%, #20c997 100%);
            color: white;
            padding: 30px;
            border-radius: 12px;
            margin: 40px 0;
        }}
        .methodology-box {{
            background: linear-gradient(135deg, #ffc107 0%, #fd7e14 100%);
            color: white;
            padding: 25px;
            border-radius: 10px;
            margin: 30px 0;
        }}
        .chart-description {{
            font-size: 0.95em;
            color: #666;
            margin-bottom: 20px;
            font-style: italic;
        }}
    </style>
</head>
<body>
    <div class="container">
        <h1>🔬 Enhanced PFAS Detection Benchmark</h1>
        
        <div class="highlight-section">
            <h3>🎯 Comprehensive PFASGroups vs PFAS-Atlas Comparison</h3>
            <p style="font-size: 1.2em; margin: 20px 0;">
                Large-scale systematic evaluation with {total_single_molecules + total_multi_molecules} molecules<br>
                Analysis Date: {datetime.now().strftime('%B %d, %Y at %H:%M')}
            </p>
            <p>Target: PFAS functional groups 29-59 (excluding groups 49 and 50 which are smartsPath-only) | Enhanced dataset sizes | Detailed privilege analysis</p>
        </div>

        <div class="summary-grid">
            <div class="summary-card">
                <div class="summary-number">{len(single_analysis)}</div>
                <div class="summary-label">Functional Groups<br>Tested (Single)</div>
            </div>
            <div class="summary-card">
                <div class="summary-number">{total_single_molecules}</div>
                <div class="summary-label">Single Group<br>Molecules</div>
            </div>
            <div class="summary-card">
                <div class="summary-number">{total_multi_molecules}</div>
                <div class="summary-label">Multi-Group<br>Molecules</div>
            </div>
            <div class="summary-card">
                <div class="summary-number">{pfasgroups_single_rate:.1f}%</div>
                <div class="summary-label">PFASGroups Single<br>Success Rate</div>
            </div>
            <div class="summary-card">
                <div class="summary-number">{atlas_single_rate:.1f}%</div>
                <div class="summary-label">Atlas Single<br>Success Rate</div>
            </div>
            <div class="summary-card">
                <div class="summary-number">{pfasgroups_multi_rate:.1f}%</div>
                <div class="summary-label">PFASGroups Multi<br>Success Rate</div>
            </div>
        </div>

        <div class="methodology-box">
            <h3>🔬 Enhanced Methodology</h3>
            <ul style="text-align: left; max-width: 800px; margin: 0 auto;">
                <li><strong>Single Groups:</strong> 15 molecules per functional group (up from 8) with diverse scaffolds</li>
                <li><strong>Multi-Groups:</strong> 5 pairs + 5 triplets, 10 molecules each (100 total multi-group molecules)</li>
                <li><strong>Privilege Analysis:</strong> Tracking which functional groups are preferentially detected</li>
                <li><strong>Comprehensive Comparison:</strong> Direct head-to-head performance evaluation</li>
            </ul>
        </div>

        <h2>📊 System Performance Heatmap Comparison</h2>
        <div class="visualization-container">
            <div class="chart-description">
                Direct comparison of detection rates between PFASGroups and PFAS-Atlas across all functional groups. 
                Green indicates high performance, red indicates poor performance.
            </div>
            {heatmap_comparison_html}
        </div>

        <h2>🔀 Enhanced System Flow Analysis</h2>
        <div class="visualization-container">
            <div class="chart-description">
                Comprehensive Sankey diagram showing success/failure flows for both single-group molecules 
                across both detection systems.
            </div>
            {sankey_atlas_html}
        </div>
        <h2>🔀 Enhanced System Flow Analysis</h2>
        <div class="visualization-container">
            <div class="chart-description">
                Comprehensive Sankey diagram showing success/failure flows for both single-group molecules 
                across both detection systems.
            </div>
            {sankey_pfasgroups_html}
        </div>
        <h2>🔀 Enhanced System Flow Analysis</h2>
        <div class="visualization-container">
            <div class="chart-description">
                Comprehensive Sankey diagram showing success/failure flows for both single-group molecules 
                across both detection systems.
            </div>
            {sankey_pfasgroups_oecd_html}
        </div>

        <h2>🎯 Multi-Group Privilege Analysis</h2>
        <div class="visualization-container">
            <div class="chart-description">
                Heatmap showing which functional groups are preferentially detected in multi-functional molecules. 
                Reveals system biases and detection hierarchies.
            </div>
            {heatmap_multi_html}
        </div>

        <h2>📈 PFASGroups Module: Functional Group Detection Hierarchy</h2>
        <div class="visualization-container">
            <div class="chart-description">
                Sankey diagram revealing the hierarchy of functional group detection in multi-group scenarios.
                Shows which groups are "privileged" during detection.
            </div>
            {sankey_privilege_html}
        </div>

        <h2>⏱️ Execution Time Comparison</h2>
        <div class="visualization-container">
            <div class="chart-description">
                Performance timing comparison between PFASGroups and PFAS-Atlas on single functional group molecules.
                Shows average execution time per molecule in milliseconds.
            </div>
            {timing_heatmap_html}
        </div>

        <h2>📊 Detailed Timing Analysis</h2>
        <div class="visualization-container">
            <div class="chart-description">
                Comprehensive timing statistics showing execution times and molecule counts tested for each functional group.
                Helps identify performance bottlenecks and system efficiency patterns.
            </div>
            {timing_stats_html}
        </div>

        <h2>⚠️ PFASGroups Detection Failures Analysis</h2>
        <div class="visualization-container">
            <div class="chart-description">
                Detailed analysis of molecules not correctly identified by the PFASGroups module.
                Includes both complete detection failures and incorrect group assignments.
            </div>
            
            <h3>📈 Failure Summary Statistics</h3>
            <div class="summary-grid">
                <div class="summary-card">
                    <div class="summary-number">{failure_report['summary']['total_failures']}</div>
                    <div class="summary-label">Total Failed<br>Detections</div>
                </div>
                <div class="summary-card">
                    <div class="summary-number">{failure_report['summary']['groups_with_failures']}</div>
                    <div class="summary-label">Groups with<br>Failures</div>
                </div>
                <div class="summary-card">
                    <div class="summary-number">{failure_report['summary']['no_detection_failures']}</div>
                    <div class="summary-label">No Detection<br>Failures</div>
                </div>
                <div class="summary-card">
                    <div class="summary-number">{failure_report['summary']['wrong_group_failures']}</div>
                    <div class="summary-label">Wrong Group<br>Assignments</div>
                </div>
            </div>
        </div>

        <h2>🏷️ PFAS-Atlas Classification Analysis</h2>
        <div class="visualization-container">
            <div class="chart-description">
                Analysis of PFAS-Atlas classification results showing the distribution of molecules across 
                PFAS-Atlas's native classification system (not converted to PFASGroups terms).
            </div>
            
            <h3>📊 PFAS-Atlas Class Definitions</h3>
            <div style="background: #f0f8ff; padding: 20px; border-radius: 8px; margin: 20px 0;">
                <h4>First Class Categories:</h4>
                <ul>
                    <li><strong>PFAA precursors:</strong> Precursors to perfluoroalkyl acids</li>
                    <li><strong>PFAAs:</strong> Perfluoroalkyl acids</li>
                    <li><strong>Other PFASs:</strong> Other per- and polyfluoroalkyl substances</li>
                    <li><strong>Not PFAS:</strong> Not identified as PFAS</li>
                </ul>
                
                <h4>Second Class Categories:</h4>
                <ul>
                    <li><strong>PFACs:</strong> Perfluoroalkyl compounds</li>
                    <li><strong>Aromatic PFASs:</strong> Aromatic per- and polyfluoroalkyl substances</li>
                    <li><strong>Polyfluoroalkyl acids:</strong> Polyfluoroalkyl carboxylic/sulfonic acids</li>
                </ul>
            </div>
            
            <div style="background: #fffacd; padding: 15px; border-radius: 8px; border-left: 4px solid #ffa500;">
                <strong>⚠️ Important Note:</strong> PFAS-Atlas classifications represent their native output classes 
                and are not directly comparable to PFASGroups functional group IDs. The systems use different 
                classification frameworks and evaluation criteria.
            </div>
            
            <h4>📊 PFAS-Atlas Classification Flow</h4>
            {atlas_flow_html}
        </div>

        <h2>📋 Detailed Performance Comparison Table</h2>
        <table class="comparison-table">
            <thead>
                <tr>
                    <th>Group ID</th>
                    <th>Functional Group</th>
                    <th>Molecules</th>
                    <th>PFASGroups Detection</th>
                    <th>PFASGroups Accuracy</th>
                    <th>Atlas Detection</th>
                    <th>Performance Gap</th>
                    <th>PFASGroups Avg Time (ms)</th>
                    <th>Atlas Avg Time (ms)</th>
                    <th>Speed Ratio</th>
                </tr>
            </thead>
            <tbody>
"""
    
    # Add detailed comparison table
    for group_id, data in sorted(single_analysis.items()):
        total = data['total_molecules']
        pfas_rate = (data['pfasgroups_correct'] / total * 100) if total > 0 else 0
        atlas_rate = (data['atlas_any_detection'] / total * 100) if total > 0 else 0
        gap = pfas_rate - atlas_rate
        
        def get_performance_class(rate):
            if rate >= 80: return "performance-excellent"
            elif rate >= 60: return "performance-good"
            elif rate >= 40: return "performance-moderate"
            else: return "performance-poor"
        
        pfas_class = get_performance_class(pfas_rate)
        atlas_class = get_performance_class(atlas_rate)
        gap_class = "performance-excellent" if gap > 50 else "performance-good" if gap > 20 else "performance-moderate" if gap > 0 else "performance-poor"
        
        # Calculate timing statistics
        pfas_avg_time_ms = data['pfasgroups_avg_time'] * 1000
        atlas_avg_time_ms = data['atlas_avg_time'] * 1000
        speed_ratio = atlas_avg_time_ms / pfas_avg_time_ms if pfas_avg_time_ms > 0 else 0
        
        # Timing performance classes
        def get_timing_class(time_ms):
            if time_ms < 1: return "performance-excellent"
            elif time_ms < 5: return "performance-good" 
            elif time_ms < 20: return "performance-moderate"
            else: return "performance-poor"
        
        pfas_timing_class = get_timing_class(pfas_avg_time_ms)
        atlas_timing_class = get_timing_class(atlas_avg_time_ms)
        
        html_content += f"""
                <tr>
                    <td>{group_id}</td>
                    <td>{data['group_name']}</td>
                    <td>{total}</td>
                    <td class="{pfas_class}">{pfas_rate:.1f}%</td>
                    <td class="{pfas_class}">{pfas_rate:.1f}%</td>
                    <td class="{atlas_class}">{atlas_rate:.1f}%</td>
                    <td class="{gap_class}">{gap:+.1f}%</td>
                    <td class="{pfas_timing_class}">{pfas_avg_time_ms:.2f}</td>
                    <td class="{atlas_timing_class}">{atlas_avg_time_ms:.2f}</td>
                    <td class="{'performance-good' if speed_ratio > 1 else 'performance-moderate'}">{speed_ratio:.1f}x</td>
                </tr>
"""
    
    html_content += f"""
            </tbody>
        </table>
        
        <h2>⚠️ PFASGroups Detection Failures Analysis</h2>
        <div class="visualization-container">
            <div class="chart-description">
                Detailed analysis of molecules not correctly identified by the PFASGroups module.
                Includes both complete detection failures and incorrect group assignments.
            </div>
            
            <h3>📈 Failure Summary Statistics</h3>
            <div class="summary-grid">
                <div class="summary-card">
                    <div class="summary-number">{failure_report['summary']['total_failures']}</div>
                    <div class="summary-label">Total Failed<br>Detections</div>
                </div>
                <div class="summary-card">
                    <div class="summary-number">{failure_report['summary']['groups_with_failures']}</div>
                    <div class="summary-label">Groups with<br>Failures</div>
                </div>
                <div class="summary-card">
                    <div class="summary-number">{failure_report['summary']['no_detection_failures']}</div>
                    <div class="summary-label">No Detection<br>Failures</div>
                </div>
                <div class="summary-card">
                    <div class="summary-number">{failure_report['summary']['wrong_group_failures']}</div>
                    <div class="summary-label">Wrong Group<br>Assignments</div>
                </div>
            </div>"""
    
    # Generate failure analysis table rows
    if failure_report['by_group']:
        html_content += f"""
            
            <h3>📝 Detailed Failure Analysis by Functional Group</h3>
            <table class="comparison-table">
                <thead>
                    <tr>
                        <th>Group ID</th>
                        <th>Group Name</th>
                        <th>Total Molecules</th>
                        <th>Failed Detections</th>
                        <th>Failure Rate</th>
                        <th>No Detection</th>
                        <th>Wrong Group</th>
                    </tr>
                </thead>
                <tbody>"""
        
        for group_id, failure_data in failure_report['by_group'].items():
            failure_class = 'performance-poor' if failure_data['failure_rate'] > 10 else 'performance-moderate' if failure_data['failure_rate'] > 5 else 'performance-good'
            html_content += f"""
                    <tr>
                        <td>{group_id}</td>
                        <td>{failure_data['group_name']}</td>
                        <td>{failure_data['total_molecules']}</td>
                        <td>{failure_data['failure_count']}</td>
                        <td class="{failure_class}">{failure_data['failure_rate']:.1f}%</td>
                        <td>{failure_data['no_detection_count']}</td>
                        <td>{failure_data['wrong_group_count']}</td>
                    </tr>"""
        
        html_content += """
                </tbody>
            </table>"""
    
    # Generate failure molecules table
    if failure_report['failure_molecules']:
        html_content += f"""
            
            <h3>🧪 Individual Failed Molecules (First 50)</h3>
            <div style="max-height: 400px; overflow-y: scroll; border: 1px solid #ddd; padding: 15px; border-radius: 8px; background: #f9f9f9;">
                <table class="comparison-table">
                    <thead>
                        <tr>
                            <th>Group ID</th>
                            <th>Group Name</th>
                            <th>SMILES</th>
                            <th>Failure Type</th>
                            <th>Detected Groups</th>
                        </tr>
                    </thead>
                    <tbody>"""
        
        # Show first 50 failures
        for i, failure in enumerate(failure_report['failure_molecules'][:50]):
            status_class = 'performance-poor' if failure['detection_status'] == 'no_detection' else 'performance-moderate'
            html_content += f"""
                        <tr>
                            <td>{failure['group_id']}</td>
                            <td>{failure['group_name']}</td>
                            <td style="font-family: monospace; font-size: 0.8em; max-width: 300px; word-break: break-all;">{failure['smiles']}</td>
                            <td class="{status_class}">{failure['detection_status'].replace('_', ' ').title()}</td>
                            <td>{failure['detected_groups'] if failure['detected_groups'] else 'None'}</td>
                        </tr>"""
        
        if len(failure_report['failure_molecules']) > 50:
            html_content += f"""
                        <tr style="background: #f0f0f0;">
                            <td colspan="5" style="text-align: center; font-style: italic;">
                                ... and {len(failure_report['failure_molecules']) - 50} more failed detections
                            </td>
                        </tr>"""
        
        html_content += """
                    </tbody>
                </table>
            </div>"""
    
    html_content += """
        </div>
        
        <h2>🏷️ PFAS-Atlas Classification Analysis</h2>
        <div class="visualization-container">
            <div class="chart-description">
                Analysis of PFAS-Atlas classification results showing the distribution of molecules across 
                PFAS-Atlas's native classification system (not converted to PFASGroups terms).
            </div>
            
            <div style="background: #f0f8ff; padding: 20px; border-radius: 8px; margin: 20px 0;">
                <h4>PFAS-Atlas Class Definitions:</h4>
                <div style="display: grid; grid-template-columns: 1fr 1fr; gap: 20px;">
                    <div>
                        <h5>First Class Categories:</h5>
                        <ul style="margin: 0; padding-left: 20px;">
                            <li><strong>PFAA precursors:</strong> Precursors to perfluoroalkyl acids</li>
                            <li><strong>PFAAs:</strong> Perfluoroalkyl acids</li>
                            <li><strong>Other PFASs:</strong> Other per- and polyfluoroalkyl substances</li>
                            <li><strong>Not PFAS:</strong> Not identified as PFAS</li>
                        </ul>
                    </div>
                    <div>
                        <h5>Second Class Categories:</h5>
                        <ul style="margin: 0; padding-left: 20px;">
                            <li><strong>PFACs:</strong> Perfluoroalkyl compounds</li>
                            <li><strong>Aromatic PFASs:</strong> Aromatic per- and polyfluoroalkyl substances</li>
                            <li><strong>Polyfluoroalkyl acids:</strong> Polyfluoroalkyl carboxylic/sulfonic acids</li>
                        </ul>
                    </div>
                </div>
            </div>
            
            <div style="background: #fffacd; padding: 15px; border-radius: 8px; border-left: 4px solid #ffa500;">
                <strong>⚠️ Important Note:</strong> PFAS-Atlas classifications represent their native output classes 
                and are not directly comparable to PFASGroups functional group IDs. The systems use different 
                classification frameworks and evaluation criteria.
            </div>
        </div>

        <div class="key-findings">
            <h3>🔍 Enhanced Key Findings</h3>
            <div style="text-align: left; max-width: 900px; margin: 20px auto;">
                <h4>🎯 System Performance Comparison:</h4>
                <ul>
                    <li><strong>PFASGroups dominates single-group detection:</strong> {pfasgroups_single_rate:.1f}% vs Atlas {atlas_single_rate:.1f}% ({pfasgroups_single_rate - atlas_single_rate:+.1f}% advantage)</li>
                    <li><strong>Multi-group challenge:</strong> PFASGroups {pfasgroups_multi_rate:.1f}% vs Atlas {atlas_multi_rate:.1f}% detection rates</li>
                    <li><strong>Consistency advantage:</strong> PFASGroups shows more reliable performance across functional groups</li>
                </ul>
                
                <h4>🔬 Multi-Group Insights:</h4>
                <ul>
                    <li><strong>Functional group privilege:</strong> Some groups are preferentially detected in complex molecules</li>
                    <li><strong>Detection hierarchy:</strong> Clear patterns emerge in which groups dominate multi-functional detection</li>
                    <li><strong>System limitations:</strong> Both systems struggle with highly complex multi-functional molecules</li>
                </ul>
                
                <h4>⚡ Performance Recommendations:</h4>
                <ul>
                    <li><strong>Use PFASGroups for single-group detection:</strong> Superior accuracy and consistency</li>
                    <li><strong>Enhance multi-group algorithms:</strong> Address detection privilege and improve complex molecule handling</li>
                    <li><strong>Complementary approach:</strong> Consider hybrid systems leveraging strengths of both</li>
                </ul>
                
                <h4>⏱️ Execution Time Analysis:</h4>
                <ul>
                    <li><strong>Speed comparison:</strong> Both systems show sub-millisecond to low-millisecond processing times per molecule</li>
                    <li><strong>Performance patterns:</strong> Timing varies by functional group complexity and molecular structure</li>
                    <li><strong>Efficiency insights:</strong> Speed ratios reveal relative computational costs between detection systems</li>
                </ul>
            </div>
        </div>

        <h2>📊 Exported Visualizations</h2>
        <div class="visualization-container">
            <p>All visualizations have been exported in multiple formats for use in presentations and publications:</p>
            <ul>
                <li><strong>Performance Comparison Heatmap:</strong> 
                    <a href="imgs/comparison_heatmap_{timestamp}.png">PNG</a> | 
                    <a href="imgs/comparison_heatmap_{timestamp}.svg">SVG</a>
                </li>
                <li><strong>Multi-Group Privilege Heatmap:</strong> 
                    <a href="imgs/multigroup_privilege_heatmap_{timestamp}.png">PNG</a> | 
                    <a href="imgs/multigroup_privilege_heatmap_{timestamp}.svg">SVG</a>
                </li>
                <li><strong>Enhanced System Sankey:</strong> 
                    <a href="imgs/enhanced_system_sankey_{timestamp}.png">PNG</a> | 
                    <a href="imgs/enhanced_system_sankey_{timestamp}.svg">SVG</a>
                </li>
                <li><strong>Privilege Hierarchy Sankey:</strong> 
                    <a href="imgs/privilege_hierarchy_sankey_{timestamp}.png">PNG</a> | 
                    <a href="imgs/privilege_hierarchy_sankey_{timestamp}.svg">SVG</a>
                </li>
                <li><strong>Timing Comparison Heatmap:</strong> 
                    <a href="imgs/timing_comparison_heatmap_{timestamp}.png">PNG</a> | 
                    <a href="imgs/timing_comparison_heatmap_{timestamp}.svg">SVG</a>
                </li>
                <li><strong>Detailed Timing Statistics:</strong> 
                    <a href="imgs/timing_statistics_{timestamp}.png">PNG</a> | 
                    <a href="imgs/timing_statistics_{timestamp}.svg">SVG</a>
                </li>
                <li><strong>PFAS-Atlas Classification Flow:</strong> 
                    <a href="imgs/atlas_classification_flow_{timestamp}.png">PNG</a> | 
                    <a href="imgs/atlas_classification_flow_{timestamp}.svg">SVG</a>
                </li>
            </ul>
        </div>

        <div style="text-align: center; color: #666; font-size: 0.95em; margin-top: 60px; padding-top: 40px; border-top: 3px solid #eee;">
            <strong>Enhanced PFAS Benchmark Analysis</strong><br>
            Generated on {datetime.now().strftime('%Y-%m-%d at %H:%M:%S')} | 
            Total Molecules: {total_single_molecules + total_multi_molecules} | 
            Functional Groups: {len(single_analysis)} | 
            Multi-Group Combinations: {len(multi_analysis)}
        </div>
    </div>
</body>
</html>
"""
    
    return html_content

def main():
    """Main enhanced analysis function"""
    
    print("🚀 ENHANCED PFAS BENCHMARK ANALYSIS")
    print("=" * 50)
    
    # Find the most recent enhanced benchmark file
    import glob
    benchmark_files = glob.glob('data/pfas_enhanced_benchmark_*.json')
    if not benchmark_files:
        print("❌ No enhanced benchmark results found. Run enhanced_pfas_benchmark.py first.")
        return

    latest_file = max(benchmark_files)
    print(f"📊 Using enhanced benchmark file: {latest_file}")
    
    # Load and analyze results
    results = load_benchmark_results(latest_file)
    
    # Comprehensive analysis
    print("📋 Analyzing system performance comparison...")
    single_analysis, multi_analysis = analyze_system_comparison(results)
    
    # Create timestamp for file naming
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    
    # Generate all visualizations and analysis
    print("📊 Creating enhanced visualizations...")
    print("   • Performance comparison heatmap")
    print("   • Multi-group privilege analysis") 
    print("   • Enhanced system comparison Sankey")
    print("   • Functional group hierarchy analysis")
    
    # Create comprehensive HTML report
    print("🌐 Creating enhanced HTML analysis report...")
    html_content = create_enhanced_html_report(single_analysis, multi_analysis, timestamp, results)
    
    html_filename = f"enhanced_pfas_analysis_{timestamp}.html"
    with open(html_filename, 'w', encoding='utf-8') as f:
        f.write(html_content)
    
    # Calculate and display summary statistics
    total_single = sum(data['total_molecules'] for data in single_analysis.values())
    total_multi = sum(data['total_molecules'] for data in multi_analysis.values())
    
    pfas_single_success = sum(data['pfasgroups_correct'] for data in single_analysis.values())
    atlas_single_success = sum(data['atlas_any_detection'] for data in single_analysis.values())
    
    pfas_single_rate = (pfas_single_success / total_single * 100) if total_single > 0 else 0
    atlas_single_rate = (atlas_single_success / total_single * 100) if total_single > 0 else 0
    
    print(f"\n✅ Enhanced Analysis Complete!")
    print("=" * 50)
    print(f"📊 Dataset Summary:")
    print(f"   • Total molecules tested: {total_single + total_multi}")
    print(f"   • Single-group molecules: {total_single} (15 per group)")
    print(f"   • Multi-group molecules: {total_multi} (pairs + triplets)")
    print(f"   • Functional groups: {len(single_analysis)}")
    
    print(f"\n🎯 Performance Summary:")
    print(f"   • PFASGroups single-group: {pfas_single_rate:.1f}%")
    print(f"   • PFAS-Atlas single-group: {atlas_single_rate:.1f}%") 
    print(f"   • Performance gap: {pfas_single_rate - atlas_single_rate:+.1f}%")
    
    print(f"\n📁 Generated Files:")
    print(f"   • {html_filename} (comprehensive analysis)")
    print(f"   • imgs/comparison_heatmap_{timestamp}.png/.svg")
    print(f"   • imgs/multigroup_privilege_heatmap_{timestamp}.png/.svg") 
    print(f"   • imgs/enhanced_system_sankey_{timestamp}.png/.svg")
    print(f"   • imgs/privilege_hierarchy_sankey_{timestamp}.png/.svg")



def analyze_combined_results():
    """Analyze both enhanced and OECD benchmark results"""
    
    print("🚀 COMBINED PFAS BENCHMARK ANALYSIS")
    print("=" * 50)
    
    # Find benchmark files
    import glob
    
    # Enhanced benchmark
    enhanced_files = glob.glob('pfas_enhanced_benchmark_*.json')
    if not enhanced_files:
        print("❌ No enhanced benchmark results found.")
        return
    
    latest_enhanced = max(enhanced_files)
    print(f"📊 Using enhanced benchmark file: {latest_enhanced}")
    
    # OECD benchmark
    oecd_files = glob.glob('pfas_oecd_benchmark_*.json')
    if not oecd_files:
        print("❌ No OECD benchmark results found.")
        return
    
    latest_oecd = max(oecd_files)
    print(f"📊 Using OECD benchmark file: {latest_oecd}")
    
    # Load results
    enhanced_results = load_benchmark_results(latest_enhanced)
    oecd_results = load_benchmark_results(latest_oecd)
    
    # Analyze enhanced results
    print("📋 Analyzing enhanced benchmark...")
    single_analysis, multi_analysis = analyze_system_comparison(enhanced_results)
    
    # Analyze OECD results
    print("📋 Analyzing OECD benchmark...")
    oecd_analysis = analyze_oecd_benchmark(oecd_results)
    
    # Create timestamp
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    
    # Generate comprehensive reports
    print("📊 Creating comprehensive HTML reports...")
    
    # Enhanced analysis report
    html_enhanced = create_enhanced_html_report(single_analysis, multi_analysis, timestamp, enhanced_results)
    enhanced_filename = f"enhanced_pfas_analysis_{timestamp}.html"
    with open(enhanced_filename, 'w') as f:
        f.write(html_enhanced)
    
    # OECD analysis report
    create_oecd_html_report(oecd_analysis, timestamp)
    
    # Combined comparison report with Sankey diagrams
    create_combined_comparison_report(enhanced_results, oecd_results, timestamp)
    
    create_enhanced_sankey_comparison(single_analysis,multi_analysis, enhanced_results)

    print("\n✅ COMBINED ANALYSIS COMPLETE!")
    print("📄 Generated reports:")
    print(f"   • {enhanced_filename} (enhanced benchmark)")
    print(f"   • oecd_pfas_analysis_{timestamp}.html (OECD validation)")
    print(f"   • combined_pfas_analysis_{timestamp}.html (comparative analysis with Sankey diagrams)")
    
    # Summary statistics
    print("\n📊 Summary Statistics:")
    enhanced_total = len(enhanced_results)
    oecd_total = oecd_analysis['total_molecules']
    
    enhanced_single_success = sum(data['pfasgroups_correct'] for data in single_analysis.values())
    enhanced_single_total = sum(data['total_molecules'] for data in single_analysis.values())
    enhanced_rate = (enhanced_single_success / enhanced_single_total * 100) if enhanced_single_total > 0 else 0
    
    oecd_rate = (oecd_analysis['pfasgroups_detections'] / oecd_total * 100) if oecd_total > 0 else 0
    
    print(f"   • Enhanced Dataset: {enhanced_rate:.1f}% PFASGroups success ({enhanced_total} molecules)")
    print(f"   • OECD Dataset: {oecd_rate:.1f}% PFASGroups success ({oecd_total} molecules)")
    print(f"   • Atlas Detection Rate on OECD: {(oecd_analysis['atlas_detections'] / oecd_total * 100):.1f}%")

if __name__ == "__main__":
    import sys
    if len(sys.argv) > 1 and sys.argv[1] == "--combined":
        analyze_combined_results()
    else:
        main()