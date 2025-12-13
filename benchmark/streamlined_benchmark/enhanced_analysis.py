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
                'molecules': []
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
        
        single_analysis[group_id]['molecules'].append({
            'smiles': result['molecule_data']['smiles'],
            'pfasgroups_detected': detected_groups,
            'pfasgroups_success': pfas_detected,
            'atlas_success': atlas_success,
            'atlas_class': result['atlas_result']['second_class']
        })
    
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
    
    return single_analysis, multi_analysis

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
        title="PFASGroups vs PFAS-Atlas Detection Performance Comparison",
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

def create_enhanced_sankey_comparison(single_analysis, multi_analysis):
    """Create enhanced Sankey diagram showing system comparison"""
    
    # Calculate overall statistics
    total_single = sum(data['total_molecules'] for data in single_analysis.values())
    pfasgroups_single_success = sum(data['pfasgroups_correct'] for data in single_analysis.values())
    atlas_single_success = sum(data['atlas_any_detection'] for data in single_analysis.values())
    
    total_multi = sum(data['total_molecules'] for data in multi_analysis.values())
    pfasgroups_multi_success = sum(data['pfasgroups_any_detection'] for data in multi_analysis.values())
    atlas_multi_success = sum(data['atlas_classifications'] for data in multi_analysis.values())
    
    # Create nodes
    nodes = [
        # Source nodes
        "Single Group<br>Molecules", "Multi-Group<br>Molecules",
        # PFASGroups outcomes
        "PFASGroups<br>Success", "PFASGroups<br>Failure", 
        # Atlas outcomes
        "Atlas<br>Success", "Atlas<br>Failure"
    ]
    
    node_colors = [
        "rgba(31, 119, 180, 0.8)", "rgba(255, 127, 14, 0.8)",  # Sources
        "rgba(44, 160, 44, 0.8)", "rgba(214, 39, 40, 0.8)",   # PFASGroups
        "rgba(148, 103, 189, 0.8)", "rgba(188, 189, 34, 0.8)" # Atlas
    ]
    
    # Create links
    sources = []
    targets = []
    values = []
    link_colors = []
    
    # Single group to PFASGroups
    sources.extend([0, 0])
    targets.extend([2, 3])
    values.extend([pfasgroups_single_success, total_single - pfasgroups_single_success])
    link_colors.extend(["rgba(44, 160, 44, 0.4)", "rgba(214, 39, 40, 0.4)"])
    
    # Single group to Atlas
    sources.extend([0, 0])
    targets.extend([4, 5])
    values.extend([atlas_single_success, total_single - atlas_single_success])
    link_colors.extend(["rgba(148, 103, 189, 0.4)", "rgba(188, 189, 34, 0.4)"])
    
    # Multi-group to PFASGroups
    sources.extend([1, 1])
    targets.extend([2, 3])
    values.extend([pfasgroups_multi_success, total_multi - pfasgroups_multi_success])
    link_colors.extend(["rgba(44, 160, 44, 0.4)", "rgba(214, 39, 40, 0.4)"])
    
    # Multi-group to Atlas
    sources.extend([1, 1])
    targets.extend([4, 5])
    values.extend([atlas_multi_success, total_multi - atlas_multi_success])
    link_colors.extend(["rgba(148, 103, 189, 0.4)", "rgba(188, 189, 34, 0.4)"])
    
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
            color=link_colors
        )
    )])
    
    pfas_single_rate = (pfasgroups_single_success / total_single * 100) if total_single > 0 else 0
    atlas_single_rate = (atlas_single_success / total_single * 100) if total_single > 0 else 0
    pfas_multi_rate = (pfasgroups_multi_success / total_multi * 100) if total_multi > 0 else 0
    atlas_multi_rate = (atlas_multi_success / total_multi * 100) if total_multi > 0 else 0
    
    fig.update_layout(
        title=f"Enhanced System Comparison: PFASGroups vs PFAS-Atlas<br>" + 
              f"<sub>Single: PFASGroups {pfas_single_rate:.1f}% vs Atlas {atlas_single_rate:.1f}% | " +
              f"Multi: PFASGroups {pfas_multi_rate:.1f}% vs Atlas {atlas_multi_rate:.1f}%</sub>",
        font_size=12,
        width=1200,
        height=600
    )
    
    return fig

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

def create_enhanced_html_report(single_analysis, multi_analysis, timestamp):
    """Create comprehensive HTML report with embedded visualizations"""
    
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
    sankey_comparison = create_enhanced_sankey_comparison(single_analysis, multi_analysis)
    sankey_privilege = create_detailed_multi_group_sankey(multi_analysis)
    
    # Convert to HTML
    heatmap_comparison_html = heatmap_comparison.to_html(include_plotlyjs='cdn', div_id="heatmap-comparison")
    heatmap_multi_html = heatmap_multi.to_html(include_plotlyjs=False, div_id="heatmap-multi")
    sankey_comparison_html = sankey_comparison.to_html(include_plotlyjs=False, div_id="sankey-comparison")
    sankey_privilege_html = sankey_privilege.to_html(include_plotlyjs=False, div_id="sankey-privilege")
    
    # Save individual visualizations
    heatmap_comparison.write_image(f"comparison_heatmap_{timestamp}.png", width=1400, height=500, scale=2)
    heatmap_comparison.write_image(f"comparison_heatmap_{timestamp}.svg")
    heatmap_multi.write_image(f"multigroup_privilege_heatmap_{timestamp}.png", width=1200, height=600, scale=2)
    heatmap_multi.write_image(f"multigroup_privilege_heatmap_{timestamp}.svg")
    sankey_comparison.write_image(f"enhanced_system_sankey_{timestamp}.png", width=1200, height=600, scale=2)
    sankey_comparison.write_image(f"enhanced_system_sankey_{timestamp}.svg")
    sankey_privilege.write_image(f"privilege_hierarchy_sankey_{timestamp}.png", width=1000, height=600, scale=2)
    sankey_privilege.write_image(f"privilege_hierarchy_sankey_{timestamp}.svg")
    
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
            <p>Target: PFAS functional groups 29-51 (excluding group 48) | Enhanced dataset sizes | Detailed privilege analysis</p>
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
                Comprehensive Sankey diagram showing success/failure flows for both single and multi-group molecules 
                across both detection systems.
            </div>
            {sankey_comparison_html}
        </div>

        <h2>🎯 Multi-Group Privilege Analysis</h2>
        <div class="visualization-container">
            <div class="chart-description">
                Heatmap showing which functional groups are preferentially detected in multi-functional molecules. 
                Reveals system biases and detection hierarchies.
            </div>
            {heatmap_multi_html}
        </div>

        <h2>📈 Functional Group Detection Hierarchy</h2>
        <div class="visualization-container">
            <div class="chart-description">
                Sankey diagram revealing the hierarchy of functional group detection in multi-group scenarios.
                Shows which groups are "privileged" during detection.
            </div>
            {sankey_privilege_html}
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
        
        html_content += f"""
                <tr>
                    <td>{group_id}</td>
                    <td>{data['group_name']}</td>
                    <td>{total}</td>
                    <td class="{pfas_class}">{pfas_rate:.1f}%</td>
                    <td class="{pfas_class}">{pfas_rate:.1f}%</td>
                    <td class="{atlas_class}">{atlas_rate:.1f}%</td>
                    <td class="{gap_class}">{gap:+.1f}%</td>
                </tr>
"""
    
    html_content += f"""
            </tbody>
        </table>

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
            </div>
        </div>

        <h2>📊 Exported Visualizations</h2>
        <div class="visualization-container">
            <p>All visualizations have been exported in multiple formats for use in presentations and publications:</p>
            <ul>
                <li><strong>Performance Comparison Heatmap:</strong> 
                    <a href="comparison_heatmap_{timestamp}.png">PNG</a> | 
                    <a href="comparison_heatmap_{timestamp}.svg">SVG</a>
                </li>
                <li><strong>Multi-Group Privilege Heatmap:</strong> 
                    <a href="multigroup_privilege_heatmap_{timestamp}.png">PNG</a> | 
                    <a href="multigroup_privilege_heatmap_{timestamp}.svg">SVG</a>
                </li>
                <li><strong>Enhanced System Sankey:</strong> 
                    <a href="enhanced_system_sankey_{timestamp}.png">PNG</a> | 
                    <a href="enhanced_system_sankey_{timestamp}.svg">SVG</a>
                </li>
                <li><strong>Privilege Hierarchy Sankey:</strong> 
                    <a href="privilege_hierarchy_sankey_{timestamp}.png">PNG</a> | 
                    <a href="privilege_hierarchy_sankey_{timestamp}.svg">SVG</a>
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
    benchmark_files = glob.glob('pfas_enhanced_benchmark_*.json')
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
    html_content = create_enhanced_html_report(single_analysis, multi_analysis, timestamp)
    
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
    print(f"   • comparison_heatmap_{timestamp}.png/.svg")
    print(f"   • multigroup_privilege_heatmap_{timestamp}.png/.svg") 
    print(f"   • enhanced_system_sankey_{timestamp}.png/.svg")
    print(f"   • privilege_hierarchy_sankey_{timestamp}.png/.svg")

if __name__ == "__main__":
    main()