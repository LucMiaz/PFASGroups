"""
Analysis and Visualization for PFAS Benchmark Results
Creates Sankey diagrams and HTML analysis
"""

import json
import pandas as pd
import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from collections import defaultdict, Counter
from datetime import datetime

def load_benchmark_results(filename):
    """Load benchmark results from JSON"""
    
    with open(filename, 'r') as f:
        results = json.load(f)
    
    print(f"📊 Loaded {len(results)} benchmark results")
    return results

def analyze_single_group_performance(results):
    """Analyze performance for single group molecules"""
    
    single_group_results = [r for r in results if r['molecule_data'].get('generation_type') == 'single_group']
    
    analysis = {}
    
    for result in single_group_results:
        group_id = result['molecule_data']['group_id']
        group_name = result['molecule_data']['group_name']
        
        if group_id not in analysis:
            analysis[group_id] = {
                'group_name': group_name,
                'total_molecules': 0,
                'pfasgroups_detections': 0,
                'atlas_classifications': 0,
                'correct_pfasgroups_detections': 0,
                'molecules': []
            }
        
        analysis[group_id]['total_molecules'] += 1
        
        # Check PFASGroups detection
        pfas_detected = result['pfasgroups_result']['success']
        if pfas_detected:
            analysis[group_id]['pfasgroups_detections'] += 1
            
            # Check if correct group was detected
            detected_groups = result['pfasgroups_result']['detected_groups']
            if group_id in detected_groups:
                analysis[group_id]['correct_pfasgroups_detections'] += 1
        
        # Check PFAS-Atlas classification
        atlas_success = result['atlas_result']['success']
        if atlas_success:
            analysis[group_id]['atlas_classifications'] += 1
        
        analysis[group_id]['molecules'].append({
            'smiles': result['molecule_data']['smiles'],
            'pfasgroups_success': pfas_detected,
            'pfasgroups_correct': group_id in result['pfasgroups_result']['detected_groups'],
            'atlas_success': atlas_success,
            'atlas_class': result['atlas_result']['second_class']
        })
    
    return analysis

def analyze_multi_group_performance(results):
    """Analyze performance for multi-group molecules"""
    
    multi_group_results = [r for r in results if r['molecule_data'].get('generation_type') == 'multi_group']
    
    analysis = {}
    
    for result in multi_group_results:
        target_groups = tuple(sorted(result['molecule_data']['target_groups']))
        
        if target_groups not in analysis:
            analysis[target_groups] = {
                'target_groups': list(target_groups),
                'total_molecules': 0,
                'pfasgroups_detections': 0,
                'atlas_classifications': 0,
                'multi_group_detections': 0,
                'molecules': []
            }
        
        analysis[target_groups]['total_molecules'] += 1
        
        # Check PFASGroups detection
        detected_groups = result['pfasgroups_result']['detected_groups']
        pfas_detected = len(detected_groups) > 0
        
        if pfas_detected:
            analysis[target_groups]['pfasgroups_detections'] += 1
            
            # Check if multiple groups detected
            target_set = set(target_groups)
            detected_set = set(detected_groups)
            overlap = len(target_set.intersection(detected_set))
            
            if overlap >= 2:
                analysis[target_groups]['multi_group_detections'] += 1
        
        # Check PFAS-Atlas classification
        atlas_success = result['atlas_result']['success']
        if atlas_success:
            analysis[target_groups]['atlas_classifications'] += 1
        
        analysis[target_groups]['molecules'].append({
            'smiles': result['molecule_data']['smiles'],
            'target_groups': list(target_groups),
            'detected_groups': detected_groups,
            'atlas_success': atlas_success,
            'atlas_class': result['atlas_result']['second_class']
        })
    
    return analysis

def create_single_group_sankey(single_group_analysis):
    """Create Sankey diagram for single group performance"""
    
    # Prepare data
    groups = list(single_group_analysis.keys())
    group_names = [single_group_analysis[g]['group_name'] for g in groups]
    
    # Node categories: Groups -> Detection -> Classification
    detection_outcomes = ['PFASGroups Detected', 'PFASGroups Not Detected']
    atlas_outcomes = ['Atlas Classified', 'Atlas Not Classified']
    
    all_nodes = group_names + detection_outcomes + atlas_outcomes
    node_dict = {node: idx for idx, node in enumerate(all_nodes)}
    
    sources = []
    targets = []
    values = []
    
    for group_id in groups:
        data = single_group_analysis[group_id]
        group_name = data['group_name']
        
        total = data['total_molecules']
        detected = data['pfasgroups_detections']
        not_detected = total - detected
        classified = data['atlas_classifications']
        not_classified = total - classified
        
        # Group -> PFASGroups Detection
        if detected > 0:
            sources.append(node_dict[group_name])
            targets.append(node_dict['PFASGroups Detected'])
            values.append(detected)
        
        if not_detected > 0:
            sources.append(node_dict[group_name])
            targets.append(node_dict['PFASGroups Not Detected'])
            values.append(not_detected)
        
        # Group -> Atlas Classification
        if classified > 0:
            sources.append(node_dict[group_name])
            targets.append(node_dict['Atlas Classified'])
            values.append(classified)
        
        if not_classified > 0:
            sources.append(node_dict[group_name])
            targets.append(node_dict['Atlas Not Classified'])
            values.append(not_classified)
    
    # Create colors
    node_colors = ['rgba(31, 119, 180, 0.8)'] * len(group_names) + \
                  ['rgba(44, 160, 44, 0.8)', 'rgba(214, 39, 40, 0.8)'] + \
                  ['rgba(255, 127, 14, 0.8)', 'rgba(148, 103, 189, 0.8)']
    
    fig = go.Figure(data=[go.Sankey(
        node=dict(
            pad=15,
            thickness=20,
            line=dict(color="black", width=0.5),
            label=all_nodes,
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
        title="Single PFAS Group Detection Performance",
        font_size=10,
        width=1200,
        height=700
    )
    
    return fig

def create_multi_group_sankey(multi_group_analysis):
    """Create Sankey diagram for multi-group performance"""
    
    # Simplify for visualization
    combo_names = []
    for combo in multi_group_analysis.keys():
        combo_name = f"Groups {'+'.join([str(g) for g in combo])}"
        combo_names.append(combo_name)
    
    outcomes = ['Multi-Group Detected', 'Single Group Detected', 'No Detection', 'Atlas Classified']
    
    all_nodes = combo_names + outcomes
    node_dict = {node: idx for idx, node in enumerate(all_nodes)}
    
    sources = []
    targets = []
    values = []
    
    for i, (combo, data) in enumerate(multi_group_analysis.items()):
        combo_name = combo_names[i]
        
        total = data['total_molecules'] 
        multi_detected = data['multi_group_detections']
        single_detected = data['pfasgroups_detections'] - multi_detected
        no_detection = total - data['pfasgroups_detections']
        atlas_classified = data['atlas_classifications']
        
        # Combo -> Detection outcomes
        if multi_detected > 0:
            sources.append(node_dict[combo_name])
            targets.append(node_dict['Multi-Group Detected'])
            values.append(multi_detected)
        
        if single_detected > 0:
            sources.append(node_dict[combo_name])
            targets.append(node_dict['Single Group Detected'])
            values.append(single_detected)
        
        if no_detection > 0:
            sources.append(node_dict[combo_name])
            targets.append(node_dict['No Detection'])
            values.append(no_detection)
        
        # Combo -> Atlas
        if atlas_classified > 0:
            sources.append(node_dict[combo_name])
            targets.append(node_dict['Atlas Classified'])
            values.append(atlas_classified)
    
    node_colors = ['rgba(31, 119, 180, 0.8)'] * len(combo_names) + \
                  ['rgba(44, 160, 44, 0.8)', 'rgba(255, 193, 7, 0.8)', 
                   'rgba(214, 39, 40, 0.8)', 'rgba(255, 127, 14, 0.8)']
    
    fig = go.Figure(data=[go.Sankey(
        node=dict(
            pad=15,
            thickness=20,
            line=dict(color="black", width=0.5),
            label=all_nodes,
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
        title="Multi-Group PFAS Detection Performance",
        font_size=10,
        width=1200,
        height=600
    )
    
    return fig

def create_html_analysis(single_group_analysis, multi_group_analysis, timestamp):
    """Create comprehensive HTML analysis report"""
    
    # Calculate overall statistics
    total_single_molecules = sum(data['total_molecules'] for data in single_group_analysis.values())
    total_multi_molecules = sum(data['total_molecules'] for data in multi_group_analysis.values())
    
    single_pfas_success = sum(data['pfasgroups_detections'] for data in single_group_analysis.values())
    single_atlas_success = sum(data['atlas_classifications'] for data in single_group_analysis.values())
    
    multi_pfas_success = sum(data['pfasgroups_detections'] for data in multi_group_analysis.values())
    multi_atlas_success = sum(data['atlas_classifications'] for data in multi_group_analysis.values())
    
    html_content = f"""
<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>PFAS Benchmark Analysis - Comprehensive Results</title>
    <style>
        body {{
            font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif;
            margin: 0;
            padding: 20px;
            background: linear-gradient(135deg, #2c3e50 0%, #3498db 100%);
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
            font-size: 2.5em;
            margin-bottom: 10px;
        }}
        h2 {{
            color: #34495e;
            border-left: 4px solid #3498db;
            padding-left: 15px;
            margin-top: 40px;
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
            padding: 20px;
            border-radius: 10px;
            text-align: center;
        }}
        .summary-number {{
            font-size: 2.5em;
            font-weight: bold;
        }}
        .summary-label {{
            font-size: 0.9em;
            opacity: 0.9;
        }}
        table {{
            width: 100%;
            border-collapse: collapse;
            margin: 20px 0;
            font-size: 0.9em;
        }}
        th, td {{
            border: 1px solid #ddd;
            padding: 10px;
            text-align: left;
        }}
        th {{
            background-color: #f8f9fa;
            font-weight: bold;
            position: sticky;
            top: 0;
        }}
        tr:nth-child(even) {{
            background-color: #f9f9f9;
        }}
        .success-rate {{
            font-weight: bold;
        }}
        .excellent {{ color: #28a745; }}
        .good {{ color: #17a2b8; }}
        .moderate {{ color: #ffc107; }}
        .poor {{ color: #dc3545; }}
        .highlight-section {{
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            color: white;
            padding: 30px;
            border-radius: 12px;
            margin: 40px 0;
            text-align: center;
        }}
        .section-container {{
            background-color: #f8f9fa;
            padding: 25px;
            border-radius: 10px;
            margin: 30px 0;
        }}
        .molecular-example {{
            font-family: monospace;
            background-color: #e9ecef;
            padding: 5px;
            border-radius: 3px;
            font-size: 0.8em;
        }}
        .performance-indicator {{
            display: inline-block;
            padding: 3px 8px;
            border-radius: 12px;
            font-size: 0.8em;
            margin-left: 10px;
        }}
    </style>
</head>
<body>
    <div class="container">
        <h1>🧪 PFAS Detection Benchmark Results</h1>
        
        <div class="highlight-section">
            <h3>Comprehensive PFASGroups vs PFAS-Atlas Comparison</h3>
            <p>Analysis Date: {datetime.now().strftime('%B %d, %Y at %H:%M')} | 
            Generated Molecules: {total_single_molecules + total_multi_molecules} | 
            Target Groups: 29-51 (excluding 48)</p>
        </div>

        <div class="summary-grid">
            <div class="summary-card">
                <div class="summary-number">{len(single_group_analysis)}</div>
                <div class="summary-label">Single Groups Tested</div>
            </div>
            <div class="summary-card">
                <div class="summary-number">{len(multi_group_analysis)}</div>
                <div class="summary-label">Multi-Group Combinations</div>
            </div>
            <div class="summary-card">
                <div class="summary-number">{single_pfas_success/total_single_molecules*100:.0f}%</div>
                <div class="summary-label">Single Group PFASGroups Success</div>
            </div>
            <div class="summary-card">
                <div class="summary-number">{multi_pfas_success/total_multi_molecules*100:.0f}%</div>
                <div class="summary-label">Multi-Group PFASGroups Success</div>
            </div>
        </div>

        <h2>📊 Single Group Performance Analysis</h2>
        <div class="section-container">
            <table>
                <thead>
                    <tr>
                        <th>Group ID</th>
                        <th>Group Name</th>
                        <th>Molecules</th>
                        <th>PFASGroups Detection</th>
                        <th>Correct Detection</th>
                        <th>Atlas Classification</th>
                        <th>Example Molecule</th>
                    </tr>
                </thead>
                <tbody>
"""
    
    # Add single group results
    for group_id, data in sorted(single_group_analysis.items()):
        total = data['total_molecules']
        pfas_detected = data['pfasgroups_detections'] 
        correct_detected = data['correct_pfasgroups_detections']
        atlas_classified = data['atlas_classifications']
        
        pfas_rate = pfas_detected / total * 100 if total > 0 else 0
        correct_rate = correct_detected / total * 100 if total > 0 else 0
        atlas_rate = atlas_classified / total * 100 if total > 0 else 0
        
        # Get performance indicators
        def get_performance_class(rate):
            if rate >= 80:
                return "excellent"
            elif rate >= 60:
                return "good" 
            elif rate >= 40:
                return "moderate"
            else:
                return "poor"
        
        pfas_class = get_performance_class(pfas_rate)
        atlas_class = get_performance_class(atlas_rate)
        
        example_mol = data['molecules'][0]['smiles'] if data['molecules'] else 'N/A'
        
        html_content += f"""
                    <tr>
                        <td>{group_id}</td>
                        <td>{data['group_name']}</td>
                        <td>{total}</td>
                        <td class="success-rate {pfas_class}">{pfas_rate:.1f}%</td>
                        <td class="success-rate {pfas_class}">{correct_rate:.1f}%</td>
                        <td class="success-rate {atlas_class}">{atlas_rate:.1f}%</td>
                        <td class="molecular-example">{example_mol[:50]}{'...' if len(example_mol) > 50 else ''}</td>
                    </tr>
"""
    
    html_content += """
                </tbody>
            </table>
        </div>

        <h2>🔬 Multi-Group Performance Analysis</h2>
        <div class="section-container">
            <table>
                <thead>
                    <tr>
                        <th>Target Groups</th>
                        <th>Molecules</th>
                        <th>PFASGroups Detection</th>
                        <th>Multi-Group Detection</th>
                        <th>Atlas Classification</th>
                        <th>Example Molecule</th>
                    </tr>
                </thead>
                <tbody>
"""
    
    # Add multi-group results
    for combo, data in multi_group_analysis.items():
        combo_str = '+'.join([str(g) for g in combo])
        total = data['total_molecules']
        pfas_detected = data['pfasgroups_detections']
        multi_detected = data['multi_group_detections']
        atlas_classified = data['atlas_classifications']
        
        pfas_rate = pfas_detected / total * 100 if total > 0 else 0
        multi_rate = multi_detected / total * 100 if total > 0 else 0
        atlas_rate = atlas_classified / total * 100 if total > 0 else 0
        
        def get_performance_class(rate):
            if rate >= 80:
                return "excellent"
            elif rate >= 60:
                return "good"
            elif rate >= 40:
                return "moderate"
            else:
                return "poor"
        
        pfas_class = get_performance_class(pfas_rate)
        multi_class = get_performance_class(multi_rate)
        atlas_class = get_performance_class(atlas_rate)
        
        example_mol = data['molecules'][0]['smiles'] if data['molecules'] else 'N/A'
        
        html_content += f"""
                    <tr>
                        <td>{combo_str}</td>
                        <td>{total}</td>
                        <td class="success-rate {pfas_class}">{pfas_rate:.1f}%</td>
                        <td class="success-rate {multi_class}">{multi_rate:.1f}%</td>
                        <td class="success-rate {atlas_class}">{atlas_rate:.1f}%</td>
                        <td class="molecular-example">{example_mol[:50]}{'...' if len(example_mol) > 50 else ''}</td>
                    </tr>
"""
    
    html_content += f"""
                </tbody>
            </table>
        </div>

        <h2>📈 Visual Analysis</h2>
        <div class="section-container">
            <p>Interactive Sankey diagrams showing detection performance:</p>
            <ul>
                <li><strong>Single Group Sankey:</strong> <a href="single_group_sankey_{timestamp}.html">single_group_sankey_{timestamp}.html</a></li>
                <li><strong>Multi-Group Sankey:</strong> <a href="multi_group_sankey_{timestamp}.html">multi_group_sankey_{timestamp}.html</a></li>
            </ul>
            <p>Static images available in PNG and SVG formats for presentations and publications.</p>
        </div>

        <div class="highlight-section">
            <h3>🔍 Key Findings</h3>
            <div style="text-align: left; max-width: 800px; margin: 20px auto;">
                <p><strong>Single Group Performance:</strong></p>
                <ul>
                    <li>PFASGroups achieved {single_pfas_success/total_single_molecules*100:.1f}% detection rate across all single groups</li>
                    <li>PFAS-Atlas achieved {single_atlas_success/total_single_molecules*100:.1f}% classification rate</li>
                    <li>Best performing groups show excellent correlation between systems</li>
                </ul>
                
                <p><strong>Multi-Group Challenge:</strong></p>
                <ul>
                    <li>Multi-functional molecules achieved {multi_pfas_success/total_multi_molecules*100:.1f}% PFASGroups detection</li>
                    <li>Complex molecules present detection challenges for both systems</li>
                    <li>Atlas classification: {multi_atlas_success/total_multi_molecules*100:.1f}% success on multi-group molecules</li>
                </ul>
            </div>
        </div>

        <h2>🎯 Recommendations</h2>
        <div class="section-container">
            <ol>
                <li><strong>Algorithm Enhancement:</strong> Focus on improving multi-functional group detection</li>
                <li><strong>Training Data:</strong> Expand datasets with complex multi-group molecules</li>
                <li><strong>Cross-Validation:</strong> Use both systems complementarily for comprehensive analysis</li>
                <li><strong>Threshold Tuning:</strong> Optimize detection thresholds for different molecular complexities</li>
            </ol>
        </div>

        <div style="text-align: center; color: #666; font-size: 0.9em; margin-top: 50px; padding-top: 30px; border-top: 2px solid #eee;">
            Generated on {datetime.now().strftime('%Y-%m-%d at %H:%M:%S')} | 
            Total Molecules Analyzed: {total_single_molecules + total_multi_molecules}
        </div>
    </div>
</body>
</html>
"""
    
    return html_content

def main():
    """Main analysis function"""
    
    print("🚀 PFAS BENCHMARK ANALYSIS AND VISUALIZATION")
    print("=" * 60)
    
    # Find the most recent benchmark file
    import glob
    benchmark_files = glob.glob('pfas_comprehensive_benchmark_*.json')
    if not benchmark_files:
        print("❌ No benchmark results found. Run comprehensive_pfas_benchmark.py first.")
        return
    
    latest_file = max(benchmark_files)
    print(f"📊 Using benchmark file: {latest_file}")
    
    # Load and analyze results
    results = load_benchmark_results(latest_file)
    
    # Analyze performance
    print("📋 Analyzing single group performance...")
    single_group_analysis = analyze_single_group_performance(results)
    
    print("🔬 Analyzing multi-group performance...")
    multi_group_analysis = analyze_multi_group_performance(results)
    
    # Create visualizations
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    
    print("📊 Creating Sankey diagrams...")
    
    # Single group Sankey
    fig1 = create_single_group_sankey(single_group_analysis)
    fig1.write_html(f"single_group_sankey_{timestamp}.html")
    fig1.write_image(f"single_group_sankey_{timestamp}.svg")
    fig1.write_image(f"single_group_sankey_{timestamp}.png", width=1200, height=700, scale=2)
    
    # Multi-group Sankey  
    fig2 = create_multi_group_sankey(multi_group_analysis)
    fig2.write_html(f"multi_group_sankey_{timestamp}.html")
    fig2.write_image(f"multi_group_sankey_{timestamp}.svg")
    fig2.write_image(f"multi_group_sankey_{timestamp}.png", width=1200, height=600, scale=2)
    
    # Create HTML analysis
    print("🌐 Creating HTML analysis report...")
    html_content = create_html_analysis(single_group_analysis, multi_group_analysis, timestamp)
    
    html_filename = f"pfas_benchmark_analysis_{timestamp}.html"
    with open(html_filename, 'w', encoding='utf-8') as f:
        f.write(html_content)
    
    print(f"\n✅ Analysis Complete!")
    print("📁 Generated files:")
    print(f"   - {html_filename} (comprehensive analysis)")
    print(f"   - single_group_sankey_{timestamp}.html/.svg/.png")
    print(f"   - multi_group_sankey_{timestamp}.html/.svg/.png")

if __name__ == "__main__":
    main()