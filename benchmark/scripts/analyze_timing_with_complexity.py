"""
Timing Performance Analysis with Graph Complexity Metrics
Analyzes timing benchmark results including new graph-based complexity metrics
"""

import json
import numpy as np
import pandas as pd
import plotly.graph_objects as go
import plotly.express as px
from plotly.subplots import make_subplots
from datetime import datetime
from scipy import stats as scipy_stats
import sys
import os

def load_timing_results(filename):
    """Load timing benchmark results from JSON"""
    
    with open(filename, 'r') as f:
        results = json.load(f)
    
    print(f"Loaded {len(results)} timing benchmark results")
    return results

def analyze_complexity_metrics(timing_results):
    """Analyze graph complexity metrics if available"""
    
    # Check if complexity metrics are present
    first_result = timing_results[0]
    has_complexity = 'complexity_score' in first_result
    
    if not has_complexity:
        print("⚠️  No graph complexity metrics found in data")
        return None
    
    # Extract complexity metrics
    complexity_data = {
        'diameter': [r.get('complexity_diameter', 0) for r in timing_results],
        'radius': [r.get('complexity_radius', 0) for r in timing_results],
        'avg_eccentricity': [r.get('complexity_avg_eccentricity', 0) for r in timing_results],
        'max_eccentricity': [r.get('complexity_max_eccentricity', 0) for r in timing_results],
        'avg_degree': [r.get('complexity_avg_degree', 0) for r in timing_results],
        'density': [r.get('complexity_density', 0) for r in timing_results],
        'num_cycles': [r.get('complexity_num_cycles', 0) for r in timing_results],
        'avg_betweenness': [r.get('complexity_avg_betweenness', 0) for r in timing_results],
        'max_betweenness': [r.get('complexity_max_betweenness', 0) for r in timing_results],
        'avg_clustering': [r.get('complexity_avg_clustering', 0) for r in timing_results],
        'complexity_score': [r.get('complexity_score', 0) for r in timing_results]
    }
    
    # Calculate statistics for each metric
    complexity_stats = {}
    for metric, values in complexity_data.items():
        complexity_stats[metric] =  {
            'mean': np.mean(values),
            'std': np.std(values),
            'min': np.min(values),
            'max': np.max(values),
            'median': np.median(values)
        }
    
    print(f"\n📊 Graph Complexity Metrics Summary:")
    print(f"   • Complexity Score: {complexity_stats['complexity_score']['mean']:.2f}±{complexity_stats['complexity_score']['std']:.2f} (range: {complexity_stats['complexity_score']['min']:.2f}-{complexity_stats['complexity_score']['max']:.2f})")
    print(f"   • Diameter: {complexity_stats['diameter']['mean']:.2f}±{complexity_stats['diameter']['std']:.2f} (range: {int(complexity_stats['diameter']['min'])}-{int(complexity_stats['diameter']['max'])})")
    print(f"   • Avg Eccentricity: {complexity_stats['avg_eccentricity']['mean']:.2f}±{complexity_stats['avg_eccentricity']['std']:.2f}")
    print(f"   • Avg Degree: {complexity_stats['avg_degree']['mean']:.2f}±{complexity_stats['avg_degree']['std']:.2f}")
    print(f"   • Number of Cycles: {complexity_stats['num_cycles']['mean']:.2f}±{complexity_stats['num_cycles']['std']:.2f}")
    
    return complexity_data, complexity_stats

def calculate_complexity_correlations(timing_results, complexity_data):
    """Calculate correlations between complexity metrics and timing"""
    
    # Extract timing data
    if 'HalogenGroups_time_avg' in timing_results[0]:
        pfas_times = np.array([r['HalogenGroups_time_avg'] * 1000 for r in timing_results])
        atlas_times = np.array([r['atlas_time_avg'] * 1000 for r in timing_results])
    else:
        pfas_times = np.array([r['HalogenGroups_time'] * 1000 for r in timing_results])
        atlas_times = np.array([r['atlas_time'] * 1000 for r in timing_results])
    
    num_atoms = np.array([r['num_atoms'] for r in timing_results])
    
    # Calculate correlations
    correlations = {
        'metric': [],
        'pfas_correlation': [],
        'pfas_pvalue': [],
        'atlas_correlation': [],
        'atlas_pvalue': []
    }
    
    for metric, values in complexity_data.items():
        values_array = np.array(values)
        
        # HalogenGroups correlation
        pfas_corr, pfas_pval = scipy_stats.pearsonr(values_array, pfas_times)
        
        # Atlas correlation
        atlas_corr, atlas_pval = scipy_stats.pearsonr(values_array, atlas_times)
        
        correlations['metric'].append(metric)
        correlations['pfas_correlation'].append(pfas_corr)
        correlations['pfas_pvalue'].append(pfas_pval)
        correlations['atlas_correlation'].append(atlas_corr)
        correlations['atlas_pvalue'].append(atlas_pval)
    
    # Add num_atoms for comparison
    atom_pfas_corr, atom_pfas_pval = scipy_stats.pearsonr(num_atoms, pfas_times)
    atom_atlas_corr, atom_atlas_pval = scipy_stats.pearsonr(num_atoms, atlas_times)
    
    correlations['metric'].append('num_atoms')
    correlations['pfas_correlation'].append(atom_pfas_corr)
    correlations['pfas_pvalue'].append(atom_pfas_pval)
    correlations['atlas_correlation'].append(atlas_corr)
    correlations['atlas_pvalue'].append(atlas_pval)
    
    corr_df = pd.DataFrame(correlations)
    
    print(f"\n🔗 Correlation Analysis (Complexity Metrics vs Timing):")
    print(f"\nTop correlations with HalogenGroups timing:")
    pfas_top = corr_df.nlargest(5, 'pfas_correlation')[['metric', 'pfas_correlation', 'pfas_pvalue']]
    for _, row in pfas_top.iterrows():
        sig = "***" if row['pfas_pvalue'] < 0.001 else "**" if row['pfas_pvalue'] < 0.01 else "*" if row['pfas_pvalue'] < 0.05 else ""
        print(f"   • {row['metric']}: r={row['pfas_correlation']:.3f} (p={row['pfas_pvalue']:.4f}){sig}")
    
    print(f"\nTop correlations with Atlas timing:")
    atlas_top = corr_df.nlargest(5, 'atlas_correlation')[['metric', 'atlas_correlation', 'atlas_pvalue']]
    for _, row in atlas_top.iterrows():
        sig = "***" if row['atlas_pvalue'] < 0.001 else "**" if row['atlas_pvalue'] < 0.01 else "*" if row['atlas_pvalue'] < 0.05 else ""
        print(f"   • {row['metric']}: r={row['atlas_correlation']:.3f} (p={row['atlas_pvalue']:.4f}){sig}")
    
    return corr_df

def create_complexity_scatter_plot(timing_results, complexity_data):
    """Create scatter plot showing timing vs complexity score"""
    
    complexity_scores = complexity_data['complexity_score']
    
    # Handle both old and new timing format
    if 'HalogenGroups_time_avg' in timing_results[0]:
        pfas_times_ms = [r['HalogenGroups_time_avg'] * 1000 for r in timing_results]
        atlas_times_ms = [r['atlas_time_avg'] * 1000 for r in timing_results]
        pfas_errors = [r.get('HalogenGroups_time_std', 0) * 1000 for r in timing_results]
        atlas_errors = [r.get('atlas_time_std', 0) * 1000 for r in timing_results]
        has_error_bars = any(pfas_errors) or any(atlas_errors)
    else:
        pfas_times_ms = [r['HalogenGroups_time'] * 1000 for r in timing_results]
        atlas_times_ms = [r['atlas_time'] * 1000 for r in timing_results]
        pfas_errors = None
        atlas_errors = None
        has_error_bars = False
    
    num_atoms = [r['num_atoms'] for r in timing_results]
    
    fig = go.Figure()
    
    # HalogenGroups scatter
    pfas_trace_args = {
        'x': complexity_scores,
        'y': pfas_times_ms,
        'mode': 'markers',
        'name': 'HalogenGroups',
        'marker': dict(
            color=num_atoms,
            colorscale='Blues',
            size=8,
            colorbar=dict(title="Num Atoms", x=1.02, len=0.5, y=0.75),
            showscale=True
        ),
        'text': [f"Complexity: {cs:.2f}, Atoms: {na}, Time: {pt:.2f}ms" + 
                (f" ±{pe:.2f}ms" if has_error_bars and pe > 0 else "")
                for cs, na, pt, pe in zip(complexity_scores, num_atoms, pfas_times_ms, pfas_errors or [0]*len(pfas_times_ms))],
        'hovertemplate': '%{text}<extra>HalogenGroups</extra>'
    }
    
    if has_error_bars and pfas_errors:
        pfas_trace_args['error_y'] = dict(type='data', array=pfas_errors, visible=True)
    
    fig.add_trace(go.Scatter(**pfas_trace_args))
    
    # Atlas scatter
    atlas_trace_args = {
        'x': complexity_scores,
        'y': atlas_times_ms,
        'mode': 'markers',
        'name': 'PFAS-Atlas',
        'marker': dict(
            color=num_atoms,
            colorscale='Reds',
            size=8,
            colorbar=dict(title="Num Atoms", x=1.12, len=0.5, y=0.25),
            showscale=True
        ),
        'text': [f"Complexity: {cs:.2f}, Atoms: {na}, Time: {at:.2f}ms" + 
                (f" ±{ae:.2f}ms" if has_error_bars and ae > 0 else "")
                for cs, na, at, ae in zip(complexity_scores, num_atoms, atlas_times_ms, atlas_errors or [0]*len(atlas_times_ms))],
        'hovertemplate': '%{text}<extra>PFAS-Atlas</extra>'
    }
    
    if has_error_bars and atlas_errors:
        atlas_trace_args['error_y'] = dict(type='data', array=atlas_errors, visible=True)
    
    fig.add_trace(go.Scatter(**atlas_trace_args))
    
    title_suffix = f" ({timing_results[0].get('iterations', 1)} iterations)" if timing_results[0].get('iterations', 1) > 1 else ""
    
    fig.update_layout(
        title=f"Timing vs Graph Complexity Score{title_suffix}<br><sub>How execution time scales with molecular graph complexity</sub>",
        xaxis_title="Graph Complexity Score",
        yaxis_title="Execution Time (ms)",
        width=1200,
        height=700,
        margin=dict(r=150)
    )
    
    return fig

def create_correlation_heatmap(corr_df):
    """Create heatmap showing correlations between metrics and timing"""
    
    # Prepare data for heatmap
    metrics = corr_df['metric'].tolist()
    pfas_corrs = corr_df['pfas_correlation'].tolist()
    atlas_corrs = corr_df['atlas_correlation'].tolist()
    
    # Create combined matrix
    correlation_matrix = np.array([pfas_corrs, atlas_corrs])
    
    fig = go.Figure(data=go.Heatmap(
        z=correlation_matrix,
        x=metrics,
        y=['HalogenGroups', 'PFAS-Atlas'],
        colorscale='RdBu',
        zmid=0,
        text=correlation_matrix,
        texttemplate='%{text:.3f}',
        textfont={"size": 10},
        colorbar=dict(title="Correlation")
    ))
    
    fig.update_layout(
        title="Correlation: Complexity Metrics vs Execution Time<br><sub>Pearson correlation coefficients</sub>",
        xaxis_title="Complexity Metric",
        yaxis_title="System",
        width=1400,
        height=400,
        xaxis=dict(tickangle=-45)
    )
    
    return fig

def create_complexity_breakdown_plot(timing_results, complexity_data):
    """Create plot showing timing by complexity quartiles"""
    
    complexity_scores = np.array(complexity_data['complexity_score'])
    
    # Calculate quartiles
    q1, q2, q3 = np.percentile(complexity_scores, [25, 50, 75])
    
    # Categorize molecules
    categories = []
    for score in complexity_scores:
        if score <= q1:
            categories.append('Low (Q1)')
        elif score <= q2:
            categories.append('Medium-Low (Q2)')
        elif score <= q3:
            categories.append('Medium-High (Q3)')
        else:
            categories.append('High (Q4)')
    
    # Extract timing data
    if 'HalogenGroups_time_avg' in timing_results[0]:
        pfas_times_ms = np.array([r['HalogenGroups_time_avg'] * 1000 for r in timing_results])
        atlas_times_ms = np.array([r['atlas_time_avg'] * 1000 for r in timing_results])
    else:
        pfas_times_ms = np.array([r['HalogenGroups_time'] * 1000 for r in timing_results])
        atlas_times_ms = np.array([r['atlas_time'] * 1000 for r in timing_results])
    
    # Calculate statistics by category
    category_order = ['Low (Q1)', 'Medium-Low (Q2)', 'Medium-High (Q3)', 'High (Q4)']
    pfas_stats = []
    atlas_stats = []
    sample_counts = []
    
    for cat in category_order:
        mask = np.array(categories) == cat
        pfas_stats.append({
            'mean': np.mean(pfas_times_ms[mask]),
            'std': np.std(pfas_times_ms[mask])
        })
        atlas_stats.append({
            'mean': np.mean(atlas_times_ms[mask]),
            'std': np.std(atlas_times_ms[mask])
        })
        sample_counts.append(np.sum(mask))
    
    fig = make_subplots(
        rows=2, cols=1,
        subplot_titles=['Average Execution Time by Complexity Quartile', 'Sample Count per Quartile'],
        vertical_spacing=0.15
    )
    
    # Timing comparison
    fig.add_trace(
        go.Bar(
            name='HalogenGroups',
            x=category_order,
            y=[s['mean'] for s in pfas_stats],
            error_y=dict(type='data', array=[s['std'] for s in pfas_stats]),
            marker_color='lightblue'
        ),
        row=1, col=1
    )
    
    fig.add_trace(
        go.Bar(
            name='PFAS-Atlas',
            x=category_order,
            y=[s['mean'] for s in atlas_stats],
            error_y=dict(type='data', array=[s['std'] for s in atlas_stats]),
            marker_color='orange'
        ),
        row=1, col=1
    )
    
    # Sample counts
    fig.add_trace(
        go.Bar(
            name='Sample Count',
            x=category_order,
            y=sample_counts,
            marker_color='lightgreen',
            showlegend=False
        ),
        row=2, col=1
    )
    
    fig.update_layout(
        title="Performance Analysis by Graph Complexity<br><sub>Molecules grouped by complexity score quartiles</sub>",
        barmode='group',
        width=1200,
        height=800
    )
    
    fig.update_xaxes(title_text="Complexity Quartile", row=2, col=1)
    fig.update_yaxes(title_text="Execution Time (ms)", row=1, col=1)
    fig.update_yaxes(title_text="Number of Molecules", row=2, col=1)
    
    return fig

def create_multi_metric_analysis(timing_results, complexity_data):
    """Create subplots showing individual complexity metrics vs timing"""
    
    metrics_to_plot = [
        ('diameter', 'Diameter'),
        ('avg_eccentricity', 'Avg Eccentricity'),
        ('num_cycles', 'Number of Cycles'),
        ('avg_betweenness', 'Avg Betweenness')
    ]
    
    if 'HalogenGroups_time_avg' in timing_results[0]:
        pfas_times_ms = [r['HalogenGroups_time_avg'] * 1000 for r in timing_results]
    else:
        pfas_times_ms = [r['HalogenGroups_time'] * 1000 for r in timing_results]
    
    fig = make_subplots(
        rows=2, cols=2,
        subplot_titles=[title for _, title in metrics_to_plot],
        vertical_spacing=0.12,
        horizontal_spacing=0.10
    )
    
    for idx, (metric_key, metric_name) in enumerate(metrics_to_plot):
        row = (idx // 2) + 1
        col = (idx % 2) + 1
        
        metric_values = complexity_data[metric_key]
        
        fig.add_trace(
            go.Scatter(
                x=metric_values,
                y=pfas_times_ms,
                mode='markers',
                name=metric_name,
                marker=dict(size=6, opacity=0.6),
                showlegend=False
            ),
            row=row, col=col
        )
        
        fig.update_xaxes(title_text=metric_name, row=row, col=col)
        fig.update_yaxes(title_text="HalogenGroups Time (ms)", row=row, col=col)
    
    fig.update_layout(
        title="Individual Complexity Metrics vs Timing<br><sub>HalogenGroups execution time vs different graph metrics</sub>",
        width=1400,
        height=900
    )
    
    return fig

def analyze_timing_performance(timing_results):
    """Analyze timing performance and scaling characteristics with iteration statistics"""
    
    # Extract data arrays - handle both old and new format
    chain_lengths = [r['chain_length'] for r in timing_results]
    mol_weights = [r['molecular_weight'] for r in timing_results]
    num_atoms = [r['num_atoms'] for r in timing_results]
    num_bonds = [r['num_bonds'] for r in timing_results]
    
    # Handle both old single-run and new multi-iteration format
    if 'HalogenGroups_time_avg' in timing_results[0]:  # New format with iterations
        pfas_times_ms = [r['HalogenGroups_time_avg'] * 1000 for r in timing_results]
        atlas_times_ms = [r['atlas_time_avg'] * 1000 for r in timing_results]
        pfas_stds_ms = [r.get('HalogenGroups_time_std', 0) * 1000 for r in timing_results]
        atlas_stds_ms = [r.get('atlas_time_std', 0) * 1000 for r in timing_results]
        iterations = timing_results[0].get('iterations', 1)
        system_specs = timing_results[0].get('system_specs', {})
    else:  # Old format for backward compatibility
        pfas_times_ms = [r['HalogenGroups_time'] * 1000 for r in timing_results]
        atlas_times_ms = [r['atlas_time'] * 1000 for r in timing_results]
        pfas_stds_ms = [0] * len(timing_results)
        atlas_stds_ms = [0] * len(timing_results)
        iterations = 1
        system_specs = {}
    
    # Calculate statistics
    stats = {
        'total_molecules': len(timing_results),
        'iterations_per_molecule': iterations,
        'HalogenGroups_avg_time': np.mean(pfas_times_ms),
        'HalogenGroups_std_time': np.std(pfas_times_ms),
        'HalogenGroups_median_time': np.median(pfas_times_ms),
        'HalogenGroups_avg_std': np.mean(pfas_stds_ms),
        'atlas_avg_time': np.mean(atlas_times_ms),
        'atlas_std_time': np.std(atlas_times_ms),
        'atlas_median_time': np.median(atlas_times_ms),
        'atlas_avg_std': np.mean(atlas_stds_ms),
        'speed_ratio': np.mean(atlas_times_ms) / np.mean(pfas_times_ms),
        'chain_length_range': (min(chain_lengths), max(chain_lengths)),
        'mol_weight_range': (min(mol_weights), max(mol_weights)),
        'atom_count_range': (min(num_atoms), max(num_atoms)),
        'system_specs': system_specs
    }
    
    print(f"\nTiming Performance Analysis:")
    print(f"   • Total molecules tested: {stats['total_molecules']}")
    print(f"   • Iterations per molecule: {stats['iterations_per_molecule']}")
    if iterations > 1:
        print(f"   • HalogenGroups: {stats['HalogenGroups_avg_time']:.2f}±{stats['HalogenGroups_std_time']:.2f}ms avg (median: {stats['HalogenGroups_median_time']:.2f}ms)")
        print(f"   • PFAS-Atlas: {stats['atlas_avg_time']:.2f}±{stats['atlas_std_time']:.2f}ms avg (median: {stats['atlas_median_time']:.2f}ms)")
    else:
        print(f"   • HalogenGroups: {stats['HalogenGroups_avg_time']:.2f}±{stats['HalogenGroups_std_time']:.2f}ms avg")
        print(f"   • PFAS-Atlas: {stats['atlas_avg_time']:.2f}±{stats['atlas_std_time']:.2f}ms avg")
    print(f"   • Speed ratio: {stats['speed_ratio']:.1f}x (Atlas/HalogenGroups)")
    print(f"   • Atom count range: {stats['atom_count_range'][0]}-{stats['atom_count_range'][1]}")
    
    return stats

def main():
    """Main timing analysis function with complexity metrics"""
    
    if len(sys.argv) < 2:
        print("Usage: python analyze_timing_with_complexity.py <timing_results.json>")
        sys.exit(1)
    
    filename = sys.argv[1]
    print(f"TIMING PERFORMANCE ANALYSIS WITH GRAPH COMPLEXITY")
    print(f"=" * 60)
    print(f"Loading results from: {filename}")
    
    # Load and analyze results
    timing_results = load_timing_results(filename)
    stats = analyze_timing_performance(timing_results)
    
    # Analyze complexity metrics
    complexity_result = analyze_complexity_metrics(timing_results)
    
    if complexity_result is None:
        print("\n⚠️  This dataset does not contain graph complexity metrics.")
        print("   Please run the timing benchmark with the updated enhanced_pfas_benchmark.py")
        print("   that includes complexity calculations.")
        sys.exit(1)
    
    complexity_data, complexity_stats = complexity_result
    
    # Calculate correlations
    corr_df = calculate_complexity_correlations(timing_results, complexity_data)
    
    # Generate visualizations
    print(f"\n📊 Creating complexity-aware visualizations...")
    
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    
    # Determine output directories
    imgs_dir = "../imgs" if os.path.exists("../imgs") else "imgs"
    reports_dir = "../reports" if os.path.exists("../reports") else "reports"
    os.makedirs(imgs_dir, exist_ok=True)
    os.makedirs(reports_dir, exist_ok=True)
    
    # Create and save complexity visualizations
    complexity_scatter = create_complexity_scatter_plot(timing_results, complexity_data)
    complexity_scatter.write_image(f"{imgs_dir}/timing_complexity_scatter_{timestamp}.png", width=1200, height=700)
    complexity_scatter.write_html(f"{reports_dir}/timing_complexity_scatter_{timestamp}.html")
    print(f"   ✓ Complexity scatter plot saved")
    
    correlation_heatmap = create_correlation_heatmap(corr_df)
    correlation_heatmap.write_image(f"{imgs_dir}/timing_correlation_heatmap_{timestamp}.png", width=1400, height=400)
    correlation_heatmap.write_html(f"{reports_dir}/timing_correlation_heatmap_{timestamp}.html")
    print(f"   ✓ Correlation heatmap saved")
    
    complexity_breakdown = create_complexity_breakdown_plot(timing_results, complexity_data)
    complexity_breakdown.write_image(f"{imgs_dir}/timing_complexity_breakdown_{timestamp}.png", width=1200, height=800)
    complexity_breakdown.write_html(f"{reports_dir}/timing_complexity_breakdown_{timestamp}.html")
    print(f"   ✓ Complexity breakdown plot saved")
    
    multi_metric = create_multi_metric_analysis(timing_results, complexity_data)
    multi_metric.write_image(f"{imgs_dir}/timing_multi_metric_{timestamp}.png", width=1400, height=900)
    multi_metric.write_html(f"{reports_dir}/timing_multi_metric_{timestamp}.html")
    print(f"   ✓ Multi-metric analysis saved")
    
    # Save correlation data
    data_dir = "../data" if os.path.exists("../data") else "data"
    os.makedirs(data_dir, exist_ok=True)
    
    corr_df.to_csv(f"{data_dir}/timing_correlations_{timestamp}.csv", index=False)
    print(f"   ✓ Correlation data saved to CSV")
    
    # Save comprehensive JSON report
    json_report = {
        "summary": {
            "total_molecules": stats['total_molecules'],
            "iterations": stats['iterations_per_molecule'],
            "HalogenGroups_avg_time": stats['HalogenGroups_avg_time'],
            "atlas_avg_time": stats['atlas_avg_time'],
            "speed_ratio": stats['speed_ratio'],
            "timestamp": timestamp
        },
        "complexity_stats": complexity_stats,
        "correlations": corr_df.to_dict('records'),
        "figures": [
            {
                "title": "Timing vs Complexity Score",
                "file": f"timing_complexity_scatter_{timestamp}.png",
                "description": "Execution time vs graph complexity score"
            },
            {
                "title": "Correlation Heatmap",
                "file": f"timing_correlation_heatmap_{timestamp}.png",
                "description": "Correlations between complexity metrics and timing"
            },
            {
                "title": "Complexity Breakdown",
                "file": f"timing_complexity_breakdown_{timestamp}.png",
                "description": "Performance by complexity quartiles"
            },
            {
                "title": "Multi-Metric Analysis",
                "file": f"timing_multi_metric_{timestamp}.png",
                "description": "Individual complexity metrics vs timing"
            }
        ]
    }
    
    json_filename = f"{data_dir}/timing_complexity_analysis_{timestamp}.json"
    with open(json_filename, 'w', encoding='utf-8') as f:
        json.dump(json_report, f, indent=2)
    
    print(f"\n✅ Complexity Analysis Complete!")
    print(f"\nGenerated Files:")
    print(f"   • timing_complexity_scatter_{timestamp}.png/.html")
    print(f"   • timing_correlation_heatmap_{timestamp}.png/.html")
    print(f"   • timing_complexity_breakdown_{timestamp}.png/.html")
    print(f"   • timing_multi_metric_{timestamp}.png/.html")
    print(f"   • {json_filename}")
    print(f"   • timing_correlations_{timestamp}.csv")
    
    # Key insights
    print(f"\n🔍 Key Complexity Insights:")
    
    # Find strongest correlations
    pfas_strongest = corr_df.loc[corr_df['pfas_correlation'].idxmax()]
    atlas_strongest = corr_df.loc[corr_df['atlas_correlation'].idxmax()]
    
    print(f"   • Strongest predictor of HalogenGroups timing: {pfas_strongest['metric']} (r={pfas_strongest['pfas_correlation']:.3f})")
    print(f"   • Strongest predictor of Atlas timing: {atlas_strongest['metric']} (r={atlas_strongest['atlas_correlation']:.3f})")
    
    # Check if complexity is better than atom count
    atom_row = corr_df[corr_df['metric'] == 'num_atoms'].iloc[0]
    complexity_row = corr_df[corr_df['metric'] == 'complexity_score'].iloc[0]
    
    if complexity_row['pfas_correlation'] > atom_row['pfas_correlation']:
        improvement = ((complexity_row['pfas_correlation'] - atom_row['pfas_correlation']) / atom_row['pfas_correlation']) * 100
        print(f"   • Complexity score is {improvement:.1f}% better predictor than atom count for HalogenGroups")
    
    if complexity_row['atlas_correlation'] > atom_row['atlas_correlation']:
        improvement = ((complexity_row['atlas_correlation'] - atom_row['atlas_correlation']) / atom_row['atlas_correlation']) * 100
        print(f"   • Complexity score is {improvement:.1f}% better predictor than atom count for Atlas")

if __name__ == "__main__":
    main()
