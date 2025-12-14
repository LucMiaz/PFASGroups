"""
Timing Performance Analysis
Analyzes timing benchmark results and creates performance visualizations
"""

import json
import numpy as np
import plotly.graph_objects as go
import plotly.express as px
from plotly.subplots import make_subplots
from datetime import datetime
import sys

def load_timing_results(filename):
    """Load timing benchmark results from JSON"""
    
    with open(filename, 'r') as f:
        results = json.load(f)
    
    print(f"📊 Loaded {len(results)} timing benchmark results")
    return results

def analyze_timing_performance(timing_results):
    """Analyze timing performance and scaling characteristics with iteration statistics"""
    
    # Extract data arrays - handle both old and new format
    chain_lengths = [r['chain_length'] for r in timing_results]
    mol_weights = [r['molecular_weight'] for r in timing_results]
    num_atoms = [r['num_atoms'] for r in timing_results]
    num_bonds = [r['num_bonds'] for r in timing_results]
    
    # Handle both old single-run and new multi-iteration format
    if 'pfasgroups_time_avg' in timing_results[0]:  # New format with iterations
        pfas_times_ms = [r['pfasgroups_time_avg'] * 1000 for r in timing_results]
        atlas_times_ms = [r['atlas_time_avg'] * 1000 for r in timing_results]
        pfas_stds_ms = [r.get('pfasgroups_time_std', 0) * 1000 for r in timing_results]
        atlas_stds_ms = [r.get('atlas_time_std', 0) * 1000 for r in timing_results]
        iterations = timing_results[0].get('iterations', 1)
        
        # Extract system specifications if available
        system_specs = timing_results[0].get('system_specs', {})
    else:  # Old format for backward compatibility
        pfas_times_ms = [r['pfasgroups_time'] * 1000 for r in timing_results]
        atlas_times_ms = [r['atlas_time'] * 1000 for r in timing_results]
        pfas_stds_ms = [0] * len(timing_results)
        atlas_stds_ms = [0] * len(timing_results)
        iterations = 1
        system_specs = {}
    
    # Calculate statistics
    stats = {
        'total_molecules': len(timing_results),
        'iterations_per_molecule': iterations,
        'pfasgroups_avg_time': np.mean(pfas_times_ms),
        'pfasgroups_std_time': np.std(pfas_times_ms),
        'pfasgroups_median_time': np.median(pfas_times_ms),
        'pfasgroups_avg_std': np.mean(pfas_stds_ms),  # Average of individual std devs
        'atlas_avg_time': np.mean(atlas_times_ms),
        'atlas_std_time': np.std(atlas_times_ms),
        'atlas_median_time': np.median(atlas_times_ms),
        'atlas_avg_std': np.mean(atlas_stds_ms),  # Average of individual std devs
        'speed_ratio': np.mean(atlas_times_ms) / np.mean(pfas_times_ms),
        'chain_length_range': (min(chain_lengths), max(chain_lengths)),
        'mol_weight_range': (min(mol_weights), max(mol_weights)),
        'atom_count_range': (min(num_atoms), max(num_atoms)),
        'system_specs': system_specs
    }
    
    print(f"\n📈 Timing Performance Analysis:")
    print(f"   • Total molecules tested: {stats['total_molecules']}")
    print(f"   • Iterations per molecule: {stats['iterations_per_molecule']}")
    if iterations > 1:
        print(f"   • PFASGroups: {stats['pfasgroups_avg_time']:.2f}±{stats['pfasgroups_std_time']:.2f}ms avg (median: {stats['pfasgroups_median_time']:.2f}ms) | Individual std: {stats['pfasgroups_avg_std']:.2f}ms")
        print(f"   • PFAS-Atlas: {stats['atlas_avg_time']:.2f}±{stats['atlas_std_time']:.2f}ms avg (median: {stats['atlas_median_time']:.2f}ms) | Individual std: {stats['atlas_avg_std']:.2f}ms")
    else:
        print(f"   • PFASGroups: {stats['pfasgroups_avg_time']:.2f}±{stats['pfasgroups_std_time']:.2f}ms avg (median: {stats['pfasgroups_median_time']:.2f}ms)")
        print(f"   • PFAS-Atlas: {stats['atlas_avg_time']:.2f}±{stats['atlas_std_time']:.2f}ms avg (median: {stats['atlas_median_time']:.2f}ms)")
    print(f"   • Speed ratio: {stats['speed_ratio']:.1f}x (Atlas/PFASGroups)")
    print(f"   • Chain length range: {stats['chain_length_range'][0]}-{stats['chain_length_range'][1]}")
    print(f"   • Molecular weight range: {stats['mol_weight_range'][0]:.1f}-{stats['mol_weight_range'][1]:.1f}")
    print(f"   • Atom count range: {stats['atom_count_range'][0]}-{stats['atom_count_range'][1]}")
    
    # Print system specifications if available
    if system_specs:
        print(f"\n💻 System Specifications:")
        print(f"   • OS: {system_specs.get('system', 'Unknown')} ({system_specs.get('architecture', 'Unknown')})")
        print(f"   • CPU: {system_specs.get('cpu_name', 'Unknown')}")
        print(f"   • Cores: {system_specs.get('cpu_cores_physical', 'Unknown')} physical, {system_specs.get('cpu_cores_logical', 'Unknown')} logical")
        print(f"   • Memory: {system_specs.get('total_memory_gb', 'Unknown')} GB total")
        print(f"   • Python: {system_specs.get('python_version', 'Unknown')}")
    
    return stats

def create_timing_scatter_plot(timing_results):
    """Create scatter plot showing timing vs molecular complexity with error bars if available"""
    
    # Prepare data - handle both old and new format
    num_atoms = [r['num_atoms'] for r in timing_results]
    chain_lengths = [r['chain_length'] for r in timing_results]
    
    # Handle both old single-run and new multi-iteration format
    if 'pfasgroups_time_avg' in timing_results[0]:
        pfas_times_ms = [r['pfasgroups_time_avg'] * 1000 for r in timing_results]
        atlas_times_ms = [r['atlas_time_avg'] * 1000 for r in timing_results]
        pfas_errors = [r.get('pfasgroups_time_std', 0) * 1000 for r in timing_results]
        atlas_errors = [r.get('atlas_time_std', 0) * 1000 for r in timing_results]
        has_error_bars = any(pfas_errors) or any(atlas_errors)
    else:
        pfas_times_ms = [r['pfasgroups_time'] * 1000 for r in timing_results]
        atlas_times_ms = [r['atlas_time'] * 1000 for r in timing_results]
        pfas_errors = None
        atlas_errors = None
        has_error_bars = False
    
    fig = go.Figure()
    
    # PFASGroups scatter with optional error bars
    pfas_trace_args = {
        'x': num_atoms,
        'y': pfas_times_ms,
        'mode': 'markers',
        'name': 'PFASGroups',
        'marker': dict(
            color=chain_lengths,
            colorscale='Blues',
            size=8,
            colorbar=dict(title="Chain Length", x=1.02, len=0.5, y=0.75),
            showscale=True
        ),
        'text': [f"Chain: {cl}, Atoms: {na}, Time: {pt:.2f}ms" + 
                (f" ±{pe:.2f}ms" if has_error_bars and pe > 0 else "")
                for cl, na, pt, pe in zip(chain_lengths, num_atoms, pfas_times_ms, pfas_errors or [0]*len(pfas_times_ms))],
        'hovertemplate': '%{text}<extra>PFASGroups</extra>'
    }
    
    if has_error_bars and pfas_errors:
        pfas_trace_args['error_y'] = dict(type='data', array=pfas_errors, visible=True)
    
    fig.add_trace(go.Scatter(**pfas_trace_args))
    
    # PFAS-Atlas scatter with optional error bars
    atlas_trace_args = {
        'x': num_atoms,
        'y': atlas_times_ms,
        'mode': 'markers',
        'name': 'PFAS-Atlas',
        'marker': dict(
            color=chain_lengths,
            colorscale='Reds',
            size=8,
            colorbar=dict(title="Chain Length", x=1.12, len=0.5, y=0.25),
            showscale=True
        ),
        'text': [f"Chain: {cl}, Atoms: {na}, Time: {at:.2f}ms" + 
                (f" ±{ae:.2f}ms" if has_error_bars and ae > 0 else "")
                for cl, na, at, ae in zip(chain_lengths, num_atoms, atlas_times_ms, atlas_errors or [0]*len(atlas_times_ms))],
        'hovertemplate': '%{text}<extra>PFAS-Atlas</extra>'
    }
    
    if has_error_bars and atlas_errors:
        atlas_trace_args['error_y'] = dict(type='data', array=atlas_errors, visible=True)
    
    fig.add_trace(go.Scatter(**atlas_trace_args))
    
    title_suffix = f" ({timing_results[0].get('iterations', 1)} iterations per point)" if timing_results[0].get('iterations', 1) > 1 else ""
    
    fig.update_layout(
        title=f"Timing Performance vs Molecular Complexity{title_suffix}<br><sub>Execution time scaling with number of atoms</sub>",
        xaxis_title="Number of Atoms",
        yaxis_title="Execution Time (ms)",
        width=1200,
        height=700,
        margin=dict(r=150)
    )
    
    return fig

def create_timing_distribution_plot(timing_results):
    """Create histogram showing timing distributions"""
    
    # Handle both old and new data formats
    if 'pfasgroups_time' in timing_results[0]:
        pfas_times_ms = [r['pfasgroups_time'] * 1000 for r in timing_results]
        atlas_times_ms = [r['atlas_time'] * 1000 for r in timing_results]
    else:
        pfas_times_ms = [r['pfasgroups_time_avg'] * 1000 for r in timing_results]
        atlas_times_ms = [r['atlas_time_avg'] * 1000 for r in timing_results]
    
    fig = make_subplots(
        rows=2, cols=1,
        subplot_titles=['PFASGroups Execution Time Distribution', 'PFAS-Atlas Execution Time Distribution'],
        vertical_spacing=0.15
    )
    
    # PFASGroups histogram
    fig.add_trace(
        go.Histogram(x=pfas_times_ms, nbinsx=30, name='PFASGroups', marker_color='lightblue'),
        row=1, col=1
    )
    
    # PFAS-Atlas histogram
    fig.add_trace(
        go.Histogram(x=atlas_times_ms, nbinsx=30, name='PFAS-Atlas', marker_color='orange'),
        row=2, col=1
    )
    
    fig.update_layout(
        title="Execution Time Distributions<br><sub>Frequency of different execution times</sub>",
        width=1000,
        height=800,
        showlegend=False
    )
    
    fig.update_xaxes(title_text="Execution Time (ms)", row=2, col=1)
    fig.update_yaxes(title_text="Frequency", row=1, col=1)
    fig.update_yaxes(title_text="Frequency", row=2, col=1)
    
    return fig

def create_scaling_analysis_plot(timing_results):
    """Create plot showing how performance scales with molecule size"""
    
    # Group by molecule size bins
    size_bins = [(0, 10), (11, 15), (16, 20), (21, 25), (26, 30), (31, 40)]
    bin_data = []
    
    for min_atoms, max_atoms in size_bins:
        molecules_in_bin = [r for r in timing_results if min_atoms <= r['num_atoms'] <= max_atoms]
        
        if molecules_in_bin:
            # Handle both old format (single timing) and new format (average from multiple iterations)
            if 'pfasgroups_time_avg' in molecules_in_bin[0]:
                pfas_times = [r['pfasgroups_time_avg'] * 1000 for r in molecules_in_bin]
                atlas_times = [r['atlas_time_avg'] * 1000 for r in molecules_in_bin]
            else:
                pfas_times = [r['pfasgroups_time'] * 1000 for r in molecules_in_bin]
                atlas_times = [r['atlas_time'] * 1000 for r in molecules_in_bin]
            
            bin_data.append({
                'size_range': f"{min_atoms}-{max_atoms}",
                'molecule_count': len(molecules_in_bin),
                'pfas_avg': np.mean(pfas_times),
                'pfas_std': np.std(pfas_times),
                'atlas_avg': np.mean(atlas_times),
                'atlas_std': np.std(atlas_times),
                'avg_atoms': np.mean([r['num_atoms'] for r in molecules_in_bin])
            })
    
    if not bin_data:
        return go.Figure().update_layout(title="No data available for scaling analysis")
    
    size_ranges = [d['size_range'] for d in bin_data]
    pfas_avgs = [d['pfas_avg'] for d in bin_data]
    pfas_stds = [d['pfas_std'] for d in bin_data]
    atlas_avgs = [d['atlas_avg'] for d in bin_data]
    atlas_stds = [d['atlas_std'] for d in bin_data]
    molecule_counts = [d['molecule_count'] for d in bin_data]
    
    fig = make_subplots(
        rows=2, cols=1,
        subplot_titles=['Average Execution Time by Molecule Size', 'Number of Molecules per Size Bin'],
        vertical_spacing=0.15
    )
    
    # Timing comparison with error bars
    fig.add_trace(
        go.Bar(name='PFASGroups', x=size_ranges, y=pfas_avgs, 
               error_y=dict(type='data', array=pfas_stds),
               marker_color='lightblue'),
        row=1, col=1
    )
    fig.add_trace(
        go.Bar(name='PFAS-Atlas', x=size_ranges, y=atlas_avgs,
               error_y=dict(type='data', array=atlas_stds),
               marker_color='orange'),
        row=1, col=1
    )
    
    # Sample counts
    fig.add_trace(
        go.Bar(name='Sample Count', x=size_ranges, y=molecule_counts,
               marker_color='lightgreen', showlegend=False),
        row=2, col=1
    )
    
    fig.update_layout(
        title="Performance Scaling Analysis<br><sub>How execution time scales with molecule size</sub>",
        barmode='group',
        width=1200,
        height=800
    )
    
    fig.update_xaxes(title_text="Molecule Size (atoms)", row=2, col=1)
    fig.update_yaxes(title_text="Execution Time (ms)", row=1, col=1)
    fig.update_yaxes(title_text="Sample Count", row=2, col=1)
    
    return fig

def create_timing_html_report(timing_results, stats, timestamp):
    """Create comprehensive HTML report for timing analysis with system specs"""
    
    # Generate plots
    scatter_plot = create_timing_scatter_plot(timing_results)
    distribution_plot = create_timing_distribution_plot(timing_results)
    scaling_plot = create_scaling_analysis_plot(timing_results)
    
    # Convert plots to HTML
    scatter_html = scatter_plot.to_html(include_plotlyjs='cdn', div_id="timing-scatter")
    distribution_html = distribution_plot.to_html(include_plotlyjs='cdn', div_id="timing-distribution")
    scaling_html = scaling_plot.to_html(include_plotlyjs='cdn', div_id="timing-scaling")
    
    # System specifications section
    system_specs = stats.get('system_specs', {})
    if system_specs:
        system_html = f"""
        <div class="system-specs">
            <h3>💻 System Specifications</h3>
            <div class="specs-grid">
                <div class="spec-item">
                    <strong>Operating System:</strong> {system_specs.get('system', 'Unknown')} ({system_specs.get('architecture', 'Unknown')})
                </div>
                <div class="spec-item">
                    <strong>Platform:</strong> {system_specs.get('platform', 'Unknown')}
                </div>
                <div class="spec-item">
                    <strong>CPU:</strong> {system_specs.get('cpu_name', 'Unknown CPU')}
                </div>
                <div class="spec-item">
                    <strong>CPU Cores:</strong> {system_specs.get('cpu_cores_physical', 'Unknown')} physical, {system_specs.get('cpu_cores_logical', 'Unknown')} logical
                </div>
                <div class="spec-item">
                    <strong>Memory:</strong> {system_specs.get('total_memory_gb', 'Unknown')} GB total, {system_specs.get('available_memory_gb', 'Unknown')} GB available
                </div>
                <div class="spec-item">
                    <strong>Python Version:</strong> {system_specs.get('python_version', 'Unknown')}
                </div>
            </div>
        </div>"""
    else:
        system_html = ""
        
    # Iteration information
    iterations = stats.get('iterations_per_molecule', 1)
    iteration_info = f"""
        <div class="iteration-info">
            <h3>🔄 Statistical Reliability</h3>
            <p><strong>Iterations per molecule:</strong> {iterations}</p>
            {'<p>Each timing measurement represents the average of multiple runs for improved statistical accuracy.</p>' if iterations > 1 else '<p>Single measurement per molecule (consider running with more iterations for better statistical reliability).</p>'}
        </div>""" if iterations > 1 else ""
    
    html_content = f"""
<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>PFAS Timing Performance Analysis - {timestamp}</title>
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
            font-size: 2.2em;
            font-weight: bold;
            margin-bottom: 10px;
        }}
        .summary-label {{
            font-size: 1.0em;
            opacity: 0.95;
            line-height: 1.3;
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
        .key-finding {{
            background: #e8f4f8;
            border-left: 4px solid #1e88e5;
            padding: 20px;
            margin: 20px 0;
            border-radius: 8px;
        }}
        .system-specs, .iteration-info {{
            background: #e8f5e8;
            border-left: 4px solid #28a745;
            padding: 20px;
            margin: 20px 0;
            border-radius: 8px;
        }}
        .specs-grid {{
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(300px, 1fr));
            gap: 15px;
            margin-top: 15px;
        }}
        .spec-item {{
            background: white;
            padding: 12px;
            border-radius: 6px;
            border: 1px solid #ddd;
        }}
        .performance-excellent {{ background-color: #4caf50; color: white; }}
        .performance-good {{ background-color: #8bc34a; color: white; }}
        .performance-moderate {{ background-color: #ff9800; color: white; }}
        .performance-poor {{ background-color: #f44336; color: white; }}
    </style>
</head>
<body>
    <div class="container">
        <h1>🚀 PFAS Timing Performance Analysis</h1>
        <p style="text-align: center; font-size: 1.2em; color: #666; margin-bottom: 50px;">
            Comprehensive timing comparison: PFASGroups vs PFAS-Atlas<br>
            Generated on {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}
        </p>

        <div class="highlight-section">
            <h2 style="border: none; color: white; margin-top: 0;">📊 Performance Summary</h2>
            <div class="summary-grid">
                <div class="summary-card">
                    <div class="summary-number">{stats['total_molecules']}</div>
                    <div class="summary-label">Molecules<br>Tested</div>
                </div>
                <div class="summary-card">
                    <div class="summary-number">{stats['pfasgroups_avg_time']:.1f}ms</div>
                    <div class="summary-label">PFASGroups<br>Average Time</div>
                </div>
                <div class="summary-card">
                    <div class="summary-number">{stats['atlas_avg_time']:.1f}ms</div>
                    <div class="summary-label">PFAS-Atlas<br>Average Time</div>
                </div>
                <div class="summary-card">
                    <div class="summary-number">{stats['speed_ratio']:.1f}x</div>
                    <div class="summary-label">Speed Ratio<br>(Atlas/PFASGroups)</div>
                </div>
                <div class="summary-card">
                    <div class="summary-number">{stats['atom_count_range'][0]}-{stats['atom_count_range'][1]}</div>
                    <div class="summary-label">Atom Count<br>Range</div>
                </div>
                <div class="summary-card">
                    <div class="summary-number">{stats['mol_weight_range'][0]:.0f}-{stats['mol_weight_range'][1]:.0f}</div>
                    <div class="summary-label">Molecular Weight<br>Range</div>
                </div>
            </div>
        </div>
        
        {system_html}
        {iteration_info}

        <h2>🎯 Performance vs Molecular Complexity</h2>
        <div class="visualization-container">
            <div style="margin-bottom: 20px;">
                <strong>Analysis:</strong> How execution time scales with molecular complexity (number of atoms).
                Points are colored by chain length to show structural diversity.
            </div>
            {scatter_html}
        </div>

        <h2>📈 Execution Time Distributions</h2>
        <div class="visualization-container">
            <div style="margin-bottom: 20px;">
                <strong>Analysis:</strong> Distribution of execution times showing the frequency of different 
                performance ranges for both detection systems.
            </div>
            {distribution_html}
        </div>

        <h2>⚖️ Scaling Performance Analysis</h2>
        <div class="visualization-container">
            <div style="margin-bottom: 20px;">
                <strong>Analysis:</strong> How performance scales across different molecule size ranges.
                Error bars show standard deviation within each size category.
            </div>
            {scaling_html}
        </div>

        <div class="key-finding">
            <h3>🔍 Key Performance Insights</h3>
            <ul>
                <li><strong>Speed Comparison:</strong> {"PFASGroups is faster" if stats['speed_ratio'] > 1 else "PFAS-Atlas is faster"} 
                    by a factor of {abs(stats['speed_ratio'] - 1):.1f}x</li>
                <li><strong>Consistency:</strong> PFASGroups std dev: {stats['pfasgroups_std_time']:.2f}ms, 
                    PFAS-Atlas std dev: {stats['atlas_std_time']:.2f}ms</li>
                <li><strong>Molecule Range:</strong> Tested molecules from {stats['atom_count_range'][0]} to {stats['atom_count_range'][1]} atoms</li>
                <li><strong>Performance Trend:</strong> {"Both systems show similar scaling patterns" if abs(stats['speed_ratio'] - 1) < 0.5 else "Significant performance difference observed"}</li>
            </ul>
        </div>

        <div style="text-align: center; margin-top: 50px; padding: 30px; background: #f8f9fa; border-radius: 10px;">
            <h3>🎯 Benchmark Methodology</h3>
            <p>This timing benchmark tested {stats['total_molecules']} molecules with carboxylic acid functional groups, 
               varying in chain length from {stats['chain_length_range'][0]} to {stats['chain_length_range'][1]} carbons. 
               Each molecule was processed by both PFASGroups and PFAS-Atlas detection systems, 
               measuring execution time with high precision timing.</p>
        </div>
    </div>
</body>
</html>
"""
    
    return html_content

def main():
    """Main timing analysis function"""
    
    if len(sys.argv) < 2:
        print("Usage: python analyze_timing.py <timing_results.json>")
        sys.exit(1)
    
    filename = sys.argv[1]
    print(f"🕒 TIMING PERFORMANCE ANALYSIS")
    print(f"=" * 50)
    print(f"📊 Loading results from: {filename}")
    
    # Load and analyze results
    timing_results = load_timing_results(filename)
    stats = analyze_timing_performance(timing_results)
    
    # Generate visualizations
    print(f"\n📊 Creating performance visualizations...")
    
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    
    # Save individual plots
    scatter_plot = create_timing_scatter_plot(timing_results)
    scatter_plot.write_image(f"timing_scatter_{timestamp}.png", width=1200, height=700)
    scatter_plot.write_html(f"timing_scatter_{timestamp}.svg")
    
    distribution_plot = create_timing_distribution_plot(timing_results)
    distribution_plot.write_image(f"timing_distribution_{timestamp}.png", width=1000, height=800)
    distribution_plot.write_html(f"timing_distribution_{timestamp}.svg")
    
    scaling_plot = create_scaling_analysis_plot(timing_results)
    scaling_plot.write_image(f"timing_scaling_{timestamp}.png", width=1200, height=800)
    scaling_plot.write_html(f"timing_scaling_{timestamp}.svg")
    
    # Create comprehensive HTML report
    html_content = create_timing_html_report(timing_results, stats, timestamp)
    html_filename = f"timing_analysis_{timestamp}.html"
    
    with open(html_filename, 'w') as f:
        f.write(html_content)
    
    print(f"\n✅ Timing Analysis Complete!")
    print(f"📁 Generated Files:")
    print(f"   • {html_filename} (comprehensive analysis)")
    print(f"   • timing_scatter_{timestamp}.png/.svg")
    print(f"   • timing_distribution_{timestamp}.png/.svg") 
    print(f"   • timing_scaling_{timestamp}.png/.svg")
    
    # Performance recommendations
    print(f"\n🎯 Performance Recommendations:")
    if stats['speed_ratio'] > 2:
        print(f"   • PFAS-Atlas is significantly slower ({stats['speed_ratio']:.1f}x) - consider optimizing for large datasets")
    elif stats['speed_ratio'] > 1.5:
        print(f"   • PFAS-Atlas is moderately slower ({stats['speed_ratio']:.1f}x) - acceptable for most use cases")
    elif stats['speed_ratio'] < 0.5:
        print(f"   • PFASGroups is significantly slower ({1/stats['speed_ratio']:.1f}x) - investigate bottlenecks")
    else:
        print(f"   • Both systems have comparable performance ({stats['speed_ratio']:.1f}x ratio)")
    
    if stats['pfasgroups_std_time'] > stats['pfasgroups_avg_time'] * 0.5:
        print(f"   • PFASGroups shows high timing variability - consider profiling edge cases")
    if stats['atlas_std_time'] > stats['atlas_avg_time'] * 0.5:
        print(f"   • PFAS-Atlas shows high timing variability - consider profiling edge cases")

if __name__ == "__main__":
    main()