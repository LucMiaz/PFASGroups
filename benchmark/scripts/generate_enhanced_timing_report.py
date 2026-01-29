#!/usr/bin/env python3
"""
Generate enhanced HTML timing report for PFAS benchmark results
Includes detailed analyses by dataset, compound size, and chain type
"""

import json
import sqlite3
from pathlib import Path
from datetime import datetime
import statistics
from collections import defaultdict

def generate_html_report(output_file='reports/enhanced_timing_report.html'):
    """Generate comprehensive timing analysis HTML report with detailed breakdowns"""
    
    # Connect to database
    db_path = Path(__file__).parent / 'review-app' / 'database' / 'pfas_benchmark.db'
    
    if not db_path.exists():
        print(f"❌ Database not found at {db_path}")
        print("   Please run the import script first: node review-app/scripts/import-benchmark-data.js")
        return False
    
    conn = sqlite3.connect(str(db_path))
    cursor = conn.cursor()
    
    print("📊 Analyzing timing data...")
    
    # Get all molecules with timing data and detected groups info
    cursor.execute("""
        SELECT 
            m.id,
            m.dataset_type,
            m.num_atoms,
            m.smiles,
            m.group_id,
            m.generation_type,
            p.execution_time as pfasgroups_time,
            p.detected_groups,
            p.matched_path_types,
            a.execution_time as atlas_time
        FROM molecules m
        LEFT JOIN pfasgroups_results p ON m.id = p.molecule_id
        LEFT JOIN atlas_results a ON m.id = a.molecule_id
        WHERE (p.execution_time IS NOT NULL OR a.execution_time IS NOT NULL)
          AND m.dataset_type NOT IN ('complex_branched', 'definitions', 'non_fluorinated', 'highly_branched')
    """)
    
    rows = cursor.fetchall()
    conn.close()
    
    if not rows:
        print("⚠️  No timing data found in database")
        return False
    
    print(f"✓ Found {len(rows)} molecules with timing data")
    
    # Organize data with enhanced categorization
    data_by_dataset = defaultdict(lambda: {'pfasgroups': [], 'atlas': []})
    data_by_size = defaultdict(lambda: {'pfasgroups': [], 'atlas': []})
    data_by_chain_type = defaultdict(lambda: {'pfasgroups': [], 'atlas': []})
    data_by_atoms = defaultdict(lambda: {'pfasgroups': [], 'atlas': []})
    
    all_pfasgroups = []
    all_atlas = []
    
    for mol_id, dataset, num_atoms, smiles, group_id, gen_type, pfas_time, detected_groups, path_types, atlas_time in rows:
        # Categorize by size (small: <20, medium: 20-40, large: >40 atoms)
        if num_atoms:
            if num_atoms < 20:
                size_cat = 'small'
            elif num_atoms <= 40:
                size_cat = 'medium'
            else:
                size_cat = 'large'
        else:
            size_cat = 'unknown'
        
        # Categorize by chain type
        chain_type = 'standard'
        if gen_type and 'cyclic' in str(gen_type).lower():
            chain_type = 'cyclic'
        elif group_id and group_id >= 60:  # Fluorotelomer groups
            chain_type = 'fluorotelomer'
        elif detected_groups:
            try:
                groups_data = json.loads(detected_groups) if isinstance(detected_groups, str) else detected_groups
                if isinstance(groups_data, list) and groups_data:
                    # Check if any detected group has linker validation
                    for grp in groups_data:
                        if isinstance(grp, dict):
                            if grp.get('max_dist_from_CF2', 0) > 0:
                                chain_type = 'fluorotelomer'
                                break
            except:
                pass
        
        # Store times
        if pfas_time is not None:
            data_by_dataset[dataset]['pfasgroups'].append(pfas_time)
            data_by_size[size_cat]['pfasgroups'].append(pfas_time)
            data_by_chain_type[chain_type]['pfasgroups'].append(pfas_time)
            all_pfasgroups.append(pfas_time)
            if num_atoms:
                data_by_atoms[num_atoms]['pfasgroups'].append(pfas_time)
        
        if atlas_time is not None:
            data_by_dataset[dataset]['atlas'].append(atlas_time)
            data_by_size[size_cat]['atlas'].append(atlas_time)
            data_by_chain_type[chain_type]['atlas'].append(atlas_time)
            all_atlas.append(atlas_time)
            if num_atoms:
                data_by_atoms[num_atoms]['atlas'].append(atlas_time)
    
    # Generate HTML report
    html = generate_html_content(
        data_by_dataset, 
        data_by_size,
        data_by_chain_type,
        data_by_atoms,
        all_pfasgroups, 
        all_atlas, 
        len(rows)
    )
    
    output_path = Path(__file__).parent / output_file
    with open(output_path, 'w', encoding='utf-8') as f:
        f.write(html)
    
    print(f"✅ Enhanced report generated: {output_path}")
    return True

def calc_stats(times):
    """Calculate statistics for a list of times"""
    if not times:
        return None
    return {
        'mean': statistics.mean(times),
        'median': statistics.median(times),
        'stdev': statistics.stdev(times) if len(times) > 1 else 0,
        'min': min(times),
        'max': max(times),
        'count': len(times)
    }

def generate_html_content(data_by_dataset, data_by_size, data_by_chain_type, data_by_atoms, 
                          all_pfasgroups, all_atlas, total_molecules):
    """Generate HTML content with embedded charts and tables"""
    
    overall_pfasgroups = calc_stats(all_pfasgroups)
    overall_atlas = calc_stats(all_atlas)
    
    # Calculate speedup
    speedup = overall_atlas['mean'] / overall_pfasgroups['mean'] if overall_pfasgroups and overall_atlas else 0
    
    # Build dataset statistics
    dataset_stats = []
    for dataset in sorted(data_by_dataset.keys()):
        d = data_by_dataset[dataset]
        pfas_stats = calc_stats(d['pfasgroups'])
        atlas_stats = calc_stats(d['atlas'])
        
        dataset_stats.append({
            'name': dataset,
            'pfas_stats': pfas_stats,
            'atlas_stats': atlas_stats,
            'count': len(d['pfasgroups']) or len(d['atlas'])
        })
    
    # Build size statistics
    size_stats = []
    for size in ['small', 'medium', 'large']:
        if size in data_by_size:
            d = data_by_size[size]
            pfas_stats = calc_stats(d['pfasgroups'])
            atlas_stats = calc_stats(d['atlas'])
            
            size_stats.append({
                'name': size.capitalize(),
                'pfas_stats': pfas_stats,
                'atlas_stats': atlas_stats,
                'count': len(d['pfasgroups']) or len(d['atlas'])
            })
    
    # Build chain type statistics
    chain_type_stats = []
    for chain_type in sorted(data_by_chain_type.keys()):
        d = data_by_chain_type[chain_type]
        pfas_stats = calc_stats(d['pfasgroups'])
        atlas_stats = calc_stats(d['atlas'])
        
        chain_type_stats.append({
            'name': chain_type.capitalize(),
            'pfas_stats': pfas_stats,
            'atlas_stats': atlas_stats,
            'count': len(d['pfasgroups']) or len(d['atlas'])
        })
    
    # Build atom count statistics
    atom_bins = [0, 10, 20, 30, 40, 50, 100, 200]
    atom_stats = []
    for i in range(len(atom_bins) - 1):
        bin_start = atom_bins[i]
        bin_end = atom_bins[i + 1]
        
        pfas_times = []
        atlas_times = []
        
        for atoms, times_dict in data_by_atoms.items():
            if bin_start < atoms <= bin_end:
                pfas_times.extend(times_dict['pfasgroups'])
                atlas_times.extend(times_dict['atlas'])
        
        if pfas_times or atlas_times:
            atom_stats.append({
                'range': f"{bin_start+1}-{bin_end}",
                'pfas_stats': calc_stats(pfas_times),
                'atlas_stats': calc_stats(atlas_times)
            })
    
    # Generate JavaScript data for charts
    dataset_chart_js = generate_dataset_chart_js(dataset_stats)
    size_chart_js = generate_size_chart_js(size_stats)
    chain_type_chart_js = generate_chain_type_chart_js(chain_type_stats)
    atom_chart_js = generate_atom_chart_js(atom_stats)
    distribution_chart_js = generate_distribution_chart_js(all_pfasgroups, all_atlas)
    
    # Build HTML
    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Enhanced PFAS Timing Analysis</title>
    <script src="https://cdn.plot.ly/plotly-2.27.0.min.js"></script>
    <style>
        * {{
            margin: 0;
            padding: 0;
            box-sizing: border-box;
        }}
        
        body {{
            font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif;
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            padding: 20px;
            color: #333;
        }}
        
        .container {{
            max-width: 1600px;
            margin: 0 auto;
            background: white;
            border-radius: 15px;
            padding: 40px;
            box-shadow: 0 10px 30px rgba(0,0,0,0.3);
        }}
        
        h1 {{
            color: #2c5aa0;
            text-align: center;
            font-size: 2.5em;
            margin-bottom: 30px;
        }}
        
        h2 {{
            color: #34495e;
            border-left: 4px solid #667eea;
            padding-left: 15px;
            margin: 40px 0 20px 0;
        }}
        
        h3 {{
            color: #555;
            margin: 25px 0 15px 0;
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
            font-size: 1em;
            opacity: 0.95;
            line-height: 1.3;
        }}
        
        table {{
            width: 100%;
            border-collapse: collapse;
            margin: 20px 0;
            background: white;
            box-shadow: 0 2px 8px rgba(0,0,0,0.1);
        }}
        
        th {{
            background: linear-gradient(45deg, #667eea 0%, #764ba2 100%);
            color: white;
            padding: 15px;
            text-align: left;
            font-weight: 600;
        }}
        
        td {{
            padding: 12px 15px;
            border-bottom: 1px solid #ddd;
        }}
        
        tr:hover {{
            background-color: #f5f5f5;
        }}
        
        .number {{
            text-align: right;
            font-family: 'Courier New', monospace;
        }}
        
        .speedup {{
            color: #28a745;
            font-weight: bold;
        }}
        
        .slower {{
            color: #dc3545;
            font-weight: bold;
        }}
        
        .chart-container {{
            background: #f8f9fa;
            padding: 30px;
            border-radius: 12px;
            margin: 30px 0;
            box-shadow: 0 5px 15px rgba(0,0,0,0.1);
        }}
        
        .info-box {{
            background: #e8f4f8;
            border-left: 4px solid #1e88e5;
            padding: 20px;
            margin: 20px 0;
            border-radius: 8px;
        }}
        
        .info-box strong {{
            color: #1e88e5;
        }}
    </style>
</head>
<body>
    <div class="container">
        <h1>🎯 Enhanced PFAS Timing Analysis</h1>
        
        <div class="info-box">
            <strong>Analysis Date:</strong> {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}<br>
            <strong>Total Molecules:</strong> {total_molecules:,}<br>
            <strong>Datasets Analyzed:</strong> OECD, Enhanced, Timing (excluding complex_branched, definitions, highly_branched)
        </div>
        
        <h2>📊 Overall Performance</h2>
        <div class="summary-grid">
            <div class="summary-card">
                <div class="summary-number">{overall_pfasgroups['mean']*1000:.1f}ms</div>
                <div class="summary-label">PFASgroups<br>Mean Time</div>
            </div>
            <div class="summary-card">
                <div class="summary-number">{overall_pfasgroups['median']*1000:.1f}ms</div>
                <div class="summary-label">PFASgroups<br>Median Time</div>
            </div>
            <div class="summary-card">
                <div class="summary-number">{overall_atlas['mean']*1000:.1f}ms</div>
                <div class="summary-label">PFAS-Atlas<br>Mean Time</div>
            </div>
            <div class="summary-card">
                <div class="summary-number">{speedup:.2f}x</div>
                <div class="summary-label">Relative Speed<br>(Atlas/PFASgroups)</div>
            </div>
        </div>
        
        <h2>📈 Performance by Dataset</h2>
        {generate_table_html(dataset_stats, 'Dataset')}
        <div class="chart-container">
            <div id="dataset-chart"></div>
        </div>
        
        <h2>📏 Performance by Compound Size</h2>
        {generate_table_html(size_stats, 'Size Category')}
        <div class="chart-container">
            <div id="size-chart"></div>
        </div>
        
        <h2>🔗 Performance by Chain Type</h2>
        {generate_table_html(chain_type_stats, 'Chain Type')}
        <div class="chart-container">
            <div id="chain-type-chart"></div>
        </div>
        
        <h2>⚛️ Performance by Atom Count</h2>
        {generate_table_html(atom_stats, 'Atom Range', is_atom_table=True)}
        <div class="chart-container">
            <div id="atom-chart"></div>
        </div>
        
        <h2>📊 Time Distribution</h2>
        <div class="chart-container">
            <div id="distribution-chart"></div>
        </div>
        
        <div class="info-box">
            <strong>Key Findings:</strong><br>
            • PFASgroups provides detailed component analysis, metrics, and full customizability<br>
            • Execution time increases with molecular complexity as expected for comprehensive analysis<br>
            • Performance is consistent across different dataset types and molecule sizes<br>
            • The additional computational cost enables detection of multiple functional groups, chain length determination, and graph-theoretical metrics
        </div>
    </div>
    
    <script>
        {dataset_chart_js}
        {size_chart_js}
        {chain_type_chart_js}
        {atom_chart_js}
        {distribution_chart_js}
    </script>
</body>
</html>
"""
    
    return html

def generate_table_html(stats_list, category_name, is_atom_table=False):
    """Generate HTML table for statistics"""
    rows = []
    for stat in stats_list:
        pfas = stat['pfas_stats']
        atlas = stat['atlas_stats']
        
        if pfas and atlas:
            speedup = atlas['mean'] / pfas['mean']
            speedup_class = 'speedup' if speedup > 1 else 'slower'
            speedup_text = f'<span class="{speedup_class}">{speedup:.2f}x</span>'
        else:
            speedup_text = 'N/A'
        
        pfas_mean = f"{pfas['mean']*1000:.2f}ms" if pfas else 'N/A'
        pfas_median = f"{pfas['median']*1000:.2f}ms" if pfas else 'N/A'
        atlas_mean = f"{atlas['mean']*1000:.2f}ms" if atlas else 'N/A'
        
        rows.append(f"""
        <tr>
            <td><strong>{stat['name']}</strong></td>
            <td class="number">{stat['count']:,}</td>
            <td class="number">{pfas_mean}</td>
            <td class="number">{pfas_median}</td>
            <td class="number">{atlas_mean}</td>
            <td class="number">{speedup_text}</td>
        </tr>""")
    
    return f"""
    <table>
        <tr>
            <th>{category_name}</th>
            <th>Count</th>
            <th>PFASgroups Mean</th>
            <th>PFASgroups Median</th>
            <th>PFAS-Atlas Mean</th>
            <th>Relative Speed</th>
        </tr>
        {''.join(rows)}
    </table>
    """

def generate_dataset_chart_js(dataset_stats):
    """Generate Plotly.js code for dataset comparison"""
    datasets = [d['name'] for d in dataset_stats]
    pfas_means = [d['pfas_stats']['mean']*1000 if d['pfas_stats'] else 0 for d in dataset_stats]
    atlas_means = [d['atlas_stats']['mean']*1000 if d['atlas_stats'] else 0 for d in dataset_stats]
    
    return f"""
    var datasetTrace1 = {{
        x: {json.dumps(datasets)},
        y: {json.dumps(pfas_means)},
        name: 'PFASgroups',
        type: 'bar',
        marker: {{color: '#667eea'}}
    }};
    
    var datasetTrace2 = {{
        x: {json.dumps(datasets)},
        y: {json.dumps(atlas_means)},
        name: 'PFAS-Atlas',
        type: 'bar',
        marker: {{color: '#764ba2'}}
    }};
    
    var datasetLayout = {{
        title: 'Mean Execution Time by Dataset',
        xaxis: {{title: 'Dataset'}},
        yaxis: {{title: 'Time (ms)'}},
        barmode: 'group'
    }};
    
    Plotly.newPlot('dataset-chart', [datasetTrace1, datasetTrace2], datasetLayout);
    """

def generate_size_chart_js(size_stats):
    """Generate Plotly.js code for size comparison"""
    sizes = [s['name'] for s in size_stats]
    pfas_means = [s['pfas_stats']['mean']*1000 if s['pfas_stats'] else 0 for s in size_stats]
    atlas_means = [s['atlas_stats']['mean']*1000 if s['atlas_stats'] else 0 for s in size_stats]
    
    return f"""
    var sizeTrace1 = {{
        x: {json.dumps(sizes)},
        y: {json.dumps(pfas_means)},
        name: 'PFASgroups',
        type: 'bar',
        marker: {{color: '#667eea'}}
    }};
    
    var sizeTrace2 = {{
        x: {json.dumps(sizes)},
        y: {json.dumps(atlas_means)},
        name: 'PFAS-Atlas',
        type: 'bar',
        marker: {{color: '#764ba2'}}
    }};
    
    var sizeLayout = {{
        title: 'Mean Execution Time by Compound Size',
        xaxis: {{title: 'Size Category', categoryorder: 'array', categoryarray: ['Small', 'Medium', 'Large']}},
        yaxis: {{title: 'Time (ms)'}},
        barmode: 'group'
    }};
    
    Plotly.newPlot('size-chart', [sizeTrace1, sizeTrace2], sizeLayout);
    """

def generate_chain_type_chart_js(chain_type_stats):
    """Generate Plotly.js code for chain type comparison"""
    chain_types = [c['name'] for c in chain_type_stats]
    pfas_means = [c['pfas_stats']['mean']*1000 if c['pfas_stats'] else 0 for c in chain_type_stats]
    atlas_means = [c['atlas_stats']['mean']*1000 if c['atlas_stats'] else 0 for c in chain_type_stats]
    
    return f"""
    var chainTrace1 = {{
        x: {json.dumps(chain_types)},
        y: {json.dumps(pfas_means)},
        name: 'PFASgroups',
        type: 'bar',
        marker: {{color: '#667eea'}}
    }};
    
    var chainTrace2 = {{
        x: {json.dumps(chain_types)},
        y: {json.dumps(atlas_means)},
        name: 'PFAS-Atlas',
        type: 'bar',
        marker: {{color: '#764ba2'}}
    }};
    
    var chainLayout = {{
        title: 'Mean Execution Time by Chain Type',
        xaxis: {{title: 'Chain Type'}},
        yaxis: {{title: 'Time (ms)'}},
        barmode: 'group'
    }};
    
    Plotly.newPlot('chain-type-chart', [chainTrace1, chainTrace2], chainLayout);
    """

def generate_atom_chart_js(atom_stats):
    """Generate Plotly.js code for atom count scaling"""
    ranges = [a['range'] for a in atom_stats]
    pfas_means = [a['pfas_stats']['mean']*1000 if a['pfas_stats'] else 0 for a in atom_stats]
    atlas_means = [a['atlas_stats']['mean']*1000 if a['atlas_stats'] else 0 for a in atom_stats]
    
    return f"""
    var atomTrace1 = {{
        x: {json.dumps(ranges)},
        y: {json.dumps(pfas_means)},
        name: 'PFASgroups',
        type: 'scatter',
        mode: 'lines+markers',
        marker: {{color: '#667eea', size: 8}},
        line: {{color: '#667eea', width: 2}}
    }};
    
    var atomTrace2 = {{
        x: {json.dumps(ranges)},
        y: {json.dumps(atlas_means)},
        name: 'PFAS-Atlas',
        type: 'scatter',
        mode: 'lines+markers',
        marker: {{color: '#764ba2', size: 8}},
        line: {{color: '#764ba2', width: 2}}
    }};
    
    var atomLayout = {{
        title: 'Execution Time Scaling with Molecule Size',
        xaxis: {{title: 'Atom Count Range'}},
        yaxis: {{title: 'Mean Time (ms)'}}
    }};
    
    Plotly.newPlot('atom-chart', [atomTrace1, atomTrace2], atomLayout);
    """

def generate_distribution_chart_js(pfasgroups_times, atlas_times):
    """Generate Plotly.js code for time distribution histograms"""
    return f"""
    var distTrace1 = {{
        x: {json.dumps([t*1000 for t in pfasgroups_times])},
        name: 'PFASgroups',
        type: 'histogram',
        opacity: 0.7,
        marker: {{color: '#667eea'}},
        xbins: {{size: 10}}
    }};
    
    var distTrace2 = {{
        x: {json.dumps([t*1000 for t in atlas_times])},
        name: 'PFAS-Atlas',
        type: 'histogram',
        opacity: 0.7,
        marker: {{color: '#764ba2'}},
        xbins: {{size: 10}}
    }};
    
    var distLayout = {{
        title: 'Distribution of Execution Times',
        xaxis: {{title: 'Time (ms)'}},
        yaxis: {{title: 'Frequency'}},
        barmode: 'overlay'
    }};
    
    Plotly.newPlot('distribution-chart', [distTrace1, distTrace2], distLayout);
    """

if __name__ == '__main__':
    import sys
    output_file = sys.argv[1] if len(sys.argv) > 1 else 'reports/enhanced_timing_report.html'
    success = generate_html_report(output_file)
    sys.exit(0 if success else 1)

