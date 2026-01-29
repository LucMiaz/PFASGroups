#!/usr/bin/env python3
"""
Generate HTML timing report for PFAS benchmark results
Analyzes execution times for PFASgroups and PFAS-Atlas across datasets
"""

import json
import sqlite3
from pathlib import Path
from datetime import datetime
import statistics
from collections import defaultdict

def generate_html_report(output_file='benchmark_timing_report.html'):
    """Generate comprehensive timing analysis HTML report"""
    
    # Connect to database
    db_path = Path(__file__).parent / 'review-app' / 'database' / 'pfas-benchmark.db'
    
    if not db_path.exists():
        print(f"❌ Database not found at {db_path}")
        print("   Please run the import script first: node review-app/scripts/import-benchmark-data.js")
        return False
    
    conn = sqlite3.connect(str(db_path))
    cursor = conn.cursor()
    
    print("📊 Analyzing timing data...")
    
    # Get all molecules with timing data
    cursor.execute("""
        SELECT 
            m.dataset_type,
            m.num_atoms,
            m.smiles,
            p.execution_time as pfasgroups_time,
            a.execution_time as atlas_time
        FROM molecules m
        LEFT JOIN pfasgroups_results p ON m.id = p.molecule_id
        LEFT JOIN atlas_results a ON m.id = a.molecule_id
        WHERE p.execution_time IS NOT NULL OR a.execution_time IS NOT NULL
    """)
    
    rows = cursor.fetchall()
    conn.close()
    
    if not rows:
        print("⚠️  No timing data found in database")
        return False
    
    print(f"✓ Found {len(rows)} molecules with timing data")
    
    # Organize data
    data_by_dataset = defaultdict(lambda: {'pfasgroups': [], 'atlas': [], 'atoms': []})
    all_pfasgroups = []
    all_atlas = []
    data_by_atoms = defaultdict(lambda: {'pfasgroups': [], 'atlas': []})
    
    for dataset, num_atoms, smiles, pfas_time, atlas_time in rows:
        if pfas_time is not None:
            data_by_dataset[dataset]['pfasgroups'].append(pfas_time)
            all_pfasgroups.append(pfas_time)
            if num_atoms:
                data_by_atoms[num_atoms]['pfasgroups'].append(pfas_time)
        
        if atlas_time is not None:
            data_by_dataset[dataset]['atlas'].append(atlas_time)
            all_atlas.append(atlas_time)
            if num_atoms:
                data_by_atoms[num_atoms]['atlas'].append(atlas_time)
        
        if num_atoms:
            data_by_dataset[dataset]['atoms'].append(num_atoms)
    
    # Generate HTML report
    html = generate_html_content(data_by_dataset, all_pfasgroups, all_atlas, data_by_atoms, len(rows))
    
    output_path = Path(__file__).parent / output_file
    with open(output_path, 'w', encoding='utf-8') as f:
        f.write(html)
    
    print(f"✅ Report generated: {output_path}")
    return True

def generate_html_content(data_by_dataset, all_pfasgroups, all_atlas, data_by_atoms, total_molecules):
    """Generate HTML content with embedded charts"""
    
    # Calculate statistics
    def calc_stats(times):
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
    
    overall_pfasgroups = calc_stats(all_pfasgroups)
    overall_atlas = calc_stats(all_atlas)
    
    # Build dataset statistics table
    dataset_rows = []
    for dataset in sorted(data_by_dataset.keys()):
        d = data_by_dataset[dataset]
        pfas_stats = calc_stats(d['pfasgroups'])
        atlas_stats = calc_stats(d['atlas'])
        
        dataset_rows.append({
            'name': dataset,
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
    
    # Generate HTML
    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>PFAS Benchmark Timing Report</title>
    <script src="https://cdn.plot.ly/plotly-2.27.0.min.js"></script>
    <style>
        * {{
            margin: 0;
            padding: 0;
            box-sizing: border-box;
        }}
        
        body {{
            font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, Oxygen, Ubuntu, Cantarell, sans-serif;
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            padding: 20px;
            min-height: 100vh;
        }}
        
        .container {{
            max-width: 1400px;
            margin: 0 auto;
            background: white;
            border-radius: 15px;
            box-shadow: 0 20px 60px rgba(0, 0, 0, 0.3);
            overflow: hidden;
        }}
        
        .header {{
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            color: white;
            padding: 40px;
            text-align: center;
        }}
        
        .header h1 {{
            font-size: 2.5rem;
            margin-bottom: 10px;
        }}
        
        .header p {{
            font-size: 1.1rem;
            opacity: 0.9;
        }}
        
        .content {{
            padding: 40px;
        }}
        
        .section {{
            margin-bottom: 50px;
        }}
        
        .section h2 {{
            font-size: 1.8rem;
            color: #333;
            margin-bottom: 20px;
            padding-bottom: 10px;
            border-bottom: 3px solid #667eea;
        }}
        
        .stats-grid {{
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(250px, 1fr));
            gap: 20px;
            margin-bottom: 30px;
        }}
        
        .stat-card {{
            background: #f8f9fa;
            border-radius: 10px;
            padding: 20px;
            border-left: 4px solid #667eea;
        }}
        
        .stat-card.atlas {{
            border-left-color: #764ba2;
        }}
        
        .stat-card h3 {{
            font-size: 1rem;
            color: #666;
            margin-bottom: 10px;
            text-transform: uppercase;
            letter-spacing: 0.5px;
        }}
        
        .stat-card .value {{
            font-size: 2rem;
            font-weight: bold;
            color: #333;
        }}
        
        .stat-card .unit {{
            font-size: 1rem;
            color: #999;
            margin-left: 5px;
        }}
        
        .stat-card .detail {{
            font-size: 0.9rem;
            color: #666;
            margin-top: 8px;
        }}
        
        table {{
            width: 100%;
            border-collapse: collapse;
            margin-bottom: 30px;
            background: white;
        }}
        
        thead {{
            background: #667eea;
            color: white;
        }}
        
        th {{
            padding: 15px;
            text-align: left;
            font-weight: 600;
        }}
        
        td {{
            padding: 12px 15px;
            border-bottom: 1px solid #eee;
        }}
        
        tbody tr:hover {{
            background: #f8f9fa;
        }}
        
        .number {{
            font-family: 'Courier New', monospace;
            text-align: right;
        }}
        
        .chart-container {{
            background: #f8f9fa;
            border-radius: 10px;
            padding: 20px;
            margin-bottom: 30px;
        }}
        
        .footer {{
            text-align: center;
            padding: 20px;
            color: #666;
            font-size: 0.9rem;
            border-top: 1px solid #eee;
        }}
        
        .speedup {{
            color: #10b981;
            font-weight: bold;
        }}
        
        .slower {{
            color: #ef4444;
            font-weight: bold;
        }}
    </style>
</head>
<body>
    <div class="container">
        <div class="header">
            <h1>⚡ PFAS Benchmark Timing Report</h1>
            <p>Performance Analysis: PFASgroups vs PFAS-Atlas</p>
            <p style="margin-top: 10px; font-size: 0.95rem;">Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}</p>
        </div>
        
        <div class="content">
            <!-- Overall Statistics -->
            <div class="section">
                <h2>📊 Overall Performance</h2>
                <div class="stats-grid">
                    <div class="stat-card">
                        <h3>PFASgroups Mean Time</h3>
                        <div class="value">{overall_pfasgroups['mean']*1000:.2f}<span class="unit">ms</span></div>
                        <div class="detail">Median: {overall_pfasgroups['median']*1000:.2f}ms</div>
                        <div class="detail">σ: {overall_pfasgroups['stdev']*1000:.2f}ms</div>
                    </div>
                    
                    <div class="stat-card">
                        <h3>PFASgroups Range</h3>
                        <div class="value">{overall_pfasgroups['min']*1000:.2f}<span class="unit">ms</span></div>
                        <div class="detail">Max: {overall_pfasgroups['max']*1000:.2f}ms</div>
                        <div class="detail">Molecules: {overall_pfasgroups['count']:,}</div>
                    </div>
                    
                    {f'''<div class="stat-card atlas">
                        <h3>PFAS-Atlas Mean Time</h3>
                        <div class="value">{overall_atlas['mean']*1000:.2f}<span class="unit">ms</span></div>
                        <div class="detail">Median: {overall_atlas['median']*1000:.2f}ms</div>
                        <div class="detail">σ: {overall_atlas['stdev']*1000:.2f}ms</div>
                    </div>
                    
                    <div class="stat-card atlas">
                        <h3>Speedup Factor</h3>
                        <div class="value {'speedup' if overall_atlas['mean'] > overall_pfasgroups['mean'] else 'slower'}">{overall_atlas['mean']/overall_pfasgroups['mean']:.2f}<span class="unit">x</span></div>
                        <div class="detail">{'PFASgroups is faster' if overall_atlas['mean'] > overall_pfasgroups['mean'] else 'PFAS-Atlas is faster'}</div>
                        <div class="detail">Molecules: {overall_atlas['count']:,}</div>
                    </div>''' if overall_atlas else ''}
                </div>
            </div>
            
            <!-- Performance by Dataset -->
            <div class="section">
                <h2>📁 Performance by Dataset</h2>
                <table>
                    <thead>
                        <tr>
                            <th>Dataset</th>
                            <th>Count</th>
                            <th>PFASgroups Mean</th>
                            <th>PFASgroups Median</th>
                            <th>Atlas Mean</th>
                            <th>Speedup</th>
                        </tr>
                    </thead>
                    <tbody>
                        {generate_dataset_rows(dataset_rows)}
                    </tbody>
                </table>
                
                <div class="chart-container">
                    <div id="datasetChart"></div>
                </div>
            </div>
            
            <!-- Performance by Molecule Size -->
            <div class="section">
                <h2>🔬 Performance by Molecule Size</h2>
                <table>
                    <thead>
                        <tr>
                            <th>Atom Count Range</th>
                            <th>PFASgroups Mean</th>
                            <th>PFASgroups Median</th>
                            <th>Atlas Mean</th>
                            <th>Speedup</th>
                        </tr>
                    </thead>
                    <tbody>
                        {generate_atom_rows(atom_stats)}
                    </tbody>
                </table>
                
                <div class="chart-container">
                    <div id="atomChart"></div>
                </div>
            </div>
            
            <!-- Distribution Charts -->
            <div class="section">
                <h2>📈 Time Distribution</h2>
                <div class="chart-container">
                    <div id="distributionChart"></div>
                </div>
            </div>
        </div>
        
        <div class="footer">
            <p>PFAS Benchmark Timing Report • Total Molecules: {total_molecules:,}</p>
            <p>Generated by PFASgroups Benchmark Suite</p>
        </div>
    </div>
    
    <script>
        // Dataset comparison chart
        {generate_dataset_chart_js(dataset_rows)}
        
        // Atom count chart
        {generate_atom_chart_js(atom_stats)}
        
        // Distribution chart
        {generate_distribution_chart_js(all_pfasgroups, all_atlas)}
    </script>
</body>
</html>"""
    
    return html

def generate_dataset_rows(dataset_rows):
    """Generate HTML table rows for dataset statistics"""
    rows = []
    for d in dataset_rows:
        pfas = d['pfas_stats']
        atlas = d['atlas_stats']
        
        if pfas and atlas:
            speedup = atlas['mean'] / pfas['mean']
            speedup_class = 'speedup' if speedup > 1 else 'slower'
            speedup_text = f'<span class="{speedup_class}">{speedup:.2f}x</span>'
        else:
            speedup_text = 'N/A'
        
        rows.append(f"""
        <tr>
            <td><strong>{d['name']}</strong></td>
            <td class="number">{d['count']:,}</td>
            <td class="number">{pfas['mean']*1000:.2f}ms</td>
            <td class="number">{pfas['median']*1000:.2f}ms</td>
            <td class="number">{atlas['mean']*1000:.2f}ms if atlas else 'N/A'}</td>
            <td class="number">{speedup_text}</td>
        </tr>""")
    
    return '\n'.join(rows)

def generate_atom_rows(atom_stats):
    """Generate HTML table rows for atom count statistics"""
    rows = []
    for a in atom_stats:
        pfas = a['pfas_stats']
        atlas = a['atlas_stats']
        
        if pfas and atlas:
            speedup = atlas['mean'] / pfas['mean']
            speedup_class = 'speedup' if speedup > 1 else 'slower'
            speedup_text = f'<span class="{speedup_class}">{speedup:.2f}x</span>'
        else:
            speedup_text = 'N/A'
        
        rows.append(f"""
        <tr>
            <td><strong>{a['range']} atoms</strong></td>
            <td class="number">{pfas['mean']*1000:.2f}ms if pfas else 'N/A'}</td>
            <td class="number">{pfas['median']*1000:.2f}ms if pfas else 'N/A'}</td>
            <td class="number">{atlas['mean']*1000:.2f}ms if atlas else 'N/A'}</td>
            <td class="number">{speedup_text}</td>
        </tr>""")
    
    return '\n'.join(rows)

def generate_dataset_chart_js(dataset_rows):
    """Generate Plotly.js code for dataset comparison chart"""
    datasets = [d['name'] for d in dataset_rows]
    pfas_means = [d['pfas_stats']['mean']*1000 if d['pfas_stats'] else 0 for d in dataset_rows]
    atlas_means = [d['atlas_stats']['mean']*1000 if d['atlas_stats'] else 0 for d in dataset_rows]
    
    return f"""
        var datasetTrace1 = {{
            x: {json.dumps(datasets)},
            y: {json.dumps(pfas_means)},
            name: 'PFASgroups',
            type: 'bar',
            marker: {{ color: '#667eea' }}
        }};
        
        var datasetTrace2 = {{
            x: {json.dumps(datasets)},
            y: {json.dumps(atlas_means)},
            name: 'PFAS-Atlas',
            type: 'bar',
            marker: {{ color: '#764ba2' }}
        }};
        
        var datasetLayout = {{
            title: 'Mean Execution Time by Dataset',
            xaxis: {{ title: 'Dataset' }},
            yaxis: {{ title: 'Time (ms)' }},
            barmode: 'group',
            height: 400
        }};
        
        Plotly.newPlot('datasetChart', [datasetTrace1, datasetTrace2], datasetLayout);
    """

def generate_atom_chart_js(atom_stats):
    """Generate Plotly.js code for atom count chart"""
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
            line: {{ color: '#667eea', width: 3 }},
            marker: {{ size: 8 }}
        }};
        
        var atomTrace2 = {{
            x: {json.dumps(ranges)},
            y: {json.dumps(atlas_means)},
            name: 'PFAS-Atlas',
            type: 'scatter',
            mode: 'lines+markers',
            line: {{ color: '#764ba2', width: 3 }},
            marker: {{ size: 8 }}
        }};
        
        var atomLayout = {{
            title: 'Execution Time vs Molecule Size',
            xaxis: {{ title: 'Number of Atoms' }},
            yaxis: {{ title: 'Time (ms)' }},
            height: 400
        }};
        
        Plotly.newPlot('atomChart', [atomTrace1, atomTrace2], atomLayout);
    """

def generate_distribution_chart_js(pfasgroups_times, atlas_times):
    """Generate Plotly.js code for distribution histogram"""
    return f"""
        var distTrace1 = {{
            x: {json.dumps([t*1000 for t in pfasgroups_times])},
            name: 'PFASgroups',
            type: 'histogram',
            opacity: 0.7,
            marker: {{ color: '#667eea' }},
            xbins: {{ size: 1 }}
        }};
        
        var distTrace2 = {{
            x: {json.dumps([t*1000 for t in atlas_times]) if atlas_times else '[]'},
            name: 'PFAS-Atlas',
            type: 'histogram',
            opacity: 0.7,
            marker: {{ color: '#764ba2' }},
            xbins: {{ size: 1 }}
        }};
        
        var distLayout = {{
            title: 'Execution Time Distribution',
            xaxis: {{ title: 'Time (ms)' }},
            yaxis: {{ title: 'Frequency' }},
            barmode: 'overlay',
            height: 400
        }};
        
        Plotly.newPlot('distributionChart', [{('distTrace1, distTrace2' if atlas_times else 'distTrace1')}], distLayout);
    """

if __name__ == '__main__':
    import sys
    output_file = sys.argv[1] if len(sys.argv) > 1 else 'benchmark_timing_report.html'
    success = generate_html_report(output_file)
    sys.exit(0 if success else 1)
