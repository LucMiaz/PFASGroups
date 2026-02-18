#!/usr/bin/env python3
"""
Export timing analysis figures from HTML report to PNG for LaTeX inclusion
Requires: plotly, kaleido
Install: pip install plotly kaleido
"""

import json
import sqlite3
from pathlib import Path
import statistics
from collections import defaultdict

try:
    import plotly.graph_objects as go
    from plotly.subplots import make_subplots
    PLOTLY_AVAILABLE = True
except ImportError:
    print("⚠️  Plotly not installed. Install with: pip install plotly kaleido")
    PLOTLY_AVAILABLE = False

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

def export_timing_figures():
    """Export timing analysis figures as PNG files"""
    
    if not PLOTLY_AVAILABLE:
        return False
    
    # Connect to database
    db_path = Path(__file__).parent / 'review-app' / 'database' / 'pfas_benchmark.db'
    
    if not db_path.exists():
        print(f"❌ Database not found at {db_path}")
        return False
    
    conn = sqlite3.connect(str(db_path))
    cursor = conn.cursor()
    
    print("📊 Loading timing data...")
    
    # Get all molecules with timing data
    cursor.execute("""
        SELECT 
            m.id,
            m.dataset_type,
            m.num_atoms,
            m.group_id,
            m.generation_type,
            p.execution_time as HalogenGroups_time,
            p.detected_groups,
            a.execution_time as atlas_time
        FROM molecules m
        LEFT JOIN HalogenGroups_results p ON m.id = p.molecule_id
        LEFT JOIN atlas_results a ON m.id = a.molecule_id
        WHERE (p.execution_time IS NOT NULL OR a.execution_time IS NOT NULL)
          AND m.dataset_type NOT IN ('complex_branched', 'definitions', 'non_fluorinated', 'highly_branched')
    """)
    
    rows = cursor.fetchall()
    conn.close()
    
    if not rows:
        print("⚠️  No timing data found")
        return False
    
    print(f"✓ Found {len(rows)} molecules")
    
    # Organize data
    data_by_dataset = defaultdict(lambda: {'HalogenGroups': [], 'atlas': []})
    data_by_size = defaultdict(lambda: {'HalogenGroups': [], 'atlas': []})
    data_by_chain_type = defaultdict(lambda: {'HalogenGroups': [], 'atlas': []})
    data_by_atoms = defaultdict(lambda: {'HalogenGroups': [], 'atlas': []})
    
    all_HalogenGroups = []
    all_atlas = []
    
    for mol_id, dataset, num_atoms, group_id, gen_type, pfas_time, detected_groups, atlas_time in rows:
        # Size category
        if num_atoms:
            if num_atoms < 20:
                size_cat = 'small'
            elif num_atoms <= 40:
                size_cat = 'medium'
            else:
                size_cat = 'large'
        else:
            size_cat = 'unknown'
        
        # Chain type
        chain_type = 'standard'
        if gen_type and 'cyclic' in str(gen_type).lower():
            chain_type = 'cyclic'
        elif group_id and group_id >= 60:
            chain_type = 'fluorotelomer'
        
        # Store times
        if pfas_time is not None:
            data_by_dataset[dataset]['HalogenGroups'].append(pfas_time)
            data_by_size[size_cat]['HalogenGroups'].append(pfas_time)
            data_by_chain_type[chain_type]['HalogenGroups'].append(pfas_time)
            all_HalogenGroups.append(pfas_time)
            if num_atoms:
                data_by_atoms[num_atoms]['HalogenGroups'].append(pfas_time)
        
        if atlas_time is not None:
            data_by_dataset[dataset]['atlas'].append(atlas_time)
            data_by_size[size_cat]['atlas'].append(atlas_time)
            data_by_chain_type[chain_type]['atlas'].append(atlas_time)
            all_atlas.append(atlas_time)
            if num_atoms:
                data_by_atoms[num_atoms]['atlas'].append(atlas_time)
    
    output_dir = Path(__file__).parent
    
    # Figure 1: Dataset comparison
    print("📈 Generating dataset comparison figure...")
    dataset_stats = []
    for dataset in sorted(data_by_dataset.keys()):
        d = data_by_dataset[dataset]
        pfas_stats = calc_stats(d['HalogenGroups'])
        atlas_stats = calc_stats(d['atlas'])
        if pfas_stats and atlas_stats:
            dataset_stats.append({
                'name': dataset,
                'pfas_mean': pfas_stats['mean'] * 1000,
                'atlas_mean': atlas_stats['mean'] * 1000
            })
    
    fig1 = go.Figure()
    fig1.add_trace(go.Bar(
        x=[d['name'] for d in dataset_stats],
        y=[d['pfas_mean'] for d in dataset_stats],
        name='HalogenGroups',
        marker_color='#667eea'
    ))
    fig1.add_trace(go.Bar(
        x=[d['name'] for d in dataset_stats],
        y=[d['atlas_mean'] for d in dataset_stats],
        name='PFAS-Atlas',
        marker_color='#764ba2'
    ))
    fig1.update_layout(
        xaxis_title='Dataset',
        yaxis_title='Mean Time (ms)',
        barmode='group',
        template='plotly_white',
        width=800,
        height=500,
        font=dict(size=14)
    )
    fig1.write_image(str(output_dir / 'enhanced_timing_report_dataset.png'), scale=2)
    print(f"✓ Saved: enhanced_timing_report_dataset.png")
    
    # Figure 2: Size comparison
    print("📈 Generating size comparison figure...")
    size_stats = []
    for size in ['small', 'medium', 'large']:
        if size in data_by_size:
            d = data_by_size[size]
            pfas_stats = calc_stats(d['HalogenGroups'])
            atlas_stats = calc_stats(d['atlas'])
            if pfas_stats and atlas_stats:
                size_stats.append({
                    'name': size.capitalize(),
                    'pfas_mean': pfas_stats['mean'] * 1000,
                    'atlas_mean': atlas_stats['mean'] * 1000
                })
    
    fig2 = go.Figure()
    fig2.add_trace(go.Bar(
        x=[s['name'] for s in size_stats],
        y=[s['pfas_mean'] for s in size_stats],
        name='HalogenGroups',
        marker_color='#667eea'
    ))
    fig2.add_trace(go.Bar(
        x=[s['name'] for s in size_stats],
        y=[s['atlas_mean'] for s in size_stats],
        name='PFAS-Atlas',
        marker_color='#764ba2'
    ))
    fig2.update_layout(
        xaxis_title='Size Category',
        yaxis_title='Mean Time (ms)',
        barmode='group',
        template='plotly_white',
        width=800,
        height=500,
        font=dict(size=14)
    )
    fig2.write_image(str(output_dir / 'enhanced_timing_report_size.png'), scale=2)
    print(f"✓ Saved: enhanced_timing_report_size.png")
    
    # Figure 3: Chain type comparison
    print("📈 Generating chain type comparison figure...")
    chain_stats = []
    for chain_type in sorted(data_by_chain_type.keys()):
        d = data_by_chain_type[chain_type]
        pfas_stats = calc_stats(d['HalogenGroups'])
        atlas_stats = calc_stats(d['atlas'])
        if pfas_stats and atlas_stats:
            chain_stats.append({
                'name': chain_type.capitalize(),
                'pfas_mean': pfas_stats['mean'] * 1000,
                'atlas_mean': atlas_stats['mean'] * 1000
            })
    
    fig3 = go.Figure()
    fig3.add_trace(go.Bar(
        x=[c['name'] for c in chain_stats],
        y=[c['pfas_mean'] for c in chain_stats],
        name='HalogenGroups',
        marker_color='#667eea'
    ))
    fig3.add_trace(go.Bar(
        x=[c['name'] for c in chain_stats],
        y=[c['atlas_mean'] for c in chain_stats],
        name='PFAS-Atlas',
        marker_color='#764ba2'
    ))
    fig3.update_layout(
        xaxis_title='Chain Type',
        yaxis_title='Mean Time (ms)',
        barmode='group',
        template='plotly_white',
        width=800,
        height=500,
        font=dict(size=14)
    )
    fig3.write_image(str(output_dir / 'enhanced_timing_report_chain.png'), scale=2)
    print(f"✓ Saved: enhanced_timing_report_chain.png")
    
    # Figure 4: Atom count scaling
    print("📈 Generating atom scaling figure...")
    atom_bins = [0, 10, 20, 30, 40, 50, 100, 200]
    atom_stats = []
    for i in range(len(atom_bins) - 1):
        bin_start = atom_bins[i]
        bin_end = atom_bins[i + 1]
        
        pfas_times = []
        atlas_times = []
        
        for atoms, times_dict in data_by_atoms.items():
            if bin_start < atoms <= bin_end:
                pfas_times.extend(times_dict['HalogenGroups'])
                atlas_times.extend(times_dict['atlas'])
        
        if pfas_times and atlas_times:
            pfas_stats = calc_stats(pfas_times)
            atlas_stats = calc_stats(atlas_times)
            atom_stats.append({
                'range': f"{bin_start+1}-{bin_end}",
                'pfas_mean': pfas_stats['mean'] * 1000,
                'atlas_mean': atlas_stats['mean'] * 1000
            })
    
    fig4 = go.Figure()
    fig4.add_trace(go.Scatter(
        x=[a['range'] for a in atom_stats],
        y=[a['pfas_mean'] for a in atom_stats],
        name='HalogenGroups',
        mode='lines+markers',
        marker=dict(color='#667eea', size=8),
        line=dict(color='#667eea', width=2)
    ))
    fig4.add_trace(go.Scatter(
        x=[a['range'] for a in atom_stats],
        y=[a['atlas_mean'] for a in atom_stats],
        name='PFAS-Atlas',
        mode='lines+markers',
        marker=dict(color='#764ba2', size=8),
        line=dict(color='#764ba2', width=2)
    ))
    fig4.update_layout(
        xaxis_title='Atom Count Range',
        yaxis_title='Mean Time (ms)',
        template='plotly_white',
        width=900,
        height=500,
        font=dict(size=14)
    )
    fig4.write_image(str(output_dir / 'enhanced_timing_report_atoms.png'), scale=2)
    print(f"✓ Saved: enhanced_timing_report_atoms.png")
    
    # Figure 5: Distribution
    print("📈 Generating distribution figure...")
    fig5 = go.Figure()
    fig5.add_trace(go.Histogram(
        x=[t*1000 for t in all_HalogenGroups],
        name='HalogenGroups',
        opacity=0.7,
        marker_color='#667eea',
        xbins=dict(size=10)
    ))
    fig5.add_trace(go.Histogram(
        x=[t*1000 for t in all_atlas],
        name='PFAS-Atlas',
        opacity=0.7,
        marker_color='#764ba2',
        xbins=dict(size=10)
    ))
    fig5.update_layout(
        xaxis_title='Time (ms)',
        yaxis_title='Frequency',
        barmode='overlay',
        template='plotly_white',
        width=900,
        height=500,
        font=dict(size=14)
    )
    fig5.write_image(str(output_dir / 'enhanced_timing_report_distribution.png'), scale=2)
    print(f"✓ Saved: enhanced_timing_report_distribution.png")
    
    print(f"\n✅ All figures exported to: {output_dir}")
    return True

if __name__ == '__main__':
    import sys
    success = export_timing_figures()
    sys.exit(0 if success else 1)

