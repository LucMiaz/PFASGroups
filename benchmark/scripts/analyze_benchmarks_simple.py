"""
Comprehensive benchmark analysis without numpy dependency.
Analyzes timing models and generates summary statistics for LaTeX documents.
Generates exponential fit visualization using Plotly with PNG/SVG export.
"""

import json
import math
import glob
import os

try:
    import plotly.graph_objects as go
    import plotly.io as pio
    PLOTLY_AVAILABLE = True
except ImportError:
    PLOTLY_AVAILABLE = False
    print("⚠️  plotly not available - plots will not be generated")

print("="*80)
print("PFASGROUPS BENCHMARK ANALYSIS")
print("="*80)

# Load all benchmark data
print("\n1. Loading benchmark data...")

# Find the latest benchmark files
import glob

def get_latest_file(pattern):
    """Get the most recent file matching the pattern."""
    files = glob.glob(pattern)
    if not files:
        return None
    return max(files, key=os.path.getmtime) if os.path.exists else max(files)

# Timing benchmark
timing_file = get_latest_file('data/pfas_timing_benchmark_*.json')
if not timing_file:
    print("  Error: No timing benchmark file found")
    exit(1)
with open(timing_file, 'r') as f:
    timing_data = json.load(f)

# Enhanced benchmark
enhanced_file = get_latest_file('data/pfas_enhanced_benchmark_*.json')
if not enhanced_file:
    print("  Error: No enhanced benchmark file found")
    exit(1)
with open(enhanced_file, 'r') as f:
    enhanced_data = json.load(f)

# Complex branched
complex_file = get_latest_file('data/pfas_complex_branched_benchmark_*.json')
if not complex_file:
    print("  Error: No complex branched benchmark file found")
    exit(1)
with open(complex_file, 'r') as f:
    complex_data = json.load(f)

# Highly branched
hb_file = get_latest_file('data/pfas_highly_branched_benchmark_*.json')
if not hb_file:
    print("  Error: No highly branched benchmark file found")
    exit(1)
with open(hb_file, 'r') as f:
    highly_branched_data = json.load(f)

# Telomer validation
with open('data/telomer_validation_results.json', 'r') as f:
    telomer_data = json.load(f)

print(f"  Timing: {timing_file}")
print(f"  Enhanced: {enhanced_file}")
print(f"  Complex branched: {complex_file}")
print(f"  Highly branched: {hb_file}")
print(f"  Telomer validation: data/telomer_validation_results.json")
print(f"  Loaded: {len(timing_data)} timing, {len(enhanced_data)} enhanced, {len(complex_data)} complex, {len(highly_branched_data)} HB, {len(telomer_data.get('results', []))} telomer")

# Analyze timing model
print("\n2. Analyzing timing model (exponential fit)...")
timing_atoms = [d['num_atoms'] for d in timing_data]
timing_avg = [d['pfasgroups_time_avg'] for d in timing_data]

# Log-linear regression for exponential model: t = a * exp(b * n)
# ln(t) = ln(a) + b*n
log_timing = [math.log(t) for t in timing_avg]

# Calculate means
n = len(timing_atoms)
mean_x = sum(timing_atoms) / n
mean_y = sum(log_timing) / n

# Calculate slope (b) and intercept (ln(a))
numerator = sum((timing_atoms[i] - mean_x) * (log_timing[i] - mean_y) for i in range(n))
denominator = sum((timing_atoms[i] - mean_x)**2 for i in range(n))
b_fit = numerator / denominator
ln_a = mean_y - b_fit * mean_x
a_fit = math.exp(ln_a)

# Calculate R²
y_pred = [a_fit * math.exp(b_fit * x) for x in timing_atoms]
ss_res = sum((timing_avg[i] - y_pred[i])**2 for i in range(n))
mean_t = sum(timing_avg) / n
ss_tot = sum((t - mean_t)**2 for t in timing_avg)
r_squared = 1 - (ss_res / ss_tot)

# Calculate Pearson correlation
mean_timing_avg = sum(timing_avg) / n
numerator = sum((timing_atoms[i] - mean_x) * (timing_avg[i] - mean_timing_avg) for i in range(n))
denominator_x = math.sqrt(sum((x - mean_x)**2 for x in timing_atoms))
denominator_y = math.sqrt(sum((y - mean_timing_avg)**2 for y in timing_avg))
corr_coef = numerator / (denominator_x * denominator_y)

# Calculate statistics
median_timing = sorted(timing_avg)[len(timing_avg)//2]
std_timing = math.sqrt(sum((t - mean_t)**2 for t in timing_avg) / n)

print(f"\n  Exponential Model: t = a x exp(b x n)")
print(f"  Parameters:")
print(f"    a = {a_fit:.6e} seconds")
print(f"    b = {b_fit:.6f} atoms^-1 (alpha in paper)")
print(f"  Quality:")
print(f"    R^2 = {r_squared:.6f}")
print(f"    Pearson r = {corr_coef:.6f}")
print(f"  Statistics:")
print(f"    Mean time: {mean_t*1000:.2f} ms")
print(f"    Median time: {median_timing*1000:.2f} ms")
print(f"    Std time: {std_timing*1000:.2f} ms")

# Analyze enhanced benchmark
print("\n3. Analyzing enhanced benchmark...")
enhanced_detected = []
enhanced_unique_groups = set()
for mol in enhanced_data:
    detected = mol['pfasgroups_result']['detected_groups']
    enhanced_detected.append(len(detected) > 0)
    enhanced_unique_groups.update(detected)

detection_rate_enhanced = 100.0 * sum(enhanced_detected) / len(enhanced_detected)
print(f"  Total molecules: {len(enhanced_data)}")
print(f"  Detected: {sum(enhanced_detected)}")
print(f"  Detection rate: {detection_rate_enhanced:.1f}%")
print(f"  Unique groups found: {len(enhanced_unique_groups)}")

# Count group frequencies
group_counts = {}
for mol in enhanced_data:
    for gid in mol['pfasgroups_result']['detected_groups']:
        group_counts[gid] = group_counts.get(gid, 0) + 1

top_groups = sorted(group_counts.items(), key=lambda x: x[1], reverse=True)[:5]
print(f"  Top 5 groups:")
for gid, count in top_groups:
    print(f"    Group {gid}: {count} detections")

# Analyze complex branched
print("\n4. Analyzing complex branched benchmark...")
complex_detected = 0
complex_total = 0
complex_unique_groups = set()
for test in complex_data:
    for mol in test['molecules']:
        complex_total += 1
        detected = mol.get('pfasgroups_groups', [])
        if len(detected) > 0:
            complex_detected += 1
        complex_unique_groups.update(detected)

detection_rate_complex = 100.0 * complex_detected / complex_total if complex_total > 0 else 0
print(f"  Total molecules: {complex_total}")
print(f"  Detected: {complex_detected}")
print(f"  Detection rate: {detection_rate_complex:.1f}%")
print(f"  Unique groups found: {len(complex_unique_groups)}")

# Analyze highly branched
print("\n5. Analyzing highly branched benchmark...")
hb_total = highly_branched_data.get('summary', {}).get('total', 0)
hb_detected = highly_branched_data.get('summary', {}).get('pfasgroups_passed', 0)
hb_unique_groups = set()
for detail in highly_branched_data.get('details', []):
    # Group ID appears directly in detail
    group_id = detail.get('group_id')
    if group_id is not None:
        hb_unique_groups.add(group_id)

detection_rate_hb = 100.0 * hb_detected / hb_total if hb_total > 0 else 0
print(f"  Total molecules: {hb_total}")
print(f"  Detected: {hb_detected}")
print(f"  Detection rate: {detection_rate_hb:.1f}%")
print(f"  Unique groups found: {len(hb_unique_groups)}")

# Analyze telomer validation
print("\n6. Analyzing telomer validation...")
telomer_total = telomer_data.get('total_molecules', 0)
telomer_detected_count = telomer_data.get('telomer_detected', 0)
detection_rate_telomer = telomer_data.get('detection_rate', 0.0)

# Get unique groups from group_counts
telomer_groups = set()
telomer_group_counts = {}
for group_info in telomer_data.get('group_counts', []):
    gid = group_info.get('id')
    count = group_info.get('count', 0)
    if gid:
        telomer_groups.add(gid)
        telomer_group_counts[gid] = count

print(f"  Total fluorotelomers: {telomer_total}")
print(f"  Detected: {telomer_detected_count}")
print(f"  Detection rate: {detection_rate_telomer:.1f}%")
print(f"  Unique groups found: {len(telomer_groups)}")

# Show top telomer groups
top_telomer_groups = sorted(telomer_group_counts.items(), key=lambda x: x[1], reverse=True)
print(f"  Telomer group detections:")
for gid, count in top_telomer_groups:
    print(f"    Group {gid}: {count} ({100.0*count/telomer_total:.1f}%)")

# Generate LaTeX table rows
print("\n" + "="*80)
print("LATEX TABLE ROWS")
print("="*80)

print("\nTiming Model Parameters:")
print(f"All metrics & ${a_fit:.2e}$ & ${b_fit:.4f}$ & ${r_squared:.4f}$ & {mean_t*1000:.1f} & {median_timing*1000:.1f} \\\\")

print("\nBenchmark Detection Statistics:")
print(f"Enhanced & {len(enhanced_data)} & {sum(enhanced_detected)} & {detection_rate_enhanced:.1f}\\% & {len(enhanced_unique_groups)} \\\\")
print(f"Complex Branched & {complex_total} & {complex_detected} & {detection_rate_complex:.1f}\\% & {len(complex_unique_groups)} \\\\")
print(f"Highly Branched & {hb_total} & {hb_detected} & {detection_rate_hb:.1f}\\% & {len(hb_unique_groups)} \\\\")
print(f"Telomer Validation & {telomer_total} & {telomer_detected_count} & {detection_rate_telomer:.1f}\\% & {len(telomer_groups)} \\\\")

# Save summary to JSON
summary = {
    'timing_model': {
        'model': 't = a × exp(b × n)',
        'a': a_fit,
        'b': b_fit,
        'r_squared': r_squared,
        'correlation': corr_coef,
        'mean_time_ms': mean_t * 1000,
        'median_time_ms': median_timing * 1000,
        'std_time_ms': std_timing * 1000,
        'n_molecules': len(timing_data),
        'atom_range': [min(timing_atoms), max(timing_atoms)]
    },
    'benchmarks': {
        'enhanced': {
            'total': len(enhanced_data),
            'detected': sum(enhanced_detected),
            'detection_rate': detection_rate_enhanced,
            'unique_groups': len(enhanced_unique_groups),
            'top_groups': dict(top_groups)
        },
        'complex_branched': {
            'total': complex_total,
            'detected': complex_detected,
            'detection_rate': detection_rate_complex,
            'unique_groups': len(complex_unique_groups)
        },
        'highly_branched': {
            'total': hb_total,
            'detected': hb_detected,
            'detection_rate': detection_rate_hb,
            'unique_groups': len(hb_unique_groups)
        },
        'telomer_validation': {
            'total': telomer_total,
            'detected': telomer_detected_count,
            'detection_rate': detection_rate_telomer,
            'unique_groups': len(telomer_groups),
            'group_detections': telomer_group_counts
        }
    }
}

with open('reports/benchmark_summary.json', 'w') as f:
    json.dump(summary, f, indent=2)

print(f"\nSummary saved to: reports/benchmark_summary.json")

# Generate exponential fit plot
if PLOTLY_AVAILABLE:
    print("\n7. Generating exponential fit visualization...")
    
    # Create fit line data
    x_fit = list(range(min(timing_atoms), max(timing_atoms) + 1))
    y_fit = [a_fit * math.exp(b_fit * x) * 1000 for x in x_fit]
    
    # Create figure
    fig = go.Figure()
    
    # Add scatter plot for measured data
    fig.add_trace(go.Scatter(
        x=timing_atoms,
        y=[t*1000 for t in timing_avg],
        mode='markers',
        name='Measured times',
        marker=dict(
            size=10,
            color='#2196F3',
            opacity=0.7,
            line=dict(color='black', width=1)
        ),
        hovertemplate='<b>Atoms:</b> %{x}<br><b>Time:</b> %{y:.2f} ms<extra></extra>'
    ))
    
    # Add exponential fit line
    fig.add_trace(go.Scatter(
        x=x_fit,
        y=y_fit,
        mode='lines',
        name=f'Exponential fit: t = {a_fit:.2e} × exp({b_fit:.4f} × n) s',
        line=dict(color='red', width=3),
        hovertemplate='<b>Atoms:</b> %{x}<br><b>Fit:</b> %{y:.2f} ms<extra></extra>'
    ))
    
    # Update layout
    fig.update_layout(
        title=dict(
            text=f'PFASGroups Exponential Scaling Model<br><sub>R² = {r_squared:.4f}, Pearson r = {corr_coef:.4f}</sub>',
            font=dict(size=20, family='Arial, sans-serif', color='#333')
        ),
        xaxis=dict(
            title=dict(text='Number of Atoms (n)', font=dict(size=16, family='Arial, sans-serif')),
            showgrid=True,
            gridcolor='rgba(0,0,0,0.1)',
            zeroline=False
        ),
        yaxis=dict(
            title=dict(text='Execution Time (ms)', font=dict(size=16, family='Arial, sans-serif')),
            showgrid=True,
            gridcolor='rgba(0,0,0,0.1)',
            zeroline=False
        ),
        plot_bgcolor='white',
        paper_bgcolor='white',
        legend=dict(
            x=0.02,
            y=0.98,
            xanchor='left',
            yanchor='top',
            bgcolor='rgba(255,255,255,0.9)',
            bordercolor='#999',
            borderwidth=1,
            font=dict(size=12)
        ),
        width=1200,
        height=800,
        margin=dict(l=80, r=200, t=100, b=80),
        hovermode='closest'
    )
    
    # Add annotation box with parameters
    annotation_text = (
        f'<b>Parameters:</b><br>'
        f'a = {a_fit:.3e} s<br>'
        f'α = {b_fit:.4f} atoms⁻¹<br><br>'
        f'<b>Statistics:</b><br>'
        f'Mean: {mean_t*1000:.1f} ms<br>'
        f'Median: {median_timing*1000:.1f} ms<br>'
        f'Std: {std_timing*1000:.1f} ms<br><br>'
        f'<b>Dataset:</b><br>'
        f'n = {len(timing_data)} molecules<br>'
        f'Atom range: {min(timing_atoms)}-{max(timing_atoms)}'
    )
    
    fig.add_annotation(
        x=0.98,
        y=0.02,
        xref='paper',
        yref='paper',
        text=annotation_text,
        showarrow=False,
        xanchor='right',
        yanchor='bottom',
        bordercolor='#999',
        borderwidth=1,
        borderpad=10,
        bgcolor='rgba(255, 245, 220, 0.9)',
        font=dict(size=11, family='Courier New, monospace'),
        align='left'
    )
    
    # Save figure in multiple formats
    os.makedirs('imgs', exist_ok=True)
    
    # Save as PNG (high resolution for publications)
    png_path = 'imgs/timing_exponential_fit.png'
    fig.write_image(png_path, width=1200, height=800, scale=2)  # 2x scale for high DPI
    print(f"  [OK] PNG saved: {png_path}")
    
    # Save as SVG (vector format for LaTeX)
    svg_path = 'imgs/timing_exponential_fit.svg'
    fig.write_image(svg_path, width=1200, height=800)
    print(f"  [OK] SVG saved: {svg_path}")
    
    # Save as PDF (alternative vector format for LaTeX)
    pdf_path = 'imgs/timing_exponential_fit.pdf'
    fig.write_image(pdf_path, width=1200, height=800)
    print(f"  [OK] PDF saved: {pdf_path}")
    
    # Also save interactive HTML version
    html_path = 'reports/timing_exponential_fit.html'
    fig.write_html(html_path)
    print(f"  [OK] HTML (interactive) saved: {html_path}")
    
else:
    print("\n[WARNING] Skipping plot generation (plotly not available)")

print("="*80)
print("="*80)
