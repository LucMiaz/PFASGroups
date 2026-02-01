"""
Comprehensive benchmark analysis without numpy dependency.
Analyzes timing models and generates summary statistics for LaTeX documents.
"""

import json
import math

print("="*80)
print("PFASGROUPS BENCHMARK ANALYSIS")
print("="*80)

# Load all benchmark data
print("\n1. Loading benchmark data...")

# Timing benchmark
with open('data/pfas_timing_benchmark_20260201_022433.json', 'r') as f:
    timing_data = json.load(f)

# Enhanced benchmark
with open('data/pfas_enhanced_benchmark_20260201_022011.json', 'r') as f:
    enhanced_data = json.load(f)

# Complex branched
with open('data/pfas_complex_branched_benchmark_20260201_022453.json', 'r') as f:
    complex_data = json.load(f)

# Highly branched
with open('data/pfas_highly_branched_benchmark_20260201_022459.json', 'r') as f:
    highly_branched_data = json.load(f)

# Telomer validation
with open('data/telomer_validation_results.json', 'r') as f:
    telomer_data = json.load(f)

print(f"  Timing: {len(timing_data)} molecules")
print(f"  Enhanced: {len(enhanced_data)} molecules")
print(f"  Complex branched: {len(complex_data)} molecules")
print(f"  Highly branched: {len(highly_branched_data)} molecules")
print(f"  Telomer validation: {len(telomer_data)} molecules")

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

print(f"\n  Exponential Model: t = a × exp(b × n)")
print(f"  Parameters:")
print(f"    a = {a_fit:.6e} seconds")
print(f"    b = {b_fit:.6f} atoms⁻¹ (α in paper)")
print(f"  Quality:")
print(f"    R² = {r_squared:.6f}")
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
hb_detected = highly_branched_data.get('summary', {}).get('passed', 0)
hb_unique_groups = set()
for detail in highly_branched_data.get('details', []):
    if detail.get('passed', False):
        hb_unique_groups.add(detail.get('group_id'))

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
print("="*80)
