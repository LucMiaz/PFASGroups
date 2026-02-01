"""
Comprehensive benchmark analysis and statistics generation.
Analyzes timing models and generates summary statistics for LaTeX documents.
"""

import json
import numpy as np

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
timing_atoms = np.array([d['num_atoms'] for d in timing_data])
timing_avg = np.array([d['pfasgroups_time_avg'] for d in timing_data])

# Log-linear regression for exponential model: t = a * exp(b * n)
log_timing = np.log(timing_avg)
coeffs = np.polyfit(timing_atoms, log_timing, 1)
b_fit = coeffs[0]
a_fit = np.exp(coeffs[1])

# Calculate R²
y_pred = a_fit * np.exp(b_fit * timing_atoms)
ss_res = np.sum((timing_avg - y_pred)**2)
ss_tot = np.sum((timing_avg - np.mean(timing_avg))**2)
r_squared = 1 - (ss_res / ss_tot)

# Calculate correlation
corr_matrix = np.corrcoef(timing_atoms, timing_avg)
corr_coef = corr_matrix[0, 1]

print(f"\n  Exponential Model: t = a × exp(b × n)")
print(f"  Parameters:")
print(f"    a = {a_fit:.6e} seconds")
print(f"    b = {b_fit:.6f} atoms⁻¹ (α in paper)")
print(f"  Quality:")
print(f"    R² = {r_squared:.6f}")
print(f"    Pearson r = {corr_coef:.6f}")
print(f"  Statistics:")
print(f"    Mean time: {np.mean(timing_avg)*1000:.2f} ms")
print(f"    Median time: {np.median(timing_avg)*1000:.2f} ms")
print(f"    Std time: {np.std(timing_avg)*1000:.2f} ms")

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
complex_detected = []
complex_unique_groups = set()
for mol in complex_data:
    detected = mol['pfasgroups_result']['detected_groups']
    complex_detected.append(len(detected) > 0)
    complex_unique_groups.update(detected)

detection_rate_complex = 100.0 * sum(complex_detected) / len(complex_detected)
print(f"  Total molecules: {len(complex_data)}")
print(f"  Detected: {sum(complex_detected)}")
print(f"  Detection rate: {detection_rate_complex:.1f}%")
print(f"  Unique groups found: {len(complex_unique_groups)}")

# Analyze highly branched
print("\n5. Analyzing highly branched benchmark...")
hb_detected = []
hb_unique_groups = set()
for mol in highly_branched_data:
    detected = mol['pfasgroups_result']['detected_groups']
    hb_detected.append(len(detected) > 0)
    hb_unique_groups.update(detected)

detection_rate_hb = 100.0 * sum(hb_detected) / len(hb_detected)
print(f"  Total molecules: {len(highly_branched_data)}")
print(f"  Detected: {sum(hb_detected)}")
print(f"  Detection rate: {detection_rate_hb:.1f}%")
print(f"  Unique groups found: {len(hb_unique_groups)}")

# Analyze telomer validation
print("\n6. Analyzing telomer validation...")
telomer_detected = []
telomer_groups = set()
for mol in telomer_data:
    detected = mol.get('detected_groups', [])
    telomer_detected.append(len(detected) > 0)
    telomer_groups.update(detected)

detection_rate_telomer = 100.0 * sum(telomer_detected) / len(telomer_detected)
print(f"  Total fluorotelomers: {len(telomer_data)}")
print(f"  Detected: {sum(telomer_detected)}")
print(f"  Detection rate: {detection_rate_telomer:.1f}%")
print(f"  Unique groups found: {len(telomer_groups)}")

# Count telomer group frequencies
telomer_group_counts = {}
for mol in telomer_data:
    for gid in mol.get('detected_groups', []):
        telomer_group_counts[gid] = telomer_group_counts.get(gid, 0) + 1

top_telomer_groups = sorted(telomer_group_counts.items(), key=lambda x: x[1], reverse=True)
print(f"  Telomer group detections:")
for gid, count in top_telomer_groups:
    print(f"    Group {gid}: {count} ({100.0*count/len(telomer_data):.1f}%)")

# Generate LaTeX table rows
print("\n" + "="*80)
print("LATEX TABLE ROWS")
print("="*80)

print("\nTiming Model Parameters:")
print(f"All metrics & ${a_fit:.2e}$ & ${b_fit:.4f}$ & ${r_squared:.4f}$ & {np.mean(timing_avg)*1000:.1f} & {np.median(timing_avg)*1000:.1f} \\\\")

print("\nBenchmark Detection Statistics:")
print(f"Enhanced & {len(enhanced_data)} & {sum(enhanced_detected)} & {detection_rate_enhanced:.1f}\\% & {len(enhanced_unique_groups)} \\\\")
print(f"Complex Branched & {len(complex_data)} & {sum(complex_detected)} & {detection_rate_complex:.1f}\\% & {len(complex_unique_groups)} \\\\")
print(f"Highly Branched & {len(highly_branched_data)} & {sum(hb_detected)} & {detection_rate_hb:.1f}\\% & {len(hb_unique_groups)} \\\\")
print(f"Telomer Validation & {len(telomer_data)} & {sum(telomer_detected)} & {detection_rate_telomer:.1f}\\% & {len(telomer_groups)} \\\\")

# Save summary to JSON
summary = {
    'timing_model': {
        'model': 't = a × exp(b × n)',
        'a': float(a_fit),
        'b': float(b_fit),
        'r_squared': float(r_squared),
        'correlation': float(corr_coef),
        'mean_time_ms': float(np.mean(timing_avg) * 1000),
        'median_time_ms': float(np.median(timing_avg) * 1000),
        'std_time_ms': float(np.std(timing_avg) * 1000),
        'n_molecules': len(timing_data),
        'atom_range': [int(min(timing_atoms)), int(max(timing_atoms))]
    },
    'benchmarks': {
        'enhanced': {
            'total': len(enhanced_data),
            'detected': sum(enhanced_detected),
            'detection_rate': detection_rate_enhanced,
            'unique_groups': len(enhanced_unique_groups)
        },
        'complex_branched': {
            'total': len(complex_data),
            'detected': sum(complex_detected),
            'detection_rate': detection_rate_complex,
            'unique_groups': len(complex_unique_groups)
        },
        'highly_branched': {
            'total': len(highly_branched_data),
            'detected': sum(hb_detected),
            'detection_rate': detection_rate_hb,
            'unique_groups': len(hb_unique_groups)
        },
        'telomer_validation': {
            'total': len(telomer_data),
            'detected': sum(telomer_detected),
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
