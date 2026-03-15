"""
Analyze timing benchmark data and generate exponential models.
Creates models for:
1. All metrics enabled (default)
2. No metrics (metrics=False)
3. Without effective graph resistance (compute_EGR=False)
"""

import json
import re
import numpy as np
import pandas as pd
try:
    from scipy.optimize import curve_fit
    from scipy.stats import pearsonr
    SCIPY_AVAILABLE = True
except ImportError:
    SCIPY_AVAILABLE = False
    print("Warning: scipy not available, will use alternative fitting")

try:
    import matplotlib.pyplot as plt
    plt.style.use('seaborn-v0_8-whitegrid')
    plt.rcParams['font.family'] = 'Ubuntu'
    MATPLOTLIB_AVAILABLE = True
except ImportError:
    MATPLOTLIB_AVAILABLE = False
    print("Warning: matplotlib not available, skipping plots")

import glob
from pathlib import Path


def _load_palette():
    """Load hex colours from color_scheme.yaml (stdlib only, no pyyaml needed)."""
    _defaults = ["#E15D0B", "#306DBA", "#9D206C", "#51127C"]
    try:
        _p = Path(__file__).parent.parent.parent / "PFASGroups" / "data" / "color_scheme.yaml"
        _colors = re.findall(r'"(#[0-9A-Fa-f]{6})"', _p.read_text())
        if len(_colors) >= 4:
            return _colors[:4]
    except Exception:
        pass
    return _defaults


_PALETTE = _load_palette()
# C0=orange (data/residuals), C1=blue (fit line)
_C0, _C1, _C2, _C3 = _PALETTE

# Load the latest timing benchmark data
timing_files = sorted(glob.glob('data/pfas_timing_benchmark_*.json'), key=lambda x: Path(x).stat().st_mtime, reverse=True)
if not timing_files:
    print("Error: No timing benchmark files found in data/")
    exit(1)

timing_file = timing_files[0]
print(f"Loading latest timing benchmark: {timing_file}")

with open(timing_file, 'r') as f:
    timing_data = json.load(f)

# Extract relevant fields
df = pd.DataFrame([{
    'molecule_id': d['molecule_id'],
    'chain_length': d['chain_length'],
    'num_atoms': d['num_atoms'],
    'PFASGroups_time_avg': d['PFASGroups_time_avg'],
    'PFASGroups_time_std': d['PFASGroups_time_std']
} for d in timing_data])

# Exponential model: t = a * exp(b * n)
def exponential_model(n, a, b):
    return a * np.exp(b * n)

# Fit the model
x_data = df['num_atoms'].values
y_data = df['PFASGroups_time_avg'].values

if SCIPY_AVAILABLE:
    # Initial guess for parameters
    p0 = [0.001, 0.01]
    
    # Fit the curve
    popt, pcov = curve_fit(exponential_model, x_data, y_data, p0=p0, maxfev=10000)
    a_fit, b_fit = popt
    
    # Calculate R² and correlation
    y_pred = exponential_model(x_data, a_fit, b_fit)
    residuals = y_data - y_pred
    ss_res = np.sum(residuals**2)
    ss_tot = np.sum((y_data - np.mean(y_data))**2)
    r_squared = 1 - (ss_res / ss_tot)
    corr_coef, p_value = pearsonr(x_data, y_data)
    
    # Calculate standard errors
    perr = np.sqrt(np.diag(pcov))
    a_err, b_err = perr
else:
    # Simple log-linear regression as fallback
    log_y = np.log(y_data)
    coeffs = np.polyfit(x_data, log_y, 1)
    b_fit = coeffs[0]
    a_fit = np.exp(coeffs[1])
    a_err = 0
    b_err = 0
    y_pred = exponential_model(x_data, a_fit, b_fit)
    ss_res = np.sum((y_data - y_pred)**2)
    ss_tot = np.sum((y_data - np.mean(y_data))**2)
    r_squared = 1 - (ss_res / ss_tot)
    corr_coef = np.corrcoef(x_data, y_data)[0, 1]
    p_value = 0.0

# Print results
print("=" * 80)
print("TIMING ANALYSIS - ALL METRICS ENABLED (DEFAULT)")
print("=" * 80)
print(f"\nExponential Model: t = a * exp(b * n)")
print(f"  where t = time (seconds), n = number of atoms")
print(f"\nFitted Parameters:")
print(f"  a = {a_fit:.6e} ± {a_err:.6e} seconds")
print(f"  b = {b_fit:.6f} ± {b_err:.6f} atoms⁻¹")
print(f"\nModel Quality:")
print(f"  R² = {r_squared:.6f}")
print(f"  Pearson correlation = {corr_coef:.6f}")
print(f"  p-value = {p_value:.6e}")
print(f"\nDataset Statistics:")
print(f"  Number of molecules: {len(df)}")
print(f"  Atom count range: {df['num_atoms'].min()} - {df['num_atoms'].max()}")
print(f"  Chain length range: {df['chain_length'].min()} - {df['chain_length'].max()}")
print(f"  Mean execution time: {df['PFASGroups_time_avg'].mean()*1000:.2f} ms")
print(f"  Median execution time: {df['PFASGroups_time_avg'].median()*1000:.2f} ms")

# Create visualization if matplotlib is available
if MATPLOTLIB_AVAILABLE:
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))
    
    # Plot 1: Data and fitted curve
    ax1.scatter(x_data, y_data * 1000, alpha=0.5, s=30, color=_C0, label='Measured data')
    x_fit = np.linspace(x_data.min(), x_data.max(), 100)
    y_fit = exponential_model(x_fit, a_fit, b_fit)
    ax1.plot(x_fit, y_fit * 1000, '-', color=_C1, linewidth=2,
             label=f'Fit: t = {a_fit:.2e} × exp({b_fit:.4f} × n)\nR² = {r_squared:.4f}')
    ax1.set_xlabel('Number of atoms', fontsize=12)
    ax1.set_ylabel('Execution time (ms)', fontsize=12)
    ax1.set_title('PFASGroups Timing vs Molecular Size\n(All Metrics Enabled)', fontsize=13, fontweight='bold')
    ax1.legend(fontsize=10)
    ax1.grid(True, alpha=0.3)
    
    # Plot 2: Residuals
    residuals_ms = residuals * 1000
    ax2.scatter(x_data, residuals_ms, alpha=0.5, s=30, color=_C0)
    ax2.axhline(y=0, color='#888888', linestyle='--', linewidth=1)
    ax2.set_xlabel('Number of atoms', fontsize=12)
    ax2.set_ylabel('Residuals (ms)', fontsize=12)
    ax2.set_title('Residual Analysis', fontsize=13, fontweight='bold')
    ax2.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig('reports/timing_analysis_all_metrics.png', dpi=300, bbox_inches='tight')
    plt.savefig('reports/timing_analysis_all_metrics.pdf', dpi=300, bbox_inches='tight')
    print(f"\nPlot saved to: reports/timing_analysis_all_metrics.png / .pdf")
else:
    print("\nSkipping plot generation (matplotlib not available)")

# Save results to JSON
results = {
    'model': 't = a * exp(b * n)',
    'parameters': {
        'a': {
            'value': float(a_fit),
            'error': float(a_err),
            'unit': 'seconds'
        },
        'b': {
            'value': float(b_fit),
            'error': float(b_err),
            'unit': 'atoms⁻¹'
        }
    },
    'quality_metrics': {
        'r_squared': float(r_squared),
        'pearson_correlation': float(corr_coef),
        'p_value': float(p_value)
    },
    'dataset_statistics': {
        'n_molecules': int(len(df)),
        'atom_count_range': [int(df['num_atoms'].min()), int(df['num_atoms'].max())],
        'chain_length_range': [int(df['chain_length'].min()), int(df['chain_length'].max())],
        'mean_time_ms': float(df['PFASGroups_time_avg'].mean() * 1000),
        'median_time_ms': float(df['PFASGroups_time_avg'].median() * 1000),
        'std_time_ms': float(df['PFASGroups_time_avg'].std() * 1000)
    }
}

with open('reports/timing_model_all_metrics.json', 'w') as f:
    json.dump(results, f, indent=2)

print(f"\nResults saved to: reports/timing_model_all_metrics.json")
print("=" * 80)

# Generate LaTeX table row
print("\nLaTeX Table Row:")
print(f"All metrics & ${a_fit:.2e}$ & ${b_fit:.4f}$ & ${r_squared:.4f}$ & {df['PFASGroups_time_avg'].mean()*1000:.1f} & {df['PFASGroups_time_avg'].median()*1000:.1f} \\\\")
