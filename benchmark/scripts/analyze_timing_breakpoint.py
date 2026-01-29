"""
Analyze timing data for breakpoint at 100 atoms
Fit two separate models and investigate the cause
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy import stats
import json

# Configure matplotlib
plt.rcParams['font.family'] = 'serif'
plt.rcParams['pdf.fonttype'] = 42

# Load data
data = np.load('data/timing_full_dataset.npz')
atom_counts = data['atom_counts']
exec_times = data['exec_times']

print(f"Total data points: {len(atom_counts)}")
print(f"Atom range: {atom_counts.min()}-{atom_counts.max()}")
print(f"Time range: {exec_times.min():.1f}-{exec_times.max():.1f} ms")

# Split data at 100 atoms
mask_low = atom_counts < 100
mask_high = atom_counts >= 100

atoms_low = atom_counts[mask_low]
times_low = exec_times[mask_low]
atoms_high = atom_counts[mask_high]
times_high = exec_times[mask_high]

print(f"\nData split:")
print(f"  < 100 atoms: {len(atoms_low)} points (range: {atoms_low.min()}-{atoms_low.max()} atoms)")
print(f"  ≥ 100 atoms: {len(atoms_high)} points (range: {atoms_high.min()}-{atoms_high.max()} atoms)")

# Analyze timing around the breakpoint
breakpoint_range = (atom_counts >= 90) & (atom_counts <= 110)
atoms_breakpoint = atom_counts[breakpoint_range]
times_breakpoint = exec_times[breakpoint_range]

print(f"\nBreakpoint region (90-110 atoms):")
for n in sorted(np.unique(atoms_breakpoint)):
    mask = atoms_breakpoint == n
    mean_time = times_breakpoint[mask].mean()
    std_time = times_breakpoint[mask].std()
    count = mask.sum()
    print(f"  {n:3d} atoms: {mean_time:6.1f} ± {std_time:5.1f} ms (n={count})")

# Calculate speed jump at 100 atoms
times_99 = exec_times[(atom_counts >= 97) & (atom_counts < 100)].mean()
times_100 = exec_times[atom_counts >= 100].mean()
jump_factor = times_100 / times_99

print(f"\nSpeed jump analysis:")
print(f"  Mean < 100 atoms: {times_99:.1f} ms")
print(f"  Mean ≥ 100 atoms: {times_100:.1f} ms")
print(f"  Jump factor: {jump_factor:.2f}×")

# Define exponential model
def exponential_model(x, a, b):
    return a * np.exp(b * x)

# Fit model for < 100 atoms
popt_low, pcov_low = curve_fit(exponential_model, atoms_low, times_low, 
                                 p0=[50, 0.015], maxfev=10000)
y_pred_low_data = exponential_model(atoms_low, *popt_low)
r2_low = 1 - np.sum((times_low - y_pred_low_data)**2) / np.sum((times_low - times_low.mean())**2)
rmse_low = np.sqrt(np.mean((times_low - y_pred_low_data)**2))

print(f"\nModel for < 100 atoms:")
print(f"  t = {popt_low[0]:.2f} × exp({popt_low[1]:.5f} × n)")
print(f"  R² = {r2_low:.4f}, RMSE = {rmse_low:.1f} ms")
print(f"  α = {popt_low[1]:.5f} atoms⁻¹")

# Fit model for ≥ 100 atoms
popt_high, pcov_high = curve_fit(exponential_model, atoms_high, times_high,
                                   p0=[200, 0.01], maxfev=10000)
y_pred_high_data = exponential_model(atoms_high, *popt_high)
r2_high = 1 - np.sum((times_high - y_pred_high_data)**2) / np.sum((times_high - times_high.mean())**2)
rmse_high = np.sqrt(np.mean((times_high - y_pred_high_data)**2))

print(f"\nModel for ≥ 100 atoms:")
print(f"  t = {popt_high[0]:.2f} × exp({popt_high[1]:.5f} × n)")
print(f"  R² = {r2_high:.4f}, RMSE = {rmse_high:.1f} ms")
print(f"  α = {popt_high[1]:.5f} atoms⁻¹")

# Compare growth rates
print(f"\nGrowth rate comparison:")
print(f"  α_low / α_high = {popt_low[1] / popt_high[1]:.2f}×")
print(f"  Growth is {popt_low[1] / popt_high[1]:.1f}× faster for < 100 atoms")

# Create visualization
fig, axes = plt.subplots(2, 2, figsize=(14, 10))

# Plot 1: Two-region fit
ax1 = axes[0, 0]
x_low = np.linspace(atoms_low.min(), atoms_low.max(), 100)
x_high = np.linspace(atoms_high.min(), atoms_high.max(), 100)
y_low = exponential_model(x_low, *popt_low)
y_high = exponential_model(x_high, *popt_high)

ax1.scatter(atoms_low, times_low, alpha=0.6, s=30, color='#667eea', label='Data < 100 atoms')
ax1.scatter(atoms_high, times_high, alpha=0.6, s=30, color='#f56565', label='Data ≥ 100 atoms')
ax1.plot(x_low, y_low, 'b-', linewidth=2, label=f'Fit < 100: α={popt_low[1]:.5f}')
ax1.plot(x_high, y_high, 'r-', linewidth=2, label=f'Fit ≥ 100: α={popt_high[1]:.5f}')
ax1.axvline(100, color='green', linestyle='--', linewidth=2, alpha=0.7, label='Breakpoint at 100')
ax1.set_xlabel('Number of Atoms', fontsize=11, fontweight='bold')
ax1.set_ylabel('Execution Time (ms)', fontsize=11, fontweight='bold')
ax1.set_title('Two-Region Exponential Fit', fontsize=12, fontweight='bold')
ax1.legend(fontsize=9)
ax1.grid(True, alpha=0.3)

# Plot 2: Residuals
ax2 = axes[0, 1]
residuals_low = times_low - y_pred_low_data
residuals_high = times_high - y_pred_high_data
ax2.scatter(atoms_low, residuals_low, alpha=0.6, s=30, color='#667eea', label='< 100 atoms')
ax2.scatter(atoms_high, residuals_high, alpha=0.6, s=30, color='#f56565', label='≥ 100 atoms')
ax2.axhline(0, color='black', linestyle='-', linewidth=1)
ax2.axvline(100, color='green', linestyle='--', linewidth=2, alpha=0.7)
ax2.set_xlabel('Number of Atoms', fontsize=11, fontweight='bold')
ax2.set_ylabel('Residuals (ms)', fontsize=11, fontweight='bold')
ax2.set_title('Model Residuals', fontsize=12, fontweight='bold')
ax2.legend(fontsize=9)
ax2.grid(True, alpha=0.3)

# Plot 3: Breakpoint region detail
ax3 = axes[1, 0]
breakpoint_atoms = np.unique(atoms_breakpoint)
breakpoint_means = [times_breakpoint[atoms_breakpoint == n].mean() for n in breakpoint_atoms]
breakpoint_stds = [times_breakpoint[atoms_breakpoint == n].std() for n in breakpoint_atoms]
ax3.errorbar(breakpoint_atoms, breakpoint_means, yerr=breakpoint_stds, 
             fmt='o-', capsize=5, markersize=8, linewidth=2, color='purple')
ax3.axvline(100, color='green', linestyle='--', linewidth=2, alpha=0.7, label='Breakpoint')
ax3.set_xlabel('Number of Atoms', fontsize=11, fontweight='bold')
ax3.set_ylabel('Mean Execution Time (ms)', fontsize=11, fontweight='bold')
ax3.set_title('Detail: 90-110 Atoms', fontsize=12, fontweight='bold')
ax3.legend(fontsize=9)
ax3.grid(True, alpha=0.3)

# Plot 4: Time per atom
ax4 = axes[1, 1]
bins = np.arange(10, 160, 10)
time_per_atom = []
bin_centers = []
for i in range(len(bins)-1):
    mask = (atom_counts >= bins[i]) & (atom_counts < bins[i+1])
    if mask.sum() > 0:
        time_per_atom.append((exec_times[mask] / atom_counts[mask]).mean())
        bin_centers.append((bins[i] + bins[i+1]) / 2)

ax4.plot(bin_centers, time_per_atom, 'o-', markersize=8, linewidth=2, color='orange')
ax4.axvline(100, color='green', linestyle='--', linewidth=2, alpha=0.7, label='Breakpoint')
ax4.set_xlabel('Number of Atoms (bin center)', fontsize=11, fontweight='bold')
ax4.set_ylabel('Time per Atom (ms/atom)', fontsize=11, fontweight='bold')
ax4.set_title('Computational Cost per Atom', fontsize=12, fontweight='bold')
ax4.legend(fontsize=9)
ax4.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig('imgs/timing_breakpoint_analysis.pdf', dpi=300, bbox_inches='tight')
plt.savefig('imgs/timing_breakpoint_analysis.png', dpi=300, bbox_inches='tight')
print(f"\n✓ Saved analysis plots to imgs/timing_breakpoint_analysis.pdf/png")

plt.close()

# Save analysis results
results = {
    'total_points': int(len(atom_counts)),
    'breakpoint': 100,
    'low_region': {
        'n_points': int(len(atoms_low)),
        'atom_range': [int(atoms_low.min()), int(atoms_low.max())],
        'time_range': [float(times_low.min()), float(times_low.max())],
        'model': f't = {popt_low[0]:.2f} × exp({popt_low[1]:.5f} × n)',
        'a': float(popt_low[0]),
        'alpha': float(popt_low[1]),
        'R2': float(r2_low),
        'RMSE': float(rmse_low)
    },
    'high_region': {
        'n_points': int(len(atoms_high)),
        'atom_range': [int(atoms_high.min()), int(atoms_high.max())],
        'time_range': [float(times_high.min()), float(times_high.max())],
        'model': f't = {popt_high[0]:.2f} × exp({popt_high[1]:.5f} × n)',
        'a': float(popt_high[0]),
        'alpha': float(popt_high[1]),
        'R2': float(r2_high),
        'RMSE': float(rmse_high)
    },
    'speed_jump': {
        'mean_below_100': float(times_99),
        'mean_above_100': float(times_100),
        'jump_factor': float(jump_factor)
    },
    'growth_rate_ratio': float(popt_low[1] / popt_high[1])
}

with open('data/timing_breakpoint_analysis.json', 'w') as f:
    json.dump(results, f, indent=2)

print(f"\n✓ Saved analysis results to data/timing_breakpoint_analysis.json")
