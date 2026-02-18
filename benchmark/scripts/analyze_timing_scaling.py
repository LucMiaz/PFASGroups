import sqlite3
import numpy as np
from scipy.optimize import curve_fit

conn = sqlite3.connect('review-app/database/pfas_benchmark.db')
c = conn.cursor()

# Get timing dataset with molecule properties
c.execute('''
    SELECT 
        m.num_atoms,
        m.molecular_weight,
        p.execution_time as pfas_time,
        a.execution_time as atlas_time
    FROM molecules m
    JOIN HalogenGroups_results p ON m.id = p.molecule_id
    JOIN atlas_results a ON m.id = a.molecule_id
    WHERE m.dataset_type = 'timing'
      AND p.execution_time IS NOT NULL 
      AND a.execution_time IS NOT NULL
    ORDER BY m.num_atoms
''')

rows = c.fetchall()
conn.close()

if not rows:
    print("No timing data found!")
    exit(1)

# Extract data
atom_counts = np.array([r[0] for r in rows])
mol_weights = np.array([r[1] for r in rows])
pfas_times = np.array([r[2] * 1000 for r in rows])  # Convert to ms
atlas_times = np.array([r[3] * 1000 for r in rows])

print(f"Timing dataset: {len(rows)} molecules")
print(f"Atom count range: {atom_counts.min()}-{atom_counts.max()}")
print(f"HalogenGroups time range: {pfas_times.min():.1f}-{pfas_times.max():.1f} ms")
print(f"Atlas time range: {atlas_times.min():.1f}-{atlas_times.max():.1f} ms")

# Define models
def exponential_model(x, a, b):
    """time = a * exp(b * x)"""
    return a * np.exp(b * x)

def power_law_model(x, a, b):
    """time = a * x^b"""
    return a * np.power(x, b)

def linear_model(x, a, b):
    """time = a * x + b"""
    return a * x + b

# Fit models for HalogenGroups
print("\n" + "="*60)
print("HalogenGroupS SCALING ANALYSIS")
print("="*60)

try:
    # Exponential fit
    popt_exp, pcov_exp = curve_fit(exponential_model, atom_counts, pfas_times, 
                                     p0=[10, 0.01], maxfev=10000)
    pred_exp = exponential_model(atom_counts, *popt_exp)
    r2_exp = 1 - (np.sum((pfas_times - pred_exp)**2) / np.sum((pfas_times - np.mean(pfas_times))**2))
    rmse_exp = np.sqrt(np.mean((pfas_times - pred_exp)**2))
    
    print(f"\nExponential Model: time = {popt_exp[0]:.3f} * exp({popt_exp[1]:.5f} * atoms)")
    print(f"  R² = {r2_exp:.4f}")
    print(f"  RMSE = {rmse_exp:.1f} ms")
except Exception as e:
    print(f"\nExponential fit failed: {e}")
    popt_exp = None
    r2_exp = -np.inf

try:
    # Power law fit
    popt_pow, pcov_pow = curve_fit(power_law_model, atom_counts, pfas_times, 
                                     p0=[1, 1.5], maxfev=10000)
    pred_pow = power_law_model(atom_counts, *popt_pow)
    r2_pow = 1 - (np.sum((pfas_times - pred_pow)**2) / np.sum((pfas_times - np.mean(pfas_times))**2))
    rmse_pow = np.sqrt(np.mean((pfas_times - pred_pow)**2))
    
    print(f"\nPower Law Model: time = {popt_pow[0]:.3f} * atoms^{popt_pow[1]:.3f}")
    print(f"  R² = {r2_pow:.4f}")
    print(f"  RMSE = {rmse_pow:.1f} ms")
except Exception as e:
    print(f"\nPower law fit failed: {e}")
    popt_pow = None
    r2_pow = -np.inf

try:
    # Linear fit
    popt_lin, pcov_lin = curve_fit(linear_model, atom_counts, pfas_times)
    pred_lin = linear_model(atom_counts, *popt_lin)
    r2_lin = 1 - (np.sum((pfas_times - pred_lin)**2) / np.sum((pfas_times - np.mean(pfas_times))**2))
    rmse_lin = np.sqrt(np.mean((pfas_times - pred_lin)**2))
    
    print(f"\nLinear Model: time = {popt_lin[0]:.3f} * atoms + {popt_lin[1]:.1f}")
    print(f"  R² = {r2_lin:.4f}")
    print(f"  RMSE = {rmse_lin:.1f} ms")
except Exception as e:
    print(f"\nLinear fit failed: {e}")
    r2_lin = -np.inf

# Determine best model
best_r2 = max(r2_exp if popt_exp is not None else -np.inf, 
              r2_pow if popt_pow is not None else -np.inf, 
              r2_lin)

if best_r2 == r2_exp:
    best_model = "exponential"
    best_params = popt_exp
elif best_r2 == r2_pow:
    best_model = "power law"
    best_params = popt_pow
else:
    best_model = "linear"
    best_params = popt_lin

print(f"\n>>> Best model: {best_model} (R² = {best_r2:.4f})")

# Atlas analysis for comparison
print("\n" + "="*60)
print("ATLAS SCALING ANALYSIS")
print("="*60)

try:
    popt_atlas_lin, _ = curve_fit(linear_model, atom_counts, atlas_times)
    pred_atlas_lin = linear_model(atom_counts, *popt_atlas_lin)
    r2_atlas_lin = 1 - (np.sum((atlas_times - pred_atlas_lin)**2) / 
                        np.sum((atlas_times - np.mean(atlas_times))**2))
    
    print(f"\nLinear Model: time = {popt_atlas_lin[0]:.3f} * atoms + {popt_atlas_lin[1]:.1f}")
    print(f"  R² = {r2_atlas_lin:.4f}")
except Exception as e:
    print(f"Linear fit failed: {e}")

# Statistics by size bins
print("\n" + "="*60)
print("STATISTICS BY SIZE BINS")
print("="*60)

bins = [(0, 20), (20, 40), (40, 60), (60, 100)]
for low, high in bins:
    mask = (atom_counts >= low) & (atom_counts < high)
    if np.sum(mask) > 0:
        print(f"\n{low}-{high} atoms ({np.sum(mask)} molecules):")
        print(f"  HalogenGroups: {pfas_times[mask].mean():.1f} ± {pfas_times[mask].std():.1f} ms")
        print(f"  Atlas:      {atlas_times[mask].mean():.1f} ± {atlas_times[mask].std():.1f} ms")
        print(f"  Ratio:      {pfas_times[mask].mean() / atlas_times[mask].mean():.2f}×")
