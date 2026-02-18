import json
import numpy as np
from scipy.optimize import curve_fit

# Data extracted from HTML timing scatter plot
# Format: atom_count -> list of execution times (ms)
data_by_atoms = {}

# Sample of extracted data (grouping by atom count)
atom_times = [
    (13, [13.82, 13.87, 15.03, 13.67, 12.97, 14.28, 13.65, 14.94]),  # Chain 3
    (16, [17.78, 17.93, 18.55, 17.64, 17.83, 18.20, 17.62, 17.11]),  # Chain 4
    (19, [22.26, 22.05, 21.59, 22.77, 22.28, 23.15, 22.88, 21.34]),  # Chain 5
    (22, [30.93, 28.39, 27.05, 27.46, 27.75, 29.23, 27.51, 27.89]),  # Chain 6
    (34, [54.49, 55.33, 56.19, 53.20, 52.15, 53.67, 53.49, 55.28, 58.22, 58.09, 56.69]),  # Chain 10
    (70, [219.68, 227.60, 249.91, 218.59]),  # Chain 22
    (73, [238.95, 230.51, 231.27, 234.15]),  # Chain 23
    (76, [352.06, 267.30, 257.18, 235.88]),  # Chain 24
    (79, [313.77, 284.71, 270.64, 271.29, 278.34, 280.55, 266.94]),  # Chain 25
    (127, [374.41, 380.80, 383.19, 379.81]),  # Chain 41
    (154, [663.66, 576.62, 626.78]),  # Chain 50
]

# Create arrays for regression
atom_counts = []
exec_times = []

for atoms, times in atom_times:
    for time in times:
        atom_counts.append(atoms)
        exec_times.append(time)

atom_counts = np.array(atom_counts)
exec_times = np.array(exec_times)

print(f"Data points: {len(atom_counts)}")
print(f"Atom range: {atom_counts.min()}-{atom_counts.max()}")
print(f"Time range: {exec_times.min():.1f}-{exec_times.max():.1f} ms")

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

def quadratic_model(x, a, b, c):
    """time = a * x^2 + b * x + c"""
    return a * x**2 + b * x + c

# Fit models
print("\n" + "="*70)
print("TIMING DATASET SCALING ANALYSIS (HalogenGroups)")
print("="*70)

models = []

try:
    popt_exp, _ = curve_fit(exponential_model, atom_counts, exec_times, 
                            p0=[1, 0.01], maxfev=10000)
    pred_exp = exponential_model(atom_counts, *popt_exp)
    r2_exp = 1 - (np.sum((exec_times - pred_exp)**2) / 
                 np.sum((exec_times - np.mean(exec_times))**2))
    rmse_exp = np.sqrt(np.mean((exec_times - pred_exp)**2))
    
    print(f"\nExponential Model: t = {popt_exp[0]:.3f} × exp({popt_exp[1]:.5f} × n)")
    print(f"  R² = {r2_exp:.4f}, RMSE = {rmse_exp:.1f} ms")
    models.append(('exponential', r2_exp, popt_exp, (popt_exp[0], popt_exp[1])))
except Exception as e:
    print(f"\nExponential fit failed: {e}")

try:
    popt_pow, _ = curve_fit(power_law_model, atom_counts, exec_times, 
                            p0=[0.1, 1.5], maxfev=10000)
    pred_pow = power_law_model(atom_counts, *popt_pow)
    r2_pow = 1 - (np.sum((exec_times - pred_pow)**2) / 
                 np.sum((exec_times - np.mean(exec_times))**2))
    rmse_pow = np.sqrt(np.mean((exec_times - pred_pow)**2))
    
    print(f"\nPower Law Model: t = {popt_pow[0]:.3f} × n^{popt_pow[1]:.3f}")
    print(f"  R² = {r2_pow:.4f}, RMSE = {rmse_pow:.1f} ms")
    print(f"  (where n = atom count, t = time in ms)")
    models.append(('power_law', r2_pow, popt_pow, (popt_pow[0], popt_pow[1])))
except Exception as e:
    print(f"\nPower law fit failed: {e}")

try:
    popt_lin, _ = curve_fit(linear_model, atom_counts, exec_times)
    pred_lin = linear_model(atom_counts, *popt_lin)
    r2_lin = 1 - (np.sum((exec_times - pred_lin)**2) / 
                 np.sum((exec_times - np.mean(exec_times))**2))
    rmse_lin = np.sqrt(np.mean((exec_times - pred_lin)**2))
    
    print(f"\nLinear Model: t = {popt_lin[0]:.3f} × n + {popt_lin[1]:.1f}")
    print(f"  R² = {r2_lin:.4f}, RMSE = {rmse_lin:.1f} ms")
    models.append(('linear', r2_lin, popt_lin, (popt_lin[0], popt_lin[1])))
except Exception as e:
    print(f"\nLinear fit failed: {e}")

try:
    popt_quad, _ = curve_fit(quadratic_model, atom_counts, exec_times, p0=[0.01, 1, 0])
    pred_quad = quadratic_model(atom_counts, *popt_quad)
    r2_quad = 1 - (np.sum((exec_times - pred_quad)**2) / 
                  np.sum((exec_times - np.mean(exec_times))**2))
    rmse_quad = np.sqrt(np.mean((exec_times - pred_quad)**2))
    
    print(f"\nQuadratic Model: t = {popt_quad[0]:.5f} × n² + {popt_quad[1]:.3f} × n + {popt_quad[2]:.1f}")
    print(f"  R² = {r2_quad:.4f}, RMSE = {rmse_quad:.1f} ms")
    models.append(('quadratic', r2_quad, popt_quad, (popt_quad[0], popt_quad[1], popt_quad[2])))
except Exception as e:
    print(f"\nQuadratic fit failed: {e}")

# Find best model
if models:
    best = max(models, key=lambda x: x[1])
    print(f"\n{'='*70}")
    print(f"BEST MODEL: {best[0].upper()} (R² = {best[1]:.4f})")
    print(f"{'='*70}")
    
    if best[0] == 'exponential':
        a, b = best[3]
        print(f"\nExecution time = {a:.3f} × exp({b:.5f} × atom_count)")
        print(f"\nInterpretation: Time grows exponentially with molecular size")
        print(f"  For small molecules (~20 atoms): ~{exponential_model(20, a, b):.0f} ms")
        print(f"  For medium molecules (~50 atoms): ~{exponential_model(50, a, b):.0f} ms")
        print(f"  For large molecules (~100 atoms): ~{exponential_model(100, a, b):.0f} ms")
    elif best[0] == 'power_law':
        a, b = best[3]
        print(f"\nExecution time = {a:.3f} × atom_count^{b:.3f}")
        print(f"\nInterpretation: Time scales as a power law (approximately {'quadratic' if b > 1.5 else 'linear'})")
        print(f"  For small molecules (~20 atoms): ~{power_law_model(20, a, b):.0f} ms")
        print(f"  For medium molecules (~50 atoms): ~{power_law_model(50, a, b):.0f} ms")
        print(f"  For large molecules (~100 atoms): ~{power_law_model(100, a, b):.0f} ms")
    elif best[0] == 'linear':
        a, b = best[3]
        print(f"\nExecution time = {a:.3f} × atom_count + {b:.1f}")
        print(f"\nInterpretation: Time scales linearly with molecular size")
        print(f"  For small molecules (~20 atoms): ~{linear_model(20, a, b):.0f} ms")
        print(f"  For medium molecules (~50 atoms): ~{linear_model(50, a, b):.0f} ms")
        print(f"  For large molecules (~100 atoms): ~{linear_model(100, a, b):.0f} ms")
    elif best[0] == 'quadratic':
        a, b, c = best[3]
        print(f"\nExecution time = {a:.5f} × atom_count² + {b:.3f} × atom_count + {c:.1f}")
        print(f"\nInterpretation: Time scales quadratically with molecular size")
        print(f"  For small molecules (~20 atoms): ~{quadratic_model(20, a, b, c):.0f} ms")
        print(f"  For medium molecules (~50 atoms): ~{quadratic_model(50, a, b, c):.0f} ms")
        print(f"  For large molecules (~100 atoms): ~{quadratic_model(100, a, b, c):.0f} ms")

# Statistics by size ranges
print(f"\n{'='*70}")
print("SIZE-BASED STATISTICS")
print(f"{'='*70}")

size_ranges = [
    ("Small (<25 atoms)", 0, 25),
    ("Medium (25-50 atoms)", 25, 50),
    ("Medium-Large (50-100 atoms)", 50, 100),
    ("Large (>100 atoms)", 100, 200)
]

for label, low, high in size_ranges:
    mask = (atom_counts >= low) & (atom_counts < high)
    if np.sum(mask) > 0:
        times = exec_times[mask]
        print(f"\n{label}: {np.sum(mask)} measurements")
        print(f"  Mean: {times.mean():.1f} ms")
        print(f"  Median: {np.median(times):.1f} ms")
        print(f"  Range: {times.min():.1f}-{times.max():.1f} ms")
        print(f"  Std dev: {times.std():.1f} ms")
