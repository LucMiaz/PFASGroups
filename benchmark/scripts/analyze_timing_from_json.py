import json
import numpy as np
from scipy.optimize import curve_fit

# Load timing analysis JSON
with open('../reports/timing_analysis.json', 'r') as f:
    data = json.load(f)

print("Timing Analysis Data:")
print(f"Total molecules: {data['summary']['total_molecules']}")
print(f"PFASgroups avg: {data['summary']['pfasgroups_avg_time']:.1f} ms")
print(f"PFASgroups std: {data['summary']['pfasgroups_std_time']:.1f} ms")
print(f"Atlas avg: {data['summary']['atlas_avg_time']:.1f} ms")
print(f"Atlas std: {data['summary']['atlas_std_time']:.1f} ms")
print(f"Speed ratio: {data['summary']['speed_ratio']:.3f}")
print(f"Iterations: {data['summary']['iterations']}")

# Try to find the raw data file
import os
import glob

# Look for timing data files
pattern = '../reports/timing_data_*.json'
files = glob.glob(pattern)

if files:
    print(f"\nFound {len(files)} timing data files")
    # Use the most recent one (matching timestamp)
    timestamp = data['summary']['timestamp']
    data_file = f'../reports/timing_data_{timestamp}.json'
    
    if os.path.exists(data_file):
        print(f"Loading: {data_file}")
        with open(data_file, 'r') as f:
            timing_data = json.load(f)
        
        # Extract molecule properties and times
        molecules = timing_data.get('molecules', [])
        print(f"\nFound {len(molecules)} molecules with timing data")
        
        if molecules:
            atom_counts = []
            pfas_times = []
            atlas_times = []
            
            for mol in molecules:
                if 'num_atoms' in mol and 'pfasgroups_time' in mol and 'atlas_time' in mol:
                    atom_counts.append(mol['num_atoms'])
                    pfas_times.append(mol['pfasgroups_time'] * 1000)  # Convert to ms
                    atlas_times.append(mol['atlas_time'] * 1000)
            
            if atom_counts:
                atom_counts = np.array(atom_counts)
                pfas_times = np.array(pfas_times)
                atlas_times = np.array(atlas_times)
                
                print(f"\nAtom count range: {atom_counts.min()}-{atom_counts.max()}")
                print(f"PFASgroups time range: {pfas_times.min():.1f}-{pfas_times.max():.1f} ms")
                
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
                
                # Fit models
                print("\n" + "="*60)
                print("PFASGROUPS SCALING ANALYSIS")
                print("="*60)
                
                try:
                    popt_exp, _ = curve_fit(exponential_model, atom_counts, pfas_times, 
                                           p0=[10, 0.01], maxfev=10000)
                    pred_exp = exponential_model(atom_counts, *popt_exp)
                    r2_exp = 1 - (np.sum((pfas_times - pred_exp)**2) / 
                                 np.sum((pfas_times - np.mean(pfas_times))**2))
                    rmse_exp = np.sqrt(np.mean((pfas_times - pred_exp)**2))
                    
                    print(f"\nExponential: time = {popt_exp[0]:.3f} * exp({popt_exp[1]:.5f} * atoms)")
                    print(f"  R² = {r2_exp:.4f}, RMSE = {rmse_exp:.1f} ms")
                except Exception as e:
                    print(f"\nExponential fit failed: {e}")
                    r2_exp = -np.inf
                
                try:
                    popt_pow, _ = curve_fit(power_law_model, atom_counts, pfas_times, 
                                           p0=[1, 1.5], maxfev=10000)
                    pred_pow = power_law_model(atom_counts, *popt_pow)
                    r2_pow = 1 - (np.sum((pfas_times - pred_pow)**2) / 
                                 np.sum((pfas_times - np.mean(pfas_times))**2))
                    rmse_pow = np.sqrt(np.mean((pfas_times - pred_pow)**2))
                    
                    print(f"\nPower Law: time = {popt_pow[0]:.3f} * atoms^{popt_pow[1]:.3f}")
                    print(f"  R² = {r2_pow:.4f}, RMSE = {rmse_pow:.1f} ms")
                except Exception as e:
                    print(f"\nPower law fit failed: {e}")
                    r2_pow = -np.inf
                
                try:
                    popt_lin, _ = curve_fit(linear_model, atom_counts, pfas_times)
                    pred_lin = linear_model(atom_counts, *popt_lin)
                    r2_lin = 1 - (np.sum((pfas_times - pred_lin)**2) / 
                                 np.sum((pfas_times - np.mean(pfas_times))**2))
                    rmse_lin = np.sqrt(np.mean((pfas_times - pred_lin)**2))
                    
                    print(f"\nLinear: time = {popt_lin[0]:.3f} * atoms + {popt_lin[1]:.1f}")
                    print(f"  R² = {r2_lin:.4f}, RMSE = {rmse_lin:.1f} ms")
                except Exception as e:
                    print(f"\nLinear fit failed: {e}")
                    r2_lin = -np.inf
                
                # Statistics by size bins
                print("\n" + "="*60)
                print("STATISTICS BY SIZE BINS")
                print("="*60)
                
                bins = [(0, 20), (20, 40), (40, 60), (60, 100)]
                for low, high in bins:
                    mask = (atom_counts >= low) & (atom_counts < high)
                    if np.sum(mask) > 0:
                        print(f"\n{low}-{high} atoms ({np.sum(mask)} molecules):")
                        print(f"  PFASgroups: {pfas_times[mask].mean():.1f} ± {pfas_times[mask].std():.1f} ms")
                        print(f"  Atlas:      {atlas_times[mask].mean():.1f} ± {atlas_times[mask].std():.1f} ms")
                        print(f"  Ratio:      {pfas_times[mask].mean() / atlas_times[mask].mean():.2f}×")
else:
    print("\nNo detailed timing data files found. Using summary statistics only.")
