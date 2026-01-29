#!/usr/bin/env python3
"""
Run timing benchmark on Windows and compare with Linux results.

This script:
1. Runs timing analysis on the current Windows machine
2. Loads the Linux timing data from timing_analysis.json
3. Compares Windows vs Linux performance
4. Generates comparative analysis plots

Author: PFASGroups Benchmark Team
Date: 2026-01-29
"""

import json
import time
import platform
import sys
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy import stats

# Add PFASGroups to path
sys.path.insert(0, str(Path(__file__).parent.parent.parent))

from PFASgroups import parse_smiles


def exponential_model(n, a, b):
    """Exponential model: t = a * exp(b * n)"""
    return a * np.exp(b * n)


def run_timing_benchmark(test_molecules, n_iterations=10):
    """
    Run timing benchmark on current machine.
    
    Parameters
    ----------
    test_molecules : list
        List of dicts with 'smiles' and 'num_atoms' keys
    n_iterations : int
        Number of iterations per molecule
        
    Returns
    -------
    list
        List of timing results with statistics
    """
    results = []
    
    print(f"\n{'='*60}")
    print(f"Running timing benchmark on {platform.system()} {platform.release()}")
    print(f"Processor: {platform.processor()}")
    print(f"Python: {platform.python_version()}")
    print(f"{'='*60}\n")
    
    for i, mol_data in enumerate(test_molecules, 1):
        smiles = mol_data['smiles']
        n_atoms = mol_data['num_atoms']
        
        print(f"[{i}/{len(test_molecules)}] Testing {n_atoms} atoms: {smiles[:50]}...")
        
        # Run iterations
        times = []
        for j in range(n_iterations):
            start_time = time.perf_counter()
            try:
                result = parse_smiles(smiles)
                elapsed = (time.perf_counter() - start_time) * 1000  # Convert to ms
                times.append(elapsed)
            except Exception as e:
                print(f"  ✗ Error on iteration {j+1}: {e}")
                continue
        
        if times:
            mean_time = np.mean(times)
            std_time = np.std(times, ddof=1)
            min_time = np.min(times)
            max_time = np.max(times)
            
            results.append({
                'smiles': smiles,
                'num_atoms': n_atoms,
                'times_ms': times,
                'mean_ms': mean_time,
                'std_ms': std_time,
                'min_ms': min_time,
                'max_ms': max_time,
                'n_iterations': len(times)
            })
            
            print(f"  ✓ Mean: {mean_time:.2f} ± {std_time:.2f} ms (n={len(times)})")
        else:
            print(f"  ✗ All iterations failed")
    
    return results


def load_linux_data(linux_file='data/timing_analysis.json'):
    """Load Linux timing data from JSON file"""
    script_dir = Path(__file__).parent.parent
    linux_path = script_dir / linux_file
    
    if not linux_path.exists():
        print(f"Warning: Linux data file not found: {linux_path}")
        return None
    
    with open(linux_path, 'r') as f:
        data = json.load(f)
    
    return data


def compare_platforms(windows_results, linux_results):
    """
    Compare Windows vs Linux timing performance.
    
    Parameters
    ----------
    windows_results : list
        Windows timing results
    linux_results : list
        Linux timing results
        
    Returns
    -------
    dict
        Comparison statistics
    """
    # Match molecules by atom count
    comparison = []
    
    for win_res in windows_results:
        n_atoms = win_res['num_atoms']
        
        # Find corresponding Linux result
        linux_match = next((lr for lr in linux_results 
                           if lr['num_atoms'] == n_atoms), None)
        
        if linux_match:
            speedup = win_res['mean_ms'] / linux_match['mean_ms']
            comparison.append({
                'num_atoms': n_atoms,
                'windows_ms': win_res['mean_ms'],
                'linux_ms': linux_match['mean_ms'],
                'speedup': speedup,  # > 1 means Linux faster
                'windows_std': win_res['std_ms'],
                'linux_std': linux_match['std_ms']
            })
    
    if not comparison:
        return None
    
    # Overall statistics
    speedups = [c['speedup'] for c in comparison]
    
    stats_dict = {
        'n_molecules': len(comparison),
        'mean_speedup': np.mean(speedups),
        'median_speedup': np.median(speedups),
        'std_speedup': np.std(speedups),
        'min_speedup': np.min(speedups),
        'max_speedup': np.max(speedups),
        'comparison': comparison
    }
    
    return stats_dict


def plot_comparison(windows_results, linux_results, comparison_stats, output_dir='imgs'):
    """Generate comparison plots"""
    script_dir = Path(__file__).parent.parent
    output_path = script_dir / output_dir
    output_path.mkdir(parents=True, exist_ok=True)
    
    # Extract data
    win_atoms = np.array([r['num_atoms'] for r in windows_results])
    win_times = np.array([r['mean_ms'] for r in windows_results])
    win_stds = np.array([r['std_ms'] for r in windows_results])
    
    linux_atoms = np.array([r['num_atoms'] for r in linux_results])
    linux_times = np.array([r['mean_ms'] for r in linux_results])
    linux_stds = np.array([r['std_ms'] for r in linux_results])
    
    # Fit exponential models
    popt_win, _ = curve_fit(exponential_model, win_atoms, win_times, 
                             p0=[50, 0.015], maxfev=10000)
    popt_linux, _ = curve_fit(exponential_model, linux_atoms, linux_times,
                               p0=[50, 0.015], maxfev=10000)
    
    # Generate smooth curves
    atoms_smooth = np.linspace(win_atoms.min(), win_atoms.max(), 300)
    win_fit = exponential_model(atoms_smooth, *popt_win)
    linux_fit = exponential_model(atoms_smooth, *popt_linux)
    
    # Create figure
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    fig.suptitle('Windows vs Linux Performance Comparison', fontsize=16, fontweight='bold')
    
    # 1. Direct timing comparison
    ax = axes[0, 0]
    ax.errorbar(win_atoms, win_times, yerr=win_stds, fmt='o', 
                label='Windows', alpha=0.6, capsize=3, markersize=4)
    ax.errorbar(linux_atoms, linux_times, yerr=linux_stds, fmt='s',
                label='Linux', alpha=0.6, capsize=3, markersize=4)
    ax.plot(atoms_smooth, win_fit, '-', color='C0', linewidth=2,
            label=f'Windows: t={popt_win[0]:.2f}×exp({popt_win[1]:.5f}×n)')
    ax.plot(atoms_smooth, linux_fit, '-', color='C1', linewidth=2,
            label=f'Linux: t={popt_linux[0]:.2f}×exp({popt_linux[1]:.5f}×n)')
    ax.axvline(100, color='red', linestyle='--', alpha=0.5, linewidth=2,
               label='100-atom threshold')
    ax.set_xlabel('Number of Atoms', fontsize=12)
    ax.set_ylabel('Execution Time (ms)', fontsize=12)
    ax.set_title('Timing Comparison', fontsize=13, fontweight='bold')
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3)
    
    # 2. Speedup ratio
    ax = axes[0, 1]
    if comparison_stats:
        comp_atoms = [c['num_atoms'] for c in comparison_stats['comparison']]
        speedups = [c['speedup'] for c in comparison_stats['comparison']]
        
        ax.scatter(comp_atoms, speedups, alpha=0.6, s=50)
        ax.axhline(1, color='black', linestyle='--', linewidth=2,
                   label='Equal performance')
        ax.axhline(comparison_stats['mean_speedup'], color='red', 
                   linestyle='-', linewidth=2,
                   label=f"Mean: {comparison_stats['mean_speedup']:.2f}×")
        ax.axvline(100, color='green', linestyle='--', alpha=0.5, linewidth=2)
        ax.set_xlabel('Number of Atoms', fontsize=12)
        ax.set_ylabel('Speedup (Windows/Linux)', fontsize=12)
        ax.set_title('Platform Speedup Ratio (>1 = Linux faster)', fontsize=13, fontweight='bold')
        ax.legend(fontsize=10)
        ax.grid(True, alpha=0.3)
    
    # 3. Time per atom
    ax = axes[1, 0]
    win_time_per_atom = win_times / win_atoms
    linux_time_per_atom = linux_times / linux_atoms
    
    ax.scatter(win_atoms, win_time_per_atom, label='Windows', alpha=0.6, s=50)
    ax.scatter(linux_atoms, linux_time_per_atom, label='Linux', alpha=0.6, s=50)
    ax.axvline(100, color='red', linestyle='--', alpha=0.5, linewidth=2,
               label='100-atom threshold')
    ax.set_xlabel('Number of Atoms', fontsize=12)
    ax.set_ylabel('Time per Atom (ms/atom)', fontsize=12)
    ax.set_title('Efficiency Comparison', fontsize=13, fontweight='bold')
    ax.legend(fontsize=10)
    ax.grid(True, alpha=0.3)
    
    # 4. Statistical summary
    ax = axes[1, 1]
    ax.axis('off')
    
    # Calculate R² for both platforms
    win_residuals = win_times - exponential_model(win_atoms, *popt_win)
    win_ss_res = np.sum(win_residuals**2)
    win_ss_tot = np.sum((win_times - np.mean(win_times))**2)
    win_r2 = 1 - (win_ss_res / win_ss_tot)
    
    linux_residuals = linux_times - exponential_model(linux_atoms, *popt_linux)
    linux_ss_res = np.sum(linux_residuals**2)
    linux_ss_tot = np.sum((linux_times - np.mean(linux_times))**2)
    linux_r2 = 1 - (linux_ss_res / linux_ss_tot)
    
    summary_text = f"""
    Platform Comparison Summary
    {'='*45}
    
    Windows ({platform.system()} {platform.release()})
      • Processor: {platform.processor()}
      • Model: t = {popt_win[0]:.2f} × exp({popt_win[1]:.5f} × n)
      • α = {popt_win[1]:.5f} atoms⁻¹
      • R² = {win_r2:.4f}
      • Mean time: {np.mean(win_times):.2f} ± {np.std(win_times):.2f} ms
    
    Linux
      • Model: t = {popt_linux[0]:.2f} × exp({popt_linux[1]:.5f} × n)
      • α = {popt_linux[1]:.5f} atoms⁻¹
      • R² = {linux_r2:.4f}
      • Mean time: {np.mean(linux_times):.2f} ± {np.std(linux_times):.2f} ms
    """
    
    if comparison_stats:
        summary_text += f"""
    Platform Speedup (Windows / Linux)
      • Mean: {comparison_stats['mean_speedup']:.2f}×
      • Median: {comparison_stats['median_speedup']:.2f}×
      • Range: {comparison_stats['min_speedup']:.2f}× - {comparison_stats['max_speedup']:.2f}×
      
    {'Linux is faster' if comparison_stats['mean_speedup'] > 1 else 'Windows is faster'}
    by an average factor of {abs(comparison_stats['mean_speedup'] - 1):.2f}×
    """
    
    ax.text(0.05, 0.95, summary_text, transform=ax.transAxes,
            fontsize=10, verticalalignment='top', family='monospace',
            bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
    
    plt.tight_layout()
    
    # Save figure
    output_pdf = output_path / 'windows_vs_linux_comparison.pdf'
    output_png = output_path / 'windows_vs_linux_comparison.png'
    
    plt.savefig(output_pdf, dpi=300, bbox_inches='tight')
    plt.savefig(output_png, dpi=300, bbox_inches='tight')
    print(f"\n✓ Saved comparison plots:")
    print(f"  {output_pdf}")
    print(f"  {output_png}")
    
    plt.close()


def main():
    """Main execution"""
    script_dir = Path(__file__).parent.parent
    
    # Load timing dataset
    print("Loading timing dataset...")
    timing_data_file = script_dir / 'data' / 'timing_full_dataset.npz'
    
    if not timing_data_file.exists():
        print(f"Error: Timing dataset not found: {timing_data_file}")
        return
    
    # Load the dataset
    data = np.load(timing_data_file, allow_pickle=True)
    atom_counts = data['atom_counts']
    smiles_list = data['smiles_list']
    
    # Get unique test molecules (one per atom count)
    unique_molecules = {}
    for smiles, n_atoms in zip(smiles_list, atom_counts):
        if n_atoms not in unique_molecules:
            unique_molecules[n_atoms] = smiles
    
    # Create test molecules list
    test_molecules = [
        {'smiles': smiles, 'num_atoms': n_atoms}
        for n_atoms, smiles in sorted(unique_molecules.items())
    ]
    
    print(f"Found {len(test_molecules)} unique molecules (atom range: {min(atom_counts)}-{max(atom_counts)})")
    
    # Run timing benchmark on Windows
    print("\nRunning Windows timing benchmark...")
    windows_results = run_timing_benchmark(test_molecules, n_iterations=10)
    
    # Save Windows results
    windows_output = {
        'platform': platform.system(),
        'platform_release': platform.release(),
        'processor': platform.processor(),
        'python_version': platform.python_version(),
        'timestamp': time.strftime('%Y-%m-%d %H:%M:%S'),
        'n_molecules': len(windows_results),
        'n_iterations': 10,
        'results': windows_results
    }
    
    output_file = script_dir / 'data' / 'timing_analysis_windows.json'
    output_file.parent.mkdir(parents=True, exist_ok=True)
    
    with open(output_file, 'w') as f:
        json.dump(windows_output, f, indent=2)
    
    print(f"\n✓ Saved Windows timing results to {output_file}")
    
    # Compare platforms - need to match with Linux data  
    print("\nLoading Linux timing data for comparison...")
    linux_data_file = script_dir / 'data' / 'timing_full_dataset.npz'
    linux_data_loaded = np.load(linux_data_file, allow_pickle=True)
    linux_atom_counts = linux_data_loaded['atom_counts']
    linux_times = linux_data_loaded['times']
    
    # Calculate Linux statistics per molecule (mean of iterations)
    linux_results = []
    for n_atoms in sorted(set(linux_atom_counts)):
        mask = linux_atom_counts == n_atoms
        times_for_atom = linux_times[mask]
        linux_results.append({
            'num_atoms': n_atoms,
            'mean_ms': np.mean(times_for_atom),
            'std_ms': np.std(times_for_atom, ddof=1),
            'min_ms': np.min(times_for_atom),
            'max_ms': np.max(times_for_atom)
        })
    
    print("\nComparing Windows vs Linux performance...")
    comparison_stats = compare_platforms(windows_results, linux_results)
    
    if comparison_stats:
        print(f"\nPlatform Performance Comparison:")
        print(f"  Mean speedup (Windows/Linux): {comparison_stats['mean_speedup']:.2f}×")
        print(f"  Median speedup: {comparison_stats['median_speedup']:.2f}×")
        print(f"  Speedup range: {comparison_stats['min_speedup']:.2f}× - {comparison_stats['max_speedup']:.2f}×")
        
        if comparison_stats['mean_speedup'] > 1.1:
            print(f"\n  → Linux is ~{(comparison_stats['mean_speedup']-1)*100:.0f}% faster on average")
        elif comparison_stats['mean_speedup'] < 0.9:
            print(f"\n  → Windows is ~{(1/comparison_stats['mean_speedup']-1)*100:.0f}% faster on average")
        else:
            print(f"\n  → Both platforms have similar performance")
        
        # Save comparison results
        comparison_file = script_dir / 'data' / 'platform_comparison.json'
        with open(comparison_file, 'w') as f:
            json.dump(comparison_stats, f, indent=2)
        
        print(f"\n✓ Saved platform comparison to {comparison_file}")
    
    # Generate comparison plots
    print("\nGenerating comparison plots...")
    plot_comparison(windows_results, linux_results, comparison_stats)
    
    print("\n" + "="*60)
    print("Windows timing benchmark complete!")
    print("="*60)


if __name__ == '__main__':
    main()
