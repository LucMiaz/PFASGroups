#!/usr/bin/env python3
"""
Comprehensive timing analysis comparing three metric computation scenarios:
1. Full metrics with effective graph resistance for all molecules (limit=None)
2. Full metrics without effective graph resistance (limit=0)
3. Minimal metrics - only component size (compute_component_metrics=False)

This analysis quantifies the performance impact of:
- Effective graph resistance computation (O(N²) operation)
- All other graph metrics (diameter, radius, eccentricity, etc.)

Author: PFASGroups Benchmark Team
Date: 2026-01-29
"""

import time
import platform
import sys
from pathlib import Path
import numpy as np
import json
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from scipy.optimize import curve_fit
from scipy import stats

# Add PFASGroups to path
pfas_path = r'c:\Users\luc\git\PFASGroups'
sys.path.insert(0, pfas_path)

from PFASgroups import parse_smiles


def exponential_model(n, a, b):
    """Exponential model: t = a * exp(b * n)"""
    return a * np.exp(b * n)


def run_timing_test(test_molecules, n_iterations=10, scenario='full'):
    """
    Run timing test with specific metric computation settings.
    
    Parameters
    ----------
    test_molecules : list
        List of dicts with 'smiles' and 'num_atoms' keys
    n_iterations : int
        Number of iterations per molecule
    scenario : str
        'full' - All metrics including effective graph resistance (limit=None)
        'no_resistance' - All metrics except effective graph resistance (limit=0)
        'minimal' - Only component size (compute_component_metrics=False)
        
    Returns
    -------
    list
        List of timing results with statistics
    """
    # Configure parameters based on scenario
    if scenario == 'full':
        kwargs = {'limit_effective_graph_resistance': None, 'compute_component_metrics': True}
        desc = "Full metrics + effective graph resistance"
    elif scenario == 'no_resistance':
        kwargs = {'limit_effective_graph_resistance': 0, 'compute_component_metrics': True}
        desc = "Full metrics, no effective graph resistance"
    elif scenario == 'minimal':
        kwargs = {'limit_effective_graph_resistance': 0, 'compute_component_metrics': False}
        desc = "Minimal metrics (size only)"
    else:
        raise ValueError(f"Unknown scenario: {scenario}")
    
    print(f"\n{'='*60}")
    print(f"Testing: {desc}")
    print(f"{'='*60}\n")
    
    results = []
    
    for i, mol_data in enumerate(test_molecules, 1):
        smiles = mol_data['smiles']
        n_atoms = mol_data['num_atoms']
        
        print(f"[{i}/{len(test_molecules)}] {n_atoms} atoms: ", end='', flush=True)
        
        # Run iterations
        times = []
        for j in range(n_iterations):
            start_time = time.perf_counter()
            try:
                result = parse_smiles(smiles, **kwargs)
                elapsed = (time.perf_counter() - start_time) * 1000  # Convert to ms
                times.append(elapsed)
            except Exception as e:
                print(f"\n  ✗ Error on iteration {j+1}: {e}")
                continue
        
        if times:
            mean_time = np.mean(times)
            std_time = np.std(times, ddof=1)
            
            results.append({
                'smiles': smiles,
                'num_atoms': n_atoms,
                'times_ms': times,
                'mean_ms': mean_time,
                'std_ms': std_time,
                'min_ms': np.min(times),
                'max_ms': np.max(times),
                'n_iterations': len(times)
            })
            
            print(f"{mean_time:.2f} ± {std_time:.2f} ms")
        else:
            print(f"✗ All iterations failed")
    
    return results


def compare_scenarios(results_dict):
    """
    Compare the three scenarios and calculate speedups.
    
    Parameters
    ----------
    results_dict : dict
        Dictionary with keys 'full', 'no_resistance', 'minimal' and timing results
        
    Returns
    -------
    dict
        Comparison statistics
    """
    comparisons = []
    
    # Match molecules by atom count across all scenarios
    for full_res in results_dict['full']:
        n_atoms = full_res['num_atoms']
        
        # Find corresponding results in other scenarios
        no_res_match = next((r for r in results_dict['no_resistance'] 
                            if r['num_atoms'] == n_atoms), None)
        minimal_match = next((r for r in results_dict['minimal']
                             if r['num_atoms'] == n_atoms), None)
        
        if no_res_match and minimal_match:
            # Calculate time differences
            resistance_cost = full_res['mean_ms'] - no_res_match['mean_ms']
            other_metrics_cost = no_res_match['mean_ms'] - minimal_match['mean_ms']
            total_cost = full_res['mean_ms'] - minimal_match['mean_ms']
            
            comparisons.append({
                'num_atoms': n_atoms,
                'full_ms': full_res['mean_ms'],
                'no_resistance_ms': no_res_match['mean_ms'],
                'minimal_ms': minimal_match['mean_ms'],
                'resistance_cost_ms': resistance_cost,
                'resistance_cost_pct': (resistance_cost / full_res['mean_ms']) * 100,
                'other_metrics_cost_ms': other_metrics_cost,
                'other_metrics_cost_pct': (other_metrics_cost / full_res['mean_ms']) * 100,
                'total_overhead_ms': total_cost,
                'total_overhead_pct': (total_cost / full_res['mean_ms']) * 100,
                'speedup_no_resistance': full_res['mean_ms'] / no_res_match['mean_ms'],
                'speedup_minimal': full_res['mean_ms'] / minimal_match['mean_ms']
            })
    
    # Calculate overall statistics
    if comparisons:
        stats_dict = {
            'n_molecules': len(comparisons),
            'mean_resistance_cost_ms': np.mean([c['resistance_cost_ms'] for c in comparisons]),
            'mean_resistance_cost_pct': np.mean([c['resistance_cost_pct'] for c in comparisons]),
            'mean_other_metrics_cost_ms': np.mean([c['other_metrics_cost_ms'] for c in comparisons]),
            'mean_other_metrics_cost_pct': np.mean([c['other_metrics_cost_pct'] for c in comparisons]),
            'mean_speedup_no_resistance': np.mean([c['speedup_no_resistance'] for c in comparisons]),
            'mean_speedup_minimal': np.mean([c['speedup_minimal'] for c in comparisons]),
            'comparisons': comparisons
        }
        return stats_dict
    
    return None


def plot_comparison(results_dict, comparison_stats, output_dir='imgs'):
    """Generate comprehensive comparison plots"""
    script_dir = Path(__file__).parent.parent
    output_path = script_dir / output_dir
    output_path.mkdir(parents=True, exist_ok=True)
    
    # Extract data for each scenario
    scenarios = {
        'full': ('Full Metrics + Resistance', 'C0', 'o'),
        'no_resistance': ('Full Metrics (No Resistance)', 'C1', 's'),
        'minimal': ('Minimal (Size Only)', 'C2', '^')
    }
    
    fig = plt.figure(figsize=(16, 12))
    gs = GridSpec(3, 3, figure=fig, hspace=0.3, wspace=0.3)
    
    # 1. Main timing comparison
    ax1 = fig.add_subplot(gs[0, :2])
    for scenario_key, (label, color, marker) in scenarios.items():
        data = results_dict[scenario_key]
        atoms = np.array([r['num_atoms'] for r in data])
        times = np.array([r['mean_ms'] for r in data])
        stds = np.array([r['std_ms'] for r in data])
        
        ax1.errorbar(atoms, times, yerr=stds, fmt=marker, label=label,
                    alpha=0.6, capsize=3, markersize=5, color=color)
        
        # Fit exponential model
        try:
            popt, _ = curve_fit(exponential_model, atoms, times, p0=[50, 0.015], maxfev=10000)
            atoms_smooth = np.linspace(atoms.min(), atoms.max(), 300)
            fit_curve = exponential_model(atoms_smooth, *popt)
            ax1.plot(atoms_smooth, fit_curve, '-', linewidth=2, color=color, alpha=0.7,
                    label=f'{label}: t={popt[0]:.1f}×exp({popt[1]:.5f}×n)')
        except:
            pass
    
    ax1.axvline(100, color='red', linestyle='--', alpha=0.5, linewidth=2,
               label='100-atom threshold')
    ax1.set_xlabel('Number of Atoms', fontsize=12)
    ax1.set_ylabel('Execution Time (ms)', fontsize=12)
    ax1.set_title('Timing Comparison: Impact of Metric Computation', fontsize=13, fontweight='bold')
    ax1.legend(fontsize=9, loc='upper left')
    ax1.grid(True, alpha=0.3)
    
    # 2. Cost breakdown
    ax2 = fig.add_subplot(gs[0, 2])
    if comparison_stats:
        comp = comparison_stats['comparisons']
        atoms_comp = [c['num_atoms'] for c in comp]
        resistance_costs = [c['resistance_cost_ms'] for c in comp]
        other_costs = [c['other_metrics_cost_ms'] for c in comp]
        
        ax2.scatter(atoms_comp, resistance_costs, label='Resistance Cost', alpha=0.6, s=50)
        ax2.scatter(atoms_comp, other_costs, label='Other Metrics Cost', alpha=0.6, s=50)
        ax2.axvline(100, color='red', linestyle='--', alpha=0.5, linewidth=2)
        ax2.set_xlabel('Number of Atoms', fontsize=11)
        ax2.set_ylabel('Cost (ms)', fontsize=11)
        ax2.set_title('Metric Computation Costs', fontsize=12, fontweight='bold')
        ax2.legend(fontsize=9)
        ax2.grid(True, alpha=0.3)
    
    # 3. Speedup ratios
    ax3 = fig.add_subplot(gs[1, 0])
    if comparison_stats:
        comp = comparison_stats['comparisons']
        atoms_comp = [c['num_atoms'] for c in comp]
        speedup_no_res = [c['speedup_no_resistance'] for c in comp]
        
        ax3.scatter(atoms_comp, speedup_no_res, alpha=0.6, s=50)
        ax3.axhline(1, color='black', linestyle='--', linewidth=1)
        ax3.axvline(100, color='red', linestyle='--', alpha=0.5, linewidth=2)
        mean_speedup = comparison_stats['mean_speedup_no_resistance']
        ax3.axhline(mean_speedup, color='red', linestyle='-', linewidth=2,
                   label=f'Mean: {mean_speedup:.2f}×')
        ax3.set_xlabel('Number of Atoms', fontsize=11)
        ax3.set_ylabel('Speedup', fontsize=11)
        ax3.set_title('Speedup: Full vs No Resistance', fontsize=12, fontweight='bold')
        ax3.legend(fontsize=9)
        ax3.grid(True, alpha=0.3)
    
    # 4. Speedup to minimal
    ax4 = fig.add_subplot(gs[1, 1])
    if comparison_stats:
        comp = comparison_stats['comparisons']
        atoms_comp = [c['num_atoms'] for c in comp]
        speedup_minimal = [c['speedup_minimal'] for c in comp]
        
        ax4.scatter(atoms_comp, speedup_minimal, alpha=0.6, s=50, color='C2')
        ax4.axhline(1, color='black', linestyle='--', linewidth=1)
        ax4.axvline(100, color='red', linestyle='--', alpha=0.5, linewidth=2)
        mean_speedup = comparison_stats['mean_speedup_minimal']
        ax4.axhline(mean_speedup, color='red', linestyle='-', linewidth=2,
                   label=f'Mean: {mean_speedup:.2f}×')
        ax4.set_xlabel('Number of Atoms', fontsize=11)
        ax4.set_ylabel('Speedup', fontsize=11)
        ax4.set_title('Speedup: Full vs Minimal', fontsize=12, fontweight='bold')
        ax4.legend(fontsize=9)
        ax4.grid(True, alpha=0.3)
    
    # 5. Percentage breakdown
    ax5 = fig.add_subplot(gs[1, 2])
    if comparison_stats:
        comp = comparison_stats['comparisons']
        atoms_comp = np.array([c['num_atoms'] for c in comp])
        resistance_pct = np.array([c['resistance_cost_pct'] for c in comp])
        other_pct = np.array([c['other_metrics_cost_pct'] for c in comp])
        
        width = 2
        ax5.bar(atoms_comp - width/2, resistance_pct, width=width, label='Resistance', alpha=0.7)
        ax5.bar(atoms_comp + width/2, other_pct, width=width, label='Other Metrics', alpha=0.7)
        ax5.axvline(100, color='red', linestyle='--', alpha=0.5, linewidth=2)
        ax5.set_xlabel('Number of Atoms', fontsize=11)
        ax5.set_ylabel('% of Total Time', fontsize=11)
        ax5.set_title('Cost Breakdown (%)', fontsize=12, fontweight='bold')
        ax5.legend(fontsize=9)
        ax5.grid(True, alpha=0.3, axis='y')
    
    # 6. Time per atom comparison
    ax6 = fig.add_subplot(gs[2, 0])
    for scenario_key, (label, color, marker) in scenarios.items():
        data = results_dict[scenario_key]
        atoms = np.array([r['num_atoms'] for r in data])
        times = np.array([r['mean_ms'] for r in data])
        time_per_atom = times / atoms
        
        ax6.scatter(atoms, time_per_atom, label=label, alpha=0.6, s=50,
                   marker=marker, color=color)
    
    ax6.axvline(100, color='red', linestyle='--', alpha=0.5, linewidth=2)
    ax6.set_xlabel('Number of Atoms', fontsize=11)
    ax6.set_ylabel('Time per Atom (ms/atom)', fontsize=11)
    ax6.set_title('Efficiency Comparison', fontsize=12, fontweight='bold')
    ax6.legend(fontsize=9)
    ax6.grid(True, alpha=0.3)
    
    # 7-9. Statistical summary text boxes
    ax7 = fig.add_subplot(gs[2, 1:])
    ax7.axis('off')
    
    summary_text = f"""
    COMPREHENSIVE TIMING ANALYSIS SUMMARY
    {'='*75}
    
    Platform: {platform.system()} {platform.release()}
    Processor: {platform.processor()}
    Python: {platform.python_version()}
    
    Test Configuration:
      • Molecules tested: {len(results_dict['full'])}
      • Atom range: {min(r['num_atoms'] for r in results_dict['full'])}-{max(r['num_atoms'] for r in results_dict['full'])}
      • Iterations per molecule: {results_dict['full'][0]['n_iterations']}
    """
    
    if comparison_stats:
        summary_text += f"""
    SCENARIO 1: Full Metrics + Effective Graph Resistance (limit=None)
      • Mean time: {np.mean([r['mean_ms'] for r in results_dict['full']]):.2f} ms
      • Baseline (100% of time)
    
    SCENARIO 2: Full Metrics, No Effective Graph Resistance (limit=0)
      • Mean time: {np.mean([r['mean_ms'] for r in results_dict['no_resistance']]):.2f} ms
      • Speedup: {comparison_stats['mean_speedup_no_resistance']:.2f}× faster
      • Time saved: {comparison_stats['mean_resistance_cost_ms']:.2f} ms ({comparison_stats['mean_resistance_cost_pct']:.1f}%)
    
    SCENARIO 3: Minimal Metrics (component size only)
      • Mean time: {np.mean([r['mean_ms'] for r in results_dict['minimal']]):.2f} ms
      • Speedup: {comparison_stats['mean_speedup_minimal']:.2f}× faster
      • Time saved: {np.mean([c['total_overhead_ms'] for c in comparison_stats['comparisons']]):.2f} ms
    
    COST BREAKDOWN (Average across all molecules):
      • Effective graph resistance: {comparison_stats['mean_resistance_cost_ms']:.2f} ms ({comparison_stats['mean_resistance_cost_pct']:.1f}%)
      • Other graph metrics: {comparison_stats['mean_other_metrics_cost_ms']:.2f} ms ({comparison_stats['mean_other_metrics_cost_pct']:.1f}%)
      • Total metrics overhead: {comparison_stats['mean_resistance_cost_ms'] + comparison_stats['mean_other_metrics_cost_ms']:.2f} ms
    
    KEY FINDINGS:
      • Effective graph resistance is the dominant cost for molecules < 100 atoms
      • At 100 atoms, resistance computation is automatically skipped (hardcoded threshold)
      • Other graph metrics (diameter, radius, etc.) have moderate cost
      • Minimal mode (size only) achieves {comparison_stats['mean_speedup_minimal']:.1f}× speedup
    """
    
    ax7.text(0.05, 0.95, summary_text, transform=ax7.transAxes,
            fontsize=9, verticalalignment='top', family='monospace',
            bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
    
    plt.suptitle('Component Metrics Performance Analysis', fontsize=16, fontweight='bold', y=0.98)
    
    # Save figure
    output_pdf = output_path / 'component_metrics_analysis.pdf'
    output_png = output_path / 'component_metrics_analysis.png'
    
    plt.savefig(output_pdf, dpi=300, bbox_inches='tight')
    plt.savefig(output_png, dpi=300, bbox_inches='tight')
    print(f"\n✓ Saved comparison plots:")
    print(f"  {output_pdf}")
    print(f"  {output_png}")
    
    plt.close()


def main():
    """Main execution"""
    print("="*60)
    print("COMPONENT METRICS PERFORMANCE ANALYSIS")
    print("="*60)
    print(f"Platform: {platform.system()} {platform.release()}")
    print(f"Python: {platform.python_version()}")
    print()
    
    script_dir = Path(__file__).parent.parent
    
    # Load timing dataset
    print("Loading timing dataset...")
    timing_data_file = script_dir / 'data' / 'timing_full_dataset.npz'
    print(f"Looking for: {timing_data_file}")
    
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
    
    # Run all three scenarios
    results_dict = {}
    
    print("\n" + "="*60)
    print("COMPREHENSIVE TIMING ANALYSIS")
    print("Testing 3 scenarios for component metrics computation")
    print("="*60)
    
    # Scenario 1: Full metrics with effective graph resistance
    results_dict['full'] = run_timing_test(test_molecules, n_iterations=10, scenario='full')
    
    # Scenario 2: Full metrics without effective graph resistance
    results_dict['no_resistance'] = run_timing_test(test_molecules, n_iterations=10, scenario='no_resistance')
    
    # Scenario 3: Minimal metrics (size only)
    results_dict['minimal'] = run_timing_test(test_molecules, n_iterations=10, scenario='minimal')
    
    # Save results
    output_data = {
        'platform': platform.system(),
        'platform_release': platform.release(),
        'processor': platform.processor(),
        'python_version': platform.python_version(),
        'timestamp': time.strftime('%Y-%m-%d %H:%M:%S'),
        'n_molecules': len(test_molecules),
        'n_iterations': 10,
        'scenarios': {
            'full': {
                'description': 'Full metrics including effective graph resistance (limit=None)',
                'results': results_dict['full']
            },
            'no_resistance': {
                'description': 'Full metrics without effective graph resistance (limit=0)',
                'results': results_dict['no_resistance']
            },
            'minimal': {
                'description': 'Minimal metrics - component size only (compute_component_metrics=False)',
                'results': results_dict['minimal']
            }
        }
    }
    
    output_file = script_dir / 'data' / 'component_metrics_timing_analysis.json'
    with open(output_file, 'w') as f:
        json.dump(output_data, f, indent=2)
    
    print(f"\n✓ Saved timing results to {output_file}")
    
    # Compare scenarios
    print("\nComparing scenarios...")
    comparison_stats = compare_scenarios(results_dict)
    
    if comparison_stats:
        print(f"\nCOMPARISON SUMMARY:")
        print(f"  Mean resistance cost: {comparison_stats['mean_resistance_cost_ms']:.2f} ms ({comparison_stats['mean_resistance_cost_pct']:.1f}%)")
        print(f"  Mean other metrics cost: {comparison_stats['mean_other_metrics_cost_ms']:.2f} ms ({comparison_stats['mean_other_metrics_cost_pct']:.1f}%)")
        print(f"  Speedup (no resistance): {comparison_stats['mean_speedup_no_resistance']:.2f}×")
        print(f"  Speedup (minimal): {comparison_stats['mean_speedup_minimal']:.2f}×")
        
        # Save comparison results
        comparison_file = script_dir / 'data' / 'component_metrics_comparison.json'
        with open(comparison_file, 'w') as f:
            json.dump(comparison_stats, f, indent=2)
        
        print(f"\n✓ Saved comparison results to {comparison_file}")
    
    # Generate comparison plots
    print("\nGenerating comparison plots...")
    plot_comparison(results_dict, comparison_stats)
    
    print("\n" + "="*60)
    print("Component metrics timing analysis complete!")
    print("="*60)


if __name__ == '__main__':
    main()
