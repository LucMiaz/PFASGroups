#!/usr/bin/env python3
"""
Simple Windows timing test - run a quick test and compare with Linux

Author: HalogenGroups Benchmark Team
Date: 2026-01-29
"""

import time
import platform
import numpy as np
from pathlib import Path
import sys
import json

print("Starting Windows timing test...")
print(f"Python: {sys.version}")

# Add HalogenGroups to path
pfas_path = r'c:\Users\luc\git\HalogenGroups'
sys.path.insert(0, pfas_path)
print(f"Added to path: {pfas_path}")

try:
    from HalogenGroups import parse_smiles
    print("✓ Successfully imported HalogenGroups")
except Exception as e:
    print(f"✗ Failed to import HalogenGroups: {e}")
    import traceback
    traceback.print_exc()
    sys.exit(1)


def test_molecule_timing(smiles, n_iterations=10):
    """Test timing for a single molecule"""
    times = []
    
    print(f"\nTesting: {smiles}")
    
    for i in range(n_iterations):
        start = time.perf_counter()
        try:
            result = parse_smiles(smiles)
            elapsed_ms = (time.perf_counter() - start) * 1000
            times.append(elapsed_ms)
            print(f"  Iteration {i+1}: {elapsed_ms:.2f} ms")
        except Exception as e:
            print(f"  ✗ Error on iteration {i+1}: {e}")
    
    if times:
        mean_time = np.mean(times)
        std_time = np.std(times, ddof=1)
        print(f"\n  Result: {mean_time:.2f} ± {std_time:.2f} ms (n={len(times)})")
        return mean_time, std_time
    else:
        return None, None


def main():
    """Run quick timing test"""
    print("="*60)
    print(f"Windows Timing Test")
    print(f"Platform: {platform.system()} {platform.release()}")
    print(f"Processor: {platform.processor()}")
    print(f"Python: {platform.python_version()}")
    print("="*60)
    
    # Load timing dataset to get test molecules
    data_file = Path(__file__).parent.parent / 'data' / 'timing_full_dataset.npz'
    
    if not data_file.exists():
        print(f"\n✗ Dataset not found: {data_file}")
        return
    
    print(f"\n✓ Loading dataset from {data_file}")
    
    data = np.load(data_file, allow_pickle=True)
    atom_counts = data['atom_counts']
    smiles_list = data['smiles_list']
    linux_times = data['times']
    
    print(f"✓ Loaded {len(atom_counts)} timing measurements")
    print(f"  Atom range: {atom_counts.min()}-{atom_counts.max()}")
    print(f"  Linux mean time: {linux_times.mean():.2f} ms")
    
    # Test a few representative molecules
    test_cases = [
        (50, "Small molecule (~50 atoms)"),
        (97, "Before breakpoint (97 atoms)"),
        (100, "At breakpoint (100 atoms)"),
        (130, "After breakpoint (130 atoms)")
    ]
    
    results = []
    
    for target_atoms, description in test_cases:
        print(f"\n{'='*60}")
        print(f"{description}")
        print(f"{'='*60}")
        
        # Find a molecule with this atom count
        mask = atom_counts == target_atoms
        if not np.any(mask):
            print(f"  ✗ No molecule with {target_atoms} atoms found")
            continue
        
        indices = np.where(mask)[0]
        idx = indices[0]
        smiles = smiles_list[idx]
        n_atoms = atom_counts[idx]
        
        # Get Linux times for this molecule
        linux_times_mol = linux_times[mask]
        linux_mean = linux_times_mol.mean()
        linux_std = linux_times_mol.std(ddof=1)
        
        print(f"  Molecule: {smiles[:60]}...")
        print(f"  Atoms: {n_atoms}")
        print(f"  Linux: {linux_mean:.2f} ± {linux_std:.2f} ms")
        
        # Test on Windows
        win_mean, win_std = test_molecule_timing(smiles, n_iterations=10)
        
        if win_mean is not None:
            speedup = win_mean / linux_mean
            print(f"\n  Windows: {win_mean:.2f} ± {win_std:.2f} ms")
            print(f"  Speedup (Win/Linux): {speedup:.2f}×")
            
            if speedup > 1.1:
                print(f"  → Linux is {(speedup-1)*100:.0f}% faster")
            elif speedup < 0.9:
                print(f"  → Windows is {(1/speedup-1)*100:.0f}% faster")
            else:
                print(f"  → Similar performance")
            
            results.append({
                'atoms': n_atoms,
                'description': description,
                'smiles': smiles,
                'windows_ms': win_mean,
                'windows_std_ms': win_std,
                'linux_ms': linux_mean,
                'linux_std_ms': linux_std,
                'speedup': speedup
            })
    
    # Save results
    if results:
        output_file = Path(__file__).parent.parent / 'data' / 'quick_windows_test.json'
        output_data = {
            'platform': platform.system(),
            'platform_release': platform.release(),
            'processor': platform.processor(),
            'python_version': platform.python_version(),
            'timestamp': time.strftime('%Y-%m-%d %H:%M:%S'),
            'results': results
        }
        
        with open(output_file, 'w') as f:
            json.dump(output_data, f, indent=2)
        
        print(f"\n{'='*60}")
        print(f"✓ Saved results to {output_file}")
        print(f"{'='*60}")
        
        # Summary
        speedups = [r['speedup'] for r in results]
        mean_speedup = np.mean(speedups)
        
        print(f"\nOverall Summary:")
        print(f"  Tested {len(results)} molecules")
        print(f"  Mean speedup (Win/Linux): {mean_speedup:.2f}×")
        
        if mean_speedup > 1.1:
            print(f"  → Linux is ~{(mean_speedup-1)*100:.0f}% faster on average")
        elif mean_speedup < 0.9:
            print(f"  → Windows is ~{(1/mean_speedup-1)*100:.0f}% faster on average")
        else:
            print(f"  → Both platforms have similar performance")


if __name__ == '__main__':
    main()
