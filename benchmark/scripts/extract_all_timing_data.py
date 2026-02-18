"""
Extract ALL timing data from HTML file - complete dataset
"""
import numpy as np
import re
import glob
from pathlib import Path
import sys

# Find the latest timing analysis HTML file
html_files = sorted(glob.glob('review-app/analysis_reports/timing_analysis_*.html'), 
                   key=lambda x: Path(x).stat().st_mtime, reverse=True)

if not html_files:
    print("Error: No timing analysis HTML files found in review-app/analysis_reports/")
    sys.exit(1)

html_file = html_files[0]
print(f"Loading latest timing analysis: {html_file}")

# Read the HTML file
with open(html_file, 'r', encoding='utf-8') as f:
    content = f.read()

# Find the x and y arrays for HalogenGroups (first occurrence)
# The x array is the atom counts
x_match = re.search(r'"x":\[([0-9,]+)\]', content)
if x_match:
    x_str = x_match.group(1)
    atom_counts = np.array([int(x) for x in x_str.split(',')])
    print(f"Found {len(atom_counts)} atom count values")
    print(f"Atom count range: {atom_counts.min()}-{atom_counts.max()}")
    print(f"Unique atom counts: {len(np.unique(atom_counts))}")

# The y array is immediately after the x array in the same object
y_match = re.search(r'"x":\[[0-9,]+\],"y":\[([0-9.,]+)\]', content)
if y_match:
    y_str = y_match.group(1)
    exec_times = np.array([float(y) for y in y_str.split(',')])
    print(f"Found {len(exec_times)} execution time values")
    print(f"Time range: {exec_times.min():.1f}-{exec_times.max():.1f} ms")

# Verify we have matching data
if len(atom_counts) == len(exec_times):
    print(f"\n✓ Data validated: {len(atom_counts)} complete measurements")
    
    # Statistics by size range
    small_mask = atom_counts < 25
    medium_mask = (atom_counts >= 25) & (atom_counts < 50)
    medium_large_mask = (atom_counts >= 50) & (atom_counts < 100)
    large_mask = atom_counts >= 100
    
    print(f"\nSize distribution:")
    print(f"  Small (<25 atoms): {small_mask.sum()} measurements")
    print(f"  Medium (25-50 atoms): {medium_mask.sum()} measurements")
    print(f"  Medium-Large (50-100 atoms): {medium_large_mask.sum()} measurements")
    print(f"  Large (≥100 atoms): {large_mask.sum()} measurements")
    
    print(f"\nMean execution times by size:")
    print(f"  Small: {exec_times[small_mask].mean():.1f} ms")
    print(f"  Medium: {exec_times[medium_mask].mean():.1f} ms")
    print(f"  Medium-Large: {exec_times[medium_large_mask].mean():.1f} ms")
    print(f"  Large: {exec_times[large_mask].mean():.1f} ms")
    
    # Save the full dataset
    np.savez('timing_full_dataset.npz', 
             atom_counts=atom_counts, 
             exec_times=exec_times)
    print(f"\n✓ Full dataset saved to timing_full_dataset.npz")
    
    # Create CSV for easy inspection
    with open('timing_full_dataset.csv', 'w') as f:
        f.write('atom_count,exec_time_ms\n')
        for atoms, time in zip(atom_counts, exec_times):
            f.write(f'{atoms},{time:.2f}\n')
    print(f"✓ CSV saved to timing_full_dataset.csv")
