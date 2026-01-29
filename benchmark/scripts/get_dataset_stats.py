import sqlite3
import statistics

conn = sqlite3.connect('review-app/database/pfas_benchmark.db')
c = conn.cursor()

# Get timing data for each of the three main datasets
# Note: Check actual case in database
datasets = ['timing', 'oecd', 'enhanced']

print("Checking datasets...")
for dataset in datasets:
    c.execute('''
        SELECT 
            p.execution_time as pfas_time,
            a.execution_time as atlas_time
        FROM molecules m
        LEFT JOIN pfasgroups_results p ON m.id = p.molecule_id
        LEFT JOIN atlas_results a ON m.id = a.molecule_id
        WHERE p.execution_time IS NOT NULL 
          AND a.execution_time IS NOT NULL
          AND m.dataset_type = ?
    ''', (dataset,))
    
    rows = c.fetchall()
    
    if rows:
        pfas_times = [r[0] for r in rows]
        atlas_times = [r[1] for r in rows]
        
        pfas_mean = statistics.mean(pfas_times) * 1000
        pfas_median = statistics.median(pfas_times) * 1000
        pfas_stdev = statistics.stdev(pfas_times) * 1000 if len(pfas_times) > 1 else 0
        pfas_min = min(pfas_times) * 1000
        pfas_max = max(pfas_times) * 1000
        
        atlas_mean = statistics.mean(atlas_times) * 1000
        atlas_median = statistics.median(atlas_times) * 1000
        atlas_stdev = statistics.stdev(atlas_times) * 1000 if len(atlas_times) > 1 else 0
        
        ratio = pfas_mean / atlas_mean if atlas_mean > 0 else 0
        
        print(f"\n{'='*60}")
        print(f"{dataset.upper()} DATASET")
        print(f"{'='*60}")
        print(f"Count: {len(rows)}")
        print(f"\nPFASgroups:")
        print(f"  Mean:   {pfas_mean:.1f}ms")
        print(f"  Median: {pfas_median:.1f}ms")
        print(f"  StdDev: {pfas_stdev:.1f}ms")
        print(f"  Range:  {pfas_min:.1f}-{pfas_max:.1f}ms")
        print(f"\nPFAS-Atlas:")
        print(f"  Mean:   {atlas_mean:.1f}ms")
        print(f"  Median: {atlas_median:.1f}ms")
        print(f"  StdDev: {atlas_stdev:.1f}ms")
        print(f"\nPerformance Ratio: {ratio:.2f}x (PFASgroups/Atlas)")
        print(f"Speed difference: {(ratio-1)*100:+.1f}%")

conn.close()
