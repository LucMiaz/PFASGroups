import sqlite3
import statistics
from collections import defaultdict

conn = sqlite3.connect('review-app/database/pfas_benchmark.db')
c = conn.cursor()

# Get all timing data by dataset
c.execute('''
    SELECT 
        m.dataset_type,
        p.execution_time as pfas_time,
        a.execution_time as atlas_time,
        m.num_atoms
    FROM molecules m
    LEFT JOIN pfasgroups_results p ON m.id = p.molecule_id
    LEFT JOIN atlas_results a ON m.id = a.molecule_id
    WHERE p.execution_time IS NOT NULL 
      AND a.execution_time IS NOT NULL
      AND m.dataset_type NOT IN ('complex_branched', 'definitions', 'non_fluorinated', 'highly_branched')
''')

rows = c.fetchall()

# Organize by dataset
by_dataset = defaultdict(lambda: {'pfas': [], 'atlas': [], 'atoms': []})
for dataset, pfas, atlas, atoms in rows:
    by_dataset[dataset]['pfas'].append(pfas)
    by_dataset[dataset]['atlas'].append(atlas)
    if atoms:
        by_dataset[dataset]['atoms'].append(atoms)

# Organize by size
by_size = {'small': {'pfas': [], 'atlas': []}, 'medium': {'pfas': [], 'atlas': []}, 'large': {'pfas': [], 'atlas': []}}
for dataset, pfas, atlas, atoms in rows:
    if atoms:
        if atoms < 20:
            by_size['small']['pfas'].append(pfas)
            by_size['small']['atlas'].append(atlas)
        elif atoms <= 40:
            by_size['medium']['pfas'].append(pfas)
            by_size['medium']['atlas'].append(atlas)
        else:
            by_size['large']['pfas'].append(pfas)
            by_size['large']['atlas'].append(atlas)

print("=" * 60)
print("BY DATASET")
print("=" * 60)
for dataset in sorted(by_dataset.keys()):
    d = by_dataset[dataset]
    pfas_mean = statistics.mean(d['pfas']) * 1000
    pfas_median = statistics.median(d['pfas']) * 1000
    atlas_mean = statistics.mean(d['atlas']) * 1000
    ratio = atlas_mean / pfas_mean if pfas_mean > 0 else 0
    count = len(d['pfas'])
    print(f"\n{dataset.upper()}:")
    print(f"  Count: {count}")
    print(f"  PFASgroups Mean: {pfas_mean:.1f}ms, Median: {pfas_median:.1f}ms")
    print(f"  Atlas Mean: {atlas_mean:.1f}ms")
    print(f"  Ratio: {ratio:.2f}x")

print("\n" + "=" * 60)
print("BY SIZE")
print("=" * 60)
for size in ['small', 'medium', 'large']:
    d = by_size[size]
    if d['pfas']:
        pfas_mean = statistics.mean(d['pfas']) * 1000
        pfas_median = statistics.median(d['pfas']) * 1000
        atlas_mean = statistics.mean(d['atlas']) * 1000
        ratio = atlas_mean / pfas_mean if pfas_mean > 0 else 0
        count = len(d['pfas'])
        print(f"\n{size.upper()}:")
        print(f"  Count: {count}")
        print(f"  PFASgroups Mean: {pfas_mean:.1f}ms, Median: {pfas_median:.1f}ms")
        print(f"  Atlas Mean: {atlas_mean:.1f}ms")
        print(f"  Ratio: {ratio:.2f}x")

conn.close()
