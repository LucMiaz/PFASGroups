import sqlite3
import statistics

conn = sqlite3.connect('review-app/database/pfas_benchmark.db')
c = conn.cursor()

# Get sample times
c.execute('''
    SELECT p.execution_time, a.execution_time 
    FROM molecules m
    LEFT JOIN HalogenGroup_results p ON m.id = p.molecule_id
    LEFT JOIN atlas_results a ON m.id = a.molecule_id
    WHERE p.execution_time IS NOT NULL 
      AND a.execution_time IS NOT NULL
      AND m.dataset_type NOT IN ('complex_branched', 'definitions', 'non_fluorinated', 'highly_branched')
    LIMIT 10
''')

rows = c.fetchall()
print('Sample times (first 10):')
for p, a in rows:
    print(f'HalogenGroup: {p:.6f}s ({p*1000:.2f}ms), Atlas: {a:.6f}s ({a*1000:.2f}ms)')

# Get averages
c.execute('''
    SELECT AVG(p.execution_time), AVG(a.execution_time), COUNT(*)
    FROM molecules m
    LEFT JOIN HalogenGroup_results p ON m.id = p.molecule_id
    LEFT JOIN atlas_results a ON m.id = a.molecule_id
    WHERE p.execution_time IS NOT NULL 
      AND a.execution_time IS NOT NULL
      AND m.dataset_type NOT IN ('complex_branched', 'definitions', 'non_fluorinated', 'highly_branched')
''')

result = c.fetchone()
print(f'\n=== AVERAGES over {result[2]} molecules ===')
print(f'HalogenGroup: {result[0]:.6f}s = {result[0]*1000:.2f}ms')
print(f'Atlas:      {result[1]:.6f}s = {result[1]*1000:.2f}ms')
print(f'Ratio:      {result[0]/result[1]:.2f}x (HalogenGroup/Atlas)')

conn.close()
