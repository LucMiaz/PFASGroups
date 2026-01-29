import sqlite3

conn = sqlite3.connect('review-app/database/pfas_benchmark.db')
c = conn.cursor()

datasets = ['timing', 'oecd', 'enhanced']

for ds in datasets:
    c.execute('SELECT COUNT(*) FROM molecules WHERE dataset_type=?', (ds,))
    count = c.fetchone()[0]
    print(f"{ds}: {count} molecules")
    
    # Check if they have timing data
    c.execute('''
        SELECT COUNT(*)
        FROM molecules m
        JOIN pfasgroups_results p ON m.id = p.molecule_id
        JOIN atlas_results a ON m.id = a.molecule_id
        WHERE m.dataset_type=? 
          AND p.execution_time IS NOT NULL 
          AND a.execution_time IS NOT NULL
    ''', (ds,))
    with_timing = c.fetchone()[0]
    print(f"  -> with timing data: {with_timing}")

conn.close()
