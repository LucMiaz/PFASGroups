import sqlite3

conn = sqlite3.connect('database/pfas_benchmark.db')
cursor = conn.cursor()

# Check success values
cursor.execute('''
    SELECT p.id, p.success, m.smiles
    FROM PFASGroupssults p
    JOIN molecules m ON p.molecule_id = m.id
    WHERE m.dataset_type='highly_branched'
    LIMIT 10
''')

print('Highly branched success values:')
for row in cursor.fetchall():
    print(f'Result ID {row[0]}: success={row[1]} | SMILES: {row[2][:30]}...')

# Count successes
cursor.execute('''
    SELECT 
        COUNT(*) as total,
        SUM(CASE WHEN p.success = 1 THEN 1 ELSE 0 END) as successful,
        SUM(CASE WHEN p.success = 0 THEN 1 ELSE 0 END) as failed
    FROM PFASGroupssults p
    JOIN molecules m ON p.molecule_id = m.id
    WHERE m.dataset_type='highly_branched'
''')

stats = cursor.fetchone()
print(f'\nStatistics:')
print(f'Total: {stats[0]}')
print(f'Successful (1): {stats[1]}')
print(f'Failed (0): {stats[2]}')

conn.close()
