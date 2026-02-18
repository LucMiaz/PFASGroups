import sqlite3

conn = sqlite3.connect('database/pfas_benchmark.db')
cursor = conn.cursor()

# Check ALL highly_branched records
cursor.execute('''
    SELECT m.id, p.success
    FROM molecules m
    LEFT JOIN HalogenGroups_results p ON m.id = p.molecule_id
    WHERE m.dataset_type='highly_branched'
    ORDER BY m.id
''')

rows = cursor.fetchall()
print(f'Total highly_branched molecules: {len(rows)}')
print(f'First 5: {rows[:5]}')
print(f'Last 5: {rows[-5:]}')

# Group by success value
success_counts = {}
for row in rows:
    success_counts[row[1]] = success_counts.get(row[1], 0) + 1

print(f'\nSuccess value distribution: {success_counts}')

conn.close()
