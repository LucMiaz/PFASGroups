import sqlite3

conn = sqlite3.connect('database/pfas_benchmark.db')
cursor = conn.cursor()

# Check for orphaned results
cursor.execute('''
    SELECT COUNT(*)
    FROM HalogenGroups_results p
    WHERE NOT EXISTS (SELECT 1 FROM molecules m WHERE m.id = p.molecule_id)
''')

orphaned = cursor.fetchone()[0]
print(f'Orphaned HalogenGroups_results records: {orphaned}')

# Check the new records
cursor.execute('''
    SELECT m.id, p.success
    FROM molecules m
    JOIN HalogenGroups_results p ON m.id = p.molecule_id
    WHERE m.id >= 6123
    ORDER BY m.id
    LIMIT 10
''')

print('\nNew records (ID >= 6123):')
for row in cursor.fetchall():
    print(f'Molecule {row[0]}: success={row[1]}')

conn.close()
