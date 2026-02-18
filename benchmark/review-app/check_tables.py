import sqlite3

conn = sqlite3.connect('database/pfas_benchmark.db')
cursor = conn.cursor()

# List all tables
cursor.execute("SELECT name FROM sqlite_master WHERE type='table'")
tables = cursor.fetchall()
print('Available tables:')
for table in tables:
    print(f'  - {table[0]}')

print()

# Check structure of HalogenGroups_results table
cursor.execute("PRAGMA table_info(HalogenGroups_results)")
columns = cursor.fetchall()
print('HalogenGroups_results columns:')
for col in columns:
    print(f'  {col[1]} ({col[2]})')

print()

# Check highly_branched results
cursor.execute('''
    SELECT m.id, m.smiles, p.detected_groups, p.success
    FROM molecules m
    LEFT JOIN HalogenGroups_results p ON m.id = p.molecule_id
    WHERE m.dataset_type='highly_branched'
    LIMIT 10
''')

print('Sample highly_branched molecules with HalogenGroups results:')
for row in cursor.fetchall():
    mol_id, smiles, detected_groups, success = row
    print(f'ID {mol_id}: {smiles[:40]}... | Success: {success} | Groups: {detected_groups}')

# Check Atlas results
cursor.execute('''
    SELECT m.id, m.smiles, a.first_class, a.success
    FROM molecules m
    LEFT JOIN atlas_results a ON m.id = a.molecule_id
    WHERE m.dataset_type='highly_branched'
    LIMIT 10
''')

print('\nSample highly_branched molecules with Atlas results:')
for row in cursor.fetchall():
    mol_id, smiles, first_class, success = row
    print(f'ID {mol_id}: {smiles[:40]}... | Success: {success} | Class: {first_class}')

conn.close()
