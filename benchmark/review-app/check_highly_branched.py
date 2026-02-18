import sqlite3
import json

conn = sqlite3.connect('database/pfas_benchmark.db')
cursor = conn.cursor()

# Check highly_branched molecules
cursor.execute('''
    SELECT id, smiles, HalogenGroups_detected, atlas_detected, 
           HalogenGroups_success, atlas_success 
    FROM molecules 
    WHERE dataset_type='highly_branched' 
    LIMIT 10
''')

print('Sample highly_branched molecules:')
print('-' * 100)
for row in cursor.fetchall():
    mol_id, smiles, pfas_det, atlas_det, pfas_success, atlas_success = row
    pfas_groups = json.loads(pfas_det) if pfas_det else []
    atlas_groups = json.loads(atlas_det) if atlas_det else []
    
    print(f'ID {mol_id}: {smiles[:40]}...')
    print(f'  HalogenGroups: success={pfas_success}, detected={len(pfas_groups)} groups: {pfas_groups[:5]}')
    print(f'  Atlas: success={atlas_success}, detected={atlas_groups}')
    print()

# Check detection statistics
cursor.execute('''
    SELECT 
        COUNT(*) as total,
        SUM(CASE WHEN HalogenGroups_success = 1 THEN 1 ELSE 0 END) as pfas_success,
        SUM(CASE WHEN atlas_success = 1 THEN 1 ELSE 0 END) as atlas_success,
        SUM(CASE WHEN json_array_length(HalogenGroups_detected) > 0 THEN 1 ELSE 0 END) as pfas_detected,
        SUM(CASE WHEN json_array_length(atlas_detected) > 0 THEN 1 ELSE 0 END) as atlas_detected
    FROM molecules 
    WHERE dataset_type='highly_branched'
''')

stats = cursor.fetchone()
print('\nStatistics for highly_branched dataset:')
print(f'Total molecules: {stats[0]}')
print(f'HalogenGroups success: {stats[1]}')
print(f'Atlas success: {stats[2]}')
print(f'HalogenGroups detected groups: {stats[3]}')
print(f'Atlas detected: {stats[4]}')

conn.close()
