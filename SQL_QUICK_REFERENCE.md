# Quick Reference: to_sql() Methods

## Installation
```bash
pip install sqlalchemy
```

## SQLite Usage
```python
from PFASgroups import parse_smiles

# Parse molecules
results = parse_smiles(["C(C(C(C(F)(F)F)(F)F)(F)F)(=O)O"])

# Export to SQLite
results.to_sql(filename='pfas.db')
```

## PostgreSQL Usage
```python
# With explicit credentials
results.to_sql(
    dbname='pfas_db',
    user='username',
    password='password',
    host='localhost',
    port=5432
)

# With environment variables
import os
os.environ['DB_USER'] = 'username'
os.environ['DB_PASSWORD'] = 'password'

results.to_sql(dbname='pfas_db')
```

## Common Parameters
```python
results.to_sql(
    filename='pfas.db',              # SQLite file path
    components_table='components',   # Component details table
    groups_table='pfas_groups',      # Group summary table  
    if_exists='append'               # 'append', 'replace', or 'fail'
)
```

## Query Results
```python
import sqlite3

conn = sqlite3.connect('pfas.db')
cursor = conn.cursor()

# Get all PFAS groups
cursor.execute("SELECT * FROM pfas_groups_in_compound")
for row in cursor.fetchall():
    print(row)

# Get components for a molecule
cursor.execute("""
    SELECT smarts_label, component_atoms 
    FROM components 
    WHERE smiles = ?
""", (smiles,))

conn.close()
```

## Tables Created

### pfas_groups_in_compound
- smiles (TEXT)
- group_id (INTEGER)
- group_name (TEXT)
- match_count (INTEGER)

### components
- smiles (TEXT)
- group_id (INTEGER)
- group_name (TEXT)
- smarts_label (TEXT)
- component_atoms (TEXT)

## See Also
- Full documentation: `SQL_EXPORT_README.md`
- Demo script: `demo_to_sql.py`
- Tests: `tests/test_results_sql.py`
