# SQL Export Functionality for PFASGroups

This document describes the new `to_sql()` methods added to the `MoleculeResult` and `ResultsModel` classes in PFASGroups.

## Overview

The `to_sql()` methods allow you to export PFAS group detection results directly to SQL databases. Both SQLite and PostgreSQL/MySQL are supported.

## Installation

The SQL export functionality requires additional dependencies:

```bash
pip install sqlalchemy
```

Or install with the optional database dependencies:

```bash
pip install PFASgroups[database]
```

## Methods

### MoleculeResult.to_sql()

Export a single molecule's PFAS group matches to a SQL database.

```python
from PFASgroups import parse_smiles

results = parse_smiles(["C(C(C(C(F)(F)F)(F)F)(F)F)(=O)O"])
molecule_result = results[0]

# Export to SQLite
molecule_result.to_sql(filename='pfas_data.db')

# Export to PostgreSQL
molecule_result.to_sql(
    dbname='pfas_database',
    user='username',
    password='password',
    host='localhost',
    port=5432
)
```

### ResultsModel.to_sql()

Export multiple molecules' data in a single batch operation (more efficient).

```python
from PFASgroups import parse_smiles

smiles_list = [
    "C(C(C(C(F)(F)F)(F)F)(F)F)(=O)O",  # PFBA
    "C(C(C(C(C(C(F)(F)F)(F)F)(F)F)(F)F)(F)F)(=O)O",  # PFHxA
    # ... more SMILES
]

results = parse_smiles(smiles_list)

# Export all to SQLite
results.to_sql(filename='pfas_results.db')

# Export all to PostgreSQL
results.to_sql(
    dbname='pfas_db',
    user='myuser',
    password='mypass',
    host='localhost',
    port=5432
)
```

## Parameters

Both methods accept the same parameters:

### Connection Parameters

**For SQLite:**
- `filename` (str, optional): Path to SQLite database file.

**For PostgreSQL/MySQL:**
- `dbname` (str, optional): Database name.
- `user` (str, optional): Database username. Defaults to `os.environ['DB_USER']` if not provided.
- `password` (str, optional): Database password. Defaults to `os.environ['DB_PASSWORD']` if not provided.
- `host` (str, optional): Database host. Defaults to `os.environ.get('DB_HOST', 'localhost')`.
- `port` (int, optional): Database port. Defaults to `os.environ.get('DB_PORT', 5432)`.

### Table Parameters

- `components_table` (str, default: "components"): Name of the table to store component-level data.
- `groups_table` (str, default: "pfas_groups_in_compound"): Name of the table to store PFAS group matches.

### Behavior Parameters

- `if_exists` (str, default: "append"): How to behave if tables exist:
  - `'fail'`: Raise an error if table exists.
  - `'replace'`: Drop the table and create a new one.
  - `'append'`: Add data to existing table.

## Database Schema

### components Table

Stores detailed component-level information:

| Column | Type | Description |
|--------|------|-------------|
| smiles | TEXT | SMILES representation of the molecule |
| group_id | INTEGER | Numeric ID of the PFAS group |
| group_name | TEXT | Name of the PFAS group |
| smarts_label | TEXT | SMARTS pattern label for the component |
| component_atoms | TEXT | Comma-separated atom indices |

### pfas_groups_in_compound Table

Stores aggregated PFAS group matches per molecule:

| Column | Type | Description |
|--------|------|-------------|
| smiles | TEXT | SMILES representation of the molecule |
| group_id | INTEGER | Numeric ID of the PFAS group |
| group_name | TEXT | Name of the PFAS group |
| match_count | INTEGER | Number of times this group matched |

## Environment Variables

To avoid hardcoding credentials, you can use environment variables:

```bash
export DB_USER="your_username"
export DB_PASSWORD="your_password"
export DB_HOST="localhost"
export DB_PORT="5432"
```

Then call without credentials:

```python
results.to_sql(dbname='pfas_database')
```

## Examples

### Basic SQLite Export

```python
from PFASgroups import parse_smiles

smiles_list = ["C(C(C(C(F)(F)F)(F)F)(F)F)(=O)O"]
results = parse_smiles(smiles_list)

# Create new database or append to existing
results.to_sql(filename='pfas_results.db', if_exists='append')
```

### Custom Table Names

```python
results.to_sql(
    filename='pfas_results.db',
    components_table='my_components',
    groups_table='my_pfas_groups'
)
```

### PostgreSQL with Environment Variables

```python
import os

# Set credentials
os.environ['DB_USER'] = 'pfas_user'
os.environ['DB_PASSWORD'] = 'secure_password'

# Export without explicit credentials
results.to_sql(
    dbname='pfas_database',
    host='db.example.com',
    port=5432
)
```

### Replace Existing Data

```python
# First export
results.to_sql(filename='pfas.db', if_exists='replace')

# Later, replace with new data
new_results.to_sql(filename='pfas.db', if_exists='replace')
```

### Batch Processing

```python
from PFASgroups import parse_smiles

# Process in batches to manage memory
batch_size = 100
all_smiles = [...]  # Large list of SMILES

for i in range(0, len(all_smiles), batch_size):
    batch = all_smiles[i:i+batch_size]
    results = parse_smiles(batch)
    
    # Append to same database
    results.to_sql(
        filename='pfas_large.db',
        if_exists='append' if i > 0 else 'replace'
    )
```

## Querying the Database

Once exported, you can query the data using SQL:

### Get all molecules with a specific PFAS group

```python
import sqlite3

conn = sqlite3.connect('pfas_results.db')
cursor = conn.cursor()

cursor.execute("""
    SELECT DISTINCT smiles, match_count
    FROM pfas_groups_in_compound
    WHERE group_name = 'Perfluoroalkyl carboxylic acid'
""")

for smiles, count in cursor.fetchall():
    print(f"{smiles}: {count} matches")

conn.close()
```

### Get component details for a molecule

```python
cursor.execute("""
    SELECT smarts_label, component_atoms
    FROM components
    WHERE smiles = ?
""", (target_smiles,))

for smarts, atoms in cursor.fetchall():
    print(f"{smarts}: atoms {atoms}")
```

## Error Handling

The methods will raise errors in the following cases:

1. **Missing dependencies**: If pandas or sqlalchemy are not installed.
2. **No connection info**: If neither `filename` nor `dbname` is provided.
3. **Missing credentials**: If using PostgreSQL without credentials (and no environment variables set).
4. **Invalid parameters**: If invalid values are provided for `if_exists` or other parameters.

```python
try:
    results.to_sql(filename='pfas.db')
except ImportError:
    print("Please install sqlalchemy: pip install sqlalchemy")
except ValueError as e:
    print(f"Configuration error: {e}")
```

## Testing

Comprehensive tests are available in `tests/test_results_sql.py`. Run with:

```bash
pytest tests/test_results_sql.py -v
```

The tests cover:
- SQLite and PostgreSQL exports
- Custom table names
- Append, replace, and fail modes
- Large datasets (hundreds of molecules)
- Data integrity verification
- Error handling

## Performance Considerations

1. **Batch operations**: Use `ResultsModel.to_sql()` instead of calling `MoleculeResult.to_sql()` multiple times.
2. **Append mode**: Use `if_exists='append'` for incremental updates rather than replacing the entire database.
3. **Indexing**: Consider adding indexes on frequently queried columns:

```sql
CREATE INDEX idx_smiles ON pfas_groups_in_compound(smiles);
CREATE INDEX idx_group_name ON pfas_groups_in_compound(group_name);
```

## Limitations

1. The methods assume standard pandas/sqlalchemy behavior for SQL operations.
2. MySQL support uses the same connection string format as PostgreSQL (may need adjustment for specific MySQL configurations).
3. Large datasets may require significant memory for the pandas DataFrame conversion.

## See Also

- [PFASGroups Documentation](https://pfasgroups.readthedocs.io)
- [SQLAlchemy Documentation](https://docs.sqlalchemy.org/)
- [Pandas to_sql Documentation](https://pandas.pydata.org/docs/reference/api/pandas.DataFrame.to_sql.html)
