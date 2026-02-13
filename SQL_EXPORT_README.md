# SQL Export Functionality for PFASGroups

This document describes the SQL database integration features in PFASGroups, including exporting results, parsing molecules from databases, and setting up PFAS groups in databases.

## Overview

PFASGroups provides three main database integration features:

1. **Export results to SQL databases** - Write PFAS detection results to databases
2. **Parse molecules from databases** - Directly analyze molecules stored in databases
3. **Setup PFAS groups metadata** - Store PFAS group definitions in databases

Both SQLite and PostgreSQL/MySQL are supported.

## Installation

The SQL features require additional dependencies:

```bash
pip install sqlalchemy
```

Or install with the optional database dependencies:

```bash
pip install PFASgroups[database]
```

## 1. Exporting Results to Databases

### Methods

#### MoleculeResult.to_sql()

Export a single molecule's PFAS group matches to a SQL database.

```python
from PFASgroups import parse_smiles

results = parse_smiles(["C(C(C(C(F)(F)F)(F)F)(F)F)(=O)O"])
molecule_result = results[0]

# Using connection string
molecule_result.to_sql(conn='postgresql://user:pass@localhost/pfas_db')

# Using SQLAlchemy engine
from sqlalchemy import create_engine
engine = create_engine('sqlite:///pfas_data.db')
molecule_result.to_sql(conn=engine)

# Legacy filename parameter (SQLite only)
molecule_result.to_sql(filename='pfas_data.db')
```

#### ResultsModel.to_sql()

Export multiple molecules' data in a single batch operation (more efficient).

```python
from PFASgroups import parse_smiles

smiles_list = [
    "C(C(C(C(F)(F)F)(F)F)(F)F)(=O)O",  # PFBA
    "C(C(C(C(C(C(F)(F)F)(F)F)(F)F)(F)F)(F)F)(=O)O",  # PFHxA
]

results = parse_smiles(smiles_list)

# Using connection string
results.to_sql(conn='postgresql://user:pass@localhost/pfas_db')

# Using SQLAlchemy engine
from sqlalchemy import create_engine
engine = create_engine('sqlite:///pfas_results.db')
results.to_sql(conn=engine)
```

### Connection Parameters

Both methods accept:

- **`conn`** (str or sqlalchemy.engine.Engine): Database connection
  - Connection string: `'postgresql://user:pass@host:port/dbname'`
  - SQLite: `'sqlite:///path/to/file.db'`
  - SQLAlchemy Engine object
- **`filename`** (str, optional): Legacy parameter for SQLite (use `conn` instead)
- **`components_table`** (str, default: "components"): Table for component data
- **`groups_table`** (str, default: "pfas_groups_in_compound"): Table for group matches
- **`if_exists`** (str, default: "append"): 'fail', 'replace', or 'append'

### Database Schema

#### components Table

| Column | Type | Description |
|--------|------|-------------|
| smiles | TEXT | SMILES representation |
| group_id | INTEGER | PFAS group numeric ID |
| group_name | TEXT | PFAS group name |
| smarts_label | TEXT | SMARTS pattern label |
| component_atoms | TEXT | Comma-separated atom indices |

#### pfas_groups_in_compound Table

| Column | Type | Description |
|--------|------|-------------|
| smiles | TEXT | SMILES representation |
| group_id | INTEGER | PFAS group numeric ID |
| group_name | TEXT | PFAS group name |
| match_count | INTEGER | Number of matches |

## 2. Parsing Molecules from Databases

The `parse_from_database()` function reads molecules directly from database tables and analyzes them for PFAS groups.

### Basic Usage

```python
from PFASgroups import parse_from_database

# Parse molecules from a table
results = parse_from_database(
    conn='postgresql://user:pass@localhost/chem_db',
    table='molecules',
    mol_column='rdkit_mol',
    smiles_column='canonical_smiles'
)

# Parse with custom SQL query
results = parse_from_database(
    conn='postgresql://user:pass@localhost/chem_db',
    query="SELECT id, mol, smiles FROM compounds WHERE mw < 1000",
    batch_size=500
)

# Parse and write results back
results = parse_from_database(
    conn='sqlite:///chemicals.db',
    table='test_compounds',
    write_results=True,
    components_table='pfas_components',
    groups_table='pfas_matches'
)
```

### Parameters

- **`conn`** (str or Engine): Database connection
- **`query`** (str, optional): Custom SQL query
- **`table`** (str, optional): Table name (if query not provided)
- **`mol_column`** (str, default: 'mol'): Column with RDKit mol objects
- **`smiles_column`** (str, default: 'smiles'): Fallback SMILES column
- **`inchi_column`** (str, default: 'inchi'): Fallback InChI column
- **`id_column`** (str, default: 'id'): Unique identifier column
- **`batch_size`** (int, default: 1000): Molecules per batch
- **`write_results`** (bool, default: True): Write to database
- **`output_table`** (str, optional): Summary results table
- **`components_table`** (str): Component details table
- **`groups_table`** (str): Group matches table

### Molecule Column Fallback Order

The function attempts to parse molecules in this order:

1. **mol_column**: Binary RDKit mol or mol block
2. **smiles_column**: SMILES string
3. **inchi_column**: InChI string

### Complete Example

```python
from sqlalchemy import create_engine
from PFASgroups import parse_from_database, setup_pfas_groups_database

# 1. Setup database with PFAS groups metadata
engine = create_engine('postgresql://user:pass@localhost/pfas_db')
setup_pfas_groups_database(conn=engine)

# 2. Parse molecules from existing table
results = parse_from_database(
    conn=engine,
    table='chemical_inventory',
    mol_column='molecule',
    smiles_column='smiles',
    id_column='compound_id',
    batch_size=1000,
    write_results=True
)

# 3. Query results
import pandas as pd
df_groups = pd.read_sql('SELECT * FROM pfas_groups_in_compound', engine)
print(f"Found {len(df_groups)} PFAS group matches")
```

## 3. Setting Up PFAS Groups in Databases

The `setup_pfas_groups_database()` function creates tables with PFAS group definitions and SMARTS patterns.

### Basic Usage

```python
from PFASgroups import setup_pfas_groups_database
from sqlalchemy import create_engine

# PostgreSQL
engine = create_engine('postgresql://user:pass@localhost/pfas_db')
stats = setup_pfas_groups_database(
    conn=engine,
    groups_info_table='pfas_groups',
    smarts_table='pfas_smarts_patterns'
)

# SQLite
stats = setup_pfas_groups_database(
    conn='sqlite:///pfas_metadata.db',
    groups_info_table='pfas_groups_info',
    smarts_table='pfas_smarts'
)

print(f"Loaded {stats['total_groups']} PFAS groups")
print(f"  - {stats['compute_groups']} compute groups")
print(f"  - {stats['aggregate_groups']} aggregate groups")
print(f"  - {stats['total_smarts']} unique SMARTS patterns")
```

### Created Tables

#### pfas_groups_info Table

| Column | Type | Description |
|--------|------|-------------|
| id | INTEGER | Group numeric ID |
| name | TEXT | Group name |
| compute | BOOLEAN | Whether to compute or aggregate |
| componentSmarts | TEXT | Component type |
| pathType | TEXT | Path finding type |
| constraints | TEXT | JSON constraints |
| smarts_patterns | TEXT | JSON SMARTS patterns |

#### pfas_smarts Table

| Column | Type | Description |
|--------|------|-------------|
| pattern_id | INTEGER | Unique pattern ID |
| smarts | TEXT | SMARTS pattern string |

## Complete Workflow Example

Here's a complete workflow showing how to set up, populate, and query a PFAS database:

```python
from sqlalchemy import create_engine
from PFASgroups import (
    parse_smiles,
    parse_from_database,
    setup_pfas_groups_database
)
import pandas as pd

# Step 1: Create database and setup PFAS groups metadata
engine = create_engine('postgresql://user:pass@localhost/pfas_analysis')

print("Setting up PFAS groups metadata...")
stats = setup_pfas_groups_database(conn=engine)

# Step 2: Parse some initial molecules and store results
print("\nParsing initial molecule set...")
smiles_list = [
    "C(C(C(C(F)(F)F)(F)F)(F)F)(=O)O",  # PFBA
    "C(C(C(C(C(C(F)(F)F)(F)F)(F)F)(F)F)(F)F)(=O)O",  # PFHxA
    "C(C(C(C(C(C(C(C(F)(F)F)(F)F)(F)F)(F)F)(F)F)(F)F)(F)F)(=O)O",  # PFOA
]

results = parse_smiles(smiles_list)
results.to_sql(conn=engine, if_exists='replace')

# Step 3: Add molecules from existing database
print("\nParsing molecules from existing table...")
results_db = parse_from_database(
    conn=engine,
    query="SELECT id, mol, smiles FROM chemical_inventory WHERE suspect_pfas = true",
    batch_size=500,
    write_results=True
)

# Step 4: Query and analyze results
print("\nAnalyzing results...")

# Get all PFAS-containing compounds
df_pfas = pd.read_sql("""
    SELECT DISTINCT smiles, COUNT(DISTINCT group_name) as num_groups
    FROM pfas_groups_in_compound
    GROUP BY smiles
    ORDER BY num_groups DESC
""", engine)

print(f"\nFound {len(df_pfas)} PFAS compounds")
print(f"Top compounds by number of PFAS groups:")
print(df_pfas.head(10))

# Get most common PFAS groups
df_groups = pd.read_sql("""
    SELECT group_name, COUNT(*) as compound_count, SUM(match_count) as total_matches
    FROM pfas_groups_in_compound
    GROUP BY group_name
    ORDER BY compound_count DESC
""", engine)

print(f"\nMost common PFAS groups:")
print(df_groups.head(10))

# Get detailed component information
df_components = pd.read_sql("""
    SELECT c.smiles, c.group_name, c.smarts_label, COUNT(*) as component_count
    FROM components c
    GROUP BY c.smiles, c.group_name, c.smarts_label
    ORDER BY component_count DESC
""", engine)

print(f"\nComponent types found:")
print(df_components.head(10))
```

## Advanced Usage

### Custom Connection with Environment Variables

```python
from PFASgroups import parse_smiles

smiles_list = ["C(C(C(C(F)(F)F)(F)F)(F)F)(=O)O"]
results = parse_smiles(smiles_list)

# Create new database or append to existing
results.to_sql(filename='pfas_results.db', if_exists='append')
```

## Advanced Usage

### Custom Connection with Environment Variables

```python
import os
from sqlalchemy import create_engine

# Set credentials as environment variables
os.environ['DB_USER'] = 'pfas_user'
os.environ['DB_PASSWORD'] = 'secure_password'

# Use environment variables in connection string
user = os.getenv('DB_USER')
password = os.getenv('DB_PASSWORD')
conn_str = f'postgresql://{user}:{password}@localhost/pfas_db'

results.to_sql(conn=conn_str)
```

### Batch Processing Large Datasets

```python
from PFASgroups import parse_from_database
from sqlalchemy import create_engine

engine = create_engine('postgresql://user:pass@localhost/large_db')

# Process in batches to manage memory
results = parse_from_database(
    conn=engine,
    query="SELECT id, mol, smiles FROM compounds WHERE processed = false",
    batch_size=500,  # Adjust based on available memory
    write_results=True
)

print(f"Processed {len(results)} compounds")
```

### Integration with Existing Chemical Databases

```python
from sqlalchemy import create_engine, text
from PFASgroups import parse_from_database
import pandas as pd

# Connect to existing database
engine = create_engine('postgresql://user:pass@localhost/chem_inventory')

# Parse molecules with custom filtering
results = parse_from_database(
    conn=engine,
    query="""
        SELECT 
            c.id,
            c.mol_object as mol,
            c.canonical_smiles as smiles,
            c.inchi
        FROM chemicals c
        LEFT JOIN pfas_analysis p ON c.id = p.chemical_id
        WHERE p.chemical_id IS NULL  -- Only unanalyzed compounds
        AND c.contains_fluorine = true
        LIMIT 10000
    """,
    mol_column='mol',
    smiles_column='smiles',
    inchi_column='inchi',
    id_column='id',
    batch_size=1000,
    write_results=True
)

# Update original table with PFAS flag
with engine.connect() as conn:
    # Get all compounds with PFAS groups
    pfas_compounds = pd.read_sql(
        "SELECT DISTINCT smiles FROM pfas_groups_in_compound",
        engine
    )
    
    # Update original table
    for smiles in pfas_compounds['smiles']:
        conn.execute(text("""
            UPDATE chemicals 
            SET is_pfas = true, 
                last_analyzed = CURRENT_TIMESTAMP
            WHERE canonical_smiles = :smiles
        """), {"smiles": smiles})
    
    conn.commit()
```

### Creating Indexed Views for Fast Queries

```python
from sqlalchemy import create_engine, text

engine = create_engine('postgresql://user:pass@localhost/pfas_db')

with engine.connect() as conn:
    # Create indexed view for fast lookups
    conn.execute(text("""
        CREATE INDEX IF NOT EXISTS idx_groups_smiles 
        ON pfas_groups_in_compound(smiles);
        
        CREATE INDEX IF NOT EXISTS idx_groups_name 
        ON pfas_groups_in_compound(group_name);
        
        CREATE INDEX IF NOT EXISTS idx_components_smiles 
        ON components(smiles);
        
        -- Create materialized view for summary statistics
        CREATE MATERIALIZED VIEW IF NOT EXISTS pfas_summary AS
        SELECT 
            g.smiles,
            COUNT(DISTINCT g.group_name) as num_pfas_groups,
            SUM(g.match_count) as total_matches,
            STRING_AGG(DISTINCT g.group_name, ', ') as pfas_groups_list,
            COUNT(c.smarts_label) as num_components
        FROM pfas_groups_in_compound g
        LEFT JOIN components c ON g.smiles = c.smiles AND g.group_id = c.group_id
        GROUP BY g.smiles;
        
        CREATE INDEX ON pfas_summary(smiles);
        CREATE INDEX ON pfas_summary(num_pfas_groups);
    """))
    
    conn.commit()
    
print("✅ Created indexes and summary view")
```

## Error Handling

```python
from PFASgroups import parse_from_database
from sqlalchemy.exc import SQLAlchemyError

try:
    results = parse_from_database(
        conn='postgresql://user:pass@localhost/pfas_db',
        table='compounds',
        write_results=True
    )
except SQLAlchemyError as e:
    print(f"Database error: {e}")
except ValueError as e:
    print(f"Configuration error: {e}")
except ImportError as e:
    print(f"Missing dependencies: {e}")
    print("Install with: pip install sqlalchemy pandas")
```

## Performance Tips

1. **Use batch processing**: Set appropriate `batch_size` based on available memory
2. **Create indexes**: Index frequently queried columns (smiles, group_name)
3. **Use connection pooling**: Reuse SQLAlchemy engines
4. **Write in batches**: Use `ResultsModel.to_sql()` instead of individual writes
5. **Use materialized views**: For complex aggregate queries
6. **Partition large tables**: For databases with millions of compounds

## Comparison with ZeroPM Database Approach

The PFASGroups database integration is inspired by the ZeroPM database module but simplified:

```python
# ZeroPM approach (Django ORM)
from elements.models import PFASGroup
from elements.scripts.load import load_pfas_groups

load_pfas_groups()  # Syncs to Django models

# PFASGroups approach (Direct SQL)
from PFASgroups import setup_pfas_groups_database

setup_pfas_groups_database(conn='postgresql://...')  # Direct to database
```

### Key Differences:

| Feature | ZeroPM | PFASGroups |
|---------|--------|------------|
| Framework | Django ORM | Direct SQL (SQLAlchemy) |
| Setup | Django models required | Standalone, no framework |
| Flexibility | Django-specific | Any SQL database |
| Complexity | Higher (Django setup) | Lower (direct connection) |
| Use Case | Web applications | Data analysis pipelines |

## Querying Results

### Get All PFAS-Positive Compounds

```python
import pandas as pd
from sqlalchemy import create_engine

engine = create_engine('sqlite:///pfas_results.db')

df = pd.read_sql("""
    SELECT DISTINCT smiles, COUNT(DISTINCT group_name) as num_groups
    FROM pfas_groups_in_compound
    GROUP BY smiles
    HAVING COUNT(DISTINCT group_name) > 0
    ORDER BY num_groups DESC
""", engine)

print(df)
```

### Get Component Details for Specific Molecule

```python
target_smiles = "C(C(C(C(F)(F)F)(F)F)(F)F)(=O)O"

df_components = pd.read_sql("""
    SELECT smarts_label, component_atoms, group_name
    FROM components
    WHERE smiles = :smiles
    ORDER BY group_name, smarts_label
""", engine, params={"smiles": target_smiles})

print(df_components)
```

### Find Molecules with Specific PFAS Groups

```python
df_matches = pd.read_sql("""
    SELECT smiles, match_count
    FROM pfas_groups_in_compound
    WHERE group_name = 'Perfluoroalkyl carboxylic acid'
    ORDER BY match_count DESC
""", engine)

print(df_matches)
```

## Troubleshooting

### Issue: "pandas and sqlalchemy are required"

```bash
pip install pandas sqlalchemy
```

### Issue: "Either 'conn' or 'filename' must be provided"

Make sure to pass either `conn` or `filename`:

```python
# Correct
results.to_sql(conn='sqlite:///pfas.db')
results.to_sql(filename='pfas.db')  # Legacy

# Wrong
results.to_sql()  # Missing connection
```

### Issue: Connection refused to PostgreSQL

Check connection string format and credentials:

```python
# Correct format
conn = 'postgresql://username:password@hostname:5432/database_name'

# Test connection
from sqlalchemy import create_engine
engine = create_engine(conn)
with engine.connect() as connection:
    print("✅ Connected successfully")
```

### Issue: Molecules not parsing from database

Check column names and data types:

```python
import pandas as pd

# Inspect table structure
df = pd.read_sql("SELECT * FROM molecules LIMIT 5", engine)
print(df.columns)
print(df.dtypes)

# Try parsing with correct column names
results = parse_from_database(
    conn=engine,
    table='molecules',
    mol_column='actual_mol_column_name',  # Adjust this
    smiles_column='actual_smiles_column',  # And this
)
```

## See Also

- [PFASGroups Documentation](https://pfasgroups.readthedocs.io)
- [SQLAlchemy Documentation](https://docs.sqlalchemy.org/)
- [Pandas to_sql Documentation](https://pandas.pydata.org/docs/reference/api/pandas.DataFrame.to_sql.html)
- [ZeroPM Database Module](https://github.com/lucmiaz/zeropmdb)

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
