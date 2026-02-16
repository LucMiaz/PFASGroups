# PFASGroups Database Schema

## Overview

The `apply_pfasgroups_all_halogens.py` script now uses the `ResultsModel.to_sql()` method from PFASGroups to save results to a PostgreSQL database with a properly normalized three-table schema.

## Database Tables

### 1. `pfasgroups`
Contains PFAS group definitions. Populated once during initialization.

**Columns:**
- `id` (INTEGER, PRIMARY KEY) - Group ID
- `name` (TEXT) - Group name
- `alias` (TEXT) - Alternative name
- `category` (TEXT) - Functional category
- `pathway_type` (TEXT) - Perfluoroalkyl, polyfluoroalkyl, etc.
- `smarts_primary` (TEXT) - Primary SMARTS pattern
- `smarts_secondary` (TEXT) - Secondary SMARTS pattern
- `created_at` (TIMESTAMP) - Creation timestamp

**Example:**
```sql
SELECT id, name FROM pfasgroups LIMIT 3;
```
| id | name |
|----|------|
| 1  | Perfluoroalkyl carboxylic acids |
| 2  | Polyfluoroalkyl carboxylic acid |
| 3  | Perfluoroalkyl dicarboxylic acids |

### 2. `pfasgroups_in_molecules`
Stores which PFAS groups were detected in which molecules.

**Columns:**
- `id` (SERIAL, PRIMARY KEY) - Auto-increment ID
- `smiles` (TEXT) - Molecule SMILES string
- `molecule_id` (INTEGER) - Foreign key to molecules table (clinventory)
- `group_id` (INTEGER) - Foreign key to pfasgroups.id
- `group_name` (TEXT) - Group name (denormalized for convenience)
- `match_count` (INTEGER) - Number of times this group matched
- `created_at` (TIMESTAMP) - Creation timestamp

**Indexes:**
- `idx_pfasgroups_in_molecules_smiles` on `smiles`
- `idx_pfasgroups_in_molecules_molecule_id` on `molecule_id`
- `idx_pfasgroups_in_molecules_group_id` on `group_id`

**Example:**
```sql
SELECT molecule_id, group_id, group_name, match_count 
FROM pfasgroups_in_molecules LIMIT 3;
```
| molecule_id | group_id | group_name | match_count |
|-------------|----------|------------|-------------|
| 9           | 52       | polyfluoroalkyl | 1 |
| 10          | 52       | polyfluoroalkyl | 1 |
| 19          | 58       | Perfluoroaryl compounds | 1 |

### 3. `components_in_molecules`
Stores detailed component information for each match.

**Columns:**
- `id` (SERIAL, PRIMARY KEY) - Auto-increment ID
- `smiles` (TEXT) - Molecule SMILES string
- `molecule_id` (INTEGER) - Foreign key to molecules table (clinventory)
- `group_id` (INTEGER) - Foreign key to pfasgroups.id
- `group_name` (TEXT) - Group name (denormalized)
- `smarts_label` (TEXT) - SMARTS pattern that matched
- `component_atoms` (TEXT) - Comma-separated atom indices
- `created_at` (TIMESTAMP) - Creation timestamp

**Indexes:**
- `idx_components_in_molecules_smiles` on `smiles`
- `idx_components_in_molecules_molecule_id` on `molecule_id`
- `idx_components_in_molecules_group_id` on `group_id`

**Example:**
```sql
SELECT molecule_id, group_id, group_name, component_atoms 
FROM components_in_molecules LIMIT 3;
```
| molecule_id | group_id | group_name | component_atoms |
|-------------|----------|------------|-----------------|
| 9           | 52       | polyfluoroalkyl | 1 |
| 10          | 52       | polyfluoroalkyl | 1 |
| 19          | 58       | Perfluoroaryl compounds | 4,8,9 |

## Usage

### Initialize and Populate Tables

The script automatically initializes all three tables on first run with `--save-to-db`:

```bash
conda run -n chem python apply_pfasgroups_all_halogens.py \
  --save-to-db \
  --limit 1000 \
  --batch-size 100
```

### Query Examples

**Count molecules with specific PFAS groups:**
```sql
SELECT g.name, COUNT(DISTINCT pm.molecule_id) as molecule_count
FROM pfasgroups_in_molecules pm
JOIN pfasgroups g ON pm.group_id = g.id
GROUP BY g.id, g.name
ORDER BY molecule_count DESC;
```

**Find molecules with multiple PFAS groups:**
```sql
SELECT molecule_id, COUNT(DISTINCT group_id) as num_groups
FROM pfasgroups_in_molecules
GROUP BY molecule_id
HAVING COUNT(DISTINCT group_id) > 1
ORDER BY num_groups DESC;
```

**Get component details for a specific molecule:**
```sql
SELECT c.group_name, c.smarts_label, c.component_atoms
FROM components_in_molecules c
WHERE c.molecule_id = 19;
```

**Join with molecules table:**
```sql
SELECT m.id, m.smiles, m.formula, 
       pm.group_name, pm.match_count
FROM molecules m
JOIN pfasgroups_in_molecules pm ON m.id = pm.molecule_id
WHERE pm.group_id IN (52, 58)
LIMIT 10;
```

## Implementation Details

### ResultsModel.to_sql() Method

The script uses PFASGroups' built-in `ResultsModel.to_sql()` method which:
- Creates DataFrames from parsed results
- Uses pandas `to_sql()` for efficient batch inserts
- Supports both SQLite and PostgreSQL
- Allows custom table names via `components_table` and `groups_table` parameters

### Molecule ID Mapping

After `ResultsModel.to_sql()` inserts data using SMILES as identifier, the script updates `molecule_id` columns:

```python
for smiles, mol_id in compound_id_map.items():
    cur.execute(
        "UPDATE pfasgroups_in_molecules SET molecule_id = %s WHERE smiles = %s AND molecule_id IS NULL",
        (mol_id, smiles)
    )
    cur.execute(
        "UPDATE components_in_molecules SET molecule_id = %s WHERE smiles = %s AND molecule_id IS NULL",
        (mol_id, smiles)
    )
```

This allows linking results back to the original `molecules` table in clinventory.

## Performance

- **Table initialization**: Once per database (checks if pfasgroups table is empty)
- **Batch inserts**: Uses pandas `to_sql()` with `if_exists='append'`
- **Indexes**: Created automatically for efficient queries
- **Throughput**: ~74 molecules/second (median processing time 10.75ms)

## Migration from Old Schema

The old `pfasgroups_results` table (if it exists) can coexist with the new schema. To migrate data or remove the old table:

```sql
-- View old table (if exists)
SELECT COUNT(*) FROM pfasgroups_results;

-- Drop old table (after verifying new tables are populated)
DROP TABLE IF EXISTS pfasgroups_results;
```
