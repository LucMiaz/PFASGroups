# PFASGroups All Halogens Processing Script

This script processes all halogenated compounds (F, Cl, Br, I) from the `clinventory` database using PFASGroups parsing and can save results to both JSON files and a PostgreSQL database table.

## Overview

The script:
1. Queries the `clinventory` database for compounds containing specified halogens
2. Optionally substitutes Cl/Br/I with F to leverage fluorine-based SMARTS patterns
3. Parses each compound using PFASGroups `parse_mol` with `bycomponent=True`
4. Saves results to:
   - **Database table** (`pfasgroups_results`) with `--save-to-db` flag
   - **JSON file** with `--output-json` flag
5. Tracks timing, substitution details, and PFAS group information

## Database Configuration

### Default Settings
- Database: `clinventory`
- Table: `elements_compound`
- Columns: `id`, `smiles`, `formula`
- Connection: `localhost:5436` as user `django`

### Configurable via Arguments
```bash
python apply_pfasgroups_all_halogens.py \
  --db-name clinventory \
  --db-user django \
  --db-password $DJANGO_USER \
  --db-host localhost \
  --db-port 5436 \
  --table-name elements_compound \
  --id-column id \
  --smiles-column smiles \
  --formula-column formula
```

### Alternative Table Structures

For SPIN database tables:
```bash
python apply_pfasgroups_all_halogens.py \
  --table-name spindb_spinstof \
  --id-column id \
  --smiles-column smiles \
  --formula-column bruttoformel
```

## Usage Examples

### Basic Usage (F/Cl/Br/I compounds)
```bash
python apply_pfasgroups_all_halogens.py \
  --output-json results.json
```

### Fluorine Only with Limit
```bash
python apply_pfasgroups_all_halogens.py \
  --halogens F \
  --limit 1000 \
  --output-json fluorine_compounds.json
```

### Chlorine/Bromine with Substitution
Process Cl/Br compounds by temporarily replacing them with F for parsing:
```bash
python apply_pfasgroups_all_halogens.py \
  --halogens Cl,Br \
  --substitute-halogens \
  --output-json chlorine_bromine_as_f.json
```

### Paginated Processing
```bash
python apply_pfasgroups_all_halogens.py \
  --offset 0 \
  --limit 5000 \
  --output-json batch_1.json

python apply_pfasgroups_all_halogens.py \
  --offset 5000 \
  --limit 5000 \
  --output-json batch_2.json
```

### Save to Database
```bash
# Save results to database table
python apply_pfasgroups_all_halogens.py \
  --halogens F \
  --save-to-db

# Save to both database and JSON
python apply_pfasgroups_all_halogens.py \
  --halogens F,Cl,Br,I \
  --save-to-db \
  --output-json results.json

# Use custom results table name
python apply_pfasgroups_all_halogens.py \
  --save-to-db \
  --results-table my_pfas_results
```

## Command-Line Arguments

| Argument | Default | Description |
|----------|---------|-------------|
| `--db-name` | `clinventory` | Database name |
| `--db-user` | `django` | Database user |
| `--db-password` | `$DJANGO_USER` | Database password |
| `--db-host` | `localhost` | Database host |
| `--db-port` | `5432` | Database port |
| `--table-name` | `molecules` | Table to query |
| `--id-column` | `id` | ID column name |
| `--smiles-column` | `smiles` | SMILES column name |
| `--formula-column` | `formula` | Formula column name |
| `--pfasgroups-path` | `.` (current dir) | Path to PFASGroups repo |
| `--run-id` | Auto-generated timestamp | Run identifier |
| `--limit` | None | Max compounds to process |
| `--offset` | `0` | Skip first N compounds |
| `--halogens` | `F,Cl,Br,I` | Comma-separated halogen list |
| `--substitute-halogens` | `False` | Replace Cl/Br/I with F |
| `--output-json` | None | JSON output file path |
| `--save-to-db` | `False` | Save results to database |
| `--results-table` | `pfasgroups_results` | Database table for results |
| `--use-halogen-columns` | `False` | Use numeric halogen columns |
| `--batch-size` | `100` | Progress reporting interval |
| `--run-id` | Auto-generated timestamp | Run identifier |
| `--limit` | None | Max compounds to process |
| `--offset` | `0` | Skip first N compounds |
| `--halogens` | `F,Cl,Br,I` | Comma-separated halogen list |
| `--substitute-halogens` | `False` | Replace Cl/Br/I with F |
| `--output-json` | None | JSON output file path |
| `--batch-size` | `100` | Progress reporting interval |

## Halogen Substitution Feature

When `--substitute-halogens` is enabled:
1. All Cl/Br/I atoms are temporarily replaced with F (atomic number 9)
2. The substituted SMILES is parsed using fluorine-based SMARTS
3. Results include `substituted: true` and `substitutions: {Cl: X, Br: Y, I: Z}`

This allows leveraging existing fluorine SMARTS patterns for other halogenated compounds.

## Output Format

The JSON output contains:

```json
{
  "run_id": "pfasgroups_halogens_20260216T123456Z",
  "timestamp": "2026-02-16T12:34:56.789012",
  "config": {
    "database": "clinventory",
    "table": "elements_compound",
    "halogens": ["F", "Cl", "Br", "I"],
    "substitute_halogens": false,
    "total_compounds": 1234
  },
  "summary": {
    "success": 1200,
    "failed": 34,
    "with_groups": 567,
    "duration_seconds": 45.678
  },
  "results": [
    {
      "compound_id": 123,
      "smiles": "C(F)(F)C(F)(F)F",
      "formula": "C2F6",
      "status": "success",
      "error": null,
      "pfas_groups": [
        {
          "group_id": 1,
          "group_name": "Perfluoroalkyl",
          "group_alias": "PFA",
          "match_count": 1,
          "chain_lengths": [2],
          "component_smarts": "[CF2]"
        }
      ],
      "substituted": false,
      "substitutions": {},
      "parse_time_seconds": 0.023
    }
  ]
}
```

## Database Output Format

When using `--save-to-db`, results are stored in the `pfasgroups_results` table with the following schema:

```sql
CREATE TABLE pfasgroups_results (
    id SERIAL PRIMARY KEY,
    molecule_id INTEGER NOT NULL,
    run_id TEXT NOT NULL,
    group_id INTEGER,
    group_name TEXT,
    match_count INTEGER,
    chain_lengths JSONB,
    component_smarts JSONB,
    num_components INTEGER,
    substituted BOOLEAN DEFAULT FALSE,
    substitutions JSONB,
    parse_time_seconds DOUBLE PRECISION,
    status TEXT,
    error_message TEXT,
    created_at TIMESTAMP DEFAULT NOW(),
    UNIQUE (molecule_id, run_id, group_id)
);
```

### Table Features

- **One row per PFAS group match** per molecule per run
- **Unique constraint** prevents duplicate entries for the same molecule/run/group combination
- **JSONB columns** for flexible storage of arrays (chain_lengths, component_smarts, substitutions)
- **Indexes** on molecule_id, run_id, group_id, and status for fast queries
- **Automatic table creation** - the table is created if it doesn't exist

### Example Database Queries

**Find molecules with multiple PFAS groups:**
```sql
SELECT molecule_id, COUNT(DISTINCT group_id) as group_count
FROM pfasgroups_results
WHERE group_name IS NOT NULL
GROUP BY molecule_id
HAVING COUNT(DISTINCT group_id) > 1
ORDER BY group_count DESC;
```

**Count molecules by PFAS group:**
```sql
SELECT group_name, COUNT(DISTINCT molecule_id) as molecule_count
FROM pfasgroups_results
WHERE group_name IS NOT NULL
GROUP BY group_name
ORDER BY molecule_count DESC;
```

**Find molecules with substituted halogens:**
```sql
SELECT molecule_id, group_name, substitutions
FROM pfasgroups_results
WHERE substituted = true;
```

**Get average chain lengths by group:**
```sql
SELECT 
    group_name,
    COUNT(*) as matches,
    AVG(jsonb_array_length(chain_lengths)) as avg_chain_length
FROM pfasgroups_results
WHERE group_name IS NOT NULL 
  AND jsonb_array_length(chain_lengths) > 0
GROUP BY group_name
ORDER BY avg_chain_length DESC;
```

**Join with molecules table for analysis:**
```sql
SELECT 
    m.smiles,
    m.formula,
    r.group_name,
    r.match_count,
    r.chain_lengths
FROM pfasgroups_results r
JOIN molecules m ON r.molecule_id = m.id
WHERE r.group_name = 'polyfluoroalkyl'
LIMIT 10;
```

## Query Details

The script executes a PostgreSQL query to find halogenated compounds:

**Default (clinventory database with formula text column):**
```sql
SELECT c.id, c.smiles, c.formula
FROM molecules c
WHERE c.smiles IS NOT NULL 
  AND c.smiles <> ''
  AND c.formula ~ '(F[0-9]*|Cl[0-9]*|Br[0-9]*|I[0-9]*)'
ORDER BY c.id
LIMIT <limit>
OFFSET <offset>
```

**With `--use-halogen-columns` flag (zeropmdbwp2 database with numeric columns):**
```sql
SELECT c.id, c.smiles, c.formula
FROM elements_compound c
WHERE c.smiles IS NOT NULL 
  AND c.smiles <> ''
  AND (c.F > 0 OR c.Cl > 0 OR c.Br > 0 OR c.I > 0)
ORDER BY c.id
LIMIT <limit>
OFFSET <offset>
```

The default query uses PostgreSQL regex matching on the formula column. Use `--use-halogen-columns` for databases with integer columns tracking element counts.

## Performance

- Typical processing rate: 20-50 compounds/second
- Memory usage: ~100-500 MB depending on SMARTS complexity
- RDKit warnings are suppressed via `RDLogger.DisableLog('rdApp.*')`

## Error Handling

- Invalid SMILES → `status: "failed"` with error message
- Parse exceptions → captured in `error` field
- Halogen substitution failures → warning printed, original SMILES used
- Missing halogens in database → zero compounds returned

## Dependencies

- Python 3.6+
- psycopg2
- RDKit
- PFASGroups (local repository)

## Environment Variables

```bash
export DATABASE_NAME=clinventory
export DATABASE_USER=django
export DJANGO_USER=your_password_here
export DATABASE_HOST=localhost
export DATABASE_PORT=5436
```

## Notes

1. **Location**: Script is in PFASGroups root directory, not zeropm_db
2. **Database**: Targets `clinventory` database by default
3. **PFASGroups Import**: Imports from current directory (no path manipulation needed)
4. **Standalone**: No Django dependencies - uses psycopg2 directly
5. **Flexible Schema**: Table and column names are configurable for different database structures
