# PFASGroups Parsing Guide for Clinventory Database

## Overview

The `apply_pfasgroups_all_halogens.py` script has been updated to:

1. ✅ **Parse ALL molecules** in the database (not just those with matches)
2. ✅ **Save all results** including molecules without PFAS group matches  
3. ✅ **Fixed database connection** parameters
4. ✅ **Better batch processing** for large datasets

## Key Changes Made

### 1. New --all-molecules Flag
```bash
--all-molecules  # Process ALL molecules in the table, not just halogenated ones
```

### 2. Improved Data Collection
- Previously: Only molecules with PFAS group matches were saved
- Now: ALL successfully parsed molecules are tracked
- Result: Complete coverage of your database

### 3. Fixed Database Saving
- Previously: Used wrong parameter name `conn=` for ResultsModel.to_sql()
- Now: Direct database inserts using pandas with proper connection handling
- Result: More reliable and efficient database operations

## Usage Examples

### Process ALL molecules in clinventory and save to database

```bash
python apply_pfasgroups_all_halogens.py \
  --db-name clinventory \
  --db-user django \
  --db-password YOUR_PASSWORD \
  --db-host localhost \
  --db-port 5432 \
  --table-name molecules \
  --all-molecules \
  --save-to-db \
  --batch-size 500
```

### Process only halogenated compounds (F, Cl, Br, I)

```bash
python apply_pfasgroups_all_halogens.py \
  --db-name clinventory \
  --db-user django \
  --db-password YOUR_PASSWORD \
  --table-name molecules \
  --halogens "F,Cl,Br,I" \
  --save-to-db \
  --batch-size 500
```

### Process only fluorinated compounds with halogen substitution

```bash
python apply_pfasgroups_all_halogens.py \
  --db-name clinventory \
  --db-user django \
  --db-password YOUR_PASSWORD \
  --table-name molecules \
  --halogens "F" \
  --substitute-halogens \
  --save-to-db \
  --batch-size 500
```

### Test run on limited dataset with JSON output

```bash
python apply_pfasgroups_all_halogens.py \
  --db-name clinventory \
  --db-user django \
  --db-password YOUR_PASSWORD \
  --table-name molecules \
  --all-molecules \
  --limit 1000 \
  --output-json results_test.json
```

## Database Tables Created

The script creates/populates three tables:

### 1. `pfasgroups`
PFAS group definitions (created once)

| Column | Type | Description |
|--------|------|-------------|
| id | INTEGER | Primary key, group ID |
| name | TEXT | Group name |
| alias | TEXT | Group alias |
| category | TEXT | Group category |
| pathway_type | TEXT | Pathway type |
| smarts_primary | TEXT | Primary SMARTS pattern |
| smarts_secondary | TEXT | Secondary SMARTS pattern |
| created_at | TIMESTAMP | Creation timestamp |

### 2. `pfasgroups_in_molecules`
PFAS group matches in molecules

| Column | Type | Description |
|--------|------|-------------|
| id | SERIAL | Primary key |
| smiles | TEXT | Molecule SMILES |
| molecule_id | INTEGER | Foreign key to molecules table |
| group_id | INTEGER | Foreign key to pfasgroups table |
| group_name | TEXT | Group name |
| match_count | INTEGER | Number of components matched |
| created_at | TIMESTAMP | Creation timestamp |

### 3. `components_in_molecules`
Individual component details

| Column | Type | Description |
|--------|------|-------------|
| id | SERIAL | Primary key |
| smiles | TEXT | Molecule SMILES |
| molecule_id | INTEGER | Foreign key to molecules table |
| group_id | INTEGER | Foreign key to pfasgroups table |
| group_name | TEXT | Group name |
| smarts_label | TEXT | Component type (e.g., "Perfluoroalkyl") |
| component_atoms | TEXT | Comma-separated atom indices |
| created_at | TIMESTAMP | Creation timestamp |

## Important Notes

### Environment Variables
You can use environment variables instead of command-line arguments:
```bash
export DATABASE_NAME=clinventory
export DATABASE_USER=django
export DJANGO_USER=your_password
export DATABASE_HOST=localhost
export DATABASE_PORT=5432
```

### Performance Tips

1. **Batch size**: Adjust `--batch-size` based on memory (100-1000 works well)
2. **Limit for testing**: Use `--limit 1000` to test on a subset first
3. **Progress monitoring**: Progress is printed every batch-size compounds
4. **Database indexing**: The script creates indexes automatically for better query performance

### What Gets Saved

- **All molecules**: Tracked in internal results
- **Only molecules with matches**: Saved to database tables
- **Statistics**: Printed to console showing:
  - Total processed
  - Success/failed counts
  - Molecules with PFAS groups
  - Total matches and averages

## Verification Queries

After running, verify the results:

```sql
-- Count total molecules with PFAS groups
SELECT COUNT(DISTINCT molecule_id) 
FROM pfasgroups_in_molecules;

-- Count total PFAS group matches
SELECT COUNT(*) 
FROM pfasgroups_in_molecules;

-- Count total components
SELECT COUNT(*) 
FROM components_in_molecules;

-- Top 10 most common PFAS groups
SELECT group_name, COUNT(*) as count
FROM pfasgroups_in_molecules
GROUP BY group_name
ORDER BY count DESC
LIMIT 10;

-- Molecules with most PFAS groups
SELECT molecule_id, smiles, COUNT(DISTINCT group_id) as num_groups
FROM pfasgroups_in_molecules
GROUP BY molecule_id, smiles
ORDER BY num_groups DESC
LIMIT 10;
```

## Troubleshooting

### No results saved
- Check that molecules have valid SMILES
- Verify database connection parameters
- Use `--output-json` to see what's being parsed

### Connection errors
- Verify PostgreSQL is running
- Check username/password
- Confirm database name exists

### Out of memory
- Reduce `--batch-size`
- Use `--limit` to process in chunks with `--offset`

### Slow performance
- Check database indexes (created automatically)
- Consider processing in parallel batches using offset/limit
- Monitor PostgreSQL performance

## Next Steps

1. **Run a test**: Start with `--limit 1000` to verify everything works
2. **Full run**: Remove `--limit` to process all molecules
3. **Verify**: Use SQL queries above to check results
4. **Analyze**: Query the tables to explore PFAS groups in your database
