# Getting Started with PFAS Benchmark Reviewer

## Problem: No Molecules in Review App

If you see an empty review app with no molecules, it means **benchmark data hasn't been generated yet**.

## Solution: Generate Benchmark Data

The review app needs JSON data files to import. Follow these steps:

### Step 1: Activate Python Environment

```powershell
conda activate chem
```

### Step 2: Generate Benchmark Data

Choose ONE of the following options:

#### Option A: OECD Benchmark (Recommended for initial testing)
```powershell
cd C:\Users\luc\git\PFASGroups\benchmark
python enhanced_pfas_benchmark.py
```

This generates:
- `data/pfas_oecd_benchmark_YYYYMMDD_HHMMSS.json` (classification results)
- `data/pfas_timing_benchmark_YYYYMMDD_HHMMSS.json` (timing data)

#### Option B: Enhanced Benchmark with Database
```powershell
cd C:\Users\luc\git\PFASGroups\benchmark
python enhanced_benchmark_with_db.py
```

This generates enhanced dataset with additional molecules.

#### Option C: Run All Benchmarks
```bash
cd C:\Users\luc\git\PFASGroups\benchmark
bash run_all_benchmarks.sh
```

### Step 3: Import Data into Review App

After generating data files:

```powershell
cd C:\Users\luc\git\PFASGroups\benchmark\review-app
node scripts\import-benchmark-data.js
```

You should see output like:
```
Found 2 JSON files to import
Importing pfas_oecd_benchmark_20251226_123456.json...
  Removed 150 duplicate(s) from pfas_oecd_benchmark_20251226_123456.json
✓ Imported 3264 records from pfas_oecd_benchmark_20251226_123456.json
...
✓ Import Statistics:
{
  "molecules": [
    {
      "dataset_type": "oecd",
      "count": 3264
    }
  ]
}
```

### Step 4: Start the Review App

```powershell
cd C:\Users\luc\git\PFASGroups\benchmark\review-app
.\start-dev.ps1
```

Open http://localhost:3000

---

## Quick Reference

### Re-import Latest Data
If you generate new benchmark data and want to refresh the database:

```powershell
# Option 1: Clear and re-import everything
cd review-app
Remove-Item database\pfas_benchmark.db
node scripts\import-benchmark-data.js

# Option 2: Use the quick-start script (includes prompts)
cd ..
.\quick-start.ps1
```

### Check What's Imported
To see import statistics without starting the server:

```powershell
cd review-app
node -e "const db = require('./database/database'); const importer = require('./scripts/import-benchmark-data'); const i = new importer(); i.getImportStats().then(s => console.log(JSON.stringify(s, null, 2)))"
```

---

## Typical Workflow

1. **First Time Setup:**
   ```powershell
   cd C:\Users\luc\git\PFASGroups\benchmark
   .\quick-start.ps1  # Sets up the app (but shows warning about no data)
   ```

2. **Generate Data:**
   ```powershell
   conda activate chem
   python enhanced_pfas_benchmark.py  # Takes 5-30 minutes depending on dataset size
   ```

3. **Import Data:**
   ```powershell
   cd review-app
   node scripts\import-benchmark-data.js
   ```

4. **Review:**
   ```powershell
   .\start-dev.ps1
   # Open http://localhost:3000
   ```

5. **After Making Changes to Algorithms:**
   - Re-run benchmark scripts
   - Clear database and re-import
   - Review changes

---

## Troubleshooting

### "Found 0 JSON files to import"
- Check that `benchmark/data/` directory has `pfas_*.json` files
- Run the Python benchmark scripts to generate data

### "Database is empty"
- Run `node scripts\import-benchmark-data.js` from the `review-app` directory
- Ensure data files exist in `benchmark/data/`

### "No molecules showing in review"
- Check browser console for errors
- Verify server is running (http://localhost:5000/api/molecules should return data)
- Check that database file exists: `review-app/database/pfas_benchmark.db`

### Python import errors
- Ensure `chem` conda environment is activated
- Check that PFASGroups is installed: `pip list | grep PFASGroups`
- Check that PFAS-Atlas is available in the expected path

---

## Data File Format

The review app expects JSON files with this structure:

```json
[
  {
    "molecule_data": {
      "smiles": "C(F)(F)F",
      "molecular_weight": 100.0,
      "num_atoms": 7
    },
    "pfasgroups_result": {
      "detected_groups": ["CF3"],
      "success": true,
      "error": null,
      "execution_time": 0.001234
    },
    "atlas_result": {
      "first_class": "PFAS",
      "second_class": "Type1",
      "success": true,
      "error": null,
      "execution_time": 0.002345
    }
  }
]
```

---

## Need Help?

- Check the main [README.md](../README.md) for overall project documentation
- See [DEDUPLICATION_FEATURE.md](DEDUPLICATION_FEATURE.md) for details on how duplicates are handled
- Review the Python benchmark scripts for data generation options
