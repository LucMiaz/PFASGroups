# PFAS Benchmark - Next Steps

## Summary of Issues Found and Fixed

### 1. ✅ OECD Benchmark - Verified but Not Run
**Status:** The `run_all_benchmarks.sh` script IS configured to run OECD when option is empty or '6'
**Issue:** The OECD benchmark was never actually executed
**Action Required:** Run OECD benchmark manually

### 2. ✅ Timing Triplicates - Fixed
**Status:** Fixed in import script
**Change:** Modified `import-benchmark-data.js` to preserve timing triplicates for statistical analysis
**Note:** Triplicates are intentional for timing measurements

### 3. ✅ Non-Fluorinated Molecules - Fixed  
**Status:** Fixed in benchmark code
**Issue:** Hardcoded SMILES contained fluorine (e.g., "FC(F)(F)CC(=O)O")
**Fix:** Replaced with truly non-fluorinated molecules (e.g., "CCC(=O)O")
**New Count:** Will generate 50+ molecules (12 variants × rounds)

### 4. ✅ Complex Branched - Fixed
**Status:** Fixed in benchmark code
**Issue:** Only 20 molecules per test (6 tests = 120 total)  
**Fix:** Increased to 50 molecules per test (6 tests = 300 total)

### 5. ✅ Analysis Reports Page - Added
**Status:** New page created in review app
**Features:**
- View timing, complex, and enhanced analysis results
- Add/edit descriptions with Markdown support
- Template insertion for each report type
- Export to Markdown with embedded figures (SVG & PNG)
- Save descriptions persistently

## Actions Required

### Step 1: Run Missing/Updated Benchmarks

```powershell
# Activate your conda environment
conda activate chem

# Navigate to benchmark directory
cd c:\Users\luc\git\PFASGroups\benchmark

# Run OECD Benchmark (4000+ molecules)
python enhanced_pfas_benchmark.py
# When prompted, enter: 2

# Run Non-Fluorinated Benchmark (50+ molecules, updated)
python enhanced_pfas_benchmark.py  
# When prompted, enter: 4

# Run Complex Branched Benchmark (300+ molecules, updated)
python enhanced_pfas_benchmark.py
# When prompted, enter: 5

# Move generated files to data directory
Move-Item pfas_*_benchmark_*.json data\
```

### Step 2: Reimport All Data (Clean Import)

```powershell
cd c:\Users\luc\git\PFASGroups\benchmark\review-app

# Clear and reimport all data with fixed deduplication
node scripts/reimport-all.js
```

### Step 3: Calculate Formulas for New Molecules

```powershell
# Calculate molecular formulas for newly imported molecules
python scripts/calculate-formulas.py
```

### Step 4: Generate Analysis Reports

```powershell
cd c:\Users\luc\git\PFASGroups\benchmark

# Generate timing analysis
python analyze_timing.py data/pfas_timing_benchmark_*.json

# Generate complex branched analysis  
python analyze_complex.py data/pfas_complex_branched_benchmark_*.json

# Generate enhanced analysis
python enhanced_analysis.py data/pfas_enhanced_benchmark_*.json data/pfas_oecd_benchmark_*.json

# Create analysis_reports directory and copy results
New-Item -ItemType Directory -Force -Path review-app/analysis_reports
Copy-Item *.json review-app/analysis_reports/
```

### Step 5: Restart Server

```powershell
cd c:\Users\luc\git\PFASGroups\benchmark\review-app

# Stop any running node processes
Get-Process | Where-Object {$_.ProcessName -eq "node"} | Stop-Process -Force

# Start server
node server.js
```

### Step 6: Access Analysis Reports

1. Open browser to: http://localhost:5000
2. Click "Analysis Reports" in navigation
3. View three tabs:
   - ⏱️ Timing Performance
   - 🧬 Complex Branched  
   - 🔬 Enhanced Analysis
4. For each report:
   - Click "📝 Insert Template" to add description template
   - Edit the description with your analysis
   - Click "💾 Save" to save your description
   - Click "📄 Export to Markdown" to download report with figures

## Expected Results After Completion

### Database Contents
- **OECD:** ~4000 molecules
- **Enhanced:** ~1675 molecules
- **Timing:** ~200 molecules (with triplicates = ~600 total)
- **Complex Branched:** ~300 molecules (6 tests × 50)
- **Non-Fluorinated:** ~50+ molecules (3 groups × 12+ variants)
- **Total:** ~6000+ molecules (unique ~5800)

### Files Generated
```
benchmark/
├── data/
│   ├── pfas_oecd_benchmark_YYYYMMDD_HHMMSS.json (NEW)
│   ├── pfas_non_fluorinated_benchmark_YYYYMMDD_HHMMSS.json (UPDATED)
│   ├── pfas_complex_branched_benchmark_YYYYMMDD_HHMMSS.json (UPDATED)
│   ├── pfas_enhanced_benchmark_YYYYMMDD_HHMMSS.json (existing)
│   └── pfas_timing_benchmark_YYYYMMDD_HHMMSS.json (existing)
├── analysis_reports/
│   ├── timing_analysis.json
│   ├── timing_description.md
│   ├── complex_analysis.json
│   ├── complex_description.md
│   ├── enhanced_analysis.json
│   ├── enhanced_description.md
│   └── figures/
│       ├── timing_fig1.svg
│       ├── timing_fig1.png
│       └── ...
└── review-app/
    └── database/
        └── pfas_benchmark.db (updated with all molecules)
```

## Description Templates

The Analysis Reports page provides three description templates:

### Timing Performance Template
- Overview
- Key Findings  
- Performance Characteristics (PFASGroups vs PFAS-Atlas)
- Scalability Analysis
- Conclusions

### Complex Branched Template
- Overview
- Test Molecules
- Detection Results
- Challenges Identified
- Conclusions

### Enhanced Analysis Template
- Overview
- Coverage Analysis
- Accuracy Results (comparison)
- Conclusions

## Notes

1. **Triplicates in Timing:** These are intentional for statistical analysis - DO NOT remove
2. **Deduplication:** Only non-timing datasets have duplicates removed
3. **Molecular Formulas:** Calculated in batch after import for performance
4. **Analysis Reports:** Figures are saved as both SVG and PNG for flexibility
5. **Markdown Export:** Includes relative paths to figures in `figures/` subdirectory

## Troubleshooting

### If Python is not found:
```powershell
conda activate chem
# Then run python commands
```

### If import fails:
```powershell
# Check data directory contents
ls c:\Users\luc\git\PFASGroups\benchmark\data\*.json

# Verify files aren't locked
Get-Process | Where-Object {$_.ProcessName -eq "python" -or $_.ProcessName -eq "node"}
```

### If server won't start:
```powershell
# Kill processes on port 5000
Get-NetTCPConnection -LocalPort 5000 | Select-Object -ExpandProperty OwningProcess | ForEach-Object { Stop-Process -Id $_ -Force }

# Then restart
node server.js
```
