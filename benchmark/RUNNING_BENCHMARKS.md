# Running the Complete PFAS Benchmark Suite

## Quick Start (Windows)

The complete benchmark suite can now be run with a single command:

```powershell
# From the benchmark directory
cd c:\Users\luc\git\PFASGroups\benchmark

# Activate conda environment (if not already active)
conda activate chem

# Run everything
.\run_all_benchmarks.ps1
```

## What This Does

The unified script performs the following steps automatically:

### 1. **Run All Benchmarks** (15-30 minutes)
- ✅ Enhanced Functional Groups Benchmark (~1675 molecules)
- ✅ OECD Validation Benchmark (~4000 molecules)  
- ✅ Timing Performance Benchmark (200 molecules × 10 iterations)
- ✅ Non-Fluorinated Exclusion Benchmark (50+ molecules)
- ✅ Complex Branched PFAS Benchmark (300 molecules)

### 2. **Generate Analysis Reports** (2-5 minutes)
- ⏱️ Timing performance analysis (`analyze_timing.py`)
- 🧬 Complex branched analysis (`analyze_complex.py`)
- 🔬 Enhanced + OECD analysis (`enhanced_analysis.py`)
- Generates JSON reports and visualization figures (SVG & PNG)

### 3. **Update Review App Database** (5-10 minutes)
- 📦 Backs up existing database
- 🗑️ Clears old data
- 📥 Imports all benchmark results
- 🧪 Calculates molecular formulas
- Organizes analysis reports in `review-app/analysis_reports/`

### 4. **File Organization**
```
benchmark/
├── data/                       # All benchmark JSON files
│   ├── pfas_enhanced_benchmark_*.json
│   ├── pfas_oecd_benchmark_*.json
│   ├── pfas_timing_benchmark_*.json
│   ├── pfas_non_fluorinated_benchmark_*.json
│   └── pfas_complex_branched_benchmark_*.json
├── html/                       # Unified HTML reports
│   └── unified_pfas_benchmark_report_*.html
├── imgs/                       # Generated plots
└── review-app/
    ├── database/
    │   ├── pfas_benchmark.db   # Updated database
    │   └── pfas_benchmark.db.backup_*  # Auto backups
    └── analysis_reports/       # Analysis results
        ├── timing_analysis.json
        ├── complex_analysis.json
        ├── enhanced_analysis.json
        └── figures/            # SVG & PNG figures
```

## Using the Review App

After running the benchmark suite:

```powershell
# Start the review app server
cd review-app
node server.js

# Open browser to http://localhost:5000
```

### Review App Features

1. **Dashboard** - Overview of all datasets and statistics
2. **Review Molecules** - Manual validation interface
3. **Accuracy Report** - Detailed accuracy metrics
4. **Analysis Reports** - NEW! View and export analysis results
   - Timing Performance tab
   - Complex Branched tab
   - Enhanced Analysis tab
   - Add descriptions with Markdown
   - Export to Markdown with embedded figures

## Linux/macOS

Use the bash script instead:

```bash
cd /path/to/PFASGroups/benchmark
chmod +x run_all_benchmarks.sh
./run_all_benchmarks.sh
```

## Expected Results

After completion, you should have:

- **~6000+ molecules** in the database:
  - OECD: ~4000
  - Enhanced: ~1675
  - Timing: ~600 (200 molecules × 3 replicates)
  - Complex Branched: ~300
  - Non-Fluorinated: ~50+

- **3 analysis reports** with figures and statistics
- **Unified HTML report** with comprehensive results
- **Review app ready** with all data imported

## Troubleshooting

### If script fails:

1. **Check conda environment:**
   ```powershell
   conda activate chem
   python --version  # Should show Python 3.x
   node --version    # Should show Node.js version
   ```

2. **Check dependencies:**
   ```powershell
   pip install rdkit pandas numpy matplotlib
   npm install  # In review-app/client directory if needed
   ```

3. **Check file permissions:**
   - Ensure write access to benchmark directory
   - Close any open database connections

4. **Manual step-by-step:**
   If the unified script fails, you can run steps manually:
   ```powershell
   # 1. Run benchmarks individually
   python enhanced_pfas_benchmark.py  # Choose option 1-5

   # 2. Run analysis scripts
   python analyze_timing.py data\pfas_timing_benchmark_*.json
   python analyze_complex.py data\pfas_complex_branched_benchmark_*.json
   python enhanced_analysis.py data\pfas_enhanced_benchmark_*.json data\pfas_oecd_benchmark_*.json

   # 3. Import to database
   cd review-app
   node scripts\import-benchmark-data.js
   python scripts\calculate-formulas.py
   ```

## Time Estimates

- **Quick test run:** 5-10 minutes (small datasets)
- **Full benchmark:** 30-45 minutes (all 6000+ molecules)
- **With analysis:** Add 10-15 minutes
- **Total:** ~1 hour for complete suite

## Notes

- **Triplicates preserved:** Timing data keeps all 3 replicates for statistical analysis
- **Auto-backup:** Database is automatically backed up before clearing
- **Deduplication:** Non-timing datasets have duplicates removed during import
- **Formula calculation:** May take 10-15 minutes for 6000+ molecules
- **Analysis reports:** JSON format compatible with the review app UI
