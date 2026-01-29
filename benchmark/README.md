# PFAS Benchmark

## Actions Required

### Step 1: Run Missing/Updated Benchmarks

```powershell
# Activate your conda environment
conda activate pfasatlas

# Navigate to benchmark directory
cd c:\Users\luc\git\PFASGroups\benchmark
cd /home/luc/git/PFASGroups/benchmark

# Run all benchmarks
python enhanced_pfas_benchmark.py
# When prompted, enter: 6 for all


```

### Step 2: Import (or reimport) All Data (Clean Import)

```powershell
cd c:\Users\luc\git\PFASGroups\benchmark\review-app

# or in linux

```bash
cd /home/luc/git/PFASGroups/benchmark/review-app

# Clear and reimport all data with fixed deduplication
node scripts/reimport-all.js
```

### Step 3: Calculate Formulas for New Molecules

```powershell
# Calculate molecular formulas for newly imported molecules
python scripts/calculate-formulas.py
```

### Step 4: Generate Analysis Reports

**Prerequisites:**
```powershell
# Install required package for image generation
pip install kaleido
```

**Windows PowerShell:**
```powershell
cd c:\Users\luc\git\PFASGroups\benchmark

# Get latest benchmark files (PowerShell doesn't expand globs like bash)
$timingFile = Get-ChildItem data\pfas_timing_benchmark_*.json | Sort-Object LastWriteTime -Descending | Select-Object -First 1 -ExpandProperty FullName
$complexFile = Get-ChildItem data\pfas_complex_branched_benchmark_*.json | Sort-Object LastWriteTime -Descending | Select-Object -First 1 -ExpandProperty FullName
$enhancedFile = Get-ChildItem data\pfas_enhanced_benchmark_*.json | Sort-Object LastWriteTime -Descending | Select-Object -First 1 -ExpandProperty FullName
$oecdFile = Get-ChildItem data\pfas_oecd_benchmark_*.json | Sort-Object LastWriteTime -Descending | Select-Object -First 1 -ExpandProperty FullName

# Generate analyses
python analyze_timing.py $timingFile
python analyze_complex.py $complexFile
python enhanced_analysis.py $enhancedFile $oecdFile

# Or use the helper script:
.\generate-analysis-reports.ps1
```

**Linux/macOS Bash:**
```bash
cd /home/luc/git/PFASGroups/benchmark

# Generate timing analysis
python analyze_timing.py data/pfas_timing_benchmark_*.json

# Generate complex branched analysis  
python analyze_complex.py data/pfas_complex_branched_benchmark_*.json

# Generate enhanced analysis
python enhanced_analysis.py data/pfas_enhanced_benchmark_*.json data/pfas_oecd_benchmark_*.json
```

**Move files to review-app:**
```powershell
# Windows
New-Item -ItemType Directory -Force -Path review-app\analysis_reports\figures
Move-Item *_analysis.json review-app\analysis_reports\ -Force
Move-Item *.png,*.svg review-app\analysis_reports\figures\ -Force

# Linux
mkdir -p review-app/analysis_reports/figures
mv *_analysis.json review-app/analysis_reports/
mv *.png *.svg review-app/analysis_reports/figures/
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
├── scripts/          # 34 Python scripts
├── reports/          # 9 HTML reports (all in one place)
├── data/             # 14 data files (JSON, CSV, NPZ, TSV)
├── imgs/             # 30 images (PNG, PDF, SVG)
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
