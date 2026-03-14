# PFASGroups Benchmark Suite

Comprehensive validation and performance benchmarking for the PFASGroups PFAS
detection algorithm. Benchmarks cover classification accuracy, timing, cross-tool
comparisons, and publication-ready report generation.

---

## Directory Structure

```
benchmark/
├── scripts/                        # All Python analysis and data-collection scripts
├── data/                           # Timestamped benchmark output JSON files
├── reports/                        # Generated LaTeX and JSON report fragments
├── imgs/                           # Generated figures (PDF + PNG)
├── results/                        # Tanimoto / cosine similarity caches and plots
├── review-app/                     # Node.js + React web app for browsing results
├── test_data/                      # Small static datasets used by certain scripts
│
├── run_all_benchmarks.sh           # Main runner — Linux/macOS (see below)
├── run_all_benchmarks.ps1          # Main runner — Windows PowerShell
├── run_clinventory_benchmark.sh    # Clinventory database classification runner
├── run-missing-benchmarks.ps1      # Re-runs the OECD, non-fluorinated, and complex benchmarks
├── run_analysis_reports.ps1        # Re-generates all figures/reports from existing data
├── generate-analysis-reports.ps1   # Legacy alias for run_analysis_reports.ps1
├── quick-start.sh                  # One-shot setup + review-app launch (Linux)
├── quick-start.ps1                 # One-shot setup + review-app launch (Windows)
├── setup-review-app.sh             # Installs Node.js deps + initialises SQLite DB (Linux)
├── setup-review-app.ps1            # Same for Windows
│
├── _rebuild.py                     # Rewrites fingerprint_structure_analysis.ipynb in-place
├── create_notebook.py              # Generates fingerprint_structure_analysis.ipynb programmatically
├── fingerprint_structure_analysis.ipynb  # Interactive Tanimoto / structure analysis notebook
└── targets_fingerprint_REACH.md   # Notes on REACH-relevant PFAS target structures
```

---

## Shell / PowerShell Scripts

### `run_all_benchmarks.sh` / `run_all_benchmarks.ps1`

**Main entry point.** Runs every benchmark type in sequence, then generates a
unified HTML report. Requires the `pfasatlas` conda environment.

```bash
# Linux
./run_all_benchmarks.sh
./run_all_benchmarks.sh --reuse-timing   # skip re-running timing benchmarks
./run_all_benchmarks.sh --quick          # reduced molecule counts for smoke testing

# Windows
.\run_all_benchmarks.ps1
```

Calls in order:
1. `scripts/enhanced_pfas_benchmark.py` — options 1, 2, 3 (×3 profiles), 4, 5
2. `scripts/test_highly_branched.py`
3. `scripts/validate_telomers.py`
4. `scripts/benchmark_pfas_definitions.py`
5. `scripts/generate_unified_report.py`
6. `review-app` database import scripts

---

### `run_clinventory_benchmark.sh`

Classifies all halogenated molecules in the PostgreSQL `clinventory` database
using **PFASGroups (HalogenGroups)**, **PFASGroups (F-only)**, and
**PFAS-Atlas**, then compares timing and PFAS-detection performance.

```bash
./run_clinventory_benchmark.sh [OPTIONS]

Options:
  --limit N          Process only first N molecules (quick test: --quick = 1 000)
  --all-molecules    Include non-halogenated molecules
  --db-password P    Database password (prompted if omitted)
  --skip-hg          Skip HalogenGroups classification
  --skip-atlas       Skip PFAS-Atlas classification
  --skip-pfasgroups  Skip PFASGroups (F-only) classification
  --skip-classify    Skip all classification steps (use existing JSON files)
```

Calls:
- `scripts/classify_halogengroups_clinventory.py`
- `scripts/classify_pfasatlas_clinventory.py`
- `scripts/classify_pfasgroups_clinventory.py`
- `scripts/compare_clinventory_classifiers.py`
- `scripts/compare_oecd_clinventory_timing.py`
- `scripts/generate_clinventory_latex.py`

---

### `run-missing-benchmarks.ps1`

Re-runs only the **OECD** (option 2), **Non-Fluorinated** (option 4), and
**Complex Branched** (option 5) benchmarks when they are missing or need
updating, without repeating the time-consuming functional-groups benchmark.
Calls `scripts/enhanced_pfas_benchmark.py` for each, then moves output JSON
files to `data/`.

---

### `run_analysis_reports.ps1` / `generate-analysis-reports.ps1`

Re-generates **all figures and report files** from previously collected
benchmark data **without re-running data collection**. Accepts `-Env` and
`-AtlasEnv` parameters to select conda environments.

```powershell
.\run_analysis_reports.ps1                        # default env: chem
.\run_analysis_reports.ps1 -Env pfasatlas -SkipAtlas
```

Runs in order (scripts 1–11):
1. `analyze_timing.py` — timing performance analysis
2. `analyze_timing_models.py` — exponential scaling model fits
3. `compare_timing_profiles.py` — multi-profile timing comparison
4. `analyze_definitions_benchmark.py` — per-definition metrics and confusion matrices
5. `enhanced_analysis.py` — combined enhanced + OECD heatmap / Sankey analysis
6. `generate_unified_report.py` — unified HTML report
7. `generate_telomer_report.py` — telomer HTML report
8. `analyze_benchmarks_simple.py` — lightweight `benchmark_summary.json`
9. `generate_latex_tables.py` — publication LaTeX tables
10. `generate_telomer_latex.py` — telomer LaTeX fragments
11. `generate_clinventory_latex.py` — clinventory LaTeX fragments

---

### `quick-start.sh` / `quick-start.ps1`

One-shot bootstrap: installs Node.js dependencies, initialises the SQLite
database, imports the latest benchmark data, and prints instructions for
starting the review app in development or production mode.

---

### `setup-review-app.sh` / `setup-review-app.ps1`

Installs `review-app` Node.js dependencies, merges React component sources,
and initialises the SQLite database. Called automatically by `quick-start.*`.

---

## `scripts/` — Python Scripts

### Data Collection

| Script | Purpose |
|--------|---------|
| `enhanced_pfas_benchmark.py` | Interactive benchmark menu (options 1–6): functional groups, OECD validation, timing (3 profiles), non-fluorinated exclusion, complex branched structures. Primary data-collection script. |
| `benchmark_pfas_definitions.py` | Benchmarks five PFAS definitions (OECD, EU Restriction, OPPT 2023, UK, PFASTRUCTv5) for accuracy, specificity, sensitivity, and inter-definition agreement. |
| `test_highly_branched.py` | Tests functional groups 29–59 at varying branching distances (0–2 bonds) from perfluorinated chains; compares PFASGroups vs PFAS-Atlas detection. |
| `validate_telomers.py` | Validates telomer detection against a PubChem fluorotelomer SDF dataset. |
| `classify_halogengroups_clinventory.py` | Classifies all halogen-containing molecules in the PostgreSQL `clinventory` database with **HalogenGroups** (full halogen set). Writes `data/PFASGroups_clinventory_<TS>.json`. |
| `classify_pfasgroups_clinventory.py` | Same as above but restricts HalogenGroups to fluorine only (`halogens='F'`). Writes `data/pfasgroups_clinventory_<TS>.json`. |
| `classify_pfasatlas_clinventory.py` | Classifies all clinventory molecules with **PFAS-Atlas** (requires `pfasatlas` conda env). Writes `data/pfasatlas_clinventory_<TS>.json`. |
| `fingerprint_tanimoto_benchmark.py` | Two-phase Tanimoto benchmark: (1) selects the top-5 most discriminating PFASGroups fingerprint configurations; (2) compares them against TxP_PFAS (Richard et al. 2023 CSRML). |
| `compare_fingerprints_vs_txppfas.py` | Detailed structural comparison of PFASGroups binary fingerprints vs the 129-bit TxP_PFAS CSRML fingerprint. |
| `compare_oecd_clinventory_timing.py` | Head-to-head timing of PFASGroups vs PFAS-Atlas on the OECD PFAS CSV dataset and the clinventory database, with molecular-complexity correlation analysis. |

### Analysis and Reporting

| Script | Purpose |
|--------|---------|
| `analyze_timing.py` | Reads a timing benchmark JSON and generates Plotly visualisations: distribution histograms, scatter plots (time vs size), scaling curves. |
| `analyze_timing_models.py` | Fits exponential scaling models to timing data for three profiles (all features / no effective graph resistance / no metrics). |
| `compare_timing_profiles.py` | Compares the three timing profiles, fits per-profile exponential curves, and writes a LaTeX summary and PNG/PDF figures. |
| `analyze_benchmarks_simple.py` | Minimal aggregation of all benchmark JSON files into `reports/benchmark_summary.json`. |
| `analyze_complex.py` | Analyses the complex-branched-structures benchmark JSON and prints group-level statistics. |
| `analyze_definitions_benchmark.py` | Generates per-definition performance tables, confusion matrices, agreement patterns, and an interactive HTML report from `benchmark_pfas_definitions.py` output. |
| `enhanced_analysis.py` | Produces heatmaps, Sankey diagrams, and detailed cross-tool statistics by combining the enhanced and OECD benchmark JSON files. |
| `compare_clinventory_classifiers.py` | Loads all three clinventory JSON files and generates timing comparison plots (by atom-count bracket, by halogen class) and PFAS-detection agreement plots. Saves a unified `data/clinventory_comparison_<TS>.json`. |
| `generate_unified_report.py` | Combines all benchmark results into a single self-contained HTML report by auto-detecting the most recent files in `data/`. |
| `generate_latex_tables.py` | Writes publication-ready LaTeX fragments (tables, timing statistics, figure captions) to `reports/`. |
| `generate_telomer_report.py` | Reads `data/telomer_validation_results.json` and produces an HTML report with group detection rates and examples for the review app. |
| `generate_telomer_latex.py` | Writes LaTeX table fragments for the telomer validation section of the article. |
| `generate_clinventory_latex.py` | Writes LaTeX table and narrative text fragments for the clinventory comparison section of the article. |

### Utilities

| Script | Purpose |
|--------|---------|
| `add_test_metadata.py` | Post-processes benchmark outputs: copies legacy files from `scripts/data/` to `benchmark/data/` and adds minimal metadata to definition-benchmark files. |
| `compare_timing_profiles.py` | Overlays exponential fits for all three timing profiles and writes a summary table. |

---

## Other Directories

### `data/`

Timestamped JSON output files produced by the data-collection scripts.
Common naming patterns:

| Pattern | Produced by |
|---------|------------|
| `pfas_enhanced_benchmark_<TS>.json` | `enhanced_pfas_benchmark.py` (option 1) |
| `pfas_oecd_benchmark_<TS>.json` | `enhanced_pfas_benchmark.py` (option 2) |
| `pfas_timing_benchmark_{profile}_<TS>.json` | `enhanced_pfas_benchmark.py` (option 3) |
| `pfas_non_fluorinated_benchmark_<TS>.json` | `enhanced_pfas_benchmark.py` (option 4) |
| `pfas_complex_branched_benchmark_<TS>.json` | `enhanced_pfas_benchmark.py` (option 5) |
| `pfas_highly_branched_benchmark_<TS>.json` | `test_highly_branched.py` |
| `pfas_definitions_benchmark_<TS>.json` | `benchmark_pfas_definitions.py` |
| `telomer_validation_results.json` | `validate_telomers.py` |
| `PFASGroups_clinventory_<TS>.json` | `classify_halogengroups_clinventory.py` |
| `pfasgroups_clinventory_<TS>.json` | `classify_pfasgroups_clinventory.py` |
| `pfasatlas_clinventory_<TS>.json` | `classify_pfasatlas_clinventory.py` |
| `clinventory_comparison_<TS>.json` | `compare_clinventory_classifiers.py` |
| `oecd_clinventory_timing_<TS>.json` | `compare_oecd_clinventory_timing.py` |
| `timing_analysis_<TS>.json` | `analyze_timing.py` |
| `PubChem_fluorotelomers.sdf` | Static — PubChem fluorotelomer reference dataset |

### `reports/`

LaTeX fragments and JSON summaries written by the `generate_*` and `analyze_*`
scripts, ready to be `\input{}`-ed into the article:

- `*.tex` — tables, figure captions, narrative text
- `benchmark_summary.json` — aggregated benchmark statistics
- `timing_analysis_*.json` — timing model parameters
- `timing_profiles_summary.{json,tex}` — multi-profile comparison

### `imgs/`

PDF and PNG figures produced by the analysis scripts, organised by benchmark
type (clinventory timing, OECD/clinventory comparison, timing scaling, etc.).

### `results/`

Cached intermediate results for the Tanimoto / cosine fingerprint benchmarks:

- `tanimoto/` — pairwise Tanimoto similarity matrices
- `cosine/` — cosine similarity matrices
- `plots/` — fingerprint discrimination figures
- `cache/` — pre-computed fingerprint vectors
- `lookups/` — cached SMARTS-match lookups

### `review-app/`

A Node.js (Express) + React web application for interactively browsing
benchmark results. Features:
- Prioritised display of misclassified molecules
- Manual review with accept / reject / unsure buttons
- RDKit.js 2-D molecular visualisation
- Real-time accuracy computation
- JSON/CSV export

Start in development mode:
```bash
cd review-app
./start-dev.sh         # Linux
.\start-dev.ps1        # Windows — serves on http://localhost:3000
```

Start in production mode:
```bash
cd review-app
./start-prod.sh        # Linux — serves on http://localhost:5000
```

### `test_data/`

Small curated molecule sets and reference files used as static inputs by
specific scripts (e.g. `benchmark_test_compounds.csv` for the definitions
benchmark, REACH target lists).

---

## Conda Environments

| Environment | Used for |
|-------------|---------|
| `pfasatlas` | `run_all_benchmarks.sh`, `classify_pfasatlas_clinventory.py`, any script that imports PFAS-Atlas |
| `chem` | All other scripts (PFASGroups, RDKit, pandas, plotly, scipy) |

Scripts that need both tools accept a `--no-atlas` flag to skip the PFAS-Atlas
classifiers and run with the `chem` environment only.

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
