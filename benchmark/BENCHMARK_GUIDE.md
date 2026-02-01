# PFAS Benchmark Suite Documentation

## Overview

The PFAS Benchmark Suite is a comprehensive testing framework for validating the accuracy and performance of PFAS (Per- and Polyfluoroalkyl Substances) detection tools, specifically PFASGroups and PFAS-Atlas. It includes multiple benchmark types, analysis scripts, and a review application for detailed result exploration.

## Quick Start

### Running All Benchmarks

**Linux/macOS:**
```bash
cd benchmark
./run_all_benchmarks.sh
```

**Windows (PowerShell):**
```powershell
cd benchmark
.\run_all_benchmarks.ps1
```

### Prerequisites

- **Python Environment**: Conda environment with required packages:
  - PFASGroups library
  - PFAS-Atlas library
  - RDKit
  - Matplotlib, Seaborn, Pandas
  - Other dependencies from requirements.txt

- **Node.js**: Required for database import and review app (optional)

- **Activation**: Activate your conda environment before running:
  ```bash
  conda activate pfasatlas  # or your environment name
  ```

## What the Benchmarks Test

### 1. Functional Groups Benchmark
**Script**: `enhanced_pfas_benchmark.py` (option 1)  
**Purpose**: Tests basic PFAS functional group detection accuracy  
**Molecules**: ~3,820 test molecules across 60+ functional groups  
**Output**: `data/pfas_enhanced_benchmark_*.json`

Tests the ability to correctly identify functional groups in PFAS molecules, including:
- Perfluoroalkyl and polyfluoroalkyl groups
- Carboxylic acids, sulfonic acids, phosphonic acids
- Ethers, esters, amides
- Aromatic and cyclic structures
- Various functional head groups

### 2. OECD Validation Benchmark
**Script**: `enhanced_pfas_benchmark.py` (option 2)  
**Purpose**: Validates detection against OECD PFAS database  
**Molecules**: OECD reference PFAS compounds  
**Output**: `data/pfas_oecd_benchmark_*.json`

Compares detection results against the authoritative OECD PFAS database to ensure compliance with international standards.

### 3. Timing Performance Benchmark
**Script**: `enhanced_pfas_benchmark.py` (option 3)  
**Purpose**: Measures execution speed and scalability  
**Molecules**: Variable sizes (10 to 500 atoms)  
**Output**: `data/pfas_timing_benchmark_*.json`

Tests computational performance by measuring:
- Processing time vs. molecule size
- Scaling behavior (linear, polynomial, exponential)
- System performance under load
- Both PFASGroups and PFAS-Atlas performance

### 4. Non-Fluorinated Exclusion Benchmark
**Script**: `enhanced_pfas_benchmark.py` (option 4)  
**Purpose**: Ensures proper exclusion of non-PFAS molecules  
**Molecules**: Non-fluorinated compounds  
**Output**: `data/pfas_non_fluorinated_benchmark_*.json`

Validates that the tools correctly reject molecules without fluorine or insufficient fluorination to be classified as PFAS.

### 5. Complex Branched Structures Benchmark
**Script**: `enhanced_pfas_benchmark.py` (option 5)  
**Purpose**: Tests detection on structurally complex PFAS  
**Molecules**: 300 complex branched molecules  
**Output**: `data/pfas_complex_branched_benchmark_*.json`

Tests challenging molecular structures including:
- Highly branched carboxylic acids
- Branched sulfonic acids
- Branched ether chains
- Multi-functional branched compounds
- Cyclic branched structures
- Aromatic branched structures

### 6. Highly Branched Compounds Test
**Script**: `test_highly_branched.py`  
**Purpose**: Tests functional groups on perfluorinated components  
**Molecules**: Functional groups 29-59 at various branching distances  
**Output**: `data/highly_branched_pfas_test_*.json`

Systematically tests how functional group detection handles branching by placing test groups at different distances from perfluorinated carbon chains.

### 7. Telomer Validation
**Script**: `validate_telomers.py`  
**Purpose**: Tests fluorotelomer detection on real PubChem data  
**Molecules**: PubChem fluorotelomer dataset  
**Output**: `data/telomer_validation_results.json`

Validates the detection of fluorotelomers (key industrial PFAS precursors) against real-world chemical database entries.

## Output Structure

### Directory Organization

```
benchmark/
├── data/               # JSON benchmark results
│   ├── pfas_enhanced_benchmark_*.json
│   ├── pfas_oecd_benchmark_*.json
│   ├── pfas_timing_benchmark_*.json
│   ├── pfas_complex_branched_benchmark_*.json
│   ├── highly_branched_pfas_test_*.json
│   └── telomer_validation_results.json
│
├── html/               # HTML reports
│   └── unified_pfas_benchmark_report_*.html
│
├── imgs/               # Generated plots and figures
│   ├── timing_*.png
│   ├── complex_benchmark_analysis_*.png
│   └── enhanced_*.png
│
├── reports/            # Analysis reports
│   └── (historical HTML reports)
│
└── review-app/         # Interactive review application
    ├── database/
    │   └── pfas_benchmark.db
    └── analysis_reports/
        ├── timing_analysis.json
        ├── complex_analysis.json
        ├── enhanced_analysis.json
        └── figures/
```

### Key Output Files

1. **Unified Report**: `html/unified_pfas_benchmark_report_*.html`
   - Combined report of all benchmarks
   - Open in any web browser
   - Contains summary statistics and visualizations

2. **Individual JSON Results**: `data/*.json`
   - Machine-readable benchmark results
   - Used by analysis scripts
   - Timestamped for version tracking

3. **Analysis Plots**: `imgs/*.png`, `imgs/*.svg`
   - Timing performance curves
   - Accuracy visualizations
   - Comparative analysis charts

4. **Database**: `review-app/database/pfas_benchmark.db`
   - SQLite database with all results
   - Enables interactive exploration
   - Powers the review web application

## Analysis Scripts

The benchmark suite includes specialized analysis scripts that run automatically:

### Timing Analysis
**Script**: `scripts/analyze_timing.py`  
**Input**: `data/pfas_timing_benchmark_*.json`  
**Output**: 
- `review-app/analysis_reports/timing_analysis.json`
- `imgs/timing_*.png` (performance curves)

Analyzes computational performance, identifies scaling behavior, and compares tools.

### Complex Branched Analysis
**Script**: `scripts/analyze_complex.py`  
**Input**: `data/pfas_complex_branched_benchmark_*.json`  
**Output**: 
- `review-app/analysis_reports/complex_analysis.json`
- `imgs/complex_benchmark_analysis_*.png`
- `complex_benchmark_report_*.html`

Analyzes detection accuracy on complex structures, accounting for perfluoro/polyfluoro chemistry distinctions.

### Enhanced Analysis
**Script**: `scripts/enhanced_analysis.py`  
**Input**: Both enhanced and OECD benchmark files  
**Output**: 
- `review-app/analysis_reports/enhanced_analysis.json`
- `imgs/enhanced_*.png`

Comprehensive analysis of functional group detection and OECD validation results.

### Analysis Scripts

After benchmarks complete, several analysis scripts process the results:

#### Timing Analysis
**Script**: `scripts/analyze_timing.py`  
**Input**: `data/pfas_timing_benchmark_*.json`  
**Output**: Interactive HTML plots and performance metrics

Analyzes timing performance comparing PFASGroups and PFAS-Atlas:
- Execution time vs. molecule size
- Scaling behavior analysis
- Performance comparison charts

#### Timing Models
**Script**: `scripts/analyze_timing_models.py`  
**Input**: `data/pfas_timing_benchmark_*.json`  
**Output**: Exponential fit models and predictions

Generates exponential models for timing data:
- Model: t = a × exp(α × n) where n is number of atoms
- Parameter estimation and goodness of fit
- Predictive models for different configurations

#### Definitions Benchmark Analysis
**Script**: `scripts/analyze_definitions_benchmark.py`  
**Input**: `data/pfas_definitions_benchmark_*.json`  
**Output**: Comprehensive HTML report with visualizations

Analyzes PFAS definition benchmark results (optional):
- Performance metrics per definition (OECD, EU, OPPT, UK, PFASTRUCTv5)
- Confusion matrices
- Agreement/disagreement patterns
- Interactive visualizations

Note: This analysis only runs if definitions benchmark data is available.

#### Complex Branched Analysis
**Script**: `scripts/analyze_complex.py`  
**Input**: `data/pfas_complex_branched_benchmark_*.json`  
**Output**: Accuracy metrics for complex structures

#### Enhanced Analysis
**Script**: `scripts/enhanced_analysis.py`  
**Input**: `data/pfas_enhanced_benchmark_*.json`, `data/pfas_oecd_benchmark_*.json`  
**Output**: Combined functional groups and OECD analysis

#### Comprehensive Statistics
**Script**: `scripts/analyze_benchmarks_simple.py`  
**Output**: LaTeX tables and benchmark summary JSON

Generates comprehensive statistics across all benchmarks:
- Exponential fit visualization
- Summary statistics
- LaTeX-formatted tables for publication
- Saved to `reports/benchmark_summary.json`

### Telomer Report
**Script**: `scripts/generate_telomer_report.py`  
**Input**: `data/telomer_validation_results.json`  
**Output**: `telomer_validation_report_*.html`

Detailed report on fluorotelomer detection performance.

## Review Application

The benchmark suite includes an interactive web application for exploring results:

### Starting the Review App

```bash
cd review-app
node server.js
```

Then open your browser to: **http://localhost:5000**

### Features

- **Molecule Browser**: Browse all tested molecules
- **Analysis Reports**: View timing, complex, and enhanced analysis
- **Manual Review**: Flag issues or add notes to specific molecules
- **Filtering**: Filter by PFAS groups, detection status, complexity
- **Visualization**: Interactive charts and molecular structure viewers

## Workflow Details

### Complete Execution Flow

1. **Benchmark Execution** (Sequential)
   - Functional Groups → OECD → Timing → Non-Fluorinated → Complex → Highly Branched → Telomer
   
2. **Report Generation**
   - Unified HTML report created from all benchmark results
   
3. **Analysis Execution**
   - Timing analysis (performance curves)
   - Timing models (exponential fit and predictions)
   - Definitions analysis (if definitions benchmark data available)
   - Complex branched analysis (accuracy assessment)
   - Enhanced analysis (functional groups + OECD)
   - Comprehensive statistics (LaTeX tables and summary)
   - Telomer report generation
   
4. **File Organization**
   - JSON files → `data/`
   - HTML reports → `html/` and `reports/`
   - Images → `imgs/`
   
5. **Database Import**
   - Backup existing database
   - Clear old data
   - Import new benchmark results
   - Calculate molecular formulas
   - Analysis reports available in Review App

### Error Handling

- Each benchmark must complete successfully before the next runs
- If a benchmark fails, the script exits with error code
- Analysis failures are non-critical (warnings shown, execution continues)
- Database import failure exits with error

## Troubleshooting

### Common Issues

**Import Errors:**
```bash
conda activate pfasatlas  # or your environment
pip install -e /path/to/PFASGroups
pip install -e /path/to/PFAS-atlas
```

**Missing Files:**
- Ensure you're in the `benchmark/` directory
- Check that `scripts/enhanced_pfas_benchmark.py` exists

**Database Errors:**
- Node.js must be installed
- Check `review-app/database/` directory exists
- Review backups in `review-app/database/*.backup_*`

**Permission Errors (Linux/macOS):**
```bash
chmod +x run_all_benchmarks.sh
```

**PowerShell Execution Policy (Windows):**
```powershell
Set-ExecutionPolicy -ExecutionPolicy RemoteSigned -Scope CurrentUser
```

## Performance Considerations

- **Full suite runtime**: ~5-15 minutes (depends on system)
- **Memory usage**: ~2-4 GB RAM
- **Disk space**: ~50-100 MB for results
- **CPU**: Benefits from multi-core processors

## Customization

### Running Individual Benchmarks

You can run individual benchmarks manually:

```bash
cd benchmark

# Functional Groups
echo "1" | python scripts/enhanced_pfas_benchmark.py

# Timing Benchmark
echo "3" | python scripts/enhanced_pfas_benchmark.py

# Highly Branched Test
python scripts/test_highly_branched.py

# Telomer Validation
python scripts/validate_telomers.py
```

### Modifying Test Parameters

Edit `scripts/enhanced_pfas_benchmark.py`:
- Change `molecules_per_test` parameter in benchmark functions
- Modify test molecule sets in benchmark definitions
- Adjust timing benchmark size range

### Custom Analysis

Create custom analysis scripts by loading the JSON results:
```python
import json

with open('data/pfas_enhanced_benchmark_20260201_120000.json', 'r') as f:
    data = json.load(f)
    
# Your custom analysis here
```

## Chemistry Notes

### Perfluoro vs. Polyfluoro

The benchmarks correctly distinguish between:
- **Perfluoro** (groups 49, 55, 57): No C-H bonds, fully fluorinated
- **Polyfluoro** (groups 51, 56, 58): Some C-H bonds present

Functional groups like -COOH, -SO₃H contain hydrogen, making molecules polyfluorinated even if the carbon chain is perfluorinated.

### Group Equivalences

The analysis scripts account for these chemistry-based equivalences:
- 49 ↔ 51 (perfluoroalkyl ↔ polyfluoroalkyl)
- 55 ↔ 56 (perfluoro cyclic ↔ polyfluoro cyclic)
- 57 ↔ 58 (perfluoroaryl ↔ polyfluoroaryl)

## Contributing

When modifying benchmarks:
1. Update expected groups based on actual chemistry
2. Run full benchmark suite to validate
3. Check analysis reports for accuracy
4. Update this documentation if adding new tests

## Support

For issues or questions:
- Check benchmark output logs for detailed error messages
- Review HTML reports for specific failures
- Use Review App to inspect individual molecule results
- Check `review-app/database/*.backup_*` for database history
