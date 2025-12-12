# PFAS Classification Benchmark: PFASGroups vs PFAS-atlas

This benchmark compares the classification results between PFASGroups and PFAS-atlas classification systems using molecules from the PFASGroups specificity test results.

## Overview

The benchmark analyzes:
- **PFASGroups**: Rule-based PFAS functional group classification system
- **PFAS-atlas**: Machine learning-based PFAS classification system

## Files

- `pfas_atlas_benchmark.py` - Main benchmark script
- `run_benchmark.sh` - Shell script to run the benchmark with proper environment
- `requirements_benchmark.txt` - Python dependencies for the benchmark
- `README.md` - This documentation

## Setup

### Prerequisites

1. **Environments**: Ensure you have the `pfasatlas` mamba environment set up with PFAS-atlas dependencies
2. **Data**: The specificity test results CSV should be available at:
   `/home/luc/git/PFASGroups/PFASgroups/tests/results/specificity_test_results.csv`
3. **PFAS-atlas**: The PFAS-atlas repository should be available at:
   `/home/luc/git/PFAS-atlas`

### Installation

```bash
# Navigate to benchmark directory
cd /home/luc/git/PFASGroups/benchmark

# Activate the pfasatlas environment
mamba activate pfasatlas

# Install any missing requirements
pip install -r requirements_benchmark.txt
```

## Usage

### Option 1: Using the shell script (Recommended)
```bash
./run_benchmark.sh
```

### Option 2: Manual execution
```bash
# Activate environment
mamba activate pfasatlas

# Run benchmark
python pfas_atlas_benchmark.py
```

## Output

The benchmark generates three output files:

1. **`pfas_atlas_benchmark_report.md`** - Detailed human-readable report with:
   - Summary statistics
   - Classification agreement analysis
   - Top classifications from both systems

2. **`pfas_atlas_benchmark_results.csv`** - Complete results table with:
   - Original molecule data
   - PFASGroups classification results
   - PFAS-atlas classification results
   - Parsed and processed classification data

3. **`pfas_atlas_benchmark_summary.json`** - Summary statistics in JSON format:
   - Coverage percentages
   - Agreement rates
   - Classification counts

## Interpretation

### Key Metrics

- **Coverage**: Percentage of molecules each system was able to classify
- **Agreement Rate**: How often both systems agree on classifying the same molecules
- **Classification Distribution**: Most common functional groups/classes identified

### Expected Results

The benchmark will show:
- How many molecules each system can classify (coverage)
- Where the systems agree or disagree
- Which functional groups/classes are most commonly identified
- Potential gaps or biases in either classification system

## Troubleshooting

### Environment Issues
```bash
# If pfasatlas environment doesn't exist, create it
mamba create -n pfasatlas python=3.9
mamba activate pfasatlas

# Install PFAS-atlas dependencies (check their README)
pip install torch rdkit-pypi mhfp
```

### Import Errors
- Ensure PFAS-atlas repository is properly installed/accessible
- Check that all required Python packages are installed in the pfasatlas environment
- Verify file paths are correct

### Data Issues
- Ensure the specificity test results CSV file exists and is readable
- Check that SMILES strings in the CSV are valid

## Extending the Benchmark

The benchmark can be extended to:
- Include additional test datasets
- Compare with other PFAS classification systems
- Analyze specific functional group performance
- Add statistical significance testing
- Include visualization of results

## Notes

- The benchmark handles edge cases like invalid SMILES, classification failures, etc.
- Progress is displayed during PFAS-atlas classification (every 100 molecules)
- The mapping of PFASGroups IDs to functional group names may need updating based on the latest PFASGroups definitions