# Extending Timing Benchmark Results

## Overview

The timing benchmark now supports **incremental extension** - you can load previous timing results and add more molecules to extend your dataset without starting from scratch.

## Why Use This Feature?

- **Preserve Previous Work**: Don't lose hours of benchmark data when you need more samples
- **Incremental Growth**: Build large datasets gradually (e.g., 200 molecules at a time)
- **Grid-Based Sampling**: The system ensures maximum coverage of parameter space
- **Statistical Power**: Easily reach 1000+ molecules for robust analysis

## Usage Methods

### Method 1: Using the Dedicated Script (Recommended)

```bash
cd benchmark
./extend_timing_benchmark.sh [num_molecules] [iterations]
```

**Examples:**
```bash
# Add 2500 more molecules with 5 iterations (default)
./extend_timing_benchmark.sh

# Add 1000 molecules with 10 iterations
./extend_timing_benchmark.sh 1000 10

# Add 500 molecules with 5 iterations (faster testing)
./extend_timing_benchmark.sh 500 5
```

### Method 2: Using run_all_benchmarks.sh

```bash
cd benchmark
./run_all_benchmarks.sh --reuse-timing
```

This runs all benchmarks, but reuses previous timing results for the timing benchmark portion.

### Method 3: Manual Python Execution

```bash
cd benchmark/scripts
python enhanced_pfas_benchmark.py
# Choose option 3 (Timing Performance Benchmark)
# When asked "Reuse previous timing results and add more? (y/N):", enter 'y'
```

## How It Works

### 1. Automatic Detection
When you enable reuse mode, the system:
- Searches `benchmark/data/` for timing result files (`pfas_timing_benchmark_*.json`)
- Loads the most recent file automatically
- Reports how many previous results were found

### 2. Seamless Merging
The benchmark then:
- Starts with your previous results
- Continues the grid-based sampling where it left off
- Adds the requested number of new molecules
- Creates a **new combined JSON file** with all results

### 3. Grid-Based Coverage
The system uses a grid of:
- **Chain lengths**: 5, 10, 15, 20, ..., 200 carbons (24 bins)
- **Functional groups**: ~80 groups (IDs 29-114)
- **Total combinations**: ~1,920 unique grid points

Each run cycles through shuffled grid points, ensuring systematic coverage of the parameter space.

## Example Workflow

```bash
# Day 1: Initial benchmark (2500 molecules, 5 iterations)
cd benchmark
printf "3\nn\n" | python scripts/enhanced_pfas_benchmark.py
# Creates: data/pfas_timing_benchmark_20260206_143022.json (2500 molecules)
# Also creates: data/pfas_timing_excluded_20260206_143022.json (excluded molecules for review)

# Day 2: Extend with 1000 more
./extend_timing_benchmark.sh 1000
# Creates: data/pfas_timing_benchmark_20260207_091544.json (3500 molecules total)

# Day 3: Extend with 1500 more
./extend_timing_benchmark.sh 1500
# Creates: data/pfas_timing_benchmark_20260208_154312.json (5000 molecules total)

# Analyze the comprehensive dataset
python scripts/analyze_timing_with_complexity.py data/pfas_timing_benchmark_20260208_154312.json
```

## Output Files

Each run creates **two timestamped JSON files**:

### 1. Main Results File
`pfas_timing_benchmark_YYYYMMDD_HHMMSS.json`
- All successfully tested molecules
- Previous molecules (with their original IDs)
- Newly generated molecules (with incremented IDs)
- Full metadata and complexity metrics

### 2. Excluded Molecules File (if any)
`pfas_timing_excluded_YYYYMMDD_HHMMSS.json`
- Molecules where target group was not detected
- Saved for review to check for:
  - Bugs in molecular generation
  - Issues with SMARTS patterns
  - Edge cases in group detection
- Contains: SMILES, target group, detected groups, exclusion reason

**Important**: Review excluded molecules to ensure they represent expected behavior and not bugs in the generation or detection logic.

**Original files are preserved** - the system never modifies existing benchmark data.

## Integration with Database & Reports

The merged results work seamlessly with all analysis tools:

```bash
# Import merged data to database
cd review-app/scripts
node import-benchmark-data.js

# Generate enhanced HTML report
cd ../../scripts
python generate_enhanced_timing_report.py

# Analyze complexity correlations
python analyze_timing_with_complexity.py ../data/pfas_timing_benchmark_*.json
```

## Benefits

✅ **No Data Loss** - Previous results are always preserved  
✅ **Flexible Growth** - Add molecules in manageable batches  
✅ **Systematic Coverage** - Grid ensures comprehensive parameter space exploration  
✅ **Time Efficient** - Avoid re-running expensive benchmarks  
✅ **Statistical Power** - Build large datasets for robust conclusions  

## Technical Details

### File Naming
- Format: `pfas_timing_benchmark_YYYYMMDD_HHMMSS.json`
- Most recent file is automatically detected by modification time

### Molecule IDs
- IDs are auto-incremented: `molecule_id: 1, 2, 3, ...`
- When extending, new molecules continue from `last_id + 1`

### Grid State
- Grid is reshuffled each run for randomization
- Previous grid position is NOT preserved (ensures fresh randomization)

### Statistics
- Progress reports show: `Generated X/Y molecules` where Y = previous + new
- Final summary includes all molecules in the dataset

## Command Reference

| Command | Purpose |
|---------|---------|
| `./extend_timing_benchmark.sh` | Add 200 molecules to latest results |
| `./extend_timing_benchmark.sh 500 15` | Add 500 molecules, 15 iterations |
| `./run_all_benchmarks.sh --reuse-timing` | Run all benchmarks, extend timing |
| `python enhanced_pfas_benchmark.py` | Manual execution with prompts |

## Tips

1. **Start Small**: Begin with 100-200 molecules to test your setup
2. **Be Patient**: Large molecules (150+ carbons) take longer to generate
3. **Monitor Progress**: Watch for "✅ Generated X/Y molecules" updates every 50 molecules
4. **Check Coverage**: Use analysis scripts to verify grid coverage statistics
5. **Backup Data**: Original files in `data/` are your backup - keep them safe!

## Troubleshooting

### "No previous timing results found"
- Check that `benchmark/data/pfas_timing_benchmark_*.json` files exist
- Run initial benchmark first without reuse option

### "Timing benchmark extension failed"
- Verify Python environment is activated (`conda activate pfasatlas`)
- Check that PFASGroups and dependencies are installed
- Review error messages for specific issues

### Unexpected molecule count
- New file contains: previous count + newly requested count
- Failed generations don't count toward the total
- Check console output for exclusion statistics
