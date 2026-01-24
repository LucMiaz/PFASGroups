# PFAS Definitions Benchmark - Quick Reference

## 🚀 Quick Start

```powershell
cd c:\Users\luc\git\PFASGroups\benchmark
.\run_definitions_benchmark.ps1
```

## 📋 What Gets Tested

### 5 PFAS Definitions
1. **OECD** - CF₃/CF₂ not bonded to H/Cl/Br/I
2. **EU** - OECD + more exclusions
3. **OPPT** - 3 inclusion criteria
4. **UK** - CF₃/CF₂ not bonded to Cl/Br/I
5. **PFASTRUCTv5** - Patterns OR F-ratio ≥0.3

### Test Categories
- ✅ **True Positives** (~30 known PFAS)
- ❌ **True Negatives** (~15 non-PFAS)
- ⚖️ **Edge Cases** (~15 borderline)
- 🎯 **Definition-Specific** (~15 targeted)
- 🔗 **Concordance** (~10 reference)
- ⏱️ **Performance** (~100 random)

## 📊 Key Metrics

### Sensitivity (True Positive Rate)
- **What**: % of PFAS correctly identified
- **Target**: >90%
- **Low means**: Missing PFAS compounds

### Specificity (True Negative Rate)
- **What**: % of non-PFAS correctly rejected
- **Target**: 100%
- **Low means**: Too many false positives

### Jaccard Similarity
- **What**: Agreement between definitions
- **Range**: 0 (no overlap) to 1 (perfect agreement)
- **High means**: Definitions are similar

## 📁 Output Files

```
benchmark/data/
├── pfas_definitions_benchmark_YYYYMMDD_HHMMSS.json  # Raw results
├── definitions_benchmark_report_YYYYMMDD_HHMMSS.html # HTML report
└── figures_definitions/                              # Visualizations
    ├── sensitivity_specificity.html
    ├── detection_heatmap.html
    ├── concordance_matrix.png
    ├── edge_case_disagreements.html
    └── performance_by_size.html
```

## 🎯 Interpreting Results

### Good Performance
```
✅ Sensitivity: 95%+ (catches most PFAS)
✅ Specificity: 100% (no false positives)
✅ Fast: <5ms per molecule
```

### Warning Signs
```
⚠️ Sensitivity: <80% (missing PFAS)
⚠️ Specificity: <95% (false positives)
⚠️ Slow: >20ms per molecule
```

## 🔧 Common Tasks

### Run Full Benchmark
```powershell
python benchmark_pfas_definitions.py
python analyze_definitions_benchmark.py
```

### Test Single Molecule
```python
from benchmark_pfas_definitions import PFASDefinitionBenchmark

benchmark = PFASDefinitionBenchmark()
result = benchmark.test_single_molecule('FC(F)(F)C(F)(F)F')
print(f"Detected by: {result['detected_ids']}")
```

### Compare Two Runs
```powershell
# Run baseline
python benchmark_pfas_definitions.py
# Note the timestamp

# Make changes to definitions

# Run new benchmark
python benchmark_pfas_definitions.py

# Compare HTML reports side-by-side
```

## 📈 Visualization Guide

### Sensitivity vs Specificity Plot
- **Top-right corner**: Ideal (high both)
- **Top-left**: Too inclusive (false positives)
- **Bottom-right**: Too restrictive (misses PFAS)

### Detection Heatmap
- **Green**: Good detection (>90%)
- **Yellow**: Moderate (70-90%)
- **Red**: Poor (<70%)

### Concordance Matrix
- **Diagonal**: Always 1.0 (self-agreement)
- **Off-diagonal**: Agreement between different definitions
- **Dark red**: High agreement
- **Light yellow**: Low agreement

## ⚡ Performance Targets

| Metric | Target | Acceptable | Needs Work |
|--------|--------|------------|------------|
| Sensitivity | >95% | 85-95% | <85% |
| Specificity | 100% | >95% | <95% |
| Execution Time | <5ms | 5-15ms | >15ms |
| Agreement (Jaccard) | N/A | >0.7 | <0.5 |

## 🐛 Troubleshooting

### No results generated
```powershell
# Check Python environment
python --version

# Check imports
python -c "from PFASgroups.core import parse_mol; print('OK')"

# Check RDKit
python -c "from rdkit import Chem; print('OK')"
```

### Slow performance
```python
# Reduce test size
benchmark.run_performance_benchmark(num_molecules=50)
```

### Visualization errors
```powershell
pip install --upgrade plotly kaleido matplotlib seaborn
```

## 📚 More Information

See full guide: `PFAS_DEFINITIONS_BENCHMARKING_GUIDE.md`

## 🔄 Update Workflow

1. **Before changes**: Run benchmark (baseline)
2. **Make changes**: Modify definitions or code
3. **After changes**: Run benchmark (updated)
4. **Compare**: Review HTML reports
5. **Document**: Note improvements/regressions

## 📞 Support

- Full documentation: `PFAS_DEFINITIONS_BENCHMARKING_GUIDE.md`
- Code: `benchmark_pfas_definitions.py`
- Analysis: `analyze_definitions_benchmark.py`

---

**Quick Command Reference**

```powershell
# Full benchmark + analysis
.\run_definitions_benchmark.ps1

# Just benchmark
python benchmark_pfas_definitions.py

# Just analysis (uses latest results)
python analyze_definitions_benchmark.py

# Analyze specific file
python analyze_definitions_benchmark.py data/pfas_definitions_benchmark_20260123_120000.json

# View report
start data/definitions_benchmark_report_YYYYMMDD_HHMMSS.html
```
