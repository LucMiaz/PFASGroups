# PFAS Algorithm Automated Test Report Generator

This directory contains scripts to automatically test the PFAS group identification algorithm and generate comprehensive reports showing accuracy, specificity, and detailed problem analysis.

## 🚀 Quick Start

To generate a test report, simply run:

```bash
python generate_test_report.py
```

This will automatically:
1. **Try to run fresh tests** using the test framework in `test_examples.py`
2. **Fall back to existing data** if there are dependency issues (like RDKit problems)
3. **Generate a comprehensive HTML report** with detailed analysis and molecular examples

## 📊 What You Get

The generated report includes:
- **Overall Performance Metrics**: Detection rate, specificity rate, average groups per test
- **Algorithm Status Assessment**: EXCELLENT/GOOD/NEEDS IMPROVEMENT/REQUIRES ATTENTION
- **Problem Case Analysis**: Detailed tables showing:
  - False negatives (expected groups not detected)
  - False positives (unexpected groups detected)
  - Low specificity cases (too many groups detected)
  - Severe specificity issues (>5 groups detected)
- **Molecular Examples**: Actual SMILES strings for each problem case
- **Group-Level Performance**: Which specific PFAS groups are most problematic
- **Actionable Recommendations**: Specific steps to improve algorithm performance

## 📁 Scripts Overview

### `generate_test_report.py` (Main Entry Point)
- **Purpose**: Smart wrapper that tries fresh testing but falls back gracefully
- **Use**: `python generate_test_report.py`
- **Output**: Comprehensive HTML report
- **Handles**: RDKit dependency issues, test failures, etc.

### `automated_test_runner.py` (Fresh Testing)
- **Purpose**: Runs fresh PFAS group detection tests using `test_examples.py`
- **Features**: 
  - Generates new test molecules
  - Runs specificity, OECD, and generic tests
  - Analyzes results in real-time
- **Requirements**: Working RDKit installation
- **Output**: `pfas_algorithm_test_report.html`

### `standalone_report_generator.py` (Existing Data Analysis)
- **Purpose**: Analyzes existing test result CSV files
- **Features**:
  - No RDKit dependency required
  - Works with any existing test data
  - Handles various CSV formats
- **Input**: `*test_results.csv` files
- **Output**: `pfas_algorithm_existing_data_report.html`

## 🔧 Troubleshooting

### "No test data found"
- Ensure you have test result CSV files in the current directory or parent directories
- Try running tests manually first: `cd ../PFASgroups/tests && python test_examples.py specificity`

### RDKit Import Errors
- The scripts automatically fall back to existing data analysis
- For fresh testing, ensure RDKit is properly installed: `conda install rdkit`

### Permission Errors
- Ensure you have write permissions in the current directory
- Try running from a different location

## 📈 Understanding the Results

### Performance Metrics
- **Detection Rate**: % of expected groups that were correctly identified
- **Specificity Rate**: % of tests where only relevant groups were detected
- **Average Groups per Test**: Lower is better (indicates specificity)

### Algorithm Status
- **EXCELLENT**: >95% detection, >90% specificity, <2 avg groups
- **GOOD**: >85% detection, >80% specificity, <3 avg groups
- **NEEDS IMPROVEMENT**: >70% detection and specificity
- **REQUIRES ATTENTION**: Below 70% on key metrics

### Problem Cases
- **False Negatives**: Algorithm missed expected PFAS groups
- **False Positives**: Algorithm detected unexpected groups
- **Low Specificity**: Too many groups detected (>3)
- **Severe Specificity**: Very poor specificity (>5 groups)

## 🔄 Regular Usage

For ongoing algorithm development:

1. **After making changes to SMARTS patterns**: `python generate_test_report.py`
2. **To track performance over time**: Save reports with timestamps
3. **To focus on specific issues**: Look at the "Most Frequently Missing/Over-detected Groups" sections
4. **For detailed debugging**: Check the molecular examples in problem case tables

## 📝 Output Files

The scripts generate several files:
- `*_report.html`: Main comprehensive report (open in web browser)
- `*test_results.csv`: Raw test data for further analysis
- `test_summary.json`: Machine-readable summary
- `false_negatives.csv`, `specificity_issues.csv`: Specific problem case data

## 🎯 Customization

To modify the analysis:
- Edit `standalone_report_generator.py` for different metrics or visualizations
- Modify `automated_test_runner.py` to change test parameters
- Update the HTML templates in either script for different report styles

## 🆘 Support

If you encounter issues:
1. Check the console output for specific error messages
2. Ensure you're in the correct directory (should contain the scripts)
3. Verify that test data files exist and are readable
4. Try running the individual scripts manually for more detailed error output

---

## Example Output

A typical successful run looks like:
```
🎯 Detection Rate: 99.3%
🎯 Specificity Rate: 100.0%
🎯 Average Groups per Test: 2.3
📊 Total Tests: 430

🔍 PROBLEM CASES:
  False Negatives: 3
  False Positives: 0
  Low Specificity: 10

📄 Report saved as: pfas_algorithm_existing_data_report.html
✅ Analysis completed successfully!
```

The generated HTML report provides detailed molecular-level analysis with specific SMILES examples for each problem case, making it easy to identify exactly which patterns need improvement.