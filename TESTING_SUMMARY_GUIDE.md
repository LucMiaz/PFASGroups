# Enhanced PFAS Groups Testing with Summary Generation

The `test_examples.py` file has been enhanced to automatically generate comprehensive summary reports when pytest is called. Here's how to use the new functionality:

## Features Added

### 1. Automatic Summary Generation
- **JSON Summary Report**: Comprehensive metrics saved to `test_summary_report.json`
- **CSV Detail Files**: Individual test results saved to separate CSV files
- **Console Summary**: Human-readable summary printed after tests complete

### 2. Metrics Collected
- **Overall Accuracy**: Percentage of correctly detected PFAS groups
- **Specificity Rate**: Percentage of tests with acceptable specificity (no false positives)
- **Detection Rate**: Percentage of expected groups that were detected
- **Group-wise Performance**: Individual performance for each PFAS group
- **Test Duration**: Time taken to complete all tests

### 3. Output Files

When tests are run, the following files are generated:

- `test_summary_report.json` - Main summary with all metrics
- `oecd_test_results.csv` - Detailed results for OECD PFAS group tests
- `generic_test_results.csv` - Detailed results for generic PFAS group tests
- `specificity_test_results.csv` - Detailed specificity test results

## Usage

### Option 1: Using pytest (Recommended)

```bash
# Run all tests with summary generation
pytest test_examples.py -v

# Run specific test methods
pytest test_examples.py::TestPFASGroups::test_oecd_pfas_groups -v
pytest test_examples.py::TestPFASGroups::test_specificity -v
```

### Option 2: Manual Test Execution

```bash
# Run all tests manually (without pytest)
python test_examples.py full

# Run specific test types
python test_examples.py specificity
python test_examples.py quick
```

### Option 3: Using the Test Runner Script

```bash
# Run with pytest
python test_runner_example.py pytest

# Run manually
python test_runner_example.py manual

# Display existing summary
python test_runner_example.py summary
```

## Summary Report Structure

The JSON summary report contains:

```json
{
  "test_metadata": {
    "timestamp": "2025-10-05T...",
    "test_duration_seconds": 45.2,
    "python_version": "3.x.x",
    "pytest_available": true
  },
  "oecd_test_results": {
    "total_tests": 150,
    "successful_detections": 135,
    "overall_detection_rate": 0.9,
    "group_wise_results": [...],
    "worst_performing_groups": [...],
    "best_performing_groups": [...]
  },
  "generic_test_results": { ... },
  "specificity_test_results": {
    "detection_rate": 0.85,
    "specificity_rate": 0.72,
    "average_detected_groups_per_test": 1.4
  },
  "overall_summary": {
    "overall_accuracy": 0.87,
    "test_status": "PASSED"
  }
}
```

## Integration with CI/CD

The enhanced test suite is designed for CI/CD integration:

1. **Exit Codes**: Tests return appropriate exit codes for automation
2. **Machine-readable Output**: JSON format for parsing by other tools
3. **Detailed Logging**: Comprehensive logs for debugging failures
4. **Thresholds**: Configurable pass/fail thresholds for different metrics

### Example CI Configuration (GitHub Actions)

```yaml
name: PFAS Groups Tests
on: [push, pull_request]

jobs:
  test:
    runs-on: ubuntu-latest
    steps:
    - uses: actions/checkout@v2
    - name: Set up Python
      uses: actions/setup-python@v2
      with:
        python-version: 3.8
    - name: Install dependencies
      run: |
        pip install rdkit-pypi pandas pytest tqdm networkx
    - name: Run PFAS Groups tests
      run: |
        cd PFASgroups/tests
        pytest test_examples.py -v
    - name: Upload test results
      uses: actions/upload-artifact@v2
      with:
        name: test-results
        path: |
          PFASgroups/tests/test_summary_report.json
          PFASgroups/tests/*_test_results.csv
```

## Understanding the Metrics

### Detection Rate
- **High (>90%)**: Excellent - Algorithm correctly identifies most PFAS groups
- **Medium (70-90%)**: Good - Some groups may be missed but generally reliable
- **Low (<70%)**: Poor - Significant improvements needed

### Specificity Rate
- **High (>80%)**: Excellent - Few false positives, high precision
- **Medium (60-80%)**: Acceptable - Some false positives but generally specific
- **Low (<60%)**: Poor - Too many false positives, low precision

### Overall Accuracy
- Combines detection and specificity metrics
- **Target**: >85% for production use
- **Minimum**: >70% for development/testing

## Troubleshooting

### Common Issues

1. **Import Errors**: Ensure all required packages are installed
   ```bash
   pip install rdkit-pypi pandas pytest tqdm networkx
   ```

2. **Missing Data Files**: Ensure PFAS group definition files are present
   - `data/PFAS_groups_smarts.json`
   - `tests/specificity_test_groups.json`

3. **Memory Issues**: For large test runs, consider:
   - Reducing test set sizes
   - Running tests in smaller batches
   - Increasing available memory

### Performance Optimization

- **Parallel Testing**: Use pytest-xdist for parallel execution
- **Test Selection**: Run specific test groups during development
- **Caching**: Results can be cached for repeated analysis

## Contributing

When adding new PFAS groups or test cases:

1. Update the group definitions in the appropriate files
2. Add corresponding test cases to the test matrices
3. Run the full test suite to ensure compatibility
4. Review the generated summary for any performance regressions

The enhanced testing framework provides comprehensive feedback to help maintain and improve the PFAS group classification algorithm.