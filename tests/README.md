# PFASGroups Testing Framework

This directory contains the testing infrastructure for validating PFAS groups and definitions against their test metadata.

## Overview

The testing framework validates that:
1. **PFAS Groups** correctly identify positive examples from test metadata
2. **PFAS Definitions** correctly classify molecules as PFAS/non-PFAS

## Test Methods

### PFASGroup.test()

Added to `PFASGroupModel.py`, this method validates a PFAS group against test examples:

```python
from PFASgroups import PFASGroup

# Load a group
group = PFASGroup(id=1, name="Perfluoroalkyl carboxylic acids", ...)

# Run tests
result = group.test()
print(f"Passed: {result['passed']}")
print(f"Failures: {result['failures']}")
```

**Returns:**
```python
{
    'passed': bool,                    # Whether all tests passed
    'total_tests': int,                # Total number of tests run
    'failures': [...],                 # List of failure details
    'category': 'OECD'|'generic'|'telomer'
}
```

### PFASDefinition.test()

Added to `PFASDefinitionModel.py`, this method validates a PFAS definition:

```python
from PFASgroups import PFASDefinition

# Load a definition
definition = PFASDefinition(
    id=1, 
    name="OECD Definition",
    smarts=["[#6X4!$([#6H1])!$([#6][#17,#35,#53])](F)F"],
    ...
)

# Run tests
result = definition.test()
print(f"True Positives: {result['stats']['true_positives']}")
print(f"True Negatives: {result['stats']['true_negatives']}")
```

**Returns:**
```python
{
    'passed': bool,
    'total_tests': int,
    'failures': [...],
    'category': 'definition',
    'stats': {
        'true_positives': int,   # Should match, did match
        'true_negatives': int,   # Should not match, did not match
        'false_positives': int,  # Should not match, but matched
        'false_negatives': int   # Should match, but did not match
    }
}
```

## Test Script Usage

### Basic Usage

```bash
# Test everything (all groups and definitions)
python tests/run_groups_definitions_tests.py

# Test only groups
python tests/run_groups_definitions_tests.py --groups-only

# Test only definitions
python tests/run_groups_definitions_tests.py --definitions-only
```

### Using pytest

For pytest users, there's a pytest-compatible test suite:

```bash
# Run all tests with pytest
pytest tests/test_pytest.py -v

# Run only groups tests
pytest tests/test_pytest.py -m groups -v

# Run only definitions tests
pytest tests/test_pytest.py -m definitions -v

# Run smoke tests (quick validation)
pytest tests/test_pytest.py -m smoke -v

# Test specific group
pytest tests/test_pytest.py::test_individual_group[1] -v

# Test specific definition
pytest tests/test_pytest.py::test_individual_definition[1] -v
```

### Verbose Output

```bash
# Show detailed output for each test
python tests/run_groups_definitions_tests.py -v

# Test specific group with verbose output
python tests/run_groups_definitions_tests.py --group-id 1 -v
```

### Specific Tests

```bash
# Test specific group by ID
python tests/run_groups_definitions_tests.py --group-id 1

# Test specific definition by ID
python tests/run_groups_definitions_tests.py --definition-id 1
```

### Output Control

```bash
# Specify custom failure report file
python tests/run_groups_definitions_tests.py --output my_failures.txt
```

## Test Metadata Structure

### PFAS Groups (PFAS_groups_smarts.json)

Each group should have a `test` key:

```json
{
  "id": 1,
  "name": "Perfluoroalkyl carboxylic acids",
  "smarts": {...},
  "test": {
    "category": "OECD",
    "examples": [
      "OC(=O)C(F)(F)C(F)(F)C(F)(F)F",
      "OC(=O)C(F)(F)C(F)(F)C(F)(F)C(F)(F)F"
    ]
  }
}
```

**Categories:**
- `OECD`: Curated examples from OECD benchmark
- `generic`: Generic functional groups with positive/negative examples
- `telomer`: Telomer groups with generation patterns

### PFAS Definitions (PFAS_definitions_smarts.json)

Each definition should have a `test` key:

```json
{
  "id": 1,
  "name": "OECD Definition",
  "smarts": ["[#6X4!$([#6H1])!$([#6][#17,#35,#53])](F)F"],
  "test": {
    "category": "definition",
    "examples": {
      "true_positives": [
        {"smiles": "C(F)(F)(F)C(F)(F)F", "category": "PFAS"}
      ],
      "true_negatives": [
        {"smiles": "CCCCCC", "category": "non-PFAS"}
      ],
      "false_positives": [],  // Known incorrect matches (for documentation)
      "false_negatives": []   // Known missed matches (for documentation)
    }
  }
}
```

## Exit Codes

- `0`: All tests passed
- `1`: One or more tests failed

## Output Files

### test_failures_report.txt

Generated when tests fail, contains:
- Detailed failure information for each failing group/definition
- SMILES strings that failed
- Expected vs actual results
- Error messages

Example:
```
================================================================================
GROUP 1: Perfluoroalkyl carboxylic acids
================================================================================

Total Tests: 4
Failures: 1
Category: OECD

Failure Details:

1. SMILES: OC(=O)C(F)(F)C(F)(F)C(F)(F)F
   Expected: True
   Got: False
   Error: Group should match but did not
```

## Integration with CI/CD

Add to your CI pipeline:

```yaml
# .github/workflows/test.yml
- name: Run PFASGroups Tests
  run: |
    python tests/run_groups_definitions_tests.py
```

## Test Metadata Generation

Test metadata can be added/updated using:

```bash
# Add test metadata from benchmarks
python benchmark/scripts/add_test_metadata.py
```

This script:
1. Loads examples from OECD benchmark
2. Loads telomer validation results
3. Adds `test` keys to all groups and definitions
4. Categorizes groups (OECD/telomer/generic)

## Test Results Summary

Current status (as of last run):
- **PFAS Groups**: 113/113 passed (100%)
- **PFAS Definitions**: Test metadata needs review

## Adding New Tests

### For a New PFAS Group

1. Add group to `PFAS_groups_smarts.json`
2. Add test metadata:
   ```json
   "test": {
     "category": "generic",
     "examples": [
       "SMILES1",  // Should match
       "SMILES2"   // Should match
     ]
   }
   ```
3. Run: `python tests/run_groups_definitions_tests.py --group-id <NEW_ID> -v`

### For a New PFAS Definition

1. Add definition to `PFAS_definitions_smarts.json`
2. Add test metadata with TP/TN/FP/FN examples
3. Run: `python tests/run_groups_definitions_tests.py --definition-id <NEW_ID> -v`

## Troubleshooting

### ImportError: No module named 'PFASgroups'

Make sure you're running from the repository root:
```bash
cd /path/to/PFASGroups
python tests/run_groups_definitions_tests.py
```

### "No test data found"

The group/definition is missing a `test` key in the JSON file. Run:
```bash
python benchmark/scripts/add_test_metadata.py
```

### Tests failing unexpectedly

1. Check test metadata is correct
2. Run with `-v` flag for detailed output
3. Verify SMILES strings are valid
4. Check if pattern/constraints are correct

## Files

- `run_groups_definitions_tests.py`: Main test script
- `__init__.py`: Test package initialization
- `test_failures_report.txt`: Generated failure report (gitignored)

## Dependencies

- RDKit: Molecule handling
- NetworkX: Graph operations (via ComponentsSolver)
- Python 3.9+: Core language features
