 # HalogenGroups Tests

This folder contains pytest-compatible tests for HalogenGroups parsing, component metrics,
and group/definition metadata validation.

## Quick Start

From the repository root:

```bash
pytest tests -v
```

Run a focused subset:

```bash
# Groups/definitions metadata checks
pytest tests/test_pytest.py -m groups -v
pytest tests/test_pytest.py -m definitions -v

# Linker SMARTS behavior
pytest tests/test_linker_smarts.py -v

# Component metrics and fractions
pytest tests/test_component_fractions.py -v
pytest tests/test_metrics.py -v
```

## What The Tests Check

### Groups and Definitions
- `tests/test_pytest.py` and `tests/test_PFASgroups_smarts.py` validate PFAS group and
	definition metadata using `PFASGroup.test()` and `PFASDefinition.test()`.
- `tests/test_run_groups_definitions_tests.py` includes pytest tests and a standalone
	runner for detailed failure reporting.

### Component Metrics and Fractions
- `tests/test_database_integration.py` verifies parse structure and presence of summary
	and component metrics (ready for database storage).
- `tests/test_component_fractions.py` and `tests/test_component_ratios.py` validate
	fraction calculations and bounds.
- `tests/test_metrics.py`, `tests/test_metrics_detailed.py`, and
	`tests/test_comprehensive_metrics.py` check that core graph metrics fields are present.

### Linker SMARTS and Telomers
- `tests/test_linker_smarts.py` validates linker-only behavior for fluorotelomer alcohols.
- `tests/test_telomer_validation.py` loads the telomer SDF if available; it skips if the
	file is missing.

### SMARTS Precompute Sanity
- `tests/test_smarts_atom_count.py` and `tests/test_smarts_precompute.py` ensure SMARTS
	precompute fields are populated and parsable.

### Results Model Helpers
- `tests/test_results_model.py` checks ResultsModel helper utilities and JSON safety.

## Standalone Validation Script

For a verbose, report-style run against the metadata files:

```bash
python tests/run_groups_definitions_tests.py -v
```

Additional options:

```bash
python tests/run_groups_definitions_tests.py --groups-only
python tests/run_groups_definitions_tests.py --definitions-only
python tests/run_groups_definitions_tests.py --group-id 1
python tests/run_groups_definitions_tests.py --definition-id 1
```

## Data Dependencies

- `tests/test_telomer_validation.py` expects `benchmark/data/PubChem_fluorotelomers.sdf`.
	The test will skip automatically if the file is missing.

## Notes

- All `tests/test_*.py` files are pytest-compatible and should not execute work at import
	time; use pytest to run them.
- If you see import errors for `PFASgroups`, ensure you run pytest from the repository
	root.
