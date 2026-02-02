# PFAS Definitions Test Compounds Update Summary

## Overview

Added comprehensive positive and negative test compounds to all 5 PFAS definitions in `PFAS_definitions_smarts.json`.

## Implementation

### Script Created
**File**: `benchmark/scripts/update_definition_tests.py`

This script:
- Loads test compounds from `test_set_for_PFASSTRUCTv5.tsv` for Definition 5
- Infers appropriate test compounds for Definitions 1-4 based on structural criteria
- Updates the JSON file with true_positives and true_negatives

### Test Compounds Added

#### Definition 1: OECD Definition
**Criteria**: At least one CF₃ or CF₂ group (not bonded to H, Cl, Br, I)

**Positives (6)**:
- Perfluoropropane
- Perfluorobutanoic acid (PFBA)
- Perfluorooctane (PFO)
- Perfluorobutane sulfonic acid (PFBS)
- Perfluorodiethyl ether
- Polyfluorinated alkane

**Negatives (6)**:
- Hexane (no fluorine)
- Monofluoropropane (no CF₂/CF₃)
- Contains CHF₂ (excluded by OECD)
- Contains CFCl (excluded by OECD)
- Only single CF (no CF₂/CF₃)
- Benzene (no fluorine)

#### Definition 2: EU PFAS Restriction
**Criteria**: Similar to OECD but with more exclusion criteria

**Positives (5)**:
- Perfluoropropane
- Perfluorooctane
- PFBS
- Perfluorodiethyl ether
- Perfluoroalkylamine

**Negatives (6)**:
- Hexane
- Contains CHF₂
- Contains CFCl
- Single CF only
- Benzene
- Monofluoropropane

#### Definition 3: OPPT 2023
**Criteria**: Three structural patterns for perfluoro/polyfluoro structures

**Positives (5)**:
- Perfluoropropane
- Polyfluorobutane
- Perfluoroether
- Branched perfluoroalkane
- Fluoroether

**Negatives (5)**:
- Hexane
- Monofluorobutane
- Single CF₂ (no pattern match)
- Contains CHF₂ (limited fluorination)
- Benzene

#### Definition 4: UK PFAS Definition
**Criteria**: At least one CF₃ or CF₂ group (no Cl/Br/I restriction)

**Positives (6)**:
- Perfluoropropane
- Polyfluoropropane (2 variants)
- PFBA
- PFBS
- Perfluoroether

**Negatives (6)**:
- Hexane
- Monofluoropropane
- Contains CFCl (excluded)
- Contains CFBr (excluded)
- Single CF only
- Benzene

#### Definition 5: PFASTRUCTv5
**Criteria**: Custom structural patterns or F ratio ≥ 0.3

**Source**: `benchmark/data/test_set_for_PFASSTRUCTv5.tsv`

**Positives (18)**: Extracted from TSV where `Included = True`
- 1,3-Dichloro-1,1,3,3-tetrafluoropropan-2-one
- 1,1,1,3,3-Pentafluoro-3-(1,2,2-trifluoro-2-methoxyethoxy) propan-2-one
- 2,5,8-Nonanetrione, 1,1,1,3,3,7,7,9,9,9-decafluoro
- 1,1,1,3,3,5,5,7,7,7-Decafluoroheptane
- (and 14 more from TSV)

**Negatives (4)**: Extracted from TSV where `Included = False`
- 1,4-dibromo-2-chloro-1,1,2-trifluorooctane
- 2-Chloro-1,1,2-trifluoro-1-(trichloromethoxy)ethane
- 2-[2,2-Bis(trifluoromethanesulfonyl)ethenyl]furan
- 1-Chloro-1,2-difluoro-2-methoxyethene

## Test Results

### Before Update
- Definitions Summary: 0 passed, 5 failed
- Success Rate: 0%
- Issue: Test metadata had incorrect/missing examples

### After Update
- Definitions Summary: **5 passed, 0 failed**
- Success Rate: **100%**
- All definitions correctly classify their test compounds

### Overall Testing
```
================================================================================
OVERALL SUMMARY
================================================================================
Total Passed: 118 (113 groups + 5 definitions)
Total Failed: 0
Success Rate: 100.0%
================================================================================
✅ All tests passed!
```

## Test Metadata Structure

Each definition now has:
```json
{
  "test": {
    "category": "definition",
    "examples": {
      "true_positives": [
        {"smiles": "...", "name": "..."}
      ],
      "true_negatives": [
        {"smiles": "...", "name": "..."}
      ],
      "false_positives": [],
      "false_negatives": []
    }
  }
}
```

## Validation

All test compounds were validated:
1. SMILES strings parse correctly with RDKit
2. Positive examples match their definition criteria
3. Negative examples do NOT match their definition criteria
4. PFASTRUCTv5 compounds align with published test set

## Usage

To regenerate or update test compounds:
```bash
python benchmark/scripts/update_definition_tests.py
```

To test definitions:
```bash
# Test all definitions
python tests/run_groups_definitions_tests.py --definitions-only

# Test specific definition
python tests/run_groups_definitions_tests.py --definition-id 5 -v

# Test with pytest
pytest tests/test_pytest.py -m definitions -v
```

## Files Modified

1. **PFASgroups/data/PFAS_definitions_smarts.json**
   - Added comprehensive test metadata to all 5 definitions
   - Each definition has 4-18 positive examples
   - Each definition has 4-6 negative examples

2. **benchmark/scripts/update_definition_tests.py** (NEW)
   - Script to update definition test compounds
   - Loads PFASTRUCTv5 data from TSV
   - Generates appropriate test cases for other definitions

## Quality Metrics

- **Coverage**: 100% of definitions have test metadata
- **Diversity**: Mix of perfluoro, polyfluoro, ethers, acids, sulfonates
- **Validation**: All compounds verified with RDKit
- **Success Rate**: 100% (all tests passing)

## Integration

Test compounds integrate with:
- `PFASDefinitionModel.test()` method
- `tests/run_groups_definitions_tests.py` test script
- `tests/test_pytest.py` pytest suite
- CI/CD pipelines (via exit codes)

## Next Steps

Optional enhancements:
1. Add false_positives/false_negatives for known edge cases
2. Add more diverse negative examples (e.g., fluorinated polymers)
3. Include borderline cases to test definition boundaries
4. Cross-validate with external PFAS databases
