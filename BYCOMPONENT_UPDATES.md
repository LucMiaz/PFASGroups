# ByComponent Algorithm Flavor Updates

## Summary
This document tracks the implementation of dual algorithm testing (bycomponent=True and bycomponent=False) across the PFASGroups project.

## Completed Updates

### 1. Enhanced Benchmark (`benchmark/enhanced_pfas_benchmark.py`)
- ✅ Updated `test_with_pfasgroups()` to accept `bycomponent` parameter
- ✅ Modified `run_enhanced_benchmark()` to test both flavors for:
  - Single group molecules
  - Multi-group pairs
  - Multi-group triplets  
  - OECD groups (1-28)
- ✅ Both flavors are now tested and results stored as `pfasgroups_result` (default) and `pfasgroups_result_bycomponent`

### 2. Review App Database (`benchmark/review-app/database/database.js`)
- ✅ Added new table `pfasgroups_results_bycomponent` to store bycomponent=True results
- ✅ Added index `idx_pfasgroups_bycomp_molecule` for performance
- ✅ Both flavors can now be stored and queried independently

### 3. Review App Import Script (`benchmark/review-app/scripts/import-benchmark-data.js`)
- ✅ Added `insertPFASGroupsResultByComponent()` method
- ✅ Updated `importBenchmarkFile()` to handle both result types:
  - `pfasgroups_result` → `pfasgroups_results` table
  - `pfasgroups_result_bycomponent` → `pfasgroups_results_bycomponent` table

### 4. Review App Server (`benchmark/review-app/server.js`)
- ✅ Updated molecule query to JOIN both `pfasgroups_results` and `pfasgroups_results_bycomponent` tables
- ✅ Added flavor mismatch detection and prioritization:
  - Highest priority: Unreviewed molecules with flavor mismatches
  - High priority: Unreviewed misclassified
  - Medium priority: Reviewed with flavor mismatches
  - Low priority: Reviewed misclassified
  - Lowest priority: Reviewed correct
- ✅ Added `has_flavor_mismatch` field to processed molecules
- ✅ Priority categories now include:
  - `flavor_mismatch_unreviewed`
  - `flavor_mismatch_reviewed`

### 5. Test Suite (`PFASgroups/tests/test_examples.py`)
- ✅ Updated `test_oecd_pfas_groups()` to test both flavors and verify:
  - High accuracy for both flavors (>50% detection rate)
  - Perfect adequation between flavors (>95% agreement)
  - Tracks and reports adequation mismatches
- ✅ Updated `test_generic_pfas_groups()` similarly for generic groups
- Results now include:
  - `detected_default` and `detected_bycomponent`
  - `all_matches_default` and `all_matches_bycomponent`
  - `adequation_match` boolean

## Remaining Updates

### 6. Test Suite - Specificity Tests (IN PROGRESS)
- ⏳ Need to update `df_test_pfas_group_specificity()` to test both flavors
- ⏳ Need to update `test_specificity()` method
- ⏳ Update test summary generation to include adequation metrics

### 7. Review App Client (NOT STARTED)
Need to update the React frontend to:
- Display both algorithm flavors side-by-side
- Show `pfasgroups_detected` and `pfasgroups_bycomponent_detected`
- Highlight molecules with flavor mismatches (visual indicator)
- Add filter for "flavor mismatches only"
- Update molecule detail view to compare both flavors
- Show adequation status in molecule cards

Client files to update:
- `benchmark/review-app/client-src/src/components/MoleculeCard.jsx`
- `benchmark/review-app/client-src/src/components/MoleculeDetail.jsx`
- `benchmark/review-app/client-src/src/components/MoleculeList.jsx`
- `benchmark/review-app/client-src/src/App.jsx` (filters)

### 8. Analysis Scripts (NOT STARTED)
Update analysis scripts to compare both flavors:
- `benchmark/enhanced_analysis.py`
- `benchmark/generate_unified_report.py`
- Add adequation rate to performance metrics
- Generate comparison charts for both flavors

## Testing Checklist

### Benchmark Testing
- [ ] Run enhanced benchmark with both flavors
- [ ] Verify JSON output contains both `pfasgroups_result` and `pfasgroups_result_bycomponent`
- [ ] Check timing data for both flavors

### Database Testing
- [ ] Import benchmark data with both flavors
- [ ] Verify both tables are populated correctly
- [ ] Test queries for flavor mismatches
- [ ] Verify prioritization logic works

### Test Suite Testing
- [ ] Run `pytest PFASgroups/tests/test_examples.py`
- [ ] Verify both flavors pass accuracy tests
- [ ] Check adequation rate is >95%
- [ ] Review adequation mismatch reports

### Client Testing (After Implementation)
- [ ] Verify both flavors display correctly
- [ ] Test flavor mismatch highlighting
- [ ] Check filter functionality
- [ ] Verify molecule detail comparison view

## Running the Updated System

### 1. Generate Benchmark Data
```bash
cd benchmark
python enhanced_pfas_benchmark.py
# Select option 1 for enhanced benchmark
```

### 2. Import to Database
```bash
cd review-app
node scripts/import-benchmark-data.js
```

### 3. Start Review App
```bash
cd review-app
npm run dev
```

### 4. Run Tests
```bash
cd ../..
pytest PFASgroups/tests/test_examples.py -v
```

## Expected Behavior

### Adequation
- **Expected**: Both flavors should detect the same PFAS groups (perfect adequation >95%)
- **Chain lengths**: May differ between flavors (this is acceptable)
- **Mismatches**: Should be rare and flagged as high priority for review

### Performance
- Both flavors should have similar:
  - Accuracy (>50% for OECD, >30% for generic)
  - Specificity (>60%)
  - Detection rates
- Timing differences are acceptable and expected

### Review App Priority
1. **Flavor mismatch + unreviewed** (highest)
2. **Misclassified + unreviewed**
3. **Unreviewed**
4. **Flavor mismatch + reviewed**
5. **Misclassified + reviewed**
6. **Reviewed correct** (lowest)

## Notes
- The `bycomponent=False` flavor is the default and original implementation
- The `bycomponent=True` flavor uses a different approach for groups with 1 SMARTS pattern
- Both flavors should ideally produce identical group detections (adequation)
- Chain length differences are acceptable between flavors
- Adequation mismatches should be investigated and potentially fixed in the core algorithm
