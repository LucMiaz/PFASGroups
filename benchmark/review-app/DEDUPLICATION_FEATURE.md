# Deduplication Feature Added to Benchmark Data Import

## Overview
Added automatic deduplication to the benchmark data import process to remove duplicate entries while preserving legitimate variations in classification results.

## Changes Made

### 1. Modified Files
- `review-app/scripts/import-benchmark-data.js`
  - Added `deduplicateRecords()` method for benchmark data
  - Added `deduplicateTimingRecords()` method for timing data
  - Integrated deduplication into `importBenchmarkFile()` and `importTimingData()`

### 2. Deduplication Logic

#### Benchmark Data (OECD, Enhanced, Complex Branched, Non-Fluorinated)
**Criteria for uniqueness:**
- SMILES string (molecular structure)
- PFASGroups detected groups
- PFASGroups bycomponent detected groups  
- Atlas first class
- Atlas second class

**Behavior:**
- Keeps the first occurrence of each unique combination
- Same molecule with different classifications → kept as separate records
- True duplicates (identical SMILES and classifications) → removed

**Example:**
```javascript
// KEPT: Different classifications
Record 1: SMILES="C(F)(F)F", PFASGroups=[CF3], Atlas="Type1"
Record 2: SMILES="C(F)(F)F", PFASGroups=[CF2], Atlas="Type2"

// REMOVED: Exact duplicate
Record 1: SMILES="C(F)(F)F", PFASGroups=[CF3], Atlas="Type1"
Record 2: SMILES="C(F)(F)F", PFASGroups=[CF3], Atlas="Type1" ← Removed
```

#### Timing Data
**Criteria for uniqueness:**
- SMILES string only

**Behavior:**
- Keeps the first occurrence of each unique SMILES
- Removes all subsequent occurrences regardless of timing values

**Example:**
```javascript
// KEPT
Record 1: SMILES="C(F)(F)F", time=0.001ms

// REMOVED
Record 2: SMILES="C(F)(F)F", time=0.002ms ← Removed
Record 3: SMILES="C(F)(F)F", time=0.003ms ← Removed
```

### 3. Output
When duplicates are found, the import process now displays:
```
Importing pfas_oecd_benchmark_20251226_002740.json...
  Removed 150 duplicate(s) from pfas_oecd_benchmark_20251226_002740.json
✓ Imported 3264 records from pfas_oecd_benchmark_20251226_002740.json
```

## Testing

### Test Script
Created `scripts/test-deduplication.js` to verify the logic:

**Test Results:**
- ✅ Benchmark data: Correctly identifies duplicates based on SMILES + classifications
- ✅ Timing data: Correctly identifies duplicates based on SMILES only
- ✅ Edge cases: Handles null/undefined values gracefully

### Running Tests
```bash
cd review-app
node scripts/test-deduplication.js
```

## Benefits

1. **Data Quality**: Eliminates truly duplicate records from datasets
2. **Storage Efficiency**: Reduces database size
3. **Analysis Accuracy**: Prevents duplicate records from skewing statistics
4. **Preserves Variation**: Keeps records where the same molecule has different classifications (important for validation)

## Usage

The deduplication is automatic and requires no configuration. When you run:
```bash
cd benchmark
.\quick-start.ps1
```

Or manually import data:
```bash
cd review-app
node scripts/import-benchmark-data.js
```

Duplicates will be automatically detected and removed with a summary showing how many were removed from each file.

## Implementation Details

### Key Functions

**`deduplicateRecords(records)`**
- Creates a composite key from SMILES and all classification results
- Uses a Map to track seen keys
- Returns deduplicated array

**`deduplicateTimingRecords(records)`**
- Creates key from SMILES only
- Uses a Set to track seen SMILES
- Returns deduplicated array

### Performance
- O(n) time complexity where n = number of records
- Minimal memory overhead (stores only keys in Map/Set)
- Processes 10,000+ records in milliseconds

## Future Enhancements

Potential improvements:
1. Add option to keep last occurrence instead of first
2. Add statistics about which specific duplicates were removed
3. Create a duplicate report file for review
4. Add command-line flags to control deduplication behavior

## Notes

- Deduplication happens **before** database insertion
- Original JSON files remain unchanged
- The "first occurrence" is kept based on the order in the JSON file
- Case-sensitive SMILES matching
- Classification groups are sorted before comparison to handle order variations
