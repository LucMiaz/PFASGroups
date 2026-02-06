# Benchmark System Changes - February 2026

## Summary

The benchmark system has been refactored to use **JSON-driven test generation**, eliminating hardcoded SMILES mappings and ensuring consistency between production code and test sets.

## Changes Made

### 1. Removed Hardcoded SMILES Mappings

**Before:**
```python
smiles_map = {
    29: {'smiles': 'O[H]', 'mode': 'attach'},  # alcohol
    30: {'smiles': 'C(=O)', 'mode': 'insert'},  # ketone
    31: {'smiles': 'O', 'mode': 'insert'},  # WRONG! This is ether, not aldehyde
    32: {'smiles': 'C(=O)OC', 'mode': 'insert'},  # WRONG! Shifted by one
    # ... 80+ more hardcoded entries
}
```

**After:**
```python
# Read directly from PFAS_groups_smarts.json
for group in self.pfas_groups:
    if 'test' in group and 'generate' in group['test']:
        generate_info = group['test']['generate']
        smiles = generate_info.get('smiles', '')
        mode = generate_info.get('mode', 'attach')
        is_telomer = generate_info.get('is_telomer', False)
        # ... use these values
```

### 2. Fixed Group Mapping Errors

The hardcoded mappings had several critical errors that caused 0% accuracy in heatmap plots:

| Group ID | Name | Was (Wrong) | Now (Correct) |
|----------|------|-------------|---------------|
| 31 | aldehyde | `O` (ether) | `C(=O)[H]` |
| 32 | ether | `C(=O)OC` (ester) | `O` |
| 33 | ester | `C(=O)O` (carboxylic acid) | `C(=O)OC` |
| 34 | carboxylic acid | `C(=O)N` (amide) | `C(=O)O` |
| ... | ... | ... | ... |

Groups 31-45 were all shifted by one position due to the initial error at group 31.

### 3. Fixed Telomer Group IDs

**Before:**
- Telomer groups incorrectly started at ID 69 (overwriting non-telomer group)
- Created duplicate definitions

**After:**
- Telomer groups correctly start at ID 74
- Groups 69-73 are preserved as non-telomer functional groups

### 4. Added Dynamic Group Loading

The system now:
- ✅ Automatically detects all groups with `test.generate` fields
- ✅ Skips groups without valid test generation patterns
- ✅ Warns when groups are missing required fields
- ✅ Supports adding new groups without code changes

## JSON Structure

Each group in `PFAS_groups_smarts.json` can include a `test.generate` field:

```json
{
  "id": 31,
  "name": "aldehyde",
  "test": {
    "category": "functional",
    "generate": {
      "smiles": "C(=O)[H]",
      "mode": "attach",
      "is_telomer": false
    },
    "examples": ["..."],
    "counter-examples": ["..."]
  }
}
```

**Fields:**
- `smiles`: SMILES pattern to attach or insert into test molecules
- `mode`: `"attach"` (terminal) or `"insert"` (internal)
- `is_telomer`: Boolean flag for fluorotelomer groups

## Benefits

### 1. Single Source of Truth
- Test patterns match production group definitions
- No risk of divergence between test and production code

### 2. Maintainability
- Add new groups by editing JSON only
- No need to update benchmark script
- Changes propagate automatically

### 3. Consistency
- All test generation uses same patterns
- Reduces bugs from manual updates
- Clear documentation in JSON file

### 4. Transparency
- Easy to see what SMILES each group uses for testing
- Can validate test patterns against group definitions
- Better debugging when tests fail

## Migration Guide

### For Adding New Groups

**Before (old way - required code changes):**
1. Add group to `PFAS_groups_smarts.json`
2. Edit `enhanced_pfas_benchmark.py`
3. Add entry to `smiles_map` dictionary
4. Remember correct group ID and position
5. Re-run benchmarks

**After (new way - JSON only):**
1. Add group to `PFAS_groups_smarts.json` with `test.generate` field
2. Re-run benchmarks (automatically includes new group)

### For Modifying Test Patterns

**Before:**
1. Find group in benchmark script's `smiles_map`
2. Update SMILES
3. Hope you got the group ID right

**After:**
1. Find group in `PFAS_groups_smarts.json`
2. Update `test.generate.smiles`
3. Changes automatically applied on next run

## Example: Adding a New Group

```json
{
  "id": 115,
  "name": "My New Functional Group",
  "componentSmarts": "Polyfluoroalkyl",
  "smarts": {
    "[custom smarts pattern]": 1
  },
  "test": {
    "category": "functional",
    "generate": {
      "smiles": "C(=O)OC(=O)",
      "mode": "insert",
      "is_telomer": false
    },
    "examples": [
      "FC(F)(F)C(F)(F)C(=O)OC(=O)C(F)(F)C(F)(F)F"
    ]
  }
}
```

This group will automatically:
- Be detected by the benchmark system
- Have test molecules generated
- Be included in accuracy reports
- Appear in comparison heatmaps

## Verification

Run the test to verify correct loading:

```bash
cd /home/luc/git/PFASGroups/benchmark
python3 -c "
from scripts.enhanced_pfas_benchmark import EnhancedPFASBenchmark
benchmark = EnhancedPFASBenchmark()
print(f'Loaded {len(benchmark.functional_smarts)} functional groups')
for gid in [29, 31, 32, 74, 77]:
    if gid in benchmark.functional_smarts:
        entry = benchmark.functional_smarts[gid]
        print(f'{gid}: {entry[\"name\"]} - {entry[\"smiles\"]}')
"
```

Expected output should show:
- Group 31 with `C(=O)[H]` (aldehyde)
- Group 32 with `O` (ether)
- Group 74 with telomer=True

## Next Steps

1. ✅ Test script loads correctly from JSON
2. 🔄 Re-run all benchmarks to generate corrected data
3. 🔄 Update analysis reports with fixed accuracy metrics
4. 🔄 Verify heatmap plots show correct accuracy values

## Files Modified

- `benchmark/scripts/enhanced_pfas_benchmark.py` - Removed hardcoded mappings
- `docs/changelog.rst` - Added version 1.2.3 section
- `docs/benchmarking.rst` - Added JSON-driven generation section
- `benchmark/README.md` - Added overview of JSON system

## Documentation

See updated documentation at:
- `/docs/benchmarking.rst` - Full guide on JSON-driven test generation
- `/docs/changelog.rst` - Version 1.2.3 changes
- `PFASgroups/data/PFAS_groups_smarts.json` - Group definitions with test patterns
