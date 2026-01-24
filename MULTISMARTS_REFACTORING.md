# Multi-SMARTS Refactoring Summary

## Overview
Successfully refactored PFASGroups to support multiple SMARTS patterns with counts per group, replacing the previous `smarts1`/`smarts2` structure.

## Changes Made

### 1. JSON Structure (`PFAS_groups_smarts.json`)
**Before:**
```json
{
  "smarts1": "[pattern]",
  "smarts2": "[pattern]" or null,
  "smarts1_size": 2,
  "smarts2_size": 2 or null
}
```

**After:**
```json
{
  "smarts": {
    "[pattern]": 1,  // count = 1 for normal groups
    "[pattern]": 2   // count = 2 for diacids/disulfonics
  }
}
```

### 2. PFASGroupModel.py
**Key Changes:**
- `__init__` now extracts SMARTS patterns and counts from dictionary:
  - `self.smarts_str`: tuple of SMARTS strings (dict keys)
  - `self.smarts_count`: tuple of required counts (dict values)
  - `self.smarts`: list of compiled RDKit Mol objects
  - `self.smarts_extra_atoms`: list of extra atom counts per SMARTS
- Fixed empty dictionary handling to avoid unpacking errors
- Removed `smarts1`, `smarts2`, `smarts1_size`, `smarts2_size` attributes

### 3. core.py
**Major Changes:**

#### New Function: `find_components_with_all_smarts`
- Handles groups requiring multiple SMARTS patterns or multiple copies
- Checks that components contain sufficient matches for each SMARTS
- Supports diacids (requiring 2 copies of same SMARTS)
- Supports groups with multiple different SMARTS patterns

#### Updated Function: `parse_groups_in_mol`
- Replaced `smarts1`/`smarts2` branching logic
- Now checks if group has multiple SMARTS or count > 1
- Routes to appropriate handler:
  - Cyclic groups: use first SMARTS with `find_aryl_components`
  - Multiple SMARTS or count > 1: use `find_components_with_all_smarts`
  - Single SMARTS with count=1: use `find_alkyl_components` with first SMARTS

#### Updated Function: `ComponentsSolver.get_matched_component_dict`
- Changed `smarts_extra_atoms` calculation from:
  ```python
  if pfas_group.smarts2 is not None:
      smarts_extra_atoms = pfas_group.smarts1_extra_atoms + pfas_group.smarts2_extra_atoms
  else:
      smarts_extra_atoms = pfas_group.smarts1_extra_atoms * len(smarts_matches)
  ```
- To:
  ```python
  if pfas_group.smarts_extra_atoms is not None:
      smarts_extra_atoms = sum(pfas_group.smarts_extra_atoms) * len(smarts_matches)
  ```

#### Updated Docstrings
- Revised notes in `parse_groups_in_mol` to reflect new logic

## Testing Results
All tests passed successfully:

### Test 1: Structure Verification ✓
- All 59 groups have new attributes (`smarts_str`, `smarts_count`, `smarts_extra_atoms`)
- Old attributes (`smarts1`, `smarts2`) successfully removed
- Diacid groups correctly have count=2

### Test 2: Simple PFCA (PFOA) ✓
- Correctly identified as "Perfluoroalkyl carboxylic acids"
- Single SMARTS with count=1 works properly

### Test 3: Diacid Detection ✓
- Perfluorooctanedicarboxylic acid correctly identified
- "Perfluoroalkyl dicarboxylic acids" group matched
- Requires 2 copies of carboxylic acid SMARTS in same component

## Backward Compatibility
⚠️ **Breaking Changes:**
- Code using `pfas_group.smarts1` or `pfas_group.smarts2` will break
- Use `pfas_group.smarts[0]`, `pfas_group.smarts[1]` etc. instead
- Code using `smarts1_extra_atoms`/`smarts2_extra_atoms` needs updating
- Use `pfas_group.smarts_extra_atoms[i]` or `sum(pfas_group.smarts_extra_atoms)`

## Benefits
1. **More flexible**: Can define groups with any number of SMARTS patterns
2. **More explicit**: Count requirements are clearly specified in the data
3. **Cleaner code**: No special-casing for smarts1 vs smarts2
4. **Future-proof**: Easy to add groups requiring 3+ different SMARTS patterns

## Files Modified
- `PFASgroups/PFASGroupModel.py`
- `PFASgroups/core.py`
- `PFASgroups/data/PFAS_groups_smarts.json` (converted by script)

## Files Created
- `convert_smarts.py` (conversion script)
- `test_multismarts.py` (test suite)
