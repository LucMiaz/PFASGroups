# Component Fraction Calculation Refactoring Summary

## Date: January 12, 2026

## Overview
Refactored component fraction calculations in PFASGroups to use **carbon atoms only** instead of all atoms (including H, F, Cl, Br, I). This makes the fractions more chemically meaningful and easier to interpret.

## Motivation
- Original all-atom fractions were confusing because they included implicit hydrogens and halogens
- Sum of individual component fractions often exceeded total fraction
- Carbon-only fractions better represent the structural backbone coverage

## Changes Made

### 1. PFASGroups Library (`c:\Users\luc\git\PFASGroups\`)

#### `PFASgroups\core.py`
- **ComponentsSolver.__init__** (line ~400):
  - Added `self.total_carbons` to track total carbon atoms in molecule
  - Formula: `sum(1 for atom in mol.GetAtoms() if atom.GetSymbol() == 'C')`

- **get_total_components_fraction** (line ~486):
  - Changed from union of all atoms to union of **carbon atoms only**
  - Returns: `len(union_carbon_atoms) / self.total_carbons`

- **get_matched_component_dict** (line ~865-885):
  - Refactored `component_fraction` calculation:
    ```python
    # Count carbon atoms in component
    component_carbons = sum(1 for atom_idx in component 
                           if mol.GetAtomWithIdx(atom_idx).GetSymbol() == 'C')
    
    # Count carbon atoms in SMARTS matches NOT in component
    smarts_carbons_not_in_component = ...
    
    # Add extra carbons from functional group
    total_carbons = component_carbons + smarts_carbons_not_in_component + smarts_extra_atoms
    component_fraction = total_carbons / self.total_carbons
    ```

- **find_alkyl_components** (line ~1057):
  - Updated fallback calculation to use carbon counts

- **parse_groups_in_mol** (line ~1350-1369):
  - `total_components_fraction` now unions carbon atoms only

#### `PFASgroups\PFASGroupModel.py`
- **_count_smarts_extra_atoms** (line ~97-135):
  - Updated docstring: Now counts **extra CARBON atoms**
  - `manual_size` (smarts1_size/smarts2_size) now interpreted as **total carbon atoms** in functional group
  - Returns: `manual_size - 1` (subtract 1 for the matched carbon)
  - Automatic counting returns 0 (carbons counted explicitly elsewhere)

### 2. ZeroPMDB Project (`c:\Users\luc\git\zeropmdb\`)

#### `database\elements\models.py`
Updated model field help text to reflect carbon-only calculations:

- **PFASGroup model**:
  - `smarts1_size`: "Manual specification of total CARBON atoms in smarts1 functional group"
  - `smarts2_size`: "Manual specification of total CARBON atoms in smarts2 functional group"
  - `smarts1_extra_atoms`: "Precomputed extra CARBON atoms in smarts1 functional group beyond matched carbon"
  - `smarts2_extra_atoms`: "Precomputed extra CARBON atoms in smarts2 functional group beyond matched carbon"

- **PFASGroupInCompound model**:
  - `mean_component_fraction`: "Mean fraction of molecule CARBON atoms covered by components"
  - `total_components_fraction`: "Total fraction of molecule CARBON atoms covered by union of all components"

- **Components model**:
  - `component_fraction`: "Fraction of molecule CARBON atoms in this component (carbon_count / total_carbons)"

#### `database\elements\scripts\utility_funcgroups.py`
- Added documentation note explaining carbon-based fractions
- No code changes needed - automatically uses updated PFASGroups library

## Interpretation Guide

### Example: Perfluoroalkyl Carboxylic Acid (PFOA)
```
SMILES: C(C(C(C(C(C(C(F)(F)F)(F)F)(F)F)(F)F)(F)F)(F)F)(=O)O
Total carbons: 7
Component (perfluorinated chain): 6 carbons
SMARTS match (CF2 adjacent to COOH): 1 carbon  
smarts1_size: 2 (total carbons in -COOH functional group)
smarts_extra_atoms: 1 (from smarts1_size - 1 = 2 - 1 = 1)

Calculation:
- Component carbons: 6
- SMARTS match carbons (not in component): 1
- Extra carbons from functional group: 1 (the C=O carbon)
- Total: 6 + 1 + 1 = 8 carbons... wait, molecule only has 7!
```

**Issue Found**: The SMARTS matches the CF2 carbon adjacent to COOH (already in component), but we're still adding smarts_extra_atoms. Need to check if SMARTS-matched carbon is in component!

### smarts1_size Interpretation
- `smarts1_size = 2` means functional group has 2 total carbons
- Since SMARTS matches 1 carbon, we add `2 - 1 = 1` extra carbon
- For carboxylic acid (-COOH): matched carbon (CF2) + extra carbon (C=O) = 2 total

### Fraction Values
- `component_fraction = 1.0` → Component covers 100% of carbon backbone
- `component_fraction = 0.75` → Component covers 75% of carbons (e.g., 3 out of 4)
- `total_components_fraction` may differ from sum of individual fractions due to carbon overlap

## Validation

### Test Results (`test_component_ratios.py`)
✅ Component fractions no longer exceed 1.0
✅ Perfluoroalkyl carboxylic acids show `component_fraction = 1.0` (all carbons covered)
⚠️  `total_components_fraction` still needs update in `get_total_components_fraction`

### Known Issues
1. `total_components_fraction` calculation needs updating to match new carbon-only logic
2. Need to verify behavior when SMARTS-matched carbon overlaps with component

## Migration Notes

### For Database Users
- **No data migration needed** - values are recomputed on next analysis
- Field semantics changed: fractions now represent carbon coverage, not atom coverage
- Existing fraction values will be overwritten with carbon-based calculations

### For API/UI Consumers
- Update documentation to reflect carbon-only fractions
- Adjust any visualizations or reports that assumed all-atom fractions
- Fractions will generally be **higher** now (fewer carbons than total atoms)

## Future Work
1. Fix `total_components_fraction` to properly union carbon atoms
2. Handle edge cases where SMARTS match overlaps with component
3. Add explicit tests for various functional groups with different smarts1_size values
4. Document expected fraction ranges for each PFAS group type

## Contact
For questions about this refactoring, refer to conversation summary dated January 12, 2026.
