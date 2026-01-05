# PFASgroups Function Renaming - December 27, 2025

## Summary of Changes

To improve API clarity and generality, the following functions have been renamed:

### Core Function Renames

| Old Name | New Name | Changes |
|----------|----------|---------|
| `parse_pfas()` | `parse_smiles()` | Now accepts single SMILES or list |
| N/A | `parse_mol()` | NEW: Parses RDKit molecules (single or list) |
| `parse_PFAS_groups()` | `parse_groups_in_mol()` | More descriptive name |
| `generate_pfas_fingerprint()` | `generate_fingerprint()` | Shorter, clearer name |

### New Features

**parse_smiles()** and **parse_mol()** now accept both:
- Single input: Returns list of matches for one molecule
- List input: Returns list of lists of matches

**Example:**
```python
from PFASgroups import parse_smiles, parse_mol
from rdkit import Chem

# Single SMILES
result = parse_smiles('FC(F)(F)C(F)(F)C(=O)O')
# Returns: list of (group, n_matches, n_CFchain, chains) tuples

# List of SMILES  
results = parse_smiles(['FC(F)(F)C(F)(F)C(=O)O', 'C(C(F)(F)F)F'])
# Returns: list of lists

# Single molecule
mol = Chem.MolFromSmiles('FC(F)(F)C(F)(F)C(=O)O')
result = parse_mol(mol)

# List of molecules
mols = [Chem.MolFromSmiles(s) for s in ['FC(F)(F)C(F)(F)C(=O)O', 'C(C(F)(F)F)F']]
results = parse_mols(mols)
```

## Files Updated

### Core Package Files
- ✅ `PFASgroups/core.py` - Function definitions and implementations
- ✅ `PFASgroups/__init__.py` - Exports
- ✅ `PFASgroups/cli.py` - Command-line interface

### Test Files
- ✅ `PFASgroups/tests/test_examples.py` - Unit tests
- ✅ `test_compile_functions.py` - Compile functions tests
- ✅ `test_python_only.py` - Python-only tests
- ✅ `test_python_js_comparison.py` - JS comparison tests

### Benchmark Files
- ✅ `benchmark/test_oecd_analysis.py`
- ✅ `benchmark/enhanced_pfas_benchmark.py`
- ✅ `benchmark/generate_unified_report.py`

### Documentation Files
- ✅ `USER_GUIDE.md` - Complete user guide
- ✅ `QUICK_REFERENCE.md` - Quick reference
- ✅ `README.md` - Main README

## Backward Compatibility

⚠️ **Breaking Changes**: Old function names are no longer available. Code using the old names will need to be updated.

### Migration Guide

```python
# OLD CODE
from PFASgroups import parse_pfas, parse_PFAS_groups, generate_pfas_fingerprint

results = parse_pfas(['smiles1', 'smiles2'])
matches = parse_PFAS_groups(mol, formula)
fps, info = generate_pfas_fingerprint(smiles)

# NEW CODE
from PFASgroups import parse_smiles, parse_mol, parse_groups_in_mol, generate_fingerprint

results = parse_smiles(['smiles1', 'smiles2'])  # or single SMILES
matches = parse_groups_in_mol(mol, formula=formula)
fps, info = generate_fingerprint(smiles)
```

## Verification

All tests pass successfully:
- ✅ Function renaming complete
- ✅ Single/list input handling works
- ✅ Compile functions test passes
- ✅ Documentation updated

## Module Name

The module name **PFASgroups** remains unchanged as requested.
