# Degradation Product Generation Functions

This document describes the new functions added to `fragmentation.py` for generating degradation products of PFAS molecules by breaking bonds that are connected to (but not part of) the fluorinated chain.

## Overview

The new functions enable systematic generation of degradation products by:
1. Identifying fluorinated chains in PFAS molecules using SMARTS patterns
2. Determining which bonds are connected to but not part of the fluorinated chain
3. Generating all possible combinations of these specific bond breaks
4. Analyzing the resulting fragments

This approach is more chemically realistic as it focuses on breaking bonds that would disconnect functional groups from the fluorinated backbone, rather than fragmenting the entire non-fluorinated portion of the molecule.

## Functions

### `find_fluorinated_chains(mol, smartsPathName='Perfluoroalkyl')`

Identifies atoms that are part of fluorinated chains in a molecule.

**Parameters:**
- `mol`: RDKit molecule object
- `smartsPathName`: SMARTS pattern name ('Perfluoroalkyl', 'Polyfluoroalkyl', etc.)

**Returns:**
- Set of atom indices that are part of fluorinated chains

### `get_bonds_connected_to_fluorinated_path(mol, fluorinated_atoms=None, smartsPathName='Perfluoroalkyl')`

Identifies bonds that are connected to but not part of the fluorinated chain.

**Parameters:**
- `mol`: RDKit molecule object
- `fluorinated_atoms`: Set of fluorinated chain atom indices (optional)
- `smartsPathName`: SMARTS pattern name

**Returns:**
- List of bond indices that connect the fluorinated chain to other molecular parts

**Key Behavior:**
- Only returns bonds where exactly one atom is in the fluorinated chain
- Excludes bonds entirely within the fluorinated chain
- Excludes bonds entirely between non-fluorinated atoms

### `get_non_fluorinated_bonds(mol, fluorinated_atoms=None, smartsPathName='Perfluoroalkyl')`

Legacy function that identifies all bonds not part of the fluorinated chain (kept for compatibility).

### `generate_degradation_products(mol, smartsPathName='Perfluoroalkyl', max_breaks=3, min_fragment_size=2, include_original=True)`

Generates degradation products by breaking bonds connected to but not part of the fluorinated chain.

**Parameters:**
- `mol`: RDKit molecule object
- `smartsPathName`: SMARTS pattern name for fluorinated chains
- `max_breaks`: Maximum number of bonds to break simultaneously
- `min_fragment_size`: Minimum number of atoms in fragments to keep
- `include_original`: Whether to include the original molecule

**Returns:**
- Dictionary with InChI keys mapping to fragment information

**Modified Behavior:**
- Now uses `get_bonds_connected_to_fluorinated_path()` instead of `get_non_fluorinated_bonds()`
- More selective bond breaking strategy
- Focuses on disconnecting functional groups from fluorinated backbone

### `generate_systematic_degradation_products(mol, smartsPathName='Perfluoroalkyl', preserve_fluorinated_chain=True, min_fragment_size=2, max_combinations=1000)`

Generates systematic degradation products exploring all possible bond breaking patterns.

**Parameters:**
- `mol`: RDKit molecule object
- `smartsPathName`: SMARTS pattern name
- `preserve_fluorinated_chain`: If True, never break bonds within fluorinated chain
- `min_fragment_size`: Minimum fragment size
- `max_combinations`: Maximum combinations to explore

**Returns:**
- Dictionary organized by number of breaks

### `analyze_degradation_pathways(mol, smartsPathName='Perfluoroalkyl', max_breaks=3)`

Analyzes degradation pathways and provides detailed information about fragments.

**Parameters:**
- `mol`: RDKit molecule object
- `smartsPathName`: SMARTS pattern name
- `max_breaks`: Maximum number of bonds to break

**Returns:**
- Dictionary with detailed pathway analysis including:
  - Original molecule information
  - Degradation summary statistics
  - Detailed product information

## Usage Examples

```python
from rdkit import Chem
from PFASgroups.fragmentation import (
    generate_degradation_products,
    analyze_degradation_pathways
)

# Create a PFAS molecule
smiles = "C(C(C(C(F)(F)F)(F)F)(F)F)(C(=O)O)(F)F"  # PFOA-like
mol = Chem.MolFromSmiles(smiles)

# Generate degradation products
products = generate_degradation_products(mol, max_breaks=2)

# Analyze pathways
analysis = analyze_degradation_pathways(mol, max_breaks=2)

print(f"Generated {len(products)} degradation products")
print(f"Products by breaks: {analysis['degradation_summary']['products_by_breaks']}")
```

## Key Features

1. **Preserves Fluorinated Chains**: The functions identify and preserve the core fluorinated chain while breaking peripheral bonds.

2. **Systematic Exploration**: Generates all possible combinations of bond breaks up to a specified limit.

3. **Fragment Analysis**: Provides detailed information about each degradation product including molecular properties and structural features.

4. **Flexible SMARTS Patterns**: Supports different types of fluorinated chains (perfluoroalkyl, polyfluoroalkyl, etc.).

5. **Chemical Validity**: Ensures generated fragments are chemically valid using RDKit sanitization.

## Applications

- Environmental fate modeling of PFAS compounds
- Predicting biotransformation products
- Understanding degradation pathways
- Identifying persistent degradation products
- Risk assessment of PFAS metabolites

## Dependencies

The functions require:
- RDKit for molecular operations
- NetworkX for graph operations (inherited from core functions)
- Access to SMARTS pattern files (fpaths.json) or fallback patterns
- Functions from core.py for molecular graph operations