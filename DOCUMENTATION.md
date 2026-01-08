# PFASgroups

[![PyPI version](https://badge.fury.io/py/PFASgroups.svg)](https://badge.fury.io/py/PFASgroups)
[![Python 3.7+](https://img.shields.io/badge/python-3.7+-blue.svg)](https://www.python.org/downloads/)
[![License: CC BY-NC 4.0](https://img.shields.io/badge/License-CC%20BY--NC%204.0-lightgrey.svg)](https://creativecommons.org/licenses/by-nc/4.0/)
[![Powered by RDKit](https://img.shields.io/badge/Powered%20by-RDKit-3838ff.svg)](https://www.rdkit.org/)

A comprehensive Python cheminformatics package for automated detection, classification, and analysis of Per- and Polyfluoroalkyl Substances (PFAS) in chemical databases.

## Overview

**PFASgroups** provides three main capabilities:

1. **PFAS Group Identification and Chain Analysis**: Automated detection of 55 functional groups directly connected to fluorinated chains, based on OECD terminology
2. **Fluorinated Chain Length Quantification**: Determination of per- or polyfluorinated alkyl chain lengths connected to functional groups
3. **Customization**: Both functional groups and pathway patterns can be customized via JSON configuration files

The algorithm combines SMARTS pattern matching, molecular formula constraints, and graph-based pathfinding using RDKit and NetworkX.

## Installation

### From PyPI

```bash
pip install PFASgroups
```

### From Conda-Forge

```bash
conda install -c conda-forge pfasgroups
```

### From Source

```bash
git clone https://github.com/yourusername/PFASGroups.git
cd PFASGroups
pip install -e .
```

### Requirements

- Python >= 3.7
- RDKit
- NumPy
- Pandas
- NetworkX
- tqdm
- svgutils

---

## Quick Start

### Python API

```python
from PFASgroups import parse_smiles, generate_fingerprint

# Parse PFAS structures
smiles_list = ["C(C(F)(F)F)F", "FC(F)(F)C(F)(F)C(=O)O"]
results = parse_smiles(smiles_list)

# Each result contains: (PFASGroup, match_count, chain_lengths, matched_chains)
for i, smiles_str in enumerate(smiles_list):
    print(f"\nAnalyzing: {smiles_str}")
    for group, n_matches, chain_lengths, chains in results[i]:
        print(f"  {group.name}: {n_matches} matches, chains: {chain_lengths}")

# Generate fingerprints for machine learning
fingerprints, group_info = generate_fingerprint(smiles_list)
```

### Command Line Interface

```bash
# Parse SMILES strings
pfasgroups parse "C(C(F)(F)F)F" "FC(F)(F)C(F)(F)C(=O)O"

# Generate fingerprints
pfasgroups fingerprint "C(C(F)(F)F)F" --format dict

# List available PFAS groups
pfasgroups list-groups

# List available pathway types
pfasgroups list-paths
```

---

## Module Reference

### core.py - Core Functions

The main module providing PFAS group identification and parsing functionality.

#### Main Functions

##### `parse_smiles(smiles, bycomponent=False, output_format='list', **kwargs)`

Parse SMILES string(s) and return PFAS group information.

**Parameters:**
- `smiles` (str or list): Single SMILES string or list of SMILES strings
- `bycomponent` (bool): Whether to use component-based analysis for cyclic structures
- `output_format` (str): Output format - 'list', 'dataframe', or 'csv'
- `**kwargs`: Additional parameters (pfas_groups, smartsPaths, etc.)

**Returns:**
- List of tuples: `(PFASGroup, match_count, chain_lengths, matched_chains)`

```python
from PFASgroups import parse_smiles

# Basic usage
results = parse_smiles("FC(F)(F)C(F)(F)C(=O)O")

# With component-based analysis
results = parse_smiles(smiles_list, bycomponent=True)

# Output as DataFrame
df = parse_smiles(smiles_list, output_format='dataframe')
```

##### `parse_groups_in_mol(mol, bycomponent=False, **kwargs)`

Parse an RDKit molecule object directly.

**Parameters:**
- `mol` (rdkit.Chem.Mol): RDKit molecule object
- `bycomponent` (bool): Whether to use component-based analysis
- `**kwargs`: Additional parameters (formula, pfas_groups, smartsPaths, etc.)

**Returns:**
- List of tuples: `(PFASGroup, match_count, chain_lengths, matched_chains)`

##### `generate_fingerprint(smiles, selected_groups=None, representation='vector', count_mode='binary', **kwargs)`

Generate PFAS group fingerprints for machine learning applications.

**Parameters:**
- `smiles` (str or list): SMILES string(s) to fingerprint
- `selected_groups` (list or range): Specific group IDs to include (default: all)
- `representation` (str): 'vector', 'dict', or 'sparse'
- `count_mode` (str): 'binary', 'count', or 'max_chain'

**Returns:**
- Tuple: `(fingerprints, group_info)`

```python
from PFASgroups import generate_fingerprint

# Binary vector (default)
fps, info = generate_fingerprint(smiles_list)

# Dictionary format
fps, info = generate_fingerprint(smiles_list, representation='dict')

# Count-based fingerprints
fps, info = generate_fingerprint(smiles_list, count_mode='count')

# Select specific groups
fps, info = generate_fingerprint(smiles_list, selected_groups=range(28, 52))
```

##### `plot_pfasgroups(smiles, **kwargs)`

Visualize PFAS group assignments on molecular structures.

**Parameters:**
- `smiles` (str or list): SMILES string(s) to visualize
- `display` (bool): Whether to display the plot
- `path` (str): Path to save the plot image
- `svg` (bool): Whether to generate SVG output
- `subwidth`, `subheight` (int): Dimensions of each sub-image
- `ncols` (int): Number of columns in grid layout
- `addAtomIndices` (bool): Whether to add atom indices

```python
from PFASgroups import plot_pfasgroups

# Basic visualization
plot_pfasgroups("FC(F)(F)C(F)(F)C(=O)O")

# Save as SVG
plot_pfasgroups(smiles_list, svg=True, path="output.svg")
```

#### Helper Functions

##### `get_smartsPaths(**kwargs)`

Get the default pathway SMARTS definitions.

```python
from PFASgroups import get_smartsPaths

paths = get_smartsPaths()
# Returns dict with keys: 'Perfluoroalkyl', 'Polyfluoroalkyl', 'Polyfluoro', 'Polyfluorobr'
```

##### `get_PFASGroups(**kwargs)`

Get the default PFAS group definitions.

```python
from PFASgroups import get_PFASGroups

groups = get_PFASGroups()
# Returns list of 55 PFASGroup objects
```

##### `compile_smartsPath(chain_smarts, end_smarts)`

Compile SMARTS patterns for custom pathway definitions.

```python
from PFASgroups import compile_smartsPath, get_smartsPaths

paths = get_smartsPaths()
paths['Perchlorinated'] = compile_smartsPath(
    "[C;X4](Cl)(Cl)!@!=!#[C;X4](Cl)(Cl)",  # chain pattern
    "[C;X4](Cl)(Cl)Cl"                      # end pattern
)
```

##### `compile_smartsPaths(paths_dict)`

Compile multiple pathway definitions from a dictionary.

```python
from PFASgroups import compile_smartsPaths

custom_paths = {
    'Perchlorinated': {
        'chain': '[C;X4](Cl)(Cl)!@!=!#[C;X4](Cl)(Cl)',
        'end': '[C;X4](Cl)(Cl)Cl'
    }
}
compiled = compile_smartsPaths(custom_paths)
```

---

### PFASGroupModel.py - PFAS Group Data Model

Defines the `PFASGroup` class representing PFAS group definitions.

#### `PFASGroup` Class

```python
from PFASgroups import PFASGroup

group = PFASGroup(
    id=999,
    name="My Custom Group",
    smarts1="[C](F)(F)F",           # Primary SMARTS pattern
    smarts2="C(=O)O",               # Secondary SMARTS pattern (optional)
    smartsPath="Perfluoroalkyl",    # Pathway type
    constraints={                    # Formula constraints
        "eq": {"O": 2},
        "gte": {"F": 1},
        "only": ["C", "F", "O", "H"],
        "rel": {"C": {"atoms": ["F"], "add": 0.5, "div": 2}}
    },
    max_dist_from_CF=2              # Max distance from fluorinated carbon (default: 2)
)
```

**Attributes:**
- `id` (int): Unique identifier
- `name` (str): Descriptive name
- `smarts1` (Mol): Primary SMARTS pattern (RDKit Mol object)
- `smarts2` (Mol): Secondary SMARTS pattern (optional)
- `smartsPath` (str): Pathway type ('Perfluoroalkyl', 'Polyfluoroalkyl', 'cyclic', or None)
- `constraints` (dict): Molecular formula constraints
- `max_dist_from_CF` (int): Maximum bond distance from fluorinated carbon terminal (applies to groups without formula constraints when `bycomponent=True`)

**Constraint Types:**
- `eq`: Exact element count (e.g., `{"O": 2}`)
- `gte`: Minimum element count (e.g., `{"F": 1}`)
- `lte`: Maximum element count
- `only`: Allowed elements only (e.g., `["C", "F", "O", "H"]`)
- `rel`: Relative element ratios (e.g., C = (F + H) / 2 - 1)

**Methods:**
- `formula_dict_satisfies_constraints(formula_dict)`: Check if a formula satisfies all constraints

---

### PFASDefinitionModel.py - PFAS Definition Data Model

Defines the `PFASDefinition` class for broader PFAS definitions (e.g., EPA, OECD definitions).

#### `PFASDefinition` Class

```python
from PFASgroups import PFASDefinition

definition = PFASDefinition(
    id=1,
    name="OECD PFAS Definition",
    smarts=["[C](F)(F)F", "[C](F)(F)[C](F)(F)"],
    fluorineRatio=0.5,
    description="OECD definition of PFAS",
    includeHydrogen=True,
    requireBoth=False
)
```

**Attributes:**
- `id` (int): Unique identifier
- `name` (str): Definition name
- `smarts` (list): List of SMARTS patterns
- `fluorineRatio` (float): Minimum fluorine ratio threshold
- `description` (str): Human-readable description
- `includeHydrogen` (bool): Include H in ratio calculation
- `requireBoth` (bool): Require both SMARTS match AND fluorine ratio

**Methods:**
- `applies_to_molecule(mol_or_smiles, formula=None, **kwargs)`: Check if definition applies

---

### generate_homologues.py - Homologue Series Generation

Generate shorter-chain analogues by iteratively removing CF₂ moieties.

#### `generate_homologues(mol, smartsPathName='Perfluoroalkyl', **kwargs)`

Generate homologue series by removing repeating units from fluorinated chains.

**Parameters:**
- `mol` (Mol): RDKit molecule object
- `smartsPathName` (str): Pathway type ('Perfluoroalkyl' or 'Polyfluoroalkyl')
- `repeating` (str): SMARTS pattern for repeating unit (default: 'C(F)F')
- `base_repeating` (list): Base atoms not to remove (default: ['C'])

**Returns:**
- Dict mapping InChIKey → {formula: Mol}

```python
from PFASgroups import generate_homologues
from rdkit import Chem

mol = Chem.MolFromSmiles("FC(F)(F)C(F)(F)C(F)(F)C(F)(F)C(=O)O")
homologues = generate_homologues(mol, smartsPathName='Perfluoroalkyl')

# Returns shorter chain analogues:
# PFOA (C8) → PFHpA (C7) → PFHxA (C6) → PFPeA (C5) → PFBA (C4)
for inchikey, formulas in homologues.items():
    for formula, mol in formulas.items():
        print(f"{formula}: {Chem.MolToSmiles(mol)}")
```

---

### fragmentation.py - Molecular Fragmentation

Fragment molecules based on bond dissociation energies (BDE).

#### `generate_fragments(mol, nb_breakingpoints=5)`

Break molecule on weakest bonds using BDE values.

**Parameters:**
- `mol` (Mol): RDKit molecule object
- `nb_breakingpoints` (int): Number of bonds to break

**Returns:**
- List of fragment Mol objects

```python
from PFASgroups.fragmentation import generate_fragments
from rdkit import Chem

mol = Chem.MolFromSmiles("FC(F)(F)C(F)(F)C(=O)O")
fragments = generate_fragments(mol, nb_breakingpoints=3)
```

#### `generate_degradation_products(mol, **kwargs)`

Generate potential degradation products by breaking weak bonds.

---

### generate_mol.py - Molecule Generation

Utilities for generating random PFAS-like molecules for testing and validation.

#### `generate_random_carbon_chain(n, cycle=False, alkene=False, alkyne=False)`

Generate a random carbon chain molecule.

**Parameters:**
- `n` (int): Number of carbon atoms
- `cycle` (bool): Generate cyclic structure
- `alkene` (bool): Include C=C bonds
- `alkyne` (bool): Include C≡C bonds

```python
from PFASgroups.generate_mol import generate_random_carbon_chain

mol = generate_random_carbon_chain(6)
```

#### `append_functional_group(mol, group_smiles, insertion='attach', m=1, **kwargs)`

Append a functional group to a molecule.

**Parameters:**
- `mol` (Mol): Base molecule
- `group_smiles` (str): Functional group as SMILES
- `insertion` (str): 'attach' or 'insert'
- `m` (int): Number of groups to attach

---

### draw_mols.py - Molecular Visualization

Utilities for drawing and combining molecular images.

#### `plot_mol(mol, **kwargs)`

Draw a single molecule.

**Parameters:**
- `mol` (Mol): RDKit molecule object
- `svg` (bool): Generate SVG output
- `subwidth`, `subheight` (int): Image dimensions
- `addAtomIndices` (bool): Add atom indices
- `addBondIndices` (bool): Add bond indices

#### `plot_mols(mols, **kwargs)`

Draw multiple molecules in a grid layout.

#### `draw_images(imgs, buffer=1, ncols=2, svg=False)`

Combine multiple images into a grid.

---

## PFAS Group Classifications

The package includes 55 PFAS group classifications:

### OECD-Defined Groups (IDs 1-28)

Strict structural requirements based on the OECD reconciling terminology report:

| ID | Name | Alias |
|----|------|-------|
| 1 | Perfluoroalkyl carboxylic acids | PFCAs |
| 2 | Polyfluoroalkyl carboxylic acid | PolyFCAs |
| 3 | Perfluoroalkyl dicarboxylic acids | PFdiCAs |
| 4 | Perfluoroalkylether carboxylic acids | PFECAs |
| 5 | Polyfluoroalkylether carboxylic acid | PolyFECAs |
| 6 | Perfluoroalkyl sulfonic acids | PFSAs |
| 7 | Polyfluoroalkyl sulfonic acid | PolyFSAs |
| 8 | Perfluoroalkyl disulfonic acids | PFdiSAs |
| 9 | Perfluoroalkylether sulfonic acids | PFESAs |
| 10 | Polyfluoroalkylether sulfonic acid | PolyFESAs |
| 11 | Perfluoroalkyl sulfinic acids | PFSiAs |
| 12 | Perfluoroalkyl phosphonic acids | PFPAs |
| 13 | Perfluoroalkyl phosphinic acids | PFPiAs |
| 14 | Perfluoroalkyl alcohols | - |
| 15 | Fluorotelomer alcohols | FTOHs |
| 16 | Perfluoropolyethers | PFPEs |
| 17 | Hydrofluoroethers | HFEs |
| 18 | Perfluoroalkene | - |
| 19 | Hydrofluoroolefins | HFOs |
| 20 | Hydrofluorocarbons | HFCs |
| 21 | Semi-fluorinated alkanes | SFAs |
| 22 | Side-chain fluorinated aromatics | - |
| 23 | Perfluoroalkane | - |
| 24 | Perfluoroalkyl-tert-amines | - |
| 25 | Perfluoroalkyl iodides | PFAIs |
| 26 | Perfluoroalkane sulfonyl fluorides | PASFs |
| 27 | Perfluoroalkyl ketones | - |
| 28 | Semi-fluoroalkyl ketones | - |

### Generic Classifications (IDs 29-55)

Looser criteria capturing fluorinated functional groups:

| ID | Name |
|----|------|
| 29 | alcohol |
| 30 | ketone |
| 31 | ether |
| 32 | ester |
| 33 | carboxylic acid |
| 34 | amide |
| 35 | acyl halide |
| 36 | sulfonic acid |
| 37 | sulfenic acid |
| 38 | sulfinic acid |
| 39 | phosphonic acid |
| 40 | phosphinic acid |
| 41 | ethene |
| 42 | iodide |
| 43 | sulfonamide |
| 44 | Heterocyclic azole |
| 45 | Heterocyclic azine |
| 46 | benzodioxole |
| 47 | amine |
| 48 | alkane |
| 49 | alkene |
| 50 | alkyne |
| 51 | Side-chain aromatics |
| 52 | Perfluoro cyclic compounds |
| 53 | Polyfluoro cyclic compounds |
| 54 | Perfluoroaryl compounds |
| 55 | Polyfluoroaryl compounds |
| 56-57 | Peroxides, Benzoyl peroxides |

---

## Custom Configuration

### Creating Custom PFAS Groups

```python
from PFASgroups import get_PFASGroups, PFASGroup, parse_smiles

# Get default groups
groups = get_PFASGroups()

# Add custom group
groups.append(PFASGroup(
    id=999,
    name="My Custom PFAS",
    smarts1="[C](F)(F)F",
    smarts2="[N+](=O)[O-]",
    smartsPath="Perfluoroalkyl",
    constraints={"gte": {"F": 3}},
    max_dist_from_CF=2
))

# Use custom groups
results = parse_smiles(["FC(F)(F)C(F)(F)[N+](=O)[O-]"], pfas_groups=groups)
```

### Creating Custom Pathway Types

```python
from PFASgroups import compile_smartsPath, get_smartsPaths, parse_smiles

paths = get_smartsPaths()

# Add chlorinated pathway
paths['Perchlorinated'] = compile_smartsPath(
    "[C;X4](Cl)(Cl)!@!=!#[C;X4](Cl)(Cl)",  # chain pattern
    "[C;X4](Cl)(Cl)Cl"                      # end pattern
)

results = parse_smiles(["ClC(Cl)(Cl)C(Cl)(Cl)C(=O)O"], smartsPaths=paths)
```

### JSON Configuration Files

#### PFAS Groups (PFAS_groups_smarts.json)

```json
[
  {
    "id": 1,
    "name": "Perfluoroalkyl carboxylic acids",
    "alias": "PFCAs",
    "base_functional_groups": ["carboxylic acid"],
    "main_group": "Perfluoroalkyl acids",
    "smarts1": "[#6$([#6][#6](=O)([OH1,O-]))]",
    "smarts2": null,
    "smartsPath": "Perfluoroalkyl",
    "constraints": {
      "eq": {"O": 2},
      "only": ["C", "F", "O", "H"],
      "gte": {"F": 1},
      "rel": {"C": {"atoms": ["F"], "add": 0.5, "div": 2}}
    },
    "max_dist_from_CF": 2
  }
]
```

#### Pathway Types (fpaths.json)

```json
{
  "Perfluoroalkyl": {
    "chain": "[C;X4;H0](F)(F)!@!=!#[C;X4;H0](F)(F)",
    "end": "[C;X4;H0](F)(F)F"
  },
  "Polyfluoroalkyl": {
    "chain": "[C;X4;H1](F)!@!=!#[C;X4](F)",
    "end": "[C;X4;H1](F)F"
  }
}
```

---

## Algorithm Overview

The detection algorithm processes molecules through five sequential steps:

### 1. Molecule Preprocessing
- Add explicit hydrogen atoms using `Chem.AddHs()`
- Sanitize molecular structure
- If sanitization fails, attempt iterative fragmentation

### 2. Formula Constraint Validation
Evaluate elemental composition against group requirements:
- **Relative ratios** (`rel`): e.g., C = (F+H+Cl+Br+I)/2 - 1
- **Minimum counts** (`gte`): e.g., F ≥ 1
- **Maximum counts** (`lte`)
- **Exact counts** (`eq`): e.g., O = 2
- **Element restrictions** (`only`): e.g., only C, F, H allowed

### 3. SMARTS Pattern Matching
Apply SMARTS pattern(s) to identify candidate functional groups.

### 4a. Fluorinated Path Finding
For non-cyclic groups:
1. Identify all substructure matches for source and target SMARTS
2. Compute shortest paths using Dijkstra's algorithm (NetworkX)
3. Validate intermediate atoms satisfy fluorinated pathway constraints
4. Return chain information (length, atom indices, pathway type)

### 4b. Fluorinated Connected Components
For cyclic groups (`smartsPath="cyclic"`):
- Extract subgraph of fluorinated atoms
- Return connected components instead of chains

### 5. Result Aggregation
- Remove redundant chains (subsets of longer chains)
- Prioritize longest valid chains
- Return all matches with chain information

---

## Testing and Benchmarking

### Unit Tests

The package includes comprehensive unit tests in the `tests/` directory:

```bash
# Run all tests
python -m pytest

# Run specific test file
python test_python_only.py
```

Test files include:
- `test_python_only.py`: Core functionality tests
- `test_compile_functions.py`: SMARTS compilation tests
- `test_python_js_comparison.py`: Python vs JavaScript parity tests

### Benchmarking Against PFAS-Atlas

The `benchmark/` directory contains tools for comparing PFASgroups against PFAS-Atlas:

```bash
cd benchmark

# Run all benchmarks
python enhanced_pfas_benchmark.py

# Generate analysis reports
python analyze_timing.py data/pfas_timing_benchmark_*.json
python enhanced_analysis.py data/pfas_enhanced_benchmark_*.json
```

Benchmark outputs include:
- Timing comparisons
- Classification accuracy metrics
- Complex branched molecule analysis
- OECD group coverage analysis

### Interactive Review App

A web-based review application for comparing results:

```bash
cd benchmark/review-app
npm install
node scripts/reimport-all.js
npm start
```

---

## Web Tool

A JavaScript implementation using RDKitJS provides browser-based PFAS analysis without installation:

- **PFAS group identification**: All 55 classifications with chain analysis
- **Real-time analysis**: Immediate feedback for single molecules
- **Batch analysis**: Process lists of SMILES

The web tool mirrors the Python module's logic and was validated against the same test set.

---

## License

This work is licensed under a [Creative Commons Attribution-NoDerivatives 4.0 International License](http://creativecommons.org/licenses/by-nd/4.0/).

Contact the author for exceptions to the No Derivatives term.

---

## Citation

If you use PFASgroups in your research, please cite:

> Miaz, L.T., Cousins, I.T. (2026). Automatic Determination and Classification of Per- and Polyfluoroalkyl Substances. *Journal of Cheminformatics* (in preparation).

---

## Acknowledgments

This project is part of the [ZeroPM project](https://zeropm.eu/) (WP2) and has received funding from the European Union's Horizon 2020 research and innovation programme under grant agreement No 101036756.

Developed at the [Department of Environmental Science](https://aces.su.se) at Stockholm University.

<img alt="EU logo" src="https://zeropm.eu/wp-content/uploads/2021/12/flag_yellow_low.jpg" width=100/>
<img alt="ZeroPM logo" src="https://zeropm.eu/wp-content/uploads/2022/01/ZeroPM-logo.png" width=250/>
<img alt="Stockholm University logo" src="https://eu01web.zoom.us/account/branding/p/5065401a-9915-4baa-9c16-665dcd743470.png" width=200/>

---

## Contributing

Contributions are welcome! Please:

1. Fork the repository
2. Create a feature branch
3. Add tests for new functionality
4. Submit a pull request

---

## Changelog

### v1.2.2 (Current)
- Added `max_dist_from_CF` parameter for functional groups
- Improved homologue generation
- Enhanced benchmarking tools
- JavaScript web tool parity

### v1.2.0
- Added CLI interface
- Custom configuration support
- Fingerprint generation
- DataFrame output format

### v1.0.0
- Initial release
- 55 PFAS group classifications
- Core parsing and visualization
