# PFASgroups User Guide

Complete guide for using the PFASgroups package to parse, analyze, and fingerprint Per- and Polyfluoroalkyl Substances (PFAS).

## Table of Contents

1. [Installation](#installation)
2. [Quick Start](#quick-start)
3. [Python API](#python-api)
4. [Custom Configuration](#custom-configuration)
5. [Command-Line Interface](#command-line-interface)
6. [Advanced Usage](#advanced-usage)
7. [API Reference](#api-reference)

---

## Installation

### From Source

```bash
cd PFASGroups
pip install -e .
```

### Requirements

- Python >= 3.7
- RDKit
- numpy
- pandas
- networkx
- tqdm
- svgutils

---

## Quick Start

### Python API

```python
from PFASgroups import parse_pfas, generate_pfas_fingerprint

# Parse PFAS structures
smiles_list = ["C(C(F)(F)F)F", "FC(F)(F)C(F)(F)C(=O)O"]
results = parse_pfas(smiles_list)

# Generate fingerprints
fingerprints, group_info = generate_pfas_fingerprint(smiles_list)
print(fingerprints)
```

### Command Line

```bash
# Parse SMILES strings
pfasgroups parse "C(C(F)(F)F)F" "FC(F)(F)C(F)(F)C(=O)O"

# Generate fingerprints
pfasgroups fingerprint "C(C(F)(F)F)F" --format dict
```

---

## Python API

### Basic Parsing

Parse SMILES strings to identify PFAS groups:

```python
from PFASgroups import parse_pfas

# Single or multiple SMILES
smiles = ["C(C(F)(F)F)F", "FC(F)(F)C(F)(F)C(=O)O"]
results = parse_pfas(smiles)

# Results contain: (PFASGroup, n_matches, n_CFchain, chains)
for i, smiles_str in enumerate(smiles):
    print(f"\nAnalyzing: {smiles_str}")
    for group, n_matches, n_CFchain, chains in results[i]:
        print(f"  {group.name}: {n_matches} matches, {len(chains)} chains")
```

### Component-Based Analysis

Use component-based analysis for certain PFAS groups:

```python
results = parse_pfas(smiles, bycomponent=True)
```

### Working with RDKit Molecules Directly

You can also parse individual molecules directly:

```python
from PFASgroups import parse_PFAS_groups
from rdkit import Chem
from rdkit.Chem.rdMolDescriptors import CalcMolFormula

mol = Chem.MolFromSmiles("C(C(F)(F)F)F")
# optional: provide molecular formula (otherwise it will be calculated using the CalcMolFormula as below)
formula = CalcMolFormula(mol)
results = parse_PFAS_groups(mol, formula=formula)
```

### Fingerprint Generation

Generate PFAS group fingerprints for machine learning:

```python
from PFASgroups import generate_pfas_fingerprint

smiles = ["C(C(F)(F)F)F", "FC(F)(F)C(F)(F)C(=O)O"]

# Binary vector (default)
fps, info = generate_pfas_fingerprint(smiles)
print(fps)  # numpy array

# Dictionary format
fps, info = generate_pfas_fingerprint(smiles, representation='dict')
print(fps)  # [{'group_name': count, ...}, ...]

# Sparse (only non-zero)
fps, info = generate_pfas_fingerprint(smiles, representation='sparse')

# Count-based instead of binary
fps, info = generate_pfas_fingerprint(
    smiles, 
    count_mode='count'  # 'binary', 'count', or 'max_chain'
)

# Select specific groups
fps, info = generate_pfas_fingerprint(
    smiles,
    selected_groups=range(28, 52)  # or [28, 29, 30, ...]
)
```

### Visualization

Plot PFAS groups on molecular structures:

```python
from PFASgroups import parse_PFAS_groups, plot_pfasgroups
from rdkit import Chem
from rdkit.Chem.rdMolDescriptors import CalcMolFormula

smiles = "FC(F)(F)C(F)(F)C(=O)O"
mol = Chem.MolFromSmiles(smiles)
formula = CalcMolFormula(mol)

# Parse groups
matches = parse_PFAS_groups(mol, formula=formula)

# Plot
plot_pfasgroups([(smiles, matches)])
```

---

## Custom Configuration

PFASgroups allows you to use custom definitions for pathtype patterns (`chain` element, used both for identifying atoms in the core chain or in connected components, and `end` element, used to identify the end of a chain) and PFAS groups by passing them as parameters to the main functions. You can either load custom files entirely, or fetch the defaults and extend/modify them.

### Compile Functions for Custom Paths

**NEW:** When creating custom path types, use the **`compile_smartsPath()`** and **`compile_smartsPaths()`** helper functions to automatically preprocess SMARTS patterns. These functions handle all required RDKit preprocessing steps automatically.

**Quick Example:**
```python
from PFASgroups import compile_smartsPath, get_smartsPaths

paths = get_smartsPaths()  # Get defaults

# Add custom path (automatically preprocessed!)
paths['Perchlorinated'] = compile_smartsPath(
    "[C;X4](Cl)(Cl)!@!=!#[C;X4](Cl)(Cl)",  # chain pattern
    "[C;X4](Cl)(Cl)Cl"                     # end pattern
)
```

See [compile_smartsPath documentation](#compile_smartspathchain_smarts-end_smarts) for more details.

### Configuration Files

Two JSON files control the behavior:

1. **fpaths.json** - Defines pathtype SMARTS patterns for chain detection
2. **PFAS_groups_smarts.json** - Defines PFAS group structures and constraints

### Method 1: Using Helper Functions

Use the built-in helper functions to load configuration. These functions allow you to:
- Load custom files entirely, OR
- Fetch defaults and extend/modify them

#### Loading Custom Files

```python
from PFASgroups import get_smartsPaths, get_PFASGroups, parse_pfas

# Load custom path definitions from file
custom_paths = get_smartsPaths()

# Load custom PFAS groups from file
custom_groups = get_PFASGroups()



# Use custom configuration
results = parse_pfas(
    ["C(C(F)(F)F)F"],
    smartsPaths=custom_paths,
    pfas_groups=custom_groups
)
```

#### Extending Default Configuration

Fetch the defaults and add your custom definitions:

```python
from PFASgroups import get_smartsPaths, get_PFASGroups, parse_pfas, PFASGroup, compile_smartsPath

# Get default paths
paths = get_smartsPaths()

# Add a custom pathtype using compile_smartsPath helper
paths['MyCustomPath'] = compile_smartsPath(
    "[C;X4](F)(F)!@!=!#[C;X4](F)",  # chain pattern
    "[C;X4](F)F"                     # end pattern
)

# Get default groups and add custom ones
groups = get_PFASGroups()
groups.append(PFASGroup(
    id=999,
    name="My Custom Group",
    smarts1="[C](F)(F)F",
    smarts2="[N+](=O)[O-]",
    smartsPath="MyCustomPath",
    constraints={"nF": [3, None]}
))

# Use extended configuration
results = parse_pfas(
    ["C(C(F)(F)F)F"],
    smartsPaths=paths,
    pfas_groups=groups
)
```

#### Filtering/Subsetting Default Configuration

Fetch defaults and filter to only what you need:

```python
from PFASgroups import get_PFASGroups, parse_pfas

# Get all default groups
all_groups = get_PFASGroups()

# Filter to only perfluorinated groups
perfluoro_groups = [
    g for g in all_groups 
    if g.smartsPath == 'Perfluoroalkyl'
]

# Use filtered configuration
results = parse_pfas(
    ["FC(F)(F)C(F)(F)C(=O)O"],
    pfas_groups=perfluoro_groups
)

# Or filter by ID range
specific_groups = [g for g in all_groups if 28 <= g.id <= 52]
```

### Method 2: Loading Directly from Files

Load and construct custom groups manually:

```python
from PFASgroups import PFASGroup, parse_pfas
import json

# Load custom PFAS groups from file
with open('path/to/custom_groups.json', 'r') as f:
    custom_groups_data = json.load(f)

custom_groups = [PFASGroup(**x) for x in custom_groups_data]

# Use with parsing
results = parse_pfas(
    ["C(C(F)(F)F)F"],
    pfas_groups=custom_groups
)
```

### Method 3: Command Line

Use custom files via CLI flags:

```bash
# Parse with custom groups
pfasgroups parse --groups-file custom_groups.json "C(C(F)(F)F)F"

# Fingerprint with custom config
pfasgroups fingerprint \
    --fpaths-file custom_fpaths.json \
    --groups-file custom_groups.json \
    --input smiles.txt
```

### Creating Custom Configuration Files

#### custom_fpaths.json

Define pathtype SMARTS patterns:

```json
{
  "Perfluoroalkyl": {
    "chain": "[C;X4;H0](F)(F)!@!=!#[C;X4;H0](F)(F)",
    "end": "[C;X4;H0](F)(F)F"
  },
  "Polyfluoroalkyl": {
    "chain": "[C;X4;H1](F)!@!=!#[C;X4](F)",
    "end": "[C;X4;H1](F)F"
  },
  "CustomPathType": {
    "chain": "[your SMARTS pattern]",
    "end": "[your SMARTS pattern]"
  }
}
```

#### custom_groups.json

Define PFAS groups:

```json
[
  {
    "id": 1,
    "name": "My Custom PFAS Group",
    "smarts1": "[C](F)(F)F",
    "smarts2": "C(=O)O",
    "smartsPath": "Perfluoroalkyl",
    "constraints": {
      "nF": [3, null],
      "nC": [2, null]
    }
  }
]
```

**Fields:**
- `id`: Unique identifier
- `name`: Descriptive name
- `smarts1`: SMARTS pattern for functional group 1
- `smarts2`: SMARTS pattern for functional group 2 (optional)
- `smartsPath`: Pathtype name from fpaths.json or "cyclic"
- `constraints`: Chemical formula constraints (min, max or null)

---

## Command-Line Interface

The CLI provides easy access to PFASgroups functionality from the terminal.

### Parse Command

Identify PFAS groups in structures:

```bash
# Parse SMILES from command line
pfasgroups parse "C(C(F)(F)F)F" "FC(F)(F)C(F)(F)C(=O)O"

# Parse from file (one SMILES per line)
pfasgroups parse --input smiles.txt --output results.json

# Component-based analysis
pfasgroups parse --bycomponent "C(C(F)(F)F)F"

# Pretty-print JSON
pfasgroups parse "C(C(F)(F)F)F" --pretty
```

**Input file format** (smiles.txt):
```
C(C(F)(F)F)F
FC(F)(F)C(F)(F)C(=O)O
FC(F)C(F)(F)F
```

**Output format** (results.json):
```json
[
  {
    "smiles": "C(C(F)(F)F)F",
    "groups": [
      {
        "name": "Perfluoroalkyl moiety (1)",
        "id": 28,
        "n_matches": 1,
        "n_CFchain": [3],
        "chains_count": 1
      }
    ]
  }
]
```

### Fingerprint Command

Generate PFAS fingerprints:

```bash
# Basic fingerprint
pfasgroups fingerprint "C(C(F)(F)F)F"

# From file
pfasgroups fingerprint --input smiles.txt --output fingerprints.json

# Select specific groups (range)
pfasgroups fingerprint "C(C(F)(F)F)F" --groups 28-52

# Select specific groups (list)
pfasgroups fingerprint "C(C(F)(F)F)F" --groups 28,29,30

# Different formats
pfasgroups fingerprint "C(C(F)(F)F)F" --format dict
pfasgroups fingerprint "C(C(F)(F)F)F" --format sparse
pfasgroups fingerprint "C(C(F)(F)F)F" --format int

# Count mode
pfasgroups fingerprint "C(C(F)(F)F)F" --count-mode count
pfasgroups fingerprint "C(C(F)(F)F)F" --count-mode max_chain
```

### List Groups Command

Display available PFAS groups:

```bash
# List all default groups
pfasgroups list-groups

# Save to file
pfasgroups list-groups --output groups.json

# List from custom file
pfasgroups list-groups --groups-file custom_groups.json
```

**Note:** In Python, use `get_PFASGroups()` to load and extend these groups.

### List Paths Command

Display available pathtype definitions:

```bash
# List all default path types
pfasgroups list-paths

# Save to file
pfasgroups list-paths --output paths.json

# List from custom file
pfasgroups list-paths --fpaths-file custom_fpaths.json
```

**Output example:**
```json
{
  "total_paths": 4,
  "paths": [
    {
      "name": "Perfluoroalkyl",
      "chain_smarts": "[C;X4;H0](F)(F)!@!=!#[C;X4;H0](F)(F)",
      "end_smarts": "[C;X4;H0](F)(F)F"
    }
  ],
  "note": "Use get_smartsPaths() in Python to load and extend these paths"
}
```

**Note:** In Python, use `get_smartsPaths()` to load and extend these paths.

### Validate Config Command

Validate configuration files:

```bash
# Validate default config
pfasgroups validate-config

# Validate custom files
pfasgroups validate-config --groups-file custom_groups.json --fpaths-file custom_fpaths.json
```

---

## Advanced Usage

### Batch Processing

Process multiple files efficiently:

```python
import pandas as pd
from PFASgroups import generate_pfas_fingerprint

# Read SMILES from CSV
df = pd.read_csv('compounds.csv')
smiles_list = df['smiles'].tolist()

# Generate fingerprints in batches
batch_size = 1000
all_fps = []

for i in range(0, len(smiles_list), batch_size):
    batch = smiles_list[i:i+batch_size]
    fps, info = generate_pfas_fingerprint(batch, representation='sparse')
    all_fps.extend(fps)

# Save results
df['pfas_fingerprint'] = all_fps
df.to_csv('compounds_with_fps.csv', index=False)
```

### Integration with Machine Learning

```python
import numpy as np
from sklearn.ensemble import RandomForestClassifier
from PFASgroups import generate_pfas_fingerprint

# Training data
train_smiles = ["C(C(F)(F)F)F", "FC(F)(F)C(F)(F)C(=O)O", ...]
train_labels = [1, 0, ...]  # Your labels

# Generate features
X_train, info = generate_pfas_fingerprint(
    train_smiles,
    representation='vector',
    count_mode='binary'
)

# Train model
clf = RandomForestClassifier()
clf.fit(X_train, train_labels)

# Predict on new data
test_smiles = ["FC(F)C(F)(F)F", ...]
X_test, _ = generate_pfas_fingerprint(
    test_smiles,
    selected_groups=info['selected_indices'],
    representation='vector',
    count_mode='binary'
)
predictions = clf.predict(X_test)
```

### Using Custom Groups for Screening

Create custom groups tailored to your screening needs:

```python
from PFASgroups import PFASGroup, parse_pfas
import json

# Define screening criteria for short-chain PFAS only
custom_groups = [{
    "id": 1,
    "name": "Short-chain perfluorinated",
    "smarts1": "[C](F)(F)F",
    "smarts2": "C(=O)O",
    "smartsPath": "Perfluoroalkyl",
    "constraints": {"nC": [2, 6]}  # C2-C6 only
}]

# Convert to PFASGroup objects
groups = [PFASGroup(**g) for g in custom_groups]

# Screen compounds
compounds = ["FC(F)(F)C(F)(F)C(=O)O",  # C3 - will match
             "FC(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(=O)O"]  # C8 - won't match

results = parse_pfas(compounds, pfas_groups=groups)

for smiles, matches in zip(compounds, results):
    if matches:
        print(f"{smiles}: MATCHES short-chain criteria")
    else:
        print(f"{smiles}: Does not match")
```

---

## API Reference

### Core Functions

#### `parse_pfas(smiles_list, bycomponent=False, **kwargs)`

Parse SMILES strings and identify PFAS groups.

**Parameters:**
- `smiles_list` (list): List of SMILES strings
- `bycomponent` (bool): Use component-based analysis
- `**kwargs`: Additional parameters:
  - `pfas_groups` (list): Custom PFASGroup objects
  - `smartsPaths` (dict): Custom path SMARTS patterns

**Returns:**
- `list`: Match results for each SMILES

**Example:**
```python
# With custom groups
results = parse_pfas(smiles_list, pfas_groups=custom_groups)
```

#### `parse_PFAS_groups(mol, bycomponent=False, **kwargs)`

Parse a single RDKit molecule and identify PFAS groups.

**Parameters:**
- `mol` (rdkit.Chem.Mol): RDKit molecule object
- `bycomponent` (bool): Use component-based analysis
- `**kwargs`: Additional parameters:
  - `formula` (str): Molecular formula (calculated if not provided)
  - `pfas_groups` (list): Custom PFASGroup objects
  - `smartsPaths` (dict): Custom path SMARTS patterns

**Returns:**
- `list`: List of (PFASGroup, n_matches, n_CFchain, chains) tuples

#### `generate_pfas_fingerprint(smiles, selected_groups=None, representation='vector', count_mode='binary', pfas_groups=None, **kwargs)`

Generate PFAS fingerprints.

**Parameters:**
- `smiles` (str or list): SMILES string(s)
- `selected_groups` (list/range): Group indices to include
- `representation` (str): 'vector', 'dict', 'sparse', 'detailed', or 'int'
- `count_mode` (str): 'binary', 'count', or 'max_chain'
- `pfas_groups` (list): Custom PFASGroup objects
- `**kwargs`: Additional parameters

**Returns:**
- `fingerprints`: Fingerprint(s) in specified format
- `group_info` (dict): Information about groups used

### Helper Functions

#### `get_smartsPaths(filename=None, **kwargs)`

Load SMARTS path patterns for chain detection.

**Parameters:**
- `filename` (str): Path to custom fpaths.json file (uses default if None)

**Returns:**
- `dict`: Dictionary mapping path names to [chain_mol, end_mol] pairs
  - Keys: Path type names (e.g., 'Perfluoroalkyl', 'Polyfluoroalkyl')
  - Values: List of two preprocessed RDKit Mol objects [chain_smarts, end_smarts]

**Examples:**

```python
# Load defaults
paths = get_smartsPaths()
print(paths.keys())  # dict_keys(['Perfluoroalkyl', 'Polyfluoroalkyl', ...])

# Load from custom file
custom_paths = get_smartsPaths(filename='my_fpaths.json')

# Extend defaults with custom path
paths = get_smartsPaths()
from rdkit import Chem
custom = Chem.MolFromSmarts("[C](F)")
custom.UpdatePropertyCache()
Chem.GetSymmSSSR(custom)
custom.GetRingInfo().NumRings()
paths['MyPath'] = [custom, custom]

# Use in parsing
results = parse_pfas(smiles, smartsPaths=paths)
```

#### `get_PFASGroups(filename=None, **kwargs)`

Load PFAS group definitions.

**Parameters:**
- `filename` (str): Path to custom PFAS_groups_smarts.json file (uses default if None)

**Returns:**
- `list`: List of PFASGroup objects with attributes:
  - `id`: Unique identifier
  - `name`: Group name
  - `smarts1`, `smarts2`: SMARTS patterns
  - `smartsPath`: Path type name
  - `constraints`: Formula constraints dict

**Examples:**

```python
# Load defaults
groups = get_PFASGroups()
print(f"Found {len(groups)} default groups")

# Load from custom file
custom_groups = get_PFASGroups(filename='my_groups.json')

# Extend defaults
from PFASgroups import PFASGroup
groups = get_PFASGroups()
groups.append(PFASGroup(
    id=999,
    name="My Custom Group",
    smarts1="[C](F)(F)F",
    smarts2="C(=O)O",
    smartsPath="Perfluoroalkyl",
    constraints={"nC": [2, 6]}
))

# Filter defaults
groups = get_PFASGroups()
short_chain = [g for g in groups if g.id in range(28, 40)]

# Use in parsing
results = parse_pfas(smiles, pfas_groups=groups)
```

#### `compile_smartsPath(chain_smarts, end_smarts)`

Compile a pair of SMARTS patterns into a ready-to-use path definition.

This function preprocesses SMARTS patterns for chain and end groups, preparing them for use in PFAS parsing functions. It handles all necessary RDKit preprocessing steps automatically.

**Parameters:**
- `chain_smarts` (str): SMARTS pattern for the repeating chain unit
- `end_smarts` (str): SMARTS pattern for the terminal group

**Returns:**
- `list`: List containing [chain_mol, end_mol] where both are preprocessed RDKit Mol objects

**Examples:**

```python
from PFASgroups import compile_smartsPath, get_smartsPaths, parse_pfas

# Get defaults and add a custom path
paths = get_smartsPaths()

# Add perchlorinated path (Cl instead of F)
paths['Perchlorinated'] = compile_smartsPath(
    "[C;X4;H0](Cl)(Cl)!@!=!#[C;X4;H0](Cl)(Cl)",  # chain
    "[C;X4;H0](Cl)(Cl)Cl"                         # end
)

# Use in parsing
results = parse_pfas(smiles, smartsPaths=paths)
```

#### `compile_smartsPaths(paths_dict)`

Compile multiple SMARTS path definitions from a dictionary.

This function takes a dictionary of path definitions and preprocesses all of them at once, making it convenient to define multiple custom path types.

**Parameters:**
- `paths_dict` (dict): Dictionary with structure:
  ```python
  {
      'PathName': {'chain': 'SMARTS', 'end': 'SMARTS'},
      ...
  }
  ```

**Returns:**
- `dict`: Dictionary mapping path names to [chain_mol, end_mol] pairs

**Examples:**

```python
from PFASgroups import compile_smartsPaths, get_smartsPaths, parse_pfas

# Define multiple custom paths
custom_paths = {
    'Perchlorinated': {
        'chain': '[C;X4;H0](Cl)(Cl)!@!=!#[C;X4;H0](Cl)(Cl)',
        'end': '[C;X4;H0](Cl)(Cl)Cl'
    },
    'MixedHalo': {
        'chain': '[C;X4]([F,Cl])!@!=!#[C;X4]([F,Cl])',
        'end': '[C;X4]([F,Cl])([F,Cl])[F,Cl]'
    }
}

# Compile all at once
compiled_paths = compile_smartsPaths(custom_paths)

# Merge with defaults
paths = get_smartsPaths()
paths.update(compiled_paths)

# Use in parsing
results = parse_pfas(smiles, smartsPaths=paths)
```

### Classes

#### `PFASGroup`

Represents a PFAS group definition.

**Attributes:**
- `id`: Unique identifier
- `name`: Group name
- `smarts1`: SMARTS pattern for functional group 1
- `smarts2`: SMARTS pattern for functional group 2 (optional)
- `smartsPath`: Pathtype name or "cyclic"
- `constraints`: Dictionary of formula constraints

---

## Examples

### Example 1: Basic Analysis

```python
from PFASgroups import parse_pfas

smiles = [
    "C(C(F)(F)F)F",                    # Perfluoropropane
    "FC(F)(F)C(F)(F)C(=O)O",          # PFPA
    "FC(F)C(F)(F)C(F)(F)C(=O)O"       # Mixed perfluoro/polyfluoro
]

results = parse_pfas(smiles)

for i, smi in enumerate(smiles):
    print(f"\n{smi}:")
    if results[i]:
        for group, n_matches, chains, _ in results[i]:
            print(f"  - {group.name}")
    else:
        print("  - No PFAS groups detected")
```

### Example 2: Custom Screening Workflow

```python
from PFASgroups import get_PFASGroups, parse_pfas
import json

# Load default groups
all_groups = get_PFASGroups()

# Filter to only perfluorinated groups (example)
perfluoro_groups = [g for g in all_groups if 'Perfluoro' in g.smartsPath]

# Screen compounds
compounds = ["FC(F)(F)C(F)(F)C(=O)O", "CCCCCC"]
results = parse_pfas(compounds, pfas_groups=perfluoro_groups)

print(f"Found {sum(len(r) for r in results)} perfluorinated groups")
```

### Example 3: Extending Default Configuration

```python
from PFASgroups import get_PFASGroups, get_smartsPaths, parse_pfas, PFASGroup
from rdkit import Chem

# Scenario: Add organization-specific PFAS groups to defaults

# Get defaults
default_groups = get_PFASGroups()
default_paths = get_smartsPaths()

# Add custom path type for partially fluorinated chains
custom_chain = Chem.MolFromSmarts("[C;X4](F)!@!=!#[C;X4]")
custom_chain.UpdatePropertyCache()
Chem.GetSymmSSSR(custom_chain)
custom_chain.GetRingInfo().NumRings()

custom_end = Chem.MolFromSmarts("[C;X4](F)([H,C])([H,C])[H,C]")
custom_end.UpdatePropertyCache() 
Chem.GetSymmSSSR(custom_end)
custom_end.GetRingInfo().NumRings()

default_paths['PartiallyFluorinated'] = [custom_chain, custom_end]

# Add custom group using new path type
default_groups.append(PFASGroup(
    id=1000,
    name="Partially fluorinated carboxylic acid",
    smarts1="[C](F)",
    smarts2="C(=O)O",
    smartsPath="PartiallyFluorinated",
    constraints={"nF": [1, None], "nC": [3, None]}
))

# Now use extended configuration
compounds = [
    "FC(H)(H)C(H)(H)C(=O)O",  # Will match new group
    "FC(F)(F)C(F)(F)C(=O)O"   # Will match default groups
]

results = parse_pfas(
    compounds,
    smartsPaths=default_paths,
    pfas_groups=default_groups
)

for smiles, matches in zip(compounds, results):
    print(f"\n{smiles}:")
    for group, _, _, _ in matches:
        print(f"  - {group.name}")
```

### Example 4: Creating Simple Custom PFAS Groups

```python
from PFASgroups import get_PFASGroups, parse_pfas, PFASGroup

# Load defaults
groups = get_PFASGroups()

# Example 1: Simple perfluorinated alcohol
groups.append(PFASGroup(
    id=2000,
    name="Perfluorinated alcohol",
    smarts1="[C](F)(F)F",          # CF3 group
    smarts2="[C][OH]",              # Alcohol group
    smartsPath="Perfluoroalkyl",    # Use existing perfluoroalkyl path
    constraints={"nF": [3, None]}   # At least 3 fluorines
))

# Example 2: Polyfluorinated ether
groups.append(PFASGroup(
    id=2001,
    name="Polyfluorinated ether",
    smarts1="[C](F)",               # Any C-F bond
    smarts2="[O][C]",               # Ether linkage
    smartsPath="Polyfluoroalkyl",   # Polyfluoroalkyl path
    constraints={"nF": [1, None]}   # At least 1 fluorine
))

# Example 3: PFAS with sulfonate group (similar to PFOS)
groups.append(PFASGroup(
    id=2002,
    name="Custom perfluoroalkyl sulfonate",
    smarts1="[C](F)(F)F",          # CF3 group
    smarts2="S(=O)(=O)[OH]",        # Sulfonic acid group
    smartsPath="Perfluoroalkyl",
    constraints={"nF": [5, None], "nS": [1, 1]}  # At least 5 F, exactly 1 S
))

# Example 4: Fluorotelomer alcohol (simpler version)
groups.append(PFASGroup(
    id=2003,
    name="Simple fluorotelomer alcohol",
    smarts1="[C](F)(F)F",          # Perfluorinated end
    smarts2="[C;X4]([H])([H])[C]([H])([H])O",  # -CH2CH2OH tail
    smartsPath="Perfluoroalkyl",
    constraints={"nF": [3, None]}
))

# Test the custom groups
test_smiles = [
    "FC(F)(F)C(F)(F)CO",           # Perfluorinated alcohol
    "FC(F)(F)C(F)COC",             # Polyfluorinated ether
    "FC(F)(F)C(F)(F)S(=O)(=O)O",  # Perfluoroalkyl sulfonate
    "FC(F)(F)C(F)(F)CCO"           # Fluorotelomer alcohol
]

results = parse_pfas(test_smiles, pfas_groups=groups)

for smiles, matches in zip(test_smiles, results):
    print(f"\n{smiles}:")
    for group, n_matches, n_CFchain, _ in matches:
        if group.id >= 2000:  # Only show our custom groups
            print(f"  ✓ {group.name}")
```

### Example 5: Creating Custom Path Types (Chlorinated Analogs)

```python
from PFASgroups import get_smartsPaths, get_PFASGroups, parse_pfas, PFASGroup, compile_smartsPath, compile_smartsPaths

# Method 1: Using compile_smartsPath for individual paths
paths = get_smartsPaths()

# Perchlorinated chains (Cl instead of F)
paths['Perchlorinated'] = compile_smartsPath(
    "[C;X4;H0](Cl)(Cl)!@!=!#[C;X4;H0](Cl)(Cl)",  # chain
    "[C;X4;H0](Cl)(Cl)Cl"                         # end
)

# Mixed fluoro/chloro chains
paths['MixedHaloalkyl'] = compile_smartsPath(
    "[C;X4]([F,Cl])([F,Cl])!@!=!#[C;X4]([F,Cl])",
    "[C;X4]([F,Cl])([F,Cl])[F,Cl,H]"
)

# Method 2: Using compile_smartsPaths for multiple paths at once
custom_paths_dict = {
    'Polychloroalkyl': {
        'chain': '[C;X4](Cl)!@!=!#[C;X4](Cl)',
        'end': '[C;X4](Cl)([H,C])([H,C])[H,C]'
    },
    'Perbrominated': {
        'chain': '[C;X4;H0](Br)(Br)!@!=!#[C;X4;H0](Br)(Br)',
        'end': '[C;X4;H0](Br)(Br)Br'
    }
}

# Compile all at once and merge with existing paths
new_paths = compile_smartsPaths(custom_paths_dict)
paths.update(new_paths)

# Now create groups that use these custom paths
groups = get_PFASGroups()

# Perchlorinated carboxylic acid
groups.append(PFASGroup(
    id=3000,
    name="Perchlorinated carboxylic acid",
    smarts1="[C](Cl)(Cl)Cl",
    smarts2="C(=O)O",
    smartsPath="Perchlorinated",
    constraints={"nCl": [3, None], "nC": [2, None]}
))

# Mixed halogenated alcohol
groups.append(PFASGroup(
    id=3001,
    name="Mixed fluoro-chloro alcohol",
    smarts1="[C]([F,Cl])",
    smarts2="CO",
    smartsPath="MixedHaloalkyl",
    constraints={"nF": [1, None]}
))

# Polychlorinated sulfonate
groups.append(PFASGroup(
    id=3002,
    name="Polychlorinated sulfonate",
    smarts1="[C](Cl)",
    smarts2="S(=O)(=O)O",
    smartsPath="Polychloroalkyl",
    constraints={"nCl": [1, None]}
))

# Test chlorinated compounds
test_chlorinated = [
    "ClC(Cl)(Cl)C(Cl)(Cl)C(=O)O",     # Perchlorinated acid
    "FC(Cl)(Cl)C(F)(Cl)CO",            # Mixed fluoro-chloro
    "ClC(Cl)(H)C(Cl)(H)S(=O)(=O)O",   # Polychlorinated sulfonate
    "BrC(Br)(Br)C(Br)(Br)C(=O)O"      # Perbrominated acid
]

results = parse_pfas(
    test_chlorinated,
    smartsPaths=paths,
    pfas_groups=groups
)

for smiles, matches in zip(test_chlorinated, results):
    print(f"\n{smiles}:")
    if matches:
        for group, n_matches, n_CFchain, _ in matches:
            if group.id >= 3000:  # Show our custom chlorinated groups
                print(f"  ✓ {group.name} (chain length: {n_CFchain})")
    else:
        print("  No custom groups matched")
```

**Key Points for Custom Paths:**
1. **Use compile helpers**: `compile_smartsPath()` for single paths, `compile_smartsPaths()` for multiple
2. **Chain pattern**: Defines repeating units with `!@!=!#` (non-ring single/double/triple bonds)
3. **End pattern**: Defines terminal groups
4. **Match halogens**: Use `[F,Cl,Br,I]` to match any halogen, or specific ones
5. **Test thoroughly**: Check both positive and negative cases

### Example 6: Large-Scale Processing

```bash
# Prepare data file
cat > compounds.txt << EOF
C(C(F)(F)F)F
FC(F)(F)C(F)(F)C(=O)O
CCCCCC
EOF

# Screen compounds
pfasgroups parse --input compounds.txt --output results.json --pretty

# Generate fingerprints
pfasgroups fingerprint --input compounds.txt --groups 28-52 --output fps.json
```

---

## Troubleshooting

### RDKit Errors

```
Chem.AtomValenceException
```

**Solution:** Check SMARTS patterns in custom configuration files. Use RDKit to validate:
```python
from rdkit import Chem
smarts = "[C](F)(F)F"
mol = Chem.MolFromSmarts(smarts)
print("Valid!" if mol else "Invalid SMARTS")
```

### Invalid JSON

```
json.JSONDecodeError: Expecting property name enclosed in double quotes
```

**Solution:** Validate JSON syntax:
```bash
python -m json.tool your_file.json
```

### Custom Configuration Not Working

Make sure you're passing the parameters correctly:
```python
# Correct
custom_groups = get_PFASGroups(filename='custom.json')
results = parse_pfas(smiles, pfas_groups=custom_groups)

# Also correct
custom_paths = get_smartsPaths(filename='custom_fpaths.json')
results = parse_pfas(smiles, smartsPaths=custom_paths)
```

---

## Support

For issues, questions, or contributions:
- GitHub: [PFASGroups repository]
- Email: luc@miaz.ch

## License

This work is licensed under a Creative Commons Attribution-NoDerivatives 4.0 International License.
Contact for exceptions to the No Derivatives term.

## Citation

If you use PFASgroups in your research, please cite:
[Citation information to be added]

---

## Acknowledgments

This project is part of the [ZeroPM project](https://zeropm.eu/) (WP2) and has received funding from the European Union's Horizon 2020 research and innovation programme under grant agreement No 101036756.

Developed at the Department of Environmental Science at Stockholm University.
