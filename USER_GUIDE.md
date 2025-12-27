# PFASgroups User Guide

Complete guide for using the PFASgroups package to parse, analyze, and fingerprint Per- and Polyfluoroalkyl Substances (PFAS).

## Table of Contents

1. [Installation](#installation)
2. [Quick Start](#quick-start)
3. [Python API](#python-api)
4. [Command-Line Interface](#command-line-interface)
5. [Custom Configuration](#custom-configuration)
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

# List available PFAS groups
pfasgroups list-groups
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
matches = parse_PFAS_groups(mol, formula)

# Plot
plot_pfasgroups([(smiles, matches)])
```

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
# List all groups
pfasgroups list-groups

# Save to file
pfasgroups list-groups --output groups.json
```

### Validate Config Command

Validate custom configuration files:

```bash
# Validate default config
pfasgroups validate-config

# Validate custom files
pfasgroups validate-config --groups-file custom_groups.json --fpaths-file custom_fpaths.json
```

---

## Custom Configuration

PFASgroups allows you to use custom definitions for pathtype patterns and PFAS groups.

### Configuration Files

Two JSON files control the behavior:

1. **fpaths.json** - Defines pathtype SMARTS patterns for chain detection
2. **PFAS_groups_smarts.json** - Defines PFAS group structures and constraints

### Method 1: Environment Variables

Set environment variables to use custom files globally:

**Windows (PowerShell):**
```powershell
$env:PFASGROUPS_FPATHS = "C:\path\to\custom_fpaths.json"
$env:PFASGROUPS_GROUPS = "C:\path\to\custom_groups.json"
```

**Linux/Mac:**
```bash
export PFASGROUPS_FPATHS="/path/to/custom_fpaths.json"
export PFASGROUPS_GROUPS="/path/to/custom_groups.json"
```

Then use PFASgroups normally:
```python
from PFASgroups import parse_pfas
results = parse_pfas(["C(C(F)(F)F)F"])  # Uses custom files
```

### Method 2: Python API - Session Default

Set custom configuration for your Python session:

```python
from PFASgroups import set_default_config, parse_pfas

# Set custom defaults
set_default_config(
    fpaths_file='custom_fpaths.json',
    groups_file='custom_groups.json'
)

# All subsequent calls use custom files
results = parse_pfas(["C(C(F)(F)F)F"])
fingerprints, info = generate_pfas_fingerprint(["C(C(F)(F)F)F"])
```

### Method 3: Python API - Per-Function

Pass custom files to individual function calls:

```python
from PFASgroups import parse_pfas, generate_pfas_fingerprint

# Use custom groups for specific call
results = parse_pfas(
    ["C(C(F)(F)F)F"],
    groups_file='custom_groups.json'
)

# Use custom configuration object
from PFASgroups import PFASConfig

config = PFASConfig(
    fpaths_file='custom_fpaths.json',
    groups_file='custom_groups.json'
)

results = parse_pfas(["C(C(F)(F)F)F"], config=config)
```

### Method 4: Command Line

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

### Validation

Always validate custom configuration files:

```bash
pfasgroups validate-config \
    --fpaths-file custom_fpaths.json \
    --groups-file custom_groups.json
```

Or in Python:

```python
from PFASgroups import PFASConfig

try:
    config = PFASConfig(
        fpaths_file='custom_fpaths.json',
        groups_file='custom_groups.json'
    )
    fpaths = config.load_fpaths()
    groups = config.load_pfas_groups()
    print(f"✓ Configuration valid: {len(groups)} groups loaded")
except Exception as e:
    print(f"✗ Configuration error: {e}")
```

---

## Advanced Usage

### Working with Configuration Objects

```python
from PFASgroups import PFASConfig, get_config

# Create custom config
config = PFASConfig(
    fpaths_file='custom_fpaths.json',
    groups_file='custom_groups.json'
)

# Load and inspect
groups = config.load_pfas_groups()
paths = config.load_smarts_paths()

print(f"Loaded {len(groups)} groups")
print(f"Available paths: {list(paths.keys())}")

# Use with specific path names
selected_paths = config.load_smarts_paths(
    path_names=['Perfluoroalkyl', 'Polyfluoroalkyl']
)

# Reset cache if files change
config.reset_cache()
```

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

---

## API Reference

### Core Functions

#### `parse_pfas(smiles_list, bycomponent=False, config=None, groups_file=None)`

Parse SMILES strings and identify PFAS groups.

**Parameters:**
- `smiles_list` (list): List of SMILES strings
- `bycomponent` (bool): Use component-based analysis
- `config` (PFASConfig): Custom configuration
- `groups_file` (str): Path to custom groups file

**Returns:**
- `list`: Match results for each SMILES

#### `generate_pfas_fingerprint(smiles, selected_groups=None, representation='vector', count_mode='binary', pfas_groups=None, config=None, groups_file=None)`

Generate PFAS fingerprints.

**Parameters:**
- `smiles` (str or list): SMILES string(s)
- `selected_groups` (list/range): Group indices to include
- `representation` (str): 'vector', 'dict', 'sparse', 'detailed', or 'int'
- `count_mode` (str): 'binary', 'count', or 'max_chain'
- `pfas_groups` (list): Custom PFASGroup objects
- `config` (PFASConfig): Custom configuration
- `groups_file` (str): Path to custom groups file

**Returns:**
- `fingerprints`: Fingerprint(s) in specified format
- `group_info` (dict): Information about groups used

### Configuration Functions

#### `set_default_config(fpaths_file=None, groups_file=None)`

Set global default configuration.

#### `get_config(fpaths_file=None, groups_file=None)`

Get a configuration instance.

#### `load_pfas_groups(groups_file=None)`

Load PFAS groups from file.

#### `load_smarts_paths(fpaths_file=None, path_names=None)`

Load SMARTS paths for chain detection.

### Classes

#### `PFASConfig(fpaths_file=None, groups_file=None)`

Configuration manager for custom files.

**Methods:**
- `load_fpaths()`: Load pathtype definitions
- `load_groups_raw()`: Load raw group data
- `load_pfas_groups()`: Load PFASGroup objects
- `load_smarts_paths(path_names=None)`: Load preprocessed SMARTS paths
- `reset_cache()`: Clear cached data

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

### Example 2: Custom Groups for Specific Use Case

```python
from PFASgroups import set_default_config, parse_pfas
import json

# Create minimal custom groups file
custom_groups = [
    {
        "id": 1,
        "name": "Short-chain perfluorinated",
        "smarts1": "[C](F)(F)F",
        "smarts2": "C(=O)O",
        "smartsPath": "Perfluoroalkyl",
        "constraints": {"nC": [2, 6]}  # C2-C6 only
    }
]

with open('short_chain_only.json', 'w') as f:
    json.dump(custom_groups, f)

# Use custom groups
set_default_config(groups_file='short_chain_only.json')
results = parse_pfas(["FC(F)(F)C(F)(F)C(=O)O"])  # Will match
results = parse_pfas(["FC(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(=O)O"])  # Won't match (C8)
```

### Example 3: Screening Large Datasets

```bash
# Prepare data file
echo "C(C(F)(F)F)F" > compounds.txt
echo "FC(F)(F)C(F)(F)C(=O)O" >> compounds.txt
echo "CCCCCC" >> compounds.txt  # Non-PFAS

# Screen compounds
pfasgroups parse --input compounds.txt --output results.json --pretty

# Generate fingerprints
pfasgroups fingerprint --input compounds.txt --groups 28-52 --output fps.json
```

---

## Troubleshooting

### Configuration File Not Found

```
FileNotFoundError: fpaths.json not found at: /path/to/file.json
```

**Solution:** Check file path and set correctly:
```python
from PFASgroups import set_default_config
set_default_config(fpaths_file='/correct/path/to/fpaths.json')
```

### Invalid JSON

```
json.JSONDecodeError: Expecting property name enclosed in double quotes
```

**Solution:** Validate JSON syntax using online validators or:
```bash
python -m json.tool your_file.json
```

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
