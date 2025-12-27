# PFASgroups Quick Reference

## Installation

```bash
pip install -e .
```

## Basic Usage

### Python API

```python
from PFASgroups import parse_pfas, generate_pfas_fingerprint

# Parse PFAS
results = parse_pfas(["C(C(F)(F)F)F", "FC(F)(F)C(F)(F)C(=O)O"])

# Generate fingerprints
fps, info = generate_pfas_fingerprint(["C(C(F)(F)F)F"])
```

### Command Line

```bash
# Parse
pfasgroups parse "C(C(F)(F)F)F"

# Fingerprint
pfasgroups fingerprint "C(C(F)(F)F)F" --format dict
```

## Custom Configuration

### Method 1: Helper Functions

```python
from PFASgroups import get_smartsPaths, get_PFASGroups, parse_pfas

# Load custom files
custom_paths = get_smartsPaths(filename='custom_fpaths.json')
custom_groups = get_PFASGroups(filename='custom_groups.json')

# Use in parsing
results = parse_pfas(
    ["C(C(F)(F)F)F"],
    smartsPaths=custom_paths,
    pfas_groups=custom_groups
)
```

### Method 1b: Extend Defaults

```python
from PFASgroups import get_PFASGroups, parse_pfas, PFASGroup

# Get defaults and add custom group
groups = get_PFASGroups()
groups.append(PFASGroup(
    id=999,
    name="Custom",
    smarts1="[C](F)(F)F",
    smarts2="[N+](=O)[O-]",
    smartsPath="Perfluoroalkyl",
    constraints={"nF": [3, None]}
))

# Use extended configuration
results = parse_pfas(["FC(F)(F)C(F)(F)[N+](=O)[O-]"], pfas_groups=groups)
```

### Method 1c: Filter Defaults

```python
from PFASgroups import get_PFASGroups, parse_pfas

# Get defaults and filter
all_groups = get_PFASGroups()
perfluoro_only = [g for g in all_groups if g.smartsPath == 'Perfluoroalkyl']

# Use filtered configuration
results = parse_pfas(["FC(F)(F)C(F)(F)C(=O)O"], pfas_groups=perfluoro_only)
```

### Method 2: Manual Loading

```python
from PFASgroups import PFASGroup, parse_pfas
import json

# Load and construct groups
with open('custom_groups.json') as f:
    data = json.load(f)
groups = [PFASGroup(**g) for g in data]

# Use directly
results = parse_pfas(["C(C(F)(F)F)F"], pfas_groups=groups)
```

### Method 3: Command Line

```bash
pfasgroups parse --groups-file custom.json "C(C(F)(F)F)F"
pfasgroups fingerprint --fpaths-file custom_fpaths.json --input smiles.txt
```

## Configuration File Formats

### fpaths.json

```json
{
  "Perfluoroalkyl": {
    "chain": "[C;X4;H0](F)(F)!@!=!#[C;X4;H0](F)(F)",
    "end": "[C;X4;H0](F)(F)F"
  }
}
```

### PFAS_groups_smarts.json

```json
[
  {
    "id": 1,
    "name": "My Group",
    "smarts1": "[C](F)(F)F",
    "smarts2": "C(=O)O",
    "smartsPath": "Perfluoroalkyl",
    "constraints": {"nC": [2, 6]}
  }
]
```

## Helper Functions

### get_smartsPaths(filename=None)

Load path SMARTS patterns. Returns dict of `{name: [chain_mol, end_mol]}`.

```python
# Get defaults
paths = get_smartsPaths()

# Get from file  
paths = get_smartsPaths(filename='custom.json')

# Extend defaults (see compile_smartsPath below for easier way)
paths = get_smartsPaths()
from rdkit import Chem
custom = Chem.MolFromSmarts("[C](Cl)")
custom.UpdatePropertyCache()
Chem.GetSymmSSSR(custom)
custom.GetRingInfo().NumRings()
paths['Custom'] = [custom, custom]
```

### get_PFASGroups(filename=None)

Load PFAS group definitions. Returns list of PFASGroup objects.

```python
# Get defaults
groups = get_PFASGroups()

# Get from file
groups = get_PFASGroups(filename='custom.json')

# Extend defaults
groups = get_PFASGroups()
from PFASgroups import PFASGroup
groups.append(PFASGroup(id=999, name="Custom", smarts1="[C](F)", smarts2="CO", 
                        smartsPath="Perfluoroalkyl", constraints={"nF": [1, None]}))
```

### compile_smartsPath(chain_smarts, end_smarts)

Compile a single path type from SMARTS strings. Returns `[chain_mol, end_mol]`.

```python
from PFASgroups import compile_smartsPath, get_smartsPaths

# Get defaults
paths = get_smartsPaths()

# Add custom path (much easier than manual preprocessing!)
paths['Perchlorinated'] = compile_smartsPath(
    "[C;X4](Cl)(Cl)!@!=!#[C;X4](Cl)(Cl)",  # chain
    "[C;X4](Cl)(Cl)Cl"                     # end
)

# Use in parsing
from PFASgroups import parse_pfas
results = parse_pfas(smiles, smartsPaths=paths)
```

### compile_smartsPaths(paths_dict)

Compile multiple path types at once. Returns dict of compiled paths.

```python
from PFASgroups import compile_smartsPaths, get_smartsPaths

# Define multiple custom paths
custom_paths = {
    'Perchlorinated': {
        'chain': '[C;X4](Cl)(Cl)!@!=!#[C;X4](Cl)(Cl)',
        'end': '[C;X4](Cl)(Cl)Cl'
    },
    'Perbrominated': {
        'chain': '[C;X4](Br)(Br)!@!=!#[C;X4](Br)(Br)',
        'end': '[C;X4](Br)(Br)Br'
    }
}

# Compile all at once
compiled = compile_smartsPaths(custom_paths)

# Merge with defaults
paths = get_smartsPaths()
paths.update(compiled)
```
custom = Chem.MolFromSmarts("[C](F)")
custom.UpdatePropertyCache()
Chem.GetSymmSSSR(custom)
custom.GetRingInfo().NumRings()
paths['MyPath'] = [custom, custom]
```

### get_PFASGroups(filename=None)

Load PFAS groups. Returns list of PFASGroup objects.

```python
# Get defaults
groups = get_PFASGroups()

# Get from file
groups = get_PFASGroups(filename='custom.json')

# Extend defaults
from PFASgroups import PFASGroup
groups = get_PFASGroups()
groups.append(PFASGroup(id=999, name="Custom", ...))

# Filter defaults
groups = get_PFASGroups()
filtered = [g for g in groups if g.id >= 28]
```

## Common Tasks

### Extend Defaults with Custom Groups

```python
from PFASgroups import get_PFASGroups, PFASGroup, parse_pfas

# Get defaults
groups = get_PFASGroups()

# Add organization-specific group
groups.append(PFASGroup(
    id=1000,
    name="Short-chain PFCA",
    smarts1="[C](F)(F)F",
    smarts2="C(=O)O",
    smartsPath="Perfluoroalkyl",
    constraints={"nC": [2, 6]}
))

results = parse_pfas(smiles_list, pfas_groups=groups)
```

```python
from PFASgroups import PFASGroup, parse_pfas

# Define short-chain only
groups = [PFASGroup(
    id=1,
    name="Short-chain",
    smarts1="[C](F)(F)F",
    smarts2="C(=O)O",
    smartsPath="Perfluoroalkyl",
    constraints={"nC": [2, 6]}
)]

results = parse_pfas(["FC(F)(F)C(F)(F)C(=O)O"], pfas_groups=groups)
```

### Batch Processing

```python
import pandas as pd
from PFASgroups import generate_pfas_fingerprint

df = pd.read_csv('compounds.csv')
fps, info = generate_pfas_fingerprint(df['smiles'].tolist())
```

### Machine Learning

```python
from sklearn.ensemble import RandomForestClassifier
from PFASgroups import generate_pfas_fingerprint

X_train, info = generate_pfas_fingerprint(train_smiles, representation='vector')
clf = RandomForestClassifier()
clf.fit(X_train, labels)
```

## API Quick Reference

| Function | Purpose | Can Extend? |
|----------|---------|-------------|
| `parse_pfas(smiles, **kwargs)` | Parse SMILES list | ✓ Pass custom kwargs |
| `parse_PFAS_groups(mol, **kwargs)` | Parse single molecule | ✓ Pass custom kwargs |
| `generate_pfas_fingerprint(smiles, **kwargs)` | Generate fingerprints | ✓ Pass custom kwargs |
| `get_smartsPaths(filename=None)` | Load path SMARTS | ✓ Returns modifiable dict |
| `get_PFASGroups(filename=None)` | Load PFAS groups | ✓ Returns modifiable list |

## Command Line Quick Reference

```bash
# Parse commands
pfasgroups parse "SMILES"
pfasgroups parse --input file.txt --output results.json
pfasgroups parse --groups-file custom.json "SMILES"

# Fingerprint commands  
pfasgroups fingerprint "SMILES"
pfasgroups fingerprint --groups 28-52 --format dict "SMILES"
pfasgroups fingerprint --input file.txt --output fp.json

# List available groups and paths
pfasgroups list-groups
pfasgroups list-paths
pfasgroups list-groups --groups-file custom.json

# Validate config
pfasgroups validate-config --groups-file custom.json
```

## Parameters

### parse_pfas / parse_PFAS_groups

- `pfas_groups`: List of PFASGroup objects
- `smartsPaths`: Dict of path SMARTS patterns
- `bycomponent`: Boolean for component-based analysis
- `formula`: Molecular formula (for parse_PFAS_groups)

### generate_pfas_fingerprint

- `selected_groups`: List or range of group indices
- `representation`: 'vector', 'dict', 'sparse', 'detailed', 'int'
- `count_mode`: 'binary', 'count', 'max_chain'
- `pfas_groups`: Custom PFASGroup objects

---

See [USER_GUIDE.md](USER_GUIDE.md) for complete documentation.
