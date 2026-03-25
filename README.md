# PFASGroups

A comprehensive cheminformatics package for automated detection, classification, and analysis of halogenated substances, in particular per- and polyfluoroalkyl substances (PFAS) in chemical databases.

## Overview

PFASGroups combines SMARTS pattern matching, molecular formula constraints, and graph-based pathfinding (using RDKit and NetworkX) to identify and classify PFAS compounds. The package enables systematic PFAS universe mapping and environmental monitoring applications.

## Key Features

### Core Capabilities
- **Halogen Group Identification**: Automated detection of 119 functional groups (114 compiled for fluorine-only embedding):
  - 27 PFAS OECD groups
  - 48 generic functional groups (IDs 29–73, 117–119; the last 3 are halogen-context-dependent or recently added)
  - 43 fluorotelomer-specific groups with CH₂ linker validation
  - 1 aggregate pattern-matching group (Group 116: Telomers, `compute=False`)
- **Atom Reference Requirement**: For non-telomer groups, SMARTS patterns must match atoms that are part of or directly connected to the fluorinated component (per/polyfluorinated carbons), respecting the `max_dist_from_comp` constraint
- **Linker Validation**: CH₂-specific validation for 40 fluorotelomer groups to distinguish from direct-attachment analogues. Telomer groups use `linker_smarts` to allow functional groups separated from perfluoro chains by non-fluorinated linkers
- **Aggregate Groups**: Pattern-matching groups that collect related PFAS groups via regex (e.g., Group 113 matches all "telomer" groups)
- **Component Length Analysis**: Quantification of per- and polyfluorinated alkyl components with CF₂ unit counting
- **Graph Metrics**: Comprehensive structural characterization (branching, eccentricity, diameter, resistance, centrality, telomer spacer length, ring size)
- **Customizable Definitions**: Easy extension to additional PFAS groups and halogenated chemical classes via JSON configuration. Component type definitions in `data/component_smarts_halogens.json` support optional per-type `constraints` (e.g. `{"gte": {"F": 2}}`) that are evaluated against the component's full atom count (backbone carbons plus attached halogens) to enforce minimum halogen counts or element exclusions — without repeating the check in every group definition

### Additional Tools
- **Homologue Series Generation**: Iterative component shortening to explore theoretical chemical space
- **Fingerprint Generation**: PFAS fingerprints for machine learning applications
- **Visualization**: Assign and visualize PFAS groupings
- **Multiple Interfaces**: Python API, command-line tool, and browser-based JavaScript version (RDKitJS)
- **Batch Processing**: Efficient analysis of large chemical databases

## Installation

### From Pypi

```sh
pip install PFASGroups
```

### From source (recommended for development)

**Prerequisites**: Python ≥ 3.7, RDKit (install via conda or pip before the steps below).

```sh
# Clone the repository
git clone https://github.com/lucmiaz/PFASGroups.git
cd PFASGroups

# Install in editable mode (development install)
pip install -e .
```

> **Note for conda users**: install RDKit via conda before running pip:
> ```sh
> conda install -c conda-forge rdkit
> ```

After installation, both `HalogenGroups` (all halogens by default) and `PFASgroups` (fluorine only by default) are importable, and the `pfasgroups` CLI command is available in your terminal.

### Verify installation

```python
from PFASGroups import parse_smiles
results = parse_smiles("FC(F)(F)C(F)(F)C(=O)O")   # → PFASEmbeddingSet
print(results)        # prints PFASEmbeddingSet summary (molecule count, matched groups, …)
print(results[0])     # prints PFASEmbedding summary for the first molecule
```

## Repository Structure

```
PFASGroups/
├── HalogenGroups/                   # Multi-halogen wrapper package
│   └── __init__.py                  #   Wraps PFASgroups; defaults halogens=['F','Cl','Br','I']
├── PFASGroups/                      # Core implementation package (also importable as PFASgroups)
│   ├── __init__.py                  #   Public API
│   ├── core.py                      #   SMARTS matching engine, component detection, decorators
│   ├── parser.py                    #   parse_smiles / parse_mols entry points
│   ├── embeddings.py                #   presets
│   ├── PFASEmbeddings.py            #   PFASEmbeddingSet, PFASEmbedding
│   ├── HalogenGroupModel.py         #   HalogenGroup data model
│   ├── PFASDefinitionModel.py       #   PFASDefinition model
│   ├── ComponentsSolverModel.py     #   Graph-based path-finding solver
│   ├── getter.py                    #   get_HalogenGroups, get_componentSMARTSs, …
│   ├── cli.py                       #   Command-line interface (halogengroups / pfasgroups)
│   ├── prioritise.py                #   prioritise_molecules
│   ├── generate_homologues.py       #   Homologue series generation
│   ├── fragmentation.py             #   Fragment utilities
│   ├── draw_mols.py                 #   Molecule and group visualisation helpers
│   └── data/                        #   JSON configuration files (bundled with the package)
│       ├── Halogen_groups_smarts.json      # 113 halogen group definitions
│       ├── component_smarts.json           # Fluorinated component SMARTS patterns
│       ├── component_smarts_halogens.json  # Multi-halogen component SMARTS patterns
│       └── PFAS_definitions_smarts.json    # PFAS regulatory definitions
├── tests/                           # Pytest test suite (25+ test modules)
│   ├── test_halogen_groups_smarts.py
│   ├── test_results_fingerprint.py
│   ├── test_results_model.py
│   ├── test_prioritise.py
│   └── …
├── examples/                        # Standalone usage example scripts
│   ├── results_fingerprint_analysis.py
│   ├── prioritization_examples.py
│   └── …
├── docs/                            # Sphinx documentation source
│   ├── quickstart.rst
│   ├── algorithm.rst
│   ├── customization.rst
│   ├── ResultsFingerprint_Guide.md
│   └── …
├── benchmark/                       # Benchmarking scripts and timing reports
│   ├── data/
│   ├── reports/
│   └── …
└── pyproject.toml                   # Package metadata and build configuration
```

**Two importable packages, one codebase:**

| Package | Default `halogens` | Typical use |
|---|---|---|
| `HalogenGroups` | `['F', 'Cl', 'Br', 'I']` (all) | Multi-halogen analysis |
| `PFASgroups` | `'F'` (fluorine only) | PFAS / fluorine-focused analysis |

## Benchmark Summary (Feb 2026, using v. 2.2.3 only on F groups)

Benchmarks were run on an Intel(R) Core(TM) i7-7700HQ CPU @ 2.80GHz (4C/8T), 15.5 GB RAM, Python 3.9.23, RDKit 2025.09.2, NetworkX 3.2.1 (the old version of Python was taken for compatibility with PFAS-atlas).

| Dataset/Profile | Count | Atom range | PFASgroups mean/median (ms) | PFAS-Atlas mean/median (ms) | Relative speed | Notes |
| --- | --- | --- | --- | --- | --- | --- |
| OECD reference (real compounds) | 3,414 | Typical small/medium | 19.2 / 14.8 | 38.8 / 37.7 | 2.02x faster | Real-world dataset representing existing compounds. |
| Timing stress-test (full metrics) | 2,500 | 11-625 | 251.8 / 24.4 | 58.5 / 34.2 | 0.23x | Synthetic stress-test with large molecules; heavy-tail runtime. |
| Timing stress-test (no resistance) | 2,500 | 11-619 | 176.7 / 24.8 | N/A | 1.43x faster vs full | Disables effective graph resistance only. |
| Timing stress-test (no metrics) | 2,500 | 11-619 | 97.7 / 19.7 | N/A | 2.58x faster vs full | Disables all component graph metrics. |

Timing profile plots (full vs no resistance vs no metrics):
- [benchmark/reports/timing_profiles_comparison.png](benchmark/reports/timing_profiles_comparison.png)
- [benchmark/reports/timing_profiles_residuals.png](benchmark/reports/timing_profiles_residuals.png)

Disable or limit graph metrics in the Python API:

```python
from HalogenGroups import parse_smiles

# Skip all component graph metrics (fastest)
parse_smiles(smiles_list, compute_component_metrics=False)

# Keep metrics but skip effective graph resistance entirely
parse_smiles(smiles_list, limit_effective_graph_resistance=0)

# Compute resistance only for components below a size threshold
parse_smiles(smiles_list, limit_effective_graph_resistance=200)
```

CLI equivalents:

```bash
# Skip all component graph metrics (fastest)
pfasgroups parse --no-component-metrics "C(C(F)(F)F)F"

# Skip effective graph resistance entirely
pfasgroups parse --limit-effective-graph-resistance 0 "C(C(F)(F)F)F"

# Compute resistance only for components below a size threshold
pfasgroups parse --limit-effective-graph-resistance 200 "C(C(F)(F)F)F"
```

## Quick Start

### Python API

```python
from PFASGroups import parse_smiles

# Parse PFAS structures — returns a PFASEmbeddingSet
smiles_list = ["C(C(F)(F)F)F", "FC(F)(F)C(F)(F)C(=O)O"]
results = parse_smiles(smiles_list)           # → PFASEmbeddingSet
print(results)                                # summary: molecule count, top groups, …
print(results[0])                             # per-molecule summary for first entry

# Embedding from a pre-parsed set (avoids re-parsing)
arr  = results.to_array()                     # default: binary, all groups
arr  = results.to_array(group_selection='oecd', component_metrics=['binary'])
cols = results.column_names()                 # matching column labels

# Filter components by halogen, form, and saturation
results_f      = parse_smiles(smiles_list, halogens='F')  # Fluorine only
results_pfa    = parse_smiles(smiles_list, halogens='F', saturation='per', form='alkyl')  # Perfluoroalkyl only
results_cyclic = parse_smiles(smiles_list, form='cyclic')  # Cyclic forms only

# Dimensionality reduction (methods on PFASEmbeddingSet)
pca_result  = results.perform_pca(n_components=5, plot=True)
tsne_result = results.perform_tsne(perplexity=30, plot=True)
umap_result = results.perform_umap(n_neighbors=15, plot=True)

# Compare two datasets using KL divergence
other_results = parse_smiles(other_smiles_list)
similarity = results.compare_kld(other_results, method='minmax')

# Save to SQL database
results.to_sql(filename='results.db')

# New in v2.2.4: Prioritization tool for screening and ranking
from HalogenGroups import prioritise_molecules

# Prioritize by similarity to reference list (e.g., known persistent PFAS)
reference = ["FC(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(=O)O"]
prioritized = prioritise_molecules(
    molecules=smiles_list,
    reference=reference,
    return_scores=True
)

# Prioritize by intrinsic fluorination properties (long chains, high F content)
prioritized = prioritise_molecules(
    molecules=smiles_list,
    a=1.0,  # Weight for total fluorination
    b=5.0,  # Weight for longest chains
    percentile=90  # Focus on 90th percentile of component sizes
)
```

### Command Line

```bash
# Parse SMILES strings
halogengroups parse "C(C(F)(F)F)F" "FC(F)(F)C(F)(F)C(=O)O"

# Generate fingerprints
halogengroups fingerprint "C(C(F)(F)F)F" --format dict

# List available PFAS groups
halogengroups list-groups
```

### Filtering Components by Halogen, Form, and Saturation

Filter component matches by specific halogens, chemical forms, or saturation levels:

```python
from HalogenGroups import parse_smiles

smiles_list = ["C(C(F)(F)F)F", "FC(F)(F)C(F)(F)C(=O)O", "C(C(Cl)(Cl)Cl)Cl"]

# Filter only fluorine components
results_f = parse_smiles(smiles_list, halogens='F')

# Filter perfluorinated alkyl compounds
results_pfa = parse_smiles(smiles_list, halogens='F', saturation='per', form='alkyl')

# Filter polyfluorinated cyclic compounds
results_polyf_cyclic = parse_smiles(smiles_list, halogens='F', saturation='poly', form='cyclic')

# Filter multiple halogens (F and Cl)
results_multi = parse_smiles(smiles_list, halogens=['F', 'Cl'])

# Valid filter options:
# - halogens: 'F', 'Cl', 'Br', 'I', or list like ['F', 'Cl']
# - saturation: 'per' or 'poly' (or list like ['per', 'poly'] for both)
# - form: 'alkyl' or 'cyclic' (or list like ['alkyl', 'cyclic'] for both)
```

## Embedding with Graph Metrics

The `to_array()` / `to_fingerprint()` methods accept a `component_metrics` list that
stacks one block of *N_G* columns per metric (default `N_G = 114` for fluorine-only).

```python
import numpy as np
from PFASGroups import parse_smiles, EMBEDDING_PRESETS

smiles = [
    "O=C(O)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)F",  # PFOA
    "O=C(O)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)F",          # PFHpA
    "OCC(F)(F)C(F)(F)C(F)(F)C(F)(F)F",                              # 4:2 FTOH
    "OCC(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)F",               # 6:2 FTOH
]
results = parse_smiles(smiles)

# --- binary embedding (default) ---
arr_bin = results.to_array()                    # shape (4, 114)

# --- preset 'best': binary + effective_graph_resistance ---
# Best discrimination (mean Tanimoto 0.184, outperforms TxP-PFAS 129-bit)
print(EMBEDDING_PRESETS['best']['description'])
arr_best = results.to_array(preset='best')      # shape (4, 228) = 114 × 2

# --- effective graph resistance directly ---
arr_egr = results.to_array(
    component_metrics=['binary', 'effective_graph_resistance']
)                                               # shape (4, 228)

# --- n_spacer: telomer CH₂ spacer length (encodes 'm' in 'm:n' notation) ---
# Zero for non-telomers; distinguishes 2:1, 4:2, 6:2 FTOHs in a single column
arr_ns = results.to_array(component_metrics=['n_spacer'])  # shape (4, 114)

# --- ring_size: smallest ring overlapping each matched component ---
# Zero for acyclic groups; 5 for azoles/furans; 6 for benzene/cyclohexane
arr_rs = results.to_array(component_metrics=['ring_size'])  # shape (4, 114)

# --- combined: EGR + n_spacer + ring_size + molecule-wide descriptors ---
arr_combined = results.to_array(
    component_metrics=['binary', 'effective_graph_resistance',
                       'n_spacer', 'ring_size'],
    molecule_metrics=['n_components', 'max_size',
                      'mean_branching', 'max_branching',
                      'mean_component_fraction', 'max_component_fraction'],
)  # shape (4, 4*114 + 6) = (4, 462)

# --- multi-halogen: one parse per halogen, then hstack ---
arrs = []
for hal in ['F', 'Cl', 'Br', 'I']:
    r_hal = parse_smiles(smiles, halogens=hal)
    arrs.append(r_hal.to_array(component_metrics=['effective_graph_resistance']))
arr_4x = np.hstack(arrs)                        # shape (4, 456) = 4 × 114

# --- column names ---
cols = results.column_names(
    component_metrics=['binary', 'effective_graph_resistance']
)
print(cols[:4])   # e.g. ['Perfluoroalkyl [binary]', ..., 'Perfluoroalkyl [EGR]', ...]
```

See [`examples/embedding_with_graph_metrics.py`](examples/embedding_with_graph_metrics.py)
for a complete runnable script covering all options.

CLI equivalents:

```bash
# Filter for fluorine components only
halogengroups parse --halogens F "C(C(F)(F)F)F"

# Filter for perfluorinated alkyl
halogengroups parse --halogens F --saturation per --form alkyl "FC(F)(F)C(F)(F)C(=O)O"

# Filter for multiple halogens
halogengroups parse --halogens F Cl "C(C(F)(F)F)Cl"

# Filter for cyclic forms only
halogengroups parse --form cyclic "your_smiles_here"

# Filter for polyfluorinated components
halogengroups parse --halogens F --saturation poly "FC(F)(F)C(C(F)(F)F)C(F)(F)F"
```

## Multi-Halogen Analysis

PFASGroups supports fluorine, chlorine, bromine, and iodine. There are two ways
to analyse all halogens at once:

### Option A – import `HalogenGroups` (all halogens by default)

```python
from HalogenGroups import parse_smiles

smiles_list = [
    "C(C(F)(F)F)F",          # fluorinated
    "ClC(Cl)(Cl)C(Cl)(Cl)Cl", # chlorinated
    "BrC(Br)(Br)CBr",         # brominated
]

# halogens defaults to ['F','Cl','Br','I'] — no extra argument needed
results = parse_smiles(smiles_list)           # → PFASEmbeddingSet

# to_array() reads the halogen info already captured during parsing
arr  = results.to_array(group_selection='oecd', component_metrics=['binary'])
cols = results.column_names(group_selection='oecd')

```

### Option B – import `PFASgroups` and specify `halogens` explicitly

```python
from PFASgroups import parse_smiles

smiles_list = [
    "C(C(F)(F)F)F",
    "ClC(Cl)(Cl)C(Cl)(Cl)Cl",
    "BrC(Br)(Br)CBr",
]

# Explicitly pass all halogens
results = parse_smiles(smiles_list, halogens=['F', 'Cl', 'Br', 'I'])  # → PFASEmbeddingSet

arr  = results.to_array(group_selection='oecd', component_metrics=['binary'])
cols = results.column_names(group_selection='oecd')
```

Both approaches produce identical results. Use `HalogenGroups` when your workflow
is always multi-halogen; use `PFASgroups` with an explicit `halogens` argument when
you want to mix fluorine-only and multi-halogen calls in the same script.

CLI equivalents (the CLI always requires an explicit `--halogens` flag):

```bash
# All four halogens
halogengroups parse --halogens F Cl Br I "C(C(F)(F)F)F" "ClC(Cl)(Cl)C(Cl)(Cl)Cl"

# Fingerprint with all halogens
halogengroups fingerprint --halogens F Cl Br I "C(C(F)(F)F)F"
```

## Custom Configuration

Use custom pathtype definitions and PFAS groups:

```python
# Load custom files entirely
from HalogenGroups import get_componentSMARTSs, get_HalogenGroups, parse_smiles

custom_paths = get_componentSMARTSs(filename='my_component_smartss.json')
custom_groups = get_HalogenGroups(filename='my_groups.json')

results = parse_smiles(
    ["C(C(F)(F)F)F"],
    componentSmartss=custom_paths,
    pfas_groups=custom_groups
)
```

```python
# Or extend defaults with your custom groups
from HalogenGroups import get_HalogenGroups, HalogenGroup, parse_smiles, compile_componentSmarts, get_componentSMARTSs

# Add custom PFAS groups
groups = get_HalogenGroups()  # Get defaults
groups.append(HalogenGroup(
    id=999,
    name="My Custom Group",
    smarts={"[C](F)(F)F": 1},
    componentSmarts=None,
    componentSaturation="per",
    componentHalogens="F",
    componentForm="alkyl",
    constraints={"gte": {"F": 1}}
))

results = parse_smiles(["FC(F)(F)C(F)(F)[N+](=O)[O-]"], pfas_groups=groups)

# Custom max_dist_from_comp parameter
# For functional groups without formula constraints,
# max_dist_from_comp limits the maximum bond distance between
# a functional group match and a fluorinated carbon terminal atom (default: 0)
groups.append(HalogenGroup(
    id=998,
    name="Extended Distance Group",
    smarts={"[#6$([#6][OH1])]": 1},
    componentSmarts=None,
    constraints={},
    max_dist_from_comp=3  # Allow up to 3 bonds from fluorinated carbon
))

# Add custom path types (e.g., chlorinated analogs)
paths = get_componentSMARTSs()
paths['Perchlorinated'] = compile_componentSmarts(
    "[C;X4](Cl)(Cl)!@!=!#[C;X4](Cl)(Cl)",  # component pattern
    "[C;X4](Cl)(Cl)Cl"                     # end pattern
)

results = parse_smiles(["ClC(Cl)(Cl)C(Cl)(Cl)C(=O)O"], componentSmartss=paths)
```

```bash
# Via command line
halogengroups parse --groups-file my_custom_groups.json "C(C(F)(F)F)F"

# List available groups and paths
halogengroups list-groups
halogengroups list-paths
```

## Documentation

- **[USER_GUIDE.md](USER_GUIDE.md)** - Complete documentation with examples
- **[QUICK_REFERENCE.md](QUICK_REFERENCE.md)** - Quick reference for common tasks

## Usage Examples

See [USER_GUIDE.md](USER_GUIDE.md) for comprehensive examples including:
- Basic PFAS parsing and analysis
- Fingerprint generation for machine learning
- Custom configuration files
- Batch processing
- Integration with pandas and scikit-learn

## Summary of changes by version

- **Version 3.2.2**: Fixed polyhalogenated alkyl components matching less than 2 halogens. Added option to pass 'halogens' to formula constraints. Added options to include component-wide formula constraints (on dist-1 neighbours from matched C-only components).

- **Version 3.2.1**: Added n-spacer metric for telomers and ring size for aryl and cyclic groups. These can be used in the embeddings.

- **Version 3.2.0**: Use of BDE for computing graph resistance.

- **Version 3.1.4**: Changed fingerprint parameters to molecule wide and component wide, with pre-configuration for best combinations.

- **Version 3.1.0**: Added support for other halogens, changed names to be more generic (with some support for backward compatibility). Added component smarts for other halogens, cyclic and alkyl components.

- **Version 2.2.4 (Feb 2026)**: Advanced fingerprint analysis with dimensionality reduction (PCA, kernel-PCA, t-SNE, UMAP), KL divergence comparison for dataset similarity assessment, and SQL persistence for results. Added molecule prioritization tool for screening applications, ranking by similarity to reference lists or by intrinsic fluorination properties. Introduced `PFASEmbeddingSet` and `PFASEmbedding` with comprehensive analysis methods, automated plot generation, and extensive documentation.

- **Version 2.2.3 (Feb 2026)**: Added `PFASEmbeddingSet` container to offer easier plotting and summarising capabilities for results.

- **Version 2.2 (Feb 2026)**: Added linked_smarts option to specify a restriction on path between smarts groups and fluorinated component. Added new PFASgroups (telomers). **v2.2.2** Fixed telomers and added examples and counter-examples to each PFASgroup. Removed boundary O in fluorinated components (for both Per and Polyalkyl components).

- **Version 2.1 (Jan 2026)**: Added support for multiple smarts, with individual minimum count, per PFASgroup.

- **Version 2.0 (Jan 2026)**: Major expansion of graph‑based component metrics, new coverage statistics, schema updates, and richer per‑component outputs.

### Version 2.2.4 (February 2026) - Advanced Fingerprint Analysis

Major enhancement adding comprehensive dimensionality reduction and statistical comparison capabilities:

**New Features:**
- **PFASEmbeddingSet / PFASEmbedding**: Unified result container with flexible embedding generation
- **Dimensionality Reduction**: PCA, kernel-PCA, t-SNE, and UMAP with automatic plotting
- **Statistical Comparison**: KL divergence for comparing dataset compositions
- **Database Persistence**: SQL save/load for results
- **Molecule Prioritization**: Screening and ranking tool for PFAS datasets
- **Comprehensive Documentation**: Complete API reference, examples, and 100+ tests

**Key Methods:**
- `PFASEmbeddingSet.to_array()`: Generate numeric embedding matrix from parsed results
- `PFASEmbeddingSet.column_names()`: Column labels matching `to_array()` output
- `PFASEmbeddingSet.perform_pca()`: Principal Component Analysis
- `PFASEmbeddingSet.perform_kernel_pca()`: Non-linear kernel PCA
- `PFASEmbeddingSet.perform_tsne()`: t-SNE visualization
- `PFASEmbeddingSet.perform_umap()`: Fast UMAP dimensionality reduction
- `PFASEmbeddingSet.compare_kld()`: Dataset similarity via KL divergence
- `prioritise_molecules()`: Rank molecules by similarity or fluorination properties
- `PFASEmbeddingSet.to_sql()`: Persist results to a SQLite/PostgreSQL database

**Prioritization Strategies:**
- **Reference-based**: Rank by distributional similarity to known PFAS (e.g., persistent chemicals)
- **Intrinsic properties**: Score by fluorination characteristics (total F content, chain length)
- **Flexible weighting**: Tune parameters for different screening objectives

**Use Cases:**
- Exploratory data analysis of PFAS inventories
- Database comparison and compositional analysis
- Cluster identification and pattern recognition
- Machine learning preprocessing
- Chemical space visualization
- Priority screening for environmental monitoring
- Regulatory watchlist generation

See `docs/ResultsFingerprint_Guide.md`, `docs/prioritization.rst`, `examples/results_fingerprint_analysis.py`, and `examples/prioritization_examples.py` for details.

### Version 2.2.3 (February 2026) - PFASEmbeddingSet Container

**New Features:**
- `PFASEmbeddingSet` container with visualization helpers
- Enhanced component plotting utilities
- Improved documentation and examples

### Version 2.0 (January 2026) - Comprehensive Graph Metrics

Major enhancement adding comprehensive NetworkX graph theory metrics for detailed component analysis:

**New Features:**
- **Component-Level Metrics**: Each fluorinated component now includes 15+ graph metrics:
  - `diameter` and `radius` - Graph eccentricity bounds
  - `center`, `periphery`, `barycenter` - Structural node sets
  - `effective_graph_resistance` - Sum of resistance distances
  - `component_fraction` - Fraction of molecule covered by component (includes all attached H, F, Cl, Br, I)
  - Distance metrics from functional groups to structural features
- **Molecular Coverage Metrics**: New fraction-based metrics quantify fluorination extent:
  - `mean_component_fraction` - Average coverage per component
  - `total_components_fraction` - Total coverage by union of all components (accounts for overlaps)
- **Summary Statistics**: Aggregated metrics across all components per PFAS group
- **Enhanced Database Models**: New `Components` model stores individual component data with all metrics
- **Improved Analysis**: Better understanding of molecular topology, branching, functional group positioning, and fluorination extent

**Breaking Changes:**
- `parse_mols` output now includes additional summary metric fields (`mean_diameter`, `mean_radius`, etc.)
- Database schema changes require migration (see `DATABASE_MIGRATION_GUIDE.md`)

**Metrics Explained:**
- `branching` (0-1): Measures linearity (1.0 = linear, 0.0 = highly branched) - renamed from "eccentricity"
- `mean_eccentricity`, `median_eccentricity`: Graph-theoretic eccentricity statistics for component nodes
- `smarts_centrality` (0-1): Functional group position (1.0 = central, 0.0 = peripheral)
- `n_spacer` (int ≥ 0): Fluorotelomer CH₂ spacer length — the "m" in "m:n" telomer notation; 0 for all non-telomeric groups
- `ring_size` (int ≥ 0): Smallest ring overlapping the matched component; 0 for acyclic chains, 5 for azoles/furans, 6 for benzene/cyclohexane derivatives
- `component_fraction` (0-1): Fraction of total molecule atoms in this component (includes all attached atoms)
- `total_components_fraction` (0-1): Fraction of molecule covered by union of all components
- `diameter`: Maximum distance between any two atoms in component
- `radius`: Minimum eccentricity across component nodes
- `barycenter`: Nodes minimizing total distance to all other nodes
- `center`: Nodes with minimum eccentricity
- `periphery`: Nodes with maximum eccentricity

See `COMPREHENSIVE_METRICS_SUMMARY.md` for complete documentation.

- **Version 1.x**: Shift to component‑based analysis with improved SMARTS matching and better handling of branched/cyclic structures.

### Version 1.x - Component-Based Analysis

- Replaced chain-finding with connected component analysis
- Added support for branched and cyclic structures
- Improved SMARTS pattern matching for diverse PFAS classes

### Version 0.x - Path-Based Analysis

- Find SMARTS match connected to either a second SMARTS or a default path-related SMARTS using networkx shortest_path.

## Tests and benchmarks

### Automated pytest suite

The module ships a comprehensive pytest suite in `tests/`. Run it from the repository root:

```bash
pytest tests/            # run all tests
pytest tests/ -v         # verbose output
pytest tests/ -k smarts  # run only SMARTS-related tests
```

Key test files:

| File | What it tests |
|------|---------------|
| `tests/test_halogen_groups_smarts.py` | Per-group positive and negative examples embedded in the JSON configuration |
| `tests/test_pfasstructv5.py` | PFASSTRUCTv5 definition against the PFASSTRUCTV5 from DSSTox (see below) |
| `tests/test_results_fingerprint.py` | Fingerprint generation, group selection, and dimensionality-reduction methods |
| `tests/test_definition_comparison.py` | Consistency of the five PFAS regulatory definitions |
| `tests/test_linker_smarts.py` | CH₂ linker validation for fluorotelomer groups |
| `tests/test_metrics.py` | Graph-theoretical component metric values |

### Benchmarking against PFAS-Atlas

Performance and coverage were benchmarked against [PFAS-Atlas](https://github.com/su-group/PFAS-atlas?tab=readme-ov-file) (Su *et al.* 2024, [DOI: 10.1016/j.scitotenv.2024.171229](https://www.sciencedirect.com/science/article/pii/S0048969724013688)).

See the table in the [Benchmark Summary](#benchmark-summary-feb-2026-using-v-223-only-on-f-groups) section above for timing numbers. Full per-dataset comparison and Sankey diagrams are in `benchmark/reports/`.

### PFASSTRUCTv5 validation against Richard *et al.* 2023

The PFASSTRUCTv5 definition (definition ID 5, `'PFASSTRUCTv5'`, fluorineRatio ≥ 0.3) was validated against the PFASSTRUCTV5 list (update August 2022) in [DSSTox database](https://comptox.epa.gov/dashboard/chemical-lists/PFASSTRUCTV5), downloaded on 2 March 2026 as SDF v3000.

The inventory (`tests/test_data/PFASSTRUCTV5_20221101.sdf`, 14,735 structures) was used as a positive test set: every molecule is expected to satisfy the PFASSTRUCTv5 definition. The automated test `tests/test_pfasstructv5.py` checks this:

```bash
# Pytest mode — a deterministic sample of 10 molecules (fast, suitable for CI):
pytest tests/test_pfasstructv5.py -v

# Direct mode — validates all 14,735 molecules and reports any failures:
python tests/test_pfasstructv5.py
```

### HalogenGroups fingerprint comparison against CSRML (Richard *et al.* 2023)

HalogenGroups binary fingerprints were compared against the 125-bit CSRML fingerprint system introduced by Richard *et al.* 2023, which was developed specifically for profiling and categorising PFAS. The comparison script `benchmark/scripts/compare_pfasgroups_vs_txppfas.py` generates:

- Group coverage statistics and information-content distribution per group
- Lorenz curves for class-richness concentration
- PCA scatter plots comparing chemical-space coverage of both fingerprint systems
- Pairwise Jaccard similarity distributions within each system
- Cross-fingerprint similarity scatter (HalogenGroups vs. CSRML Jaccard)

> Richard, A. M. *et al.* (2023). *A New CSRML Structure-Based Fingerprint Method for Profiling and Categorizing Per- and Polyfluoroalkyl Substances (PFAS).* Chemical Research in Toxicology, 36(3), 318–338. [DOI: 10.1021/acs.chemrestox.2c00403](https://doi.org/10.1021/acs.chemrestox.2c00403)

```bash
cd benchmark/scripts

# HalogenGroups fingerprint analysis only (no CSRML CSV needed):
python compare_pfasgroups_vs_txppfas.py --smiles ../data/test_set_for_PFASSTRUCTv5.tsv

# Full comparison (requires Richard 2023 SI Table S2 exported as CSV):
python compare_pfasgroups_vs_txppfas.py \
  --smiles ../data/test_set_for_PFASSTRUCTv5.tsv \
  --txppfas_csv ../data/richard2023_SI_table_S2.csv
```


## Licence
<a rel="license" href="http://creativecommons.org/licenses/by-nd/4.0/"><img alt="Creative Commons License" style="border-width:0" src="https://i.creativecommons.org/l/by-nd/4.0/88x31.png" /></a><br />This work is licensed under a <a rel="license" href="http://creativecommons.org/licenses/by-nd/4.0/">Creative Commons Attribution-NoDerivatives 4.0 International License</a>.

Contact me in case you want an exception to the No Derivatives term.

## Acknowledgments
This project is part of the [ZeroPM project](https://zeropm.eu/) (WP2) and has received funding from the European Union’s Horizon 2020 research and innovation programme under grant agreement No 101036756. This work was developed at the [Department of Environmental Science](https://aces.su.se) at Stockholm University.<br />


<img alt="EU logo" src="https://zeropm.eu/wp-content/uploads/2021/12/flag_yellow_low.jpg" width=100/>     <a rel='zeropm_web' href="https://zeropm.eu/"/><img alt="zeropm logo" src="https://zeropm.eu/wp-content/uploads/2022/01/ZeroPM-logo.png" width=250 /></a><a rel='zeropm_web' href="https://su.se/"/><img alt="zeropm logo" src="https://eu01web.zoom.us/account/branding/p/5065401a-9915-4baa-9c16-665dcd743470.png" width=200 /></a>

[![Powered by RDKit](https://img.shields.io/badge/Powered%20by-RDKit-3838ff.svg?logo=data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAABAAAAAQBAMAAADt3eJSAAAABGdBTUEAALGPC/xhBQAAACBjSFJNAAB6JgAAgIQAAPoAAACA6AAAdTAAAOpgAAA6mAAAF3CculE8AAAAFVBMVEXc3NwUFP8UPP9kZP+MjP+0tP////9ZXZotAAAAAXRSTlMAQObYZgAAAAFiS0dEBmFmuH0AAAAHdElNRQfmAwsPGi+MyC9RAAAAQElEQVQI12NgQABGQUEBMENISUkRLKBsbGwEEhIyBgJFsICLC0iIUdnExcUZwnANQWfApKCK4doRBsKtQFgKAQC5Ww1JEHSEkAAAACV0RVh0ZGF0ZTpjcmVhdGUAMjAyMi0wMy0xMVQxNToyNjo0NyswMDowMDzr2J4AAAAldEVYdGRhdGU6bW9kaWZ5ADIwMjItMDMtMTFUMTU6MjY6NDcrMDA6MDBNtmAiAAAAAElFTkSuQmCC)](https://www.rdkit.org/)
    
