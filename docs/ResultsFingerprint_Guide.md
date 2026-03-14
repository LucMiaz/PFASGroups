:orphan:

# PFASEmbeddingSet: Embeddings and Analysis

## Overview

`PFASEmbeddingSet` is the container returned by `parse_smiles()`. Each element is a
`PFASEmbedding` — a dict-like object holding the parsed PFAS group results for one molecule.
Call `.to_array()` on either to generate numeric embedding vectors from the stored data.

- **`PFASEmbedding`** — single-molecule result + embedding generator (`to_array()`)
- **`PFASEmbeddingSet`** — list of `PFASEmbedding` objects; batch `to_array()` → `(n_mols, n_cols)` matrix

Backward-compatible aliases: `MoleculeResult = PFASEmbedding`, `ResultsModel = PFASEmbeddingSet`.

## Quick Start

```python
from PFASGroups import parse_smiles

smiles_list = [
    "C(C(C(C(C(C(C(C(F)(F)F)(F)F)(F)F)(F)F)(F)F)(F)F)(F)F)(C(=O)O)F",  # PFOA
    "FC(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)S(=O)(=O)O",  # PFOS
]

results = parse_smiles(smiles_list)   # → PFASEmbeddingSet

# Binary fingerprint (default): shape (2, 116)
arr = results.to_array()

# Best-performing preset: binary + effective_graph_resistance
arr, cols = results.to_array(preset='best'), results.column_names(preset='best')

# Custom metrics with median aggregation
arr = results.to_array(
    component_metrics=['binary', 'effective_graph_resistance', 'effective_graph_resistance_BDE'],
    aggregation='median',
)

# Single-molecule embedding
emb = results[0]                       # PFASEmbedding (dict subclass)
vec = emb.to_array(preset='best')      # 1-D numpy array

# Column names without computing values
cols = results.column_names(preset='best')
```

## API Reference

### `PFASEmbeddingSet.to_array()`

Stacks per-molecule rows into a `(n_mols, n_cols)` matrix.

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `component_metrics` | `list[str]` | `['binary']` | Per-component metrics — see table below |
| `molecule_metrics` | `list[str]` | `None` | Molecule-wide scalars appended as trailing columns |
| `group_selection` | `str` | `'all'` | Which groups to include (see table below) |
| `selected_group_ids` | `list[int]` | `None` | Explicit group IDs — overrides `group_selection` |
| `aggregation` | `str` | `'mean'` | Aggregate multiple components per group: `'mean'` or `'median'` |
| `preset` | `str` | `None` | Named configuration — overrides `component_metrics` / `molecule_metrics` |
| `pfas_groups` | `list` | `None` | Custom group list (loaded from defaults when `None`) |

**Returns:** `numpy.ndarray` of shape `(n_mols, n_cols)`.

### `PFASEmbedding.to_array()`

Same parameters as `PFASEmbeddingSet.to_array()`. **Returns:** 1-D `numpy.ndarray` of length `n_cols`.

### `PFASEmbeddingSet.column_names()` / `PFASEmbedding.column_names()`

Returns the list of column labels without computing values. Parameters match `to_array()`
(the `aggregation` parameter is not relevant for column names).

---

### Group selection

| Value | Groups included | Count |
|-------|----------------|-------|
| `'all'` | All computable groups | 116 |
| `'oecd'` | OECD 2021 list | 28 |
| `'generic'` | Generic/structural groups | 27 |
| `'telomers'` | Telomer-type groups | 43 |
| `'generic+telomers'` | Combined | 70 |

---

### Component metrics

| Metric | Type | Description |
|--------|------|-------------|
| `'binary'` | count mode | 1 if the group is present, 0 otherwise |
| `'count'` | count mode | Number of matched components for this group |
| `'max_component'` | count mode | Maximum component size (C-atom count) |
| `'total_component'` | count mode | Sum of component sizes |
| `'effective_graph_resistance'` | graph metric | Kirchhoff index — uniform edge weights, original C-skeleton component |
| `'effective_graph_resistance_BDE'` | graph metric | Kirchhoff index — BDE-calibrated weights, component expanded 1 hop (includes F/H/Cl…) |
| `'branching'` | graph metric | Branching metric (1.0 = linear, 0.0 = highly branched) |
| `'mean_eccentricity'` | graph metric | Mean node eccentricity |
| `'median_eccentricity'` | graph metric | Median node eccentricity |
| `'diameter'` | graph metric | Maximum shortest-path length |
| `'radius'` | graph metric | Minimum eccentricity |
| `'component_fraction'` | graph metric | Fraction of molecule C-atoms in this component |
| `'min_dist_to_barycenter'` | graph metric | Min topological distance from component to barycenter |
| `'min_dist_to_center'` | graph metric | Min topological distance from component to graph center |
| `'max_dist_to_periphery'` | graph metric | Max topological distance from component to periphery |
| `'size'` | graph metric | Number of C-atoms in the component |

When a group has multiple matched components the graph metric is aggregated across them
using the `aggregation` parameter (`'mean'` by default, or `'median'`).

### Molecule metrics

Appended as trailing columns (prefixed `mol:` in column names).

| Metric | Description |
|--------|-------------|
| `'n_components'` | Total number of matched components |
| `'total_size'` | Sum of component sizes (C atoms) |
| `'mean_size'` | Mean component size |
| `'max_size'` | Maximum component size |
| `'mean_branching'` | Mean branching across components |
| `'max_branching'` | Maximum branching value |
| `'mean_eccentricity'` | Mean eccentricity aggregated across components |
| `'max_diameter'` | Largest component diameter |
| `'mean_component_fraction'` | Average component fraction |
| `'max_component_fraction'` | Maximum component fraction |

---

### Presets (`EMBEDDING_PRESETS`)

Named configurations benchmarked on inter-group Tanimoto discrimination.

| Preset | Metrics | Description |
|--------|---------|-------------|
| `'best'` | `['binary', 'effective_graph_resistance']` | Rank 1 — best discrimination |
| `'best_2'` | `['binary', 'branching', 'effective_graph_resistance']` | Rank 2 |
| `'best_3'` | `['binary', 'branching', 'mean_eccentricity', 'effective_graph_resistance']` | Rank 3 |
| `'binary'` | `['binary']` | Plain binary fingerprint |
| `'count'` | `['count']` | Count fingerprint |
| `'max_component'` | `['max_component']` | Max-component size fingerprint |

Access via `from PFASGroups import EMBEDDING_PRESETS` (alias: `FINGERPRINT_PRESETS`).

---

## `generate_embedding()` convenience function

Parses SMILES and returns `(array, column_names)` in one call.

```python
from PFASGroups import generate_embedding

arr, cols = generate_embedding(
    ["FC(F)(F)C(=O)O", "OCCS"],
    preset='best',
)
# arr.shape == (2, 232)   # 116 groups × 2 metrics
# cols[:3] == ['Perfluoroalkyl [binary]', 'Perfluoroalkyl [effective_graph_resistance]', ...]
```

| Parameter | Default | Description |
|-----------|---------|-------------|
| `smiles` | — | `str` or `list[str]` |
| `component_metrics` | `['binary']` | Per-component metrics |
| `molecule_metrics` | `None` | Molecule-wide scalars |
| `group_selection` | `'all'` | Group subset |
| `selected_group_ids` | `None` | Explicit group IDs |
| `aggregation` | `'mean'` | `'mean'` or `'median'` |
| `preset` | `None` | Named preset |
| `halogens` | `'F'` | Halogen(s) |
| `saturation` | `'per'` | Saturation filter |
| `progress` | `False` | tqdm progress bar |

**Returns:** `(numpy.ndarray, list[str])` — embedding array + column names.

---

## Column naming

Column names follow the pattern `"{group_name} [{metric}]"` for component metrics
and `"mol:{metric}"` for molecule metrics. Examples:

```
Perfluoroalkyl [binary]
Perfluoroalkyl [effective_graph_resistance]
Perfluoroalkyl [effective_graph_resistance_BDE]
mol:n_components
mol:mean_branching
```
