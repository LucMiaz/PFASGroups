# PFAS embedding options

This table lists all valid `component_metrics` and `molecule_metrics` options for
`PFASEmbedding.to_array()` / `PFASEmbeddingSet.to_array()`, their level, description,
definition, and a rough importance ranking.

## Component metrics

Each component metric produces one block of columns — one per PFAS group (`n_groups` ≈ 116).
When a group has multiple matched components the graph metric is aggregated with
`aggregation='mean'` (default) or `'median'`.

| Option | Level | Description | Formula / definition | Importance |
|---|---:|---|---|---|
| `binary` | Count mode per PFAS group | Presence indicator for a PFAS group in the molecule | 1 if any component matching the group exists, else 0 | High |
| `count` | Count mode per PFAS group | Number of matched components for the PFAS group | Count of components matched to the group | Medium |
| `max_component` | Count mode per PFAS group | Maximum component size (C-atom count) among matched components | max(comp.size) | Medium |
| `total_component` | Count mode per PFAS group | Sum of component sizes for the group | sum(comp.size) | Medium |
| `size` | Graph metric (same across groups) | Number of C-atoms in the component | integer atom count | Medium |
| `branching` | Graph metric (same across groups) | Branching metric (1.0 = linear → 0 = highly branched) | computed graph branching metric per component | Medium |
| `mean_eccentricity` | Graph metric (same across groups) | Average eccentricity of nodes in component's graph | mean(eccentricity(node)) | Medium |
| `median_eccentricity` | Graph metric (same across groups) | Median eccentricity across component nodes | median(eccentricity) | Low–Medium |
| `diameter` | Graph metric (same across groups) | Maximum shortest-path distance within component graph | max shortest-path length | Medium |
| `radius` | Graph metric (same across groups) | Minimum eccentricity of nodes (graph radius) | min eccentricity | Low |
| `effective_graph_resistance` | Graph metric (same across groups) | Kirchhoff index — **uniform** edge weights (=1), **original** C-skeleton component | sum of pairwise resistance distances via Laplacian pseudoinverse, all conductances = 1 | High |
| `effective_graph_resistance_BDE` | Graph metric (same across groups) | Kirchhoff index — **BDE-calibrated** edge weights, component **expanded 1 hop** (includes F/H/Cl/Br neighbours) | same formula but conductance = BDE(z1,z2,order)/BDE(C,C,1) on the expanded subgraph | High |
| `component_fraction` | Graph metric (same across groups) | Fraction of molecule C-atoms covered by this component | C-atoms in component / total C-atoms in molecule | High |
| `min_dist_to_center` | Graph metric per PFAS group | Minimum topological distance from component atoms to molecular graph center | min(shortest_path(atom, center)) | Medium |
| `max_dist_to_periphery` | Graph metric per PFAS group | Maximum topological distance from component atoms to periphery nodes | max(shortest_path(atom, periphery)) | Medium |
| `min_dist_to_barycenter` | Graph metric per PFAS group | Minimum topological distance from component to barycenter | min(shortest_path(atom, barycenter)) | Medium |

## Molecule metrics

Molecule-wide scalars appended as trailing columns (prefixed `mol:` in column names).
These are independent of group selection — they summarise all matched components across
the whole molecule.

| Option | Description | Formula / definition | Importance |
|---|---|---|---|
| `n_components` | Total matched components across the molecule | count of components | High |
| `total_size` | Sum of component sizes (C-atoms) across molecule | sum(comp.size) | High |
| `mean_size` | Mean component size | mean(comp.size) | Medium |
| `max_size` | Maximum component size | max(comp.size) | Medium |
| `mean_branching` | Mean branching across components | mean(comp.branching) | Medium |
| `max_branching` | Maximum branching value | max(comp.branching) | Low |
| `mean_eccentricity` | Mean eccentricity aggregated across components | mean(component mean_eccentricity) | Medium |
| `max_diameter` | Largest component diameter | max(component.diameter) | Medium |
| `mean_component_fraction` | Average fraction coverage per component | mean(component_fraction) | High |
| `max_component_fraction` | Maximum component fraction | max(component_fraction) | Medium |

## Notes

- **Level** indicates whether the metric contributes columns per PFAS group (graph metrics
  per group), molecule-level scalar columns appended once, or component-level values
  intrinsic to components regardless of grouping.
- **Formula** entries are concise descriptions; see `PFASGroups/ComponentsSolverModel.py`
  and `PFASGroups/results_model.py` for exact implementations.
- **Importance** is expert judgement; the `EMBEDDING_PRESETS` dict (alias
  `FINGERPRINT_PRESETS`) gives benchmarked presets (e.g. `'best'`) optimised for
  inter-group discrimination.
- The two `effective_graph_resistance` variants are always computed together:
  - `effective_graph_resistance` — topological (unweighted), original C-skeleton
  - `effective_graph_resistance_BDE` — bond-strength-weighted, 1-hop expanded component
    (adds all atoms directly bonded to any C-skeleton atom, including F, H, Cl, Br)
