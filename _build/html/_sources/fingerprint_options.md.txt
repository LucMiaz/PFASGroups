# PFAS embedding options

## All metrics

Per-group metrics add $N_G$ columns (one per selected PFAS group; $N_G$ = 114 for `group_selection='all'` with the default fluorine-only embedding, 119 groups defined in total,
28 for `'oecd'`, etc.), aggregated with `aggregation='mean'` (default) or `'median'` when a group
has multiple matched components. Molecule-wide metrics add a single trailing column (prefixed `mol:`)
independent of group selection. The number of halogens does not affect column count.

Notation: $V$ = C-atoms of the matched component; $|V|$ = their count; $d(u,v)$ = shortest-path
distance; $\epsilon(v) = \max_u d(v,u)$ = eccentricity; $r_{ij}$ = effective resistance
(pseudo-inverse of Laplacian); $|M|$ = heavy-atom count of the molecule; $N_C$ = total matched
components across all groups.

| Scope | Option | Bits added | Type | Formula | Agg.† | Structural insight |
|---|---|:---:|---|---|:---:|---|
| per group | `binary` | $N_G$ | binary | $b_g = \mathbf{1}[\exists\,C_i \text{ matching } g]$ | — | Presence/absence indicator. Most information-dense metric per bit; default. Sufficient alone for OECD category assignment and group-level classification. |
| per group | `count` | $N_G$ | int ≥ 0 | $c_g = \lvert\{C_i : C_i \text{ matches } g\}\rvert$ | — | Quantifies structural repetition. Distinguishes monomeric from oligomeric/polymeric structures where the same group motif recurs. |
| per group | `max_component` | $N_G$ | int ≥ 0 | $x_g = \max_i \lvert C_i \rvert$ | — | Size (C-atoms) of the dominant matched fragment. Proxy for longest chain length; correlated with persistence and bioaccumulation potential. |
| per group | `total_component` | $N_G$ | int ≥ 0 | $x_g = \sum_i \lvert C_i \rvert$ | — | Total halogenated carbon burden for this group. Measures cumulative fluorination load; useful for hazard scoring and mass-fraction estimation. |
| per group | `size` | $N_G$ | int ≥ 0 | $s = \lvert V \rvert$ | yes | Raw fragment size independent of group identity. Redundant with `max_component` for single-component matches, but distinct when a group has multiple components. |
| per group | `n_spacer` | $N_G$ | int ≥ 0 | $n = \lvert\text{CH}_2\text{ linker atoms}\rvert$ | yes | Telomer CH₂ spacer length — the "m" in "m:n" fluorotelomer notation, e.g. 2 for 4:2 FTOH, 4 for 6:2 FTOH. Zero for all non-telomer groups. Encodes the linker chain length validated by the `linker_smarts` constraint; enables one feature to distinguish 2:1, 4:2, 6:2, and 8:2 telomers. |
| per group | `ring_size` | $N_G$ | int ≥ 0 | $r = \min\{\lvert R\rvert : R \in \mathcal{R}(G),\, R \cap C \neq \emptyset\}$ | yes | Size of the smallest ring overlapping the matched component. Zero for acyclic chains; 5 for five-membered heterocycles (azoles, furans, oxazolines); 6 for benzene rings and cyclohexane derivatives. Only non-zero for cyclic and aryl component forms. |
| per group | `branching` | $N_G$ | float [0, 1] | $\beta = d(G)\,/\,(\lvert V\rvert - 1)$; $\beta = 1$ for a linear path, $\to 0$ for star topology | yes | Topological linearity of the halogenated chain. Linear perfluoroalkyl chains score near 1; highly branched or cyclic structures score lower. Complementary to $\Omega$ in multi-metric presets. |
| per group | `mean_eccentricity` | $N_G$ | float ≥ 0 | $\bar\epsilon = \tfrac{1}{\lvert V\rvert}\sum_{v \in V}\epsilon(v)$ | yes | Average topological reach per atom. Large for long linear chains, small for compact rings. Correlated with chain-length-driven persistence. |
| per group | `median_eccentricity` | $N_G$ | float ≥ 0 | $\tilde\epsilon = \mathrm{median}_{v \in V}\,\epsilon(v)$ | yes | Robust version of mean eccentricity; less sensitive to a single terminal atom. Useful for structures with heterogeneous fragment lengths. |
| per group | `diameter` | $N_G$ | int ≥ 0 | $d = \max_{u,v \in V} d(u,v)$ | yes | Longest shortest path in the fragment. Directly tracks chain length for linear structures ($d = \lvert V\rvert - 1$); key predictor of passive diffusion and membrane permeability. |
| per group | `radius` | $N_G$ | int ≥ 0 | $r = \min_{v \in V}\epsilon(v)$ | yes | Minimum eccentricity; how "central" the most central atom is. Low for compact symmetric fragments, high for highly asymmetric chains. Complements diameter. |
| per group | `effective_graph_resistance` | $N_G$ | float ≥ 0 | $\Omega = \sum_{i < j \in V} r_{ij}$, uniform $w_{ij} = 1$, C-skeleton | yes | Kirchhoff index (resistance-distance sum). Captures chain topology globally: large for long linear chains, smaller for branched or cyclic structures. **Best single graph metric for group discrimination** (used in presets `best`–`best_3`). |
| per group | `effective_graph_resistance_BDE` | $N_G$ | float ≥ 0 | $\Omega_\mathrm{BDE} = \sum_{i < j \in V'} r_{ij}^{(\mathrm{BDE})}$, $w_{ij} = \mathrm{BDE}_{ij}/\mathrm{BDE}_{CC}$, 1-hop expanded $V'$ | yes | BDE-weighted Kirchhoff index on the expanded subgraph (C-skeleton + bonded F/H/Cl…). Encodes bond-strength distribution; distinguishes weak C–F bonds in highly fluorinated chains from stronger C–C backbones. |
| per group | `component_fraction` | $N_G$ | float [0, 1] | $f = \lvert V'\rvert / \lvert M\rvert$, $V'$ = expanded component | yes | Fraction of the molecule covered by this fragment. Approaches 1 for highly fluorinated molecules (e.g. PFOA); indicates whether the matched group dominates the molecular structure. |
| per group | `min_dist_to_centre` | $N_G$ | int ≥ 0 | $\min_{v \in S,\, c \in \mathrm{center}(G)} d(v,c)$ | yes | Topological proximity of the SMARTS-matched atoms to the molecular graph centre. Zero when the functional group is at the centre; large when it is a chain terminus. |
| per group | `max_dist_to_periphery` | $N_G$ | int ≥ 0 | $\max_{v \in S,\, p \in \mathrm{periphery}(G)} d(v,p)$ | yes | How far the functional group extends towards the molecular periphery. Encodes whether the matched group is internally embedded vs. terminally exposed. |
| per group | `min_dist_to_barycentre` | $N_G$ | int ≥ 0 | $\min_{v \in S,\, b \in \mathrm{bary}(G)} d(v,b)$, bary = node(s) minimising $\sum_u d(v,u)$ | yes | Distance of the functional group to the topological centre of mass. Complements `min_dist_to_centre` for asymmetric chains and multi-functional molecules. |
| molecule-wide | `n_components` | 1 | int ≥ 0 | $N_C = \lvert\{(g,i)\}\rvert$ across all groups | — | Total structural complexity indicator. Distinguishes simple mono-functional PFAS from complex oligomeric or multi-head structures. |
| molecule-wide | `total_size` | 1 | int ≥ 0 | $S = \sum_{g,i}\lvert C_{g,i}\rvert$ | — | Cumulative halogenated carbon count. Strongly correlated with molecular weight of the halogenated fraction and overall fluorine payload. |
| molecule-wide | `mean_size` | 1 | float ≥ 0 | $\bar s = S\,/\,N_C$ | — | Average fragment size. Distinguishes many short fragments from few long chains at equal total halogenation. |
| molecule-wide | `max_size` | 1 | int ≥ 0 | $s_\max = \max_{g,i}\lvert C_{g,i}\rvert$ | — | Size of the single largest fragment. Dominant predictor of the longest halogenated chain; primary driver of persistence and bioaccumulation. |
| molecule-wide | `mean_branching` | 1 | float [0, 1] | $\bar\beta = N_C^{-1}\sum_i \beta_i$ | — | Average chain linearity across all components. Useful for ranking molecules by their overall chain vs. branched character. |
| molecule-wide | `max_branching` | 1 | float [0, 1] | $\beta_\max = \max_i \beta_i$ | — | Linearity of the most linear component. Useful when a molecule contains a mix of chain types (e.g. one linear PFAS chain + one branched fragment). |
| molecule-wide | `mean_eccentricity` | 1 | float ≥ 0 | $N_C^{-1}\sum_i \bar\epsilon_i$ | — | Mean eccentricity averaged over all components; molecule-level chain-length proxy independent of which specific groups were matched. |
| molecule-wide | `max_diameter` | 1 | int ≥ 0 | $d_\max = \max_i d_i$ | — | Diameter of the longest fragment in the molecule. Directly tracks the longest halogenated chain; correlated with slow elimination rates. |
| molecule-wide | `mean_component_fraction` | 1 | float [0, 1] | $\bar f = N_C^{-1}\sum_i f_i$ | — | Average fraction of the molecule each fragment spans. Distinguishes molecules where halogenation is localised vs. spread over the whole backbone. |
| molecule-wide | `max_component_fraction` | 1 | float [0, 1] | $f_\max = \max_i f_i$ | — | Fraction of the molecule covered by the dominant fragment. Approaches 1 for highly fluorinated structures like PFOA; strong indicator of overall fluorination dominance. |

† **Agg.** applies only to per-group graph metrics when a group has **multiple matched components**.
`yes` = the per-component values are reduced to a single number via `aggregation='mean'` (default)
or `aggregation='median'`. Count-mode metrics (`binary`, `count`, `max_component`, `total_component`)
and all molecule-wide metrics are never aggregated this way.

## Notes

- See `PFASGroups/ComponentsSolverModel.py` and `PFASGroups/results_model.py` for exact implementations.
- Benchmarked presets are available via `EMBEDDING_PRESETS` (alias `FINGERPRINT_PRESETS`);
  see `ResultsFingerprint_Guide.md` for the full preset table with mean Tanimoto scores.
- The two `effective_graph_resistance` variants are always computed together since they share the
  Laplacian pseudo-inverse step:
  - `effective_graph_resistance` — topological (unweighted), original C-skeleton
  - `effective_graph_resistance_BDE` — bond-strength-weighted, 1-hop expanded component
    (adds all atoms directly bonded to any C-skeleton atom, including F, H, Cl, Br)


