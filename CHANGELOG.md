# Changelog

All notable changes to PFASGroups (and its `HalogenGroups` entry point) are
listed here, newest first. Versions follow [semantic versioning](https://semver.org/).

## Unreleased

**Documentation**

- Examples updated to the 3.4 API: custom groups are passed as `pfas_groups=` (`halogen_groups=` was silently ignored), `HalogenGroup` takes `id` and a `{SMARTS: count}` dict, definition matches are read from `mol.matches` (`m.is_definition`).
- PFAS definitions page rewritten for the six built-in definitions and the `ratioBasis` option; prioritization page rewritten to describe the actual scores (intrinsic fluorination, or cosine similarity to a reference).
- Group and column counts corrected (118 compiled groups, 117 per halogen block); version shown in the docs read from the package metadata.

## Version 3.4.6

**Released:** 2026-11-12

**Added**

- `ratioBasis` key for PFAS definitions (`ratio_basis` in Python): the fluorine ratio is computed by mass (`"mass"`, the default, as in 3.4.x) or by atom count (`"count"`). `applies_to_molecule` and `match_details` accept a `ratio_basis` argument for one call, and `match_details` and the parser results report `fluorine_ratio_basis`.

**Changed**

- PFASSTRUCTv5 uses `"ratioBasis": "count"`: at least 30% of the non-hydrogen atoms are fluorine, as Gaines et al. (2023, <doi:10.1002/ieam.4735>) define it ("strictly on atom counting without considering the weight of the atoms"). With the mass ratio of 3.4.x, molecules such as 4-(trifluoromethyl)benzoic acid (0.31 by mass, 0.23 by count) matched PFASSTRUCTv5 although they are outside the published definition.

## Version 3.4.5

**Released:** 2026-09-11

- Fixed difluoromethane in the PFAS definitions.

## Version 3.4.4

- PFAS definition matches report the fluorine ratio and whether the structural (SMARTS) criterion matched (`fluorine_ratio`, `fluorine_ratio_threshold`, `smarts_match`, `ratio_match`), for debugging or display.

## Version 3.4.3

- Fixed the fluorine-ratio computation of PFASSTRUCTv5.

## Version 3.4.2

- Added OPPT 2021 next to OPPT 2023 to remove the ambiguity between the proposed (2021) and final (2023) EPA definitions; OPPT 2023 leaves the R groups unconstrained, as in 40 CFR 705.3.

## Version 3.4.1

- Generic embedding for non-halogenated compounds, with the wildcard functional groups and the halogen component size.
- Fixed the H components and added cyclic and aryl H components.

## Version 3.4.0

- H components (`halogens='H'`), used to generate hydrocarbon homologue series, and wildcard groups (`halogens='*'`) matching generic functional groups beyond PFAS. See `advanced_features`.

## Version 3.3.x

- GUI prototype.

## Version 3.2.2

- Fixed polyhalogenated alkyl components matching fewer than 2 halogens.
- `halogens` can be used in formula constraints, and component-wide formula constraints apply to the distance-1 neighbours of C-only components.

## Version 3.2.1

- n-spacer metric for telomers and ring size for aryl and cyclic groups, usable in the embeddings.

## Version 3.2.0

**New features:**

- **\`\`get_compiled_HalogenGroups()\`\`**: Returns compiled `HalogenGroup` instances (`compute=True` groups only), suitable for extending with custom groups and passing directly to `parse_smiles` or `generate_fingerprint`. The older `get_HalogenGroups` continues to return raw JSON dicts.

  ``` python
  from HalogenGroups import get_compiled_HalogenGroups, HalogenGroup, parse_smiles

  groups = get_compiled_HalogenGroups()
  groups.append(HalogenGroup(
      id=200, name="my_custom_group",
      smarts={"[CX4](F)(F)(F)": 1},
  ))
  results = parse_smiles(["FC(F)(F)C(F)(F)F"], pfas_groups=groups)
  ```

**API changes:**

- `generate_fingerprint` return value is now a 2-D `numpy.ndarray` of shape `(n_molecules, n_groups)` instead of a Python list.

## Version 3.1.4

- Fingerprint parameters split into molecule-wide and component-wide, with presets for the best combinations.

## Version 3.1.3

**Released:** 2026

- Documentation rewrite for ReadTheDocs publication
- Both `HalogenGroups` and `PFASGroups` entry points documented
- All code examples updated to current API (`ResultsModel` pattern)
- CLI reference page added
- API reference pages updated with autodoc

## Version 3.1.0

**New features:**

- **Multi-halogen fingerprinting**: `generate_fingerprint` and `to_fingerprint` now accept a `halogens` parameter. Multiple halogens produce horizontally-stacked vectors (116 × n_halogens columns).
- **Saturation filter**: `saturation` parameter (`'saturated'` \| `'unsaturated'` \| `None`) for both parsing and fingerprinting.
- **Corrected group-selection indexing**: selections resolve by canonical group IDs rather than list positions.
- **\`\`ResultsFingerprint\`\` metadata**: `halogens` and `saturation` stored in the object and shown in `__repr__` / `summary()`.

## Version 2.2.4

**New features:**

- **ResultsFingerprint class** with group selection, multiple encoding modes (binary / count / max_component), PCA, t-SNE, UMAP, and KL divergence.
- **Database persistence**: `to_sql()` / `from_sql()` for ResultsFingerprint and ResultsModel (SQLite and PostgreSQL via SQLAlchemy).
- 100+ new unit tests.

## Version 2.2.3

**New features:**

- **ResultsModel container**: wraps `parse_smiles` / `parse_mols` output with helpers `show()`, `summarise()`, `table()`, and `to_dataframe()`.
- Extended visualisation utilities (`plot_pfasgroups`).
- Updated PFAS group definitions to 116 groups.

## Version 2.2

**Released:** February 2026

- `linked_smarts` option to restrict the path between the SMARTS groups and the fluorinated component.
- New telomer groups.
- **2.2.2:** fixed the telomer groups, added examples and counter-examples to each group, and removed the boundary O from fluorinated components (per- and polyfluoroalkyl).

## Version 2.1

**Released:** January 2026

- Several SMARTS per group, each with its own minimum count.

## Version 2.0

- Complete rewrite of the SMARTS-based matcher.
- New JSON-based group definition format.
- Regulatory PFAS definition support (OECD, REACH, OPPT, UK EA, PFASSTRUCTv5).
- CLI entry points `halogengroups` / `pfasgroups` added.

## Version 1.x — component-based analysis

- Chain finding replaced by connected-component analysis.
- Branched and cyclic structures supported.
- Improved SMARTS matching for diverse PFAS classes.

## Version 0.x — path-based analysis

- A SMARTS match connected to a second SMARTS, or to a default path SMARTS, found with networkx `shortest_path`.
