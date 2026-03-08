from rdkit import Chem
from rdkit.Chem.rdMolDescriptors import CalcMolFormula
from .parser import load_HalogenGroups, parse_groups_in_mol
from typing import Union, List, Optional, Dict, Any
import numpy as np


# ---------------------------------------------------------------------------
# Supported per-group graph-metric keys (present in matched_components dicts)
# ---------------------------------------------------------------------------
_COMPONENT_GRAPH_METRICS = {
    'branching', 'mean_eccentricity', 'median_eccentricity', 'diameter',
    'radius', 'effective_graph_resistance', 'component_fraction',
    'min_dist_to_center', 'max_dist_to_periphery',
    'min_dist_to_barycenter', 'min_resistance_dist_to_barycenter',
    'min_resistance_dist_to_center', 'max_resistance_dist_to_periphery',
    'size',
}

# Aggregation function used to collapse multiple components → one scalar per group
_GRAPH_METRIC_AGG: Dict[str, str] = {}   # empty = default 'mean' for all

# Supported molecule-wide metric names and how they are computed
# (from the pool of *all* matched_components across all groups for a molecule)
_MOL_METRIC_KEYS = {
    'n_components', 'total_size', 'mean_size', 'max_size',
    'mean_branching', 'max_branching',
    'mean_eccentricity', 'max_diameter',
    'mean_component_fraction', 'max_component_fraction',
}

# ---------------------------------------------------------------------------
# Benchmark-validated fingerprint presets (PFASGroups Tanimoto benchmark 2026)
# Sorted by mean inter-group Tanimoto (lower = more discriminating).
# ---------------------------------------------------------------------------
FINGERPRINT_PRESETS: Dict[str, Any] = {
    # ── Top-5 PFASGroups configs from Phase 1 (all-group pairwise benchmark) ─
    'best': {
        'count_mode': 'binary',
        'graph_metrics': ['effective_graph_resistance'],
        'molecule_metrics': None,
        'description': (
            'Rank 1 — binary + effective_graph_resistance.\n'
            'Mean inter-group Tanimoto = 0.184 (best discrimination, '
            'outperforms TxP-PFAS 129-bit fingerprint).'
        ),
    },
    'best_2': {
        'count_mode': 'binary',
        'graph_metrics': ['branching', 'effective_graph_resistance'],
        'molecule_metrics': None,
        'description': (
            'Rank 2 — binary + branching + effective_graph_resistance.\n'
            'Mean inter-group Tanimoto = 0.189.'
        ),
    },
    'best_3': {
        'count_mode': 'binary',
        'graph_metrics': ['branching', 'mean_eccentricity', 'effective_graph_resistance'],
        'molecule_metrics': None,
        'description': (
            'Rank 3 — binary + branching + mean_eccentricity + effective_graph_resistance.\n'
            'Mean inter-group Tanimoto = 0.196.'
        ),
    },
    'best_4': {
        'count_mode': 'total_component',
        'graph_metrics': ['branching', 'mean_eccentricity', 'effective_graph_resistance'],
        'molecule_metrics': None,
        'description': (
            'Rank 4 — total_component + branching + mean_eccentricity + effective_graph_resistance.\n'
            'Mean inter-group Tanimoto = 0.206.'
        ),
    },
    'best_5': {
        'count_mode': 'binary',
        'graph_metrics': ['branching', 'mean_eccentricity', 'diameter', 'radius',
                          'effective_graph_resistance'],
        'molecule_metrics': None,
        'description': (
            'Rank 5 — binary + 5 graph metrics.\n'
            'Mean inter-group Tanimoto = 0.207.'
        ),
    },
    # ── Convenience presets ────────────────────────────────────────────────
    'binary': {
        'count_mode': 'binary',
        'graph_metrics': [],
        'molecule_metrics': None,
        'description': 'Plain binary fingerprint — 1 if group present, 0 absent. No graph metrics.',
    },
    'count': {
        'count_mode': 'count',
        'graph_metrics': [],
        'molecule_metrics': None,
        'description': 'Count fingerprint — number of matched components per group.',
    },
    'max_component': {
        'count_mode': 'max_component',
        'graph_metrics': [],
        'molecule_metrics': None,
        'description': 'Max-component fingerprint — maximum component size (C-atom count) per group.',
    },
}


@load_HalogenGroups()
def generate_fingerprint(smiles: Union[str, List[str]],
                             selected_groups: Union[List[int], range, None] = None,
                             representation: str = 'vector',
                             count_mode: str = 'binary',
                             halogens: Union[str, List[str]] = 'F',
                             saturation: Optional[str] = 'per',
                             graph_metrics: Optional[List[str]] = None,
                             molecule_metrics: Optional[List[str]] = None,
                             preset: Optional[str] = None,
                             **kwargs):
    """
    Generate PFAS group fingerprints from SMILES strings.

    Parameters:
    -----------
    smiles : str or list of str
        SMILES string(s) to generate fingerprints for
    selected_groups : list of int, range, or None
        Indices of PFAS groups to include in fingerprint.
        If None, all groups are used. Examples: [28, 29, 30], range(28, 52)
    representation : str, default 'vector'
        Type of fingerprint representation:
        - 'vector': Binary/count vector (numpy array)
        - 'dict': Dictionary mapping group names to counts
        - 'sparse': Dictionary with only non-zero group counts
        - 'detailed': Full match information including chain details
        - 'int': Convert binary to int using int(x,2)
    count_mode : str, default 'binary'
        How to count matches:
        - 'binary': 1 if group present, 0 if absent
        - 'count': Number of matched components found
        - 'max_component': Maximum component size (C-atom count) among all
          components matching the group
        - 'total_component': Sum of all component sizes matching the group
        Ignored when ``preset`` is given.
    halogens : str or list of str, default 'F'
        Which halogen(s) to use for component SMARTS matching.
        - Single halogen (e.g. 'F'): standard fingerprint of length n_groups
        - List of halogens (e.g. ['F', 'Cl', 'Br', 'I']): fingerprints are generated per
          halogen and stacked, yielding length n_groups * n_halogens.
          Group names are suffixed with '[F]', '[Cl]', etc.
        Available: 'F', 'Cl', 'Br', 'I'
    saturation : str or None, default 'per'
        Saturation filter for component SMARTS matching:
        - 'per': perfluorinated / perhalogenated components only
        - 'poly': polyfluorinated / polyhalogenated components only
        - None: no saturation filter (use all)
        Only applies to groups that use component SMARTS (e.g. OECD groups).
        Groups without component SMARTS (generic/telomer) are unaffected.
    graph_metrics : list of str, optional
        Per-group component-level graph metrics to append as additional columns.
        For each metric, *n_groups* extra columns are added (mean aggregated over
        all matched components for that group), yielding a vector of length
        n_groups * (1 + len(graph_metrics)).  Column names are suffixed with
        the metric name, e.g. ``"Perfluoroalkyl [branching]"``.
        Supported: 'branching', 'mean_eccentricity', 'diameter', 'radius',
        'component_fraction', 'min_dist_to_center', 'max_dist_to_periphery',
        'min_dist_to_barycenter', 'effective_graph_resistance', etc.
    molecule_metrics : list of str, optional
        Molecule-wide scalar metrics appended as extra columns *after* the
        per-group columns.  These are computed from the pool of ALL matched
        components across all PFAS groups in the molecule.
        Supported values: 'n_components', 'total_size', 'mean_size',
        'max_size', 'mean_branching', 'max_branching', 'mean_eccentricity',
        'max_diameter', 'mean_component_fraction', 'max_component_fraction'.
    preset : str, optional
        Named benchmark-validated configuration.  When given, overrides
        ``count_mode``, ``graph_metrics``, and ``molecule_metrics`` with the
        preset values.  Available presets (see ``FINGERPRINT_PRESETS``):
        - ``'best'``: binary + effective_graph_resistance (rank 1, mean T=0.184)
        - ``'best_2'``: binary + branching + EGR (rank 2, mean T=0.189)
        - ``'best_3'``: binary + branching + mean_eccentricity + EGR (rank 3)
        - ``'best_4'``: total_component + branching + mean_eccentricity + EGR (rank 4)
        - ``'best_5'``: binary + 5 graph metrics (rank 5)
        - ``'binary'``: plain binary, no graph metrics
        - ``'count'``: component count per group
        - ``'max_component'``: max component size per group
    pfas_groups : list, optional
        Custom list of PFAS groups. If None, uses default groups.

    Returns:
    --------
    fingerprints : numpy.ndarray, dict, or list
        Fingerprint representation(s) based on 'representation' parameter.
        For single SMILES input, returns single fingerprint.
        For list input, returns list of fingerprints.
        When multiple halogens are given and representation='vector', the
        vectors are stacked horizontally (shape: n_molecules x (n_groups*n_hal)).
        When graph_metrics are requested the vector is extended to
        n_groups * (1 + len(graph_metrics)) [+ len(molecule_metrics)].
    group_info : dict
        Information aboutthe groups used in fingerprinting:
        - 'group_names': List of column names in fingerprint order
        - 'group_ids': List of group IDs (None for extra metric columns)
        - 'selected_indices': Indices of selected groups
        - 'halogens': Halogen(s) used
        - 'saturation': Saturation filter used
        - 'graph_metrics': Per-group metric names (may be empty list)
        - 'molecule_metrics': Molecule-level metric names (may be empty list)

    Examples:
    ---------
    >>> # Binary vector for all groups, fluorine only (default)
    >>> fp, info = generate_fingerprint('CC(F)(F)C(=O)O')

    >>> # Stacked fingerprint for F and Cl
    >>> fp, info = generate_fingerprint('CC(F)(F)C(=O)O', halogens=['F', 'Cl'])

    >>> # With per-group branching appended (vector length = 2 * n_groups)
    >>> fp, info = generate_fingerprint('CC(F)(F)C(=O)O',
    ...                                  graph_metrics=['branching'])

    >>> # With molecule-wide metrics appended
    >>> fp, info = generate_fingerprint(['CCF', 'CCFF'],
    ...                                  molecule_metrics=['n_components', 'mean_branching'])

    >>> # Best benchmark-validated fingerprint (binary + effective_graph_resistance)
    >>> fp, info = generate_fingerprint('CC(F)(F)C(=O)O', preset='best')

    >>> # Count-based dictionary representation
    >>> fp, info = generate_fingerprint(['CCF', 'CCFF'],
    ...                                  representation='dict',
    ...                                  count_mode='count')
    """
    # Apply preset if requested (overrides count_mode / graph_metrics / molecule_metrics)
    if preset is not None:
        if preset not in FINGERPRINT_PRESETS:
            raise ValueError(
                f"Unknown preset: {preset!r}. "
                f"Available: {sorted(FINGERPRINT_PRESETS)}"
            )
        _p = FINGERPRINT_PRESETS[preset]
        count_mode = _p.get('count_mode', count_mode)
        if _p.get('graph_metrics') is not None:
            graph_metrics = _p['graph_metrics']
        if _p.get('molecule_metrics') is not None:
            molecule_metrics = _p['molecule_metrics']

    pfas_groups = kwargs.get('pfas_groups')
    graph_metrics = list(graph_metrics) if graph_metrics else []
    molecule_metrics = list(molecule_metrics) if molecule_metrics else []

    # Validate graph_metrics
    unsupported_gm = [m for m in graph_metrics if m not in _COMPONENT_GRAPH_METRICS]
    if unsupported_gm:
        raise ValueError(
            f"Unsupported graph_metrics: {unsupported_gm}. "
            f"Supported: {sorted(_COMPONENT_GRAPH_METRICS)}"
        )
    unsupported_mm = [m for m in molecule_metrics if m not in _MOL_METRIC_KEYS]
    if unsupported_mm:
        raise ValueError(
            f"Unsupported molecule_metrics: {unsupported_mm}. "
            f"Supported: {sorted(_MOL_METRIC_KEYS)}"
        )

    # Normalise halogens to a list
    if isinstance(halogens, str):
        halogens_list = [halogens]
    else:
        halogens_list = list(halogens)

    # Handle single SMILES input
    if isinstance(smiles, str):
        smiles_list = [smiles]
        single_input = True
    else:
        smiles_list = smiles
        single_input = False

    # Determine which groups to use
    if selected_groups is None:
        selected_indices = list(range(len(pfas_groups)))
        selected_pfas_groups = pfas_groups
    else:
        if isinstance(selected_groups, range):
            selected_indices = list(selected_groups)
        else:
            selected_indices = selected_groups

        max_index = len(pfas_groups) - 1
        invalid_indices = [i for i in selected_indices if i < 0 or i > max_index]
        if invalid_indices:
            raise ValueError(f"Invalid group indices {invalid_indices}. "
                           f"Valid range is 0-{max_index}")

        selected_pfas_groups = [pfas_groups[i] for i in selected_indices]

    # Build base group names/ids (with halogen suffix when multi-halogen)
    multi_halogen = len(halogens_list) > 1
    if multi_halogen:
        base_group_names = [f"{g.name} [{h}]" for h in halogens_list for g in selected_pfas_groups]
        base_group_ids: List[Optional[int]] = [g.id for h in halogens_list for g in selected_pfas_groups]
    else:
        base_group_names = [g.name for g in selected_pfas_groups]
        base_group_ids = [g.id for g in selected_pfas_groups]

    # Extend names/ids with per-group graph metric columns
    extended_group_names = list(base_group_names)
    extended_group_ids: List[Optional[int]] = list(base_group_ids)
    for gm in graph_metrics:
        for name, gid in zip(base_group_names, base_group_ids):
            extended_group_names.append(f"{name} [{gm}]")
            extended_group_ids.append(gid)

    # Append molecule-metric column names (one scalar each)
    mol_metric_names = [f"mol:{m}" for m in molecule_metrics]
    all_column_names = extended_group_names + mol_metric_names
    all_column_ids = extended_group_ids + [None] * len(molecule_metrics)

    group_info: Dict[str, Any] = {
        'group_names': all_column_names,
        'group_ids': all_column_ids,
        'selected_indices': selected_indices,
        'total_groups': len(pfas_groups),
        'halogens': halogens_list,
        'saturation': saturation,
        'graph_metrics': graph_metrics,
        'molecule_metrics': molecule_metrics,
    }

    # Base kwargs without injected keys
    kwargs_base = {k: v for k, v in kwargs.items() if k != 'pfas_groups'}

    def _parse_for_halogen(mol, formula, halogen):
        """Run parse_groups_in_mol for one halogen."""
        kw = dict(kwargs_base)
        kw['halogens'] = halogen
        if saturation is not None:
            kw['saturation'] = saturation
        all_matches, _ = parse_groups_in_mol(
            mol, formula=formula, pfas_groups=pfas_groups,
            include_PFAS_definitions=False, **kw
        )
        match_dict = {}
        for group, match_count, component_sizes, matched_components in all_matches:
            match_dict[group.id] = {
                'group': group,
                'match_count': match_count,
                'component_sizes': component_sizes,
                'matched_components': matched_components,
            }
        return match_dict

    def _agg_comp_metric(comps, metric):
        """Return mean of *metric* over a list of component dicts (0.0 if empty)."""
        values = [comp.get(metric) for comp in comps if comp.get(metric) is not None]
        return float(np.mean(values)) if values else 0.0

    def _build_graph_metric_vec(md, groups, metric):
        """Per-group mean of *metric* across matched components (float vector)."""
        vec = np.zeros(len(groups), dtype=float)
        for i, group in enumerate(groups):
            if group.id in md:
                comps = md[group.id].get('matched_components', [])
                vec[i] = _agg_comp_metric(comps, metric)
        return vec

    def _build_mol_metrics(all_match_dicts):
        """Compute molecule-wide scalar metrics from the union of all match dicts."""
        # Gather ALL matched_components across all halogens and all groups
        all_comps = []
        for md in all_match_dicts.values():
            for mi in md.values():
                all_comps.extend(mi.get('matched_components', []))
        result = []
        for m in molecule_metrics:
            if m == 'n_components':
                result.append(float(len(all_comps)))
            elif m == 'total_size':
                result.append(float(sum(c.get('size', 0) or 0 for c in all_comps)))
            elif m == 'mean_size':
                sizes = [c.get('size', 0) or 0 for c in all_comps]
                result.append(float(np.mean(sizes)) if sizes else 0.0)
            elif m == 'max_size':
                sizes = [c.get('size', 0) or 0 for c in all_comps]
                result.append(float(max(sizes)) if sizes else 0.0)
            elif m == 'mean_branching':
                result.append(_agg_comp_metric(all_comps, 'branching'))
            elif m == 'max_branching':
                vals = [c.get('branching') for c in all_comps if c.get('branching') is not None]
                result.append(float(max(vals)) if vals else 0.0)
            elif m == 'mean_eccentricity':
                result.append(_agg_comp_metric(all_comps, 'mean_eccentricity'))
            elif m == 'max_diameter':
                vals = [c.get('diameter') for c in all_comps if c.get('diameter') is not None]
                result.append(float(max(vals)) if vals else 0.0)
            elif m == 'mean_component_fraction':
                result.append(_agg_comp_metric(all_comps, 'component_fraction'))
            elif m == 'max_component_fraction':
                vals = [c.get('component_fraction') for c in all_comps
                        if c.get('component_fraction') is not None]
                result.append(float(max(vals)) if vals else 0.0)
            else:
                result.append(0.0)
        return np.array(result, dtype=float)

    fingerprints = []

    for smiles_str in smiles_list:
        try:
            mol = Chem.MolFromSmiles(smiles_str)
            if mol is None:
                raise ValueError(f"Invalid SMILES: {smiles_str}")

            formula = CalcMolFormula(mol)

            # Collect match dicts per halogen
            match_dicts = {h: _parse_for_halogen(mol, formula, h) for h in halogens_list}
            # First halogen's match_dict used for non-stacked repr (dict/sparse/detailed)
            match_dict = match_dicts[halogens_list[0]]

            def _build_vector(md, groups):
                """Build a count/binary vector for one match_dict + group list."""
                vec = np.zeros(len(groups), dtype=float if graph_metrics else int)
                for i, group in enumerate(groups):
                    if group.id in md:
                        match_info = md[group.id]
                        if count_mode == 'binary' or representation == 'int':
                            vec[i] = 1
                        elif count_mode == 'count':
                            vec[i] = match_info['match_count']
                        elif count_mode == 'max_component':
                            if match_info['matched_components']:
                                vec[i] = max(comp.get('size', 0) or 0
                                             for comp in match_info['matched_components'])
                            else:
                                vec[i] = match_info['match_count'] if match_info['match_count'] > 0 else 0
                        elif count_mode == 'total_component':
                            if match_info['matched_components']:
                                vec[i] = sum(comp.get('size', 0) or 0
                                             for comp in match_info['matched_components'])
                            else:
                                vec[i] = match_info['match_count'] if match_info['match_count'] > 0 else 0
                        else:
                            raise ValueError(f"Unknown count_mode: {count_mode}")
                return vec

            if representation == 'vector' or representation == 'int':
                if multi_halogen:
                    parts = [_build_vector(match_dicts[h], selected_pfas_groups)
                             for h in halogens_list]
                    base_vec = np.concatenate(parts)
                else:
                    base_vec = _build_vector(match_dict, selected_pfas_groups)

                # Append per-group graph metric columns
                extra_parts = []
                for gm in graph_metrics:
                    if multi_halogen:
                        gm_parts = [_build_graph_metric_vec(match_dicts[h], selected_pfas_groups, gm)
                                    for h in halogens_list]
                        extra_parts.append(np.concatenate(gm_parts))
                    else:
                        extra_parts.append(_build_graph_metric_vec(match_dict, selected_pfas_groups, gm))

                # Append molecule-wide metric scalars
                mol_vec = _build_mol_metrics(match_dicts) if molecule_metrics else np.zeros(0)

                if extra_parts or len(mol_vec):
                    fingerprint = np.concatenate(
                        [base_vec.astype(float)] + [p.astype(float) for p in extra_parts]
                        + ([mol_vec] if len(mol_vec) else [])
                    )
                else:
                    fingerprint = base_vec

                if representation == 'int':
                    fingerprints.append(int(''.join([str(x) for x in fingerprint.astype(int)]), 2))
                else:
                    fingerprints.append(fingerprint)

            elif representation == 'dict':
                fingerprint = {}
                for group in selected_pfas_groups:
                    if group.id in match_dict:
                        match_info = match_dict[group.id]
                        if count_mode == 'binary':
                            fingerprint[group.name] = 1
                        elif count_mode == 'count':
                            fingerprint[group.name] = match_info['match_count']
                        elif count_mode == 'max_component':
                            if match_info['matched_components']:
                                fingerprint[group.name] = max(
                                    comp.get('size', 0) or 0
                                    for comp in match_info['matched_components'])
                            else:
                                fingerprint[group.name] = match_info['match_count'] if match_info['match_count'] > 0 else 0
                        elif count_mode == 'total_component':
                            if match_info['matched_components']:
                                fingerprint[group.name] = sum(
                                    comp.get('size', 0) or 0
                                    for comp in match_info['matched_components'])
                            else:
                                fingerprint[group.name] = match_info['match_count'] if match_info['match_count'] > 0 else 0
                    else:
                        fingerprint[group.name] = 0
                # Append per-group graph metrics to dict
                for gm in graph_metrics:
                    for group in selected_pfas_groups:
                        col = f"{group.name} [{gm}]"
                        if group.id in match_dict:
                            r = _agg_comp_metric(match_dict[group.id].get('matched_components', []), gm)
                        else:
                            r = 0.0
                        fingerprint[col] = r
                # Append molecule metrics to dict
                if molecule_metrics:
                    mv = _build_mol_metrics(match_dicts)
                    for name, val in zip(mol_metric_names, mv):
                        fingerprint[name] = val
                fingerprints.append(fingerprint)

            elif representation == 'sparse':
                fingerprint = {}
                for group in selected_pfas_groups:
                    if group.id in match_dict:
                        match_info = match_dict[group.id]
                        value: float = 0.0
                        if count_mode == 'binary':
                            value = 1.0
                        elif count_mode == 'count':
                            value = float(match_info['match_count'])
                        elif count_mode == 'max_component':
                            if match_info['matched_components']:
                                value = float(max(comp.get('size', 0) or 0
                                                  for comp in match_info['matched_components']))
                            else:
                                value = float(match_info['match_count']) if match_info['match_count'] > 0 else 0.0
                        elif count_mode == 'total_component':
                            if match_info['matched_components']:
                                value = float(sum(comp.get('size', 0) or 0
                                                  for comp in match_info['matched_components']))
                            else:
                                value = float(match_info['match_count']) if match_info['match_count'] > 0 else 0.0
                        if value > 0:
                            fingerprint[group.name] = value
                # Sparse: only non-zero graph metric entries
                for gm in graph_metrics:
                    for group in selected_pfas_groups:
                        if group.id in match_dict:
                            r = _agg_comp_metric(match_dict[group.id].get('matched_components', []), gm)
                            if r != 0.0:
                                fingerprint[f"{group.name} [{gm}]"] = r
                if molecule_metrics:
                    mv = _build_mol_metrics(match_dicts)
                    for name, val in zip(mol_metric_names, mv):
                        if val != 0.0:
                            fingerprint[name] = val
                fingerprints.append(fingerprint)

            elif representation == 'detailed':
                fingerprint = {}
                for group in selected_pfas_groups:
                    if group.id in match_dict:
                        match_info = match_dict[group.id]
                        fingerprint[group.name] = {
                            'group_id': group.id,
                            'match_count': match_info['match_count'],
                            'component_sizes': match_info['component_sizes'],
                            'matched_components': match_info['matched_components'],
                            'smarts1': [Chem.MolToSmarts(x) for x in group.smarts if x is not None],
                        }
                fingerprints.append(fingerprint)
            else:
                raise ValueError(f"Unknown representation: {representation}")

        except Exception as e:
            raise ValueError(f"Error processing SMILES '{smiles_str}': {str(e)}")

    # Return single fingerprint for single input, stacked matrix for multiple
    if single_input:
        return fingerprints[0], group_info
    else:
        if representation == 'vector':
            return np.vstack(fingerprints), group_info
        return fingerprints, group_info
