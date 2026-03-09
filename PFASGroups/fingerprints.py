"""PFAS group fingerprints.

:class:`PFASFingerprint` is a :class:`numpy.ndarray` subclass that encapsulates
fingerprint configuration, computation, analysis, and IO. It behaves as an
ordinary NumPy array for indexing, slicing, and arithmetic while carrying
chemical metadata and offering ML helper methods.

Shape convention
----------------
- Single molecule  →  1-D array  ``(n_columns,)``
- Multiple molecules  →  2-D array  ``(n_molecules, n_columns)``

Typical usage
-------------
>>> fp = PFASFingerprint("OC(=O)C(F)(F)F", preset='best')
>>> fp.shape                       # (230,)  — 1-D for single molecule
>>> fp.group_names[:2]             # ['Perfluoroalkyl', ...]

>>> fps = PFASFingerprint(["OC(=O)C(F)(F)F", "OCCS"], preset='best')
>>> fps.shape                      # (2, 230)  — 2-D for multiple molecules

>>> results = parse_smiles(["OC(=O)C(F)(F)F"])
>>> fps = PFASFingerprint(results, preset='best')      # reuses parsed matches
>>> fps = PFASFingerprint(results[0], preset='best')   # single MoleculeResult
"""

from __future__ import annotations

import warnings
from typing import TYPE_CHECKING, Any, Dict, List, Optional, Tuple, Union

import numpy as np
from rdkit import Chem
from rdkit.Chem.rdMolDescriptors import CalcMolFormula

from .parser import parse_groups_in_mol

if TYPE_CHECKING:
    from .results_model import MoleculeResult, ResultsModel

# ---------------------------------------------------------------------------
# Supported per-group component graph metric keys
# ---------------------------------------------------------------------------

_COMPONENT_GRAPH_METRICS = {
    'branching', 'mean_eccentricity', 'median_eccentricity', 'diameter',
    'radius', 'effective_graph_resistance', 'component_fraction',
    'min_dist_to_center', 'max_dist_to_periphery',
    'min_dist_to_barycenter', 'min_resistance_dist_to_barycenter',
    'min_resistance_dist_to_center', 'max_resistance_dist_to_periphery',
    'size',
}

# Supported molecule-wide metric names
_MOL_METRIC_KEYS = {
    'n_components', 'total_size', 'mean_size', 'max_size',
    'mean_branching', 'max_branching',
    'mean_eccentricity', 'max_diameter',
    'mean_component_fraction', 'max_component_fraction',
}

# ---------------------------------------------------------------------------
# Benchmark-validated fingerprint presets (PFASGroups Tanimoto benchmark 2026)
# Sorted ascending by mean inter-group Tanimoto (lower = more discriminating).
# ---------------------------------------------------------------------------

FINGERPRINT_PRESETS: Dict[str, Any] = {
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

# ---------------------------------------------------------------------------
# PFASFingerprint
# ---------------------------------------------------------------------------


class PFASFingerprint(np.ndarray):
    """PFAS group fingerprint — a NumPy array with chemical metadata and ML methods.

    Behaves exactly like a :class:`numpy.ndarray` for all array operations
    (indexing, arithmetic, broadcasting) while also carrying named group metadata
    and offering dimensionality-reduction, comparison, and IO helper methods.

    Shape convention
    ----------------
    - ``(n_columns,)`` for a single molecule (``str`` or
      :class:`~PFASGroups.results_model.MoleculeResult` source).
    - ``(n_molecules, n_columns)`` for multiple molecules (``list[str]`` or
      :class:`~PFASGroups.results_model.ResultsModel` source).

    Parameters
    ----------
    source : str | list[str] | MoleculeResult | ResultsModel
        Input molecules.  When a :class:`~PFASGroups.results_model.MoleculeResult`
        or :class:`~PFASGroups.results_model.ResultsModel` is given the
        pre-computed group matches are reused—``parse_groups_in_mol`` is *not*
        called again.
    preset : str, optional
        Named benchmark-validated configuration.  Overrides *count_mode*,
        *graph_metrics*, and *molecule_metrics*.  See :data:`FINGERPRINT_PRESETS`.
    count_mode : str, default ``'binary'``
        ``'binary'``, ``'count'``, ``'max_component'``, or ``'total_component'``.
    group_selection : str, default ``'all'``
        Which groups to include: ``'all'``, ``'oecd'``, ``'generic'``,
        ``'telomers'``, or ``'generic+telomers'``.
    selected_group_ids : list[int], optional
        Explicit group IDs (overrides *group_selection*).
    halogens : str | list[str], default ``'F'``
        Which halogen(s) to use for component SMARTS matching.  Multiple
        halogens produce stacked column blocks suffixed with ``[F]``, ``[Cl]``, etc.
    saturation : str | None, default ``'per'``
        ``'per'``, ``'poly'``, or ``None`` (no saturation filter).
    graph_metrics : list[str], optional
        Per-group component graph metrics appended as extra column blocks after
        the base columns.  Each metric adds ``n_base_cols`` more columns.
        Supported: ``'branching'``, ``'mean_eccentricity'``, ``'diameter'``,
        ``'radius'``, ``'effective_graph_resistance'``, and more.
    molecule_metrics : list[str], optional
        Molecule-wide scalar metrics appended as final columns.
        Supported: ``'n_components'``, ``'total_size'``, ``'mean_branching'``, etc.
    pfas_groups : list, optional
        Custom compiled :class:`~PFASGroups.HalogenGroupModel.HalogenGroup`
        instances.  When *None* the default groups are loaded automatically.

    Attributes
    ----------
    smiles : list[str]
        SMILES of each molecule (same order as rows).
    group_names : list[str]
        Column names (includes metric suffixes and ``mol:`` prefix for
        molecule metrics).
    group_selection : str
    count_mode : str
    halogens : list[str]
    saturation : str | None
    graph_metrics : list[str]
    molecule_metrics : list[str]
    preset : str | None

    Examples
    --------
    >>> from PFASGroups import PFASFingerprint, parse_smiles

    >>> # Single molecule — 1-D result
    >>> fp = PFASFingerprint("OC(=O)C(F)(F)C(F)(F)F", preset='best')
    >>> fp.shape                            # (230,)
    >>> fp.group_names[:2]                  # ['Perfluoroalkyl', ...]

    >>> # Batch — 2-D result
    >>> fps = PFASFingerprint(["OC(=O)C(F)(F)F", "OCCS"], preset='best')
    >>> fps.shape                           # (2, 230)

    >>> # From pre-parsed ResultsModel  — skips re-parsing
    >>> results = parse_smiles(["OC(=O)C(F)(F)F"])
    >>> fps = PFASFingerprint(results, preset='best')

    >>> # OECD groups, count mode
    >>> fp_oecd = PFASFingerprint("OC(=O)C(F)(F)F",
    ...                           group_selection='oecd', count_mode='count')
    """

    # Metadata attribute names — propagated by __array_finalize__
    _METADATA_ATTRS: Tuple[str, ...] = (
        'smiles', 'group_names', 'group_selection', 'count_mode',
        'halogens', 'saturation', 'graph_metrics', 'molecule_metrics', 'preset',
    )
    _METADATA_DEFAULTS: Dict[str, Any] = {
        'smiles': [],
        'group_names': [],
        'group_selection': 'all',
        'count_mode': 'binary',
        'halogens': [],
        'saturation': 'per',
        'graph_metrics': [],
        'molecule_metrics': [],
        'preset': None,
    }

    # ------------------------------------------------------------------
    # Construction
    # ------------------------------------------------------------------

    def __new__(
        cls,
        source: Union[str, List[str], 'MoleculeResult', 'ResultsModel'],
        preset: Optional[str] = None,
        count_mode: str = 'binary',
        group_selection: str = 'all',
        selected_group_ids: Optional[List[int]] = None,
        halogens: Union[str, List[str]] = 'F',
        saturation: Optional[str] = 'per',
        graph_metrics: Optional[List[str]] = None,
        molecule_metrics: Optional[List[str]] = None,
        pfas_groups: Optional[list] = None,
    ) -> 'PFASFingerprint':
        # Lazy import avoids circular dependency at module load time.
        from .results_model import MoleculeResult, ResultsModel  # noqa: F401

        # ── Resolve preset ────────────────────────────────────────────
        resolved_cm = count_mode
        resolved_gm: List[str] = list(graph_metrics) if graph_metrics else []
        resolved_mm: List[str] = list(molecule_metrics) if molecule_metrics else []

        if preset is not None:
            if preset not in FINGERPRINT_PRESETS:
                raise ValueError(
                    f"Unknown preset: {preset!r}. "
                    f"Available: {sorted(FINGERPRINT_PRESETS)}"
                )
            _p = FINGERPRINT_PRESETS[preset]
            resolved_cm = _p.get('count_mode', count_mode)
            if _p.get('graph_metrics') is not None:
                resolved_gm = list(_p['graph_metrics'])
            if _p.get('molecule_metrics') is not None:
                resolved_mm = list(_p['molecule_metrics'])

        # ── Validate metric names ─────────────────────────────────────
        bad_gm = [m for m in resolved_gm if m not in _COMPONENT_GRAPH_METRICS]
        if bad_gm:
            raise ValueError(
                f"Unsupported graph_metrics: {bad_gm}. "
                f"Supported: {sorted(_COMPONENT_GRAPH_METRICS)}"
            )
        bad_mm = [m for m in resolved_mm if m not in _MOL_METRIC_KEYS]
        if bad_mm:
            raise ValueError(
                f"Unsupported molecule_metrics: {bad_mm}. "
                f"Supported: {sorted(_MOL_METRIC_KEYS)}"
            )

        # ── Normalise halogens ────────────────────────────────────────
        halogens_list = [halogens] if isinstance(halogens, str) else list(halogens)

        # ── Load groups if not supplied ───────────────────────────────
        if pfas_groups is None:
            from .getter import get_compiled_HalogenGroups
            pfas_groups = get_compiled_HalogenGroups()

        # ── Resolve group selection → 0-based index list ──────────────
        sel_indices = cls._resolve_indices(
            group_selection, selected_group_ids, halogens_list, pfas_groups
        )
        sel_groups = (
            pfas_groups if sel_indices is None
            else [pfas_groups[i] for i in sel_indices]
        )

        # ── Dispatch on source type ───────────────────────────────────
        # NOTE: ResultsModel subclasses list — must be checked before list.
        single = False
        if isinstance(source, str):
            single = True
            smiles_list = [source]
            precomputed: Optional[Dict[str, Any]] = None
        elif isinstance(source, MoleculeResult):
            single = True
            smiles_list = [source.smiles]
            md = cls._match_dict_from_result(source)
            precomputed = {source.smiles: {h: md for h in halogens_list}}
        elif isinstance(source, ResultsModel):
            smiles_list = [mol.smiles for mol in source]
            precomputed = {
                mol.smiles: {h: cls._match_dict_from_result(mol) for h in halogens_list}
                for mol in source
            }
        elif isinstance(source, list):
            smiles_list = source
            precomputed = None
        else:
            raise TypeError(
                f"source must be str, list[str], MoleculeResult, or ResultsModel; "
                f"got {type(source).__name__!r}"
            )

        # ── Compute fingerprint matrix ────────────────────────────────
        matrix, col_names = cls._compute_matrix(
            smiles_list, sel_groups, halogens_list,
            resolved_cm, saturation,
            resolved_gm, resolved_mm,
            precomputed, pfas_groups,
        )

        # ── Build numpy array (1-D for single, 2-D for batch) ─────────
        arr = matrix[0] if single else matrix
        obj = np.asarray(arr, dtype=float).view(cls)

        # ── Attach metadata ───────────────────────────────────────────
        obj.smiles = smiles_list
        obj.group_names = col_names
        obj.group_selection = group_selection
        obj.count_mode = resolved_cm
        obj.halogens = halogens_list
        obj.saturation = saturation
        obj.graph_metrics = resolved_gm
        obj.molecule_metrics = resolved_mm
        obj.preset = preset

        return obj

    def __array_finalize__(self, obj: Optional[np.ndarray]) -> None:
        """Copy metadata when NumPy creates a view/slice/ufunc result."""
        if obj is None:
            return
        for attr, default in self._METADATA_DEFAULTS.items():
            setattr(self, attr, getattr(obj, attr, default))

    @classmethod
    def _from_array(
        cls,
        array: np.ndarray,
        smiles: List[str],
        group_names: List[str],
        group_selection: str = 'all',
        count_mode: str = 'binary',
        halogens: Optional[List[str]] = None,
        saturation: Optional[str] = None,
        graph_metrics: Optional[List[str]] = None,
        molecule_metrics: Optional[List[str]] = None,
    ) -> 'PFASFingerprint':
        """Internal factory — wraps a pre-computed array with metadata.

        Used by :meth:`from_sql` and other deserialisation paths where the
        fingerprint data already exists as a NumPy array.
        """
        _arr = np.asarray(array, dtype=float)
        if _arr.ndim == 2 and _arr.shape[0] != len(smiles):
            raise ValueError(
                f"length mismatch: array has {_arr.shape[0]} rows but "
                f"{len(smiles)} SMILES provided"
            )
        obj = _arr.view(cls)
        obj.smiles = list(smiles)
        obj.group_names = list(group_names)
        obj.group_selection = group_selection
        obj.count_mode = count_mode
        obj.halogens = list(halogens) if halogens else []
        obj.saturation = saturation
        obj.graph_metrics = list(graph_metrics) if graph_metrics else []
        obj.molecule_metrics = list(molecule_metrics) if molecule_metrics else []
        obj.preset = None
        return obj

    # ------------------------------------------------------------------
    # Private computation
    # ------------------------------------------------------------------

    @staticmethod
    def _match_dict_from_result(mol_result: 'MoleculeResult') -> Dict[int, Dict]:
        """Build ``{group_id: {match_count, component_sizes, matched_components}}``
        from a pre-computed :class:`~PFASGroups.results_model.MoleculeResult`.
        """
        result: Dict[int, Dict] = {}
        for match in mol_result.matches:
            if not match.is_group:
                continue
            gid = match.group_id
            if gid is None:
                continue
            components = match.components
            result[gid] = {
                'match_count': match.get('match_count', len(components)),
                'component_sizes': [c.size for c in components],
                'matched_components': [c.data for c in components],
            }
        return result

    @staticmethod
    def _resolve_indices(
        group_selection: str,
        selected_group_ids: Optional[List[int]],
        halogens_list: List[str],
        pfas_groups: list,
    ) -> Optional[List[int]]:
        """Translate *group_selection* / *selected_group_ids* to 0-based
        indices into *pfas_groups*.  Returns ``None`` for ``'all'``.
        """
        id_to_idx = {g.id: i for i, g in enumerate(pfas_groups)}

        if selected_group_ids is not None:
            return [id_to_idx[gid] for gid in selected_group_ids if gid in id_to_idx]

        if group_selection is None or group_selection == 'all':
            return None
        if group_selection == 'oecd':
            return [id_to_idx[gid] for gid in range(1, 29) if gid in id_to_idx]
        if group_selection == 'generic':
            return [id_to_idx[gid] for gid in range(29, 56) if gid in id_to_idx]
        if group_selection == 'telomers':
            return [id_to_idx[gid] for gid in range(74, 117) if gid in id_to_idx]
        if group_selection == 'generic+telomers':
            ids = list(range(29, 56)) + list(range(74, 117))
            return [id_to_idx[gid] for gid in ids if gid in id_to_idx]

        # Dynamic category lookup via raw group JSON
        from .getter import get_HalogenGroups
        raw = get_HalogenGroups()
        compute_raw = [g for g in raw if g.get('compute', True)]
        cats = {g.get('test', {}).get('category', 'other') for g in compute_raw}
        if group_selection in cats:
            matching_ids = {
                g['id'] for g in compute_raw
                if g.get('test', {}).get('category', 'other') == group_selection
            }
            return [id_to_idx[gid] for gid in matching_ids if gid in id_to_idx]

        raise ValueError(
            f"Unknown group_selection: {group_selection!r}. "
            f"Choose from: 'all', 'oecd', 'generic', 'telomers', 'generic+telomers'"
        )

    @staticmethod
    def _agg_metric(comps: List[Dict], metric: str) -> float:
        """Mean of *metric* over matched components (0.0 when none found)."""
        vals = [c[metric] for c in comps if c.get(metric) is not None]
        return float(np.mean(vals)) if vals else 0.0

    @staticmethod
    def _encode_count(match_info: Dict, count_mode: str) -> float:
        """Scalar encoding for one matched group given *count_mode*."""
        if count_mode == 'binary':
            return 1.0
        if count_mode == 'count':
            return float(match_info['match_count'])
        comps = match_info.get('matched_components', [])
        sizes = [c.get('size', 0) or 0 for c in comps]
        if count_mode == 'max_component':
            return float(max(sizes)) if sizes else float(match_info['match_count'])
        if count_mode == 'total_component':
            return float(sum(sizes)) if sizes else float(match_info['match_count'])
        raise ValueError(f"Unknown count_mode: {count_mode!r}")

    @classmethod
    def _mol_metric_vec(
        cls,
        match_dicts: Dict[str, Dict],
        mol_metrics: List[str],
    ) -> np.ndarray:
        """Molecule-wide metric vector from all matched components across halogens."""
        all_comps = [
            c
            for md in match_dicts.values()
            for mi in md.values()
            for c in mi.get('matched_components', [])
        ]
        result: List[float] = []
        for m in mol_metrics:
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
                result.append(cls._agg_metric(all_comps, 'branching'))
            elif m == 'max_branching':
                vals = [c.get('branching') for c in all_comps if c.get('branching') is not None]
                result.append(float(max(vals)) if vals else 0.0)
            elif m == 'mean_eccentricity':
                result.append(cls._agg_metric(all_comps, 'mean_eccentricity'))
            elif m == 'max_diameter':
                vals = [c.get('diameter') for c in all_comps if c.get('diameter') is not None]
                result.append(float(max(vals)) if vals else 0.0)
            elif m == 'mean_component_fraction':
                result.append(cls._agg_metric(all_comps, 'component_fraction'))
            elif m == 'max_component_fraction':
                vals = [c.get('component_fraction') for c in all_comps if c.get('component_fraction') is not None]
                result.append(float(max(vals)) if vals else 0.0)
            else:
                result.append(0.0)
        return np.array(result, dtype=float)

    @classmethod
    def _compute_matrix(
        cls,
        smiles_list: List[str],
        sel_groups: list,
        halogens_list: List[str],
        count_mode: str,
        saturation: Optional[str],
        graph_metrics: List[str],
        mol_metrics: List[str],
        precomputed: Optional[Dict[str, Any]],
        pfas_groups: list,
    ) -> Tuple[np.ndarray, List[str]]:
        """Compute the ``(n_molecules, n_columns)`` matrix and column names."""
        multi_hal = len(halogens_list) > 1

        # ── Column names ──────────────────────────────────────────────
        if multi_hal:
            base_names = [f"{g.name} [{h}]" for h in halogens_list for g in sel_groups]
        else:
            base_names = [g.name for g in sel_groups]

        n_base = len(base_names)
        gm_names: List[str] = [
            f"{name} [{gm}]"
            for gm in graph_metrics
            for name in base_names
        ]
        mm_names = [f"mol:{m}" for m in mol_metrics]
        col_names = base_names + gm_names + mm_names
        n_cols = len(col_names)

        matrix = np.zeros((len(smiles_list), n_cols), dtype=float)

        for row_idx, smi in enumerate(smiles_list):
            try:
                mol = Chem.MolFromSmiles(smi)
                if mol is None:
                    raise ValueError(f"Invalid SMILES: {smi!r}")
                formula = CalcMolFormula(mol)

                # Collect match dicts per halogen (precomputed or fresh)
                precomp_smi = (precomputed or {}).get(smi, {})
                match_dicts: Dict[str, Dict[int, Dict]] = {}
                for h in halogens_list:
                    if h in precomp_smi:
                        match_dicts[h] = precomp_smi[h]
                    else:
                        kw: Dict[str, Any] = {'halogens': h}
                        if saturation is not None:
                            kw['saturation'] = saturation
                        all_matches, _ = parse_groups_in_mol(
                            mol, formula=formula, pfas_groups=pfas_groups,
                            include_PFAS_definitions=False, **kw
                        )
                        match_dicts[h] = {
                            g.id: {
                                'match_count': mc,
                                'component_sizes': sizes,
                                'matched_components': comps,
                            }
                            for g, mc, sizes, comps in all_matches
                        }

                # ── Base columns ──────────────────────────────────────
                if multi_hal:
                    for hi, h in enumerate(halogens_list):
                        off = hi * len(sel_groups)
                        for gi, g in enumerate(sel_groups):
                            if g.id in match_dicts[h]:
                                matrix[row_idx, off + gi] = cls._encode_count(
                                    match_dicts[h][g.id], count_mode
                                )
                else:
                    md0 = match_dicts[halogens_list[0]]
                    for gi, g in enumerate(sel_groups):
                        if g.id in md0:
                            matrix[row_idx, gi] = cls._encode_count(md0[g.id], count_mode)

                # ── Per-group graph metric blocks ─────────────────────
                for bi, gm in enumerate(graph_metrics):
                    col_off = n_base + bi * n_base
                    for ci in range(n_base):
                        if multi_hal:
                            hi = ci // len(sel_groups)
                            gi = ci % len(sel_groups)
                            g = sel_groups[gi]
                            md = match_dicts[halogens_list[hi]]
                        else:
                            gi = ci
                            g = sel_groups[gi]
                            md = match_dicts[halogens_list[0]]
                        if g.id in md:
                            comps = md[g.id].get('matched_components', [])
                            matrix[row_idx, col_off + ci] = cls._agg_metric(comps, gm)

                # ── Molecule-level metric columns ─────────────────────
                if mol_metrics:
                    mm_start = n_base + len(graph_metrics) * n_base
                    matrix[row_idx, mm_start:] = cls._mol_metric_vec(match_dicts, mol_metrics)

            except Exception as exc:
                raise ValueError(f"Error processing SMILES {smi!r}: {exc}") from exc

        return matrix, col_names

    # ------------------------------------------------------------------
    # Public properties
    # ------------------------------------------------------------------

    @property
    def n_molecules(self) -> int:
        """Number of molecules (always 1 for a 1-D fingerprint)."""
        return 1 if self.ndim == 1 else self.shape[0]

    @property
    def n_columns(self) -> int:
        """Number of fingerprint columns."""
        return self.shape[-1]

    @property
    def fingerprints(self) -> np.ndarray:
        """Backward-compatible 2-D view of the underlying NumPy data.

        Always has shape ``(n_molecules, n_columns)`` regardless of whether
        this instance is 1-D or 2-D, so that code written against the old
        ``ResultsFingerprint.fingerprints`` attribute continues to work.
        """
        return np.atleast_2d(np.asarray(self))

    # ------------------------------------------------------------------
    # Display
    # ------------------------------------------------------------------

    def __repr__(self) -> str:
        hal_str = '+'.join(self.halogens) if self.halogens else 'all'
        parts: List[str] = []
        if self.preset:
            parts.append(f"preset={self.preset!r}")
        parts += [
            f"shape={self.shape}",
            f"count_mode={self.count_mode!r}",
            f"group_selection={self.group_selection!r}",
            f"halogens={hal_str!r}",
        ]
        if self.graph_metrics:
            parts.append(f"graph_metrics={self.graph_metrics!r}")
        if self.molecule_metrics:
            parts.append(f"molecule_metrics={self.molecule_metrics!r}")
        return f"PFASFingerprint({', '.join(parts)})"

    def _compute_group_distribution(self) -> np.ndarray:
        """Return per-group presence count (number of molecules that have each group > 0).

        Returns
        -------
        np.ndarray of shape (n_columns,)
        """
        mat = np.atleast_2d(np.asarray(self))
        return np.sum(mat > 0, axis=0).astype(float)

    def summary(self) -> str:
        """Return a human-readable text summary of fingerprint statistics."""
        hal_str = '+'.join(self.halogens) if self.halogens else 'all'
        mat = np.atleast_2d(np.asarray(self))
        lines = [
            "PFASFingerprint Summary",
            "=" * 50,
            f"Molecules  : {self.n_molecules}",
            f"Columns    : {self.n_columns}",
            f"Group sel. : {self.group_selection}",
            f"Halogens   : {hal_str}",
            f"Saturation : {self.saturation}",
            f"Count mode : {self.count_mode}",
        ]
        if self.preset:
            lines.append(f"Preset     : {self.preset}")
        if self.graph_metrics:
            lines.append(f"Graph metrics    : {', '.join(self.graph_metrics)}")
        if self.molecule_metrics:
            lines.append(f"Molecule metrics : {', '.join(self.molecule_metrics)}")
        lines += [
            f"Shape      : {mat.shape}",
            f"Non-zero   : {int(np.count_nonzero(mat))}",
            f"Sparsity   : {1 - np.count_nonzero(mat) / max(mat.size, 1):.2%}",
        ]

        # Most frequent base groups (skip metric-extension columns)
        base_idx = [
            i for i, n in enumerate(self.group_names)
            if '[' not in n and not n.startswith('mol:')
        ] or list(range(self.n_columns))

        if mat.shape[0] > 0 and base_idx:
            freq = np.sum(mat[:, base_idx] > 0, axis=0)
            top = np.argsort(freq)[::-1][:10]
            lines.append("\nMost common groups:")
            for rel in top:
                if freq[rel] > 0:
                    abs_i = base_idx[rel]
                    lines.append(
                        f"  {self.group_names[abs_i]}: {int(freq[rel])} molecules "
                        f"({100 * freq[rel] / self.n_molecules:.1f}%)"
                    )
        return "\n".join(lines)

    # ------------------------------------------------------------------
    # ML analysis
    # ------------------------------------------------------------------

    def perform_pca(
        self,
        n_components: int = 2,
        plot: bool = True,
        output_file: Optional[str] = None,
    ) -> Dict[str, Any]:
        """Perform PCA on the fingerprint matrix.

        Parameters
        ----------
        n_components : int, default 2
        plot : bool, default True
        output_file : str, optional

        Returns
        -------
        dict
            Keys: ``'transformed'``, ``'explained_variance'``, ``'components'``,
            ``'pca_model'``, ``'scaler'``.
        """
        try:
            from sklearn.decomposition import PCA
            from sklearn.preprocessing import StandardScaler
            import matplotlib.pyplot as plt
        except ImportError as exc:
            raise ImportError(
                "scikit-learn and matplotlib required: pip install scikit-learn matplotlib"
            ) from exc

        mat = np.atleast_2d(np.asarray(self))
        scaler = StandardScaler()
        X = scaler.fit_transform(mat)
        pca = PCA(n_components=n_components)
        X_pca = pca.fit_transform(X)

        if plot and n_components >= 2:
            fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))
            ax1.scatter(X_pca[:, 0], X_pca[:, 1], alpha=0.6, s=50)
            ax1.set_xlabel(f'PC1 ({pca.explained_variance_ratio_[0]:.1%} variance)')
            ax1.set_ylabel(f'PC2 ({pca.explained_variance_ratio_[1]:.1%} variance)')
            ax1.set_title('PCA of PFAS Group Fingerprints')
            ax1.grid(True, alpha=0.3)
            ax2.bar(range(1, n_components + 1), pca.explained_variance_ratio_)
            ax2.set_xlabel('Principal Component')
            ax2.set_ylabel('Explained Variance Ratio')
            ax2.set_title('Scree Plot')
            ax2.grid(True, alpha=0.3)
            plt.tight_layout()
            if output_file:
                plt.savefig(output_file, dpi=300, bbox_inches='tight')
            else:
                plt.show()
            plt.close()

        return {
            'transformed': X_pca,
            'explained_variance': pca.explained_variance_ratio_,
            'components': pca.components_,
            'pca_model': pca,
            'scaler': scaler,
        }

    def perform_kernel_pca(
        self,
        n_components: int = 2,
        kernel: str = 'rbf',
        gamma: Optional[float] = None,
        plot: bool = True,
        output_file: Optional[str] = None,
    ) -> Dict[str, Any]:
        """Perform kernel PCA.

        Parameters
        ----------
        n_components : int, default 2
        kernel : str, default ``'rbf'``
            ``'linear'``, ``'poly'``, ``'rbf'``, ``'sigmoid'``, ``'cosine'``.
        gamma : float, optional
        plot : bool, default True
        output_file : str, optional

        Returns
        -------
        dict
            Keys: ``'transformed'``, ``'kpca_model'``, ``'scaler'``,
            ``'kernel'``, ``'gamma'``.
        """
        try:
            from sklearn.decomposition import KernelPCA
            from sklearn.preprocessing import StandardScaler
            import matplotlib.pyplot as plt
        except ImportError as exc:
            raise ImportError(
                "scikit-learn and matplotlib required: pip install scikit-learn matplotlib"
            ) from exc

        mat = np.atleast_2d(np.asarray(self))
        scaler = StandardScaler()
        X = scaler.fit_transform(mat)
        gamma = gamma if gamma is not None else 1.0 / X.shape[1]
        kpca = KernelPCA(n_components=n_components, kernel=kernel, gamma=gamma)
        X_kpca = kpca.fit_transform(X)

        if plot and n_components >= 2:
            fig, ax = plt.subplots(figsize=(8, 6))
            ax.scatter(X_kpca[:, 0], X_kpca[:, 1], alpha=0.6, s=50)
            ax.set_xlabel('Kernel PC1')
            ax.set_ylabel('Kernel PC2')
            ax.set_title(f'Kernel PCA ({kernel}) of PFAS Group Fingerprints')
            ax.grid(True, alpha=0.3)
            plt.tight_layout()
            if output_file:
                plt.savefig(output_file, dpi=300, bbox_inches='tight')
            else:
                plt.show()
            plt.close()

        return {
            'transformed': X_kpca,
            'kpca_model': kpca,
            'scaler': scaler,
            'kernel': kernel,
            'gamma': gamma,
        }

    def perform_tsne(
        self,
        n_components: int = 2,
        perplexity: float = 30.0,
        learning_rate: float = 200.0,
        max_iter: int = 1000,
        plot: bool = True,
        output_file: Optional[str] = None,
    ) -> Dict[str, Any]:
        """Perform t-SNE dimensionality reduction.

        Parameters
        ----------
        n_components : int, default 2
        perplexity : float, default 30.0
        learning_rate : float, default 200.0
        max_iter : int, default 1000
        plot : bool, default True
        output_file : str, optional

        Returns
        -------
        dict
            Keys: ``'transformed'``, ``'tsne_model'``, ``'scaler'``, ``'perplexity'``.
        """
        try:
            from sklearn.manifold import TSNE
            from sklearn.preprocessing import StandardScaler
            import matplotlib.pyplot as plt
        except ImportError as exc:
            raise ImportError(
                "scikit-learn and matplotlib required: pip install scikit-learn matplotlib"
            ) from exc

        mat = np.atleast_2d(np.asarray(self))
        scaler = StandardScaler()
        X = scaler.fit_transform(mat)
        tsne = TSNE(
            n_components=n_components,
            perplexity=perplexity,
            learning_rate=learning_rate,
            max_iter=max_iter,
            random_state=42,
        )
        X_tsne = tsne.fit_transform(X)

        if plot and n_components >= 2:
            fig, ax = plt.subplots(figsize=(8, 6))
            ax.scatter(X_tsne[:, 0], X_tsne[:, 1], alpha=0.6, s=50)
            ax.set_xlabel('t-SNE 1')
            ax.set_ylabel('t-SNE 2')
            ax.set_title(f't-SNE (perplexity={perplexity}) of PFAS Group Fingerprints')
            ax.grid(True, alpha=0.3)
            plt.tight_layout()
            if output_file:
                plt.savefig(output_file, dpi=300, bbox_inches='tight')
            else:
                plt.show()
            plt.close()

        return {
            'transformed': X_tsne,
            'tsne_model': tsne,
            'scaler': scaler,
            'perplexity': perplexity,
        }

    def perform_umap(
        self,
        n_components: int = 2,
        n_neighbors: int = 15,
        min_dist: float = 0.1,
        metric: str = 'euclidean',
        plot: bool = True,
        output_file: Optional[str] = None,
    ) -> Dict[str, Any]:
        """Perform UMAP dimensionality reduction.

        Parameters
        ----------
        n_components : int, default 2
        n_neighbors : int, default 15
        min_dist : float, default 0.1
        metric : str, default ``'euclidean'``
        plot : bool, default True
        output_file : str, optional

        Returns
        -------
        dict
            Keys: ``'transformed'``, ``'umap_model'``, ``'scaler'``,
            ``'n_neighbors'``, ``'min_dist'``.
        """
        try:
            import umap
            from sklearn.preprocessing import StandardScaler
            import matplotlib.pyplot as plt
        except ImportError as exc:
            raise ImportError(
                "umap-learn and matplotlib required: pip install umap-learn matplotlib"
            ) from exc

        import os as _os
        _os.environ.setdefault('KMP_DUPLICATE_LIB_OK', 'TRUE')

        mat = np.atleast_2d(np.asarray(self))
        scaler = StandardScaler()
        X = scaler.fit_transform(mat)

        with warnings.catch_warnings():
            warnings.filterwarnings('ignore', message=r'.*n_jobs.*overridden.*', category=UserWarning)
            warnings.filterwarnings('ignore', message=r'.*Intel OpenMP.*LLVM OpenMP.*', category=RuntimeWarning)
            reducer = umap.UMAP(
                n_components=n_components,
                n_neighbors=n_neighbors,
                min_dist=min_dist,
                metric=metric,
                random_state=42,
            )
            X_umap = reducer.fit_transform(X)

        if plot and n_components >= 2:
            fig, ax = plt.subplots(figsize=(8, 6))
            ax.scatter(X_umap[:, 0], X_umap[:, 1], alpha=0.6, s=50)
            ax.set_xlabel('UMAP 1')
            ax.set_ylabel('UMAP 2')
            ax.set_title(f'UMAP (n_neighbors={n_neighbors}) of PFAS Group Fingerprints')
            ax.grid(True, alpha=0.3)
            plt.tight_layout()
            if output_file:
                plt.savefig(output_file, dpi=300, bbox_inches='tight')
            else:
                plt.show()
            plt.close()

        return {
            'transformed': X_umap,
            'umap_model': reducer,
            'scaler': scaler,
            'n_neighbors': n_neighbors,
            'min_dist': min_dist,
        }

    def compare_kld(
        self,
        other: 'PFASFingerprint',
        method: str = 'minmax',
    ) -> float:
        """Compare two fingerprint sets using KL divergence.

        Implements the minmaxKLd method for comparing PFAS group fingerprint
        distributions between datasets.

        Parameters
        ----------
        other : PFASFingerprint
        method : str, default ``'minmax'``
            ``'minmax'``, ``'forward'``, ``'reverse'``, or ``'symmetric'``.

        Returns
        -------
        float
            KL divergence value (lower = more similar).
        """
        from scipy.stats import entropy as _entropy

        if self.n_columns != other.n_columns:
            raise ValueError(
                f"Column count mismatch: {self.n_columns} vs {other.n_columns}"
            )

        eps = 1e-10
        p = np.sum(np.atleast_2d(np.asarray(self)) > 0, axis=0).astype(float) + eps
        q = np.sum(np.atleast_2d(np.asarray(other)) > 0, axis=0).astype(float) + eps
        p /= p.sum()
        q /= q.sum()

        if method == 'forward':
            return float(_entropy(p, q))
        if method == 'reverse':
            return float(_entropy(q, p))
        if method == 'symmetric':
            return float((_entropy(p, q) + _entropy(q, p)) / 2)
        if method == 'minmax':
            kl_fwd = _entropy(p, q)
            kl_rev = _entropy(q, p)
            kl_sym = (kl_fwd + kl_rev) / 2
            max_kl = np.log(len(p))
            return float(kl_sym / max_kl) if max_kl > 0 else 0.0
        raise ValueError(f"Unknown method: {method!r}")

    # ------------------------------------------------------------------
    # IO
    # ------------------------------------------------------------------

    def to_sql(
        self,
        conn: Optional[Any] = None,
        filename: Optional[str] = None,
        table_name: str = "fingerprints",
        metadata_table: str = "fingerprint_metadata",
        if_exists: str = "append",
    ) -> None:
        """Save fingerprints to a SQL database (only non-zero entries stored).

        Parameters
        ----------
        conn : str | SQLAlchemy Engine, optional
        filename : str, optional
            SQLite file path (alternative to *conn*).
        table_name : str, default ``"fingerprints"``
        metadata_table : str, default ``"fingerprint_metadata"``
        if_exists : str, default ``"append"``
        """
        try:
            import sqlalchemy
            import pandas as pd
        except ImportError as exc:
            raise ImportError(
                "sqlalchemy and pandas required: pip install sqlalchemy pandas"
            ) from exc

        if conn is not None:
            engine = sqlalchemy.create_engine(conn) if isinstance(conn, str) else conn
        elif filename is not None:
            engine = sqlalchemy.create_engine(f"sqlite:///{filename}")
        else:
            raise ValueError("Either conn or filename must be provided")

        mat = np.atleast_2d(np.asarray(self))
        rows = [
            {'smiles': smi, 'group_name': self.group_names[j],
             'group_index': j, 'value': float(mat[i, j])}
            for i, smi in enumerate(self.smiles)
            for j in range(self.n_columns)
            if mat[i, j] > 0
        ]
        pd.DataFrame(rows).to_sql(table_name, engine, if_exists=if_exists, index=False)
        pd.DataFrame([{
            'group_selection': self.group_selection,
            'count_mode': self.count_mode,
            'n_molecules': self.n_molecules,
            'n_groups': self.n_columns,
        }]).to_sql(metadata_table, engine, if_exists=if_exists, index=False)

    @classmethod
    def from_sql(
        cls,
        conn: Optional[Any] = None,
        filename: Optional[str] = None,
        table_name: str = "fingerprints",
        metadata_table: str = "fingerprint_metadata",
        limit: Optional[int] = None,
    ) -> 'PFASFingerprint':
        """Load fingerprints from a SQL database.

        Parameters
        ----------
        conn : str | SQLAlchemy Engine, optional
        filename : str, optional
        table_name : str, default ``"fingerprints"``
        metadata_table : str, default ``"fingerprint_metadata"``
        limit : int, optional

        Returns
        -------
        PFASFingerprint
        """
        try:
            import sqlalchemy
            import pandas as pd
        except ImportError as exc:
            raise ImportError(
                "sqlalchemy and pandas required: pip install sqlalchemy pandas"
            ) from exc

        if conn is not None:
            engine = sqlalchemy.create_engine(conn) if isinstance(conn, str) else conn
        elif filename is not None:
            engine = sqlalchemy.create_engine(f"sqlite:///{filename}")
        else:
            raise ValueError("Either conn or filename must be provided")

        meta = pd.read_sql(f"SELECT * FROM {metadata_table} LIMIT 1", engine)
        group_selection = meta['group_selection'].iloc[0]
        count_mode = meta['count_mode'].iloc[0]

        query = f"SELECT * FROM {table_name}"
        if limit is not None:
            query += f" LIMIT {limit}"
        df = pd.read_sql(query, engine)

        smiles_list = sorted(df['smiles'].unique())
        group_names = sorted(df['group_name'].unique())
        matrix = np.zeros((len(smiles_list), len(group_names)), dtype=float)
        for _, row in df.iterrows():
            i = smiles_list.index(row['smiles'])
            j = group_names.index(row['group_name'])
            matrix[i, j] = row['value']

        return cls._from_array(
            matrix,
            smiles=smiles_list,
            group_names=group_names,
            group_selection=group_selection,
            count_mode=count_mode,
        )


# ---------------------------------------------------------------------------
# Backward-compatible generate_fingerprint() function
# ---------------------------------------------------------------------------

def generate_fingerprint(
    smiles: Union[str, List[str]],
    selected_groups: Union[List[int], range, None] = None,
    representation: str = 'vector',
    count_mode: str = 'binary',
    halogens: Union[str, List[str]] = 'F',
    saturation: Optional[str] = 'per',
    graph_metrics: Optional[List[str]] = None,
    molecule_metrics: Optional[List[str]] = None,
    preset: Optional[str] = None,
    **kwargs,
) -> Tuple[np.ndarray, Dict[str, Any]]:
    """Generate PFAS group fingerprints from SMILES strings.

    Backward-compatible wrapper around :class:`PFASFingerprint`.
    Returns ``(fp_array, group_info)`` where ``group_info`` is a dict
    containing ``group_names``, ``count_mode``, ``halogens``, and
    ``saturation``.

    Parameters
    ----------
    smiles : str or list of str
    selected_groups : list of int or None
        Group IDs to include (``None`` = all groups).
    representation : str, default ``'vector'``
        Only ``'vector'`` is supported; other values are accepted but
        ignored (the result is always a NumPy array).
    count_mode : str, default ``'binary'``
    halogens : str or list of str, default ``'F'``
    saturation : str or None, default ``'per'``
    graph_metrics : list of str or None
    molecule_metrics : list of str or None
    preset : str or None
    **kwargs
        Forwarded to :class:`PFASFingerprint`.

    Returns
    -------
    fp_array : numpy.ndarray
        Shape ``(n_columns,)`` for a single SMILES or
        ``(n_molecules, n_columns)`` for a list.
    group_info : dict
        ``{'group_names': [...], 'count_mode': ..., 'halogens': [...],
        'saturation': ..., 'selected_group_ids': [...] | None}``
    """
    fp = PFASFingerprint(
        smiles,
        preset=preset,
        count_mode=count_mode,
        selected_group_ids=list(selected_groups) if selected_groups is not None else None,
        halogens=halogens,
        saturation=saturation,
        graph_metrics=graph_metrics,
        molecule_metrics=molecule_metrics,
        **kwargs,
    )
    group_info: Dict[str, Any] = {
        'group_names': fp.group_names,
        'count_mode': fp.count_mode,
        'halogens': fp.halogens,
        'saturation': fp.saturation,
        'selected_group_ids': list(selected_groups) if selected_groups is not None else None,
    }
    return np.asarray(fp), group_info
