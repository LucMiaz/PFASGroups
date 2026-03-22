"""PFAS group fingerprints.

:class:`PFASEmbedding` is a :class:`numpy.ndarray` subclass that encapsulates
fingerprint configuration, computation, analysis, and IO. It behaves as an
ordinary NumPy array for indexing, slicing, and arithmetic while carrying
chemical metadata and offering ML helper methods.

Shape convention
----------------
- Single molecule  →  1-D array  ``(n_columns,)``
- Multiple molecules  →  2-D array  ``(n_molecules, n_columns)``

Typical usage
-------------
>>> fp = PFASEmbedding("OC(=O)C(F)(F)F", preset='best')
>>> fp.shape                       # (230,)  — 1-D for single molecule
>>> fp.group_names[:2]             # ['Perfluoroalkyl', ...]

>>> fps = PFASEmbedding(["OC(=O)C(F)(F)F", "OCCS"], preset='best')
>>> fps.shape                      # (2, 230)  — 2-D for multiple molecules

>>> results = parse_smiles(["OC(=O)C(F)(F)F"])
>>> fps = PFASEmbedding(results, preset='best')      # reuses parsed matches
>>> fps = PFASEmbedding(results[0], preset='best')   # single MoleculeResult
"""

from __future__ import annotations

import warnings
from typing import TYPE_CHECKING, Any, Dict, List, Optional, Tuple, Union

import numpy as np
from rdkit import Chem
from rdkit.Chem.rdMolDescriptors import CalcMolFormula

from .parser import parse_groups_in_mol

if TYPE_CHECKING:
    from .PFASEmbeddings import MoleculeResult, ResultsModel

# ---------------------------------------------------------------------------
# Valid values for the component_metrics parameter
# ---------------------------------------------------------------------------

# Count modes: each produces an integer-like scalar per group (presence/count/size).
_COUNT_MODES = {'binary', 'count', 'max_component', 'total_component'}

# Graph metrics: real-valued graph-theoretic descriptors computed per matched
# component.  These are group-dependent because the component subgraph used for
# computation depends on which atoms belong to the matched group.
_GRAPH_METRICS = {
    'branching', 'mean_eccentricity', 'median_eccentricity', 'diameter',
    'radius',
    'effective_graph_resistance',      # uniform weights, original C-skeleton
    'effective_graph_resistance_BDE',  # BDE weights, 1-hop expanded component
    'component_fraction',
    'min_dist_to_center', 'max_dist_to_periphery',
    'min_dist_to_barycenter',
    'size',
    'min_resistance_dist_to_center', 'max_resistance_dist_to_periphery',
    'min_resistance_dist_to_barycenter',
    'n_spacer',    # telomer CH₂ linker chain length (0 for non-telomers)
    'ring_size',   # smallest ring overlapping the component (0 for acyclic groups)
}

# Unified set of all valid component_metrics entries.
_COMPONENT_METRICS = _COUNT_MODES | _GRAPH_METRICS

# Backward-compatible alias (private; kept so that any internal or third-party
# code that imported this name continues to work).
_COMPONENT_GRAPH_METRICS = _GRAPH_METRICS

# ---------------------------------------------------------------------------
# Valid values for the molecule_metrics parameter (molecule-wide aggregates)
# ---------------------------------------------------------------------------

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
        'component_metrics': ['binary', 'effective_graph_resistance'],
        'molecule_metrics': None,
        'description': (
            'Rank 1 — binary + effective_graph_resistance.\n'
            'Mean inter-group Tanimoto = 0.184 (best discrimination, '
            'outperforms TxP-PFAS 129-bit fingerprint).'
        ),
    },
    'best_2': {
        'component_metrics': ['binary', 'branching', 'effective_graph_resistance'],
        'molecule_metrics': None,
        'description': (
            'Rank 2 — binary + branching + effective_graph_resistance.\n'
            'Mean inter-group Tanimoto = 0.189.'
        ),
    },
    'best_3': {
        'component_metrics': ['binary', 'branching', 'mean_eccentricity', 'effective_graph_resistance'],
        'molecule_metrics': None,
        'description': (
            'Rank 3 — binary + branching + mean_eccentricity + effective_graph_resistance.\n'
            'Mean inter-group Tanimoto = 0.196.'
        ),
    },
    'best_4': {
        'component_metrics': ['total_component', 'branching', 'mean_eccentricity', 'effective_graph_resistance'],
        'molecule_metrics': None,
        'description': (
            'Rank 4 — total_component + branching + mean_eccentricity + effective_graph_resistance.\n'
            'Mean inter-group Tanimoto = 0.206.'
        ),
    },
    'best_5': {
        'component_metrics': ['binary', 'branching', 'mean_eccentricity', 'diameter', 'radius',
                              'effective_graph_resistance'],
        'molecule_metrics': None,
        'description': (
            'Rank 5 — binary + 5 graph metrics.\n'
            'Mean inter-group Tanimoto = 0.207.'
        ),
    },
    'binary': {
        'component_metrics': ['binary'],
        'molecule_metrics': None,
        'description': 'Plain binary fingerprint — 1 if group present, 0 absent.',
    },
    'count': {
        'component_metrics': ['count'],
        'molecule_metrics': None,
        'description': 'Count fingerprint — number of matched components per group.',
    },
    'max_component': {
        'component_metrics': ['max_component'],
        'molecule_metrics': None,
        'description': 'Max-component fingerprint — maximum component size (C-atom count) per group.',
    },
}

# Alias: EMBEDDING_PRESETS is the canonical name; FINGERPRINT_PRESETS kept for backward compat.
EMBEDDING_PRESETS = FINGERPRINT_PRESETS
