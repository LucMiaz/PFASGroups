"""Structured feature extraction from PFASGroups halogen group detection.

This module provides :func:`extract_group_features`, a single public entry
point that runs the necessary :func:`~PFASGroups.parser.parse_mol` calls and
returns a :class:`GroupFeatureResult` with four named feature dictionaries
covering the following four group categories:

- **Polyhalogenated groups** (ids 35, 38, 45): aggregate match counts
- **Perhalogenated groups** (ids 34, 37, 44): per-halogen max component sizes
- **Alkyl chain (H pseudo-halogen)** aggregates: convenience view of the H
  column of perhalogenated results
- **Generic functional groups** (ids 29-76): wildcard bag-of-groups counts

The six component groups (34/35/37/38/44/45) are excluded from wildcard
matching by design, so they never appear in the generic-group counts.

Example
-------
>>> from rdkit import Chem
>>> from PFASGroups import extract_group_features
>>> mol = Chem.MolFromSmiles("FC(F)(F)C(F)(F)C(F)(F)C(=O)O")  # PFOA-like
>>> r = extract_group_features(mol)
>>> r.poly_counts
{'poly_alkyl': 1.0, 'poly_aryl': 0.0, 'poly_cyclic': 0.0}
>>> r.per_halogen_sizes['g34_F']  # longest perfluoroalkyl component (carbons)
4.0
>>> r.h_chain_sizes
{'alkyl_H': 0.0, 'aryl_H': 0.0, 'cyclic_H': 0.0}
>>> len(r.to_array())  # 3 + 15 + 48
66
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Union

import numpy as np
from rdkit import Chem

from .parser import parse_mol as _parse_mol  # type: ignore[attr-defined]
from .getter import get_compiled_HalogenGroups as _get_compiled_HalogenGroups

# ---------------------------------------------------------------------------
# Public constants
# ---------------------------------------------------------------------------

#: Group IDs for perhalogenated groups (alkyl, aryl, cyclic).
#: These groups require *all* halogen-bearing carbons to carry only one halogen
#: type.  Each match is attributed to a single halogen (or H for alkyl chains).
PER_GROUP_IDS: list[int] = [34, 37, 44]

#: Group IDs for polyhalogenated groups (alkyl, aryl, cyclic).
#: These groups allow a *mixture* of halogens and/or partial halogenation.
#: Match counts are aggregated across all halogens.
POLY_GROUP_IDS: list[int] = [35, 38, 45]

#: Halogen types recognised in :attr:`GroupFeatureResult.per_halogen_sizes`.
#: ``'H'`` is the pseudo-halogen used to detect un-substituted alkyl chains
#: (perhalogenated-alkyl, group 34) via :class:`~PFASGroups.ComponentsSolverModel.ComponentsSolver`.
HALOGENS_ORDER: list[str] = ["F", "Cl", "Br", "I", "H"]

#: Ordered list of group IDs covered by :attr:`GroupFeatureResult.generic_groups`.
#: Ids {34, 35, 37, 38, 44, 45} are excluded from wildcard matching and will
#: always be zero in this vector.
GENERIC_GROUP_VOCAB: list[int] = list(range(29, 77))  # 48 entries

#: Human-readable names for every group id in :data:`GENERIC_GROUP_VOCAB`.
GENERIC_GROUP_NAMES: dict[int, str] = {
    29: "acrylate",
    30: "acyl halide",
    31: "alcohol",
    32: "aldehyde",
    33: "alkene",
    34: "perhalogenated alkyl",
    35: "polyhalogenated alkyl",
    36: "alkyne",
    37: "perhalogenated aryl compounds",
    38: "polyhalogenated aryl compounds",
    39: "benzodioxole",
    40: "benzoyl peroxydes",
    41: "bromide",
    42: "carboxylic acid",
    43: "chloride",
    44: "perhalogenated cyclic compounds",
    45: "polyhalogenated cyclic compounds",
    46: "ester",
    47: "ether",
    48: "fluoride",
    49: "glucuronate",
    50: "iodide",
    51: "ketone",
    52: "methacrylate",
    53: "peroxydes",
    54: "side-chain aromatics",
    55: "sulfenic acid",
    56: "sulfenyl halide",
    57: "sulfinic acid",
    58: "sulfinyl amido sulfonic acid",
    59: "sulfonamide",
    60: "sulfonamidoethanol",
    61: "sulfonic acid",
    62: "sulfonyl halide",
    63: "sulfonyl propanoic acid",
    64: "sulfuric acid",
    65: "thioester keto dicarboxylic acid",
    66: "thiocyanic acid",
    67: "phosphinic acid",
    68: "phosphonic acid",
    69: "amide",
    70: "amine",
    71: "heterocyclic azine",
    72: "heterocyclic azole",
    73: "betaine",
    74: "glycine",
    75: "trichlorosilane",
    76: "silane",
}

# ---------------------------------------------------------------------------
# Internal state
# ---------------------------------------------------------------------------

_COMPONENT_GROUPS: list | None = None  # cached list of HalogenGroups {34..45}
_GENERIC_GROUP_INDEX: dict[int, int] = {gid: i for i, gid in enumerate(GENERIC_GROUP_VOCAB)}
_HAL_TO_COL: dict[str, int] = {h: i for i, h in enumerate(HALOGENS_ORDER)}
_PER_GID_TO_ROW: dict[int, int] = {34: 0, 37: 1, 44: 2}
_POLY_GID_TO_IDX: dict[int, int] = {35: 0, 38: 1, 45: 2}

# Short labels used in h_chain_sizes
_H_CHAIN_KEYS: list[str] = ["alkyl_H", "aryl_H", "cyclic_H"]
# Corresponding per_halogen_sizes keys (H column of each per group)
_H_CHAIN_PHS_KEYS: list[str] = ["g34_H", "g37_H", "g44_H"]


def _get_component_groups() -> list:
    """Return a cached list of :class:`~PFASGroups.HalogenGroup` objects with ids in ``{34,35,37,38,44,45}``."""
    global _COMPONENT_GROUPS
    if _COMPONENT_GROUPS is None:
        _ids = {34, 35, 37, 38, 44, 45}
        _COMPONENT_GROUPS = [g for g in _get_compiled_HalogenGroups() if g.id in _ids]
    return _COMPONENT_GROUPS


# ---------------------------------------------------------------------------
# GroupFeatureResult
# ---------------------------------------------------------------------------

@dataclass
class GroupFeatureResult:
    """Structured feature extraction result for a single molecule.

    This dataclass groups features from two :func:`~PFASGroups.parser.parse_mol`
    calls into four semantically distinct dictionaries.  It is returned by
    :func:`extract_group_features`.

    Attributes
    ----------
    poly_counts : dict[str, float]
        Match counts for **polyhalogenated** groups (ids 35, 38, 45), aggregated
        across all halogens (F, Cl, Br, I, H).  Keys: ``'poly_alkyl'``,
        ``'poly_aryl'``, ``'poly_cyclic'``.

    per_halogen_sizes : dict[str, float]
        Maximum carbon-component size for **perhalogenated** groups (ids 34, 37,
        44) resolved per halogen.  Keys follow the pattern ``'g{id}_{hal}'``
        where *id* is in ``{34, 37, 44}`` and *hal* is in
        ``['F', 'Cl', 'Br', 'I', 'H']``.  Zero when no match is found.

    h_chain_sizes : dict[str, float]
        Convenience view of the ``'H'`` column of :attr:`per_halogen_sizes`:
        the largest un-substituted alkyl / aryl / cyclic component detected via
        the H pseudo-halogen mechanism.  Keys: ``'alkyl_H'``, ``'aryl_H'``,
        ``'cyclic_H'``.  **Not included in** :meth:`to_array`.

    generic_groups : dict[str, float]
        Wildcard functional-group match counts for group ids 29-76.  Keys:
        ``'g29'`` … ``'g76'``.  Ids ``{34, 35, 37, 38, 44, 45}`` are excluded
        from wildcard matching and will always be zero.
    """

    poly_counts: dict[str, float] = field(default_factory=dict)
    per_halogen_sizes: dict[str, float] = field(default_factory=dict)
    h_chain_sizes: dict[str, float] = field(default_factory=dict)
    generic_groups: dict[str, float] = field(default_factory=dict)

    def to_array(self) -> np.ndarray:
        """Return a fixed-length float32 array of shape ``(66,)``.

        Layout::

            [0:3]   poly_counts    (poly_alkyl, poly_aryl, poly_cyclic)
            [3:18]  per_halogen_sizes  (g34_F…g44_H, 3 groups × 5 halogens)
            [18:66] generic_groups  (g29…g76, 48 entries)

        Note
        ----
        :attr:`h_chain_sizes` is excluded from this array (it is a strict
        subset of :attr:`per_halogen_sizes`).
        """
        poly = np.array(
            [self.poly_counts[k] for k in ("poly_alkyl", "poly_aryl", "poly_cyclic")],
            dtype=np.float32,
        )
        per_keys = [
            f"g{gid}_{h}"
            for gid in PER_GROUP_IDS
            for h in HALOGENS_ORDER
        ]
        per = np.array([self.per_halogen_sizes[k] for k in per_keys], dtype=np.float32)
        gen = np.array(
            [self.generic_groups[f"g{gid}"] for gid in GENERIC_GROUP_VOCAB],
            dtype=np.float32,
        )
        return np.concatenate([poly, per, gen])

    def feature_names(self) -> list[str]:
        """Return the 66 feature names in the same order as :meth:`to_array`.

        Returns
        -------
        list[str]
            Labels: ``['poly_alkyl', 'poly_aryl', 'poly_cyclic',
            'g34_F', 'g34_Cl', …, 'g44_H', 'g29', 'g30', …, 'g76']``.
        """
        poly_names = ["poly_alkyl", "poly_aryl", "poly_cyclic"]
        per_names = [
            f"g{gid}_{h}"
            for gid in PER_GROUP_IDS
            for h in HALOGENS_ORDER
        ]
        gen_names = [f"g{gid}" for gid in GENERIC_GROUP_VOCAB]
        return poly_names + per_names + gen_names

    def __repr__(self) -> str:
        nz_poly = sum(1 for v in self.poly_counts.values() if v)
        nz_per = sum(1 for v in self.per_halogen_sizes.values() if v)
        nz_gen = sum(1 for v in self.generic_groups.values() if v)
        return (
            f"GroupFeatureResult("
            f"poly_counts={nz_poly} nonzero, "
            f"per_halogen_sizes={nz_per} nonzero, "
            f"generic_groups={nz_gen} nonzero)"
        )


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def extract_group_features(
    mol: Union[Chem.Mol, str],
) -> GroupFeatureResult:
    """Extract structured halogen-group features for a single molecule.

    Two :func:`~PFASGroups.parser.parse_mol` calls are made internally:

    1. ``halogens=['H','F','Cl','Br','I']`` with the six component groups
       (ids 34, 35, 37, 38, 44, 45) — populates :attr:`~GroupFeatureResult.poly_counts`
       and :attr:`~GroupFeatureResult.per_halogen_sizes`.
    2. ``halogens=['*']`` (wildcard) — populates
       :attr:`~GroupFeatureResult.generic_groups`.

    Parameters
    ----------
    mol : rdkit.Chem.Mol or str
        RDKit molecule object or a SMILES string.

    Returns
    -------
    GroupFeatureResult
        Populated result object.  Call :meth:`~GroupFeatureResult.to_array` to
        obtain a fixed-length float32 array of shape ``(66,)``.

    Raises
    ------
    ValueError
        If a SMILES string cannot be parsed by RDKit.
    RuntimeError
        If an unexpected error occurs during :func:`parse_mol`.

    Examples
    --------
    >>> from rdkit import Chem
    >>> from PFASGroups import extract_group_features
    >>> pfoa = Chem.MolFromSmiles("FC(F)(F)C(F)(F)C(F)(F)C(=O)O")
    >>> r = extract_group_features(pfoa)
    >>> r.poly_counts
    {'poly_alkyl': 1.0, 'poly_aryl': 0.0, 'poly_cyclic': 0.0}
    >>> r.per_halogen_sizes['g34_F']
    4.0
    >>> r.h_chain_sizes
    {'alkyl_H': 0.0, 'aryl_H': 0.0, 'cyclic_H': 0.0}
    >>> arr = r.to_array()
    >>> arr.shape
    (66,)
    >>> # Octane — only H alkyl chain detected
    >>> octane = Chem.MolFromSmiles("CCCCCCCC")
    >>> r2 = extract_group_features(octane)
    >>> r2.h_chain_sizes['alkyl_H']
    6.0
    """
    if isinstance(mol, str):
        mol = Chem.MolFromSmiles(mol)
        if mol is None:
            raise ValueError(f"RDKit could not parse SMILES string.")

    # --- Call 1: component groups with all halogens including H ---
    emb_chain = _parse_mol(
        mol,
        halogens=["H", "F", "Cl", "Br", "I"],
        pfas_groups=_get_component_groups(),
        compute_component_metrics=False,
    )

    # --- Call 2: wildcard generic functional groups ---
    emb_wc = _parse_mol(
        mol,
        halogens=["*"],
        compute_component_metrics=False,
    )

    # ------------------------------------------------------------------
    # Build poly_counts
    # ------------------------------------------------------------------
    poly_vals: dict[str, float] = {"poly_alkyl": 0.0, "poly_aryl": 0.0, "poly_cyclic": 0.0}
    _poly_key = {35: "poly_alkyl", 38: "poly_aryl", 45: "poly_cyclic"}
    for match in emb_chain.matches:
        if not match.is_group:
            continue
        key = _poly_key.get(match.group_id)
        if key is not None:
            poly_vals[key] += 1.0

    # ------------------------------------------------------------------
    # Build per_halogen_sizes (15 entries: 3 groups × 5 halogens)
    # ------------------------------------------------------------------
    phs: dict[str, float] = {
        f"g{gid}_{h}": 0.0
        for gid in PER_GROUP_IDS
        for h in HALOGENS_ORDER
    }
    for match in emb_chain.matches:
        if not match.is_group:
            continue
        row = _PER_GID_TO_ROW.get(match.group_id)
        if row is None:
            continue
        gid = PER_GROUP_IDS[row]
        h_val = match.get("halogen")
        hal_hits: list[str]
        if isinstance(h_val, list):
            hal_hits = [x for x in h_val if x in _HAL_TO_COL]
        elif h_val in _HAL_TO_COL:
            hal_hits = [h_val]
        else:
            continue
        max_size = 0.0
        for comp in match.components:
            s = comp.data.get("size", 0) or 0
            if s > max_size:
                max_size = float(s)
        for hal in hal_hits:
            k = f"g{gid}_{hal}"
            if max_size > phs[k]:
                phs[k] = max_size

    # ------------------------------------------------------------------
    # Build h_chain_sizes (convenience subset: H column of per_halogen_sizes)
    # ------------------------------------------------------------------
    h_chain: dict[str, float] = {
        "alkyl_H": phs["g34_H"],
        "aryl_H": phs["g37_H"],
        "cyclic_H": phs["g44_H"],
    }

    # ------------------------------------------------------------------
    # Build generic_groups (48 entries, ids 29-76)
    # ------------------------------------------------------------------
    gen: dict[str, float] = {f"g{gid}": 0.0 for gid in GENERIC_GROUP_VOCAB}
    for match in emb_wc.matches:
        if match.get("halogen") != "*":
            continue
        gid = match.group_id
        k = f"g{gid}"
        if k in gen:
            gen[k] += 1.0

    return GroupFeatureResult(
        poly_counts=poly_vals,
        per_halogen_sizes=phs,
        h_chain_sizes=h_chain,
        generic_groups=gen,
    )
