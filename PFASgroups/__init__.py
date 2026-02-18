"""Compatibility wrapper for legacy PFASgroups imports.

This module provides the legacy PFASgroups API while delegating implementation
to HalogenGroups. It also defaults to fluorine-only components (halogens='F')
for backwards compatibility with PFAS-focused workflows.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Dict, Iterable, List, Optional
import json

from HalogenGroups import (
    HalogenGroup,
    PFASDefinition,
    ComponentsSolver,
    ResultsModel,
    rdkit_disable_log,
    HALOGEN_GROUPS_FILE,
    parse_groups_in_mol as _hg_parse_groups_in_mol,
    parse_smiles as _hg_parse_smiles,
    parse_mols as _hg_parse_mols,
    parse_mol as _hg_parse_mol,
    parse_from_database as _hg_parse_from_database,
    setup_halogen_groups_database as _hg_setup_halogen_groups_database,
    compile_componentSmarts,
    compile_componentSmartss,
    get_componentSmartss,
    get_HalogenGroups,
    get_PFASDefinitions,
    generate_fingerprint,
    generate_homologues,
    generate_degradation_products,
    plot_mol,
    plot_mols,
    plot_HalogenGroups,
    prioritise_molecules,
    prioritize_molecules,
    get_priority_statistics,
)


__version__ = "3.0.0"


@dataclass
class _FallbackGroup:
    id: Optional[int]
    name: str

    def __str__(self) -> str:  # pragma: no cover - trivial
        return self.name


_GROUP_CACHE: Optional[Dict[Optional[int], Any]] = None


def _ensure_default_halogens(kwargs: Dict[str, Any]) -> Dict[str, Any]:
    if kwargs.get("halogens") is None:
        kwargs["halogens"] = "F"
    return kwargs


def _get_group_map(pfas_groups: Optional[Iterable[Any]] = None) -> Dict[Optional[int], Any]:
    if pfas_groups:
        return {getattr(g, "id", None): g for g in pfas_groups if getattr(g, "id", None) is not None}

    global _GROUP_CACHE
    if _GROUP_CACHE is None:
        with open(HALOGEN_GROUPS_FILE, "r", encoding="utf-8") as f:
            data = json.load(f)
        _GROUP_CACHE = {g.id: g for g in (HalogenGroup(**x) for x in data)}
    return _GROUP_CACHE


def _to_legacy(results: Iterable[dict], group_map: Dict[Optional[int], Any]) -> List[list]:
    legacy: List[list] = []
    for mol_res in results:
        matches: List[tuple] = []
        for match in mol_res.get("matches", []):
            if match.get("type") != "HalogenGroup":
                continue
            group_id = match.get("id")
            group = group_map.get(group_id)
            if group is None:
                group = _FallbackGroup(group_id, match.get("group_name") or f"Group {group_id}")
            matches.append(
                (
                    group,
                    match.get("match_count", 0),
                    match.get("components_sizes", []),
                    match.get("components", []),
                )
            )
        legacy.append(matches)
    return legacy


def parse_groups_in_mol(mol, *args, **kwargs):
    _ensure_default_halogens(kwargs)
    return _hg_parse_groups_in_mol(mol, *args, **kwargs)


def parse_mol(mol, *args, include_PFAS_definitions: bool = False, **kwargs):
    _ensure_default_halogens(kwargs)
    result = _hg_parse_mol(mol, *args, include_PFAS_definitions=include_PFAS_definitions, **kwargs)
    if include_PFAS_definitions:
        return result
    group_map = _get_group_map(kwargs.get("pfas_groups"))
    return _to_legacy([result], group_map)[0]


def parse_mols(
    mols,
    *args,
    output_format: str = "list",
    include_PFAS_definitions: bool = False,
    **kwargs,
):
    _ensure_default_halogens(kwargs)
    result = _hg_parse_mols(
        mols,
        *args,
        output_format=output_format,
        include_PFAS_definitions=include_PFAS_definitions,
        **kwargs,
    )
    if output_format != "list" or include_PFAS_definitions:
        return result
    group_map = _get_group_map(kwargs.get("pfas_groups"))
    return _to_legacy(result, group_map)


def parse_smiles(
    smiles,
    *args,
    output_format: str = "list",
    include_PFAS_definitions: bool = False,
    **kwargs,
):
    _ensure_default_halogens(kwargs)
    result = _hg_parse_smiles(
        smiles,
        *args,
        output_format=output_format,
        include_PFAS_definitions=include_PFAS_definitions,
        **kwargs,
    )
    if output_format != "list" or include_PFAS_definitions:
        return result
    group_map = _get_group_map(kwargs.get("pfas_groups"))
    return _to_legacy(result, group_map)


def parse_from_database(*args, **kwargs):
    _ensure_default_halogens(kwargs)
    return _hg_parse_from_database(*args, **kwargs)


def setup_halogen_groups_database(*args, **kwargs):
    return _hg_setup_halogen_groups_database(*args, **kwargs)


__all__ = [
    "HalogenGroup",
    "PFASDefinition",
    "ComponentsSolver",
    "ResultsModel",
    "rdkit_disable_log",
    "HALOGEN_GROUPS_FILE",
    "parse_groups_in_mol",
    "parse_smiles",
    "parse_mols",
    "parse_mol",
    "parse_from_database",
    "setup_halogen_groups_database",
    "compile_componentSmarts",
    "compile_componentSmartss",
    "get_componentSmartss",
    "get_HalogenGroups",
    "get_PFASDefinitions",
    "generate_fingerprint",
    "generate_homologues",
    "generate_degradation_products",
    "plot_mol",
    "plot_mols",
    "plot_HalogenGroups",
    "prioritise_molecules",
    "prioritize_molecules",
    "get_priority_statistics",
]
