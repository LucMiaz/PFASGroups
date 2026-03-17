"""HalogenGroups - Multi-halogen analysis package.

This package extends PFASGroups to support all halogens (F, Cl, Br, I) by
default. When imported as ``HalogenGroups``, all public functions default to
``halogens=['F', 'Cl', 'Br', 'I']``, whereas ``PFASGroups`` defaults to
``halogens='F'`` (fluorine only).

Examples
--------
>>> from HalogenGroups import parse_smiles          # defaults to all halogens
>>> from PFASGroups import parse_smiles             # defaults to F only
"""
import functools
from typing import Union, List, Optional, Any, Iterable, Dict

# ---------------------------------------------------------------------------
# Re-export everything from PFASGroups that does not need default overriding.
# ---------------------------------------------------------------------------
from PFASGroups import (
    HalogenGroup,
    PFASDefinition,
    ComponentsSolver,
    rdkit_disable_log,
    HALOGEN_GROUPS_FILE,
    parse_mol,
    parse_groups_in_mol,
    parse_from_database,
    setup_halogen_groups_database,
    compile_componentSmarts,
    compile_componentSmartss,
    load_HalogenGroups,
    plot_mol,
    plot_mols,
    plot_HalogenGroups,
    get_componentSMARTSs,
    get_HalogenGroups,
    get_compiled_HalogenGroups,
    get_compiled_PFASGroups,
    get_PFASDefinitions,
    generate_homologues,
    HomologueSeries,
    HomologueEntry,
    generate_degradation_products,
    MoleculeResult,
    PFASEmbedding,
    prioritise_molecules,
    prioritize_molecules,
    get_priority_statistics,
    FINGERPRINT_PRESETS,
    EMBEDDING_PRESETS,
)
from PFASGroups import (
    parse_smiles as _parse_smiles_base,
    parse_mols as _parse_mols_base,
    PFASEmbeddingSet as _PFASEmbeddingSet,
    generate_fingerprint as _generate_fingerprint_base,
)

__version__ = "3.1.0"

_ALL_HALOGENS = ['F', 'Cl', 'Br', 'I']


# ---------------------------------------------------------------------------
# PFASEmbeddingSet subclass — kept for backward compatibility.
# halogens are already baked in at parse time; to_array() does not re-parse.
# ---------------------------------------------------------------------------

class PFASEmbeddingSet(_PFASEmbeddingSet):
    """PFASEmbeddingSet with all-halogens default in to_fingerprint().

    Identical to ``PFASGroups.PFASEmbeddingSet`` except that the deprecated
    :meth:`to_fingerprint` records that all four halogens were used.
    """

    def to_fingerprint(
        self,
        *,
        group_selection: str = 'all',
        component_metrics: Optional[List[str]] = None,
        selected_group_ids: Optional[List[int]] = None,
        halogens: Optional[Union[str, List[str]]] = None,
        saturation: Optional[str] = 'per',
        molecule_metrics: Optional[List[str]] = None,
        pfas_groups: Optional[List[Dict]] = None,
        preset: Optional[str] = None,
        count_mode: Optional[str] = None,
        graph_metrics: Optional[List[str]] = None,
        **kwargs,
    ):
        """Deprecated — use :meth:`to_array` instead.

        ``halogens`` is accepted for backward compatibility but has no effect:
        the halogen filter is applied at parse time (in :func:`parse_smiles`),
        not at embedding time.
        """
        return super().to_fingerprint(
            group_selection=group_selection,
            component_metrics=component_metrics,
            selected_group_ids=selected_group_ids,
            halogens=halogens or _ALL_HALOGENS,
            saturation=saturation,
            molecule_metrics=molecule_metrics,
            pfas_groups=pfas_groups,
            preset=preset,
            count_mode=count_mode,
            graph_metrics=graph_metrics,
            **kwargs,
        )


# Backward-compatible alias
ResultsModel = PFASEmbeddingSet


# ---------------------------------------------------------------------------
# Wrapper functions with all-halogens defaults
# ---------------------------------------------------------------------------

def parse_smiles(smiles, *, bycomponent=False, output_format='list',
                 limit_effective_graph_resistance=None,
                 compute_component_metrics=True,
                 halogens=None, form=None, saturation=None, **kwargs):
    """Parse SMILES string(s), defaulting to all halogens.

    Identical to ``PFASGroups.parse_smiles`` except the default for
    ``halogens`` is ``None`` (no filter = all halogens: F, Cl, Br, I)
    instead of ``'F'``.

    Parameters
    ----------
    smiles : str or list of str
    halogens : str or list of str or None, default None
        ``None`` → all halogens (F, Cl, Br, I).
        Pass ``'F'`` to restrict to fluorine only.

    See Also
    --------
    PFASGroups.parse_smiles : fluorine-only default version.
    """
    result = _parse_smiles_base(
        smiles,
        bycomponent=bycomponent,
        output_format=output_format,
        limit_effective_graph_resistance=limit_effective_graph_resistance,
        compute_component_metrics=compute_component_metrics,
        halogens=halogens,
        form=form,
        saturation=saturation,
        **kwargs,
    )
    # Wrap result in HalogenGroups PFASEmbeddingSet so that to_fingerprint()
    # also defaults to all halogens.
    if isinstance(result, _PFASEmbeddingSet):
        return PFASEmbeddingSet(result)
    return result


def parse_mols(mols, *, output_format='list', include_PFAS_definitions=True,
               limit_effective_graph_resistance=None,
               compute_component_metrics=True,
               halogens=None, form=None, saturation=None, **kwargs):
    """Parse RDKit molecule(s), defaulting to all halogens.

    Identical to ``PFASGroups.parse_mols`` except the default for
    ``halogens`` is ``None`` (no filter = all halogens: F, Cl, Br, I)
    instead of ``'F'``.

    Parameters
    ----------
    mols : list of rdkit.Chem.Mol
    halogens : str or list of str or None, default None
        ``None`` → all halogens (F, Cl, Br, I).

    See Also
    --------
    PFASGroups.parse_mols : fluorine-only default version.
    """
    result = _parse_mols_base(
        mols,
        output_format=output_format,
        include_PFAS_definitions=include_PFAS_definitions,
        limit_effective_graph_resistance=limit_effective_graph_resistance,
        compute_component_metrics=compute_component_metrics,
        halogens=halogens,
        form=form,
        saturation=saturation,
        **kwargs,
    )
    if isinstance(result, _PFASEmbeddingSet):
        return PFASEmbeddingSet(result)
    return result


def generate_fingerprint(smiles, *, selected_groups=None, representation='vector',
                         component_metrics=None,
                         halogens: Optional[Union[str, List[str]]] = None,
                         saturation: Optional[str] = 'per',
                         count_mode=None,
                         **kwargs):
    """Generate halogen-group fingerprints, defaulting to all halogens.

    Identical to ``PFASGroups.generate_fingerprint`` except the default for
    ``halogens`` is ``['F', 'Cl', 'Br', 'I']`` (stacked vector) instead of
    ``'F'``.

    Parameters
    ----------
    smiles : str or list of str
    halogens : str or list of str, default ['F', 'Cl', 'Br', 'I']
        Halogens to include.

    See Also
    --------
    PFASGroups.generate_fingerprint : fluorine-only default version.
    """
    if halogens is None:
        halogens = _ALL_HALOGENS
    return _generate_fingerprint_base(
        smiles,
        selected_groups=selected_groups,
        representation=representation,
        component_metrics=component_metrics,
        halogens=halogens,
        saturation=saturation,
        count_mode=count_mode,
        **kwargs,
    )


# ---------------------------------------------------------------------------
# __all__
# ---------------------------------------------------------------------------
__all__ = [
    'HalogenGroup', 'PFASDefinition', 'ComponentsSolver',
    'rdkit_disable_log', 'HALOGEN_GROUPS_FILE',
    'parse_smiles', 'parse_mols', 'parse_mol', 'parse_groups_in_mol',
    'parse_from_database', 'setup_halogen_groups_database',
    'compile_componentSmarts', 'compile_componentSmartss', 'load_HalogenGroups',
    'plot_HalogenGroups', 'plot_mol', 'plot_mols',
    'get_componentSMARTSs', 'get_HalogenGroups', 'get_compiled_HalogenGroups',
    'get_compiled_PFASGroups', 'get_PFASDefinitions',
    'generate_fingerprint', 
    'FINGERPRINT_PRESETS', 'EMBEDDING_PRESETS',
    'generate_homologues', 'HomologueSeries', 'HomologueEntry',
    'generate_degradation_products',
    'PFASEmbedding', 'PFASEmbeddingSet', 'ResultsModel', 'MoleculeResult',
    'prioritise_molecules', 'prioritize_molecules', 'get_priority_statistics',
]
