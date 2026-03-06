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
from typing import Union, List, Optional, Any, Iterable

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
    get_componentSmartss,
    get_HalogenGroups,
    get_compiled_HalogenGroups,
    get_compiled_PFASGroups,
    get_PFASDefinitions,
    generate_homologues,
    HomologueSeries,
    HomologueEntry,
    generate_degradation_products,
    MoleculeResult,
    ResultsFingerprint,
    prioritise_molecules,
    prioritize_molecules,
    get_priority_statistics,
)
from PFASGroups import (
    parse_smiles as _parse_smiles_base,
    parse_mols as _parse_mols_base,
    generate_fingerprint as _generate_fingerprint_base,
    ResultsModel as _ResultsModel,
)

__version__ = "3.1.0"

_ALL_HALOGENS = ['F', 'Cl', 'Br', 'I']


# ---------------------------------------------------------------------------
# ResultsModel subclass with all-halogens default in to_fingerprint()
# ---------------------------------------------------------------------------

class ResultsModel(_ResultsModel):
    """ResultsModel with all-halogens default for fingerprint generation.

    Identical to ``PFASGroups.ResultsModel`` except that
    :meth:`to_fingerprint` defaults to ``halogens=['F', 'Cl', 'Br', 'I']``
    instead of ``halogens='F'``.
    """

    def to_fingerprint(
        self,
        group_selection: str = 'all',
        count_mode: str = 'binary',
        selected_group_ids: Optional[List[int]] = None,
        halogens: Union[str, List[str]] = _ALL_HALOGENS,
        saturation: Optional[str] = 'per',
    ) -> 'ResultsFingerprint':
        """Convert ResultsModel to ResultsFingerprint for dimensionality reduction.

        Parameters
        ----------
        group_selection : str, default 'all'
            Which groups to include in fingerprint:
            - 'all': All available groups
            - 'oecd': OECD groups (IDs 1–28)
            - 'generic': Generic functional-group groups (IDs 29–55)
            - 'telomers': Telomer-related groups (IDs 74–116)
            - 'generic+telomers': Combination of generic and telomer groups
        count_mode : str, default 'binary'
            How to encode group matches:
            - 'binary': 1 if present, 0 if absent
            - 'count': Number of matches
            - 'max_component': Maximum component size
        selected_group_ids : list of int, optional
            Explicit list of group IDs to include (overrides group_selection)
        halogens : str or list of str, default ['F', 'Cl', 'Br', 'I']
            Which halogen(s) to match component SMARTS against.
            Multiple values produce one fingerprint per halogen,
            concatenated into a vector of length n_groups × n_halogens.
            Group names are suffixed with ``[F]``, ``[Cl]``, etc.
            Available: ``'F'``, ``'Cl'``, ``'Br'``, ``'I'``
        saturation : str or None, default 'per'
            Saturation filter applied to component SMARTS groups:
            - ``'per'``: perhalogenated only
            - ``'poly'``: polyhalogenated only
            - ``None``: no filter
        """
        return super().to_fingerprint(
            group_selection=group_selection,
            count_mode=count_mode,
            selected_group_ids=selected_group_ids,
            halogens=halogens,
            saturation=saturation,
        )


# ---------------------------------------------------------------------------
# Wrapper functions with all-halogens defaults
# ---------------------------------------------------------------------------

def parse_smiles(smiles, bycomponent=False, output_format='list',
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
    # Wrap result in HalogenGroups ResultsModel so that to_fingerprint()
    # also defaults to all halogens.
    if isinstance(result, _ResultsModel):
        return ResultsModel(result)
    return result


def parse_mols(mols, output_format='list', include_PFAS_definitions=True,
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
    if isinstance(result, _ResultsModel):
        return ResultsModel(result)
    return result


def generate_fingerprint(smiles, selected_groups=None, representation='vector',
                         count_mode='binary',
                         halogens: Union[str, List[str]] = _ALL_HALOGENS,
                         saturation: Optional[str] = 'per',
                         **kwargs):
    """Generate halogen-group fingerprints, defaulting to all halogens.

    Identical to ``PFASGroups.generate_fingerprint`` except the default for
    ``halogens`` is ``['F', 'Cl', 'Br', 'I']`` (stacked vector) instead of
    ``'F'``.

    Parameters
    ----------
    smiles : str or list of str
    halogens : str or list of str, default ['F', 'Cl', 'Br', 'I']
        Halogens to include. Multiple values produce a stacked fingerprint
        of length n_groups × n_halogens.

    See Also
    --------
    PFASGroups.generate_fingerprint : fluorine-only default version.
    """
    return _generate_fingerprint_base(
        smiles,
        selected_groups=selected_groups,
        representation=representation,
        count_mode=count_mode,
        halogens=halogens,
        saturation=saturation,
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
    'get_componentSmartss', 'get_HalogenGroups', 'get_compiled_HalogenGroups', 'get_compiled_PFASGroups', 'get_PFASDefinitions',
    'generate_fingerprint',
    'generate_homologues', 'HomologueSeries', 'HomologueEntry',
    'generate_degradation_products',
    'ResultsModel', 'MoleculeResult', 'ResultsFingerprint',
    'prioritise_molecules', 'prioritize_molecules', 'get_priority_statistics',
]
