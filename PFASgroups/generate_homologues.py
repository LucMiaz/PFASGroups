from .core import add_componentSmarts, get_substruct, remove_atoms
from .homologue_series import HomologueSeries
from rdkit import Chem
from rdkit.Chem.rdMolDescriptors import CalcMolFormula
import networkx as nx
from itertools import combinations


def _mol_from_input(mol_or_string):
    """Accept a Mol, SMILES string, or InChI string; always return a Chem.Mol.

    Parameters
    ----------
    mol_or_string :
        ``rdkit.Chem.Mol``, a SMILES string, or an InChI string
        (must start with ``'InChI='``).

    Returns
    -------
    rdkit.Chem.Mol

    Raises
    ------
    ValueError
        If the string cannot be parsed as a valid molecule.
    TypeError
        If the input type is not supported.
    """
    if isinstance(mol_or_string, Chem.Mol):
        return mol_or_string
    if isinstance(mol_or_string, str):
        s = mol_or_string.strip()
        if s.startswith('InChI='):
            mol = Chem.MolFromInchi(s)
        else:
            mol = Chem.MolFromSmiles(s)
        if mol is None:
            raise ValueError(
                f"Could not parse input as a molecule: {s!r}"
            )
        return mol
    raise TypeError(
        f"Expected rdkit.Chem.Mol or str, got {type(mol_or_string).__name__}"
    )

# Halogen element symbols and their atomic numbers used to build per-halogen SMARTS
_HALOGEN_ATOMIC_NUM = {'F': 9, 'Cl': 17, 'Br': 35, 'I': 53}

# Default componentSmartsName per halogen (per-halo, alkyl chain)
_DEFAULT_COMPONENT_NAME = {
    'F': 'Perfluoroalkyl',
    'Cl': 'Perchloroalkyl',
    'Br': 'Perbromoalkyl',
    'I': 'Periodoalkyl',
}


def find_halogenated_components(mol, component_smarts, halogen='F'):
    """Find connected halogenated subgraphs in *mol* using a component SMARTS.

    Unlike the old ``find_chain`` approach, this function does **not** require
    ``end`` SMARTS or shortest-path traversal.  Instead it:

    1. Matches every atom in *mol* against *component_smarts* (the ``component``
       entry from ``component_smarts_halogens.json``).
    2. Builds a subgraph restricted to the matched atoms.
    3. Returns the connected components of that subgraph together with the set
       of halogen-bearing carbon indices (candidate repeating-unit carbons)
       within each component.

    Parameters
    ----------
    mol : rdkit.Chem.Mol
        Input molecule.
    component_smarts : rdkit.Chem.Mol
        Pre-compiled SMARTS mol object for the halogenated component pattern.
    halogen : str, default ``'F'``
        Halogen element symbol (``'F'``, ``'Cl'``, ``'Br'``, ``'I'``).

    Returns
    -------
    list of dict
        Each entry describes one halogenated component:

        - ``'component'``  (``frozenset[int]``) – atom indices of all component atoms.
        - ``'cx2_carbons'`` (``list[int]``)      – backbone C atoms bearing ≥ 2 halogen
          substituents (i.e., candidate ``CX2`` units for homologue generation).
    """
    component_atom_idxs = get_substruct(mol, component_smarts)
    if not component_atom_idxs:
        return []

    # Build subgraph restricted to component atoms
    G = nx.Graph()
    G.add_nodes_from(component_atom_idxs)
    for bond in mol.GetBonds():
        a, b = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        if a in component_atom_idxs and b in component_atom_idxs:
            G.add_edge(a, b)

    results = []
    for cc in nx.connected_components(G):
        # Identify CX2 backbone carbons: C in component bearing *exactly* 2 halogen
        # neighbours.  Terminal CF3 (3 halogens) and mono-halo atoms are excluded —
        # only true -CX2- repeating units qualify.
        cx2 = [
            idx for idx in cc
            if mol.GetAtomWithIdx(idx).GetSymbol() == 'C'
            and sum(
                1 for nb in mol.GetAtomWithIdx(idx).GetNeighbors()
                if nb.GetSymbol() == halogen
            ) == 2
        ]
        results.append({'component': frozenset(cc), 'cx2_carbons': cx2})

    return results


@add_componentSmarts()
def generate_homologues(mol_input, componentSmartsName=None, componentSmartss=None,
                        halogen='F', base_repeating=None):
    """Generate a homologous series by removing halogenated repeating units.

    This implementation uses a **component-based** approach:

    1. The ``component`` SMARTS from *componentSmartss* is used to identify all
       atoms belonging to halogenated sub-structures (no ``end`` SMARTS needed).
    2. Within each connected halogenated component, every carbon bearing ≥ 2
       halogen atoms (a ``CX2`` unit) is a candidate for removal.
    3. All non-empty subsets of those ``CX2`` units are tried; each subset is
       removed via :func:`~HalogenGroups.core.remove_atoms`, rebuilding
       connectivity across the gap.
    4. Valid, connected products are collected and deduplicated by InChIKey.

    Parameters
    ----------
    mol : rdkit.Chem.Mol
        Parent molecule.
    componentSmartsName : str, optional
        Key in *componentSmartss* selecting which component pattern to use.
        Defaults to the per-halogen value in ``_DEFAULT_COMPONENT_NAME``
        (e.g. ``'Perfluoroalkyl'`` for ``halogen='F'``).
    componentSmartss : dict, optional
        Mapping ``{name: {'component': <Mol>, ...}}``.
        Populated automatically by the ``@add_componentSmarts`` decorator.
    halogen : str, default ``'F'``
        Halogen element symbol.  Controls which atoms count as the repeating
        substituents and which component pattern is selected by default.
        Accepted values: ``'F'``, ``'Cl'``, ``'Br'``, ``'I'``.
    base_repeating : list of str, optional
        Backbone element symbols that are *kept* when a CX2 unit is removed
        (their halogen substituents are stripped, not the backbone C itself).
        Defaults to ``['C']``.

    Returns
    -------
    dict
        ``{InChIKey: {formula: rdkit.Chem.Mol}}`` – all unique shorter
        homologues, keyed first by InChIKey and then by molecular formula.

    Raises
    ------
    ValueError
        If *halogen* is not one of ``'F'``, ``'Cl'``, ``'Br'``, ``'I'``.
    ValueError
        If atom removal yields a fragmented molecule (diagnostic aid).

    Examples
    --------
    >>> from rdkit import Chem
    >>> from HalogenGroups.generate_homologues import generate_homologues
    >>> # PFOA – generate all shorter perfluoroalkyl chain homologues
    >>> pfoa = Chem.MolFromSmiles('OC(=O)' + 'C(F)(F)' * 7 + 'F')
    >>> homologues = generate_homologues(pfoa)
    >>> print(len(homologues))  # 6  (C2-C7)

    >>> # Chlorinated analogue
    >>> pca = Chem.MolFromSmiles('OC(=O)' + 'C(Cl)(Cl)' * 5 + 'Cl')
    >>> homologues = generate_homologues(pca, halogen='Cl')

    >>> # From SMILES string directly
    >>> homologues = generate_homologues('OC(=O)C(F)(F)C(F)(F)C(F)(F)F')

    >>> # From InChI string
    >>> homologues = generate_homologues('InChI=1S/C3HF5O2/...')

    Notes
    -----
    * Removing a CX2 unit strips the backbone carbon **and** all directly
      attached halogen atoms; connectivity to the rest of the molecule is
      preserved by :func:`~HalogenGroups.core.remove_atoms`.
    * Branched structures are handled correctly because component detection
      is graph-based, not path-based.
    * The original molecule is *not* included in the output.
    """
    # Accept SMILES / InChI strings in addition to Chem.Mol
    mol = _mol_from_input(mol_input)

    if halogen not in _HALOGEN_ATOMIC_NUM:
        raise ValueError(
            f"halogen must be one of {list(_HALOGEN_ATOMIC_NUM)}, got {halogen!r}"
        )
    if base_repeating is None:
        base_repeating = ['C']

    # Resolve component SMARTS name
    if componentSmartsName is None:
        componentSmartsName = _DEFAULT_COMPONENT_NAME[halogen]

    # Look up the component SMARTS mol object
    entry = componentSmartss.get(componentSmartsName)
    if entry is None:
        raise ValueError(
            f"Component '{componentSmartsName}' not found in componentSmartss. "
            f"Available: {list(componentSmartss)}"
        )
    component_smarts = entry['component'] if isinstance(entry, dict) else entry

    # Atoms that are removed together with a CX2 backbone carbon
    # (only the halogen neighbours are stripped alongside the backbone C)

    # Discover halogenated components and their CX2 carbons
    components = find_halogenated_components(mol, component_smarts, halogen=halogen)

    # Collect all CX2 carbon indices across all components
    all_cx2 = [idx for comp in components for idx in comp['cx2_carbons']]

    series = HomologueSeries()

    if not all_cx2:
        series._set_metadata(mol, halogen, componentSmartsName)
        return series

    # Enumerate every non-empty proper subset of CX2 units
    for r in range(1, len(all_cx2) + 1):
        for subset in combinations(all_cx2, r):
            flat_idx = list(subset)
            try:
                h = remove_atoms(mol, flat_idx, removable=[halogen])
            except Exception as e:
                print(
                    f"Error removing atoms {flat_idx} from "
                    f"{Chem.MolToSmiles(mol)}: {e}"
                )
                raise e
            if len(Chem.GetMolFrags(h)) > 1:
                raise ValueError(
                    f"Fragmented molecule {Chem.MolToSmiles(h)} after removing "
                    f"atoms {flat_idx} from {Chem.MolToSmiles(mol)}"
                )
            inchikey = Chem.MolToInchiKey(h)
            formula = CalcMolFormula(h)
            series.setdefault(inchikey, {})[formula] = h

    series._set_metadata(mol, halogen, componentSmartsName)
    return series
