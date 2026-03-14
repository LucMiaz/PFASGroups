"""
HalogenGroups core functions.

This module provides functions for parsing and plotting halogen groups.
"""

import json
import functools
from typing import Optional

from rdkit import Chem
from rdkit.Chem.rdMolDescriptors import CalcMolFormula

from .HalogenGroupModel import HalogenGroup
from .PFASDefinitionModel import PFASDefinition
from .ComponentsSolverModel import ComponentsSolver
from .core import (
    fragment_until_valence_is_correct,
    n_from_formula,
    add_componentSmarts,
    PFAS_DEFINITIONS_FILE,
    HALOGEN_GROUPS_FILE,
    rdkit_disable_log,
)



# --- Load halogen groups from PFAS_groups_smarts.json ---
def load_HalogenGroups():
    """
    Adds default HalogenGroups to function
    """
    with open(HALOGEN_GROUPS_FILE,'r') as f:
        _pfg = json.load(f)
    pfg = [HalogenGroup(**x) for x in _pfg if x.get('compute',True)]
    agg_pfg = [HalogenGroup(**x) for x in _pfg if not x.get('compute',True)]
    # list of HalogenGroup names
    pfg_names =  {pf.name:pf.id for pf in pfg}
    # list groups aggregated by groups with compute=FALSE
    agg_pfg = {ppf:list(map(pfg_names.get,list(filter(ppf.re_search.search,pfg_names.keys())))) for ppf in agg_pfg}
    def inner(func):
        @functools.wraps(func)
        def wrapper(*args,**kwargs):
            kwargs['pfas_groups'] = kwargs.get('pfas_groups',[p for p in pfg if p.excludeHalogens is None or set(p.excludeHalogens).isdisjoint(kwargs.get('halogens', ['F','Cl','Br','I']))])
            kwargs['agg_pfas_groups'] = kwargs.get('agg_pfas_groups',agg_pfg)
            return func(*args, **kwargs)
        return wrapper
    return inner

def load_PFASDefinitions():
    """
    Adds default PFAS definitions to function
    """
    with open(PFAS_DEFINITIONS_FILE,'r') as f:
        pfg = json.load(f)
    pfg = [PFASDefinition(**x) for x in pfg]
    def inner(func):
        @functools.wraps(func)
        def wrapper(*args,**kwargs):
            kwargs['pfas_definitions'] = kwargs.get('pfas_definitions',pfg)
            return func(*args, **kwargs)
        return wrapper
    return inner


def load_componentsSolver(**kwargs):
    """
    Adds componentsSolver to function (creates it per call with the molecule)
    """
    def inner(func):
        @functools.wraps(func)
        def wrapper(*args,**kwargs):
            mol = args[0]
            # Add hydrogens to molecule before creating ComponentsSolver
            # This ensures atom indices are consistent throughout the analysis
            mol = Chem.AddHs(mol)
            args = list(args)  # Convert to mutable list
            args[0] = mol
            args = tuple(args)  # Convert back to tuple
            halogens = kwargs.get('halogens', ['F','Cl','Br','I'])
            _smarts = '[{}]'.format(','.join(halogens)) if isinstance(halogens, list) else f'[{halogens}]'
            # check organic halogen:
            if not mol.GetSubstructMatch(Chem.MolFromSmarts(_smarts)):
                #logger.debug("No organic halogens found, skipping componentsSolver")
                return [], mol
            # Pass through component metric options
            solver_kwargs = {}
            if 'limit_effective_graph_resistance' in kwargs:
                solver_kwargs['limit_effective_graph_resistance'] = kwargs['limit_effective_graph_resistance']
            if 'compute_component_metrics' in kwargs:
                solver_kwargs['compute_component_metrics'] = kwargs['compute_component_metrics']
            if 'resistance_weights' in kwargs:
                solver_kwargs['resistance_weights'] = kwargs['resistance_weights']
            if 'component_expansion' in kwargs:
                solver_kwargs['component_expansion'] = kwargs['component_expansion']
            with ComponentsSolver(mol, halogens=kwargs.get('halogens'), **solver_kwargs) as fluorinated_components_dict:
                kwargs['fluorinated_components_dict'] = kwargs.get('fluorinated_components_dict',fluorinated_components_dict)
                return func(*args, **kwargs)
        return wrapper
    return inner


@load_PFASDefinitions()
def parse_definitions_in_mol(mol, **kwargs):
    mol = Chem.AddHs(mol)
    formula = kwargs.get("formula", CalcMolFormula(mol))
    try:
        Chem.SanitizeMol(mol)
    except Chem.AtomValenceException:
        #logger.debug("failed sanitisation, fragmenting")
        frags = fragment_until_valence_is_correct(mol, [])
    else:
        frags = [mol]
    definition_matches = []
    pfas_definitions = kwargs.get('pfas_definitions')
    for pdef in pfas_definitions:
        matched = False
        for mol in frags:
            if pdef.applies_to_molecule(mol_or_smiles=mol, **kwargs) is True:
                matched = True
                break
        if matched is True:
            definition_matches.append(pdef)
    return definition_matches


# --- Main halogen group parsing functions ---
@add_componentSmarts()
@load_HalogenGroups()
@load_componentsSolver()
@rdkit_disable_log(level='warning')
def parse_groups_in_mol(mol, fluorinated_components_dict=None, pfas_groups = None, **kwargs):
    """Iterates over halogen groups and finds the ones that match the molecule.

    Parameters
    ----------
    mol : rdkit.Chem.Mol
        RDKit molecule object
    bycomponent : bool, optional
        Whether to look for fluorinated components or for chains between functional groups
    **kwargs : dict
        Additional parameters (formula, pfas_groups, componentSmartss, etc.)

    Returns
    -------
    list of tuples
        List of (HalogenGroup, match_count, component_sizes, matched_components) tuples where:
        - HalogenGroup: The matched halogen group object
        - match_count: Number of times this group pattern was matched (int)
        - component_sizes: List of carbon component sizes found (list of int)
        - matched_components: Detailed information about matched components (list of dicts)
            Each dict contains:
            - 'component': list of atom indices
            - 'size': number of atoms in component
            - 'SMARTS': component type (e.g., 'Perfluoroalkyl', 'alkyl', 'cyclic')
            - 'branching': float [0-1], 1.0 = linear, 0.0 = highly branched
            - 'smarts_centrality': float [0-1], 1.0 = functional group at center, 0.0 = at periphery
            - 'mean_eccentricity': float, mean graph eccentricity across nodes in component
            - 'median_eccentricity': float, median graph eccentricity across nodes in component

    Notes
    -----
    For HalogenGroups
    1. with componentSmarts = 'cyclic', search for connected component matching first smarts
    2. with multiple smarts patterns or counts > 1: Find components containing all required SMARTS matches with minimum counts
    3. with single smarts pattern: Search for connected components of fluorinated atoms (for each pathType) where smarts match is in the component
    4. with smarts defined and formula constraints: search for substructure matches of smarts, and for given componentSmarts for the HalogenGroup, search for connected components of fluorinated atoms where smarts match is in the component

    """
    mol = Chem.AddHs(mol)
    formula = kwargs.get("formula", CalcMolFormula(mol))
    try:
        Chem.SanitizeMol(mol)
    except Chem.AtomValenceException:
        #logger.debug("failed sanitisation, fragmenting")
        frags = fragment_until_valence_is_correct(mol, [])
        formulas = [n_from_formula(CalcMolFormula(frag)) for frag in frags]
    else:
        frags = [mol]
        formulas = [n_from_formula(formula)]# formula as a dictionary
    agg_pfas_groups = kwargs.get('agg_pfas_groups',{})
    # halogen groups
    group_matches = []
    # map of group_id -> list of matches for quick lookup
    group_id_to_matches = {}
    for pf in pfas_groups:
        #logger.debug(f"{pf.name}")
        matched1_len = 0
        for fd,mol in zip(formulas,frags):
            # match = pfgroup_obj, match_count, component_sizes, matched_components
            # matched_components = list(dicts) with entries: 'component','size','component_fraction','smarts_matches','smarts_extra_atoms','SMARTS','branching','smarts_centrality','diameter','radius','effective_graph_resistance','eccentricity_values','mean_eccentricity','median_eccentricity','center','periphery','barycenter','min_dist_to_barycenter','min_resistance_dist_to_barycenter','min_dist_to_center','min_resistance_dist_to_center','max_dist_to_periphery','max_resistance_dist_to_periphery'
            match = pf.find_components(mol, fd, fluorinated_components_dict, **kwargs)
            if match is not None and len(match)>0:
                group_matches.extend(match)
                group_id_to_matches.setdefault(pf.id, []).extend(match)

    # Process aggregate halogen groups efficiently
    if agg_pfas_groups:

        # For each aggregate group, collect and deduplicate components
        for agg_group, component_group_ids in agg_pfas_groups.items():
            # Check if any component groups were matched
            matched_component_ids = [gid for gid in component_group_ids if gid in group_id_to_matches]

            if matched_component_ids:
                # Collect all components from matched groups
                all_components = []
                for gid in matched_component_ids:
                    for _, match_count, component_sizes, matched_components in group_id_to_matches[gid]:
                        all_components.extend(matched_components)

                # Deduplicate components by atom set while keeping different SMARTS types
                # Filter by componentSmarts if aggregate has one
                unique_components = []
                seen_keys = set()  # Track (atom_set, SMARTS_type) combinations

                for comp in all_components:
                    atoms_key = frozenset(comp.get('component', []))
                    smarts_type = comp.get('SMARTS')

                    # Filter by componentSmarts if aggregate has one (not None)
                    if agg_group.componentSmarts is not None:
                        allowed = agg_group.componentSmarts
                        if isinstance(allowed, (list, tuple, set)):
                            if smarts_type not in allowed:
                                continue
                        elif smarts_type != allowed:
                            continue

                    # Deduplicate by (atoms, SMARTS_type) key
                    key = (atoms_key, smarts_type)
                    if key not in seen_keys:
                        seen_keys.add(key)
                        unique_components.append(comp)

                # Create match entry for aggregate group if we have unique components
                if unique_components:
                    component_sizes = [comp.get('size', 0) for comp in unique_components]
                    match_count = len(unique_components)
                    group_matches.append((agg_group, match_count, component_sizes, unique_components))

    return group_matches, mol


@rdkit_disable_log(level='warning')
def parse_smiles(smiles, bycomponent=False, output_format='list',
                  limit_effective_graph_resistance=None, compute_component_metrics=True,
                  halogens='F', form=None, saturation=None, progress=False,
                  **kwargs):
    """
    Parse SMILES string(s) and return halogen group information.

    Parameters:
    -----------
    smiles : str or list of str
        Single SMILES string or list of SMILES strings
    bycomponent : bool
        Whether to use component-based analysis
    output_format : str, default 'list'
        Output format: 'list' (default), 'dataframe', or 'csv'
        - 'list': Returns nested lists of tuples (default behavior)
        - 'dataframe': Returns pandas DataFrame with one row per match
        - 'csv': Returns CSV string
    limit_effective_graph_resistance : int or None, default None
        Maximum component size for computing effective graph resistance.
        - None: Compute for all components (default, may be slow for large molecules)
        - int > 0: Only compute for components with fewer atoms than this limit
        - 0: Skip computation for all components (set to NaN)
    compute_component_metrics : bool, default True
        Whether to compute graph metrics (diameter, radius, etc.) for components.
        - True: Compute all metrics (default)
        - False: Only compute component size, skip all other metrics
    halogens : str or list of str or None, default 'F'
        Filter components by halogen element.
        - 'F' (default): fluorine only
        - str (e.g. 'Cl'): restrict to that single halogen
        - list (e.g. ['F', 'Cl']): restrict to those halogens
        - None: no filter (include all halogens)
    form : str or list of str, optional
        Filter components by form type (e.g., 'alkyl', ['alkyl', 'cyclic'], or None for all)
    saturation : str or list of str, optional
        Filter components by saturation (e.g., 'per', 'poly', or None for all)
    progress : bool, default False
        If True, display a tqdm progress bar during parsing.
    **kwargs : dict
        Additional parameters (pfas_groups, componentSmartss, etc.)

    Returns:
    --------
    list, pandas.DataFrame, or str
        Depends on output_format parameter
    """
    # Convert single input to list for uniform processing
    single_input = isinstance(smiles, str)
    smiles_list = [smiles] if single_input else smiles
    mol_list = [Chem.MolFromSmiles(smi) for smi in smiles_list]
    # Parse all molecules
    kwargs['limit_effective_graph_resistance'] = limit_effective_graph_resistance
    kwargs['compute_component_metrics'] = compute_component_metrics
    kwargs['halogens'] = halogens
    kwargs['form'] = form
    kwargs['saturation'] = saturation
    return parse_mols(mol_list, bycomponent=bycomponent, output_format=output_format, progress=progress, **kwargs)

from .results_model import PFASEmbeddingSet


def parse_from_database(
    conn,
    query: Optional[str] = None,
    table: Optional[str] = None,
    mol_column: str = 'mol',
    smiles_column: Optional[str] = 'smiles',
    inchi_column: Optional[str] = 'inchi',
    id_column: Optional[str] = 'id',
    batch_size: int = 1000,
    output_table: Optional[str] = None,
    components_table: str = "components",
    groups_table: str = "pfas_groups_in_compound",
    write_results: bool = True,
    progress: bool = False,
    **kwargs
):
    """Parse halogen groups from molecules stored in a database.

    This function reads molecules from a database table, parses them for halogen groups,
    and optionally writes the results back to the database.

    Parameters
    ----------
    conn : str or sqlalchemy.engine.Engine
        Database connection. Can be:
        - SQLAlchemy Engine object
        - Connection string (e.g., 'postgresql://user:pass@host:port/db')
    query : str, optional
        SQL query to select molecules. If not provided, selects all from table.
    table : str, optional
        Table name to query (used if query is not provided).
    mol_column : str, default 'mol'
        Column name containing RDKit mol objects (binary format).
    smiles_column : str, optional, default 'smiles'
        Fallback column with SMILES strings if mol parsing fails.
    inchi_column : str, optional, default 'inchi'
        Fallback column with InChI strings if both mol and SMILES fail.
    id_column : str, optional, default 'id'
        Column name for unique identifier (for tracking results).
    batch_size : int, default 1000
        Number of molecules to process in each batch.
    output_table : str, optional
        If provided, creates/updates this table with summary results.
    components_table : str, default 'components'
        Table name for detailed component data.
    groups_table : str, default 'pfas_groups_in_compound'
        Table name for halogen group matches.
    write_results : bool, default True
        Whether to write results back to database.
    progress : bool, default False
        If True, display a tqdm progress bar during parsing.
    **kwargs : dict
        Additional parameters passed to parse_mols (e.g., pfas_groups, componentSmartss).

    Returns
    -------
    PFASEmbeddingSet
        Parsed results for all molecules.

    Examples
    --------
    >>> # Using connection string
    >>> results = parse_from_database(
    ...     conn='postgresql://user:pass@localhost/chem_db',
    ...     table='molecules',
    ...     mol_column='rdkit_mol',
    ...     smiles_column='canonical_smiles'
    ... )
    >>> \n>>> # Using SQLAlchemy engine with custom query
    >>> from sqlalchemy import create_engine
    >>> engine = create_engine('sqlite:///chemicals.db')
    >>> results = parse_from_database(
    ...     conn=engine,
    ...     query=\"SELECT id, mol, smiles FROM compounds WHERE molecular_weight < 1000\",
    ...     batch_size=500
    ... )
    >>> \n>>> # Process and write back to database
    >>> results = parse_from_database(
    ...     conn='postgresql://user:pass@localhost/pfas_db',
    ...     table='test_compounds',
    ...     write_results=True,
    ...     components_table='pfas_components',
    ...     groups_table='pfas_groups'
    ... )
    """
    try:
        import pandas as pd
        import sqlalchemy
        from sqlalchemy import text
    except ImportError as exc:
        raise ImportError("pandas and sqlalchemy are required. Install with: pip install pandas sqlalchemy") from exc
    # Create engine if conn is a string
    if isinstance(conn, str):
        engine = sqlalchemy.create_engine(conn)
    else:
        engine = conn

    # Build query if not provided
    if query is None:
        if table is None:
            raise ValueError("Either 'query' or 'table' must be provided.")
        query = f"SELECT * FROM {table}"

    # Read data in batches
    print(f"Reading molecules from database...")
    df = pd.read_sql(query, engine)
    total_rows = len(df)
    print(f"Found {total_rows} molecules to process")

    all_results = []

    # Process in batches
    for batch_start in range(0, total_rows, batch_size):
        batch_end = min(batch_start + batch_size, total_rows)
        batch_df = df.iloc[batch_start:batch_end]

        print(f"Processing batch {batch_start+1}-{batch_end} of {total_rows}...")

        mols = []
        mol_ids = []

        for idx, row in batch_df.iterrows():
            mol = None
            mol_id = row.get(id_column) if id_column and id_column in row else idx

            # Try to parse mol from binary column
            if mol_column in row and row[mol_column] is not None:
                try:
                    # Assuming mol is stored as binary or mol block
                    if isinstance(row[mol_column], bytes):
                        mol = Chem.Mol(row[mol_column])
                    elif isinstance(row[mol_column], str):
                        mol = Chem.MolFromMolBlock(row[mol_column])
                except Exception as e:
                    print(f"  Warning: Failed to parse mol for {mol_id}: {e}")

            # Fallback to SMILES
            if mol is None and smiles_column and smiles_column in row and row[smiles_column]:
                try:
                    mol = Chem.MolFromSmiles(row[smiles_column])
                except Exception as e:
                    print(f"  Warning: Failed to parse SMILES for {mol_id}: {e}")

            # Fallback to InChI
            if mol is None and inchi_column and inchi_column in row and row[inchi_column]:
                try:
                    mol = Chem.MolFromInchi(row[inchi_column])
                except Exception as e:
                    print(f"  Warning: Failed to parse InChI for {mol_id}: {e}")

            if mol is not None:
                mols.append(mol)
                mol_ids.append(mol_id)
            else:
                print(f"  Error: Could not parse molecule {mol_id}")

        # Parse this batch
        if mols:
            batch_results = parse_mols(mols, progress=progress, **kwargs)
            all_results.extend(batch_results)

    # Combine all results
    results = PFASEmbeddingSet(all_results)

    # Write results to database if requested
    if write_results and len(results) > 0:
        print("Writing results to database...")
        results.to_sql(
            conn=engine,
            components_table=components_table,
            groups_table=groups_table,
            if_exists='append'
        )
        print(f"✅ Successfully wrote {len(results)} results to database")

    return results


def setup_halogen_groups_database(
    conn,
    groups_info_table: str = 'halogen_groups_info',
    smarts_table: str = 'halogen_smarts',
    if_exists: str = 'replace'
):
    """Set up halogen groups metadata tables in a database.

    This function creates tables to store halogen group definitions and SMARTS patterns,
    similar to the load_pfas_groups function in zeropmdb.

    Parameters
    ----------
    conn : str or sqlalchemy.engine.Engine
        Database connection.
    groups_info_table : str, default 'halogen_groups_info'
        Table name for halogen group metadata.
    smarts_table : str, default 'halogen_smarts'
        Table name for SMARTS patterns used in groups.
    if_exists : str, default 'replace'
        How to behave if tables exist: 'fail', 'replace', or 'append'.

    Returns
    -------
    dict
        Statistics about loaded groups.

    Examples
    --------
    >>> # Set up in PostgreSQL
    >>> setup_halogen_groups_database(
    ...     conn='postgresql://user:pass@localhost/halogen_db',
    ...     groups_info_table='halogen_groups',
    ...     smarts_table='halogen_smarts_patterns'
    ... )
    >>>
    >>> # Set up in SQLite
    >>> from sqlalchemy import create_engine
    >>> engine = create_engine('sqlite:///halogen_database.db')
    >>> stats = setup_halogen_groups_database(engine)
    >>> print(f"Loaded {stats['total_groups']} halogen groups")
    """
    try:
        import pandas as pd
        import sqlalchemy
    except ImportError as exc:
        raise ImportError("pandas and sqlalchemy required. Install with: pip install pandas sqlalchemy") from exc
    # Create engine if conn is a string
    if isinstance(conn, str):
        engine = sqlalchemy.create_engine(conn)
    else:
        engine = conn

    # Load halogen groups from the module
    from .HalogenGroupModel import HalogenGroup
    import json

    # Load groups from JSON file
    from .core import HALOGEN_GROUPS_FILE_GROUPS_FILE
    with open(HALOGEN_GROUPS_FILE, 'r') as f:
        groups_data = json.load(f)

    # Prepare groups info data
    groups_info = []
    all_smarts = set()

    for group_data in groups_data:
        group_info = {
            'id': group_data['id'],
            'name': group_data['name'],
            'compute': group_data.get('compute', True),
            'componentSmarts': group_data.get('componentSmarts'),
            'pathType': group_data.get('pathType'),
            'constraints': json.dumps(group_data.get('_constraints', {})),
            'smarts_patterns': json.dumps(group_data.get('smarts', {})),
        }
        groups_info.append(group_info)

        # Collect unique SMARTS patterns
        if 'smarts' in group_data:
            for pattern in group_data['smarts'].keys():
                all_smarts.add(pattern)

    # Create DataFrames
    df_groups = pd.DataFrame(groups_info)
    df_smarts = pd.DataFrame([
        {'smarts': pattern, 'pattern_id': idx}
        for idx, pattern in enumerate(sorted(all_smarts))
    ])

    # Write to database
    print(f"Writing {len(df_groups)} halogen groups to table '{groups_info_table}'...")
    df_groups.to_sql(groups_info_table, engine, if_exists=if_exists, index=False)

    print(f"Writing {len(df_smarts)} SMARTS patterns to table '{smarts_table}'...")
    df_smarts.to_sql(smarts_table, engine, if_exists=if_exists, index=False)

    stats = {
        'total_groups': len(df_groups),
        'compute_groups': len(df_groups[df_groups['compute'].astype(bool)]),
        'aggregate_groups': len(df_groups[~df_groups['compute'].astype(bool)]),
        'total_smarts': len(df_smarts),
    }

    print("✅ Successfully set up halogen groups database")
    print(f"   - {stats['total_groups']} total groups ({stats['compute_groups']} compute, {stats['aggregate_groups']} aggregate)")
    print(f"   - {stats['total_smarts']} unique SMARTS patterns")

    return stats


def parse_mol(mol, progress=False, **kwargs):
    """Wrapper for parse_mols to handle single molecule input.

    Returns a single ``PFASEmbedding`` when ``output_format='list'``
    (default), preserving backwards-compatible dict behaviour while
    enabling richer navigation helpers.
    """
    return parse_mols([mol], progress=progress, **kwargs)[0]

def parse_mols(mols, output_format='list', include_PFAS_definitions=True,
               limit_effective_graph_resistance=None, compute_component_metrics=True,
               halogens='F', form=None, saturation=None, progress=False,
               **kwargs):
    """
    Parse RDKit molecule(s) and return halogen group information.

    Parameters:
    -----------
    mols : list of rdkit.Chem.Mol
        Single RDKit molecule or list of molecules
    bycomponent : bool
        Whether to use component-based analysis
    output_format : str, default 'list'
        Output format: 'list' (default), 'dataframe', or 'csv'
        - 'list': Returns nested lists of tuples (default behavior)
        - 'dataframe': Returns pandas DataFrame with one row per match
        - 'csv': Returns CSV string
    limit_effective_graph_resistance : int or None, default None
        Maximum component size for computing effective graph resistance.
        - None: Compute for all components (default, may be slow for large molecules)
        - int > 0: Only compute for components with fewer atoms than this limit
        - 0: Skip computation for all components (set to NaN)
    compute_component_metrics : bool, default True
        Whether to compute graph metrics (diameter, radius, etc.) for components.
        - True: Compute all metrics (default)
        - False: Only compute component size, skip all other metrics
    halogens : str or list of str or None, default 'F'
        Filter components by halogen element.
        - 'F' (default): fluorine only
        - str (e.g. 'Cl'): restrict to that single halogen
        - list (e.g. ['F', 'Cl']): restrict to those halogens
        - None: no filter (include all halogens)
    form : str or list of str, optional
        Filter components by form type (e.g., 'alkyl', ['alkyl', 'cyclic'], or None for all)
    saturation : str or list of str, optional
        Filter components by saturation (e.g., 'per', 'poly', or None for all)
    progress : bool, default False
        If True, display a tqdm progress bar during parsing.
    **kwargs : dict
        Additional parameters (halogen_groups, componentSmartss, etc.)

    Returns:
    --------
    list, pandas.DataFrame, or str
        Depends on output_format parameter
    """

    # Parse all molecules
    results = {}
    _iter = mols
    if progress:
        try:
            from tqdm.auto import tqdm as _tqdm
        except ImportError:
            from tqdm import tqdm as _tqdm
        _iter = _tqdm(mols, desc='parse_mols', total=len(mols))
    for mol in _iter:
        # Add hydrogens to ensure consistent atom indexing
        mol_with_h = Chem.AddHs(mol)
        formula = CalcMolFormula(mol)
        bycomponent = kwargs.pop('bycomponent', False)
        # Pass through component metric options and filters
        kwargs['limit_effective_graph_resistance'] = limit_effective_graph_resistance
        kwargs['compute_component_metrics'] = compute_component_metrics
        kwargs['halogens'] = halogens
        kwargs['form'] = form
        kwargs['saturation'] = saturation
        matches, mol_with_h = parse_groups_in_mol(mol, formula=formula, bycomponent=bycomponent, **kwargs)
        inchikey = Chem.MolToInchiKey(mol)
        inchi = Chem.MolToInchi(mol)
        smi = Chem.MolToSmiles(mol)
        # Store molblock with explicit H to preserve atom ordering for visualization
        # (MolToSmiles canonicalises atom order, breaking component atom indices)
        molblock_with_h = Chem.MolToMolBlock(mol_with_h)
        results.setdefault(inchikey,{}).update({"smiles": smi,
                        "inchikey": inchikey,
                        "inchi": inchi,
                        "formula": formula,
                        "smiles_with_h": Chem.MolToSmiles(mol_with_h),  # kept for backward compat
                        "molblock_with_h": molblock_with_h})

        # Build match results with comprehensive summary metrics
        match_results = []
        for group, match_count, components_sizes, matched_components in matches:
            # Calculate summary metrics across all components for this group
            if len(matched_components) > 0:
                # Basic metrics summaries
                mean_branching = sum([c['branching'] for c in matched_components])/len(matched_components)
                mean_smarts_centrality = sum([c['smarts_centrality'] for c in matched_components])/len(matched_components)
                mean_mean_eccentricity = sum([c['mean_eccentricity'] for c in matched_components])/len(matched_components)
                mean_median_eccentricity = sum([c['median_eccentricity'] for c in matched_components])/len(matched_components)
                mean_component_fraction = sum([c['component_fraction'] for c in matched_components])/len(matched_components)
                total_branching = matched_components[0].get('total_branching', 0.0)
                sum_component_branching = sum([c['branching'] for c in matched_components])
                sum_component_branching_ratio = sum_component_branching / total_branching if total_branching > 0 else 0.0

                # Calculate total fraction covered by union of all carbon atoms in components
                union_carbon_atoms = set()
                # Use molecule with hydrogens for consistent atom indexing
                total_carbons = sum(1 for atom in mol_with_h.GetAtoms() if atom.GetSymbol() == 'C')
                extra_carbons_count = 0

                for comp_dict in matched_components:
                    component = set(comp_dict.get('component', []))
                    smarts_matches = comp_dict.get('smarts_matches')

                    # Add carbon atoms from component
                    for atom_idx in component:
                        if mol_with_h.GetAtomWithIdx(atom_idx).GetSymbol() == 'C':
                            union_carbon_atoms.add(atom_idx)

                    # Add carbon atoms from SMARTS matches
                    if smarts_matches is not None:
                        for atom_idx in smarts_matches:
                            if mol_with_h.GetAtomWithIdx(atom_idx).GetSymbol() == 'C':
                                union_carbon_atoms.add(atom_idx)

                    # Add extra carbons from functional groups
                    if 'smarts_extra_atoms' in comp_dict:
                        extra_carbons_count += comp_dict['smarts_extra_atoms']

                # Total = union of matched carbons + extra functional group carbons
                # Cap at total_carbons to avoid exceeding 1.0
                total_carbon_count = min(len(union_carbon_atoms) + extra_carbons_count, total_carbons)
                total_components_fraction = total_carbon_count / total_carbons if total_carbons > 0 else 0.0

                # Graph structure metrics summaries
                diameters = [c['diameter'] for c in matched_components if not (isinstance(c['diameter'], float) and (c['diameter'] != c['diameter'] or c['diameter'] == float('inf')))]
                radii = [c['radius'] for c in matched_components if not (isinstance(c['radius'], float) and (c['radius'] != c['radius'] or c['radius'] == float('inf')))]
                resistances = [c['effective_graph_resistance'] for c in matched_components if not (isinstance(c['effective_graph_resistance'], float) and (c['effective_graph_resistance'] != c['effective_graph_resistance'] or c['effective_graph_resistance'] == float('inf')))]

                mean_diameter = sum(diameters)/len(diameters) if len(diameters) > 0 else float('nan')
                mean_radius = sum(radii)/len(radii) if len(radii) > 0 else float('nan')
                mean_resistance = sum(resistances)/len(resistances) if len(resistances) > 0 else float('nan')

                def _nanmean(vals):
                    clean = [v for v in vals if v == v and v != float('inf')]
                    return sum(clean) / len(clean) if clean else float('nan')

                # BDE-weighted EGR on expanded component
                egr_bde_vals = [c.get('effective_graph_resistance_BDE', float('nan')) for c in matched_components]
                mean_effective_graph_resistance_BDE = _nanmean(egr_bde_vals)

                # Distance metrics summaries
                min_dists_bc = [c['min_dist_to_barycenter'] for c in matched_components if c['min_dist_to_barycenter'] < float('inf')]
                min_dists_center = [c['min_dist_to_center'] for c in matched_components if c['min_dist_to_center'] < float('inf')]
                max_dists_periph = [c['max_dist_to_periphery'] for c in matched_components if c['max_dist_to_periphery'] > 0]

                mean_dist_to_barycenter = sum(min_dists_bc)/len(min_dists_bc) if len(min_dists_bc) > 0 else 0
                mean_dist_to_center = sum(min_dists_center)/len(min_dists_center) if len(min_dists_center) > 0 else 0
                mean_dist_to_periphery = sum(max_dists_periph)/len(max_dists_periph) if len(max_dists_periph) > 0 else 0

                summary_metrics = {
                    'mean_branching': mean_branching,
                    'total_branching': total_branching,
                    'sum_component_branching_ratio': sum_component_branching_ratio,
                    'mean_smarts_centrality': mean_smarts_centrality,
                    'mean_component_fraction': mean_component_fraction,
                    'total_components_fraction': total_components_fraction,
                    'mean_eccentricity': mean_mean_eccentricity,
                    'median_eccentricity': mean_median_eccentricity,
                    'mean_diameter': mean_diameter,
                    'mean_radius': mean_radius,
                    'mean_effective_graph_resistance': mean_resistance,
                    'mean_effective_graph_resistance_BDE': mean_effective_graph_resistance_BDE,
                    'mean_dist_to_barycenter': mean_dist_to_barycenter,
                    'mean_dist_to_center': mean_dist_to_center,
                    'mean_dist_to_periphery': mean_dist_to_periphery,
                }
            else:
                summary_metrics = {
                    'mean_branching': 0.0,
                    'total_branching': 0.0,
                    'sum_component_branching_ratio': 0.0,
                    'mean_smarts_centrality': 0.0,
                    'mean_component_fraction': 0.0,
                    'total_components_fraction': 0.0,
                    'mean_eccentricity': 0.0,
                    'median_eccentricity': 0.0,
                    'mean_diameter': float('nan'),
                    'mean_radius': float('nan'),
                    'mean_effective_graph_resistance': float('nan'),
                    'mean_effective_graph_resistance_BDE': float('nan'),
                    'mean_dist_to_barycenter': 0,
                    'mean_dist_to_center': 0,
                    'mean_dist_to_periphery': 0,
                }

            match_results.append({
                'match_id': f"G{group.id}",
                'id': group.id,
                'group_name': group.name,
                'match_count': match_count,
                'components_sizes': components_sizes,
                'num_components': len(matched_components),
                'components': matched_components,
                'components_types': list(set(str(c['SMARTS']) if isinstance(c['SMARTS'], list) else c['SMARTS'] for c in matched_components)),
                'type':'HalogenGroup',
                **summary_metrics
            })

        results[inchikey].setdefault('matches',[]).extend(match_results)
    if include_PFAS_definitions is True:
        for i, mol in enumerate(mols):
            formula = CalcMolFormula(mol)
            formula_dict = n_from_formula(formula)
            definitions = parse_definitions_in_mol(mol, formula=formula_dict, **kwargs)
            inchikey = Chem.MolToInchiKey(mol)
            results.setdefault(inchikey,{
                        "smiles": Chem.MolToSmiles(mol),
                        "inchikey": inchikey,
                        "inchi": Chem.MolToInchi(mol),
                        "formula": formula}).setdefault("matches",[]).extend([
                            {'match_id': f"D{definition.id}",
                            'id': definition.id,
                            'definition_name': definition.name,
                            'type':'PFASdefinition'} for definition in definitions])
    # Convert results to list format (one entry per molecule)
    results_list = [r for r in results.values()]
    # Format output based on requested format
    if output_format in ['dataframe', 'csv']:
        import pandas as pd
        rows = []
        for entry in results_list:
            for match in entry['matches']:
                if match['type'] == 'HalogenGroup':
                    rows.append({
                        'smiles': smi,
                        'match_id': match['match_id'],
                        'match_name': match['group_name'],
                        'match_count': match['match_count'],
                        'components_sizes': match['components_sizes'],
                        'num_chains': match['num_chains'],
                        'match_type': match['type']
                    })
                elif match['type'] == 'PFASdefinition':
                    rows.append({
                        'smiles': smi,
                        'match_id': match['match_id'],
                        'match_name': match['definition_name'],
                        'match_type': match['type']
                    })
        df = pd.DataFrame(rows)
        return df.to_csv(index=False) if output_format == 'csv' else df

    # Default: list-like container with navigation helpers that
    # remains JSON-serialisable and behaves like a normal list of
    # dicts for existing consumers.
    return PFASEmbeddingSet(results_list)

def compile_componentSmarts(chain_smarts, end_smarts):
    """
    Compile a pair of SMARTS patterns into a ready-to-use path definition.

    This function preprocesses SMARTS patterns for chain and end groups,
    preparing them for use in PFAS parsing functions.

    Parameters
    ----------
    chain_smarts : str
        SMARTS pattern for the repeating chain unit
    end_smarts : str
        SMARTS pattern for the terminal group

    Returns
    -------
    list
        List containing [chain_mol, end_mol] where both are preprocessed RDKit Mol objects

    Examples
    --------
    >>> chain = compile_componentSmarts(
    ...     "[C;X4](Cl)(Cl)!@!=!#[C;X4](Cl)(Cl)",
    ...     "[C;X4](Cl)(Cl)Cl"
    ... )
    >>> paths = {'Perchlorinated': chain}
    """
    chain_mol = Chem.MolFromSmarts(chain_smarts)
    chain_mol.UpdatePropertyCache()
    Chem.GetSymmSSSR(chain_mol)
    chain_mol.GetRingInfo().NumRings()

    end_mol = Chem.MolFromSmarts(end_smarts)
    end_mol.UpdatePropertyCache()
    Chem.GetSymmSSSR(end_mol)
    end_mol.GetRingInfo().NumRings()

    return [chain_mol, end_mol]

def compile_componentSmartss(paths_dict):
    """
    Compile multiple SMARTS path definitions from a dictionary.

    This function takes a dictionary of path definitions (with 'component' and 'end' keys)
    and preprocesses them for use in PFAS parsing functions.

    Parameters
    ----------
    paths_dict : dict
        Dictionary with structure::

            {
                'PathName': {'component': 'SMARTS', 'end': 'SMARTS'},
                ...
            }

    Returns
    -------
    dict
        Dictionary mapping path names to [chain_mol, end_mol] pairs

    Examples
    --------
    >>> custom_paths = {
    ...     'Perchlorinated': {
    ...         'component': '[C;X4](Cl)(Cl)!@!=!#[C;X4](Cl)(Cl)',
    ...         'end': '[C;X4](Cl)(Cl)Cl'
    ...     },
    ...     'MixedHalo': {
    ...         'component': '[C;X4]([F,Cl])([F,Cl])!@!=!#[C;X4]([F,Cl])',
    ...         'end': '[C;X4]([F,Cl])([F,Cl])[F,Cl]'
    ...     }
    ... }
    >>> compiled = compile_componentSmartss(custom_paths)
    >>> results = parse_smiles(smiles, componentSmartss=compiled)
    """
    compiled = {}
    for name, patterns in paths_dict.items():
        if isinstance(patterns, dict):
            if 'smarts' in patterns:
                smarts_val = patterns['smarts']
                if isinstance(smarts_val, str):
                    smarts_val = Chem.MolFromSmarts(smarts_val)
                    if smarts_val is None:
                        raise ValueError(f"Invalid SMARTS for path '{name}'")
                    smarts_val.UpdatePropertyCache()
                    Chem.GetSymmSSSR(smarts_val)
                    smarts_val.GetRingInfo().NumRings()
                compiled[name] = {
                    **patterns,
                    'smarts': smarts_val,
                    'component': smarts_val
                }
            elif 'component' in patterns and 'end' in patterns:
                compiled[name] = compile_componentSmarts(patterns['component'], patterns['end'])
            elif 'component' in patterns:
                smarts_val = patterns['component']
                if isinstance(smarts_val, str):
                    smarts_val = Chem.MolFromSmarts(smarts_val)
                    if smarts_val is None:
                        raise ValueError(f"Invalid SMARTS for path '{name}'")
                    smarts_val.UpdatePropertyCache()
                    Chem.GetSymmSSSR(smarts_val)
                    smarts_val.GetRingInfo().NumRings()
                compiled[name] = {
                    **patterns,
                    'smarts': smarts_val,
                    'component': smarts_val
                }
            else:
                raise ValueError(f"Path '{name}' must define 'smarts' or 'component' (and optionally 'end')")
        else:
            raise ValueError(f"Path '{name}' must be a dict")
    return compiled
