"""
PFASgroups core functions.

This module provides functions for parsing and plotting PFAS groups.
"""

import numpy as np

from rdkit import Chem
from rdkit.Chem.rdMolDescriptors import CalcMolFormula
from .PFASGroupModel import PFASGroup
from .PFASDefinitionModel import PFASDefinition
from .ComponentsSolverModel import ComponentsSolver
from .core import fragment_until_valence_is_correct, n_from_formula, add_componentSmarts, PFAS_DEFINITIONS_FILE, PFAS_GROUPS_FILE, rdkit_disable_log
import json



# --- Load PFAS groups from PFAS_groups_smarts.json ---
def load_PFASGroups():
    """
    Adds default PFASGroups to function
    """
    with open(PFAS_GROUPS_FILE,'r') as f:
        _pfg = json.load(f)
    pfg = [PFASGroup(**x) for x in _pfg if x.get('compute',True)]
    agg_pfg = [PFASGroup(**x) for x in _pfg if not x.get('compute',True)]
    # list of PFASgroup names
    pfg_names =  {pf.name:pf.id for pf in pfg}
    # list groups aggregated by groups with compute=FALSE
    agg_pfg = {ppf:list(map(pfg_names.get,list(filter(ppf.re_search.search,pfg_names.keys())))) for ppf in agg_pfg}
    def inner(func):
        def wrapper(*args,**kwargs):
            kwargs['pfas_groups'] = kwargs.get('pfas_groups',pfg)
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
        def wrapper(*args,**kwargs):
            mol = args[0]
            # Add hydrogens to molecule before creating ComponentsSolver
            # This ensures atom indices are consistent throughout the analysis
            mol = Chem.AddHs(mol)
            args = list(args)  # Convert to mutable list
            args[0] = mol
            args = tuple(args)  # Convert back to tuple
            # Pass through component metric options
            solver_kwargs = {}
            if 'limit_effective_graph_resistance' in kwargs:
                solver_kwargs['limit_effective_graph_resistance'] = kwargs['limit_effective_graph_resistance']
            if 'compute_component_metrics' in kwargs:
                solver_kwargs['compute_component_metrics'] = kwargs['compute_component_metrics']
            with ComponentsSolver(mol, **solver_kwargs) as fluorinated_components_dict:
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


# --- Main PFAS group parsing functions ---
@add_componentSmarts()
@load_PFASGroups()
@load_componentsSolver()
def parse_groups_in_mol(mol, fluorinated_components_dict=None, pfas_groups = None, **kwargs):
    """Iterates over PFAS groups and finds the ones that match the molecule.
    
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
        List of (PFASGroup, match_count, component_sizes, matched_components) tuples where:
        - PFASGroup: The matched PFAS group object
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
    For PFASgroups 
    1. with componentSmarts = 'cyclic', search for connected component matching first smarts
    2. with multiple smarts patterns or counts > 1: Find components containing all required SMARTS matches with minimum counts
    3. with single smarts pattern: Search for connected components of fluorinated atoms (for each pathType) where smarts match is in the component
    4. with smarts defined and formula constraints: search for substructure matches of smarts, and for given componentSmarts for the PFASgroup, search for connected components of fluorinated atoms where smarts match is in the component

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
    # PFAS groups
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
    
    # Process aggregate PFAS groups efficiently
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
                    if agg_group.componentSmarts is not None and smarts_type != agg_group.componentSmarts:
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

    return group_matches


@rdkit_disable_log(level='warning')
def parse_smiles(smiles, bycomponent=False, output_format='list', 
                  limit_effective_graph_resistance=None, compute_component_metrics=True, **kwargs):
    """
    Parse SMILES string(s) and return PFAS group information.
    
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
    return parse_mols(mol_list, bycomponent=bycomponent, output_format=output_format, **kwargs)

def parse_mol(mol, **kwargs):
    """Wrapper for parse_mols to handle single molecule input."""
    return parse_mols([mol], **kwargs)[0]

def parse_mols(mols, output_format='list', include_PFAS_definitions=True, 
               limit_effective_graph_resistance=None, compute_component_metrics=True, **kwargs):
    """
    Parse RDKit molecule(s) and return PFAS group information.
    
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
    **kwargs : dict
        Additional parameters (pfas_groups, componentSmartss, etc.)
    
    Returns:
    --------
    list, pandas.DataFrame, or str
        Depends on output_format parameter
    """
    
    # Parse all molecules
    results = {}
    for mol in mols:
        # Add hydrogens to ensure consistent atom indexing
        mol_with_h = Chem.AddHs(mol)
        formula = CalcMolFormula(mol)
        bycomponent = kwargs.pop('bycomponent', False)
        # Pass through component metric options
        kwargs['limit_effective_graph_resistance'] = limit_effective_graph_resistance
        kwargs['compute_component_metrics'] = compute_component_metrics
        matches = parse_groups_in_mol(mol, formula=formula, bycomponent=bycomponent, **kwargs)
        inchikey = Chem.MolToInchiKey(mol)
        inchi = Chem.MolToInchi(mol)
        smi = Chem.MolToSmiles(mol)
        results.setdefault(inchikey,{}).update({"smiles": smi,
                        "inchikey": inchikey,
                        "inchi": inchi,
                        "formula": formula})
        
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
                
                # Distance metrics summaries
                min_dists_bc = [c['min_dist_to_barycenter'] for c in matched_components if c['min_dist_to_barycenter'] < float('inf')]
                min_dists_center = [c['min_dist_to_center'] for c in matched_components if c['min_dist_to_center'] < float('inf')]
                max_dists_periph = [c['max_dist_to_periphery'] for c in matched_components if c['max_dist_to_periphery'] > 0]
                
                mean_dist_to_barycenter = sum(min_dists_bc)/len(min_dists_bc) if len(min_dists_bc) > 0 else 0
                mean_dist_to_center = sum(min_dists_center)/len(min_dists_center) if len(min_dists_center) > 0 else 0
                mean_dist_to_periphery = sum(max_dists_periph)/len(max_dists_periph) if len(max_dists_periph) > 0 else 0
                
                summary_metrics = {
                    'mean_branching': mean_branching,
                    'mean_smarts_centrality': mean_smarts_centrality,
                    'mean_component_fraction': mean_component_fraction,
                    'total_components_fraction': total_components_fraction,
                    'mean_eccentricity': mean_mean_eccentricity,
                    'median_eccentricity': mean_median_eccentricity,
                    'mean_diameter': mean_diameter,
                    'mean_radius': mean_radius,
                    'mean_effective_graph_resistance': mean_resistance,
                    'mean_dist_to_barycenter': mean_dist_to_barycenter,
                    'mean_dist_to_center': mean_dist_to_center,
                    'mean_dist_to_periphery': mean_dist_to_periphery
                }
            else:
                summary_metrics = {
                    'mean_branching': 0.0,
                    'mean_smarts_centrality': 0.0,
                    'mean_component_fraction': 0.0,
                    'total_components_fraction': 0.0,
                    'mean_eccentricity': 0.0,
                    'median_eccentricity': 0.0,
                    'mean_diameter': float('nan'),
                    'mean_radius': float('nan'),
                    'mean_effective_graph_resistance': float('nan'),
                    'mean_dist_to_barycenter': 0,
                    'mean_dist_to_center': 0,
                    'mean_dist_to_periphery': 0
                }
            
            match_results.append({
                'match_id': f"G{group.id}",
                'id': group.id,
                'group_name': group.name,
                'match_count': match_count,
                'components_sizes': components_sizes,
                'num_components': len(matched_components),
                'components': matched_components,
                'components_types': [x for x in set([c['SMARTS'] for c in matched_components])],
                'type':'PFASgroup',
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
    # Convert results to list format
    results = [r for r in results.values()]
    # Format output based on requested format
    if output_format in ['dataframe', 'csv']:
        import pandas as pd
        rows = []
        for entry in results:
            for match in entry['matches']:
                if match['type'] == 'PFASgroup':
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
    return results

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
        if isinstance(patterns, dict) and 'component' in patterns and 'end' in patterns:
            compiled[name] = compile_componentSmarts(patterns['component'], patterns['end'])
        else:
            raise ValueError(f"Path '{name}' must have 'component' and 'end' keys")
    return compiled