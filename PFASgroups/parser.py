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
from .core import fragment_until_valence_is_correct, n_from_formula, add_smartsPath, add_smarts, PFAS_DEFINITIONS_FILE, PFAS_GROUPS_FILE
from typing import Union, List, Dict
from PIL import Image
from io import BytesIO
import re
import json
import os
from itertools import product, groupby
import svgutils.transform as sg
from .draw_mols import draw_images, plot_mols




# --- Load PFAS groups from PFAS_groups_smarts.json ---
def load_PFASGroups():
    """
    Adds default PFASGroups to function
    """
    with open(PFAS_GROUPS_FILE,'r') as f:
        pfg = json.load(f)
    pfg = [PFASGroup(**x) for x in pfg]
    def inner(func):
        def wrapper(*args,**kwargs):
            kwargs['pfas_groups'] = kwargs.get('pfas_groups',pfg)
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
            with ComponentsSolver(mol) as fluorinated_components_dict:
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
@add_smartsPath()
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
        Additional parameters (formula, pfas_groups, smartsPaths, etc.)
    
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
    1. with smartsPath = 'cyclic', search for connected component matching first smarts
    2. with multiple smarts patterns or counts > 1: Find components containing all required SMARTS matches with minimum counts
    3. with single smarts pattern: Search for connected components of fluorinated atoms (for each pathType) where smarts match is in the component
    4. with smarts defined and formula constraints: search for substructure matches of smarts, and for given smartsPath for the PFASgroup, search for connected components of fluorinated atoms where smarts match is in the component

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

    # PFAS groups
    group_matches = []
    for pf in pfas_groups:
        #logger.debug(f"{pf.name}")
        matched1_len = 0
        for fd,mol in zip(formulas,frags):
            match = pf.find_components(mol, fd, fluorinated_components_dict, **kwargs)
            if match is not None and len(match)>0:
                group_matches.extend(match)
    return group_matches

def parse_smiles(smiles, bycomponent=False, output_format='list', **kwargs):
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
    **kwargs : dict
        Additional parameters (pfas_groups, smartsPaths, etc.)
    
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
    return parse_mols(mol_list, bycomponent=bycomponent, output_format=output_format, **kwargs)

def parse_mol(mol, **kwargs):
    """Wrapper for parse_mols to handle single molecule input."""
    return parse_mols([mol], **kwargs)[0]

def parse_mols(mols, output_format='list', include_PFAS_definitions=True, **kwargs):
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
    **kwargs : dict
        Additional parameters (pfas_groups, smartsPaths, etc.)
    
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

@load_PFASGroups()
def generate_fingerprint(smiles: Union[str, List[str]], 
                             selected_groups: Union[List[int], range, None] = None,
                             representation: str = 'vector',
                             count_mode: str = 'binary',
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
        - 'count': Number of matches found
        - 'max_component': Maximum component size for groups with components
    pfas_groups : list, optional
        Custom list of PFAS groups. If None, uses default groups.
    
    Returns:
    --------
    fingerprints : numpy.ndarray, dict, or list
        Fingerprint representation(s) based on 'representation' parameter.
        For single SMILES input, returns single fingerprint.
        For list input, returns list of fingerprints.
    group_info : dict
        Information about the groups used in fingerprinting:
        - 'group_names': List of group names in fingerprint order
        - 'group_ids': List of group IDs in fingerprint order
        - 'selected_indices': Indices of selected groups
    
    Examples:
    ---------
    >>> # Binary vector for specific groups
    >>> fp, info = generate_pfas_fingerprint('CC(F)(F)C(=O)O', 
    ...                                       selected_groups=range(28, 52))
    
    >>> # Count-based dictionary representation
    >>> fp, info = generate_pfas_fingerprint(['CCF', 'CCFF'], 
    ...                                       representation='dict',
    ...                                       count_mode='count')
    
    >>> # Detailed match information
    >>> fp, info = generate_pfas_fingerprint('PFOA_SMILES', 
    ...                                       representation='detailed')
    """
    pfas_groups = kwargs.get('pfas_groups')
    # Handle single SMILES input
    if isinstance(smiles, str):
        smiles_list = [smiles]
        single_input = True
    else:
        smiles_list = smiles
        single_input = False
    
    # Determine which groups to use
    if selected_groups is None:
        # Use all groups
        selected_indices = list(range(len(pfas_groups)))
        selected_pfas_groups = pfas_groups
    else:
        # Convert range to list if needed
        if isinstance(selected_groups, range):
            selected_indices = list(selected_groups)
        else:
            selected_indices = selected_groups
            
        # Validate indices
        max_index = len(pfas_groups) - 1
        invalid_indices = [i for i in selected_indices if i < 0 or i > max_index]
        if invalid_indices:
            raise ValueError(f"Invalid group indices {invalid_indices}. "
                           f"Valid range is 0-{max_index}")
        
        selected_pfas_groups = [pfas_groups[i] for i in selected_indices]
    
    # Create group info
    group_info = {
        'group_names': [group.name for group in selected_pfas_groups],
        'group_ids': [group.id for group in selected_pfas_groups],
        'selected_indices': selected_indices,
        'total_groups': len(pfas_groups)
    }
    
    fingerprints = []
    
    for smiles_str in smiles_list:
        try:
            mol = Chem.MolFromSmiles(smiles_str)
            if mol is None:
                raise ValueError(f"Invalid SMILES: {smiles_str}")
            
            formula = CalcMolFormula(mol)
            all_matches = parse_groups_in_mol(mol, formula=formula, pfas_groups=pfas_groups,include_PFAS_definitions=False, **kwargs)
            
            # Create mapping from group ID to match information
            match_dict = {}
            for group, match_count, component_sizes, matched_components in all_matches:
                match_dict[group.id] = {
                    'group': group,
                    'match_count': match_count,
                    'component_sizes': component_sizes,
                    'matched_components': matched_components
                }
            
            if representation == 'vector' or representation == 'int':
                # Create binary or count vector
                fingerprint = np.zeros(len(selected_pfas_groups), dtype=int)
                for i, group in enumerate(selected_pfas_groups):
                    if group.id in match_dict:
                        match_info = match_dict[group.id]
                        if count_mode == 'binary' or representation == 'int':
                            fingerprint[i] = 1
                        elif count_mode == 'count':
                            fingerprint[i] = match_info['match_count']
                        elif count_mode == 'max_component':
                            if match_info['matched_components']:
                                fingerprint[i] = max([comp['size'] for comp in match_info['matched_components']])
                            else:
                                fingerprint[i] = match_info['match_count'] if match_info['match_count'] > 0 else 0
                        else:
                            raise ValueError(f"Unknown count_mode: {count_mode}")
                if representation == 'int':
                    fingerprints.append([int(''.join([str(x) for x in f]),2) for f in fingerprint])
                else:
                    fingerprints.append(fingerprint)
                
            elif representation == 'dict':
                # Create dictionary with all groups (including zeros)
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
                                fingerprint[group.name] = max([comp['size'] for comp in match_info['matched_components']])
                            else:
                                fingerprint[group.name] = match_info['match_count'] if match_info['match_count'] > 0 else 0
                    else:
                        fingerprint[group.name] = 0
                fingerprints.append(fingerprint)
                
            elif representation == 'sparse':
                # Create dictionary with only non-zero entries
                fingerprint = {}
                for group in selected_pfas_groups:
                    if group.id in match_dict:
                        match_info = match_dict[group.id]
                        value = 0
                        if count_mode == 'binary':
                            value = 1
                        elif count_mode == 'count':
                            value = match_info['match_count']
                        elif count_mode == 'max_component':
                            if match_info['matched_components']:
                                value = max([comp['size'] for comp in match_info['matched_components']])
                            else:
                                value = match_info['match_count'] if match_info['match_count'] > 0 else 0
                        
                        if value > 0:
                            fingerprint[group.name] = value
                fingerprints.append(fingerprint)
                
            elif representation == 'detailed':
                # Return full match information for selected groups
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
    
    # Return single fingerprint for single input, list for multiple inputs
    if single_input:
        return fingerprints[0], group_info
    else:
        return fingerprints, group_info

@add_smartsPath()
def plot_pfasgroups(smiles: Union[list, str], display=True, path=None, svg=False, ipython=False, subwidth=300, subheight=300, ncols=2, addAtomIndices=True, addBondIndices=False, paths=[0, 1, 2, 3], split_matches = False, SMARTS=None, **kwargs):
    """
    Plot PFAS group assignments for a list of SMILES strings.

    :params smiles: List of SMILES strings or a single SMILES string.
    :params display: Whether to display the plot.
    :params path: Path to save the plot image.
    :params svg: Whether to generate SVG images.
    :params ipython: Whether to display in an IPython environment.
    :params subwidth: Width of each sub-image.
    :params subheight: Height of each sub-image.
    :params ncols: Number of columns in the grid layout.
    :params addAtomIndices: Whether to add atom indices to the plot.
    :params addBondIndices: Whether to add bond indices to the plot.
    :params paths: List of PFAS group indices or names to include in the plot.
    :params split_matches: Whether to create separate images for each match.
    :params SMARTS: Optional SMARTS pattern to highlight in the plots.
    :params kwargs: Additional keyword arguments for customization.
    """
    from rdkit.Chem import Draw
    if isinstance(smiles, str):
        smiles = [smiles]
    imgs = []
    path_names = list(kwargs.get('smartsPaths',{'Perfluoroalkyl':'Perfluoroalkyl','Polyfluoroalkyl':'Polyfluoroalkyl'}).keys())
    for i, s in enumerate(paths):
        if isinstance(s, int):
            paths[i] = path_names[s]
    def draw_subfig(legend, atoms=[]):
        if svg is True:
            d2d = Draw.MolDraw2DSVG(subwidth, subheight)
        else:
            d2d = Draw.MolDraw2DCairo(subwidth, subheight)
        dopts = d2d.drawOptions()
        dopts.useBWAtomPalette()
        dopts.fixedBondLength = 20
        dopts.addAtomIndices = addAtomIndices
        dopts.addBondIndices = addBondIndices
        dopts.maxFontSize = 16
        dopts.minFontSize = 13
        d2d.DrawMolecule(mol, legend=legend, highlightAtoms=atoms)
        d2d.FinishDrawing()
        return d2d.GetDrawingText()
    for i, s in enumerate(smiles):
        mol = Chem.MolFromSmiles(s)
        mol = Chem.AddHs(mol)
        matches = parse_groups_in_mol(mol, CalcMolFormula(mol), bycomponent=kwargs.get('bycomponent',False))
        highlight_atoms = []
        for pf, n, n_cfchains, match_indices in matches:
            for match in match_indices:
                if SMARTS is None or match['SMARTS'] in SMARTS:
                    highlight_atoms.extend(match['chain'])
                    if split_matches is True:
                        new_img = draw_subfig(f"{pf.name}, {match['SMARTS']}", atoms=match['chain'])
                        imgs.append(new_img)
        if split_matches is False:
            new_img = draw_subfig(f"{SMARTS if SMARTS is not None else ''}", atoms=highlight_atoms)
            imgs.append(new_img)
    if len(imgs) == 0:
        if svg is True:
            d2d = Draw.MolDraw2DSVG(subwidth, subheight)
        else:
            d2d = Draw.MolDraw2DCairo(subwidth, subheight)
        dopts = d2d.drawOptions()
        dopts.useBWAtomPalette()
        dopts.fixedBondLength = 20
        dopts.addAtomIndices = addAtomIndices
        dopts.addBondIndices = addBondIndices
        dopts.maxFontSize = 16
        dopts.minFontSize = 13
        d2d.DrawMolecule(mol)
        d2d.FinishDrawing()
        imgs.append(d2d.GetDrawingText())
    # For now, just return the images as PIL Images
    if svg is True:
        imgs = [sg.fromstring(img) for img in imgs]
    else:
        imgs = [Image.open(BytesIO(img)) for img in imgs]
    # Simple grid
    return draw_images(imgs, buffer = kwargs.get('buffer',2), ncols = ncols, svg = svg)
    width = subwidth * min(ncols, len(imgs))
    height = subheight * ((len(imgs) + ncols - 1) // ncols)
    grid = Image.new('RGBA', (width, height), (255, 255, 255, 0))
    for idx, img in enumerate(imgs):
        x = (idx % ncols) * subwidth
        y = (idx // ncols) * subheight
        grid.paste(img, (x, y))
    if path is not None:
        grid.save(path)
    if display:
        grid.show()
    return grid, width, height

@add_smartsPath()
def get_smartsPaths(**kwargs):
    return kwargs.get('smartsPaths')

@load_PFASGroups()
def get_PFASGroups(**kwargs):
    return kwargs.get('pfas_groups')

@load_PFASDefinitions()
def get_PFASDefinitions(**kwargs):
    return kwargs.get('pfas_definitions')

def compile_smartsPath(chain_smarts, end_smarts):
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
    >>> chain = compile_smartsPath(
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

def compile_smartsPaths(paths_dict):
    """
    Compile multiple SMARTS path definitions from a dictionary.
    
    This function takes a dictionary of path definitions (with 'chain' and 'end' keys)
    and preprocesses them for use in PFAS parsing functions.
    
    Parameters
    ----------
    paths_dict : dict
        Dictionary with structure:
        {
            'PathName': {'chain': 'SMARTS', 'end': 'SMARTS'},
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
    ...         'chain': '[C;X4](Cl)(Cl)!@!=!#[C;X4](Cl)(Cl)',
    ...         'end': '[C;X4](Cl)(Cl)Cl'
    ...     },
    ...     'MixedHalo': {
    ...         'chain': '[C;X4]([F,Cl])([F,Cl])!@!=!#[C;X4]([F,Cl])',
    ...         'end': '[C;X4]([F,Cl])([F,Cl])[F,Cl]'
    ...     }
    ... }
    >>> compiled = compile_smartsPaths(custom_paths)
    >>> results = parse_smiles(smiles, smartsPaths=compiled)
    """
    compiled = {}
    for name, patterns in paths_dict.items():
        if isinstance(patterns, dict) and 'chain' in patterns and 'end' in patterns:
            compiled[name] = compile_smartsPath(patterns['chain'], patterns['end'])
        else:
            raise ValueError(f"Path '{name}' must have 'chain' and 'end' keys")
    return compiled