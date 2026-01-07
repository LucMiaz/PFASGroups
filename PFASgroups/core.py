"""
PFASgroups core functions.

This module provides functions for parsing and plotting PFAS groups.
"""

import numpy as np
import networkx as nx
from rdkit import Chem
from rdkit.Chem.rdMolDescriptors import CalcMolFormula
from .PFASGroupModel import PFASGroup
from .PFASDefinitionModel import PFASDefinition
from typing import Union, List, Dict
from PIL import Image
from io import BytesIO
import re
import json
import os
from itertools import product, groupby
import svgutils.transform as sg
from .draw_mols import draw_images, plot_mols

PATH_NAMES = ['Perfluoroalkyl','Polyfluoroalkyl','Polyfluoro','Polyfluorobr']

# --- Load SMARTS paths from fpaths.json ---
MODULE_DIR = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = os.path.join(MODULE_DIR, 'data')
FPATHS_FILE = os.path.join(DATA_DIR, 'fpaths.json')
PFAS_GROUPS_FILE = os.path.join(DATA_DIR, 'PFAS_groups_smarts.json')
PFAS_DEFINITIONS_FILE = os.path.join(DATA_DIR, 'PFAS_definitions_smarts.json')


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

# --- Add SMARTS paths to function ---
def add_smartsPath(names =PATH_NAMES, filename = FPATHS_FILE):
    """Yields SMARTS for chains"""
    paths = {}
    pathsEnd = {}
    with open(filename,'r') as f:
        fpaths = json.load(f) 
    for n in names:
        s = fpaths[n]['chain']
        e = fpaths[n]['end']
        smol = Chem.MolFromSmarts(s)
        smol.UpdatePropertyCache()
        Chem.GetSymmSSSR(smol)
        smol.GetRingInfo().NumRings()
        emol = Chem.MolFromSmarts(e)
        emol.UpdatePropertyCache()
        Chem.GetSymmSSSR(emol)
        emol.GetRingInfo().NumRings()
        paths[n] = [smol,emol]
    def inner(func):
        def wrapper(*args,**kwargs):
            kwargs["smartsPaths"] = kwargs.get("smartsPaths",paths)
            return func(*args, **kwargs)
        return wrapper
    return inner

def add_smarts(name = 'smarts'):
    """Yields preprocessed SMARTS to decorated function"""
    smarts = {}
    def add(s):
        smol = Chem.MolFromSmarts(s)
        smol.UpdatePropertyCache()
        Chem.GetSymmSSSR(smol)
        smol.GetRingInfo().NumRings()
        smarts[s]=smol
        return smol
    def get(s):
        return smarts.get(s,add(s))
    def inner(func):
        def wrapper(*args,**kwargs):
            kwargs[name] = get(kwargs.get(name))
            return func(*args, **kwargs)
        return wrapper
    return inner

# --- Helper functions ---

def remove_atoms(mol, idxs, removable = ['H','F','Cl','Br','I'], show_on_error = False):
    """Remove atoms by indices and maintain connectivity.
    
    This function removes the specified atoms and their removable neighbors,
    then reconnects the remaining structure to maintain molecular integrity.
    """
    if not idxs:
        return mol
    
    to_remove = set()
    # Map each removed atom to its non-removable neighbors
    removed_to_neighbors = {}
    
    # First pass: identify all atoms to remove and their connections
    for idx in idxs:
        atom = mol.GetAtomWithIdx(idx)
        neighbors_r = [x.GetIdx() for x in atom.GetNeighbors() if x.GetSymbol() in removable]
        neighbors_c = [x.GetIdx() for x in atom.GetNeighbors() if x.GetSymbol() not in removable]
        
        # Add the atom and its removable neighbors to removal list
        to_remove.add(idx)
        to_remove.update(neighbors_r)
        
        # Store non-removable neighbors for reconnection
        if neighbors_c:
            removed_to_neighbors[idx] = neighbors_c
    
    # Build a graph of connectivity between removed atoms and their neighbors
    # to determine how to reconnect the structure
    new_bonds = []
    processed_chains = set()
    
    # Process each chain of consecutive removed atoms
    for start_idx in idxs:
        if start_idx in processed_chains:
            continue
            
        # Find the chain of consecutive removed atoms containing start_idx
        chain = [start_idx]
        processed_chains.add(start_idx)
        
        # Extend chain in both directions
        queue = [start_idx]
        while queue:
            current = queue.pop(0)
            for neighbor_idx in [n.GetIdx() for n in mol.GetAtomWithIdx(current).GetNeighbors()]:
                if neighbor_idx in idxs and neighbor_idx not in processed_chains:
                    chain.append(neighbor_idx)
                    processed_chains.add(neighbor_idx)
                    queue.append(neighbor_idx)
        
        # Find the endpoints of this chain (atoms that connect to non-removable parts)
        chain_endpoints = []
        for atom_idx in chain:
            if atom_idx in removed_to_neighbors:
                non_removable_neighbors = [n for n in removed_to_neighbors[atom_idx] if n not in to_remove]
                if non_removable_neighbors:
                    chain_endpoints.extend(non_removable_neighbors)
        
        # Connect the endpoints if there are exactly 2
        if len(chain_endpoints) == 2 and chain_endpoints[0] != chain_endpoints[1]:
            new_bonds.append((chain_endpoints[0], chain_endpoints[1]))
        elif len(chain_endpoints) > 2:
            # For branched structures, don't create connections that would change topology
            # This prevents fragmentation but may not be chemically meaningful
            pass
    
    # Create new molecule without removed atoms
    rwm = Chem.RWMol()
    _rwm = Chem.RWMol(mol)
    Chem.Kekulize(_rwm)
    
    # Map from old atom indices to new atom indices
    old_to_new = {}
    charged_atoms = []
    
    # Add all atoms except those to be removed
    for i, atom in enumerate(_rwm.GetAtoms()):
        if i not in to_remove:
            new_atom = Chem.Atom(atom.GetAtomicNum())
            new_idx = rwm.AddAtom(new_atom)
            # Copy formal charge if present
            if atom.GetFormalCharge() != 0:
                charged_atoms.append(atom.GetIdx())
            old_to_new[i] = new_idx

    # Copy existing bonds that don't involve removed atoms
    for bond in _rwm.GetBonds():
        a1 = bond.GetBeginAtomIdx()
        a2 = bond.GetEndAtomIdx()
        if a1 in old_to_new and a2 in old_to_new:
            rwm.AddBond(old_to_new[a1], old_to_new[a2], bond.GetBondType())

    # Add new bonds to maintain connectivity
    for a, b in new_bonds:
        if a in old_to_new and b in old_to_new and old_to_new[a] != old_to_new[b]:
            # Check if bond already exists to avoid duplication
            existing_bond = rwm.GetBondBetweenAtoms(old_to_new[a], old_to_new[b])
            if existing_bond is None:
                rwm.AddBond(old_to_new[a], old_to_new[b], Chem.BondType.SINGLE)
    
    # Restore formal charges
    for idx in charged_atoms:
        if idx in old_to_new:
            atom = rwm.GetAtomWithIdx(old_to_new[idx])
            atom.SetFormalCharge(_rwm.GetAtomWithIdx(idx).GetFormalCharge())
    
    try:
        Chem.SanitizeMol(rwm)
    except Exception as e:
        if show_on_error is True:
            _mol = rwm.GetMol()
            try:
                img,_,_ = plot_mols([Chem.MolToSmiles(_mol)], subwidth=600, subheight=600, svg=False, addAtomIndices=True, bondLineWidth=0.5, fixedBondLength=15, minFontSize=12)
                img.show()
            except:
                pass
        raise e
    
    return rwm.GetMol()


def n_from_formula(formula:str, element=None)->Union[int,dict]:
    """
    Compute the number of elements (any or one specific)

    :params formula: Formula to parse
    :params element: Element's symbol to find
    :return: Number of elements in the formula
    """
    if element is not None:
        PAT = f"([{element}])"+r"(\d*)"
    else:
        PAT = r"([A-Z][a-z]?)(\d*)"
    mat = re.findall(PAT,formula)
    formula_dict = {}
    for sym,nb in mat:
        if nb != '':
            formula_dict[sym] = formula_dict.setdefault(sym,0) + int(nb)
        else:
            formula_dict[sym] =formula_dict.setdefault(sym,0) + 1
    if element is not None:
        return formula_dict[element]
    return formula_dict

def fragment_on_bond(mol, a1, a2):
    """Fragment a molecule on a bond between atoms a1 and a2 (indices)"""
    bond = mol.GetBondBetweenAtoms(a1, a2)
    mms =  Chem.FragmentOnBonds(mol, [bond.GetIdx()], addDummies=False)
    return [x for x in Chem.GetMolFrags(mms, asMols=True, sanitizeFrags = False)]

def fragment_until_valence_is_correct(mol, frags):
    """Iterate over the molecule and fragment it until  valence is corrected"""
    try:    
        Chem.SanitizeMol(mol)
    except Chem.AtomValenceException as e:
        all = [int(x) for x in re.findall(r"(?<=#\s)(\d)",str(e))]
        if len(all)==0:
            raise e
        neighbours = mol.GetAtomWithIdx(all[0]).GetNeighbors()
        atom = sorted([(mol.GetBondBetweenAtoms(all[0], x.GetIdx()).GetBondType(),x.GetIdx()) for x in neighbours], reverse = True)[0][1] # return index of neighbour with bond of highest degree
        mols = fragment_on_bond(mol, all[0], atom)
        for m in mols:
            try:
                frags = fragment_until_valence_is_correct(m, frags)
            except IndexError:
                # assume that element in frag is actually isolated
                return frags 
        return frags
    else:
        return frags + [mol]


def mol_to_nx(mol):
    """Construct a networkx graph from a molecule"""
    G = nx.Graph()
    for n, atom in enumerate(mol.GetAtoms()):
        element_Z = atom.GetAtomicNum()
        node_params = {"element" : element_Z,
                       "symbol" : atom.GetSymbol()}
        G.add_node(atom.GetIdx(),
                   **node_params
                   )
    for bond in mol.GetBonds():
        edgeOrder=bond.GetBondTypeAsDouble()
        node1 = bond.GetBeginAtom()
        node1 = node1.GetAtomicNum()
        node2 = bond.GetEndAtom()
        node2 = node2.GetAtomicNum()
        G.add_edge(bond.GetBeginAtomIdx(),
                   bond.GetEndAtomIdx(),
                   order=edgeOrder)
    return G

def get_substruct(_mol:Chem.Mol,struct:Chem.Mol):
    """Returns the indices of the atoms in the molecule that match the substructure"""
    return set([x[0] for x in _mol.GetSubstructMatches(struct)])

@add_smartsPath()
def find_path_between_smarts(mol,smarts1,smarts2,G, **kwargs):
    """Iterates over substructures matched by smarts1 and smarts2 and finds fluorinated paths between them.
    Fluorinated atoms are defined by the SMARTS yield by the decorator add_smartsPath.
    if smarts2 is None, uses SMARTS corresponding to smartsPath."""
    chains = {}
    # Extract shortest_path_method from kwargs (default to 'dijkstra' if not provided)
    shortest_path_method = kwargs.pop('shortest_path_method', 'dijkstra')
    smartsPaths = kwargs.get('smartsPaths')
    smartsMatches1 = get_substruct(mol, smarts1)
    if smarts2 is not None:
        digroup_smarts = (Chem.MolToSmarts(smarts1) != Chem.MolToSmarts(smarts2))
        pairs = {smarts2:[(path,get_substruct(mol,smartsPath)) for path, (smartsPath,_)  in smartsPaths.items()]}
    else:
        digroup_smarts = False
        pairs = {_smarts2:[(path,get_substruct(mol, smartsPath))] for path, (smartsPath,_smarts2) in smartsPaths.items()}
    def match_intersection(_pair,_setp):
        for path,pmatch in _pair:
            if len(_setp)>0 and _setp == _setp.intersection(pmatch):
                chains.setdefault(path,[]).append(_setp)
                return path
        return None
    for mi, match1 in enumerate(smartsMatches1):
        for _smarts2, smartsPath in pairs.items():
            smartsMatches2 = get_substruct(mol,_smarts2)
            for mj, match2 in enumerate(smartsMatches2):
                # Fix: Compare SMARTS strings, not object identity, to avoid self-matching
                # For disulfonic acids, allow same-atom matching when we have identical SMARTS patterns
                if (match1!=match2) or digroup_smarts is False or (digroup_smarts is True and mi != mj):# avoid double matching for diacids
                    if match1 == match2:
                        # Special case: both functional groups on same atom (e.g., disulfonic acids)
                        # Create a minimal path containing just the shared atom
                        path_idx = [match1]
                    else:
                        try:
                            path_idx = nx.shortest_path(G, match1, match2, method=shortest_path_method)
                        except nx.NetworkXNoPath:
                            path_idx = []
                    setp = set(path_idx)
                    match_intersection(smartsPath,setp)
    chains = {k:sorted(v, key=len, reverse=True) for k,v in chains.items()}
    for path, chain in chains.items():
        remove = []
        for i in range(0,len(chain)):
            if sum([chain[i].issubset(x) for x in chain[0:i]])>0:
                remove.append(i)
        chains[path] = [x for i,x in enumerate(chain) if i not in remove]
    
    if len(chains.keys())==0:
        n = 0
        chain_lengths = []
    else:
        n = max([0]+[len(x) for x in chains.values()])#{k:len(chain) for k,chain in chains.items()}
        chain_lengths = [len(x) for x in chains.get('Perfluoroalkyl',[[]])]#{k:[len(x) for x in chain] for k,chain in chains.items()}
    matched_chains = [{'chain':list(chain), 'length':len(chain),'SMARTS':path} for path,all_chains in chains.items() for chain in all_chains]
    return n, chain_lengths, len(smartsMatches1), matched_chains  


def path_between_smarts(mol:Chem.Mol, smarts1,smarts2,**kwargs):
    """
    :params compound: Molecule
    :params smarts1, smarts2: SMARTS as string
    :params pathsmarts: SMARTS (parsed as Chem.MolFromSmarts) with constraint for the atoms on the path
    :return: number of chains, and max_length in chain for path between an occurrence of smarts1 and smarts2 with all atoms satisfying pathsmarts."""
    try:
        G = mol_to_nx(mol)
    except ValueError as e:
        raise e
    return find_path_between_smarts(mol,
                                    smarts1,
                                    smarts2,
                                    G, **kwargs)

def connected_components(mol, subset):
    """Find connected components in a molecule."""
    G = mol_to_nx(mol)
    G = G.subgraph(subset)
    components = list(nx.connected_components(G))
    return components

@add_smartsPath()
def get_fluorinated_subgraph(mol, **kwargs):
    """Get the fluorinated indices by connected components of a molecule based on path SMARTS."""
    subsets = {}
    smartsPaths = kwargs.get('smartsPaths')
    for pathName, d in smartsPaths.items():
        path_smarts = d[0] # chain
        matches = mol.GetSubstructMatches(path_smarts)
        subset = [y for x in matches for y in x]
        if len(subset)==0:
            subsets[pathName] = []
            continue
        components = connected_components(mol, subset)
        subsets[pathName] = components
    return subsets

def find_alkyl_components(mol, smarts, components):
    """Find alkyl components in a molecule."""
    matches = mol.GetSubstructMatches(smarts)
    subset = set(y for x in matches for y in x)
    if len(subset)==0:
        return 0, [], 0, []
    # Filter components connected to the smarts
    relevant_components = [comp for comp in components if len(subset.intersection(comp))>0]
    chain_lengths = [len(x) for x in relevant_components]
    # Convert components to dictionary format
    matched_chains = [{'chain': list(comp), 'length': len(comp), 'SMARTS': 'alkyl'} for comp in relevant_components]
    return max([0]+[len(x) for x in relevant_components]), chain_lengths, len(relevant_components), matched_chains

def find_aryl_components(mol, aryl_smarts):
    """Find aryl components in a molecule."""
    matches = mol.GetSubstructMatches(aryl_smarts)
    subset = [y for x in matches for y in x]
    if len(subset)==0:
        return 0, [], 0, []
    components = connected_components(mol, subset)
    chain_lengths = [len(x) for x in components]
    # Convert components to the same format as other functions return
    matched_chains = [{'chain': list(comp), 'length': len(comp), 'SMARTS': 'cyclic'} for comp in components]
    return max([0]+[len(x) for x in components]), chain_lengths, len(components), matched_chains


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
def parse_groups_in_mol(mol, bycomponent=False, **kwargs):
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
        List of (PFASGroup, match_count, chain_lengths, matched_chains) tuples where:
        - PFASGroup: The matched PFAS group object
        - match_count: Number of times this group pattern was matched (int)
        - chain_lengths: List of carbon chain lengths found (list of int)
        - matched_chains: Detailed information about matched chains (list of dicts)
    
    Notes
    -----
    For PFASgroups 
    1. with smartsPath = 'cyclic', search for connected component matching smarts1
    2. with smarts1 and smarts2 defined: Find shortest path between each pairs of matches
    3. with smarts1 defined and smarts2 null, no formula constraints: either a. (if bycomponent is False) use default smarts2, do like 2. for each default pathtype. b. (if bycomponent is True) search for connected components of fluorinated atoms (for each pathType) where smarts1 C match is in the component
    4. with smarts1 defined and formula constraints: search for substructure matches of smarts1, and for given smartsPath for the PFASgroup, search for connected components of fluorinated atoms where smarts1 C match is in the component

    """
    pfas_groups = kwargs.get('pfas_groups')
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
    fluorinated_components_dict = get_fluorinated_subgraph(mol, **kwargs)
    for pf in pfas_groups:
        #logger.debug(f"{pf.name}")
        matched1_len = 0
        for fd,mol in zip(formulas,frags):
            if pf.formula_dict_satisfies_constraints(fd) is True:
                chain_lengths = []
                if pf.smartsPath =='cyclic':
                    # treat cyclic groups separately
                    match_count, chain_lengths, matched1_len, matched_chains = find_aryl_components(mol, pf.smarts1)
                elif pf.smarts2 is not None and pf.smarts1 is not None:
                    # treat cases with both smarts1 and smarts2, finding path
                    try:
                        match_count, chain_lengths, matched1_len, matched_chains = path_between_smarts(mol,
                                                                            pf.smarts1,
                                                                            pf.smarts2, **kwargs)
                    except ValueError as e:
                        raise e
                elif pf.smarts1 is not None and len(pf.constraints.keys())==0:
                    # treat cases with only smarts1, using default smarts2 (if looking for chains) else 
                    # if bycomponent is True, find fluorinated components for each path type
                    if bycomponent is True:
                        match_count, chain_lengths, matched1_len, matched_chains = find_alkyl_components(mol, pf.smarts1, path_components)
                    else:
                        
                        try:
                            match_count, chain_lengths, matched1_len, matched_chains = path_between_smarts(mol,
                                                                                pf.smarts1,
                                                                                None, **kwargs)
                            
                            # For groups with no constraints and smarts2=null, verify that smarts1 matches
                            # are within 1-2 bonds of a fluorinated terminal atom
                            # This prevents matching non-fluorinated functional groups that are >2 bonds
                            # away from the fluorinated region (since chain SMARTS allow $(CF) and $(CCF))
                            if match_count > 0 and matched_chains:
                                smartsPaths = kwargs.get('smartsPaths')
                                # Get smarts1 matches
                                smarts1_matches = get_substruct(mol, pf.smarts1)
                                # Filter matched_chains where smarts1 atoms are ≤2 bonds from an END atom
                                valid_chains = []
                                for chain_info in matched_chains:
                                    path_type = chain_info['SMARTS']
                                    chain_atoms = set(chain_info['chain'])
                                    # Get the end SMARTS for this path type
                                    if smartsPaths and path_type in smartsPaths:
                                        end_smarts = smartsPaths[path_type][1]  # [chain, end]
                                        end_matches = get_substruct(mol, end_smarts)
                                        # Get smarts1 atoms that are in this chain
                                        smarts1_in_chain = chain_atoms.intersection(smarts1_matches)
                                        # Check if any smarts1 atom is within 2 bonds of an END atom
                                        chain_is_valid = False
                                        for s1_atom in smarts1_in_chain:
                                            # Check if smarts1 atom itself is an END
                                            if s1_atom in end_matches:
                                                chain_is_valid = True
                                                break
                                            # Otherwise, check the actual molecular bond distance to nearest END atom
                                            # (not position in chain array, since path can jump around in molecule)
                                            G = mol_to_nx(mol)
                                            for end_atom in end_matches:
                                                if end_atom in chain_atoms:
                                                    try:
                                                        # Get shortest path distance in molecular graph
                                                        # Uses Breadth-First Search by default (unweighted graph)
                                                        path = nx.shortest_path(G, s1_atom, end_atom)
                                                        distance = len(path) - 1  # Number of bonds
                                                        if distance <= pf.max_dist_from_CF:
                                                            chain_is_valid = True
                                                            break
                                                    except nx.NetworkXNoPath:
                                                        continue
                                            if chain_is_valid:
                                                break
                                        
                                        if chain_is_valid:
                                            valid_chains.append(chain_info)
                                
                                # Update results with only valid chains
                                if valid_chains:
                                    matched_chains = valid_chains
                                    match_count = len(valid_chains)
                                    chain_lengths = [c['length'] for c in valid_chains]
                                    matched1_len = len(set([atom for c in valid_chains for atom in c['chain']]).intersection(smarts1_matches))
                                else:
                                    # No valid chains - smarts1 matches are >2 bonds from fluorinated terminals
                                    matched_chains = []
                                    match_count = 0
                                    chain_lengths = []
                                    matched1_len = 0
                                
                        except ValueError as e:
                            raise e
                elif pf.smarts1 is not None:
                    # treat cases with only smarts1 as simple substructure match, no path finding
                    # if bycomponent is True, find fluorinated components for relevant path type
                    #  (given by pf parameter smartsPath)
                    try:
                        n = len(mol.GetSubstructMatch(pf.smarts1))
                    except:
                        try:
                            Chem.SanitizeMol(mol)
                        except Chem.AtomValenceException:
                            frags = fragment_until_valence_is_correct(mol, [])
                        n = 0
                        for frag in frags:
                            n += len(frag.GetSubstructMatches(pf.smarts1))
                    path_components = fluorinated_components_dict.get(pf.smartsPath,"Perfluoroalkyl")
                    match_count, chain_lengths, matched1_len, matched_chains = find_alkyl_components(mol, pf.smarts1, path_components)
                else:
                    # teat cases with no smarts1, just formula constraints
                    match_count = 1
                    matched_chains = []
                    matched1_len = 1
                if match_count > 0 and matched1_len > 0:
                    # add to matches if functional group was found
                    group_matches.append((pf, match_count, chain_lengths, matched_chains))
    
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

def parse_mols(mols, bycomponent=False, output_format='list', include_PFAS_definitions=True, **kwargs):
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
        formula = CalcMolFormula(mol)
        matches = parse_groups_in_mol(mol, formula=formula, bycomponent=bycomponent, **kwargs)
        inchikey = Chem.MolToInchiKey(mol)
        inchi = Chem.MolToInchi(mol)
        smi = Chem.MolToSmiles(mol)
        results.setdefault(inchikey,{}).update({"smiles": smi,
                        "inchikey": inchikey,
                        "inchi": inchi,
                        "formula": formula})
        results[inchikey].setdefault('matches',[]).extend([{'match_id': f"G{group.id}",
                        'id': group.id,
                        'group_name': group.name,
                        'match_count': match_count,
                        'chain_lengths': chain_lengths,
                        'num_chains': len(matched_chains),
                        'chains': matched_chains,
                        'chain_types': [x for x in set([c['SMARTS'] for c in matched_chains])],
                        'type':'PFASgroup'} for group, match_count, chain_lengths, matched_chains in matches])
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
                        'chain_lengths': match['chain_lengths'],
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
                             include_PFAS_definitions=False,
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
        - 'max_chain': Maximum chain length for groups with chains
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
            all_matches = parse_groups_in_mol(mol, formula=formula, pfas_groups=pfas_groups,include_PFAS_definitions=include_PFAS_definitions, **kwargs)
            
            # Create mapping from group ID to match information
            match_dict = {}
            for group, match_count, chain_lengths, matched_chains in all_matches:
                match_dict[group.id] = {
                    'group': group,
                    'match_count': match_count,
                    'chain_lengths': chain_lengths,
                    'matched_chains': matched_chains
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
                        elif count_mode == 'max_chain':
                            if match_info['matched_chains']:
                                fingerprint[i] = max([chain['length'] for chain in match_info['matched_chains']])
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
                        elif count_mode == 'max_chain':
                            if match_info['matched_chains']:
                                fingerprint[group.name] = max([chain['length'] for chain in match_info['matched_chains']])
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
                            value = match_info['n_matches']
                        elif count_mode == 'max_chain':
                            if match_info['chains']:
                                value = max([chain['length'] for chain in match_info['chains']])
                            else:
                                value = match_info['n_matches'] if match_info['n_matches'] > 0 else 0
                        
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
                            'n_matches': match_info['n_matches'],
                            'n_cfchains': match_info['n_cfchains'],
                            'chains': match_info['chains'],
                            'smarts1': Chem.MolToSmarts(group.smarts1) if group.smarts1 else None,
                            'smarts2': Chem.MolToSmarts(group.smarts2) if group.smarts2 else None
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
    for i, s in enumerate(paths):
        if isinstance(s, int):
            paths[i] = PATH_NAMES[s]
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