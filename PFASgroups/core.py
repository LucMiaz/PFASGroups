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

PATH_NAMES = ['Perfluoroalkyl','Polyfluoroalkyl']

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
def add_smartsPath(filename = FPATHS_FILE):
    """Yields SMARTS for chains"""
    paths = {}
    with open(filename,'r') as f:
        fpaths = json.load(f)
    names = fpaths.keys()
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
def find_components_between_smarts(mol, smarts1, smarts2, component_solver, **kwargs):
    """Find fluorinated components containing both smarts1 and smarts2 matches.
    Returns components that connect the two functional groups.
    If smarts2 is None, returns components containing smarts1 matches.
    
    When max_dist > 0, matches against extended components but returns original + connecting atoms.
    """
    smartsMatches1 = get_substruct(mol, smarts1)
    if len(smartsMatches1) == 0:
        return 0, [], 0, []
    
    # Get all component types
    max_dist = kwargs.get('max_dist', 0)
    components_by_type = {}
    smartsPaths = kwargs.get('smartsPaths', {})
    
    for path_type in smartsPaths.keys():
        components_by_type[path_type] = component_solver.get(path_type, max_dist=max_dist, default=[])
    
    matched_components = {}
    
    if smarts2 is not None:
        # Both smarts1 and smarts2 defined: find components containing both
        smartsMatches2 = get_substruct(mol, smarts2)
        if len(smartsMatches2) == 0:
            return 0, [], 0, []
        
        smarts_matches = smartsMatches1.union(smartsMatches2)
        for path_type, components in components_by_type.items():
            for i, comp in enumerate(components):
                # Check if component contains at least one match from each SMARTS
                has_smarts1 = any(match in comp for match in smartsMatches1)
                has_smarts2 = any(match in comp for match in smartsMatches2)
                
                if has_smarts1 and has_smarts2:
                    # Get augmented component (handles max_dist > 0 automatically)
                    augmented = component_solver.get_augmented_component(
                        path_type, max_dist, i, smarts_matches
                    )
                    matched_components.setdefault(path_type, []).append(augmented)
    else:
        # Only smarts1 defined: find components containing smarts1 matches
        for path_type, components in components_by_type.items():
            for i, comp in enumerate(components):
                if any(match in comp for match in smartsMatches1):
                    # Get augmented component (handles max_dist > 0 automatically)
                    augmented = component_solver.get_augmented_component(
                        path_type, max_dist, i, smartsMatches1
                    )
                    matched_components.setdefault(path_type, []).append(augmented)
    
    if len(matched_components.keys()) == 0:
        n = 0
        component_lengths = []
    else:
        n = max([0] + [len(comp) for comps in matched_components.values() for comp in comps])
        component_lengths = [len(comp) for comp in matched_components.get('Perfluoroalkyl', [])]
    
    # Get molecular graph for metrics calculation
    G = kwargs.get('G', mol_to_nx(mol))
    pfas_group = kwargs.get('pfas_group', None)
    
    # Add comprehensive metrics to matched components
    matched_components_list = []
    for path_type, comps in matched_components.items():
        for comp in comps:
            matched_components_list.append(
                component_solver.get_matched_component_dict(comp, smartsMatches1, path_type, pfas_group)
            )
    
    return n, component_lengths, len(smartsMatches1), matched_components_list


def components_between_smarts(mol: Chem.Mol, smarts1, smarts2, component_solver, **kwargs):
    """
    Find components containing smarts1 and smarts2 matches.
    
    :param mol: Molecule
    :param smarts1: Primary SMARTS pattern
    :param smarts2: Secondary SMARTS pattern (or None)
    :param component_solver: ComponentsSolver instance
    :return: tuple of (max_component_length, component_lengths, num_smarts1_matches, matched_components)
    """
    return find_components_between_smarts(mol, smarts1, smarts2, component_solver, **kwargs)

@add_smartsPath()
def find_components_with_all_smarts(mol, smarts_list, smarts_counts, component_solver, **kwargs):
    """Find fluorinated components containing all required SMARTS patterns with minimum counts.
    
    This function handles groups that require multiple SMARTS patterns to match simultaneously,
    such as diacids (which need 2 copies of the carboxylic acid SMARTS) or groups with 
    multiple different functional groups.
    
    Parameters
    ----------
    mol : rdkit.Chem.Mol
        Molecule to search
    smarts_list : list of rdkit.Chem.Mol
        List of compiled SMARTS patterns to match
    smarts_counts : tuple of int
        Minimum count required for each SMARTS pattern (e.g., (2,) for diacids, (1, 1) for two different groups)
    component_solver : ComponentsSolver
        Component solver instance
    **kwargs : dict
        Additional parameters including max_dist, smartsPaths, etc.
        
    Returns
    -------
    tuple
        (max_component_length, component_lengths, total_matches, matched_components_list)
    """
    if not smarts_list or len(smarts_list) == 0:
        return 0, [], 0, []
    
    # Get all matches for each SMARTS pattern
    all_matches_by_smarts = []
    for smarts in smarts_list:
        matches = get_substruct(mol, smarts)
        if len(matches) == 0:
            return 0, [], 0, []  # If any required SMARTS is missing, no match
        all_matches_by_smarts.append(matches)
    
    # Combine all matched atoms for component searching
    all_smarts_matches = set()
    for matches in all_matches_by_smarts:
        all_smarts_matches.update(matches)
    
    # Get all component types
    max_dist = kwargs.get('max_dist', 0)
    components_by_type = {}
    smartsPaths = kwargs.get('smartsPaths', {})
    
    for path_type in smartsPaths.keys():
        components_by_type[path_type] = component_solver.get(path_type, max_dist=max_dist, default=[])
    
    matched_components = {}
    
    # Find components that contain sufficient copies of each SMARTS
    for path_type, components in components_by_type.items():
        for i, comp in enumerate(components):
            # Check if this component has enough matches for each SMARTS pattern
            satisfies_all = True
            for smarts_idx, (matches, required_count) in enumerate(zip(all_matches_by_smarts, smarts_counts)):
                # Count how many matches of this SMARTS are in the component
                matches_in_comp = sum(1 for match in matches if match in comp)
                if matches_in_comp < required_count:
                    satisfies_all = False
                    break
            
            if satisfies_all:
                # Get augmented component (handles max_dist > 0 automatically)
                augmented = component_solver.get_augmented_component(
                    path_type, max_dist, i, all_smarts_matches
                )
                matched_components.setdefault(path_type, []).append(augmented)
    
    if len(matched_components.keys()) == 0:
        n = 0
        component_lengths = []
    else:
        n = max([0] + [len(comp) for comps in matched_components.values() for comp in comps])
        component_lengths = [len(comp) for comp in matched_components.get('Perfluoroalkyl', [])]
    
    # Get molecular graph for metrics calculation
    G = kwargs.get('G', mol_to_nx(mol))
    pfas_group = kwargs.get('pfas_group', None)
    
    # Add comprehensive metrics to matched components
    matched_components_list = []
    for path_type, comps in matched_components.items():
        for comp in comps:
            matched_components_list.append(
                component_solver.get_matched_component_dict(comp, all_matches_by_smarts[0], path_type, pfas_group)
            )
    
    return n, component_lengths, len(all_matches_by_smarts[0]), matched_components_list

class ComponentsSolver:
    """Class to hold components information with comprehensive graph metrics."""
    @add_smartsPath()
    def __init__(self, mol, **kwargs):
        self.smartsPaths = kwargs.get('smartsPaths')
        self.mol = mol
        self.mol_size = mol.GetNumAtoms()  # Total atoms in molecule for fraction calculation
        self.total_carbons = sum(1 for atom in mol.GetAtoms() if atom.GetSymbol() == 'C')  # Total carbon atoms
        self.G = mol_to_nx(mol)
        self.components = self.get_fluorinated_subgraph()
        self.extended_components = {k:{0:v} for k,v in self.components.items()}
        # Mapping from (pathType, max_dist, extended_component_index) -> original_component_index
        self.component_to_original_index = {}
        self.levels = {0}
        # Cache for component metrics
        self._component_metrics_cache = {}
        # Precompute full component sizes (including all attached atoms: H, F, Cl, Br, I)
        self.component_full_sizes = {}
        for path_type, components_list in self.components.items():
            if path_type not in self.component_full_sizes:
                self.component_full_sizes[path_type] = {}
            for i, comp in enumerate(components_list):
                full_comp = self.get_full_component_atoms(comp)
                self.component_full_sizes[path_type][i] = len(full_comp)
        # Compute and store metrics for all components on creation
        self._precompute_component_metrics()
        # Initialize mapping for level 0 (original components)
        self._init_component_mapping()
        
    def __enter__(self):
        return self
    def __exit__(self, exc_type, exc_value, traceback):
        self.components = None
        self.extended_components = None
        self.mol = None
        self.G = None
        self._component_metrics_cache = None
        
    def __len__(self):
        return len(self.components)
        
    def get(self, pathType, max_dist=0, default = []):
        if max_dist not in self.levels:
            self.extend_components(max_dist)
        return self.extended_components.get(pathType, {}).get(max_dist, default)
        
    def max_size(self):
        return max([len(x) for x in self.components]) if len(self.components)>0 else 0
        
    def sizes(self):
        return [len(x) for x in self.components]
        
    def _connected_components(self, subset):
        """Find connected components in a molecule."""
        G = self.G.subgraph(subset)
        components = list(nx.connected_components(G))
        return components
        
    def get_fluorinated_subgraph(self, **kwargs):
        """Get the fluorinated indices by connected components of a molecule based on path SMARTS."""
        subsets = {}
        for pathName, d in self.smartsPaths.items():
            path_smarts = d[0] # chain
            matches = self.mol.GetSubstructMatches(path_smarts)
            subset = [y for x in matches for y in x]
            if len(subset)==0:
                subsets[pathName]= []
                continue
            components = self._connected_components(subset)
            subsets[pathName] = components
        return subsets
    
    def get_full_component_atoms(self, component):
        """Get all atoms in a component including those attached (H, F, halogens).
        
        Parameters
        ----------
        component : set
            Set of atom indices representing the carbon backbone
            
        Returns
        -------
        set
            Set of all atom indices including the component and all directly attached atoms
        """
        full_component = set(component)
        # Add all neighbors of component atoms that are not already in the component
        # This includes H, F, and other halogens attached to the carbon backbone
        for atom_idx in component:
            for neighbor_idx in self.G.neighbors(atom_idx):
                if self.mol.GetAtomWithIdx(neighbor_idx).GetSymbol() in ['H', 'F', 'Cl', 'Br', 'I']:
                    full_component.add(neighbor_idx)
        return full_component
    
    def get_total_components_fraction(self, matched_components_list):
        """Calculate the fraction of carbon atoms in the molecule covered by the union of all components.
        
        Parameters
        ----------
        matched_components_list : list of dict
            List of matched component dictionaries, each with 'component' and 'smarts_matches' keys
            
        Returns
        -------
        float
            Fraction of carbon atoms covered by the union of all components (0.0 to 1.0)
        """
        if len(matched_components_list) == 0 or self.total_carbons == 0:
            return 0.0
        
        # Union all carbon atoms from components and SMARTS matches
        union_carbon_atoms = set()
        for comp_dict in matched_components_list:
            component = set(comp_dict.get('component', []))
            smarts_matches = comp_dict.get('smarts_matches')
            
            # Add carbon atoms from component
            for atom_idx in component:
                if self.mol.GetAtomWithIdx(atom_idx).GetSymbol() == 'C':
                    union_carbon_atoms.add(atom_idx)
            
            # Add carbon atoms from SMARTS matches
            if smarts_matches is not None:
                for atom_idx in smarts_matches:
                    if self.mol.GetAtomWithIdx(atom_idx).GetSymbol() == 'C':
                        union_carbon_atoms.add(atom_idx)
        
        # Add extra carbons from functional groups
        # This accounts for carbons in functional groups that are not matched by SMARTS
        # For example, in -COOH, if SMARTS matches the adjacent carbon, we need to add
        # the carbonyl carbon from smarts_extra_atoms
        extra_carbons_count = 0
        for comp_dict in matched_components_list:
            # Get smarts_extra_atoms from pfas_group if available
            if 'smarts_extra_atoms' in comp_dict:
                extra_carbons_count += comp_dict['smarts_extra_atoms']
        
        # Total = union of matched carbons + extra functional group carbons
        # Cap at total_carbons to avoid exceeding 1.0
        total_carbon_count = min(len(union_carbon_atoms) + extra_carbons_count, self.total_carbons)
        total_fraction = total_carbon_count / self.total_carbons
        return total_fraction
    
    def _precompute_component_metrics(self):
        """Precompute metrics for all initial components."""
        for path_type, components_list in self.components.items():
            for comp in components_list:
                # This will populate the cache
                self.compute_component_metrics(comp)
    
    def _init_component_mapping(self):
        """Initialize mapping for original components (max_dist=0)."""
        for pathType in self.components.keys():
            for i in range(len(self.components[pathType])):
                self.component_to_original_index[(pathType, 0, i)] = i
        
    def extend_components(self, max_dist):
        """Extend a component in a graph by a maximum distance.
        This is used to match functional groups that are not directly connected to the component. Different components that overlap are not merged.
        """
        if max_dist>0:
            for pathType, components in self.components.items():
                extended_components = []
                for i, component in enumerate(components):
                    extended = component.copy()
                    for node in component:
                        lengths = nx.single_source_shortest_path_length(self.G, node, cutoff=max_dist)
                        extended.update([n for n,d in lengths.items() if d<=max_dist])
                    extended_components.append(extended)
                    # Map this extended component back to its original component
                    self.component_to_original_index[(pathType, max_dist, i)] = i
                self.extended_components.setdefault(pathType, {})[max_dist] = extended_components
            self.levels.add(max_dist)
    
    def get_augmented_component(self, pathType, max_dist, component_index, smarts_matches):
        """Get original component augmented with connecting atoms to SMARTS matches.
        
        Parameters
        ----------
        pathType : str
            Type of component path (e.g., 'Perfluoroalkyl')
        max_dist : int
            Distance used for extension
        component_index : int
            Index of the component in the extended components list
        smarts_matches : set
            Set of atom indices that matched the SMARTS pattern
            
        Returns
        -------
        set
            Original component augmented with shortest path atoms connecting SMARTS matches
        """
        if max_dist == 0:
            # No augmentation needed, return original component
            return self.components[pathType][component_index]
        
        # Get original and extended components
        orig_index = self.component_to_original_index.get((pathType, max_dist, component_index), component_index)
        orig_comp = self.components[pathType][orig_index]
        ext_comp = self.extended_components[pathType][max_dist][component_index]
        
        # Start with original component
        augmented = set(orig_comp)
        
        # Add shortest paths from SMARTS matches to original component
        for smarts_atom in smarts_matches:
            if smarts_atom in ext_comp and smarts_atom not in orig_comp:
                min_path = None
                min_len = float('inf')
                for base_atom in orig_comp:
                    try:
                        path = nx.shortest_path(self.G, smarts_atom, base_atom)
                        if len(path) < min_len:
                            min_path = path
                            min_len = len(path)
                    except nx.NetworkXNoPath:
                        continue
                if min_path:
                    augmented.update(min_path)
            elif smarts_atom in orig_comp:
                # SMARTS atom already in original component
                augmented.add(smarts_atom)
        
        return augmented
    
    def compute_component_metrics(self, component):
        """Compute comprehensive graph metrics for a component.
        
        Parameters
        ----------
        component : set or frozenset
            Set of atom indices in the component
            
        Returns
        -------
        dict
            Dictionary with graph metrics including:
            - diameter: maximum eccentricity
            - radius: minimum eccentricity
            - eccentricity_values: dict mapping node to its eccentricity
            - center: nodes with minimum eccentricity
            - periphery: nodes with maximum eccentricity
            - barycenter: nodes minimizing sum of distances
            - effective_graph_resistance: sum of resistance distances
        """
        # Use frozenset for caching
        comp_key = frozenset(component)
        if comp_key in self._component_metrics_cache:
            return self._component_metrics_cache[comp_key]
        
        if len(component) <= 1:
            metrics = {
                'diameter': 0,
                'radius': 0,
                'eccentricity_values': {list(component)[0]: 0} if len(component) == 1 else {},
                'center': list(component),
                'periphery': list(component),
                'barycenter': list(component),
                'effective_graph_resistance': 0.0
            }
            self._component_metrics_cache[comp_key] = metrics
            return metrics
        
        # Create subgraph for this component
        subG = self.G.subgraph(component)
        
        # Check if connected
        if not nx.is_connected(subG):
            # For disconnected components, compute metrics separately
            metrics = {
                'diameter': float('inf'),
                'radius': 0,
                'eccentricity_values': {},
                'center': [],
                'periphery': [],
                'barycenter': [],
                'effective_graph_resistance': float('inf')
            }
            self._component_metrics_cache[comp_key] = metrics
            return metrics
        
        try:
            # Compute eccentricity for each node
            eccentricity_values = nx.eccentricity(subG)
            
            # Diameter and radius
            diameter = nx.diameter(subG)
            radius = nx.radius(subG)
            
            # Center and periphery
            center = nx.center(subG)
            periphery = nx.periphery(subG)
            
            # Barycenter: nodes that minimize total distance to all other nodes
            total_distances = {}
            for node in subG.nodes():
                lengths = nx.single_source_shortest_path_length(subG, node)
                total_distances[node] = sum(lengths.values())
            
            min_total_dist = min(total_distances.values())
            barycenter = [node for node, dist in total_distances.items() if dist == min_total_dist]
            
            # Effective graph resistance (requires matrix operations)
            try:
                # Compute resistance distance for small graphs
                if len(component) < 100:  # Limit for performance
                    resistance_sum = 0.0
                    nodes = list(subG.nodes())
                    for i, u in enumerate(nodes):
                        for v in nodes[i+1:]:
                            try:
                                # Resistance distance approximation using shortest path
                                # For more accurate computation, would need Laplacian pseudoinverse
                                sp_length = nx.shortest_path_length(subG, u, v)
                                resistance_sum += sp_length
                            except:
                                resistance_sum += float('inf')
                    effective_graph_resistance = resistance_sum
                else:
                    effective_graph_resistance = float('nan')
            except:
                effective_graph_resistance = float('nan')
            
            metrics = {
                'diameter': diameter,
                'radius': radius,
                'eccentricity_values': eccentricity_values,
                'center': center,
                'periphery': periphery,
                'barycenter': barycenter,
                'effective_graph_resistance': effective_graph_resistance
            }
            
        except Exception as e:
            # Fallback for any computation errors
            metrics = {
                'diameter': float('nan'),
                'radius': float('nan'),
                'eccentricity_values': {},
                'center': [],
                'periphery': [],
                'barycenter': [],
                'effective_graph_resistance': float('nan')
            }
        
        self._component_metrics_cache[comp_key] = metrics
        return metrics
    
    def compute_smarts_component_metrics(self, component, smarts_matches):
        """Compute metrics relating SMARTS matches to component structural features.
        
        Parameters
        ----------
        component : set
            Set of atom indices in the component
        smarts_matches : set
            Set of atom indices matching the SMARTS pattern
            
        Returns
        -------
        dict
            Dictionary with SMARTS-specific metrics, or None if no smarts_matches
        """
        if smarts_matches is None or len(smarts_matches) == 0:
            return None
            
        comp_metrics = self.compute_component_metrics(component)
        smarts_in_comp = smarts_matches.intersection(component)
        
        if len(smarts_in_comp) == 0 or len(component) <= 1:
            return {
                'min_dist_to_barycenter': 0,
                'min_resistance_dist_to_barycenter': 0.0,
                'min_dist_to_center': 0,
                'min_resistance_dist_to_center': 0.0,
                'max_dist_to_periphery': 0,
                'max_resistance_dist_to_periphery': 0.0
            }
        
        subG = self.G.subgraph(component)
        
        # Initialize metrics
        min_dist_to_barycenter = float('inf')
        min_resistance_to_barycenter = float('inf')
        min_dist_to_center = float('inf')
        min_resistance_to_center = float('inf')
        max_dist_to_periphery = 0
        max_resistance_to_periphery = 0.0
        
        try:
            # Compute distances from SMARTS matches to structural features
            for smarts_node in smarts_in_comp:
                if smarts_node not in subG:
                    continue
                    
                # Distances to barycenter nodes
                for bc_node in comp_metrics['barycenter']:
                    try:
                        dist = nx.shortest_path_length(subG, smarts_node, bc_node)
                        min_dist_to_barycenter = min(min_dist_to_barycenter, dist)
                        # Resistance distance approximation
                        min_resistance_to_barycenter = min(min_resistance_to_barycenter, float(dist))
                    except:
                        pass
                
                # Distances to center nodes
                for center_node in comp_metrics['center']:
                    try:
                        dist = nx.shortest_path_length(subG, smarts_node, center_node)
                        min_dist_to_center = min(min_dist_to_center, dist)
                        min_resistance_to_center = min(min_resistance_to_center, float(dist))
                    except:
                        pass
                
                # Distances to periphery nodes
                for periph_node in comp_metrics['periphery']:
                    try:
                        dist = nx.shortest_path_length(subG, smarts_node, periph_node)
                        max_dist_to_periphery = max(max_dist_to_periphery, dist)
                        max_resistance_to_periphery = max(max_resistance_to_periphery, float(dist))
                    except:
                        pass
            
            # Handle cases where no valid distances were found
            if min_dist_to_barycenter == float('inf'):
                min_dist_to_barycenter = 0
                min_resistance_to_barycenter = 0.0
            if min_dist_to_center == float('inf'):
                min_dist_to_center = 0
                min_resistance_to_center = 0.0
                
        except Exception as e:
            # Fallback values
            min_dist_to_barycenter = 0
            min_resistance_to_barycenter = 0.0
            min_dist_to_center = 0
            min_resistance_to_center = 0.0
            max_dist_to_periphery = 0
            max_resistance_to_periphery = 0.0
        
        return {
            'min_dist_to_barycenter': min_dist_to_barycenter,
            'min_resistance_dist_to_barycenter': min_resistance_to_barycenter,
            'min_dist_to_center': min_dist_to_center,
            'min_resistance_dist_to_center': min_resistance_to_center,
            'max_dist_to_periphery': max_dist_to_periphery,
            'max_resistance_dist_to_periphery': max_resistance_to_periphery
        }
    
    def get_matched_component_dict(self, component, smarts_matches=None, smarts_type='unknown', pfas_group=None):
        """Get a dictionary with all metrics for a matched component.
        
        Parameters
        ----------
        component : set
            Set of atom indices in the component
        smarts_matches : set or None
            Set of atom indices matching the SMARTS pattern (None if no SMARTS)
        smarts_type : str
            Type identifier for the SMARTS pattern
        pfas_group : PFASGroup or None
            PFASGroup object with precomputed SMARTS atom counts
            
        Returns
        -------
        dict
            Complete dictionary with all component metrics
        """
        # Basic branching metric
        smarts_set = smarts_matches if smarts_matches is not None else set()
        basic_metrics = calculate_component_metrics(self.G, component, smarts_set)
        
        # Comprehensive graph metrics (cached)
        comp_metrics = self.compute_component_metrics(component)
        
        # SMARTS-specific metrics (computed on the fly, None if no smarts)
        smarts_metrics = self.compute_smarts_component_metrics(component, smarts_matches)
        
        # Get precomputed SMARTS extra atoms count if pfas_group is available
        smarts_extra_atoms = 0
        if pfas_group is not None and smarts_matches is not None and len(smarts_matches) > 0:
            if pfas_group.smarts_extra_atoms is not None:
                # Sum the extra atoms from all SMARTS patterns
                # For groups with multiple matches, we count each match
                smarts_extra_atoms = sum(pfas_group.smarts_extra_atoms) * len(smarts_matches)
        
        # Calculate mean and median eccentricity from eccentricity_values
        eccentricity_values = comp_metrics.get('eccentricity_values', {})
        if len(eccentricity_values) > 0:
            ecc_list = list(eccentricity_values.values())
            mean_eccentricity = sum(ecc_list) / len(ecc_list)
            sorted_ecc = sorted(ecc_list)
            n = len(sorted_ecc)
            if n % 2 == 0:
                median_eccentricity = (sorted_ecc[n//2 - 1] + sorted_ecc[n//2]) / 2.0
            else:
                median_eccentricity = sorted_ecc[n//2]
        else:
            mean_eccentricity = 0.0
            median_eccentricity = 0.0
        
        # Calculate component fraction based on carbon atoms only
        # Count carbon atoms in component
        component_carbons = sum(1 for atom_idx in component if self.mol.GetAtomWithIdx(atom_idx).GetSymbol() == 'C')
        
        # Count carbon atoms in SMARTS matches that are NOT already in the component
        smarts_carbons_not_in_component = 0
        if smarts_matches is not None:
            for atom_idx in smarts_matches:
                if self.mol.GetAtomWithIdx(atom_idx).GetSymbol() == 'C' and atom_idx not in component:
                    smarts_carbons_not_in_component += 1
        
        # smarts_extra_atoms represents additional carbons in the functional group beyond the matched carbon
        # For each SMARTS pattern, extra_atoms indicates carbons beyond the primary matched atom
        # These are summed across all SMARTS matches to get the total extra carbons
        # So we always add smarts_extra_atoms regardless of whether matched carbon is in component
        
        total_carbons_in_component = component_carbons + smarts_carbons_not_in_component + smarts_extra_atoms
        
        component_fraction = total_carbons_in_component / self.total_carbons if self.total_carbons > 0 else 0.0
        
        result = {
            'component': list(component), 
            'size': len(component), 
            'component_fraction': component_fraction,  # Fraction of molecule covered by component
            'smarts_matches': list(smarts_matches) if smarts_matches is not None else None,  # Store for union calculation
            'smarts_extra_atoms': smarts_extra_atoms,  # Extra carbons from functional group
            'SMARTS': smarts_type,
            # Basic metrics
            'branching': basic_metrics['branching'],
            'smarts_centrality': basic_metrics['smarts_centrality'],
            # Graph structure metrics
            'diameter': comp_metrics.get('diameter', float('nan')),
            'radius': comp_metrics.get('radius', float('nan')),
            'effective_graph_resistance': comp_metrics.get('effective_graph_resistance', float('nan')),
            'eccentricity_values': comp_metrics.get('eccentricity_values', {}),
            'mean_eccentricity': mean_eccentricity,
            'median_eccentricity': median_eccentricity,
            'center': comp_metrics.get('center', []),
            'periphery': comp_metrics.get('periphery', []),
            'barycenter': comp_metrics.get('barycenter', []),
            # Distance metrics (with defaults)
            'min_dist_to_barycenter': 0,
            'min_resistance_dist_to_barycenter': 0.0,
            'min_dist_to_center': 0,
            'min_resistance_dist_to_center': 0.0,
            'max_dist_to_periphery': 0,
            'max_resistance_dist_to_periphery': 0.0
        }
        
        # Override distance metrics with SMARTS-specific values if available
        if smarts_metrics is not None:
            result.update({
                'min_dist_to_barycenter': smarts_metrics.get('min_dist_to_barycenter', 0),
                'min_resistance_dist_to_barycenter': smarts_metrics.get('min_resistance_dist_to_barycenter', 0.0),
                'min_dist_to_center': smarts_metrics.get('min_dist_to_center', 0),
                'min_resistance_dist_to_center': smarts_metrics.get('min_resistance_dist_to_center', 0.0),
                'max_dist_to_periphery': smarts_metrics.get('max_dist_to_periphery', 0),
                'max_resistance_dist_to_periphery': smarts_metrics.get('max_resistance_dist_to_periphery', 0.0)
            })
        
        return result

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

def calculate_component_metrics(G, component, smarts_matches):
    """Calculate branching and centrality metrics for a component.
    
    Parameters
    ----------
    G : networkx.Graph
        Molecular graph
    component : set
        Set of atom indices in the component
    smarts_matches : set
        Set of atom indices matching the SMARTS pattern
    
    Returns
    -------
    dict
        Dictionary with 'branching' (float) and 'smarts_centrality' (float)
    """
    if len(component) <= 1:
        return {'branching': 0.0, 'smarts_centrality': 1.0}
    
    # Create subgraph for this component
    subG = G.subgraph(component)
    
    # Calculate branching: measure of branching vs linearity
    # For linear chains: branching → 1.0
    # For highly branched: branching → 0.0
    try:
        # Count branch points (degree > 2)
        branch_points = sum(1 for node in subG.nodes() if subG.degree(node) > 2)
        # Normalize by component size
        branching = 1.0 - (branch_points / max(1, len(component) - 2))  # -2 to account for endpoints
        branching = max(0.0, min(1.0, branching))  # Clamp to [0, 1]
    except:
        branching = 0.0
    
    # Calculate SMARTS centrality: how central the matched atoms are
    smarts_in_component = smarts_matches.intersection(component)
    if len(smarts_in_component) == 0:
        smarts_centrality = 0.0
    else:
        try:
            # Calculate average shortest path distance from SMARTS matches to all other nodes
            total_distance = 0
            count = 0
            for smarts_node in smarts_in_component:
                if smarts_node in subG:
                    lengths = nx.single_source_shortest_path_length(subG, smarts_node)
                    for node, dist in lengths.items():
                        if node != smarts_node:
                            total_distance += dist
                            count += 1
            
            if count > 0:
                avg_distance = total_distance / count
                # Calculate maximum possible average distance (for peripheral node)
                # For a linear chain of n nodes, max avg distance is ~n/3
                max_possible_distance = len(component) / 3.0
                # Centrality: 1.0 = central, 0.0 = peripheral
                smarts_centrality = 1.0 - min(1.0, avg_distance / max(1.0, max_possible_distance))
            else:
                smarts_centrality = 1.0
        except:
            smarts_centrality = 0.5
    
    return {'branching': branching, 'smarts_centrality': smarts_centrality}

def find_alkyl_components(mol, smarts, components, pathType, component_solver, max_dist=0, **kwargs):
    """Find alkyl components in a molecule with comprehensive metrics.
    
    When max_dist > 0:
    - Matches SMARTS against extended components
    - Returns original component augmented with connecting atoms from SMARTS to original component
    """
    matches = mol.GetSubstructMatches(smarts)
    subset = set(y for x in matches for y in x)
    if len(subset)==0:
        return 0, [], 0, []
    
    if pathType is not None:
        pathTypes = [pathType]
    else:
        pathTypes = component_solver.smartsPaths.keys()
    
    # Filter components connected to the smarts and get augmented versions
    pfas_group = kwargs.get('pfas_group', None)
    matched_components = []
    for _pathType in pathTypes:
        extended_components = component_solver.get(_pathType, max_dist, [])
        for i, comp in enumerate(extended_components):
            # Check if this component is connected to SMARTS matches
            augmented = component_solver.get_augmented_component(
                    _pathType, max_dist, i, subset
                )
            if len(subset.intersection(augmented)) > 0:
                matched_components.append(
                component_solver.get_matched_component_dict(comp, subset, _pathType, pfas_group)
                )
    
    if len(matched_components) == 0:
        return 0, [], 0, []
    
    # Get all component sizes from all path types
    all_components = list(set([comp for comps in matched_components for comp in comps]))
    component_sizes = [len(x) for x in all_components]
    
    # Convert components to dictionary format with comprehensive metrics
    
    
    return max([0] + component_sizes), component_sizes, len(all_components), matched_components

def find_aryl_components(mol, aryl_smarts, component_solver=None, **kwargs):
    """Find aryl components in a molecule with comprehensive metrics."""
    matches = mol.GetSubstructMatches(aryl_smarts)
    subset = [y for x in matches for y in x]
    if len(subset)==0:
        return 0, [], 0, []
    
    components = component_solver._connected_components(subset)
    component_sizes = [len(x) for x in components]
    
    # Get molecular graph for metrics calculation
    subset_set = set(subset)
    pfas_group = kwargs.get('pfas_group', None)
    
    # Convert components to the same format as other functions return with comprehensive metrics
    matched_components = []
    for comp in components:
        matched_components.append(
            component_solver.get_matched_component_dict(comp, subset_set, 'cyclic', pfas_group)
        )
    
    return max([0]+[len(x) for x in components]), component_sizes, len(components), matched_components


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
        kwargs['max_dist'] = pf.max_dist_from_CF if pf.max_dist_from_CF is not None else 0
        for fd,mol in zip(formulas,frags):
            if pf.formula_dict_satisfies_constraints(fd) is True:
                component_sizes = []
                # Create molecular graph once for this fragment
                G = mol_to_nx(mol)
                kwargs['G'] = G
                
                if pf.smartsPath =='cyclic':
                    # treat cyclic groups separately
                    kwargs['pfas_group'] = pf
                    # Use first SMARTS pattern for cyclic
                    match_count, component_sizes, matched1_len, matched_components = find_aryl_components(mol, pf.smarts[0] if pf.smarts else None, component_solver=fluorinated_components_dict, **kwargs)
                elif pf.smarts is not None and len(pf.smarts) > 0:
                    # Handle groups with SMARTS patterns
                    kwargs['pfas_group'] = pf
                    
                    # Check if we need to match multiple SMARTS simultaneously
                    # (e.g., diacids need 2 copies of the same SMARTS)
                    has_multi_count = any(count > 1 for count in pf.smarts_count)
                    has_multiple_smarts = len(pf.smarts) > 1
                    
                    if has_multiple_smarts or has_multi_count:
                        # Multiple SMARTS or multiple copies required
                        # Use components_between_smarts to find components containing all required patterns
                        try:
                            match_count, component_sizes, matched1_len, matched_components = find_components_with_all_smarts(
                                mol, pf.smarts, pf.smarts_count, fluorinated_components_dict, **kwargs
                            )
                        except ValueError as e:
                            raise e
                    else:
                        # Single SMARTS with count=1
                        if pf.smartsPath is None:
                            # If no smartsPath specified, check all available component types
                            all_components = []
                            max_dist_val = pf.max_dist_from_CF if pf.max_dist_from_CF is not None else 0
                            for path_type in kwargs.get('smartsPaths',{'Perfluoroalkyl','Polyfluoroalkyl'}).keys():
                                path_comps = fluorinated_components_dict.get(path_type, max_dist = max_dist_val, default=[])
                                all_components.extend(path_comps)
                            path_components = all_components
                            used_pathType = None  # Mixed types, can't use get_augmented_component properly
                        else:
                            max_dist_val = pf.max_dist_from_CF if pf.max_dist_from_CF is not None else 0
                            path_components = fluorinated_components_dict.get(pf.smartsPath, max_dist = max_dist_val, default = fluorinated_components_dict.get("Polyfluoroalkyl", max_dist = max_dist_val, default=[]))
                            used_pathType = pf.smartsPath
                        kwargs['max_dist'] = pf.max_dist_from_CF if pf.max_dist_from_CF is not None else 0
                        match_count, component_sizes, matched1_len, matched_components = find_alkyl_components(
                            mol, pf.smarts[0], path_components,
                            pathType=used_pathType,
                            component_solver=fluorinated_components_dict,
                            **kwargs)
                elif pf.smartsPath is not None:
                    # treat cases with only smartsPath defined (no SMARTS patterns), find all components of that path type
                    max_dist_val = pf.max_dist_from_CF if pf.max_dist_from_CF is not None else 0
                    path_components = fluorinated_components_dict.get(pf.smartsPath, max_dist = max_dist_val, default = fluorinated_components_dict.get("Polyfluoroalkyl", max_dist = max_dist_val, default=[]))
                    match_count = len(path_components)
                    component_sizes = [len(x) for x in path_components]
                    matched1_len = len(path_components)  # Set matched1_len to enable group matching
                    matched_components = []
                    for comp in path_components:
                        # Use get_matched_component_dict with no SMARTS matches
                        matched_components.append(
                            fluorinated_components_dict.get_matched_component_dict(comp, None, pf.smartsPath, pf)
                        )
                else:
                    # treat cases with no SMARTS patterns, just formula constraints
                    match_count = 1
                    matched_components = []
                    matched1_len = 1
                if match_count > 0 and matched1_len > 0:
                    # add to matches if functional group was found
                    group_matches.append((pf, match_count, component_sizes, matched_components))
    
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