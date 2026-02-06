import os
import json
import networkx as nx
import re
from rdkit import Chem
from rdkit.Chem.rdMolDescriptors import CalcMolFormula
from typing import Union, List, Dict
from rdkit import rdBase

PATH_NAMES = ['Perfluoroalkyl','Polyfluoroalkyl']
# --- Load SMARTS paths from fpaths.json ---
MODULE_DIR = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = os.path.join(MODULE_DIR, 'data')
PFAS_GROUPS_FILE = os.path.join(DATA_DIR, 'PFAS_groups_smarts.json')
PFAS_DEFINITIONS_FILE = os.path.join(DATA_DIR, 'PFAS_definitions_smarts.json')
FPATHS_FILE = os.path.join(DATA_DIR, 'fpaths.json')

def rdkit_disable_log(level='warning'):
    """Disable RDKit warnings and errors logging to stderr"""
    def disable_logs():
        if level == 'error':
            rdBase.DisableLog('rdApp.error')
            rdBase.DisableLog('rdApp.warning')
        elif level == 'warning':
            rdBase.DisableLog('rdApp.warning')
        else:
            rdBase.DisableLog('rdApp.*')
    def enable_logs():
        rdBase.EnableLog('rdApp.*')
    def inner(func):
        def wrapper(*args,**kwargs):
            disable_logs()
            func_ret = func(*args, **kwargs)
            enable_logs()
            return func_ret
        return wrapper
    return inner
    rdBase.DisableLog('rdApp.error')
    rdBase.DisableLog('rdApp.warning')

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



# --- Add SMARTS paths to function ---
def add_componentSmarts(filename = FPATHS_FILE):
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
            kwargs["componentSmartss"] = kwargs.get("componentSmartss",paths)
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
