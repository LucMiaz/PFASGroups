"""
PFASgroups core functions.

This module provides functions for parsing and plotting PFAS groups.
"""

import numpy as np
import networkx as nx
from rdkit import Chem
from rdkit.Chem.rdMolDescriptors import CalcMolFormula
from PFASGroupModel import PFASGroup
from typing import Union, List, Dict
from PIL import Image
from io import BytesIO
import re
import json
import os

PATH_NAMES = ['Perfluoroalkyl','Polyfluoroalkyl','Polyfluoro','Polyfluorobr']

# --- Load SMARTS paths from fpaths.json ---
MODULE_DIR = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = os.path.join(MODULE_DIR, 'data')
FPATHS_FILE = os.path.join(DATA_DIR, 'fpaths.json')
PFAS_GROUPS_FILE = os.path.join(DATA_DIR, 'PFAS_groups_smarts.json')


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

# --- Add SMARTS paths to function ---
def add_smartsPath(names =PATH_NAMES):
    """Yields SMARTS for chains"""
    paths = {}
    pathsEnd = {}
    with open(FPATHS_FILE,'r') as f:
        fpaths = json.load(f) 
    for n in names:
        s = fpaths.get(n)
        e = fpaths.get(n+"_end")
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
            kwargs["smartsPaths"] = paths
            return func(*args, **kwargs)
        return wrapper
    return inner

# --- Helper functions ---


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
def find_path_between_smarts(mol,smarts1,smarts2,G, smartsPaths):
    """Iterates over substructures matched by smarts1 and smarts2 and finds fluorinated paths between them.
    Fluorinated atoms are defined by the SMARTS yield by the decorator add_smartsPath.
    if smarts2 is None, uses SMARTS corresponding to smartsPath."""
    chains = {}
    smartsMatches1 = get_substruct(mol, smarts1)
    if smarts2 is not None:
        pairs = {smarts2:[(path,get_substruct(mol,smartsPath)) for path, (smartsPath,_)  in smartsPaths.items()]}
    else:
        pairs = {_smarts2:[(path,get_substruct(mol, smartsPath))] for path, (smartsPath,_smarts2) in smartsPaths.items()}
    def match_intersection(_pair,_setp):
        for path,pmatch in _pair:
            if _setp == _setp.intersection(pmatch):
                chains.setdefault(path,[]).append(_setp)
                return path
        return None
    for match1 in smartsMatches1:
        for _smarts2, smartsPath in pairs.items():
            smartsMatches2 = get_substruct(mol,_smarts2)
            for match2 in smartsMatches2:
                if (smarts1 != _smarts2) or (match1!=match2):# avoid double matching for diacids
                    path_idx = nx.shortest_path(G, match1, match2, method='dijkstra')
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
        n_CFchain = []
    else:
        n = min([len(x) for x in chains.values()])#{k:len(chain) for k,chain in chains.items()}
        n_CFchain = [len(x) for x in chains.get('Perfluoroalkyl',[[]])]#{k:[len(x) for x in chain] for k,chain in chains.items()}
    chains = [{'chain':list(chain), 'length':len(chain),'SMARTS':path} for path,all_chains in chains.items() for chain in all_chains]
    return n,n_CFchain, len(smartsMatches1), chains  


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
                                    G)

# --- Main PFAS group parsing functions ---
@load_PFASGroups()
def parse_PFAS_groups(mol, formula, pfas_groups=None):
    """Iterates over PFAS groups and finds the ones that match the molecule."""
    mol = Chem.AddHs(mol)
    try:
        Chem.SanitizeMol(mol)
    except Chem.AtomValenceException:
        #logger.debug("failed sanitisation, fragmenting")
        frags = fragment_until_valence_is_correct(mol, [])
        formulas = [n_from_formula(CalcMolFormula(frag)) for frag in frags]
    else:
        frags = [mol]
        formulas = [n_from_formula(formula)]# formula as a dictionary
    group_matches = []
    for pf in pfas_groups:
        #logger.debug(f"{pf.name}")
        matched1_len = 0
        for fd,mol in zip(formulas,frags):
            if pf.formula_dict_satisfies_constraints(fd) is True:
                n_CFchain = []
                if pf.smarts2 is not None and pf.smarts1 is not None:
                    try:
                        n,n_CFchain,matched1_len, chains = path_between_smarts(mol,
                                                                            pf.smarts1,
                                                                            pf.smarts2)
                    except ValueError as e:
                        raise e
                elif pf.smarts1 is not None and len(pf.constraints.keys())==0:
                    try:
                        n,n_CFchain,matched1_len, chains = path_between_smarts(mol,
                                                                            pf.smarts1,
                                                                            None)
                    except ValueError as e:
                        raise e
                elif pf.smarts1 is not None:
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
                    chains = []
                    matched1_len = 1
                else:
                    n=1
                    chains = []
                    matched1_len = 1
                if n>0 and matched1_len>0:
                    group_matches.append((pf,n,n_CFchain, chains))
    return group_matches

def parse_pfas(smiles_list):
    """
    Parse a list of SMILES strings and return PFAS group information.
    """
    results = []
    for smiles in smiles_list:
        mol = Chem.MolFromSmiles(smiles)
        formula = CalcMolFormula(mol)
        matches = parse_PFAS_groups(mol, formula)
        results.append(matches)
    return results

def plot_pfasgroups(smiles: Union[list, str], display=True, path=None, svg=False, ipython=False, subwidth=300, subheight=300, buffer=1, ncols=2, addAtomIndices=True, addBondIndices=False, paths=[0, 1, 2, 3]):
    """
    Plot PFAS group assignments for a list of SMILES strings.
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
        matches = parse_PFAS_groups(mol, CalcMolFormula(mol))
        highlight_atoms = []
        for pf, n, match_indices in matches:
            for match in match_indices:
                highlight_atoms.extend(match)
        new_img = draw_subfig(f"{i}", atoms=highlight_atoms)
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
        imgs = [Image.open(BytesIO(img.encode('utf-8'))) for img in imgs]
    else:
        imgs = [Image.open(BytesIO(img)) for img in imgs]
    # Simple grid
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
