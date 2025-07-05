"""
PFASgroups core functions.

This module provides functions for parsing and plotting PFAS groups.
"""

import numpy as np
import networkx as nx
from rdkit import Chem
from rdkit.Chem.rdMolDescriptors import CalcMolFormula
from typing import Union, List, Dict
from PIL import Image
from io import BytesIO
import re
import json
import os

PATH_NAMES = ['Perfluoroalkyl','Polyfluoroalkyl','Polyfluoro','Polyfluorobr']

class PFASGroup:
    def __init__(self, id, name, smarts1, smarts2, smartsPath, constraints, **kwargs):
        self.id = id
        self.name = name
        self.smarts1 = smarts1
        self.smarts2 = smarts2
        self.smartsPath = smartsPath
        if self.smarts1 != "" and self.smarts1 is not None:
            self.smarts1 = Chem.MolFromSmarts(self.smarts1)
            self.smarts1.UpdatePropertyCache()
            Chem.GetSymmSSSR(self.smarts1)
            self.smarts1.GetRingInfo().NumRings()
            if self.smarts2 != "" and self.smarts2 is not None:
                self.smarts2 = Chem.MolFromSmarts(self.smarts2)
                self.smarts2.UpdatePropertyCache()
                Chem.GetSymmSSSR(self.smarts2)
                self.smarts2.GetRingInfo().NumRings()
            else:
                self.smarts2 = None
        else:
            self.smarts1 = None
            self.smarts2 = None
        self.constraints = constraints
    def __str__(self):
        return self.name
    def constraint_gte(self, formula_dict):
        success = True
        for e, n in self.constraints.get('gte', {}).items():
            success = success and formula_dict.get(e, 0) >= n
        return success
    def constraint_lte(self, formula_dict):
        success = True
        for e, n in self.constraints.get('lte', {}).items():
            success = success and formula_dict.get(e, 0) <= n
        return success
    def constraint_eq(self, formula_dict):
        success = True
        for e, n in self.constraints.get('eq', {}).items():
            success = success and formula_dict.get(e, 0) == n
        return success
    def constraint_only(self, formula_dict):
        success = True
        if 'only' in self.constraints.keys():
            tot = sum(formula_dict.values())
            nn = 0
            for e in self.constraints['only']:
                nn += formula_dict.get(e, 0)
            success = success and tot == nn
        return success
    def constraint_rel(self, formula_dict):
        success = True
        for e, v in self.constraints.get('rel', {}).items():
            n = sum([formula_dict.get(x, 0) for x in v.get('atoms', [])])
            success = success and formula_dict.get(e, 0) == n / v.get('div', 1) + v.get('add', 0) + sum([formula_dict.get(x, 0) for x in v.get('add_atoms', [])])
        return success
    def formula_dict_satisfies_constraints(self, formula_dict):
        if len(self.constraints.keys()) == 0:
            return True
        success = True
        process = [None, self.constraint_rel, self.constraint_only, self.constraint_eq, self.constraint_lte, self.constraint_gte]
        k = process.pop()
        while success and k is not None:
            success = k(formula_dict)
            k = process.pop()
        return success

# --- Load SMARTS paths from fpaths.json ---
MODULE_DIR = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = os.path.join(MODULE_DIR, 'data')
FPATHS_FILE = os.path.join(DATA_DIR, 'fpaths.json')
PFAS_GROUPS_FILE = os.path.join(DATA_DIR, 'PFAS_groups_smarts.json')

with open(FPATHS_FILE, 'r') as f:
    SMARTS_PATHS = json.load(f)

def get_smarts_path(name):
    """Return Chem.MolFromSmarts for a named path from fpaths.json."""
    s = SMARTS_PATHS.get(name)
    if s:
        return Chem.MolFromSmarts(s)
    return None

# --- Load PFAS groups from PFAS_groups_smarts.json ---
def load_pfas_groups():
    """Load PFAS groups from the JSON file and return a list of PFASGroup objects."""
    with open(PFAS_GROUPS_FILE, 'r') as f:
        group_data = json.load(f)
    groups = []
    for g in group_data:
        smarts1 = g.get('smarts1')
        smarts2 = g.get('smarts2')
        smartsPath = g.get('smartsPath')
        constraints = g.get('constraints', {})
        # If smartsPath is set, use fpaths.json
        if smartsPath:
            smartsPathMol = get_smarts_path(smartsPath)
        else:
            smartsPathMol = None
        groups.append(PFASGroup(
            name=g.get('name'),
            smarts1=smarts1,
            smarts2=smarts2,
            constraints=constraints
        ))
    return groups

PFAS_GROUPS = load_pfas_groups()

# --- Helper functions ---

def n_from_formula(formula: str, element=None) -> Union[int, dict]:
    """
    Compute the number of elements (any or one specific) in a formula string.
    """
    if element is not None:
        PAT = f"([{element}])"+r"(\d*)"
    else:
        PAT = r"([A-Z][a-z]?)(\d*)"
    mat = re.findall(PAT, formula)
    formula_dict = {}
    for sym, nb in mat:
        if nb != '':
            formula_dict[sym] = formula_dict.setdefault(sym, 0) + int(nb)
        else:
            formula_dict[sym] = formula_dict.setdefault(sym, 0) + 1
    if element is not None:
        return formula_dict[element]
    return formula_dict

def dry_mol_to_nx(mol):
    """
    Construct a networkx graph from a molecule.
    """
    G = nx.Graph()
    for n, atom in enumerate(mol.GetAtoms()):
        element_Z = atom.GetAtomicNum()
        node_params = {"element": element_Z, "symbol": atom.GetSymbol()}
        G.add_node(atom.GetIdx(), **node_params)
    for bond in mol.GetBonds():
        edgeOrder = bond.GetBondTypeAsDouble()
        G.add_edge(bond.GetBeginAtomIdx(), bond.GetEndAtomIdx(), order=edgeOrder)
    return G

def get_substruct(_mol: Chem.Mol, struct: Chem.Mol):
    """
    Returns the indices of the atoms in the molecule that match the substructure.
    """
    return set([x[0] for x in _mol.GetSubstructMatches(struct)])

# --- Main PFAS group parsing functions ---

def parse_PFAS_groups(mol, formula, pfas_groups: List[PFASGroup] = None):
    """
    Iterates over PFAS groups and finds the ones that match the molecule.
    Returns a list of (PFASGroup, n_matches, match_indices).
    """
    mol = Chem.AddHs(mol)
    try:
        Chem.SanitizeMol(mol)
    except Chem.AtomValenceException:
        # Fallback: just return empty for now
        return []
    else:
        frags = [mol]
        formulas = [n_from_formula(formula)]
    if pfas_groups is None:
        pfas_groups = PFAS_GROUPS
    group_matches = []
    for pf in pfas_groups:
        for fd, frag in zip(formulas, frags):
            if pf.formula_dict_satisfies_constraints(fd):
                if pf.smarts1 is not None:
                    matches = frag.GetSubstructMatches(pf.smarts1)
                    n = len(matches)
                    if n > 0:
                        group_matches.append((pf, n, matches))
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
    if paths is None:
        paths = [0, 1, 2, 3]
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
