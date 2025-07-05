"""
Molecule generation utilities for PFASgroups.

This module provides functions to generate random PFAS-like molecules and attach functional groups.
"""

import numpy as np
from rdkit import Chem
from rdkit.Chem.rdMolDescriptors import CalcMolFormula
from rdkit.Chem.Descriptors import ExactMolWt
from typing import Union
from itertools import product, groupby

np.random.seed(2025)  # For reproducibility

def generate_random_carbon_chain(n, cycle=False, alkene=False, alkyne=False):
    """
    Generate a random carbon chain molecule with n carbons.
    """
    if cycle:
        _smiles = 'c1ccccc1'
    else:
        _smiles = 'C'
    m = Chem.MolFromSmiles(_smiles)
    m = Chem.RemoveHs(m)
    rwm = Chem.RWMol(m)
    for _ in range(n - 1):
        i = rwm.AddAtom(Chem.Atom(6))
        j = None
        while j is None:
            j = np.random.randint(0, len(rwm.GetAtoms())-1)
            if j == i:
                j = None
            elif rwm.GetAtoms()[j].GetExplicitValence() == 4:
                j = None
        bondtype = Chem.BondType.SINGLE
        if alkene and np.random.uniform(0,1) < 0.5 and rwm.GetAtoms()[j].GetExplicitValence() < 3:
            bondtype = Chem.BondType.DOUBLE
        elif alkyne and bondtype != Chem.BondType.DOUBLE and np.random.uniform(0,1) < 0.5 and rwm.GetAtoms()[j].GetExplicitValence() < 2:
            bondtype = Chem.BondType.TRIPLE
        rwm.AddBond(i, j, bondtype)
        Chem.SanitizeMol(rwm)
    m2 = rwm.GetMol()
    Chem.SanitizeMol(m2)
    return m2

def fluorinate_mol(mol, perfluorinated=True, p=0.3, phigh=1):
    """
    Fluorinate a molecule by replacing hydrogens with fluorines.
    """
    mol = Chem.AddHs(mol)
    rwm = Chem.RWMol(mol)
    nF = 0
    atomsH_index = [x.GetIdx() for x in rwm.GetAtoms() if x.GetSymbol() == 'H']
    nH = len(atomsH_index)
    np.random.shuffle(atomsH_index)
    for i, atom_index in enumerate(atomsH_index):
        atom = rwm.GetAtomWithIdx(atom_index)
        if perfluorinated or np.random.uniform(0, 1) < max(p, min(phigh, 1-nF/nH)):
            atom.SetAtomicNum(9)  # Replace H with F
            nF += 1
        else:
            atom.SetAtomicNum(1)  # Keep H
    Chem.SanitizeMol(rwm)
    return rwm.GetMol()

def generate_random_mol(n, functional_groups, perfluorinated=True, cycle=False, alkene=False, alkyne=False, **kwargs):
    """
    Generate a random molecule with n carbon atoms and attach functional groups.
    """
    if isinstance(functional_groups, dict):
        functional_groups = [functional_groups]
    elif isinstance(functional_groups, str):
        functional_groups = [{'group_smiles': functional_groups, 'n': kwargs.get('m', 1), 'mode': kwargs.get('mode', 'attach')}]
    elif isinstance(functional_groups, list) and isinstance(functional_groups[0], str):
        functional_groups = [{'group_smiles': functional_group[0], 'n': kwargs.get('m', 1), 'mode': kwargs.get('mode', 'attach')} for functional_group in functional_groups]
    mol = generate_random_carbon_chain(n, cycle, alkene, alkyne)
    mol = fluorinate_mol(mol, perfluorinated=perfluorinated)
    # For brevity, functional group attachment is omitted in this stub
    return mol
