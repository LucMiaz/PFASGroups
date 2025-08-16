"""
Fragmentation utilities for PFASgroups.

This module provides functions for fragmenting molecules based on bond dissociation energies.
"""

import numpy as np
from rdkit import Chem
from rdkit.Chem.rdMolDescriptors import CalcMolFormula
from rdkit.Chem.Descriptors import ExactMolWt
import json
import os

def keysToInt(x):
    """Convert string keys to integers in a dictionary"""
    return {int(k): v for k, v in x.items()}

def yield_scheme(name):
    """Yields SMARTS for chains"""
    if name == 'BDE':
        # For standalone use, we would need to provide the bond dissociation energy data
        # This is a simplified version that creates a basic BDE dictionary
        # In a full implementation, this would load from a data file
        data = {
            6: {6: 346, 1: 414, 9: 485, 17: 339, 35: 285},  # Carbon bonds
            1: {6: 414, 1: 436, 9: 568, 17: 431, 35: 366},  # Hydrogen bonds  
            9: {6: 485, 1: 568, 9: 159},  # Fluorine bonds
            17: {6: 339, 1: 431, 17: 243}, # Chlorine bonds
            35: {6: 285, 1: 366, 35: 193}, # Bromine bonds
        }
        # Make symmetric
        for atom1 in data:
            for atom2 in data[atom1]:
                if atom2 not in data:
                    data[atom2] = {}
                data[atom2][atom1] = data[atom1][atom2]
    else:
        raise ValueError(f"Yield Scheme not implemented for {name}")
    
    def inner(func):
        def wrapper(*args,**kwargs):
            kwargs["bondWeights"] = data
            return func(*args, **kwargs)
        return wrapper
    return inner

def get_fragments(mol, idx, all=False):
    """Fragment a molecule by breaking specified bonds"""
    frags = []
    if all is True:
        with Chem.RWMol(mol) as rwmol:
            for b_idx in idx:
                b = rwmol.GetBondWithIdx(b_idx)
                rwmol.RemoveBond(b.GetBeginAtomIdx(), b.GetEndAtomIdx())
        frags = frags + [x for x in Chem.GetMolFrags(rwmol, asMols=True)]
    else:
        for b_idx in idx:
            with Chem.RWMol(mol) as rwmol:
                b = rwmol.GetBondWithIdx(b_idx)
                rwmol.RemoveBond(b.GetBeginAtomIdx(), b.GetEndAtomIdx())
            frags = frags + [x for x in Chem.GetMolFrags(rwmol, asMols=True)]
    return frags

@yield_scheme('BDE')
def generate_fragments(mol, bondWeights=None, nb_breakingpoints=5):
    """Break molecule on weakest bonds using BDE"""
    bonds = []
    strengths = []
    for i, bond in enumerate(mol.GetBonds()):
        node1 = bond.GetBeginAtom().GetAtomicNum()
        node2 = bond.GetEndAtom().GetAtomicNum()
        try:
            value = bondWeights[node1][node2]
        except:
            try:
                value = bondWeights[node2][node1]
            except:
                value = 300  # Default bond strength
        bonds.append(i)
        strengths.append(value)
    
    # Invert strengths for probability calculation (weaker bonds more likely to break)
    inv_strengths = [1.0/s for s in strengths]
    prob = np.array(inv_strengths) / np.sum(inv_strengths)
    
    fragments = []
    breaking_points_idx = np.random.choice(range(len(bonds)), 
                                         size=min(nb_breakingpoints, len(bonds)-1), 
                                         replace=False, p=prob)
    frags = get_fragments(mol, [bonds[i] for i in breaking_points_idx])
    fragments += frags
    return fragments

def fragment(mol, bde_dict, max_depth=4, nb_breakingpoints=5):
    """Fragment a molecule into smaller fragments by breaking the weakest bonds.
    mol: rdkit.Chem.Mol object
    bde_dict: dictionary with bond dissociation energies"""
    if not isinstance(mol, Chem.Mol):
        raise TypeError("Input mol must be an rdkit.Chem.Mol object")
    
    fragments = []
    process_fragments = [mol]
    
    for i in range(max_depth):
        next_fragments = []
        for frag in process_fragments:
            if len(frag.GetBonds()) > 1:
                frags = generate_fragments(frag, bondWeights=bde_dict, nb_breakingpoints=nb_breakingpoints)
                fragments = fragments + frags
                next_fragments = next_fragments + frags
        process_fragments = next_fragments
    
    fragments_d = {}
    for frag in fragments:
        try:
            inchi_key = Chem.MolToInchiKey(frag)
            formula = CalcMolFormula(frag)
            if inchi_key not in fragments_d:
                fragments_d[inchi_key] = {}
            if formula not in fragments_d[inchi_key]:
                fragments_d[inchi_key][formula] = {'n': 0, 'mol': frag}
            fragments_d[inchi_key][formula]['n'] += 1
        except:
            # Skip problematic fragments
            continue
    
    return fragments_d

def fragment_to_distribution(mol, bde_dict, max_depth=4, nb_breakingpoints=5):
    """Fragment a molecule into smaller fragments and return mass spectrum-like distribution.
    mol: rdkit.Chem.Mol object
    bde_dict: dictionary with bond dissociation energies"""
    if not isinstance(mol, Chem.Mol):
        raise TypeError("Input mol must be an rdkit.Chem.Mol object")
    
    fragments = fragment(mol, bde_dict, max_depth=max_depth, nb_breakingpoints=nb_breakingpoints)
    intensities = {}
    
    for key, value in fragments.items():
        for formula, data in value.items():
            try:
                molw = ExactMolWt(data['mol'])
                if molw not in intensities:
                    intensities[molw] = 0
                intensities[molw] += data['n']
            except:
                continue
    
    mz = list(intensities.keys())
    intensity = list(intensities.values())
    
    if len(intensity) > 0:
        intensity_distribution = np.array(intensity) / np.sum(intensity)
        sorted_indices = np.argsort(mz)
        mz = np.array(mz)[sorted_indices]
        intensity_distribution = intensity_distribution[sorted_indices]
        return mz, intensity_distribution
    else:
        return np.array([]), np.array([])
