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
from itertools import combinations, product
import networkx as nx
from core import mol_to_nx

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


def find_fluorinated_chains(mol, smartsPathName='Perfluoroalkyl'):
    """
    Find fluorinated chains in a molecule using SMARTS patterns.
    Returns atom indices that are part of fluorinated chains.
    """
    # Load SMARTS paths from fpaths.json
    MODULE_DIR = os.path.dirname(os.path.abspath(__file__))
    DATA_DIR = os.path.join(MODULE_DIR, 'data')
    FPATHS_FILE = os.path.join(DATA_DIR, 'fpaths.json')
    
    try:
        with open(FPATHS_FILE, 'r') as f:
            fpaths = json.load(f)
    except FileNotFoundError:
        # Fallback SMARTS patterns if file not found
        fpaths = {
            'Perfluoroalkyl': '[#6]([#9])([#9])[#6]([#9])([#9])',
            'Perfluoroalkyl_end': '[#6]([#9])([#9])([#6,#9,#17,#35,#53])',
            'Polyfluoroalkyl': '[#6]([#9,#1,#17,#35,#53])([#9,#1,#17,#35,#53])[#6]([#9,#1,#17,#35,#53])([#9,#1,#17,#35,#53])',
            'Polyfluoroalkyl_end': '[#6]([#9,#1,#17,#35,#53])([#9,#1,#17,#35,#53])([#6,#9,#1,#17,#35,#53])'
        }
    
    # Get SMARTS patterns
    path_smarts = fpaths.get(smartsPathName)
    end_smarts = fpaths.get(smartsPathName + '_end')
    
    if not path_smarts or not end_smarts:
        return set()
    
    # Convert SMARTS to molecules
    try:
        path_mol = Chem.MolFromSmarts(path_smarts)
        end_mol = Chem.MolFromSmarts(end_smarts)
    except:
        return set()
    
    if not path_mol or not end_mol:
        return set()
    
    # Find all atoms that are part of fluorinated chains
    fluorinated_atoms = set()
    
    # Get all matches for path and end patterns
    path_matches = mol.GetSubstructMatches(path_mol)
    end_matches = mol.GetSubstructMatches(end_mol)
    
    # Add all atoms from path matches
    for match in path_matches:
        fluorinated_atoms.update(match)
    
    # Add all atoms from end matches
    for match in end_matches:
        fluorinated_atoms.update(match)
    
    return fluorinated_atoms


def get_non_fluorinated_bonds(mol, fluorinated_atoms=None, smartsPathName='Perfluoroalkyl'):
    """
    Get bond indices that are not part of the fluorinated chain.
    
    Args:
        mol: RDKit molecule object
        fluorinated_atoms: Set of atom indices that are part of fluorinated chains
        smartsPathName: Name of the SMARTS pattern to use for identifying fluorinated chains
    
    Returns:
        List of bond indices that are not in the fluorinated chain
    """
    if fluorinated_atoms is None:
        fluorinated_atoms = find_fluorinated_chains(mol, smartsPathName)
    
    non_fluorinated_bonds = []
    
    for bond in mol.GetBonds():
        begin_idx = bond.GetBeginAtomIdx()
        end_idx = bond.GetEndAtomIdx()
        
        # If both atoms are in the fluorinated chain, keep the bond intact
        if begin_idx in fluorinated_atoms and end_idx in fluorinated_atoms:
            continue
        
        # If bond connects fluorinated chain to other parts, it can be broken
        # If bond is between non-fluorinated atoms, it can be broken
        non_fluorinated_bonds.append(bond.GetIdx())
    
    return non_fluorinated_bonds


def get_bonds_connected_to_fluorinated_path(mol, fluorinated_atoms=None, smartsPathName='Perfluoroalkyl'):
    """
    Get bond indices that are connected to (but not part of) the fluorinated path.
    
    Args:
        mol: RDKit molecule object
        fluorinated_atoms: Set of atom indices that are part of fluorinated chains
        smartsPathName: Name of the SMARTS pattern to use for identifying fluorinated chains
    
    Returns:
        List of bond indices that are connected to but not part of the fluorinated chain
    """
    if fluorinated_atoms is None:
        fluorinated_atoms = find_fluorinated_chains(mol, smartsPathName)
    
    connected_bonds = []
    
    for bond in mol.GetBonds():
        begin_idx = bond.GetBeginAtomIdx()
        end_idx = bond.GetEndAtomIdx()
        
        # Skip bonds that are entirely within the fluorinated chain
        if begin_idx in fluorinated_atoms and end_idx in fluorinated_atoms:
            continue
        
        # Only include bonds where exactly one atom is in the fluorinated chain
        # This ensures we're breaking bonds that connect the fluorinated path to other parts
        if (begin_idx in fluorinated_atoms and end_idx not in fluorinated_atoms) or \
           (begin_idx not in fluorinated_atoms and end_idx in fluorinated_atoms):
            connected_bonds.append(bond.GetIdx())
    
    return connected_bonds


def generate_degradation_products(mol, smartsPathName='Perfluoroalkyl', max_breaks=3, 
                                min_fragment_size=2, include_original=True):
    """
    Generate degradation products by breaking bonds that are connected to but not part of the fluorinated chain.
    
    Args:
        mol: RDKit molecule object
        smartsPathName: Name of the SMARTS pattern to use for identifying fluorinated chains
        max_breaks: Maximum number of bonds to break simultaneously
        min_fragment_size: Minimum number of atoms in a fragment to keep it
        include_original: Whether to include the original molecule in results
    
    Returns:
        Dictionary with InChI keys mapping to fragment information
    """
    if not isinstance(mol, Chem.Mol):
        raise TypeError("Input mol must be an rdkit.Chem.Mol object")
    
    # Find fluorinated chain atoms
    fluorinated_atoms = find_fluorinated_chains(mol, smartsPathName)
    
    # Get bonds that can be broken (connected to but not part of fluorinated chain)
    breakable_bonds = get_bonds_connected_to_fluorinated_path(mol, fluorinated_atoms, smartsPathName)
    
    if not breakable_bonds:
        # No breakable bonds found, return original molecule
        if include_original:
            inchi_key = Chem.MolToInchiKey(mol)
            formula = CalcMolFormula(mol)
            return {inchi_key: {formula: {'mol': mol, 'n_breaks': 0, 'broken_bonds': []}}}
        else:
            return {}
    
    degradation_products = {}
    
    # Include original molecule if requested
    if include_original:
        try:
            inchi_key = Chem.MolToInchiKey(mol)
            formula = CalcMolFormula(mol)
            degradation_products[inchi_key] = {
                formula: {'mol': mol, 'n_breaks': 0, 'broken_bonds': []}
            }
        except:
            pass
    
    # Generate all possible combinations of bond breaks
    for n_breaks in range(1, min(max_breaks + 1, len(breakable_bonds) + 1)):
        for bond_combination in combinations(breakable_bonds, n_breaks):
            try:
                # Break the bonds
                rwmol = Chem.RWMol(mol)
                
                # Sort bonds in descending order to avoid index shifting issues
                sorted_bonds = sorted(bond_combination, reverse=True)
                
                for bond_idx in sorted_bonds:
                    bond = rwmol.GetBondWithIdx(bond_idx)
                    rwmol.RemoveBond(bond.GetBeginAtomIdx(), bond.GetEndAtomIdx())
                
                # Get fragments
                fragments = Chem.GetMolFrags(rwmol, asMols=True, sanitizeFrags=True)
                
                # Process each fragment
                for frag in fragments:
                    # Skip fragments that are too small
                    if frag.GetNumAtoms() < min_fragment_size:
                        continue
                    
                    try:
                        inchi_key = Chem.MolToInchiKey(frag)
                        formula = CalcMolFormula(frag)
                        
                        if inchi_key not in degradation_products:
                            degradation_products[inchi_key] = {}
                        
                        if formula not in degradation_products[inchi_key]:
                            degradation_products[inchi_key][formula] = {
                                'mol': frag, 
                                'n_breaks': n_breaks, 
                                'broken_bonds': list(bond_combination)
                            }
                    except:
                        # Skip problematic fragments
                        continue
                        
            except Exception as e:
                # Skip problematic bond breaking combinations
                continue
    
    return degradation_products


def generate_systematic_degradation_products(mol, smartsPathName='Perfluoroalkyl', 
                                           preserve_fluorinated_chain=True, 
                                           min_fragment_size=2, 
                                           max_combinations=1000):
    """
    Generate systematic degradation products by exploring all possible bond breaking patterns.
    
    Args:
        mol: RDKit molecule object
        smartsPathName: Name of the SMARTS pattern to use for identifying fluorinated chains
        preserve_fluorinated_chain: If True, only break bonds connected to (not part of) fluorinated chain
        min_fragment_size: Minimum number of atoms in a fragment to keep it
        max_combinations: Maximum number of bond combinations to explore
    
    Returns:
        Dictionary with degradation products organized by number of breaks
    """
    if not isinstance(mol, Chem.Mol):
        raise TypeError("Input mol must be an rdkit.Chem.Mol object")
    
    # Find fluorinated chain atoms
    fluorinated_atoms = find_fluorinated_chains(mol, smartsPathName)
    
    # Get all bonds in the molecule
    all_bonds = list(range(mol.GetNumBonds()))
    
    if preserve_fluorinated_chain:
        # Only consider bonds connected to but not part of the fluorinated chain
        breakable_bonds = get_bonds_connected_to_fluorinated_path(mol, fluorinated_atoms, smartsPathName)
    else:
        # Consider all bonds
        breakable_bonds = all_bonds
    
    if not breakable_bonds:
        return {"0_breaks": {Chem.MolToInchiKey(mol): {CalcMolFormula(mol): mol}}}
    
    degradation_products = {}
    combinations_explored = 0
    
    # Explore different numbers of breaks
    for n_breaks in range(len(breakable_bonds) + 1):
        degradation_products[f"{n_breaks}_breaks"] = {}
        
        if n_breaks == 0:
            # Original molecule
            try:
                inchi_key = Chem.MolToInchiKey(mol)
                formula = CalcMolFormula(mol)
                degradation_products[f"{n_breaks}_breaks"][inchi_key] = {formula: mol}
            except:
                pass
            continue
        
        # Generate combinations of bonds to break
        for bond_combination in combinations(breakable_bonds, n_breaks):
            if combinations_explored >= max_combinations:
                break
                
            combinations_explored += 1
            
            try:
                # Create a copy and break the bonds
                rwmol = Chem.RWMol(mol)
                
                # Sort bonds in descending order to avoid index shifting
                sorted_bonds = sorted(bond_combination, reverse=True)
                
                for bond_idx in sorted_bonds:
                    bond = rwmol.GetBondWithIdx(bond_idx)
                    rwmol.RemoveBond(bond.GetBeginAtomIdx(), bond.GetEndAtomIdx())
                
                # Get fragments
                fragments = Chem.GetMolFrags(rwmol, asMols=True, sanitizeFrags=True)
                
                # Store each fragment
                for frag in fragments:
                    if frag.GetNumAtoms() < min_fragment_size:
                        continue
                    
                    try:
                        inchi_key = Chem.MolToInchiKey(frag)
                        formula = CalcMolFormula(frag)
                        
                        if inchi_key not in degradation_products[f"{n_breaks}_breaks"]:
                            degradation_products[f"{n_breaks}_breaks"][inchi_key] = {}
                        
                        degradation_products[f"{n_breaks}_breaks"][inchi_key][formula] = frag
                        
                    except:
                        continue
                        
            except:
                continue
        
        if combinations_explored >= max_combinations:
            break
    
    return degradation_products


def analyze_degradation_pathways(mol, smartsPathName='Perfluoroalkyl', max_breaks=3):
    """
    Analyze possible degradation pathways and return detailed information about fragments.
    
    Args:
        mol: RDKit molecule object
        smartsPathName: Name of the SMARTS pattern to use for identifying fluorinated chains
        max_breaks: Maximum number of bonds to break
    
    Returns:
        Dictionary with detailed pathway analysis
    """
    if not isinstance(mol, Chem.Mol):
        raise TypeError("Input mol must be an rdkit.Chem.Mol object")
    
    # Get degradation products
    degradation_products = generate_degradation_products(
        mol, smartsPathName=smartsPathName, max_breaks=max_breaks
    )
    
    # Analyze the results
    analysis = {
        'original_molecule': {
            'smiles': Chem.MolToSmiles(mol),
            'formula': CalcMolFormula(mol),
            'molecular_weight': ExactMolWt(mol),
            'fluorinated_atoms': len(find_fluorinated_chains(mol, smartsPathName))
        },
        'degradation_summary': {
            'total_products': len(degradation_products),
            'products_by_breaks': {}
        },
        'products': []
    }
    
    # Organize products by number of breaks
    for inchi_key, formulas in degradation_products.items():
        for formula, data in formulas.items():
            frag_mol = data['mol']
            n_breaks = data['n_breaks']
            
            if n_breaks not in analysis['degradation_summary']['products_by_breaks']:
                analysis['degradation_summary']['products_by_breaks'][n_breaks] = 0
            analysis['degradation_summary']['products_by_breaks'][n_breaks] += 1
            
            # Analyze fragment
            frag_fluorinated_atoms = find_fluorinated_chains(frag_mol, smartsPathName)
            
            product_info = {
                'inchi_key': inchi_key,
                'smiles': Chem.MolToSmiles(frag_mol),
                'formula': formula,
                'molecular_weight': ExactMolWt(frag_mol),
                'num_atoms': frag_mol.GetNumAtoms(),
                'num_bonds': frag_mol.GetNumBonds(),
                'fluorinated_atoms': len(frag_fluorinated_atoms),
                'contains_fluorinated_chain': len(frag_fluorinated_atoms) > 0,
                'n_breaks': n_breaks,
                'broken_bonds': data['broken_bonds']
            }
            
            analysis['products'].append(product_info)
    
    # Sort products by number of breaks and molecular weight
    analysis['products'].sort(key=lambda x: (x['n_breaks'], x['molecular_weight']))
    
    return analysis
