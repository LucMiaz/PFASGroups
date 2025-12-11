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

def get_attachment(mol, m, atom_symbols = ['C'], neighbors_symbols = {'C':['F','H']}):
    """Find the first atom in the main molecule to attach the functional group"""
    candidates = []
    for atom in mol.GetAtoms():
        # Check if the atom is 4 bonds and has H or F neighbors
        if atom.GetSymbol() in atom_symbols:
            neighbors = [(atom.GetIdx(),x.GetIdx()) for x in atom.GetNeighbors() if x.GetSymbol() in neighbors_symbols.get(atom.GetSymbol(),['F','H'])]
            candidates.extend(neighbors)
    # check that no bond appears two times:
    candidates = [x for i,x in enumerate(candidates) if i == len(candidates) or ((x[1],x[0]) not in candidates[i+1:] and (x[0],x[1]) not in candidates[i+1:])]
    return [candidates[x] for x in np.random.choice(range(len(candidates)),size = min(len(candidates),m), replace = False)]

def append_functional_group(mol, group_smiles, insertion = 'attach', m=1, atom_indices:list = None, neighbor_atoms = ['C'],sanitize = False):
    """Append a functional group to a molecule.
    Args: insertion: 'attach' to attach to a specific atom, 'insert' place between two bonded atoms, replacing their bond.
    atom_indices: list of lists atom indices + neighbor_indices (atom_index,neighbor_index)
    """
    group_mol = Chem.MolFromSmiles(group_smiles)
    if group_mol is None:
        return mol
    if atom_indices is None:
        atom_indices = get_attachment(mol, m, atom_symbols=neighbor_atoms, neighbors_symbols = {neighbor_atoms[0]:['F','H']})
    for atom_index, neighbor_atom_index in atom_indices:
        if insertion == 'attach':
            mol = attach_mol(mol, group_mol, atom_index)
        elif insertion == 'insert':
            mol = insert_mol(mol, group_mol, atom_index, neighbor_atom_index)
    if sanitize is True:
        Chem.SanitizeMol(mol)
    return mol


def attach_mol(mol, submol, atom_index):
    """Insert a submolecule into a molecule at a specific atom index on mol
    attach at index 0 on submol."""
    rwm = Chem.RWMol(mol)
    submol_index = rwm.GetNumAtoms() # index of the first atom in the submol
    rwm.InsertMol(submol)
    rwm.BeginBatchEdit() # start a batch
    atom = rwm.GetAtomWithIdx(atom_index)
    # choose random neighbor atom (either 'F' or 'H') to remove
    neighbors = [x.GetIdx() for x in atom.GetNeighbors() if x.GetSymbol() in ['F', 'H']]
    atom_to_remove = int(np.random.choice(neighbors))
    rwm.RemoveBond(atom_index, atom_to_remove)
    rwm.RemoveAtom(atom_to_remove)  # remove the atom
    rwm.AddBond(atom_index, submol_index, Chem.BondType.SINGLE)
    rwm.CommitBatchEdit()  # finish the batch
    Chem.SanitizeMol(rwm)
    return rwm.GetMol()


def insert_mol(mol, group_mol, atom_index, neighbor_index):
    """Insert a submolecule into a molecule between two bonded atoms."""
    rwm = Chem.RWMol(mol)
    submol_index = rwm.GetNumAtoms() # index of the first atom in the submol
    rwm.BeginBatchEdit() # start a batch
    # Get the atoms to be replaced
    atom1 = rwm.GetAtomWithIdx(atom_index)
    atom2 = rwm.GetAtomWithIdx(neighbor_index)
    # Remove the existing bond between the two atoms
    rwm.RemoveBond(atom1.GetIdx(), atom2.GetIdx())
    # Insert the submolecule
    rwm.InsertMol(group_mol)
    # Add bonds from the functional group to the original atoms
    rwm.AddBond(atom1.GetIdx(), submol_index, Chem.BondType.SINGLE)
    end_atom_index = rwm.GetNumAtoms()-1
    try:
        while rwm.GetAtoms()[end_atom_index].GetImplicitValence() == 0:
            end_atom_index -= 1
    except:
        pass
    rwm.AddBond(atom2.GetIdx(), end_atom_index, Chem.BondType.SINGLE)
    rwm.CommitBatchEdit()  # finish the batch
    Chem.SanitizeMol(rwm)
    return rwm.GetMol()

def remove_atoms(mol, idxs, removable = ['H','F','Cl','Br','I'], show_on_error = False):
    #print(f'{idxs}')
    #rwm = Chem.RWMol(mol)
    #rwm.BeginBatchEdit() # start a batch
    # Get the atoms to be replaced
    to_remove = []
    bonds = []
    for idx in idxs:
        atom = mol.GetAtomWithIdx(idx)
        neighbors_r = [x.GetIdx() for x in atom.GetNeighbors() if x.GetSymbol() in removable]
        neighbors_c = [x.GetIdx() for x in atom.GetNeighbors() if x.GetSymbol() not in removable]
        if len(neighbors_c)>2:
            raise ValueError(f"There are too many neighbors not in {removable} for atom {idx} in {Chem.MolToSmiles(mol)}, found {neighbors_c=}")
        to_remove = to_remove + neighbors_r + [idx]
        bonds.append([neighbors_c[0],idx])
        bonds.append([idx,neighbors_c[-1]])
    bonds.sort()
    bonds = list(k for k,_ in groupby(bonds))
    new_bonds = []
    last_items = []
    for i,(a,b) in enumerate(bonds):
        if a in last_items:
            for j, (x,y) in enumerate(new_bonds):
                if a == y:
                    last_items.pop(last_items.index(y))
                    last_items.append(b)
                    new_bonds[j][1]=b
        else:
            new_bonds.append([a,b])
            last_items.append(b)
    last_items = []
    bonds = new_bonds
    bonds.sort(reverse = True)
    new_bonds=[]
    for i,(a,b) in enumerate(bonds):
        if a in last_items:
            for j, (x,y) in enumerate(new_bonds):
                if a == y:
                    last_items.pop(last_items.index(y))
                    last_items.append(b)
                    new_bonds[j][1]=b
        else:
            new_bonds.append([a,b])
            last_items.append(b)
    #print(f'{to_remove=},{new_bonds}')
    # for idx in sorted(set(to_remove),reverse=True):
    #     rwm.RemoveAtom(idx)
    # for a,b in new_bonds:
    #     # Add bonds from the functional group to the original atoms
    #     rwm.AddBond(a, b, Chem.BondType.SINGLE)
    # rwm.CommitBatchEdit()  # finish the batch
    # Create a new editable molecule
    rwm = Chem.RWMol()
    _rwm = Chem.RWMol(mol)
    Chem.Kekulize(_rwm)
    # Map from old atom indices to new atom indices
    old_to_new = {}
    charged_atoms = []
    for i, atom in enumerate(_rwm.GetAtoms()):
        if i not in set(to_remove):
            new_atom = Chem.Atom(atom.GetAtomicNum())
            new_idx = rwm.AddAtom(new_atom)
            # Copy formal charge if present
            if atom.GetFormalCharge() != 0:
                charged_atoms.append(atom.GetIdx())
            old_to_new[i] = new_idx

    # Copy bonds that do not involve removed atoms
    for bond in _rwm.GetBonds():
        a1 = bond.GetBeginAtomIdx()
        a2 = bond.GetEndAtomIdx()
        if a1 in old_to_new and a2 in old_to_new:
            rwm.AddBond(old_to_new[a1], old_to_new[a2], bond.GetBondType())

    # Add new bonds
    for a, b in new_bonds:
        if a in old_to_new and b in old_to_new and old_to_new[a]!=old_to_new[b]:
            rwm.AddBond(old_to_new[a], old_to_new[b], Chem.BondType.SINGLE)
    for idx in charged_atoms:
        atom = rwm.GetAtomWithIdx(old_to_new[idx])
        atom.SetFormalCharge(_rwm.GetAtomWithIdx(idx).GetFormalCharge())
    try:
        Chem.SanitizeMol(rwm)
    except Exception as e:
        if show_on_error is True:
            _mol = rwm.GetMol()
            from elements.scripts.utility_image import plot_mol
            img,_,_ = plot_mol(_mol, subwidth=600, subheight=600, svg=False, addAtomIndices=True, bondLineWidth=0.5, fixedBondLength=15, minFontSize=12)
            img.show()
        raise e
    return rwm.GetMol()


def get_attachment(mol, m, atom_symbols = ['C'], neighbors_symbols = {'C':['F','H']}):
    # Find the first atom in the main molecule to attach the functional group
    candidates = []
    for atom in mol.GetAtoms():
        # Check if the atom is 4 bonds and has H or F neighbors
        if atom.GetSymbol() in atom_symbols:
            neighbors = [(atom.GetIdx(),x.GetIdx()) for x in atom.GetNeighbors() if x.GetSymbol() in neighbors_symbols.get(atom.GetSymbol(),['F','H'])]
            candidates+=neighbors
    # check that no bond appears two times:
    candidates = [x for i,x in enumerate(candidates) if i == len(candidates) or ((x[1],x[0]) not in candidates[i+1:] and (x[0],x[1]) not in candidates[i+1:])]
    return [candidates[x] for x in np.random.choice(range(len(candidates)),size = min(len(candidates),m), replace = False)]

def append_functional_group(mol, group_smiles, insertion = 'attach', m=1, atom_indices:list = None, neighbor_atoms = ['C'],sanitize = False):
    """Append a functional group to a molecule.
    Args: insertion: 'attach' to attach to a specific atom, 'insert' place between two bonded atoms, replacing their bond.
    atom_indices: list of lists atom indices + neighbor_indices (atom_index,neighbor_index)
    """
    group_mol = Chem.MolFromSmiles(group_smiles)
    if group_mol is None:
        raise ValueError(f"Invalid SMILES for functional group: {group_smiles}")
    if atom_indices is None:
        atom_indices = get_attachment(mol,m, neighbors_symbols=neighbor_atoms)
    for atom_index, neighbor_atom_index in atom_indices:
        main_atom = mol.GetAtomWithIdx(atom_index)
        # Add the functional group to the main molecule
        if insertion == 'attach':
            mol = attach_mol(mol, group_mol, main_atom.GetIdx())
        elif insertion == 'insert':
            neighbor_atom = mol.GetAtomWithIdx(neighbor_atom_index)
            # Remove the existing bond and add the functional group
            mol = insert_mol(mol, group_mol, main_atom.GetIdx(), neighbor_atom.GetIdx())
    if sanitize is True:
        Chem.SanitizeMol(mol)
    return mol

def fluorinate_mol(mol, perfluorinated=True, p=0.3, phigh=1):
    """
    Fluorinate a molecule by replacing hydrogens with fluorines.
    p is the lowest probability of replacing a H with a F
    phigh is the highest probability of replacing a H with a F, by default it is 1 and will be computed as nb of H / nb of H+F, this is one at the beginning, unless phigh is set below 1
    """
    mol = Chem.AddHs(mol)
    rwm = Chem.RWMol(mol)
    nF = 0
    atomsH_index = [x.GetIdx() for x in rwm.GetAtoms() if x.GetSymbol() == 'H']
    nH = len(atomsH_index)
    np.random.shuffle(atomsH_index)
    prob = []
    for i, atom_index in enumerate(atomsH_index):
        prob.append(max(p, min(phigh, 1-nF/nH)))
        atom = rwm.GetAtomWithIdx(atom_index)
        if perfluorinated or np.random.uniform(0, 1) < prob[-1]:
            atom.SetAtomicNum(9)  # Replace H with F
            nF += 1
        else:
            atom.SetAtomicNum(1)  # Keep H
    Chem.SanitizeMol(rwm)
    return rwm.GetMol()


def append_functional_groups(mol, functional_groups:list, **kwargs):
    """Append a functional group to a molecule.
    Args: insertion: 'attach' to attach to a specific atom, 'insert' place between two bonded atoms, replacing their bond.
    functional_groups as dictionary to include options 'n'= number of occurrence, 'mode'= mode of insertion ('attach' or 'insert') and 'neighbours'=type of neighbours atoms
    """
    total_m = {'attach':0,'insert':0}
    for i,params in enumerate(functional_groups):
        n = params.get('n',1)
        if isinstance(n,str):
            low,high = [int(x) for x in n[1:-1].split(',')]
            n = int(np.random.randint(low,min(kwargs.get('chain_n',high),high),size=1)[0])
            functional_groups[i]['n']=n
        total_m[params.get('mode','attach')] += n
    atom_indices = {}
    atom_indices['attach'] = get_attachment(mol,total_m['attach'], atom_symbols=['C'],neighbors_symbols = {'C':['F','H']})
    atom_indices['insert'] = get_attachment(mol,total_m['insert'], atom_symbols=['C'],neighbors_symbols = {'C':['C']})
    for functional_group in functional_groups:
        group_smiles = functional_group.get('group_smiles','')
        if 'items' in functional_group.keys():
            for item in functional_group['items']:
                ni = item.get('n',1)
                if isinstance(ni,str):
                    low,high = [int(x) for x in ni[1:-1].split(',')]
                    ni = int(np.random.randint(low,high,size=1)[0])
                group_smiles+=item.get('smiles','')*ni
        atom_indices_popped = [atom_indices[functional_group.get('mode','attach')].pop() for _ in range(functional_group.get('n',1))]
        mol = append_functional_group(mol,group_smiles,
                                      insertion = functional_group.get('mode','attach'),
                                      m = functional_group.get('n',1),
                                      atom_indices = atom_indices_popped,
                                      neighbor_atoms=functional_group.get('neighbours',['C']),
                                      sanitize = False)
        try:
            Chem.SanitizeMol(mol)
        except Exception as e:
            print(f"{Chem.MolToSmiles(mol)=}, {functional_group=}, {atom_indices_popped=}")
            raise e
    return mol

def generate_random_mol(n, functional_groups, perfluorinated=True, cycle=False, alkene=False, alkyne=False, **kwargs):
    """Generate a random molecule with n carbon atoms.
    pass functional_groups as dictionary to include options 'n'= number of occurence, 'mode'= mode of insertion ('attach' or 'insert') and 'neighbours'=type of neighbours atoms"""
    if isinstance(functional_groups, dict):
        functional_groups = [functional_groups]
    elif isinstance(functional_groups, str) :
        functional_groups = [{'group_smiles':functional_groups,'n':kwargs.get('m',1),'mode':kwargs.get('mode','attach')}]
    elif isinstance(functional_groups,list) and isinstance(functional_groups[0],str):
        functional_groups = [{'group_smiles':functional_group[0],'n':kwargs.get('m',1),'mode':kwargs.get('mode','attach')} for functional_group in functional_groups]
    mol = generate_random_carbon_chain(n, cycle, alkene, alkyne)
    # Randomly fluorinate the molecule
    mol = fluorinate_mol(mol, perfluorinated=perfluorinated)
    mol = append_functional_groups(mol,functional_groups,chain_n=n,**kwargs)
    return mol
