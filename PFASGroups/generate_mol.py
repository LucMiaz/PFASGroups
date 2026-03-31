"""
Molecule generation utilities for HalogenGroups.

This module provides functions to generate random PFAS-like molecules by:
- Creating carbon chain backbones (linear, cyclic, with multiple bonds)
- Fluorinating the chains (full or partial fluorination)
- Attaching functional groups (carboxylic acids, sulfonic acids, etc.)

Typical workflow:
1. generate_random_carbon_chain() - Create backbone
2. fluorinate_mol() - Add fluorine atoms
3. append_functional_groups() - Add functional groups

Or use generate_random_mol() which combines all steps.
"""

import numpy as np
from rdkit import Chem
from itertools import groupby

np.random.seed(2025)  # For reproducibility

def generate_random_carbon_chain(n, cycle=False, alkene=False, alkyne=False, branching_range = None):
    """Generate a random carbon chain with n carbons and optional unsaturation.

    Creates molecular backbones by randomly connecting carbon atoms with
    single, double, or triple bonds based on specified probabilities.

    Parameters
    ----------
    n : int
        Number of carbon atoms in the chain
    cycle : bool, default=False
        If True, starts with benzene ring instead of linear chain
    alkene : bool, default=False
        If True, allows C=C double bonds (50% probability per bond)
    alkyne : bool, default=False
        If True, allows C≡C triple bonds (50% probability per bond if no double bond)
    branching_range : tuple of float, optional
        If specified, the allowed range for the branching index of the generated molecule (1.0 = fully linear). After generating a molecule with n carbons, the branching index is calculated and if it falls outside the range, additional carbon atoms are added until the branching index is within the range or a maximum number of attempts is reached.
    Returns
    -------
    rdkit.Chem.Mol
        RDKit molecule object with hydrogen atoms implicit

    Examples
    --------
    >>> # Linear saturated chain (C5H12)
    >>> mol = generate_random_carbon_chain(5)

    >>> # Chain with possible double bonds
    >>> mol = generate_random_carbon_chain(6, alkene=True)

    >>> # Aromatic/cyclic starting point
    >>> mol = generate_random_carbon_chain(8, cycle=True)

    Notes
    -----
    - Carbons are added one at a time with random attachment points
    - Valence rules are enforced (max 4 bonds per carbon)
    - Double/triple bonds are added probabilistically during construction
    - Final molecule is sanitized to ensure valid chemistry
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
            elif rwm.GetAtoms()[j].GetValence(which=Chem.rdchem.ValenceType.EXPLICIT) == 4:
                j = None
        bondtype = Chem.BondType.SINGLE
        if alkene and np.random.uniform(0,1) < 0.5 and rwm.GetAtoms()[j].GetValence(which=Chem.rdchem.ValenceType.EXPLICIT) < 3:
            bondtype = Chem.BondType.DOUBLE
        elif alkyne and bondtype != Chem.BondType.DOUBLE and np.random.uniform(0,1) < 0.5 and rwm.GetAtoms()[j].GetValence(which=Chem.rdchem.ValenceType.EXPLICIT) < 2:
            bondtype = Chem.BondType.TRIPLE
        rwm.AddBond(i, j, bondtype)
        Chem.SanitizeMol(rwm)
    m2 = rwm.GetMol()
    if branching_range is not None:
        branching_index = get_branching_index(m2)
        attempts = 0
        while (branching_index < branching_range[0] or branching_index > branching_range[1]) and attempts < n*5:
            i = rwm.AddAtom(Chem.Atom(6))
            j = None
            while j is None:
                j = np.random.randint(0, len(rwm.GetAtoms())-1)
                if j == i:
                    j = None
                elif rwm.GetAtoms()[j].GetValence(which=Chem.rdchem.ValenceType.EXPLICIT) == 4:
                    j = None
            rwm.AddBond(i, j, Chem.BondType.SINGLE)
            Chem.SanitizeMol(rwm)
            m2 = rwm.GetMol()
            branching_index = get_branching_index(m2)
            attempts += 1
    Chem.SanitizeMol(m2)
    return m2

def get_branching_index(mol):
    """
    Calculate branching of carbon atoms in a molecules
    Measure of branching vs linearity
    For linear chains: branching → 1.0
    For highly branched: branching → 0.0

    1. extract carbon atoms and their connectivity
    2. Count branch points (degree > 2 in the carbon subgraph)
    3. Normalize by component size
    """
    # Count branch points: carbons with more than 2 carbon-carbon bonds
    # (only C neighbours are counted so that functional-group heteroatoms, e.g.
    #  the two oxygens on a -COOH terminus, do not create false branch points)
    carbon_nodes = [atom.GetIdx() for atom in mol.GetAtoms() if atom.GetSymbol() == 'C']
    branch_points = sum(
        max(0, sum(1 for nb in mol.GetAtomWithIdx(node).GetNeighbors() if nb.GetSymbol() == 'C') - 2) for node in carbon_nodes
    )
    # Normalize by component size
    branching_index = 1.0 - 2 * (branch_points / max(1, len(carbon_nodes)))
    branching_index = max(0.0, min(1.0, branching_index))  # Clamp to [0, 1]
    return branching_index

def get_attachment(mol, m, atom_symbols = ['C'], neighbors_symbols = {'C':['F','H']}):
    """Find attachment points in a molecule for adding functional groups.

    Searches for atoms that meet criteria for functional group attachment,
    typically carbon atoms with F or H neighbors that can be replaced.

    Parameters
    ----------
    mol : rdkit.Chem.Mol
        Molecule to search for attachment points
    m : int
        Number of attachment points to return (randomly selected)
    atom_symbols : list of str, default=['C']
        Atom types to consider as attachment points (e.g., ['C'], ['C', 'N'])
    neighbors_symbols : dict, default={'C': ['F', 'H']}
        For each atom type, which neighbor types are valid for replacement

    Returns
    -------
    list of tuple
        List of (atom_index, neighbor_index) tuples indicating where functional
        groups can be attached. Length is min(m, available_candidates).

    Examples
    --------
    >>> mol = Chem.MolFromSmiles("FC(F)C(F)F")
    >>> sites = get_attachment(mol, m=2, atom_symbols=['C'], neighbors_symbols={'C': ['F']})
    >>> # Returns 2 random (C_idx, F_neighbor_idx) pairs

    Notes
    -----
    - Removes duplicate bonds (a,b) and (b,a) from candidates
    - Random selection without replacement if more candidates than m
    - Used internally by append_functional_group for automatic site selection
    """
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
    """Attach a submolecule to a main molecule at a specific atom.

    Removes a random F or H neighbor from the attachment atom and bonds
    the submolecule at that position.

    Parameters
    ----------
    mol : rdkit.Chem.Mol
        Main molecule
    submol : rdkit.Chem.Mol
        Submolecule to attach (functional group)
    atom_index : int
        Index of atom in main molecule to attach to

    Returns
    -------
    rdkit.Chem.Mol
        Modified molecule with submolecule attached

    Notes
    -----
    - Attaches submolecule at its atom index 0
    - Randomly selects F or H neighbor of attachment atom to replace
    - Creates single bond between main molecule and submolecule
    """
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
    """Insert a submolecule between two bonded atoms.

    Breaks the bond between two atoms and inserts the functional group,
    creating new bonds from the functional group to each original atom.

    Parameters
    ----------
    mol : rdkit.Chem.Mol
        Main molecule
    group_mol : rdkit.Chem.Mol
        Functional group to insert
    atom_index : int
        Index of first atom in bond to break
    neighbor_index : int
        Index of second atom in bond to break

    Returns
    -------
    rdkit.Chem.Mol
        Modified molecule with functional group inserted

    Notes
    -----
    - Bonds first atom of functional group to atom_index
    - Bonds last atom of functional group (with available valence) to neighbor_index
    - Original bond is removed before insertion
    - Useful for adding linkers or bridging groups
    """
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
        while rwm.GetAtoms()[end_atom_index].GetValence(which=Chem.rdchem.ValenceType.IMPLICIT) == 0:
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


def get_attachment(mol, m, atom_symbols = ['C'], neighbors_symbols = {'C':['F','H']}):  # pylint: disable=function-redefined
    candidates = []
    for atom in mol.GetAtoms():
        # Check if the atom is 4 bonds and has H or F neighbors
        if atom.GetSymbol() in atom_symbols:
            neighbors = [(atom.GetIdx(),x.GetIdx()) for x in atom.GetNeighbors() if x.GetSymbol() in neighbors_symbols.get(atom.GetSymbol(),['F','H'])]
            candidates+=neighbors
    # check that no bond appears two times:
    candidates = [x for i,x in enumerate(candidates) if i == len(candidates) or ((x[1],x[0]) not in candidates[i+1:] and (x[0],x[1]) not in candidates[i+1:])]
    return [candidates[x] for x in np.random.choice(range(len(candidates)),size = min(len(candidates),m), replace = False)]

def append_functional_group(mol, group_smiles, insertion = 'attach', m=1, atom_indices:list = None, neighbor_atoms = ['C'],sanitize = False):  # pylint: disable=function-redefined
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
    """Fluorinate a molecule by replacing hydrogen atoms with fluorine.

    Converts H atoms to F atoms based on perfluorination setting and
    probability parameters that adapt during fluorination.

    Parameters
    ----------
    mol : rdkit.Chem.Mol
        Molecule to fluorinate (H atoms added automatically)
    perfluorinated : bool, default=True
        If True, replaces all H with F. If False, uses probabilistic replacement.
    p : float, default=0.3
        Minimum probability of H→F replacement (when polyfluorinated)
    phigh : float, default=1
        Maximum probability of H→F replacement. Actual maximum is min(phigh, 1-nF/nH)
        where nF is current F count and nH is total initial H count.

    Returns
    -------
    rdkit.Chem.Mol
        Fluorinated molecule

    Examples
    --------
    >>> mol = Chem.MolFromSmiles("CCCC")
    >>> # Fully fluorinated (C4F10)
    >>> perfluoro = fluorinate_mol(mol, perfluorinated=True)
    >>> # Partially fluorinated (probabilistic)
    >>> polyfluoro = fluorinate_mol(mol, perfluorinated=False, p=0.5)

    Notes
    -----
    - Explicit H atoms are added before fluorination
    - For perfluorinated=False, probability adapts as fluorination progresses
    - Probability formula: max(p, min(phigh, 1 - nF/nH))
    - H atoms are processed in random order
    """
    mol = Chem.AddHs(mol)
    rwm = Chem.RWMol(mol)
    nF = 0
    # Only fluorinate hydrogens attached to carbon atoms.
    # This avoids converting heteroatom-H bonds (e.g., Si-H, N-H) into Si-F/N-F.
    atomsH_index = []
    for atom in rwm.GetAtoms():
        if atom.GetSymbol() == 'H':
            neighbors = atom.GetNeighbors()
            if any(n.GetSymbol() == 'C' for n in neighbors):
                atomsH_index.append(atom.GetIdx())
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
    """Append multiple functional groups to a molecule with complex specifications.

    Handles functional groups with variable occurrence counts, composite groups
    (multiple items), and different insertion modes.

    Parameters
    ----------
    mol : rdkit.Chem.Mol
        Base molecule
    functional_groups : list of dict
        List of functional group specifications, each with:
        - 'group_smiles' (str): SMILES of functional group
        - 'n' (int or str): Number of occurrences. String format "[min,max]" for random count
        - 'mode' (str): 'attach' or 'insert'
        - 'neighbours' (list): Neighbor atom types for attachment
        - 'items' (list, optional): Sub-components with their own 'smiles' and 'n'
    **kwargs : dict
        Additional parameters:
        - chain_n (int): Carbon chain length (used to cap maximum functional group count)

    Returns
    -------
    rdkit.Chem.Mol
        Molecule with all functional groups attached

    Examples
    --------
    >>> mol = Chem.MolFromSmiles("FC(F)C(F)C(F)F")
    >>> groups = [
    ...     {'group_smiles': 'C(=O)O', 'n': 1, 'mode': 'attach'},
    ...     {'group_smiles': 'O', 'n': '[1,3]', 'mode': 'insert'}
    ... ]
    >>> mol = append_functional_groups(mol, groups, chain_n=4)

    Notes
    -----
    - Attachment sites are pre-computed for all groups to avoid conflicts
    - Groups with 'items' have their SMILES concatenated (e.g., polyethylene oxide)
    - Random counts are drawn once at the start and constrained by chain_n
    - Each group is added sequentially with intermediate sanitization
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

def generate_random_mol(n, functional_groups=None, perfluorinated=True, cycle=False, alkene=False, alkyne=False, **kwargs):
    """Generate a random PFAS-like molecule with functional groups.

    Complete pipeline combining chain generation, fluorination, and functional
    group attachment to create synthetic PFAS molecules.

    Parameters
    ----------
    n : int
        Number of carbon atoms in the base chain
    functional_groups : dict, str, or list
        Functional group specifications:
        - dict: Single group with 'group_smiles', 'n', 'mode', 'neighbours'
        - str: SMILES string (uses defaults: m=kwargs['m'] or 1, mode='attach')
        - list: Multiple groups (dicts or strings)
    perfluorinated : bool, default=True
        If True, fully fluorinate. If False, partial fluorination.
    cycle : bool, default=False
        Start with cyclic structure (benzene)
    alkene : bool, default=False
        Allow C=C double bonds
    alkyne : bool, default=False
        Allow C≡C triple bonds
    **kwargs : dict
        Additional parameters passed to fluorination and functional group attachment

    Returns
    -------
    rdkit.Chem.Mol
        Complete PFAS-like molecule

    Examples
    --------
    >>> # Perfluoroalkyl carboxylic acid (PFAA)
    >>> mol = generate_random_mol(8, "C(=O)O", perfluorinated=True)

    >>> # Complex molecule with multiple functional groups
    >>> groups = [
    ...     {'group_smiles': 'C(=O)O', 'n': 1, 'mode': 'attach'},
    ...     {'group_smiles': 'S(=O)(=O)O', 'n': '[0,2]', 'mode': 'insert'}
    ... ]
    >>> mol = generate_random_mol(10, groups, perfluorinated=False, alkene=True)

    Notes
    -----
    Pipeline order:
    1. Generate carbon chain backbone
    2. Fluorinate the chain
    3. Attach functional groups

    This is the primary high-level function for PFAS molecule generation.
    """
    if functional_groups is None:
        functional_groups = []
    elif isinstance(functional_groups, dict):
        functional_groups = [functional_groups]
    elif isinstance(functional_groups, str) :
        functional_groups = [{'group_smiles':functional_groups,'n':kwargs.get('m',1),'mode':kwargs.get('mode','attach')}]
    elif isinstance(functional_groups,list):
        if len(functional_groups) >0 and isinstance(functional_groups[0],str):
            functional_groups = [{'group_smiles':functional_group[0],'n':kwargs.get('m',1),'mode':kwargs.get('mode','attach')} for functional_group in functional_groups]
    _p     = kwargs.get('p', 0.3)
    _phigh = kwargs.get('phigh', 1)
    _max_defs     = kwargs.get('max_matched_definitions')
    _skip_keys    = {'branching_range', 'max_matched_definitions', 'max_generation_attempts', 'p', 'phigh'}
    _max_attempts = kwargs.get('max_generation_attempts', n * 10) if _max_defs is not None else 1
    for _ in range(_max_attempts):
        _fg = [dict(g) for g in functional_groups]  # copy so 'n' mutations don't bleed across retries
        mol = generate_random_carbon_chain(n, cycle, alkene, alkyne, branching_range=kwargs.get('branching_range'))
        mol = fluorinate_mol(mol, perfluorinated=perfluorinated, p=_p, phigh=_phigh)
        if _fg:
            mol = append_functional_groups(mol, _fg, chain_n=n,
                                           **{k: v for k, v in kwargs.items() if k not in _skip_keys})
        if _max_defs is not None:
            from PFASGroups import parse_mol as _parse_mol  # lazy import – avoids circular dependency
            _result = _parse_mol(mol, include_PFAS_definitions=True)
            _detected = {m['id'] for m in _result.get('matches', []) if m.get('type') == 'PFASdefinition'}
            if len(_detected) <= _max_defs:
                return mol
        else:
            return mol
    return mol  # return last attempt if constraint was not satisfied within max_attempts
