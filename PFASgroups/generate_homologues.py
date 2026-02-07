from .core import add_smarts, add_componentSmarts, get_substruct, remove_atoms, mol_to_nx
from rdkit import Chem
from rdkit.Chem.rdMolDescriptors import CalcMolFormula
import networkx as nx
import re
from itertools import product, groupby

@add_smarts(name = 'repeating')
def find_chain(mol,pathsmarts,endsmarts, repeating = 'C(F)(F)'):
    """Find fluorinated chains with repeating units between terminal groups.
    
    Identifies shortest paths between matching end groups that consist entirely
    of fluorinated atoms, then isolates the repeating unit sections within those chains.
    
    Parameters
    ----------
    mol : rdkit.Chem.Mol
        Molecule to search for chains
    pathsmarts : rdkit.Chem.Mol
        SMARTS pattern defining fluorinated path atoms
    endsmarts : rdkit.Chem.Mol
        SMARTS pattern defining chain terminal groups
    repeating : str, default='C(F)(F)'
        SMARTS pattern for repeating unit (e.g., -CF2- for perfluoroalkyl)
    
    Returns
    -------
    list of dict
        Each chain dictionary contains:
        - 'component' (set): All atom indices in the complete chain
        - 'start' (int): Atom index of one terminal group
        - 'end' (int): Atom index of other terminal group
        - 'belly' (set): Atom indices of repeating unit sections
        - 'consecutive_parts' (list): Lists of consecutive repeating unit atoms
    
    Notes
    -----
    - Uses shortest path algorithm to find connections between terminals
    - Filters to only fully fluorinated paths (all atoms match pathsmarts)
    - Removes chains that are subsets of longer chains
    - Chains sorted by length (longest first)
    - Consecutive parts allow systematic removal of CF2 units
    
    Examples
    --------
    >>> # Find perfluoroalkyl chains in PFOA
    >>> chains = find_chain(pfoa_mol, perfluoro_smarts, end_smarts, 'C(F)(F)')
    >>> print(f"Found {len(chains)} chains with belly units: {chains[0]['belly']}")
    """
    chains = []
    try:
        G = mol_to_nx(mol)
    except ValueError as e:
        raise e
    #logger.debug("====================")
    targetMatches = get_substruct(mol, endsmarts)
    pathMatches = get_substruct(mol,pathsmarts)
    pairs = [(x,y) for x in targetMatches for y in targetMatches if x!=y]
    repeated_struct = get_substruct(mol,repeating)
    for match1,match2 in pairs:
        path_idx = nx.shortest_path(G, match1, match2, method='dijkstra')
        setp = set(path_idx)
        if setp == setp.intersection(pathMatches):
            belly = setp.difference([match1,match2]).intersection(repeated_struct)
            chain = {"chain":setp,"start":match1, "end":match2, "belly":belly}
            parts = []
            part = []
            for i, id in enumerate(path_idx):
                if id in belly:
                    part.append(id)
                else:
                    if len(part)>0:
                        parts.append(part)
                        part = []
            chain["consecutive_parts"]=parts
            chains.append(chain)
    chains = sorted(chains, key=lambda x: len(x['component']), reverse=True)
    remove = []
    for i in range(0,len(chains)):
        if sum([chains[i]['component'].issubset(x['component']) for x in chains[0:i]])>0:
            remove.append(i)
    chains = [x for i,x in enumerate(chains) if i not in remove]
    return chains


@add_componentSmarts()
def generate_homologues(mol, componentSmartsName = 'Perfluoroalkyl', componentSmartss=None, repeating = 'C(F)F', base_repeating = ['C']):
    """Generate homologous series by systematically removing repeating units from fluorinated chains.
    
    Finds fluorinated chains with repeating units (e.g., -CF2-) and generates all possible
    shorter homologues by removing combinations of these units while maintaining connectivity.
    
    Parameters
    ----------
    mol : rdkit.Chem.Mol
        Parent molecule to generate homologues from
    componentSmartsName : str, default='Perfluoroalkyl'
        Type of fluorinated chain to search for:
        - 'Perfluoroalkyl': Fully fluorinated chains (-CF2-CF2-)
        - 'Polyfluoroalkyl': Partially fluorinated chains
        - Custom: User-defined in componentSmartss
    componentSmartss : dict, optional
        Dictionary of {name: (path_smarts, end_smarts)} patterns.
        Provided by @add_componentSmarts decorator if not specified.
    repeating : str, default='C(F)F'
        SMARTS pattern for repeating unit to remove:
        - 'C(F)F' for perfluoroalkyl (-CF2-)
        - 'C([F,H,I,Br,Cl])[F,H,I,Br,Cl]' for polyfluoroalkyl
    base_repeating : list of str, default=['C']
        Elements in repeating unit that form the backbone (not removed as neighbors)
    
    Returns
    -------
    dict
        Homologues organized as {InChIKey: {formula: mol}}
        Each unique molecule is stored by its InChIKey with its molecular formula
    
    Examples
    --------
    >>> from rdkit import Chem
    >>> # Generate C2-C7 homologues from C8 PFOA
    >>> pfoa = Chem.MolFromSmiles(\"C(=O)O\" + \"C(F)(F)\" * 7 + \"F\")
    >>> homologues = generate_homologues(pfoa, repeating='C(F)F')
    >>> print(f\"Generated {len(homologues)} unique homologues\")
    
    >>> # For polyfluoroalkyl chains
    >>> homologues = generate_homologues(
    ...     mol, 
    ...     componentSmartsName='Polyfluoroalkyl',
    ...     repeating='C([F,H])[F,H]'
    ... )
    
    Algorithm
    ---------
    1. Find fluorinated chains with terminal groups
    2. Identify consecutive sections of repeating units
    3. Generate all possible sub-chains (partial removals)
    4. Remove atoms for each combination using remove_atoms()
    5. Collect unique products by InChIKey
    
    Notes
    -----
    - Maintains molecular connectivity when removing units
    - Removes entire repeating unit + attached F/Cl/Br/I atoms
    - Filters atoms with <3 non-removable neighbors (prevents over-branching removal)
    - Generates all possible chain lengths from original down to shortest
    - Useful for degradation product prediction and homologous series enumeration
    
    Raises
    ------
    ValueError
        If removal causes molecular fragmentation (indicates issue with chain detection)
    Exception
        If atom removal fails (prints diagnostic information)
    
    See Also
    --------
    find_chain : Identifies chains before homologue generation
    remove_atoms : Performs the actual atom removal with connectivity preservation
    """
    path,endSmarts = componentSmartss[componentSmartsName]
    removable = [x for x in set(re.findall(r'[A-Z][a-z]?',repeating)) if x not in base_repeating]
    subchains = lambda x: [x[0:i] for i in range(1,len(x)+1)]
    chains = find_chain(mol, path, endSmarts, repeating = repeating)
    chains = [part for chain in chains for part in chain['consecutive_parts']]
    chains.sort()
    chains = [[x for x in chain if len([y.GetIdx() for y in mol.GetAtomWithIdx(x).GetNeighbors() if y.GetSymbol() not in removable])<3] for chain in chains]
    chains = [chain for chain in chains if len(chain)>0]
    chains = list(k for k,_ in groupby(chains))
    removable_items = [subchains(chain) for chain in chains]
    homologues = {}
    if len(removable_items) >0:
        for removable_idxs in product(*removable_items):
            flat_idx = [x for y in removable_idxs for x in y]
            try:
                h = remove_atoms(mol, flat_idx, removable = removable)
            except Exception as e:
                print(f"Error removing atoms {flat_idx} from {Chem.MolToSmiles(mol)}, had {chains=}, {Chem.MolToSmarts(path)=},{Chem.MolToSmarts(endSmarts)=}: {e}")
                raise e
            else:
                if len(Chem.GetMolFrags(h))>1:
                    raise ValueError(f"Fragmented molecule {Chem.MolToSmiles(h)} after removing atoms {flat_idx} from {Chem.MolToSmiles(mol)}, had {chains=}, {Chem.MolToSmarts(path)=},{Chem.MolToSmarts(endSmarts)=}")
                inchikey = Chem.MolToInchiKey(h)
                formula = CalcMolFormula(h)
                homologues.setdefault(inchikey,{})[formula] = h
    return homologues