from .core import add_smarts, add_smartsPath, get_substruct, remove_atoms, mol_to_nx
from rdkit import Chem
from rdkit.Chem.rdMolDescriptors import CalcMolFormula
import networkx as nx
import re
from itertools import product, groupby

@add_smarts(name = 'repeating')
def find_chain(mol,pathsmarts,endsmarts, repeating = 'C(F)(F)'):
    """Iterates over substructures matched by smarts1 and smarts2 and finds fluorinated paths between them.
    Fluorinated atoms are defined by the SMARTS yield by the decorator add_smartsPath.
    if smarts2 is None, uses SMARTS corresponding to smartsPath."""
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
    chains = sorted(chains, key=lambda x: len(x['chain']), reverse=True)
    remove = []
    for i in range(0,len(chains)):
        if sum([chains[i]['chain'].issubset(x['chain']) for x in chains[0:i]])>0:
            remove.append(i)
    chains = [x for i,x in enumerate(chains) if i not in remove]
    return chains


@add_smartsPath()
def generate_homologues(mol, smartsPathName = 'Perfluoroalkyl', smartsPaths=None, repeating = 'C(F)F', base_repeating = ['C']):
    """finds path, then cuts iteratively the repeating part in the chain
    for perfluoroalkyl use smartsPathName = 'Perfluoroalkyl' and repeating = 'C(F)F',
    for polyfluoroalkyl use smartsPathName = 'Polyfluoroalkyl' and repeating='C([F,H,I,Br,Cl])[F,H,I,Br,Cl]'
    for other non-fluorinated, define a new set of SMARTS for paths"""
    path,endSmarts = smartsPaths[smartsPathName]
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