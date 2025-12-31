"""Debug why amine is detected with 9 non-F carbons"""
from PFASgroups import parse_mol
from PFASgroups.core import get_substruct, get_smartsPaths, path_between_smarts
from rdkit import Chem

smiles = 'FC(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)CCCCCCCCCN=[N+]=[N-]'
mol = Chem.MolFromSmiles(smiles)

amine_smarts = '[#6$([#6;!$([#6]=O)][N;!$(*[#6]=O)])]'
amine_pattern = Chem.MolFromSmarts(amine_smarts)
smartsPaths = get_smartsPaths()

print('Testing 9 non-F carbon case:')
print()

# Get smarts1 matches
smarts1_matches = get_substruct(mol, amine_pattern)
print(f'smarts1 matches (C attached to N): {smarts1_matches}')

# Find paths
match_count, chain_lengths, matched1_len, matched_chains = path_between_smarts(
    mol, amine_pattern, None, smartsPaths=smartsPaths
)

print(f'Initial match_count: {match_count}')
print(f'Number of chains: {len(matched_chains)}')
print()

for i, chain_info in enumerate(matched_chains):
    print(f'Chain {i+1}:')
    print(f'  Type: {chain_info["SMARTS"]}')
    print(f'  Length: {chain_info["length"]}')
    print(f'  Atoms in path: {chain_info["chain"]}')
    
    path_type = chain_info['SMARTS']
    if path_type in smartsPaths:
        end_smarts = smartsPaths[path_type][1]
        end_matches = get_substruct(mol, end_smarts)
        print(f'  END pattern matches: {end_matches}')
        
        chain_atoms = set(chain_info['chain'])
        smarts1_in_chain = chain_atoms.intersection(smarts1_matches)
        print(f'  smarts1 atoms in this chain: {smarts1_in_chain}')
        
        # Calculate distances
        chain_list = chain_info['chain']
        for s1_atom in smarts1_in_chain:
            print(f'  For smarts1 atom {s1_atom}:')
            if s1_atom in end_matches:
                print(f'    - Is an END atom itself')
            else:
                s1_idx = chain_list.index(s1_atom)
                print(f'    - Position in chain: {s1_idx}')
                for end_atom in end_matches:
                    if end_atom in chain_atoms:
                        end_idx = chain_list.index(end_atom)
                        distance = abs(s1_idx - end_idx)
                        print(f'    - Distance to END atom {end_atom} (pos {end_idx}): {distance}')
