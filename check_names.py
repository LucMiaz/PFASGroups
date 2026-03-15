import sys
sys.path.insert(0, r'C:\Users\luc\git\PFASGroups')
from PFASGroups.parser import parse_smiles

pfoa_smiles = "C(C(C(C(C(C(C(C(=O)O)(F)F)(F)F)(F)F)(F)F)(F)F)(F)F)(F)(F)F"
r = parse_smiles([pfoa_smiles], bycomponent=True)
with open(r'C:\Users\luc\git\PFASGroups\names_out.txt', 'w') as f:
    if r and 'matches' in r[0]:
        for m in r[0]['matches']:
            f.write(m.get('group_name', '') + '\n')
    else:
        f.write('NO MATCHES\n')
        f.write(repr(r) + '\n')
