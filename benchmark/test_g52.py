import sys
sys.path.insert(0, r'c:\Users\luc\git\PFASGroups')
from PFASGroups import parse_mol
from rdkit import Chem, RDLogger
RDLogger.DisableLog('rdApp.*')

tests = [
    ('CC(C)(C)OC(=O)N1C[C@@H](F)C[C@H]1CC(=O)O', False, 'mono-F ring (old false positive)'),
    ('O=C(O)CCF',                                  False, '3-fluoropropanoic acid (1 F)'),
    ('FC(F)(F)C(F)(F)C(F)(F)C(F)(F)F',            True,  'perfluoroalkyl (many F)'),
    ('O=Cc1cc(OC(F)F)ccc1Br',                     True,  '5-bromo-2-(difluoromethoxy)benzaldehyde (2 F on ether)'),
    ('FC(F)(Cl)CCl',                               True,  'CF2Cl-CCl (2 F on one C)'),
    ('ClCCl',                                      False, 'DCM-like (2 Cl, no F, check Cl path)'),
]

ok = True
for smi, expect_g52, name in tests:
    mol = Chem.MolFromSmiles(smi)
    r = parse_mol(mol, include_PFAS_definitions=False)
    groups = [m['group_name'] for m in r.get('matches', []) if m.get('type') == 'HalogenGroup'] if r else []
    g52 = 'polyhalogenated alkyl' in groups
    status = 'OK' if g52 == expect_g52 else 'FAIL'
    if status == 'FAIL':
        ok = False
    print(f'[{status}] {name}: g52={g52} (expected {expect_g52}) | groups={groups}')

print('\nAll OK' if ok else '\nSome tests FAILED')
