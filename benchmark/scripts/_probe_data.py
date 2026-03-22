"""Probe the pre-computed fingerprint TSV files."""
import sys
sys.path.insert(0, r'C:\Users\luc\git\PFASGroups')
import pandas as pd
import numpy as np
from rdkit import Chem

DATA = r'C:\Users\luc\git\PFASGroups\benchmark\test_data'
txp = pd.read_csv(f'{DATA}/TxP_PFAS_v1.tsv', sep='\t', index_col=0)
tp  = pd.read_csv(f'{DATA}/toxprint_V2.tsv',  sep='\t', index_col=0)
smi = [l.strip() for l in open(f'{DATA}/toxcast_smiles.smi')]

print('TxP_PFAS shape:', txp.shape)
print('ToxPrint  shape:', tp.shape)
print('SMILES count:', len(smi))
print('TxP first col:', repr(txp.columns[0]))
print('TP  first col:', repr(tp.columns[0]))
print('TxP first col sum:', txp.iloc[:, 0].sum())
print('TP  first col sum:', tp.iloc[:, 0].sum())
print('Non-zero TxP rows:', int((txp.values > 0).any(axis=1).sum()))
print('Non-zero TP  rows:', int((tp.values  > 0).any(axis=1).sum()))

# CF-containing
cf_idx = []
for i, s in enumerate(smi):
    m = Chem.MolFromSmiles(s)
    if m and any(a.GetAtomicNum() == 9 for a in m.GetAtoms()):
        cf_idx.append(i)
print('F-containing molecules:', len(cf_idx))

# check M_READ_ERROR in TxP
if 'M_READ_ERROR' in txp.columns[0]:
    bad = (txp.iloc[:, 0] == 1).sum()
    print('M_READ_ERROR (parse failures):', bad)
