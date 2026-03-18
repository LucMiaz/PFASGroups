import warnings
warnings.filterwarnings('always')
import sys
sys.path.insert(0, r'c:\Users\luc\git\PFASGroups')
from PFASGroups import parse_smiles

# Test 1: valid SMILES
result = parse_smiles('OC(=O)C(F)(F)C(F)(F)F')
print(f'Test 1 (valid): {len(result)} result(s), error={result[0].get("error")}')

# Test 2: invalid SMILES in a batch
with warnings.catch_warnings(record=True) as w:
    warnings.simplefilter('always')
    result2 = parse_smiles(['OC(=O)C(F)(F)C(F)(F)F', 'INVALID_SMILES_!!', 'FC(F)(F)F'])
    print(f'Test 2 (mixed batch): {len(result2)} result(s)')
    for i, r in enumerate(result2):
        err = r.get("error")
        smi = r.get("smiles")
        print(f'  [{i}] smiles={smi!r:30s} error={err}')
    if w:
        print(f'Warnings: {[str(x.message) for x in w]}')
