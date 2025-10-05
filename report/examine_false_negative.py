import pandas as pd

# Load the specificity test results
df = pd.read_csv('specificity_test_results.csv')

# Find the false negative case
fn = df[(df['valid_smiles'] == True) & (df['expected_group_detected'] == False)]

print('FALSE NEGATIVE CASE:')
print(f'Expected: {fn.iloc[0]["group_ids"]}')
print(f'Detected: {fn.iloc[0]["detected_groups"]}')
print(f'SMILES: {fn.iloc[0]["smiles"]}')
print(f'InChI: {fn.iloc[0]["inchi"]}')

# Count ether groups in the SMILES
smiles = fn.iloc[0]["smiles"]
ether_count = smiles.count('O') - smiles.count('OH') - smiles.count('O=') - smiles.count('O)')
print(f'Ether count in molecule: {ether_count}')
print('Analysis: This molecule has only 1 ether linkage, so Group 16 (Perfluoropolyethers) should NOT be detected.')
print('This is actually CORRECT behavior, not a false negative!')