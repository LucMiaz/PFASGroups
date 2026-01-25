"""Simple test to inspect the output structure."""

from PFASgroups.parser import parse_smiles
import json

smiles = 'C(=O)(O)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)F'

print("Testing with bycomponent=True, output_format='dict'")
result = parse_smiles(smiles, bycomponent=True, output_format='dict')
print(f"Type: {type(result)}")

if isinstance(result, list) and len(result) > 0:
    mol_result = result[0]
    print(f"\nKeys: {list(mol_result.keys())}")
    print(f"\nmatches type: {type(mol_result['matches'])}")
    print(f"matches length: {len(mol_result['matches'])}")
    
    if mol_result['matches']:
        print(f"\nFirst match:")
        first_match = mol_result['matches'][0]
        print(json.dumps(first_match, indent=2, default=str))
