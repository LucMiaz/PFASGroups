from PFASgroups import parse_smiles
import json

smiles_list = [
    'c1c(C(C(F)(F)F)(C(F)(F)F)C(F)(F)F)c(F)c(C(F)(F)C(=O)O)c(F)c1F',
    'C1(C(F)(F)F)(C(F)(F)F)C(C(F)(F)C(F)(F)F)C(C(F)(F)F)C(C(=O)O)C1(F)F'
]

print("Testing molecules for classification errors:\n")

for i, smiles in enumerate(smiles_list, 1):
    print(f"Molecule {i}: {smiles}")
    try:
        results = parse_smiles(smiles)
        print(f"Type: {type(results)}")
        print(f"Results: {results}")
        if isinstance(results, list) and len(results) > 0:
            result = results[0]
            if hasattr(result, 'groups'):
                print(f"Groups: {result.groups}")
            if hasattr(result, 'success'):
                print(f"Success: {result.success}")
            if hasattr(result, 'error'):
                print(f"Error: {result.error}")
    except Exception as e:
        print(f"Exception: {type(e).__name__}: {e}")
        import traceback
        traceback.print_exc()
    print("-" * 80)
