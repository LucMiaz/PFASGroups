"""Debug index out of bounds error"""
from rdkit import Chem
from PFASgroups import parse_smiles

# Problematic SMILES
smiles = "FC(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)CCOC(=O)C1=CC=CN=C1"

print(f"Testing SMILES: {smiles}")

# Parse molecule to check basic properties
mol = Chem.MolFromSmiles(smiles)
if mol is None:
    print("ERROR: Could not parse SMILES")
else:
    print(f"Number of atoms (without H): {mol.GetNumAtoms()}")
    
    # Add hydrogens to see total atom count
    mol_with_h = Chem.AddHs(mol)
    print(f"Number of atoms (with H): {mol_with_h.GetNumAtoms()}")
    
    print("\nAtom details:")
    for i, atom in enumerate(mol_with_h.GetAtoms()):
        print(f"  Atom {i}: {atom.GetSymbol()}")
    
    print("\nNow testing with PFASgroups...")
    try:
        results = parse_smiles([smiles], bycomponent=True, output_format='dict')
        print(f"Success! Found {len(results)} results")
        if results:
            result = results[0]
            matches = result.get('matches', [])
            print(f"Found {len(matches)} group matches")
            for match in matches:
                if isinstance(match, dict):
                    group_name = match.get('group_name', 'Unknown')
                    match_count = match.get('match_count', 0)
                    if match_count > 0:
                        print(f"  - {group_name}: {match_count} matches")
    except Exception as e:
        print(f"ERROR: {type(e).__name__}: {e}")
        import traceback
        traceback.print_exc()
