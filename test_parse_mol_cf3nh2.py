#!/usr/bin/env python3

from rdkit import Chem
from PFASgroups.core import parse_mol

def test_both_smiles():
    """Test both SMILES variations of trifluoromethylamine."""
    
    smiles_variations = [
        'C(F)(F)(F)N',   # Used in EU exclusions test
        'FC(F)(F)N',     # Used in EU specific test
    ]
    
    print("Testing trifluoromethylamine against EU PFAS Restriction")
    print("="*70)
    
    for smiles in smiles_variations:
        print(f"\nTesting SMILES: {smiles}")
        mol = Chem.MolFromSmiles(smiles)
        canonical = Chem.MolToSmiles(mol)
        print(f"Canonical form: {canonical}")
        
        try:
            result = parse_mol(mol, include_PFAS_definitions=True)
            
            if 'matches' in result:
                eu_matched = False
                for match in result['matches']:
                    if match.get('type') == 'PFASdefinition':
                        if match.get('id') == 2:  # EU PFAS Restriction is ID 2
                            eu_matched = True
                            print(f"  EU PFAS Restriction (ID=2): MATCHED")
                            print(f"    Match details: {match}")
                        elif match.get('type') == 'PFASdefinition':
                            print(f"  Definition {match.get('id')} ({match.get('name', 'Unknown')}): matched")
                
                if not eu_matched:
                    print(f"  EU PFAS Restriction (ID=2): NOT MATCHED (correct)")
            else:
                print(f"  Result format unexpected: {result.keys()}")
                
        except Exception as e:
            print(f"  ERROR: {e}")
            import traceback
            traceback.print_exc()

if __name__ == "__main__":
    test_both_smiles()