#!/usr/bin/env python3

from PFASgroups.core import parse_mol

def test_trifluoromethylamine():
    """Test FC(F)(F)N against all PFAS definitions."""
    
    mol_smiles = 'FC(F)(F)N'
    print(f"Testing {mol_smiles} (trifluoromethylamine):")
    print("="*50)
    
    try:
        results = parse_mol(mol_smiles)
        
        print("PFAS Definition Results:")
        for def_name, result in results['definitions'].items():
            status = "MATCHES" if result else "no match"
            print(f"  {def_name:20} : {status}")
            
        print(f"\nOverall PFAS: {results['pfas']}")
        
        if results['pfas']:
            print("\nDetected by groups:")
            for group_result in results['groups']:
                if group_result.get('detected', False):
                    print(f"  - {group_result['name']}")
    
    except Exception as e:
        print(f"Error: {e}")

if __name__ == "__main__":
    test_trifluoromethylamine()