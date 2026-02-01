#!/usr/bin/env python3

from rdkit import Chem

def test_fixed_eu_patterns():
    """Test that the fixed EU patterns now correctly exclude FC(F)(F)N."""
    
    # Test molecules
    test_cases = [
        ('FC(F)(F)N', 'Trifluoromethylamine (CF3-NH2)', False),
        ('C(F)(F)(F)N', 'Trifluoromethylamine alt SMILES', False),
        ('FC(F)(F)O', 'Trifluoromethanol (CF3-OH)', False),
        ('FC(F)(F)OC', 'Trifluoromethyl methyl ether (CF3-OCH3)', False),
        ('FC(F)(F)F', 'Perfluoromethane (CF4)', True),
        ('FC(F)(F)C(F)(F)F', 'Perfluoroethane', True),
        ('FC(F)(F)C(F)(F)O', 'CF3-CF2-OH (ether at distance)', True),
    ]
    
    # Corrected EU patterns
    eu_cf3 = "[#6X4!$([#6][#8H1,#8$([#8][#6H3,#6H2X4,a,#6X3$([#6]=[#8])])])!$([#6][#7H2,#7H1$([#7][#6H3,#6H2X4,a,#6X3$([#6]=[#8])]),#7$([#7]([#6H3,#6H2X4,a,#6X3$([#6]=[#8])])[#6H3,#6H2X4,a,#6X3$([#6]=[#8])])])!$([#6H1])!$([#6][#17,#35,#53])](F)(F)F"
    
    eu_cf2 = "[#6X4!$([#6](F)(F)(F))!$([#6]([#8H1,#8$([#8][#6H3,#6H2X4,a,#6X3$([#6]=[#8])])])[#6H3,#6H2X4,#6X3$([#6]=[#8]),a,#8H1,#8$([#8][#6H3,#6H2X4,a,#6X3$([#6]=[#8])]),#16$([#16H1,#16$([#16][#6H3,#6H2X4,a,#6X3$([#6]=[#8])])]),#7$([#7H2,#7H1$([#7#6H3,#7#6H2X4,#7a,#6X3$([#6]=[#8])]),#7$([#7]([#6H3,#6H2X4,a,#6X3$([#6]=[#8])])[#6H3,#6H2X4,a,#6X3$([#6]=[#8])])])])!$([#6]([#7H2,#7H1$([#7][#6H3,#6H2X4,a,#6X3$([#6]=[#8])]),#7$([#7]([#6H3,#6H2X4,a,#6X3$([#6]=[#8])])[#6H3,#6H2X4,a,#6X3$([#6]=[#8])])])[#6H3,#6H2X4,#6X3$([#6]=[#8]),a,#8H1,#8$([#8][#6H3,#6H2X4,a,#6X3$([#6]=[#8])]),#16$([#16H1,#16$([#16][#6H3,#6H2X4,a,#6X3$([#6]=[#8])])]),#7$([#7H2,#7H1$([#7][#6H3,#6H2X4,a,#6X3$([#6]=[#8])]),#7$([#7]([#6H3,#6H2X4,a,#6X3$([#6]=[#8])])[#6H3,#6H2X4,a,#6X3$([#6]=[#8])])])])!$([#6H1])!$([#6][#17,#35,#53])](F)F"
    
    print("Testing Fixed EU PFAS Restriction Patterns")
    print("="*70)
    
    patt_cf3 = Chem.MolFromSmarts(eu_cf3)
    patt_cf2 = Chem.MolFromSmarts(eu_cf2)
    
    if not patt_cf3:
        print("ERROR: CF3 pattern failed to compile!")
        return
    if not patt_cf2:
        print("ERROR: CF2 pattern failed to compile!")
        return
    
    print("\nResults:")
    passed = 0
    failed = 0
    
    for smiles, name, expected_match in test_cases:
        mol = Chem.MolFromSmiles(smiles)
        canonical = Chem.MolToSmiles(mol)
        
        match_cf3 = mol.HasSubstructMatch(patt_cf3)
        match_cf2 = mol.HasSubstructMatch(patt_cf2)
        actual_match = match_cf3 or match_cf2
        
        status = "PASS" if actual_match == expected_match else "FAIL"
        if status == "PASS":
            passed += 1
        else:
            failed += 1
        
        print(f"  [{status}] {name}")
        print(f"        SMILES: {smiles} -> {canonical}")
        print(f"        Expected: {expected_match}, Actual: {actual_match} (CF3={match_cf3}, CF2={match_cf2})")
        print()
    
    print("="*70)
    print(f"Summary: {passed} passed, {failed} failed out of {len(test_cases)} tests")
    
    if failed == 0:
        print("\n✓ All tests passed! The EU patterns now correctly exclude amines.")
    else:
        print(f"\n✗ {failed} test(s) failed. Further investigation needed.")

if __name__ == "__main__":
    test_fixed_eu_patterns()
