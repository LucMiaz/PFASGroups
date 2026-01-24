#!/usr/bin/env python3

from rdkit import Chem

def test_smarts_patterns():
    """Test various SMARTS patterns against FC(F)(F)N to understand matching."""
    
    # Test molecule: FC(F)(F)N (trifluoromethylamine)
    mol = Chem.MolFromSmiles('FC(F)(F)N')
    print(f"Testing molecule: {Chem.MolToSmiles(mol)} (trifluoromethylamine)")
    print("="*60)
    
    # EU PFAS Restriction CF3 pattern (Pattern 1)
    eu_cf3_pattern = "[#6X4!$([#6][#8H1,#8$([#8][#6H3,#6H2X4,a,#6X3$([#6]=[#8])])]);#6X4!$([#6][#7H2,#7H1$([#7][#6H3,#6H2X4,a,#6X3$([#6]=[#8])]),#7$([#7]([#6H3,#6H2X4,a,#6X3$([#6]=[#8])])[#6H3,#6H2X4,a,#6X3$([#6]=[#8])])]);#6X4!$([#6H1]);#6X4!$([#6][#17,#35,#53])](F)(F)F"
    
    try:
        eu_patt = Chem.MolFromSmarts(eu_cf3_pattern)
        if eu_patt:
            matches = mol.HasSubstructMatch(eu_patt)
            print(f"EU CF3 pattern matches: {matches}")
        else:
            print("EU CF3 pattern failed to compile")
    except Exception as e:
        print(f"Error with EU pattern: {e}")
    
    print()
    
    # Test individual components
    patterns_to_test = [
        ("Simple CF3", "[#6](F)(F)F"),
        ("CF3 not bonded to N", "[#6!$([#6]N)](F)(F)F"),
        ("CF3 not bonded to NH2", "[#6!$([#6][#7H2])](F)(F)F"),
        ("C bonded to N", "[#6][#7]"),
        ("Primary amine", "[#7H2]"),
        ("OECD Pattern", "[#6X4!$([#6H1]);#6X4!$([#6][#17,#35,#53])](F)F"),
    ]
    
    for name, pattern in patterns_to_test:
        try:
            patt = Chem.MolFromSmarts(pattern)
            if patt:
                matches = mol.HasSubstructMatch(patt)
                print(f"{name:25} | {pattern:35} | Matches: {matches}")
            else:
                print(f"{name:25} | {pattern:35} | Failed to compile")
        except Exception as e:
            print(f"{name:25} | {pattern:35} | Error: {e}")
    
    print()
    print("EU CF3 Pattern Analysis:")
    print("="*60)
    
    # Break down the EU pattern
    eu_parts = [
        ("#6X4 (carbon with 4 bonds)", "[#6X4]"),
        ("Not bonded to OH/ethers", "!$([#6][#8H1,#8$([#8][#6H3,#6H2X4,a,#6X3$([#6]=[#8])])])"),
        ("Not bonded to amines", "!$([#6][#7H2,#7H1$([#7][#6H3,#6H2X4,a,#6X3$([#6]=[#8])]),#7$([#7]([#6H3,#6H2X4,a,#6X3$([#6]=[#8])])[#6H3,#6H2X4,a,#6X3$([#6]=[#8])])])"),
        ("Not bonded to H", "!$([#6H1])"),
        ("Not bonded to Cl/Br/I", "!$([#6][#17,#35,#53])"),
        ("Has three F", "(F)(F)F")
    ]
    
    print("The EU pattern should exclude CF3 bonded to:")
    print("- OH groups and ethers")
    print("- Primary and secondary amines")  
    print("- H atoms")
    print("- Cl, Br, I")
    print()
    
    # Test the amine exclusion part specifically
    amine_exclusions = [
        ("Primary amine [#7H2]", "[#6][#7H2]"),
        ("Secondary amine", "[#6][#7H1$([#7][#6H3,#6H2X4,a,#6X3$([#6]=[#8])])]"),
        ("Complete amine exclusion", "[#6][#7H2,#7H1$([#7][#6H3,#6H2X4,a,#6X3$([#6]=[#8])]),#7$([#7]([#6H3,#6H2X4,a,#6X3$([#6]=[#8])])[#6H3,#6H2X4,a,#6X3$([#6]=[#8])])]")
    ]
    
    for name, pattern in amine_exclusions:
        try:
            patt = Chem.MolFromSmarts(pattern)
            if patt:
                matches = mol.HasSubstructMatch(patt)
                print(f"{name:25} | Matches: {matches}")
                if matches:
                    match_atoms = mol.GetSubstructMatch(patt)
                    print(f"  Matched atoms: {match_atoms}")
            else:
                print(f"{name:25} | Failed to compile")
        except Exception as e:
            print(f"{name:25} | Error: {e}")

if __name__ == "__main__":
    test_smarts_patterns()