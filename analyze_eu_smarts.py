#!/usr/bin/env python3

from rdkit import Chem

def analyze_eu_cf3_pattern():
    """Analyze why EU CF3 pattern matches FC(F)(F)N when it shouldn't."""
    
    mol = Chem.MolFromSmiles('FC(F)(F)N')
    print(f"Testing molecule: {Chem.MolToSmiles(mol)} (trifluoromethylamine)")
    print("="*70)
    
    # The EU CF3 pattern
    eu_cf3 = "[#6X4!$([#6][#8H1,#8$([#8][#6H3,#6H2X4,a,#6X3$([#6]=[#8])])]);#6X4!$([#6][#7H2,#7H1$([#7][#6H3,#6H2X4,a,#6X3$([#6]=[#8])]),#7$([#7]([#6H3,#6H2X4,a,#6X3$([#6]=[#8])])[#6H3,#6H2X4,a,#6X3$([#6]=[#8])])]);#6X4!$([#6H1]);#6X4!$([#6][#17,#35,#53])](F)(F)F"
    
    print("\nEU CF3 SMARTS pattern breakdown:")
    print("- #6X4 = carbon with 4 bonds")
    print("- !$([#6][#8H1,...]) = NOT bonded to OH/ethers")
    print("- !$([#6][#7H2,...]) = NOT bonded to amines")
    print("- !$([#6H1]) = NOT bonded to H")
    print("- !$([#6][#17,#35,#53]) = NOT bonded to Cl/Br/I")
    print("- (F)(F)F = has three fluorines")
    print()
    
    patt = Chem.MolFromSmarts(eu_cf3)
    matches = mol.HasSubstructMatch(patt)
    print(f"Full EU CF3 pattern matches: {matches}")
    
    if matches:
        match_atoms = mol.GetSubstructMatch(patt)
        print(f"Matched atoms: {match_atoms}")
    
    print("\n" + "="*70)
    print("Testing individual exclusions:")
    print("="*70)
    
    # Test each exclusion separately
    tests = [
        ("Base CF3", "[#6X4](F)(F)F"),
        ("Exclude H", "[#6X4!$([#6H1])](F)(F)F"),
        ("Exclude halogens", "[#6X4!$([#6][#17,#35,#53])](F)(F)F"),
        ("Exclude OH", "[#6X4!$([#6][#8H1])](F)(F)F"),
        ("Exclude primary amine", "[#6X4!$([#6][#7H2])](F)(F)F"),
        ("Exclude amine (NH2 only)", "[#6X4!$([#6][NH2])](F)(F)F"),
        ("Exclude nitrogen", "[#6X4!$([#6][#7])](F)(F)F"),
        ("Exclude N (simple)", "[#6X4!$([#6]N)](F)(F)F"),
    ]
    
    for name, pattern in tests:
        try:
            patt = Chem.MolFromSmarts(pattern)
            if patt:
                matches = mol.HasSubstructMatch(patt)
                status = "MATCHES" if matches else "no match"
                print(f"{name:30} | {status}")
        except Exception as e:
            print(f"{name:30} | ERROR: {e}")
    
    print("\n" + "="*70)
    print("Analyzing the amine exclusion part in detail:")
    print("="*70)
    
    # The amine exclusion from the EU pattern
    amine_exclusion = "!$([#6][#7H2,#7H1$([#7][#6H3,#6H2X4,a,#6X3$([#6]=[#8])]),#7$([#7]([#6H3,#6H2X4,a,#6X3$([#6]=[#8])])[#6H3,#6H2X4,a,#6X3$([#6]=[#8])])])"
    
    print("\nThe amine exclusion pattern attempts to exclude:")
    print("1. #7H2 = primary amines (NH2)")
    print("2. #7H1$(...) = secondary amines (NHR)")
    print("3. #7$(...) = tertiary amines (NR2)")
    print()
    
    # Test if carbon is bonded to nitrogen
    carbon_n = "[#6][#7]"
    patt_cn = Chem.MolFromSmarts(carbon_n)
    print(f"Is carbon bonded to nitrogen: {mol.HasSubstructMatch(patt_cn)}")
    
    # Test if nitrogen is NH2
    nh2 = "[NH2]"
    patt_nh2 = Chem.MolFromSmarts(nh2)
    print(f"Does molecule have NH2: {mol.HasSubstructMatch(patt_nh2)}")
    
    # Test if carbon is bonded to NH2
    c_nh2 = "[#6][NH2]"
    patt_c_nh2 = Chem.MolFromSmarts(c_nh2)
    print(f"Is carbon bonded to NH2: {mol.HasSubstructMatch(patt_c_nh2)}")
    
    # Test the exact amine exclusion from EU pattern
    c_nh2_full = "[#6][#7H2]"
    patt_full = Chem.MolFromSmarts(c_nh2_full)
    print(f"EU pattern's C-NH2 test: {mol.HasSubstructMatch(patt_full)}")
    
    print("\n" + "="*70)
    print("THE PROBLEM:")
    print("="*70)
    
    # The issue is likely in how the semicolons work
    print("\nThe EU pattern uses semicolons (;) to chain multiple conditions.")
    print("In SMARTS, ';' means AND for recursive queries.")
    print()
    print("Let's test if the amine exclusion is even being evaluated:")
    
    # Simplified version that should work
    working_pattern = "[#6!$([#6][#7])](F)(F)F"
    patt_work = Chem.MolFromSmarts(working_pattern)
    print(f"\nSimplified exclusion [#6!$([#6][#7])](F)(F)F: {mol.HasSubstructMatch(patt_work)}")
    
    # Test with explicit X4
    working_pattern2 = "[#6X4!$([#6][#7])](F)(F)F"
    patt_work2 = Chem.MolFromSmarts(working_pattern2)
    print(f"With X4 [#6X4!$([#6][#7])](F)(F)F: {mol.HasSubstructMatch(patt_work2)}")

if __name__ == "__main__":
    analyze_eu_cf3_pattern()