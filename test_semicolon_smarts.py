#!/usr/bin/env python3

from rdkit import Chem

def test_semicolon_smarts():
    """Test if the semicolon chaining in EU SMARTS is the issue."""
    
    mol = Chem.MolFromSmiles('FC(F)(F)N')
    print(f"Testing: {Chem.MolToSmiles(mol)}")
    print("="*70)
    
    # The EU CF3 pattern uses semicolons to chain recursive queries
    # This should mean "carbon that satisfies ALL these conditions"
    
    # Let's test step by step
    
    print("\nStep 1: Base pattern")
    patt1 = "[#6X4](F)(F)F"
    mol1 = Chem.MolFromSmarts(patt1)
    print(f"[#6X4](F)(F)F: {mol.HasSubstructMatch(mol1)}")
    
    print("\nStep 2: Add hydrogen exclusion with semicolon")
    patt2 = "[#6X4;#6X4!$([#6H1])](F)(F)F"
    mol2 = Chem.MolFromSmarts(patt2)
    print(f"[#6X4;#6X4!$([#6H1])](F)(F)F: {mol.HasSubstructMatch(mol2)}")
    
    print("\nStep 3: Add amine exclusion")
    patt3 = "[#6X4;#6X4!$([#6][#7H2]))](F)(F)F"
    try:
        mol3 = Chem.MolFromSmarts(patt3)
        if mol3:
            print(f"[#6X4;#6X4!$([#6][#7H2]))](F)(F)F: {mol.HasSubstructMatch(mol3)}")
        else:
            print("Pattern failed to compile!")
    except Exception as e:
        print(f"ERROR: {e}")
    
    print("\n" + "="*70)
    print("POTENTIAL ISSUE: Semicolon usage in recursive SMARTS")
    print("="*70)
    
    # The EU pattern has this structure:
    # [#6X4!$(...);#6X4!$(...);#6X4!$(...);#6X4!$(...)](F)(F)F
    # 
    # This means the carbon must match #6X4 AND NOT match pattern1 AND NOT match pattern2, etc.
    # But the way it's written might not work correctly.
    
    print("\nTesting simplified amine exclusion:")
    simple_exclude = "[#6X4!$([#6][#7])](F)(F)F"
    patt_simple = Chem.MolFromSmarts(simple_exclude)
    result = mol.HasSubstructMatch(patt_simple)
    print(f"[#6X4!$([#6][#7])](F)(F)F: {result}")
    print("This correctly returns False (no match) because C is bonded to N")
    
    print("\n" + "="*70)
    print("Testing the ACTUAL EU pattern structure:")
    print("="*70)
    
    # Let's try to understand if the multiple semicolons are the problem
    # The pattern is: [A;B;C;D](F)(F)F
    # Where A, B, C, D are all conditions that must be true
    
    test_patterns = [
        ("Single condition", "[#6X4](F)(F)F"),
        ("With one negation", "[#6X4!$([#6H1])](F)(F)F"),
        ("Two conditions with ;", "[#6X4;!$([#6H1])](F)(F)F"),
        ("Two negations with ;", "[#6X4!$([#6H1]);#6X4!$([#6][#7])](F)(F)F"),
    ]
    
    for name, pattern in test_patterns:
        try:
            patt = Chem.MolFromSmarts(pattern)
            if patt:
                result = mol.HasSubstructMatch(patt)
                print(f"\n{name}:")
                print(f"  Pattern: {pattern}")
                print(f"  Result: {result}")
            else:
                print(f"\n{name}: FAILED TO COMPILE")
        except Exception as e:
            print(f"\n{name}: ERROR - {e}")

if __name__ == "__main__":
    test_semicolon_smarts()