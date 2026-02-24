#!/usr/bin/env python3
"""Diagnose SMARTS pattern matching issues"""

from rdkit import Chem
import json

# Test molecules that aren't being detected
test_cases = {
    67: {
        'smarts': '[C$(C(=O)SC(=O)C([OH1,Oh1,O-])C(=O)[OH1,Oh1,O-])]',
        'generated_smiles': 'O=C(O)C(O)C(=O)SC(=O)C(F)(C(F)(F)F)C(F)(F)F',
        'constraints': {'gte': {'O': 5, 'S': 1}},
        'max_dist_from_comp': 0
    },
    104: {
        'smarts': '[CH2$(CS(=O)(=O)[CH2][CH2]C(=O)[OH1,Oh1,O-])]',
        'generated_smiles': 'O=C(O)CCS(=O)(=O)CC(C(F)(F)F)(C(F)(F)F)C(F)(F)C(F)(F)C(F)(F)F',
        'constraints': {'gte': {'O': 4, 'S': 1}},
        'max_dist_from_comp': 0
    },
    106: {
        'smarts': '[CH2$(CS(=O)[CH2][CH2]C(=O)[NH][CH2]([CH3])([CH3])[CH2]S(=O)(=O)[OH1,Oh1,O-])]',
        'generated_smiles': 'CC(C)(CNC(=O)CCS(=O)CC(F)(F)C(F)(F)C(F)(F)F)CS(=O)(=O)O',
        'constraints': {'gte': {'O': 5, 'S': 2, 'N': 1}},
        'max_dist_from_comp': 0
    }
}

print("🔬 SMARTS PATTERN DIAGNOSIS\n")
print("=" * 80)

for group_id, test_data in test_cases.items():
    print(f"\n📊 Group {group_id}:")
    print(f"   SMARTS: {test_data['smarts']}")
    print(f"   Test molecule: {test_data['generated_smiles']}")
    print(f"   Constraints: {test_data['constraints']}")
    print(f"   max_dist_from_comp: {test_data['max_dist_from_comp']}")
    
    mol = Chem.MolFromSmiles(test_data['generated_smiles'])
    if mol is None:
        print("   ❌ ERROR: Invalid SMILES")
        continue
    
    # Check SMARTS match
    pattern = Chem.MolFromSmarts(test_data['smarts'])
    if pattern is None:
        print("   ❌ ERROR: Invalid SMARTS pattern")
        continue
    
    matches = mol.GetSubstructMatches(pattern)
    print(f"   SMARTS matches: {len(matches)}")
    
    if matches:
        print(f"   ✅ SMARTS pattern DOES match!")
        for i, match in enumerate(matches[:3], 1):
            print(f"      Match {i}: atoms {match}")
    else:
        print(f"   ❌ SMARTS pattern DOES NOT match")
        
        # Let's break down the SMARTS pattern
        print("\n   🔍 Breaking down the SMARTS pattern:")
        
        if group_id == 67:
            # [C$(C(=O)SC(=O)C([OH1,Oh1,O-])C(=O)[OH1,Oh1,O-])]
            print("      Pattern expects: C bonded to C(=O)SC(=O)C(OH)C(=O)OH")
            print("      Generated molecule has: C(=O)SC(=O)C(O)C(=O)O attached to perfluoro chain")
            
            # Test simpler patterns
            print("\n      Testing sub-patterns:")
            sub_patterns = [
                ('C(=O)SC(=O)', 'Thioester linkage'),
                ('C(O)C(=O)O', 'Hydroxy carboxylic acid'),
                ('C(=O)SC(=O)C', 'Full thioester + carbon'),
            ]
            
            for smarts, desc in sub_patterns:
                pat = Chem.MolFromSmarts(smarts)
                if pat and mol.HasSubstructMatch(pat):
                    print(f"         ✅ {desc}: {smarts}")
                else:
                    print(f"         ❌ {desc}: {smarts}")
            
            # Check if the issue is with the anchor carbon
            print("\n      Issue: The SMARTS expects a CH2 carbon [C$(...)] directly attached")
            print("      But generated molecule has C(=O) as first atom, not a CH2")
            print("      The functional group itself is present, but missing the CH2 anchor!")
            
        elif group_id == 104:
            # [CH2$(CS(=O)(=O)[CH2][CH2]C(=O)[OH1,Oh1,O-])]
            print("      Pattern expects: CH2 bonded to CS(=O)(=O)CH2CH2C(=O)OH")
            print("      Generated molecule has: CS(=O)(=O)C where C is branched")
            
            print("\n      Testing sub-patterns:")
            sub_patterns = [
                ('S(=O)(=O)', 'Sulfonyl group'),
                ('CCS(=O)(=O)', 'Ethyl sulfonyl'),
                ('C(=O)O', 'Carboxylic acid'),
                ('CS(=O)(=O)C', 'Simple sulfonyl connection'),
                ('[CH2]S(=O)(=O)[CH2]', 'CH2-sulfonyl-CH2'),
            ]
            
            for smarts, desc in sub_patterns:
                pat = Chem.MolFromSmarts(smarts)
                if pat and mol.HasSubstructMatch(pat):
                    matches_sub = mol.GetSubstructMatches(pat)
                    print(f"         ✅ {desc}: {smarts} ({len(matches_sub)} matches)")
                else:
                    print(f"         ❌ {desc}: {smarts}")
            
            # Look at the carbon attached to S
            print("\n      🔍 Analyzing carbon attached to sulfonyl:")
            s_pattern = Chem.MolFromSmarts('[S$(S(=O)(=O)C)]')
            if s_pattern:
                s_matches = mol.GetSubstructMatches(s_pattern)
                if s_matches:
                    for s_idx in s_matches[0]:
                        s_atom = mol.GetAtomWithIdx(s_idx)
                        print(f"         Sulfur at index {s_idx}")
                        for neighbor in s_atom.GetNeighbors():
                            if neighbor.GetSymbol() == 'C':
                                print(f"         Carbon neighbor: index {neighbor.GetIdx()}, degree={neighbor.GetDegree()}, hybridization={neighbor.GetHybridization()}")
                                print(f"         Total H count: {neighbor.GetTotalNumHs()}")
                                if neighbor.GetDegree() > 2:
                                    print(f"         ❌ This carbon is branched (degree {neighbor.GetDegree()}), not CH2!")
            
            print("\n      Issue: SMARTS expects [CH2] (exactly 2 hydrogens) but molecule has branched carbon")
            
        elif group_id == 106:
            # [CH2$(CS(=O)[CH2][CH2]C(=O)[NH][CH2]([CH3])([CH3])[CH2]S(=O)(=O)[OH1,Oh1,O-])]
            print("      Pattern expects: CH2 bonded to very complex sulfinyl amido sulfonic acid")
            print("      Generated molecule has: CS(=O)CCC(=O)... structure")
            
            print("\n      Testing sub-patterns:")
            sub_patterns = [
                ('S(=O)', 'Sulfinyl group (not sulfonyl)'),
                ('CS(=O)C', 'C-sulfinyl-C'),
                ('[S$(S(=O))]', 'Sulfur with exactly one =O'),
                ('C(=O)N', 'Amide'),
                ('S(=O)(=O)O', 'Sulfonic acid'),
                ('[CH2]S(=O)[CH2]', 'CH2-sulfinyl-CH2'),
                ('NC(C)(C)CS(=O)(=O)O', 'N-dimethyl-sulfonic acid part'),
            ]
            
            for smarts, desc in sub_patterns:
                pat = Chem.MolFromSmarts(smarts)
                if pat and mol.HasSubstructMatch(pat):
                    matches_sub = mol.GetSubstructMatches(pat)
                    print(f"         ✅ {desc}: {smarts} ({len(matches_sub)} matches)")
                else:
                    print(f"         ❌ {desc}: {smarts}")
            
            print("\n      Issue: Similar to group 104 - needs [CH2] anchor but has branched carbon")

print("\n\n" + "=" * 80)
print("💡 SUMMARY OF ISSUES:\n")
print("Group 67: The functional group is present BUT the SMARTS expects a CH2 anchor carbon")
print("          that's directly attached to the C(=O)SC(=O)... structure. The generated")
print("          molecule doesn't have this anchor - it starts directly with C(=O).")
print()
print("Group 104: The SMARTS pattern requires [CH2] (exactly 2 hydrogens) attached to the")
print("           sulfonyl group, but the generated molecules have branched carbons due to")
print("           perfluoro branching. Need molecules with linear -CH2-S(=O)(=O)- structure.")
print()
print("Group 106: Same issue as 104 - needs [CH2] anchor but molecules have branched carbons.")
print("\n" + "=" * 80)
print("🔧 FIXES NEEDED:\n")
print("1. Group 67: Generation SMILES should create an anchor: 'CC(=O)SC(=O)C(O)C(=O)O'")
print("             (add a CH2 before the thioester)")
print()
print("2. Group 104: Ensure linear structure: 'O=C(O)CCS(=O)(=O)C' without branching on S-C")
print("              The current generation creates branches which prevent [CH2] match")
print()
print("3. Group 106: Similar - need linear 'CS(=O)CCC(=O)NCC(C)(C)CS(=O)(=O)O' structure")
print("              without branching at the sulfinyl carbon")
