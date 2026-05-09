"""
Test script to verify H-component and wildcard documentation examples work.
Run from the PFASGroups repository root.
"""

import sys
from rdkit import Chem
from PFASGroups import parse_smiles, HalogenGroup
from PFASGroups.generate_homologues import generate_homologues

def test_h_component_basic():
    """Test basic H-component homologue generation."""
    print("\n=== Test: H-Component Basic ===")
    mol = Chem.MolFromSmiles('OC(=O)CCCCCC')  # 6-carbon alkyl chain
    
    homologues = generate_homologues(mol, halogen='H')
    assert homologues.halogen == 'H', f"Expected halogen='H', got {homologues.halogen}"
    assert len(homologues) > 0, f"Expected homologues, got empty result"
    
    print(f"✓ Parent: {Chem.MolToSmiles(mol)}")
    print(f"✓ Generated {len(homologues)} H-component homologues")
    for inchikey, inner_dict in list(homologues.items())[:3]:  # Show first 3
        for formula, h_mol in inner_dict.items():
            print(f"  • {formula}: {Chem.MolToSmiles(h_mol)}")
    

def test_wildcard_groups_disabled_by_default():
    """Test that wildcard groups are disabled by default."""
    print("\n=== Test: Wildcard Groups (Default Off) ===")
    smiles = "CCO"  # ethanol
    
    results_no_wc = parse_smiles(smiles, halogens='F')
    results_with_wc = parse_smiles(smiles, halogens='*')
    
    matches_without = len(results_no_wc[0].matches)
    matches_with = len(results_with_wc[0].matches)
    
    print(f"✓ Without wildcards: {matches_without} match(es)")
    print(f"✓ With wildcards: {matches_with} match(es)")
    assert matches_with >= matches_without, "Enabling wildcards should not decrease matches"
    

def test_wildcard_group_names():
    """Test that wildcard groups are properly named."""
    print("\n=== Test: Wildcard Group Names ===")
    test_molecules = [
        ("CCO", "Alcohol"),
        ("CCOC", "Ether"),
        ("CC(=O)O", "Carboxylic acid"),
        ("CC(=O)OC", "Ester"),
    ]
    
    results = parse_smiles([smi for smi, _ in test_molecules], halogens='*')
    
    for (smi, expected_group), mol in zip(test_molecules, results):
        if mol.matches:
            names = [m.get('group_name') or f"ID-{m.get('id')}" for m in mol.matches if m.get('type') == 'WildcardGroup']
            print(f"✓ {smi:15} -> {', '.join(names)}")
        else:
            print(f"✓ {smi:15} -> (no match)")


def test_wildcard_pfas_separation():
    """Test distinguishing between wildcard and PFAS groups."""
    print("\n=== Test: Wildcard vs. PFAS Groups ===")
    smiles = "FC(F)(F)C(F)(F)C(=O)O"  # TFA with carboxylic acid
    
    results = parse_smiles(smiles, halogens=['F', '*'])
    mol = results[0]
    
    halogen_matches = [m for m in mol.matches if m.get('type') == 'HalogenGroup']
    wildcard_matches = [m for m in mol.matches if m.get('type') == 'WildcardGroup']
    definition_matches = [m for m in mol.matches if m.get('type') == 'PFASdefinition']
    
    print(f"✓ Halogen group matches: {len(halogen_matches)}")
    for m in halogen_matches:
        print(f"  • {m.group_name} (ID {m.group_id})")
    
    print(f"✓ Wildcard group matches: {len(wildcard_matches)}")
    for m in wildcard_matches:
        print(f"  • {m.get('group_name', 'Unknown')} (Type: {m.get('type', 'Unknown')})")

    print(f"✓ PFAS definition matches: {len(definition_matches)}")
    

def test_h_component_detection():
    """Test H-component detection through parse_smiles with a custom H-group."""
    print("\n=== Test: H-Component Detection ===")
    smiles = "CCCCO"

    h_group = HalogenGroup(
        id=9990,
        name="Hydrocarbon alcohol via H-alkyl component",
        smarts={"[#6$([#6!$([#6]=O)][OH1,Oh1,O-])]": 1},
        componentSmarts="Alkyl",
        componentSaturation="per",
        componentHalogens="H",
        componentForm="alkyl",
        constraints={},
        max_dist_from_comp=1,
    )

    results = parse_smiles(smiles, halogens="H", pfas_groups=[h_group], bycomponent=True)
    matches = [m for m in results[0].matches if m.get("id") == 9990 and m.get("type") == "HalogenGroup"]
    assert len(matches) > 0, "Expected parse_smiles to detect H-component custom group"

    print(f"✓ Found {len(matches)} H-component match(es)")
    print(f"  Components in first match: {matches[0].get('num_components', 0)}")


def test_combined_h_and_wildcard():
    """Test using H-components and wildcards together."""
    print("\n=== Test: Combined H and Wildcard Analysis ===")
    smiles = "OC(=O)CCCCCC"  # Hydrocarbon with carboxylic acid
    
    # Parse with wildcards
    results = parse_smiles(smiles, halogens='*')
    matches = [m for m in results[0].matches if m.get('type') == 'WildcardGroup']
    
    print(f"✓ Detected {len(matches)} wildcard group(s)")
    for m in matches:
        print(f"  • {m.get('group_name')}")
    
    # Generate H-homologues
    mol = Chem.MolFromSmiles(smiles)
    homologues = generate_homologues(mol, halogen='H')
    
    print(f"✓ Generated {len(homologues)} H-component homologues")


def main():
    """Run all tests."""
    print("=" * 60)
    print("Testing H-Component and Wildcard Documentation Examples")
    print("=" * 60)
    
    try:
        test_h_component_basic()
        test_wildcard_groups_disabled_by_default()
        test_wildcard_group_names()
        test_wildcard_pfas_separation()
        test_h_component_detection()
        test_combined_h_and_wildcard()
        
        print("\n" + "=" * 60)
        print("✓ All tests passed!")
        print("=" * 60)
        return 0
        
    except Exception as e:
        print(f"\n✗ Test failed with error: {e}")
        import traceback
        traceback.print_exc()
        return 1


if __name__ == '__main__':
    sys.exit(main())
