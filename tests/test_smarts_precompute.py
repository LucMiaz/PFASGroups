"""Test SMARTS precomputation in HalogenGroup."""
from PFASGroups.parser import parse_smiles


def test_smarts_precompute_pfoa():
    pfoa_smiles = "C(C(C(C(C(C(C(C(=O)O)(F)F)(F)F)(F)F)(F)F)(F)F)(F)F)(F)(F)F"
    results = parse_smiles([pfoa_smiles], bycomponent=True)
    assert results
    result = results[0]
    assert 'formula' in result
    assert 'matches' in result
