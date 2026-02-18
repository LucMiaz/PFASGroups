"""Test comprehensive graph metrics in ComponentsSolver."""

from rdkit import Chem
from HalogenGroups.parser import parse_smiles


def test_comprehensive_metrics():
    test_cases = [
        {
            'name': 'PFOA (linear perfluoroalkyl)',
            'smiles': 'C(=O)(O)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)F',
        },
        {
            'name': 'Branched PFAS',
            'smiles': 'C(C(F)(F)F)(C(F)(F)F)C(F)(F)C(F)(F)C(F)(F)C(=O)O',
        },
        {
            'name': 'PFOS (sulfonic acid)',
            'smiles': 'C(C(C(C(C(C(C(C(F)(F)S(=O)(=O)O)(F)F)(F)F)(F)F)(F)F)(F)F)(F)F)(F)F',
        },
        {
            'name': '6:2 FTOH (terminal OH)',
            'smiles': 'C(CCO)(F)(C(C(C(C(C(F)(F)F)(F)F)(F)F)(F)F)(F)F)F',
        },
    ]

    for case in test_cases:
        result = parse_smiles(case['smiles'], bycomponent=True, output_format='dict')
        assert result and isinstance(result, list), f"No results for {case['name']}"
        mol_result = result[0]
        matches = mol_result.get('matches', [])
        assert isinstance(matches, list)

        found_component = False
        for match in matches:
            components = match.get('components', [])
            for comp_data in components:
                found_component = True
                assert comp_data.get('size', 0) > 0
                assert 'smarts_centrality' in comp_data
                assert 'diameter' in comp_data
                assert 'radius' in comp_data
                assert 'barycenter' in comp_data
                assert 'min_dist_to_center' in comp_data
            if found_component:
                break
        assert found_component, f"No components found for {case['name']}"
