"""Tests for ResultsModel helpers."""

import json

from PFASGroups.PFASEmbeddings import ResultsModel


def test_results_model_helpers():
    raw_results = [
        {
            "smiles": "C(F)(F)C(F)(F)C(=O)O",
            "matches": [
                {
                    "match_id": "G1",
                    "id": 1,
                    "group_name": "Test Group",
                    "match_count": 1,
                    "type": "HalogenGroup",
                    "components": [
                        {"component": [1, 2, 3], "SMARTS": "Perfluoroalkyl"},
                        {"component": [3, 4], "SMARTS": "Perfluoroalkyl"},
                    ],
                },
                {
                    "match_id": "D1",
                    "id": 1,
                    "definition_name": "Definition",
                    "type": "PFASdefinition",
                },
            ],
        }
    ]

    results = ResultsModel.from_raw(raw_results)
    assert isinstance(results, ResultsModel)
    assert results[0].smiles == "C(F)(F)C(F)(F)C(=O)O"

    matches = list(results.iter_group_matches(group_id=1))
    assert len(matches) == 1
    mol_result, match = matches[0]
    assert mol_result.smiles
    assert match.group_name == "Test Group"
    assert match.is_group
    assert not match.is_definition

    atoms = mol_result.collect_component_atoms(group_id=1)
    assert atoms == [1, 2, 3, 4]

    # Ensure dict-like results remain JSON serializable
    json.dumps(results)
