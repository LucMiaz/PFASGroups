"""Understand what SMARTS patterns match."""
from rdkit import Chem


def test_smarts_matching_patterns():
    smarts_patterns = {
        "Carboxylic acid (current)": "[#6$([#6][#6](=O)([OH1,Oh1,O-]))]",
        "Carboxylic acid (full)": "C(=O)O",
        "Sulfonic acid (current)": "[#6$([#6][#16](=O)(=O)[OH1,Oh1,O-])]",
        "Sulfonic acid (full)": "S(=O)(=O)O",
    }

    test_mol = Chem.MolFromSmiles("CC(=O)O")
    assert test_mol is not None
    test_mol = Chem.AddHs(test_mol)

    for name, smarts_str in smarts_patterns.items():
        smarts = Chem.MolFromSmarts(smarts_str)
        assert smarts is not None, f"Failed to parse SMARTS: {name}"
        matches = test_mol.GetSubstructMatches(smarts)
        assert isinstance(matches, tuple)

    full_matches = test_mol.GetSubstructMatches(Chem.MolFromSmarts("C(=O)O"))
    assert len(full_matches) > 0
