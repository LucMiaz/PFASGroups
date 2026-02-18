"""Test SMARTS atom count precomputation."""
import json
from pathlib import Path

from HalogenGroups import HalogenGroup


def test_smarts_atom_count_precompute():
    data_dir = Path(__file__).parent.parent / 'HalogenGroups' / 'data'
    groups_file = data_dir / 'Halogen_groups_smarts.json'
    with open(groups_file, 'r') as f:
        groups_data = json.load(f)

    assert groups_data, "No halogen groups loaded"

    for group_data in groups_data[:10]:
        group = HalogenGroup(**group_data)
        if not group.smarts:
            continue
        assert len(group.smarts) == len(group.smarts_count)
        assert len(group.smarts) == len(group.smarts_str)
        for i, smarts_mol in enumerate(group.smarts):
            assert smarts_mol.GetNumAtoms() > 0
            assert group.smarts_count[i] >= 1
