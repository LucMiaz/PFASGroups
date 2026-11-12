"""Fluorine ratio of PFAS definitions by mass and by atom count.

Checks :class:`PFASGroups.PFASDefinitionModel.PFASDefinition` against ratios
worked out by hand from molecular formulas, and checks that PFASSTRUCTv5 uses
the atom count of Gaines et al. (2023, doi:10.1002/ieam.4735): at least 30% of
the non-hydrogen atoms are fluorine.
"""
from __future__ import annotations

import json
from pathlib import Path

import pytest

from PFASGroups.PFASDefinitionModel import DEFAULT_RATIO_BASIS, PFASDefinition

DEFINITIONS_JSON = Path(__file__).parent.parent / "PFASGroups" / "data" / "PFAS_definitions_smarts.json"

#: (name, SMILES, count ratio without H, mass ratio without H). The count ratio is
#: F atoms / non-hydrogen atoms; the mass ratio uses C 12.011, O 15.999, F 18.998.
MOLECULES = [
    # C3H3F3O2: 3 F of 8 heavy atoms; 56.99 / (36.03 + 32.00 + 56.99)
    ("3,3,3-trifluoropropanoic acid", "OC(=O)CC(F)(F)F", 3 / 8, 0.4559),
    # C9H15F3O2: 3 F of 14 heavy atoms; 56.99 / (108.10 + 32.00 + 56.99)
    ("2-propyl-2-(trifluoromethyl)pentanoic acid", "CCCC(CCC)(C(=O)O)C(F)(F)F", 3 / 14, 0.2892),
    # C8H5F3O2: 3 F of 13 heavy atoms; 56.99 / (96.09 + 32.00 + 56.99)
    ("4-(trifluoromethyl)benzoic acid", "OC(=O)c1ccc(cc1)C(F)(F)F", 3 / 13, 0.3079),
]


def _ratio_only_definition(**kwargs) -> PFASDefinition:
    """A definition with no SMARTS, so that only the fluorine ratio decides."""
    return PFASDefinition(id=99, name="ratio only", smarts=[], fluorineRatio=0.3,
                          description="test", includeHydrogen=False, **kwargs)


def _pfasstructv5() -> PFASDefinition:
    """PFASSTRUCTv5 as shipped in PFAS_definitions_smarts.json."""
    entries = json.loads(DEFINITIONS_JSON.read_text(encoding="utf-8"))
    return PFASDefinition(**next(e for e in entries if e["name"] == "PFASSTRUCTv5"))


@pytest.mark.parametrize("name, smiles, count_ratio, mass_ratio", MOLECULES)
def test_ratio_matches_hand_computed_values(name, smiles, count_ratio, mass_ratio):
    """Both bases reproduce the ratios computed from the formula."""
    definition = _ratio_only_definition()
    by_count = definition.match_details(smiles, ratio_basis="count")
    by_mass = definition.match_details(smiles, ratio_basis="mass")
    assert by_count["fluorine_ratio"] == pytest.approx(count_ratio, abs=1e-4)
    assert by_mass["fluorine_ratio"] == pytest.approx(mass_ratio, abs=1e-3)
    assert by_count["fluorine_ratio_basis"] == "count"
    assert by_mass["fluorine_ratio_basis"] == "mass"


def test_default_basis_is_mass():
    """A definition without ratioBasis keeps the mass ratio of versions 3.4.x."""
    assert DEFAULT_RATIO_BASIS == "mass"
    assert _ratio_only_definition().ratio_basis == "mass"


def test_ratio_basis_is_read_from_both_spellings():
    """The JSON key (ratioBasis) and the Python keyword (ratio_basis) set the basis."""
    assert _ratio_only_definition(ratioBasis="count").ratio_basis == "count"
    assert _ratio_only_definition(ratio_basis="COUNT").ratio_basis == "count"


def test_unknown_basis_raises():
    """A basis other than mass or count is refused at construction and at call time."""
    with pytest.raises(ValueError, match="ratio_basis"):
        _ratio_only_definition(ratioBasis="volume")
    with pytest.raises(ValueError, match="ratio_basis"):
        _ratio_only_definition().applies_to_molecule("FC(F)(F)C(F)(F)F", ratio_basis="volume")


def test_mass_and_count_disagree_between_the_two_thresholds():
    """4-(Trifluoromethyl)benzoic acid passes 0.3 by mass (0.31) but not by count (0.23)."""
    smiles = "OC(=O)c1ccc(cc1)C(F)(F)F"
    assert _ratio_only_definition(ratioBasis="mass").applies_to_molecule(smiles)
    assert not _ratio_only_definition(ratioBasis="count").applies_to_molecule(smiles)


def test_pfasstructv5_counts_atoms():
    """PFASSTRUCTv5 follows Gaines et al. (2023): 30% of the non-hydrogen atoms."""
    definition = _pfasstructv5()
    assert definition.ratio_basis == "count"
    assert definition.includeHydrogen is False
    # 3 of 8 non-hydrogen atoms are fluorine: inside the definition by ratio alone.
    assert definition.applies_to_molecule("OC(=O)CC(F)(F)F")
    # No PFASSTRUCTv5 substructure and 3 of 13 non-hydrogen atoms: outside, although
    # the mass ratio (0.31) is above 0.3.
    details = definition.match_details("OC(=O)c1ccc(cc1)C(F)(F)F")
    assert not details["smarts_match"]
    assert not details["applies"]
    # Benzotrifluoride, C7H5F3: exactly 3 of 10, and the threshold is inclusive.
    assert definition.applies_to_molecule("FC(F)(F)c1ccccc1")
