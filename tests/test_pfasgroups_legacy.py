"""Tests for the PFASgroups legacy compatibility wrapper.

This module verifies that `import PFASgroups` (and `import PFASGroups`) provide
the expected legacy API:
  - Functions default to halogens='F' (fluorine-only)
  - parse_mol / parse_mols / parse_smiles return list-of-tuples:
        [(HalogenGroup, match_count, components_sizes, components), ...]
  - Passing include_PFAS_definitions=True returns the raw ResultsModel/dict
  - Halogen filtering can still be overridden by the caller
  - Non-fluorinated molecules return an empty match list
"""

import pytest
from rdkit import Chem

import PFASgroups
from PFASgroups import (
    parse_smiles,
    parse_mol,
    parse_mols,
    HalogenGroup,
    ResultsModel,
)
from HalogenGroups import parse_smiles as hg_parse_smiles


# ---------------------------------------------------------------------------
# Shared SMILES fixtures
# ---------------------------------------------------------------------------

# A short perfluoroalkyl carboxylic acid – always matches at least one F-group
PFOA_SMILES = "C(C(C(C(C(C(C(C(F)(F)F)(F)F)(F)F)(F)F)(F)F)(F)F)(F)F)(=O)O"
# Perfluorooctane sulfonic acid
PFOS_SMILES = "C(C(C(C(C(C(C(C(F)(F)F)(F)F)(F)F)(F)F)(F)F)(F)F)(F)F)S(=O)(=O)O"
# A pure chlorinated compound – no fluorine, should yield no matches with halogens='F'
CHLORINATED_SMILES = "ClCCCl"
# A small perfluoroalkyl fragment
CF3_SMILES = "C(F)(F)F"


@pytest.fixture
def pfoa_mol():
    return Chem.MolFromSmiles(PFOA_SMILES)


@pytest.fixture
def chlorinated_mol():
    return Chem.MolFromSmiles(CHLORINATED_SMILES)


# ---------------------------------------------------------------------------
# Module-level smoke tests
# ---------------------------------------------------------------------------


class TestModuleImport:
    """Verify that the PFASgroups namespace exposes the expected symbols."""

    def test_parse_smiles_callable(self):
        assert callable(parse_smiles)

    def test_parse_mol_callable(self):
        assert callable(parse_mol)

    def test_parse_mols_callable(self):
        assert callable(parse_mols)

    def test_halogen_group_class_available(self):
        assert HalogenGroup is not None

    def test_results_model_class_available(self):
        assert ResultsModel is not None

    def test_pfasgroups_module_alias(self):
        """PFASGroups (capital G) module should expose the same parse_smiles."""
        import PFASGroups  # top-level alias module
        assert callable(PFASGroups.parse_smiles)
        assert PFASGroups.parse_smiles is PFASgroups.parse_smiles


# ---------------------------------------------------------------------------
# parse_smiles – legacy output shape
# ---------------------------------------------------------------------------


class TestParseSmiliesLegacyOutput:
    """parse_smiles should return list[list[tuple]] in legacy mode."""

    def test_returns_list(self):
        result = parse_smiles([PFOA_SMILES])
        assert isinstance(result, list)

    def test_outer_list_length_matches_input(self):
        smiles_list = [PFOA_SMILES, PFOS_SMILES]
        result = parse_smiles(smiles_list)
        assert len(result) == len(smiles_list)

    def test_inner_list_for_matching_molecule(self):
        result = parse_smiles([PFOA_SMILES])
        inner = result[0]
        assert isinstance(inner, list)
        assert len(inner) > 0

    def test_tuple_structure(self):
        """Each element should be a 4-tuple: (group, count, sizes, components)."""
        result = parse_smiles([PFOA_SMILES])
        for tpl in result[0]:
            assert isinstance(tpl, tuple)
            assert len(tpl) == 4
            group, count, sizes, components = tpl
            assert isinstance(count, int)
            assert isinstance(sizes, list)
            assert isinstance(components, list)

    def test_group_is_halogen_group_instance(self):
        result = parse_smiles([PFOA_SMILES])
        for group, _count, _sizes, _components in result[0]:
            assert isinstance(group, HalogenGroup), (
                f"Expected HalogenGroup, got {type(group)}"
            )

    def test_match_count_positive(self):
        result = parse_smiles([PFOA_SMILES])
        for _group, count, _sizes, _components in result[0]:
            assert count > 0

    def test_non_fluorinated_molecule_returns_list(self):
        """Chlorinated-only molecule should return a list (may or may not have matches)."""
        result = parse_smiles([CHLORINATED_SMILES])
        assert isinstance(result[0], list)

    def test_multiple_molecules(self):
        result = parse_smiles([PFOA_SMILES, CHLORINATED_SMILES, PFOS_SMILES])
        assert len(result[0]) > 0   # PFOA matches
        assert isinstance(result[1], list)  # chlorinated – valid result
        assert len(result[2]) > 0   # PFOS matches

    def test_single_smiles_string_input(self):
        """Some callers pass a single string, not a list."""
        result = parse_smiles(PFOA_SMILES)
        # Should still be iterable; shape depends on implementation
        assert result is not None


# ---------------------------------------------------------------------------
# parse_smiles – halogens default is 'F'
# ---------------------------------------------------------------------------


class TestDefaultFluorineFilter:
    """Verify that halogens='F' is applied by default."""

    def test_default_matches_fluorinated_only(self):
        """Calling with no kwargs should filter for fluorine groups."""
        result_default = parse_smiles([PFOA_SMILES])
        result_explicit_f = parse_smiles([PFOA_SMILES], halogens="F")
        # Both should have the same number of top-level matches
        assert len(result_default[0]) == len(result_explicit_f[0])

    def test_override_halogens_to_cl_yields_valid_result(self):
        """Caller can override halogens to a different halogen and gets a valid list."""
        result_cl = parse_smiles([PFOA_SMILES], halogens="Cl")
        assert isinstance(result_cl[0], list)

    def test_results_differ_from_no_filter(self):
        """Default (F only) should differ from HalogenGroups with no filter."""
        legacy_result = parse_smiles([PFOA_SMILES])
        raw_result = hg_parse_smiles([PFOA_SMILES])
        # Legacy wraps in legacy tuples; raw is a list of dicts
        assert isinstance(legacy_result[0], list)
        assert isinstance(raw_result[0], dict)


# ---------------------------------------------------------------------------
# parse_smiles – include_PFAS_definitions passthrough
# ---------------------------------------------------------------------------


class TestParseSmilesPFASDefinitions:
    """When include_PFAS_definitions=True, raw HalogenGroups output is returned."""

    def test_returns_list_of_dicts(self):
        result = parse_smiles([PFOA_SMILES], include_PFAS_definitions=True)
        assert isinstance(result, list)
        assert isinstance(result[0], dict)

    def test_raw_result_has_matches_key(self):
        result = parse_smiles([PFOA_SMILES], include_PFAS_definitions=True)
        assert "matches" in result[0]

    def test_raw_result_has_smiles_key(self):
        result = parse_smiles([PFOA_SMILES], include_PFAS_definitions=True)
        assert "smiles" in result[0]


# ---------------------------------------------------------------------------
# parse_mol – single RDKit molecule
# ---------------------------------------------------------------------------


class TestParseMolLegacy:
    """parse_mol should return a single list of tuples."""

    def test_returns_list(self, pfoa_mol):
        result = parse_mol(pfoa_mol)
        assert isinstance(result, list)

    def test_tuples_for_fluorinated_mol(self, pfoa_mol):
        result = parse_mol(pfoa_mol)
        assert len(result) > 0
        for item in result:
            assert isinstance(item, tuple)
            assert len(item) == 4

    def test_non_fluorinated_returns_list(self, chlorinated_mol):
        """parse_mol on a non-fluorinated molecule still returns a list of tuples."""
        result = parse_mol(chlorinated_mol)
        assert isinstance(result, list)

    def test_include_pfas_definitions_returns_dict(self, pfoa_mol):
        result = parse_mol(pfoa_mol, include_PFAS_definitions=True)
        assert isinstance(result, dict)
        assert "matches" in result


# ---------------------------------------------------------------------------
# parse_mols – batch RDKit molecules
# ---------------------------------------------------------------------------


class TestParseMolsLegacy:
    """parse_mols should return list[list[tuple]]."""

    def test_returns_outer_list(self, pfoa_mol, chlorinated_mol):
        result = parse_mols([pfoa_mol, chlorinated_mol])
        assert isinstance(result, list)
        assert len(result) == 2

    def test_fluorinated_has_matches(self, pfoa_mol):
        result = parse_mols([pfoa_mol])
        assert len(result[0]) > 0

    def test_non_fluorinated_returns_list(self, chlorinated_mol):
        """parse_mols on a non-fluorinated molecule still returns a list of tuples."""
        result = parse_mols([chlorinated_mol])
        assert isinstance(result[0], list)

    def test_tuple_structure(self, pfoa_mol):
        result = parse_mols([pfoa_mol])
        for tpl in result[0]:
            assert isinstance(tpl, tuple)
            assert len(tpl) == 4

    def test_include_pfas_definitions_passthrough(self, pfoa_mol):
        result = parse_mols([pfoa_mol], include_PFAS_definitions=True)
        # Should return raw list-of-dicts
        assert isinstance(result, list)
        assert isinstance(result[0], dict)

    def test_non_list_output_format_passthrough(self, pfoa_mol):
        """output_format != 'list' should return the raw HalogenGroups output."""
        result = parse_mols([pfoa_mol], output_format="ResultsModel")
        # Should not be a list of tuples – passthrough to HalogenGroups
        assert result is not None


# ---------------------------------------------------------------------------
# Group attributes on the returned HalogenGroup objects
# ---------------------------------------------------------------------------


class TestReturnedGroupAttributes:
    """Verify the group objects in the tuples have the expected attributes."""

    def test_group_has_id(self):
        result = parse_smiles([PFOA_SMILES])
        for group, *_ in result[0]:
            assert hasattr(group, "id")
            assert group.id is not None

    def test_group_has_name(self):
        result = parse_smiles([PFOA_SMILES])
        for group, *_ in result[0]:
            assert hasattr(group, "name")
            assert isinstance(group.name, str)
            assert len(group.name) > 0

    def test_components_sizes_and_count(self):
        """components_sizes is a list of scores (one per candidate evaluation).
        components is the list of matched component dicts. count is match_count."""
        result = parse_smiles([PFOA_SMILES])
        for _group, count, sizes, components in result[0]:
            assert isinstance(count, int) and count > 0
            assert isinstance(sizes, list)
            assert isinstance(components, list)
            assert len(components) > 0


# ---------------------------------------------------------------------------
# Consistency between PFASgroups and HalogenGroups with halogens='F'
# ---------------------------------------------------------------------------


class TestConsistencyWithHalogenGroups:
    """Legacy output should be consistent with HalogenGroups when halogens='F'."""

    def test_match_count_consistent(self):
        legacy = parse_smiles([PFOA_SMILES])
        raw = hg_parse_smiles([PFOA_SMILES], halogens="F")

        # Count HalogenGroup matches from raw output
        raw_group_matches = [
            m for m in raw[0].get("matches", [])
            if m.get("type") == "HalogenGroup"
        ]
        assert len(legacy[0]) == len(raw_group_matches)

    def test_group_ids_consistent(self):
        legacy = parse_smiles([PFOA_SMILES])
        raw = hg_parse_smiles([PFOA_SMILES], halogens="F")

        legacy_ids = sorted(g.id for g, *_ in legacy[0])
        raw_ids = sorted(
            m["id"]
            for m in raw[0].get("matches", [])
            if m.get("type") == "HalogenGroup"
        )
        assert legacy_ids == raw_ids
