import pytest
from rdkit import Chem

from PFASGroups import HalogenGroup, parse_smiles
from PFASGroups.parser import _match_wildcard_groups_in_mol


_EXCLUDED_WILDCARD_IDS = {34, 35, 37, 38, 44, 45}


def _entry_matches(result_entry):
    return result_entry.get("matches", [])


def _type_matches(result_entry, match_type):
    return [m for m in _entry_matches(result_entry) if m.get("type") == match_type]


def _build_h_custom_group(group_id):
    return HalogenGroup(
        id=group_id,
        name="Hydrocarbon alcohol via H-alkyl component",
        smarts={"[#6$([#6!$([#6]=O)][OH1,Oh1,O-])]": 1},
        componentSmarts="Alkyl",
        componentSaturation="per",
        componentHalogens="H",
        componentForm="alkyl",
        constraints={},
        max_dist_from_comp=1,
    )


class TestWildcardHelper:
    def test_match_wildcard_groups_none_input_returns_empty(self):
        assert _match_wildcard_groups_in_mol(None) == []

    def test_match_wildcard_groups_schema_and_id_filters(self):
        mol = Chem.MolFromSmiles("CCO")
        assert mol is not None

        hits = _match_wildcard_groups_in_mol(mol)
        assert isinstance(hits, list)
        assert hits, "Expected at least one wildcard match for ethanol"

        ids = []
        for hit in hits:
            assert isinstance(hit, dict)
            assert "id" in hit and "name" in hit
            assert isinstance(hit["id"], int)
            assert isinstance(hit["name"], str)
            assert 29 <= hit["id"] <= 76
            assert hit["id"] not in _EXCLUDED_WILDCARD_IDS
            ids.append(hit["id"])

        assert len(ids) == len(set(ids)), "Wildcard IDs should not repeat for one molecule"


class TestWildcardParseIntegration:
    def test_parse_smiles_wildcard_groups_adds_wildcard_matches(self):
        smiles = "CCO"
        res_off = parse_smiles(smiles, halogens='F')
        res_on = parse_smiles(smiles, halogens='*')

        assert len(res_off) == 1
        assert len(res_on) == 1

        w_off = _type_matches(res_off[0], "WildcardGroup")
        w_on = _type_matches(res_on[0], "WildcardGroup")

        assert w_off == []
        assert len(w_on) >= 1

        for w in w_on:
            assert str(w.get("match_id", "")).startswith("W")
            assert w.get("components") == []
            assert w.get("components_sizes") == []
            assert w.get("num_components") == 0

    def test_parse_smiles_wildcard_toggle_does_not_change_halogen_group_hits(self):
        # PFOA-like input should keep the same HalogenGroup matches irrespective of wildcard toggle.
        smiles = "OC(=O)C(F)(F)C(F)(F)C(F)(F)C(F)(F)F"
        res_off = parse_smiles(smiles, halogens='F')
        res_on = parse_smiles(smiles, halogens=['F', '*'])

        assert len(res_off) == 1
        assert len(res_on) == 1

        g_off = {
            m.get("match_id")
            for m in _type_matches(res_off[0], "HalogenGroup")
        }
        g_on = {
            m.get("match_id")
            for m in _type_matches(res_on[0], "HalogenGroup")
        }

        assert g_off == g_on

    def test_wildcard_expected_group_names_for_reference_panel(self):
        panel = [
            ("CCO", "alcohol"),
            ("CCOC", "ether"),
            ("CC(=O)O", "carboxylic acid"),
            ("CC(=O)OC", "ester"),
        ]

        results = parse_smiles([s for s, _ in panel], halogens='*')
        assert len(results) == len(panel)

        for (smiles, expected_name), entry in zip(panel, results):
            wildcard_names = {
                str(m.get("group_name", "")).lower()
                for m in _type_matches(entry, "WildcardGroup")
            }
            assert expected_name in wildcard_names, (
                f"Expected wildcard group '{expected_name}' for {smiles}, "
                f"got {sorted(wildcard_names)}"
            )

    def test_telomer_groups_not_returned_for_h_or_wildcard_only(self):
        smiles = "OCCCCCCC(F)(F)C(F)(F)C(F)(F)C(F)(F)F"

        res_h = parse_smiles(smiles, halogens='H')
        res_wc = parse_smiles(smiles, halogens='*')

        names_h = [str(m.get('group_name', '')).lower() for m in _type_matches(res_h[0], "HalogenGroup")]
        names_wc = [str(m.get('group_name', '')).lower() for m in _type_matches(res_wc[0], "HalogenGroup")]

        assert not any('telomer' in n for n in names_h)
        assert not any('telomer' in n for n in names_wc)

    def test_expected_results_with_and_without_h_and_wildcard(self):
        smiles = "CCCCO"
        custom_groups = [_build_h_custom_group(997)]

        res_base = parse_smiles(smiles, pfas_groups=custom_groups, halogens='F')
        res_h = parse_smiles(smiles, pfas_groups=custom_groups, halogens='H')
        res_wc = parse_smiles(smiles, pfas_groups=custom_groups, halogens='*')
        res_both = parse_smiles(smiles, pfas_groups=custom_groups, halogens=['H', '*'])

        base_h_hits = [m for m in _type_matches(res_base[0], "HalogenGroup") if m.get("id") == 997]
        h_hits = [m for m in _type_matches(res_h[0], "HalogenGroup") if m.get("id") == 997]
        wc_hits = _type_matches(res_wc[0], "WildcardGroup")
        both_h_hits = [m for m in _type_matches(res_both[0], "HalogenGroup") if m.get("id") == 997]
        both_wc_hits = _type_matches(res_both[0], "WildcardGroup")

        # Baseline F-only: no custom H group and no wildcard groups.
        assert base_h_hits == []
        assert _type_matches(res_base[0], "WildcardGroup") == []

        # H-only: custom H group appears.
        assert len(h_hits) >= 1

        # Wildcard-only: wildcard groups appear, H custom group does not.
        assert len(wc_hits) >= 1
        assert [m for m in _type_matches(res_wc[0], "HalogenGroup") if m.get("id") == 997] == []

        # H + wildcard together: both expected outputs appear.
        assert len(both_h_hits) >= 1
        assert len(both_wc_hits) >= 1


class TestHComponentParseIntegration:
    def test_parse_smiles_h_custom_group_matches_ch2_connected_functional_group(self):
        # Custom group: alcohol SMARTS constrained to Alkyl component with H halogen path.
        custom_groups = [_build_h_custom_group(995)]

        res = parse_smiles("CCCCO", pfas_groups=custom_groups, halogens="H", bycomponent=True)
        assert len(res) == 1

        custom_hits = [
            m for m in _type_matches(res[0], "HalogenGroup")
            if m.get("id") == 995
        ]
        assert custom_hits, "Expected custom H-path group to match but got no match"
        assert custom_hits[0].get("components"), "Expected non-empty matched components for H-alkyl path"

    def test_parse_smiles_h_custom_group_no_match_without_ch2_component(self):
        custom_groups = [_build_h_custom_group(996)]

        # Methanol has no CH2 unit, so it should not satisfy the Alkyl=[CH2] component path.
        res = parse_smiles("CO", pfas_groups=custom_groups, halogens="H", bycomponent=True)
        assert len(res) == 1
        custom_hits = [
            m for m in _type_matches(res[0], "HalogenGroup")
            if m.get("id") == 996
        ]
        assert custom_hits == []
