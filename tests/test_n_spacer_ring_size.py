"""Tests for n_spacer (telomer CH₂ linker length) and ring_size metrics."""

import pytest
from rdkit import Chem
from PFASGroups.parser import parse_groups_in_mol


def _components_by_group(mol, group_name_substr):
    """Return all components whose group name contains *group_name_substr*."""
    matches, *_ = parse_groups_in_mol(mol)
    found = []
    for grp, _idx, _smarts, comps in matches:
        if group_name_substr.lower() in grp.name.lower():
            found.extend(comps)
    return found


# ─────────────────────────────────────────────────────────────────────────────
# n_spacer  – telomer CH₂ linker chain length
# ─────────────────────────────────────────────────────────────────────────────

class TestNSpacer:
    """n_spacer increases by 1 for each additional CH₂ in the spacer arm."""

    def test_telomer_iodide_n1(self):
        # perfluorobutyl iodide: CF₃(CF₂)₃–CH₂–I  →  1 CH₂ spacer
        mol = Chem.MolFromSmiles("ICC(F)(F)C(F)(F)C(F)(F)C(F)(F)F")
        assert mol is not None
        comps = _components_by_group(mol, "telomer iodide")
        assert comps, "Expected a 'telomer iodide' match"
        assert comps[0]["n_spacer"] == 1

    def test_telomer_iodide_n2(self):
        # CF₃(CF₂)₃–CH₂–CH₂–I  →  2 CH₂ spacers
        mol = Chem.MolFromSmiles("ICCC(F)(F)C(F)(F)C(F)(F)C(F)(F)F")
        assert mol is not None
        comps = _components_by_group(mol, "telomer iodide")
        assert comps, "Expected a 'telomer iodide' match"
        assert comps[0]["n_spacer"] == 2

    def test_telomer_alcohol_n2(self):
        # 4:2 FTOH  CF₃(CF₂)₃–CH₂–CH₂–OH  →  2 CH₂ spacers
        mol = Chem.MolFromSmiles("OCCC(F)(F)C(F)(F)C(F)(F)C(F)(F)F")
        assert mol is not None
        comps = _components_by_group(mol, "telomer alcohol")
        assert comps, "Expected a 'telomer alcohol' match"
        spacers = [c["n_spacer"] for c in comps]
        assert 2 in spacers, f"Expected n_spacer=2 among {spacers}"

    def test_telomer_alcohol_n3(self):
        # CF₃(CF₂)₃–CH₂–CH₂–CH₂–OH  →  3 CH₂ spacers
        mol = Chem.MolFromSmiles("OCCCC(F)(F)C(F)(F)C(F)(F)C(F)(F)F")
        assert mol is not None
        comps = _components_by_group(mol, "telomer alcohol")
        assert comps, "Expected a 'telomer alcohol' match"
        spacers = [c["n_spacer"] for c in comps]
        assert 3 in spacers, f"Expected n_spacer=3 among {spacers}"

    def test_telomer_carboxylic_acid_n1(self):
        # CF₃(CF₂)₃–CH₂–COOH  →  1 CH₂ spacer
        mol = Chem.MolFromSmiles("OC(=O)CC(F)(F)C(F)(F)C(F)(F)C(F)(F)F")
        assert mol is not None
        comps = _components_by_group(mol, "telomer carboxylic acids")
        assert comps, "Expected a 'telomer carboxylic acids' match"
        assert comps[0]["n_spacer"] == 1

    def test_telomer_carboxylic_acid_n2(self):
        # CF₃(CF₂)₃–CH₂–CH₂–COOH  →  2 CH₂ spacers
        mol = Chem.MolFromSmiles("OC(=O)CCC(F)(F)C(F)(F)C(F)(F)C(F)(F)F")
        assert mol is not None
        comps = _components_by_group(mol, "telomer carboxylic acids")
        assert comps, "Expected a 'telomer carboxylic acids' match"
        assert comps[0]["n_spacer"] == 2

    def test_n_spacer_zero_for_non_telomer(self):
        # PFOA (perfluorooctanoic acid) – no CH₂ spacer; non-telomer groups → 0
        mol = Chem.MolFromSmiles("OC(=O)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)F")
        assert mol is not None
        matches, *_ = parse_groups_in_mol(mol)
        for grp, _idx, _smarts, comps in matches:
            if "telomer" not in grp.name.lower():
                for comp in comps:
                    assert comp["n_spacer"] == 0, (
                        f"Non-telomer group '{grp.name}' should have n_spacer=0, "
                        f"got {comp['n_spacer']}"
                    )


# ─────────────────────────────────────────────────────────────────────────────
# ring_size – smallest ring overlapping with matched component
# ─────────────────────────────────────────────────────────────────────────────

class TestRingSize:
    """ring_size reflects the actual ring size of the halogenated ring fragment."""

    def test_cyclic_6_perfluorocyclohexane(self):
        # Perfluorocyclohexane – 6-membered ring
        mol = Chem.MolFromSmiles("FC1(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C1(F)F")
        assert mol is not None
        comps = _components_by_group(mol, "perhalogenated cyclic")
        assert comps, "Expected a 'perhalogenated cyclic compounds' match"
        assert comps[0]["ring_size"] == 6

    def test_cyclic_5_perfluorocyclopentane(self):
        # Perfluorocyclopentane – 5-membered ring
        mol = Chem.MolFromSmiles("FC1(F)C(F)(F)C(F)(F)C(F)(F)C1(F)F")
        assert mol is not None
        comps = _components_by_group(mol, "perhalogenated cyclic")
        assert comps, "Expected a 'perhalogenated cyclic compounds' match"
        assert comps[0]["ring_size"] == 5

    def test_aryl_6_perfluorobenzene(self):
        # Perfluorobenzene – 6-membered aromatic ring
        mol = Chem.MolFromSmiles("Fc1c(F)c(F)c(F)c(F)c1F")
        assert mol is not None
        comps = _components_by_group(mol, "perhalogenated aryl")
        assert comps, "Expected a 'perhalogenated aryl compounds' match"
        assert comps[0]["ring_size"] == 6

    def test_aryl_6_polyfluorobenzene(self):
        # 1,2-difluorobenzene – matched as polyhalogenated aryl, 6-membered ring
        mol = Chem.MolFromSmiles("Fc1cccc(F)c1")
        assert mol is not None
        comps = _components_by_group(mol, "polyhalogenated aryl")
        assert comps, "Expected a 'polyhalogenated aryl compounds' match"
        assert comps[0]["ring_size"] == 6

    def test_heterocyclic_azole_5(self):
        # 4-(heptadecafluorooctyl)-1H-pyrazole – pyrazole is a 5-membered ring
        mol = Chem.MolFromSmiles(
            "FC(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C1=CNN=C1"
        )
        assert mol is not None
        comps = _components_by_group(mol, "heterocyclic azole")
        assert comps, "Expected a 'heterocyclic azole' match"
        ring_sizes = [c["ring_size"] for c in comps if c["ring_size"] > 0]
        assert ring_sizes, "At least one 'heterocyclic azole' component should have ring_size > 0"
        assert all(rs == 5 for rs in ring_sizes), (
            f"All heterocyclic azole ring sizes should be 5, got {ring_sizes}"
        )

    def test_ring_size_zero_for_linear(self):
        # Linear PFAS (PFOS) – no ring; all components must have ring_size=0
        mol = Chem.MolFromSmiles(
            "C(C(C(C(C(C(C(C(F)(F)S(=O)(=O)O)(F)F)(F)F)(F)F)(F)F)(F)F)(F)F)(F)F"
        )
        assert mol is not None
        matches, *_ = parse_groups_in_mol(mol)
        for grp, _idx, _smarts, comps in matches:
            for comp in comps:
                assert comp["ring_size"] == 0, (
                    f"Linear molecule group '{grp.name}' should have ring_size=0, "
                    f"got {comp['ring_size']}"
                )
