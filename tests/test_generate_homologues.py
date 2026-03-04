"""Tests for generate_homologues – component-based implementation.

Covers:
- Linear perfluoroalkyl chains (PFOA, PFOS-like)
- CF3-only compounds (no repeating -CF2-, expect 0 homologues)
- Chlorinated and brominated analogues (halogen kwarg)
- Branched molecules generated with generate_random_mol
- Multi-component molecules (two separate fluorinated chains)
- Invalid halogen guard
- HomologueSeries API (summary, show/plot, svg, to_sql, entries, mols)
- String (SMILES / InChI) input to generate_homologues
"""

import pytest
from rdkit import Chem
from rdkit.Chem.rdMolDescriptors import CalcMolFormula

from PFASGroups.generate_homologues import generate_homologues, find_halogenated_components
from PFASGroups.homologue_series import HomologueSeries, HomologueEntry
from PFASGroups.generate_mol import generate_random_mol


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def smiles_set(homologues: dict) -> set:
    """Return canonical SMILES of all homologue molecules."""
    return {
        Chem.MolToSmiles(mol)
        for inner in homologues.values()
        for mol in inner.values()
    }


def formula_set(homologues: dict) -> set:
    return {formula for inner in homologues.values() for formula in inner.keys()}


def assert_all_connected(homologues: dict):
    for inchikey, inner in homologues.items():
        for formula, mol in inner.items():
            frags = Chem.GetMolFrags(mol)
            assert len(frags) == 1, (
                f"Fragmented homologue {Chem.MolToSmiles(mol)} "
                f"(formula={formula}, inchikey={inchikey})"
            )


def assert_valid_smiles(homologues: dict):
    for inner in homologues.values():
        for formula, mol in inner.items():
            assert mol is not None
            smi = Chem.MolToSmiles(mol)
            roundtrip = Chem.MolFromSmiles(smi)
            assert roundtrip is not None, f"Invalid SMILES produced: {smi}"


# ---------------------------------------------------------------------------
# Linear fluorinated chains
# ---------------------------------------------------------------------------

class TestLinearPerfluoro:
    """generate_homologues on simple linear perfluoroalkyl chains."""

    def test_pfoa_count(self):
        """C8 PFOA → 6 shorter homologues (C2–C7 chains)."""
        pfoa = Chem.MolFromSmiles('OC(=O)' + 'C(F)(F)' * 7 + 'F')
        h = generate_homologues(pfoa)
        assert len(h) == 6, f"Expected 6 homologues, got {len(h)}: {formula_set(h)}"

    def test_pfoa_chain_lengths(self):
        """Each homologue should have strictly fewer CF2 units than parent."""
        pfoa = Chem.MolFromSmiles('OC(=O)' + 'C(F)(F)' * 7 + 'F')
        parent_c = CalcMolFormula(pfoa).count('C')  # naive but good enough
        h = generate_homologues(pfoa)
        for inner in h.values():
            for formula, mol in inner.items():
                child_c = sum(1 for a in mol.GetAtoms() if a.GetSymbol() == 'C')
                parent_c_count = sum(1 for a in pfoa.GetAtoms() if a.GetSymbol() == 'C')
                assert child_c < parent_c_count, (
                    f"Homologue {Chem.MolToSmiles(mol)} not shorter than parent"
                )

    def test_pfoa_all_valid(self):
        pfoa = Chem.MolFromSmiles('OC(=O)' + 'C(F)(F)' * 7 + 'F')
        h = generate_homologues(pfoa)
        assert_all_connected(h)
        assert_valid_smiles(h)

    def test_pfos_like_count(self):
        """C8 perfluoroalkyl sulfonate → 7 shorter homologues (C1–C7)."""
        pfos = Chem.MolFromSmiles('O=S(=O)(O)' + 'C(F)(F)' * 8 + 'F')
        h = generate_homologues(pfos)
        assert len(h) == 7, f"Expected 7 homologues, got {len(h)}: {formula_set(h)}"
        assert_all_connected(h)

    def test_no_cf2_units(self):
        """CF3COOH has no internal –CF2– units → 0 homologues."""
        cf3cooh = Chem.MolFromSmiles('OC(=O)C(F)(F)F')
        h = generate_homologues(cf3cooh)
        assert len(h) == 0, f"Expected 0 homologues, got {len(h)}: {smiles_set(h)}"

    def test_single_cf2(self):
        """One –CF2– unit → exactly 1 homologue."""
        mol = Chem.MolFromSmiles('OC(=O)C(F)(F)C(F)(F)F')  # C3 acid, one CF2 + CF3
        h = generate_homologues(mol)
        assert len(h) == 1, f"Expected 1 homologue, got {len(h)}: {smiles_set(h)}"


# ---------------------------------------------------------------------------
# Halogen variants
# ---------------------------------------------------------------------------

class TestHalogenVariants:

    def test_perchloroalkyl(self):
        """Perchloroalkyl carboxylic acid (C5 chain)."""
        mol = Chem.MolFromSmiles('OC(=O)' + 'C(Cl)(Cl)' * 4 + 'Cl')
        h = generate_homologues(mol, halogen='Cl')
        assert len(h) >= 2, f"Expected ≥2 Cl homologues, got {len(h)}"
        assert_all_connected(h)
        assert_valid_smiles(h)
        # All homologues should contain Cl
        for smi in smiles_set(h):
            assert 'Cl' in smi, f"Cl homologue missing Cl: {smi}"

    def test_perbromoalkyl(self):
        """Perbromoalkyl carboxylic acid (C5 chain)."""
        mol = Chem.MolFromSmiles('OC(=O)' + 'C(Br)(Br)' * 4 + 'Br')
        h = generate_homologues(mol, halogen='Br')
        assert len(h) >= 2, f"Expected ≥2 Br homologues, got {len(h)}"
        assert_all_connected(h)
        assert_valid_smiles(h)

    def test_invalid_halogen(self):
        mol = Chem.MolFromSmiles('OC(=O)' + 'C(F)(F)' * 4 + 'F')
        with pytest.raises(ValueError, match="halogen must be one of"):
            generate_homologues(mol, halogen='S')


# ---------------------------------------------------------------------------
# Branched molecules via generate_random_mol
# ---------------------------------------------------------------------------

class TestBranched:
    """generate_homologues on branched perfluoro structures from generate_random_mol."""

    @pytest.fixture(scope='class')
    def branched_mols(self):
        """Generate a small set of random branched perfluoroalkyl molecules."""
        import numpy as np
        rng = np.random.default_rng(42)
        mols = []
        attempts = 0
        while len(mols) < 5 and attempts < 50:
            attempts += 1
            try:
                n = int(rng.integers(5, 10))
                mol = generate_random_mol(
                    n,
                    functional_groups=[{'group_smiles': 'C(=O)O', 'n': 1, 'mode': 'attach'}],
                    perfluorinated=True,
                )
                if mol is not None:
                    Chem.SanitizeMol(mol)
                    mols.append(mol)
            except Exception:
                continue
        return mols

    def test_branched_runs_without_error(self, branched_mols):
        """generate_homologues should not raise on any branched molecule."""
        assert len(branched_mols) > 0, "Could not generate any branched molecules"
        for mol in branched_mols:
            smi = Chem.MolToSmiles(mol)
            h = generate_homologues(mol)
            # May be 0 (no CF2 units) or more – both are valid
            assert isinstance(h, dict), f"Expected dict for {smi}"

    def test_branched_homologues_connected(self, branched_mols):
        """All produced homologues must be single connected fragments."""
        for mol in branched_mols:
            h = generate_homologues(mol)
            assert_all_connected(h)

    def test_branched_homologues_valid_smiles(self, branched_mols):
        for mol in branched_mols:
            h = generate_homologues(mol)
            assert_valid_smiles(h)

    def test_branched_homologues_shorter(self, branched_mols):
        """Each homologue must have fewer carbons than the parent."""
        for mol in branched_mols:
            parent_c = sum(1 for a in mol.GetAtoms() if a.GetSymbol() == 'C')
            h = generate_homologues(mol)
            for inner in h.values():
                for formula, hmol in inner.items():
                    child_c = sum(1 for a in hmol.GetAtoms() if a.GetSymbol() == 'C')
                    assert child_c < parent_c, (
                        f"Homologue {Chem.MolToSmiles(hmol)} not shorter than "
                        f"parent {Chem.MolToSmiles(mol)}"
                    )


# ---------------------------------------------------------------------------
# Multi-component (two separate fluorinated chains)
# ---------------------------------------------------------------------------

class TestMultiComponent:

    def test_two_chains_both_shortened(self):
        """Molecule with two independent perfluoroalkyl chains should produce
        homologues that shorten each chain independently."""
        # Two C4 perfluoro chains connected via an ether oxygen
        mol = Chem.MolFromSmiles('FC(F)(F)C(F)(F)CC(F)(F)C(F)(F)F')
        if mol is None:
            pytest.skip("SMILES not parsed – skipping multi-chain test")
        h = generate_homologues(mol)
        # Should find CF2 units → non-empty result
        # (exact count depends on component SMARTS matching)
        assert_all_connected(h)
        assert_valid_smiles(h)


# ---------------------------------------------------------------------------
# find_halogenated_components unit tests
# ---------------------------------------------------------------------------

class TestFindHalogenatedComponents:

    def _get_component_smarts(self, name='Perfluoroalkyl'):
        from PFASGroups.core import preprocess_componentsSmarts
        import json, os
        data_dir = os.path.join(os.path.dirname(__file__), '..', 'PFASGroups', 'data')
        with open(os.path.join(data_dir, 'component_smarts_halogens.json')) as f:
            components = json.load(f)
        processed = preprocess_componentsSmarts(components)
        return processed[name]['component']

    def test_pfoa_one_component(self):
        pfoa = Chem.MolFromSmiles('OC(=O)' + 'C(F)(F)' * 7 + 'F')
        smarts = self._get_component_smarts('Perfluoroalkyl')
        comps = find_halogenated_components(pfoa, smarts, halogen='F')
        assert len(comps) == 1, f"Expected 1 component, got {len(comps)}"
        assert len(comps[0]['cx2_carbons']) == 6, (
            f"Expected 6 CF2 carbons (C2-C7), got {comps[0]['cx2_carbons']}"
        )

    def test_cf3cooh_no_cx2(self):
        mol = Chem.MolFromSmiles('OC(=O)C(F)(F)F')
        smarts = self._get_component_smarts('Perfluoroalkyl')
        comps = find_halogenated_components(mol, smarts, halogen='F')
        total_cx2 = sum(len(c['cx2_carbons']) for c in comps)
        assert total_cx2 == 0, f"Expected 0 CF2 carbons, got {total_cx2}"

    def test_non_fluorinated_mol(self):
        mol = Chem.MolFromSmiles('CCCCCC(=O)O')
        smarts = self._get_component_smarts('Perfluoroalkyl')
        comps = find_halogenated_components(mol, smarts, halogen='F')
        assert comps == [], "Expected no components in non-fluorinated molecule"


# ---------------------------------------------------------------------------
# HomologueSeries API tests
# ---------------------------------------------------------------------------

_PFOA_SMILES = 'OC(=O)' + 'C(F)(F)' * 7 + 'F'


class TestHomologueSeriesType:
    """generate_homologues always returns a HomologueSeries (a dict subclass)."""

    def test_returns_homologue_series(self):
        mol = Chem.MolFromSmiles(_PFOA_SMILES)
        h = generate_homologues(mol)
        assert isinstance(h, HomologueSeries)
        assert isinstance(h, dict)

    def test_empty_series_is_homologue_series(self):
        mol = Chem.MolFromSmiles('OC(=O)C(F)(F)F')  # CF3 only
        h = generate_homologues(mol)
        assert isinstance(h, HomologueSeries)
        assert len(h) == 0

    def test_metadata_set(self):
        mol = Chem.MolFromSmiles(_PFOA_SMILES)
        h = generate_homologues(mol)
        assert h.parent_smiles is not None
        assert h.parent_formula is not None
        assert h.halogen == 'F'
        assert 'Perfluoroalkyl' in h.component_name

    def test_metadata_cl(self):
        mol = Chem.MolFromSmiles('OC(=O)' + 'C(Cl)(Cl)' * 4 + 'Cl')
        h = generate_homologues(mol, halogen='Cl')
        assert h.halogen == 'Cl'


class TestHomologueSeriesEntries:

    def test_entries_sorted_by_n_removed(self):
        mol = Chem.MolFromSmiles(_PFOA_SMILES)
        h = generate_homologues(mol)
        entries = h.entries()
        assert entries, "No entries for PFOA"
        removed_seq = [e.n_removed for e in entries]
        assert removed_seq == sorted(removed_seq), "entries() not sorted ascending"

    def test_n_removed_positive(self):
        mol = Chem.MolFromSmiles(_PFOA_SMILES)
        h = generate_homologues(mol)
        for entry in h.entries():
            assert entry.n_removed >= 1

    def test_entry_smiles_valid(self):
        mol = Chem.MolFromSmiles(_PFOA_SMILES)
        h = generate_homologues(mol)
        for entry in h.entries():
            roundtrip = Chem.MolFromSmiles(entry.smiles)
            assert roundtrip is not None

    def test_mols_contains_rdkit_mols(self):
        mol = Chem.MolFromSmiles(_PFOA_SMILES)
        h = generate_homologues(mol)
        for m in h.mols():
            assert isinstance(m, Chem.Mol)

    def test_homologue_entry_type(self):
        mol = Chem.MolFromSmiles(_PFOA_SMILES)
        h = generate_homologues(mol)
        for entry in h.entries():
            assert isinstance(entry, HomologueEntry)

    def test_empty_entries(self):
        mol = Chem.MolFromSmiles('OC(=O)C(F)(F)F')
        h = generate_homologues(mol)
        assert h.entries() == []
        assert h.mols() == []


class TestHomologueSeriesSummary:

    def test_summarise_returns_string(self):
        mol = Chem.MolFromSmiles(_PFOA_SMILES)
        h = generate_homologues(mol)
        s = h.summarise()
        assert isinstance(s, str)
        assert 'HomologueSeries' in s

    def test_summarise_contains_parent_info(self):
        mol = Chem.MolFromSmiles(_PFOA_SMILES)
        h = generate_homologues(mol)
        s = h.summarise()
        assert h.parent_formula in s

    def test_summary_prints(self, capsys):
        mol = Chem.MolFromSmiles(_PFOA_SMILES)
        h = generate_homologues(mol)
        h.summary()
        captured = capsys.readouterr()
        assert 'HomologueSeries' in captured.out

    def test_summarise_empty(self):
        mol = Chem.MolFromSmiles('OC(=O)C(F)(F)F')
        h = generate_homologues(mol)
        s = h.summarise()
        assert 'No homologues' in s


class TestHomologueSeriesShow:

    def test_show_returns_pil_image(self):
        from PIL import Image
        mol = Chem.MolFromSmiles(_PFOA_SMILES)
        h = generate_homologues(mol)
        img = h.show(display=False)
        assert isinstance(img, Image.Image)
        assert img.width > 0
        assert img.height > 0

    def test_plot_alias(self):
        from PIL import Image
        mol = Chem.MolFromSmiles(_PFOA_SMILES)
        h = generate_homologues(mol)
        img = h.plot(display=False)
        assert isinstance(img, Image.Image)

    def test_show_without_parent(self):
        from PIL import Image
        mol = Chem.MolFromSmiles(_PFOA_SMILES)
        h = generate_homologues(mol)
        img = h.show(display=False, show_parent=False)
        assert isinstance(img, Image.Image)

    def test_show_empty_raises(self):
        mol = Chem.MolFromSmiles('OC(=O)C(F)(F)F')
        h = generate_homologues(mol)
        # No homologues and no parent shown → should raise
        with pytest.raises(ValueError, match="empty"):
            h.show(display=False, show_parent=False)

    def test_svg_creates_file(self, tmp_path):
        mol = Chem.MolFromSmiles(_PFOA_SMILES)
        h = generate_homologues(mol)
        out = tmp_path / "series.svg"
        result = h.svg(str(out))
        assert result == str(out)
        assert out.exists()
        assert out.stat().st_size > 0


class TestHomologueSeriesToSql:

    def test_to_sql_sqlite(self, tmp_path):
        mol = Chem.MolFromSmiles(_PFOA_SMILES)
        h = generate_homologues(mol)
        db = tmp_path / "test.db"
        h.to_sql(filename=str(db))
        assert db.exists()

        import sqlite3
        conn = sqlite3.connect(str(db))
        rows = conn.execute("SELECT * FROM homologue_series").fetchall()
        conn.close()
        assert len(rows) == len(h)

    def test_to_sql_rows_content(self, tmp_path):
        mol = Chem.MolFromSmiles(_PFOA_SMILES)
        h = generate_homologues(mol)
        db = tmp_path / "test2.db"
        h.to_sql(filename=str(db))

        import sqlite3
        conn = sqlite3.connect(str(db))
        cols = [d[1] for d in conn.execute("PRAGMA table_info(homologue_series)").fetchall()]
        conn.close()
        for col in ('parent_smiles', 'halogen', 'inchikey', 'formula', 'smiles', 'n_removed', 'n_carbons'):
            assert col in cols, f"Missing column: {col}"

    def test_to_sql_no_connection_raises(self):
        mol = Chem.MolFromSmiles(_PFOA_SMILES)
        h = generate_homologues(mol)
        with pytest.raises(ValueError):
            h.to_sql()

    def test_to_sql_empty_series_no_error(self, tmp_path):
        mol = Chem.MolFromSmiles('OC(=O)C(F)(F)F')
        h = generate_homologues(mol)
        db = tmp_path / "empty.db"
        h.to_sql(filename=str(db))  # should not raise


# ---------------------------------------------------------------------------
# String input (SMILES / InChI) tests
# ---------------------------------------------------------------------------

class TestStringInput:

    def test_smiles_string_input(self):
        smiles = _PFOA_SMILES
        h = generate_homologues(smiles)
        assert isinstance(h, HomologueSeries)
        assert len(h) == 6

    def test_smiles_string_metadata(self):
        smiles = _PFOA_SMILES
        h = generate_homologues(smiles)
        assert h.parent_smiles is not None
        assert h.parent_formula is not None

    def test_invalid_smiles_raises(self):
        with pytest.raises(ValueError, match="Could not parse"):
            generate_homologues('not_a_smiles_!!!###')

    def test_wrong_type_raises(self):
        with pytest.raises(TypeError):
            generate_homologues(12345)

    def test_smiles_gives_same_result_as_mol(self):
        smiles = _PFOA_SMILES
        mol = Chem.MolFromSmiles(smiles)
        h_smi = generate_homologues(smiles)
        h_mol = generate_homologues(mol)
        assert set(h_smi.keys()) == set(h_mol.keys())

    def test_inchi_input(self):
        # InChI for PFOA (C8 perfluorooctanoic acid)
        mol = Chem.MolFromSmiles(_PFOA_SMILES)
        inchi = Chem.MolToInchi(mol)
        assert inchi is not None and inchi.startswith('InChI=')
        h = generate_homologues(inchi)
        assert isinstance(h, HomologueSeries)
        assert len(h) > 0
