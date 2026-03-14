"""Tests for PFASEmbedding.to_array() and PFASEmbeddingSet.to_array()."""

import pytest
import numpy as np
from HalogenGroups import parse_smiles
from PFASGroups.results_model import PFASEmbedding, PFASEmbeddingSet
from PFASGroups import generate_embedding

TEST_SMILES = [
    "C(C(C(C(C(C(C(C(F)(F)F)(F)F)(F)F)(F)F)(F)F)(F)F)(F)F)(C(=O)O)F",  # PFOA
    "C(C(C(C(C(F)(F)F)(F)F)(F)F)(F)F)(C(=O)O)F",  # PFBA
    "C(C(C(C(C(C(F)(F)F)(F)F)(F)F)(F)F)(F)F)(S(=O)(=O)O)F",  # PFHxS
    "FC(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)S(=O)(=O)O",  # PFOS
    "FC(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(=O)O",  # PFHxA
]


@pytest.fixture(scope="module")
def results():
    """Parse TEST_SMILES once for the whole module."""
    return parse_smiles(TEST_SMILES)


# ── TestToArray ────────────────────────────────────────────────────────────────

class TestToArray:
    """Tests for PFASEmbeddingSet.to_array()."""

    def test_default_binary_shape(self, results):
        """to_array() returns shape (n_mols, n_groups) with n_groups > 0."""
        arr = results.to_array()
        assert isinstance(arr, np.ndarray)
        assert arr.ndim == 2
        assert arr.shape[0] == len(TEST_SMILES)
        assert arr.shape[1] > 0

    def test_binary_values(self, results):
        """Default binary mode contains only 0s and 1s."""
        arr = results.to_array()
        assert np.all((arr == 0) | (arr == 1)), (
            'Binary array should contain only 0 and 1'
        )

    def test_count_mode(self, results):
        """count mode values are non-negative; positions that are 0 in binary stay 0 in count."""
        arr_bin = results.to_array(component_metrics=['binary'])
        arr_cnt = results.to_array(component_metrics=['count'])
        assert np.all(arr_cnt >= 0), 'Count values must be non-negative'
        # Wherever binary is 0 the count must also be 0
        assert np.all(arr_cnt[arr_bin == 0] == 0), (
            'Count should be 0 wherever binary is 0'
        )

    def test_multiple_metrics(self, results):
        """['binary', 'effective_graph_resistance'] gives shape (n_mols, 2*n_groups)."""
        arr_bin = results.to_array(component_metrics=['binary'])
        n_groups = arr_bin.shape[1]
        arr_multi = results.to_array(
            component_metrics=['binary', 'effective_graph_resistance']
        )
        assert arr_multi.shape == (len(TEST_SMILES), 2 * n_groups)

    def test_aggregation_mean_median(self, results):
        """mean and median aggregations give the same shape; may differ for multi-component mols."""
        arr_mean = results.to_array(
            component_metrics=['effective_graph_resistance'], aggregation='mean'
        )
        arr_median = results.to_array(
            component_metrics=['effective_graph_resistance'], aggregation='median'
        )
        assert arr_mean.shape == arr_median.shape

    def test_group_selection_oecd(self, results):
        """group_selection='oecd' gives fewer columns than 'all'."""
        arr_all  = results.to_array(group_selection='all')
        arr_oecd = results.to_array(group_selection='oecd')
        assert arr_oecd.shape[1] < arr_all.shape[1], (
            f"oecd ({arr_oecd.shape[1]}) should have fewer columns than all ({arr_all.shape[1]})"
        )

    def test_molecule_metrics(self, results):
        """molecule_metrics=['n_components'] appends one extra column."""
        arr_base = results.to_array()
        n_base = arr_base.shape[1]
        arr_mm = results.to_array(molecule_metrics=['n_components'])
        assert arr_mm.shape[1] == n_base + 1

    def test_preset_best(self, results):
        """preset='best' gives a non-empty array."""
        arr = results.to_array(preset='best')
        assert isinstance(arr, np.ndarray)
        assert arr.size > 0

    def test_egr_bde_metric(self, results):
        """component_metrics=['binary', 'effective_graph_resistance_BDE'] works and is non-negative."""
        arr = results.to_array(
            component_metrics=['binary', 'effective_graph_resistance_BDE']
        )
        assert arr.ndim == 2
        assert arr.shape[0] == len(TEST_SMILES)
        # EGR values (second block) should be non-negative where defined
        n_groups = arr.shape[1] // 2
        egr_block = arr[:, n_groups:]
        assert np.all(egr_block >= 0), 'effective_graph_resistance_BDE values must be non-negative'


# ── TestColumnNames ────────────────────────────────────────────────────────────

class TestColumnNames:
    """Tests for column_names()."""

    def test_format(self, results):
        """Column names follow pattern '{group} [metric]', e.g. ' [binary]' substring."""
        cols = results.column_names()
        assert len(cols) > 0
        assert all('[binary]' in c for c in cols), (
            'Default (binary) column names should all contain "[binary]"'
        )

    def test_count_matches_array(self, results):
        """len(column_names()) matches to_array().shape[1]."""
        arr = results.to_array()
        cols = results.column_names()
        assert len(cols) == arr.shape[1]

    def test_molecule_metrics_prefix(self, results):
        """Molecule metric columns carry a 'mol:' prefix."""
        cols = results.column_names(molecule_metrics=['n_components', 'mean_branching'])
        mol_cols = [c for c in cols if c.startswith('mol:')]
        assert 'mol:n_components' in mol_cols
        assert 'mol:mean_branching' in mol_cols

    def test_preset_consistency(self, results):
        """Same preset gives identical column names from both to_array and column_names."""
        arr = results.to_array(preset='best')
        cols = results.column_names(preset='best')
        assert len(cols) == arr.shape[1]


# ── TestSingleMolecule ────────────────────────────────────────────────────────

class TestSingleMolecule:
    """Tests for PFASEmbedding.to_array() on a single molecule."""

    @pytest.fixture(scope="class")
    def single_result(self):
        res = parse_smiles(TEST_SMILES)
        return res[0]  # PFASEmbedding dict subclass

    def test_single_returns_1d(self, single_result):
        """results[0].to_array() returns a 1-D array."""
        arr = single_result.to_array()
        assert isinstance(arr, np.ndarray)
        assert arr.ndim == 1
        assert arr.shape[0] > 0

    def test_single_column_names(self, single_result):
        """results[0].column_names() length matches to_array() length."""
        arr = single_result.to_array()
        cols = single_result.column_names()
        assert len(cols) == arr.shape[0]

    def test_single_consistent_with_batch(self, results, single_result):
        """Single-molecule to_array() matches row 0 of batch to_array()."""
        arr_single = single_result.to_array()
        arr_batch = results.to_array()
        np.testing.assert_array_equal(arr_single, arr_batch[0])


# ── TestGenerateEmbedding ─────────────────────────────────────────────────────

class TestGenerateEmbedding:
    """Tests for the generate_embedding() convenience function."""

    def test_returns_tuple(self):
        """generate_embedding returns a (array, column_names) tuple."""
        result = generate_embedding(TEST_SMILES[0])
        assert isinstance(result, tuple)
        assert len(result) == 2

    def test_single_smiles_1d(self):
        """Single SMILES gives a 1-D embedding array."""
        arr, cols = generate_embedding(TEST_SMILES[0])
        assert arr.ndim == 1
        assert arr.shape[0] > 0

    def test_list_smiles_2d(self):
        """List of SMILES gives a 2-D embedding array."""
        arr, cols = generate_embedding(TEST_SMILES)
        assert arr.ndim == 2
        assert arr.shape[0] == len(TEST_SMILES)
        assert arr.shape[1] > 0

    def test_column_names_length(self):
        """len(cols) == arr.shape[-1] for both 1-D and 2-D outputs."""
        arr1, cols1 = generate_embedding(TEST_SMILES[0])
        assert len(cols1) == arr1.shape[-1]

        arr2, cols2 = generate_embedding(TEST_SMILES)
        assert len(cols2) == arr2.shape[-1]

    def test_preset(self):
        """preset='best' works and produces a non-empty result."""
        arr, cols = generate_embedding(TEST_SMILES, preset='best')
        assert arr.size > 0
        assert len(cols) > 0


# ── TestBackwardCompat ────────────────────────────────────────────────────────

class TestBackwardCompat:
    """Tests for deprecated to_fingerprint()."""

    def test_to_fingerprint_deprecated(self, results):
        """Calling to_fingerprint() raises DeprecationWarning."""
        with pytest.warns(DeprecationWarning):
            results.to_fingerprint()

    def test_to_fingerprint_returns_array(self, results):
        """to_fingerprint() returns a numpy ndarray."""
        with pytest.warns(DeprecationWarning):
            result = results.to_fingerprint()
        assert isinstance(result, np.ndarray)

    def test_to_fingerprint_shape(self, results):
        """to_fingerprint() shape matches to_array() shape for same parameters."""
        with pytest.warns(DeprecationWarning):
            fp = results.to_fingerprint()
        arr = results.to_array()
        assert fp.shape == arr.shape


# ── TestEdgeCases ─────────────────────────────────────────────────────────────

class TestEdgeCases:

    def test_empty_results(self):
        """PFASEmbeddingSet([]).to_array() returns shape (0, 0)."""
        empty = PFASEmbeddingSet([])
        arr = empty.to_array()
        assert arr.shape == (0, 0)

    def test_invalid_preset(self, results):
        """to_array(preset='nonexistent') raises ValueError."""
        with pytest.raises(ValueError):
            results.to_array(preset='nonexistent')
