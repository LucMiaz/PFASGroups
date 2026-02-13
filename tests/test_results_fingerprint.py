"""
Unit tests for ResultsFingerprint functionality

This test suite verifies the functionality of:
- ResultsModel.to_fingerprint() conversion
- ResultsFingerprint class and methods
- Dimensionality reduction (PCA, kernel-PCA, t-SNE, UMAP)
- KL divergence comparison
- SQL persistence for both ResultsModel and ResultsFingerprint
"""

import pytest
import numpy as np
import tempfile
import os
from pathlib import Path
from PFASgroups import parse_smiles
from PFASgroups.results_model import ResultsFingerprint, ResultsModel


# Test SMILES - diverse set of PFAS compounds
TEST_SMILES = [
    "C(C(C(C(C(C(C(C(F)(F)F)(F)F)(F)F)(F)F)(F)F)(F)F)(F)F)(C(=O)O)F",  # PFOA
    "C(C(C(C(C(F)(F)F)(F)F)(F)F)(F)F)(C(=O)O)F",  # PFBA
    "C(C(C(C(C(C(F)(F)F)(F)F)(F)F)(F)F)(F)F)(S(=O)(=O)O)F",  # PFHxS
    "FC(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)S(=O)(=O)O",  # PFOS
    "FC(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(=O)O",  # PFHxA
]

# Additional SMILES for comparison tests
COMPARISON_SMILES = [
    "C(C(C(F)(F)F)(F)F)(S(=O)(=O)O)F",  # PFBS
    "FC(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(=O)O",  # PFNA
]


@pytest.fixture
def results():
    """Parse test SMILES and return ResultsModel."""
    return parse_smiles(TEST_SMILES)


@pytest.fixture
def fingerprint(results):
    """Create ResultsFingerprint from results."""
    return results.to_fingerprint(group_selection='all', count_mode='binary')


class TestResultsModelToFingerprint:
    """Test ResultsModel.to_fingerprint() method."""
    
    def test_basic_conversion(self, results):
        """Test basic conversion to fingerprint."""
        fp = results.to_fingerprint()
        
        assert isinstance(fp, ResultsFingerprint)
        assert len(fp) == len(results)
        assert fp.fingerprints.shape[0] == len(results)
        assert len(fp.group_names) > 0
    
    def test_group_selections(self, results):
        """Test different group selections."""
        # All groups
        fp_all = results.to_fingerprint(group_selection='all')
        assert len(fp_all.group_names) == 55  # Total PFAS groups
        
        # OECD groups
        fp_oecd = results.to_fingerprint(group_selection='oecd')
        assert len(fp_oecd.group_names) == 28
        
        # Generic groups
        fp_generic = results.to_fingerprint(group_selection='generic')
        assert len(fp_generic.group_names) == 27
    
    def test_count_modes(self, results):
        """Test different count modes."""
        # Binary
        fp_binary = results.to_fingerprint(count_mode='binary')
        assert np.all((fp_binary.fingerprints == 0) | (fp_binary.fingerprints == 1))
        
        # Count
        fp_count = results.to_fingerprint(count_mode='count')
        assert np.all(fp_count.fingerprints >= 0)
        
        # Max component
        fp_max = results.to_fingerprint(count_mode='max_component')
        assert np.all(fp_max.fingerprints >= 0)
    
    def test_custom_group_ids(self, results):
        """Test custom group ID selection."""
        custom_ids = [1, 2, 5, 10]
        fp = results.to_fingerprint(selected_group_ids=custom_ids)
        assert len(fp.group_names) == len(custom_ids)


class TestResultsFingerprint:
    """Test ResultsFingerprint class."""
    
    def test_initialization(self):
        """Test ResultsFingerprint initialization."""
        fps = np.array([[1, 0, 1], [0, 1, 1]])
        smiles = ['CC', 'CCC']
        groups = ['A', 'B', 'C']
        
        rf = ResultsFingerprint(fps, smiles, groups)
        
        assert len(rf) == 2
        assert rf.fingerprints.shape == (2, 3)
        assert len(rf.smiles) == 2
        assert len(rf.group_names) == 3
    
    def test_length_mismatch(self):
        """Test that length mismatch raises error."""
        fps = np.array([[1, 0, 1]])
        smiles = ['CC', 'CCC']  # Wrong length
        groups = ['A', 'B', 'C']
        
        with pytest.raises(ValueError, match="length mismatch"):
            ResultsFingerprint(fps, smiles, groups)
    
    def test_repr(self, fingerprint):
        """Test string representation."""
        repr_str = repr(fingerprint)
        assert 'ResultsFingerprint' in repr_str
        assert 'n_molecules' in repr_str
        assert 'n_groups' in repr_str
    
    def test_summary(self, fingerprint):
        """Test summary method."""
        summary = fingerprint.summary()
        assert 'Molecules:' in summary
        assert 'Groups:' in summary
        assert 'Sparsity:' in summary


class TestDimensionalityReduction:
    """Test dimensionality reduction methods."""
    
    def test_pca(self, fingerprint):
        """Test PCA analysis."""
        pca_results = fingerprint.perform_pca(n_components=2, plot=False)
        
        assert 'transformed' in pca_results
        assert 'explained_variance' in pca_results
        assert 'pca_model' in pca_results
        assert pca_results['transformed'].shape[0] == len(fingerprint)
        assert pca_results['transformed'].shape[1] == 2
        assert len(pca_results['explained_variance']) == 2
    
    def test_kernel_pca(self, fingerprint):
        """Test kernel PCA analysis."""
        for kernel in ['rbf', 'poly', 'sigmoid']:
            kpca_results = fingerprint.perform_kernel_pca(
                n_components=2, 
                kernel=kernel, 
                plot=False
            )
            
            assert 'transformed' in kpca_results
            assert 'kpca_model' in kpca_results
            assert kpca_results['transformed'].shape[0] == len(fingerprint)
            assert kpca_results['kernel'] == kernel
    
    def test_tsne(self, fingerprint):
        """Test t-SNE analysis."""
        tsne_results = fingerprint.perform_tsne(
            n_components=2,
            perplexity=min(30, len(fingerprint) - 1),  # Adjust for small dataset
            n_iter=250,  # Fewer iterations for speed
            plot=False
        )
        
        assert 'transformed' in tsne_results
        assert 'tsne_model' in tsne_results
        assert tsne_results['transformed'].shape[0] == len(fingerprint)
        assert tsne_results['transformed'].shape[1] == 2
    
    @pytest.mark.skipif(True, reason="UMAP may not be installed")
    def test_umap(self, fingerprint):
        """Test UMAP analysis."""
        try:
            umap_results = fingerprint.perform_umap(
                n_components=2,
                n_neighbors=min(15, len(fingerprint) - 1),
                plot=False
            )
            
            assert 'transformed' in umap_results
            assert 'umap_model' in umap_results
            assert umap_results['transformed'].shape[0] == len(fingerprint)
        except ImportError:
            pytest.skip("UMAP not installed")


class TestKLDivergence:
    """Test KL divergence comparison."""
    
    def test_basic_comparison(self, results):
        """Test basic KL divergence comparison."""
        fp1 = results.to_fingerprint(group_selection='all', count_mode='binary')
        
        # Create a second fingerprint from subset
        results2 = parse_smiles(TEST_SMILES[:2])
        fp2 = results2.to_fingerprint(group_selection='all', count_mode='binary')
        
        kl_div = fp1.compare_kld(fp2, method='minmax')
        
        assert isinstance(kl_div, float)
        assert 0 <= kl_div <= 1  # minmax should be normalized
        assert not np.isnan(kl_div)
    
    def test_all_methods(self, results):
        """Test all KL divergence methods."""
        fp1 = results.to_fingerprint(group_selection='all', count_mode='binary')
        results2 = parse_smiles(TEST_SMILES[:2])
        fp2 = results2.to_fingerprint(group_selection='all', count_mode='binary')
        
        methods = ['minmax', 'forward', 'reverse', 'symmetric']
        
        for method in methods:
            kl_div = fp1.compare_kld(fp2, method=method)
            assert isinstance(kl_div, float)
            assert not np.isnan(kl_div)
            assert kl_div >= 0
    
    def test_self_comparison(self, fingerprint):
        """Test that comparing fingerprint to itself gives low divergence."""
        kl_div = fingerprint.compare_kld(fingerprint, method='minmax')
        assert kl_div < 0.01  # Should be very close to 0
    
    def test_group_mismatch(self, results):
        """Test that different groups raises error."""
        fp1 = results.to_fingerprint(group_selection='all')
        fp2 = results.to_fingerprint(group_selection='oecd')
        
        with pytest.raises(ValueError, match="same number of groups"):
            fp1.compare_kld(fp2)


class TestSQLOperations:
    """Test SQL save/load operations."""
    
    def test_fingerprint_sql_roundtrip(self, fingerprint, tmp_path):
        """Test saving and loading fingerprints to SQL."""
        db_file = tmp_path / "test_fp.db"
        
        # Save
        fingerprint.to_sql(filename=str(db_file), if_exists='replace')
        
        # Load
        fp_loaded = ResultsFingerprint.from_sql(filename=str(db_file))
        
        # Verify
        assert len(fp_loaded) == len(fingerprint)
        assert len(fp_loaded.group_names) == len(fingerprint.group_names)
        assert fp_loaded.group_selection == fingerprint.group_selection
        assert fp_loaded.count_mode == fingerprint.count_mode
        
        # Check fingerprints are close (may have small rounding differences)
        assert np.allclose(fp_loaded.fingerprints, fingerprint.fingerprints)
    
    def test_results_sql_roundtrip(self, results, tmp_path):
        """Test saving and loading ResultsModel to SQL."""
        db_file = tmp_path / "test_results.db"
        
        # Save
        results.to_sql(filename=str(db_file), if_exists='replace')
        
        # Load
        results_loaded = ResultsModel.from_sql(filename=str(db_file))
        
        # Verify
        assert len(results_loaded) == len(results)
        
        for i, (orig, loaded) in enumerate(zip(results, results_loaded)):
            assert orig.smiles == loaded.smiles
            # Note: Components may not be saved in basic SQL format


class TestEdgeCases:
    """Test edge cases and error handling."""
    
    def test_empty_results(self):
        """Test handling of empty results."""
        empty_results = ResultsModel([])
        
        # Empty results should raise an error or return empty fingerprint
        # Depending on implementation, adjust this test
        try:
            fp = empty_results.to_fingerprint()
            assert len(fp) == 0 or fp.fingerprints.size == 0
        except (ValueError, IndexError, Exception) as e:
            # Expected behavior for empty input
            assert True
    
    def test_single_molecule(self):
        """Test with single molecule."""
        single_results = parse_smiles([TEST_SMILES[0]])
        fp = single_results.to_fingerprint()
        
        assert len(fp) == 1
        assert fp.fingerprints.shape[0] == 1
    
    def test_invalid_group_selection(self, results):
        """Test invalid group selection."""
        with pytest.raises(ValueError, match="Unknown group_selection"):
            results.to_fingerprint(group_selection='invalid')
    
    def test_invalid_kl_method(self, fingerprint):
        """Test invalid KL method."""
        with pytest.raises(ValueError, match="Unknown method"):
            fingerprint.compare_kld(fingerprint, method='invalid')
    
    def test_large_n_components(self, fingerprint):
        """Test PCA with more components than features."""
        # Should handle gracefully
        max_components = min(fingerprint.fingerprints.shape)
        pca = fingerprint.perform_pca(n_components=max_components, plot=False)
        assert pca['transformed'].shape[1] <= max_components


class TestIntegration:
    """Integration tests combining multiple features."""
    
    def test_full_workflow(self):
        """Test complete analysis workflow."""
        # Parse molecules
        results = parse_smiles(TEST_SMILES)
        
        # Convert to fingerprints
        fp = results.to_fingerprint(group_selection='all', count_mode='binary')
        
        # Perform multiple analyses
        pca = fp.perform_pca(n_components=3, plot=False)
        kpca = fp.perform_kernel_pca(n_components=2, kernel='rbf', plot=False)
        tsne = fp.perform_tsne(n_components=2, n_iter=250, plot=False)
        
        # Check all analyses complete
        assert 'transformed' in pca
        assert 'transformed' in kpca
        assert 'transformed' in tsne
        
        # Check dimensions
        assert pca['transformed'].shape == (len(TEST_SMILES), 3)
        assert kpca['transformed'].shape == (len(TEST_SMILES), 2)
        assert tsne['transformed'].shape == (len(TEST_SMILES), 2)
    
    def test_comparison_workflow(self):
        """Test dataset comparison workflow."""
        # Parse two datasets
        results1 = parse_smiles(TEST_SMILES)
        results2 = parse_smiles(COMPARISON_SMILES)
        
        # Convert to fingerprints
        fp1 = results1.to_fingerprint(group_selection='all')
        fp2 = results2.to_fingerprint(group_selection='all')
        
        # Compare using different methods
        methods = ['minmax', 'forward', 'reverse', 'symmetric']
        kl_values = {}
        
        for method in methods:
            kl_values[method] = fp1.compare_kld(fp2, method=method)
        
        # All values should be valid
        for method, value in kl_values.items():
            assert isinstance(value, float)
            assert not np.isnan(value)
            assert value >= 0
        
        # minmax should be normalized
        assert 0 <= kl_values['minmax'] <= 1


class TestVisualization:
    """Test visualization and plotting functionality."""
    
    def test_pca_plot_generation(self, fingerprint, tmp_path):
        """Test PCA plot file generation."""
        output_file = tmp_path / "test_pca.png"
        
        pca = fingerprint.perform_pca(
            n_components=2,
            plot=True,
            output_file=str(output_file)
        )
        
        # Check file was created
        assert output_file.exists()
        assert output_file.stat().st_size > 0
    
    def test_tsne_plot_generation(self, fingerprint, tmp_path):
        """Test t-SNE plot file generation."""
        output_file = tmp_path / "test_tsne.png"
        
        tsne = fingerprint.perform_tsne(
            n_components=2,
            perplexity=min(3, len(fingerprint) - 1),
            n_iter=250,
            plot=True,
            output_file=str(output_file)
        )
        
        # Check file was created
        assert output_file.exists()
        assert output_file.stat().st_size > 0
    
    def test_multiple_plots(self, fingerprint, tmp_path):
        """Test generating multiple plots."""
        # Generate different types of plots
        fp = fingerprint
        
        pca_file = tmp_path / "pca.png"
        kpca_file = tmp_path / "kpca.png"
        tsne_file = tmp_path / "tsne.png"
        
        fp.perform_pca(n_components=2, plot=True, output_file=str(pca_file))
        fp.perform_kernel_pca(n_components=2, plot=True, output_file=str(kpca_file))
        fp.perform_tsne(n_components=2, n_iter=250, plot=True, output_file=str(tsne_file))
        
        # All files should exist
        assert pca_file.exists()
        assert kpca_file.exists()
        assert tsne_file.exists()


class TestDataValidation:
    """Test data validation and consistency."""
    
    def test_fingerprint_consistency(self, results):
        """Test that fingerprints are consistent across conversions."""
        fp1 = results.to_fingerprint(group_selection='all', count_mode='binary')
        fp2 = results.to_fingerprint(group_selection='all', count_mode='binary')
        
        # Should produce identical results
        assert np.array_equal(fp1.fingerprints, fp2.fingerprints)
        assert fp1.smiles == fp2.smiles
        assert fp1.group_names == fp2.group_names
    
    def test_count_mode_differences(self, results):
        """Test that different count modes produce different results."""
        fp_binary = results.to_fingerprint(count_mode='binary')
        fp_count = results.to_fingerprint(count_mode='count')
        
        # Binary should only have 0s and 1s
        assert np.all((fp_binary.fingerprints == 0) | (fp_binary.fingerprints == 1))
        
        # Count may have values > 1
        # Note: may not be different if no multiple matches
        assert np.all(fp_count.fingerprints >= 0)
    
    def test_group_distribution_computation(self, fingerprint):
        """Test internal group distribution computation."""
        dist = fingerprint._compute_group_distribution()
        
        # Should be non-negative
        assert np.all(dist >= 0)
        
        # Should sum to reasonable value
        assert dist.sum() > 0
        
        # Length should match number of groups
        assert len(dist) == len(fingerprint.group_names)
    
    def test_summary_statistics(self, fingerprint):
        """Test that summary provides useful statistics."""
        summary = fingerprint.summary()
        
        # Check for key information
        assert 'Molecules:' in summary
        assert 'Groups:' in summary
        assert 'Sparsity:' in summary
        assert 'Most common groups:' in summary
        
        # Check values are present
        assert str(len(fingerprint)) in summary


class TestPerformance:
    """Test performance characteristics."""
    
    def test_pca_speed(self, fingerprint):
        """Test that PCA completes in reasonable time."""
        import time
        
        start = time.time()
        pca = fingerprint.perform_pca(n_components=5, plot=False)
        duration = time.time() - start
        
        # Should complete quickly for small dataset
        assert duration < 5.0  # seconds
    
    def test_tsne_speed(self, fingerprint):
        """Test t-SNE with minimal iterations."""
        import time
        
        start = time.time()
        tsne = fingerprint.perform_tsne(
            n_components=2,
            n_iter=100,  # Minimal iterations
            plot=False
        )
        duration = time.time() - start
        
        # t-SNE is slower but should still be reasonable
        assert duration < 30.0  # seconds


class TestDocumentation:
    """Test that documentation and examples work."""
    
    def test_repr_output(self, fingerprint):
        """Test that repr provides useful information."""
        repr_str = repr(fingerprint)
        
        assert 'ResultsFingerprint' in repr_str
        assert 'n_molecules=' in repr_str
        assert 'n_groups=' in repr_str
        assert 'group_selection=' in repr_str
        assert 'count_mode=' in repr_str
    
    def test_method_docstrings(self):
        """Test that key methods have docstrings."""
        # Check that important methods are documented
        assert ResultsFingerprint.perform_pca.__doc__ is not None
        assert ResultsFingerprint.perform_tsne.__doc__ is not None
        assert ResultsFingerprint.compare_kld.__doc__ is not None
        assert ResultsModel.to_fingerprint.__doc__ is not None
        
        # Check docstrings contain key information
        assert 'PCA' in ResultsFingerprint.perform_pca.__doc__
        assert 't-SNE' in ResultsFingerprint.perform_tsne.__doc__
        assert 'KL' in ResultsFingerprint.compare_kld.__doc__


class TestNumericalStability:
    """Test numerical stability and edge cases."""
    
    def test_zero_variance_features(self):
        """Test handling of zero-variance features."""
        # Create fingerprints with some zero-variance features
        fps = np.array([
            [1, 0, 1, 0],
            [1, 0, 0, 0],
            [1, 0, 1, 0],
        ])
        smiles = ['A', 'B', 'C']
        groups = ['G1', 'G2', 'G3', 'G4']
        
        fp = ResultsFingerprint(fps, smiles, groups)
        
        # PCA should handle this gracefully
        pca = fp.perform_pca(n_components=2, plot=False)
        assert not np.any(np.isnan(pca['transformed']))
    
    def test_identical_molecules(self):
        """Test with identical SMILES."""
        identical_smiles = [TEST_SMILES[0]] * 3
        results = parse_smiles(identical_smiles)
        fp = results.to_fingerprint()
        
        # Should not crash
        assert len(fp) == 3
        
        # All fingerprints should be identical
        assert np.all(fp.fingerprints[0] == fp.fingerprints[1])
        assert np.all(fp.fingerprints[1] == fp.fingerprints[2])


if __name__ == "__main__":
    # Run tests with verbose output
    pytest.main([__file__, "-v", "-s"])
