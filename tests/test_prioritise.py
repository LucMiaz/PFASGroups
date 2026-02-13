"""
Tests for PFAS molecule prioritization functionality.
"""

import pytest
import numpy as np
from rdkit import Chem

from PFASgroups import parse_smiles, prioritise_molecules, prioritize_molecules, get_priority_statistics
from PFASgroups.results_model import ResultsModel


# Test molecules with varying fluorination
TEST_SMILES = [
    "FC(F)(F)C(F)(F)C(=O)O",  # PFPA (C3)
    "FC(F)(F)C(F)(F)C(F)(F)C(F)(F)C(=O)O",  # PFBA (C5)
    "FC(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(=O)O",  # PFOA (C8)
    "CC(F)(F)C(F)(F)F",  # Short chain (C3, one H)
    "FC(F)(F)C(F)(F)C(F)(F)F",  # Medium chain (C4)
]

REFERENCE_SMILES = [
    "FC(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(=O)O",  # PFHxA (C6)
    "FC(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(=O)O",  # PFHpA (C7)
]


class TestPrioritiseMolecules:
    """Test prioritise_molecules function."""
    
    def test_basic_reference_prioritization(self):
        """Test basic prioritization with reference list."""
        results, scores = prioritise_molecules(
            TEST_SMILES,
            reference=REFERENCE_SMILES,
            return_scores=True
        )
        
        assert isinstance(results, ResultsModel)
        assert len(results) == len(TEST_SMILES)
        assert len(scores) == len(TEST_SMILES)
        assert np.all(scores >= 0) and np.all(scores <= 1)
        
        # PFOA (C8) should be most similar to reference (C6-C7)
        top_smiles = results[0]['smiles']
        assert 'C(=O)O' in top_smiles  # Should be a carboxylic acid
    
    def test_intrinsic_prioritization(self):
        """Test prioritization without reference (intrinsic properties)."""
        results, scores = prioritise_molecules(
            TEST_SMILES,
            reference=None,
            a=1.0,
            b=1.0,
            percentile=90,
            return_scores=True
        )
        
        assert isinstance(results, ResultsModel)
        assert len(results) == len(TEST_SMILES)
        assert len(scores) == len(TEST_SMILES)
        assert np.all(scores >= 0)
        
        # PFOA (C8) should have highest score (longest chain)
        top_smiles = results[0]['smiles']
        # Should be the longest molecule
        assert len(top_smiles) > 40  # PFOA is long
    
    def test_weight_parameters(self):
        """Test that weight parameters affect prioritization."""
        # Emphasize total fluorination
        results_a, scores_a = prioritise_molecules(
            TEST_SMILES,
            a=2.0,
            b=0.5,
            return_scores=True
        )
        
        # Emphasize longest components
        results_b, scores_b = prioritise_molecules(
            TEST_SMILES,
            a=0.5,
            b=2.0,
            percentile=95,
            return_scores=True
        )
        
        # Scores should be different
        assert not np.allclose(scores_a, scores_b)
    
    def test_percentile_parameter(self):
        """Test percentile parameter affects results."""
        results_90, scores_90 = prioritise_molecules(
            TEST_SMILES,
            a=1.0,
            b=1.0,
            percentile=90,
            return_scores=True
        )
        
        results_50, scores_50 = prioritise_molecules(
            TEST_SMILES,
            a=1.0,
            b=1.0,
            percentile=50,
            return_scores=True
        )
        
        # Scores should be different
        assert not np.allclose(scores_90, scores_50)
    
    def test_ascending_order(self):
        """Test ascending vs descending sort order."""
        results_desc, scores_desc = prioritise_molecules(
            TEST_SMILES,
            ascending=False,
            return_scores=True
        )
        
        results_asc, scores_asc = prioritise_molecules(
            TEST_SMILES,
            ascending=True,
            return_scores=True
        )
        
        # Scores should be in reverse order
        assert np.allclose(scores_desc, scores_asc[::-1])
        
        # First molecule should be different
        assert results_desc[0]['smiles'] != results_asc[0]['smiles']
    
    def test_return_scores_flag(self):
        """Test return_scores parameter."""
        # With scores
        result_with = prioritise_molecules(
            TEST_SMILES,
            return_scores=True
        )
        assert isinstance(result_with, tuple)
        assert len(result_with) == 2
        assert isinstance(result_with[0], ResultsModel)
        assert isinstance(result_with[1], np.ndarray)
        
        # Without scores
        result_without = prioritise_molecules(
            TEST_SMILES,
            return_scores=False
        )
        assert isinstance(result_without, ResultsModel)
    
    def test_results_model_input(self):
        """Test using ResultsModel as input."""
        results = parse_smiles(TEST_SMILES)
        
        prioritized, scores = prioritise_molecules(
            results,
            return_scores=True
        )
        
        assert isinstance(prioritized, ResultsModel)
        assert len(prioritized) == len(results)
    
    def test_mol_object_input(self):
        """Test using RDKit Mol objects as input."""
        mols = [Chem.MolFromSmiles(smi) for smi in TEST_SMILES]
        
        results, scores = prioritise_molecules(
            mols,
            return_scores=True
        )
        
        assert isinstance(results, ResultsModel)
        assert len(results) == len(mols)
    
    def test_group_selection(self):
        """Test different group selections for fingerprint."""
        results_all, scores_all = prioritise_molecules(
            TEST_SMILES,
            reference=REFERENCE_SMILES,
            group_selection='all',
            return_scores=True
        )
        
        results_oecd, scores_oecd = prioritise_molecules(
            TEST_SMILES,
            reference=REFERENCE_SMILES,
            group_selection='oecd',
            return_scores=True
        )
        
        # Results might be similar but not identical
        assert len(results_all) == len(results_oecd)
    
    def test_count_mode(self):
        """Test different count modes for fingerprint."""
        results_binary, scores_binary = prioritise_molecules(
            TEST_SMILES,
            reference=REFERENCE_SMILES,
            count_mode='binary',
            return_scores=True
        )
        
        results_count, scores_count = prioritise_molecules(
            TEST_SMILES,
            reference=REFERENCE_SMILES,
            count_mode='count',
            return_scores=True
        )
        
        assert len(results_binary) == len(results_count)
    
    def test_empty_input(self):
        """Test error handling for empty input."""
        with pytest.raises(ValueError):
            prioritise_molecules([])
    
    def test_invalid_input_type(self):
        """Test error handling for invalid input types."""
        with pytest.raises(TypeError):
            prioritise_molecules("not a list or ResultsModel")
    
    def test_invalid_list_contents(self):
        """Test error handling for invalid list contents."""
        with pytest.raises(TypeError):
            prioritise_molecules([123, 456])  # Numbers, not SMILES or Mols
    
    def test_american_spelling_alias(self):
        """Test that American spelling works (prioritize vs prioritise)."""
        results1, scores1 = prioritise_molecules(TEST_SMILES, return_scores=True)
        results2, scores2 = prioritize_molecules(TEST_SMILES, return_scores=True)
        
        # Should give identical results
        assert len(results1) == len(results2)
        assert np.allclose(scores1, scores2)


class TestGetPriorityStatistics:
    """Test get_priority_statistics function."""
    
    def test_basic_statistics(self):
        """Test basic statistics computation."""
        results, scores = prioritise_molecules(TEST_SMILES, return_scores=True)
        stats = get_priority_statistics(results, scores, top_n=3)
        
        assert 'n_molecules' in stats
        assert stats['n_molecules'] == len(TEST_SMILES)
        
        assert 'score_mean' in stats
        assert 'score_std' in stats
        assert 'score_min' in stats
        assert 'score_max' in stats
        assert 'score_median' in stats
        
        assert stats['score_min'] <= stats['score_mean'] <= stats['score_max']
    
    def test_top_n_information(self):
        """Test top N molecule information."""
        results, scores = prioritise_molecules(TEST_SMILES, return_scores=True)
        stats = get_priority_statistics(results, scores, top_n=3)
        
        assert 'top_n_smiles' in stats
        assert len(stats['top_n_smiles']) == 3
        
        assert 'top_n_scores' in stats
        assert len(stats['top_n_scores']) == 3
        
        # First score should be highest
        assert stats['top_n_scores'][0] >= stats['top_n_scores'][1]
        assert stats['top_n_scores'][1] >= stats['top_n_scores'][2]
    
    def test_group_analysis(self):
        """Test group frequency analysis."""
        results, scores = prioritise_molecules(TEST_SMILES, return_scores=True)
        stats = get_priority_statistics(results, scores, top_n=5)
        
        assert 'top_n_groups' in stats
        assert isinstance(stats['top_n_groups'], list)
        
        # Should contain (group_name, count) tuples
        if len(stats['top_n_groups']) > 0:
            assert isinstance(stats['top_n_groups'][0], tuple)
            assert len(stats['top_n_groups'][0]) == 2
    
    def test_top_n_larger_than_dataset(self):
        """Test when top_n is larger than dataset."""
        results, scores = prioritise_molecules(TEST_SMILES, return_scores=True)
        stats = get_priority_statistics(results, scores, top_n=100)
        
        # Should only return available molecules
        assert len(stats['top_n_smiles']) == len(TEST_SMILES)


class TestIntegrationScenarios:
    """Integration tests for real-world scenarios."""
    
    def test_environmental_persistence_ranking(self):
        """Test ranking for environmental persistence (long chains)."""
        # Emphasize longest chains
        results, scores = prioritise_molecules(
            TEST_SMILES,
            a=0.5,
            b=2.0,
            percentile=95,
            return_scores=True
        )
        
        # Longer chains should score higher
        assert scores[0] > scores[-1]
        
        stats = get_priority_statistics(results, scores, top_n=2)
        assert stats['score_max'] > stats['score_min']
    
    def test_bioaccumulation_screening(self):
        """Test ranking for bioaccumulation potential."""
        # Balance total fluorination and chain length
        results, scores = prioritise_molecules(
            TEST_SMILES,
            a=1.0,
            b=1.0,
            percentile=75,
            return_scores=True
        )
        
        # Should prioritize well-fluorinated molecules
        assert np.all(scores >= 0)
        
        stats = get_priority_statistics(results, scores)
        assert stats['n_molecules'] == len(TEST_SMILES)
    
    def test_similarity_to_regulated_compounds(self):
        """Test similarity-based ranking to known regulated compounds."""
        # Use PFOA-like molecules as reference
        reference = ["FC(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(=O)O"]
        
        results, scores = prioritise_molecules(
            TEST_SMILES,
            reference=reference,
            return_scores=True
        )
        
        # Molecules similar to PFOA should score high
        assert np.all(scores >= 0) and np.all(scores <= 1)
        
        stats = get_priority_statistics(results, scores, top_n=3)
        # Top molecules should have reasonable similarity
        assert stats['top_n_scores'][0] >= 0
    
    def test_mixed_inventory_prioritization(self):
        """Test prioritization of mixed PFAS inventory."""
        # Mix of different chain lengths and functional groups
        mixed_inventory = [
            "FC(F)(F)C(F)(F)S(=O)(=O)O",  # Sulfonic acid
            "FC(F)(F)C(F)(F)C(F)(F)C(F)(F)C(=O)O",  # Carboxylic acid
            "OCC(F)(F)C(F)(F)C(F)(F)F",  # Alcohol
            "FC(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)F",  # Alkane
        ]
        
        results, scores = prioritise_molecules(
            mixed_inventory,
            a=1.0,
            b=1.0,
            percentile=90,
            return_scores=True
        )
        
        assert len(results) == len(mixed_inventory)
        assert len(scores) == len(mixed_inventory)
        
        # Check that prioritization is stable
        stats = get_priority_statistics(results, scores)
        assert stats['n_molecules'] == len(mixed_inventory)


class TestEdgeCases:
    """Test edge cases and error conditions."""
    
    def test_single_molecule(self):
        """Test with single molecule."""
        results, scores = prioritise_molecules(
            [TEST_SMILES[0]],
            return_scores=True
        )
        
        assert len(results) == 1
        assert len(scores) == 1
    
    def test_molecules_with_no_fluorination(self):
        """Test with non-fluorinated molecules."""
        non_pfas = ["CC(=O)O", "CCCCCC"]  # Acetic acid, hexane
        
        results, scores = prioritise_molecules(
            non_pfas,
            return_scores=True
        )
        
        # Should handle gracefully (score = 0)
        assert len(results) == len(non_pfas)
        # Scores might be zero or very low
        assert np.all(scores >= 0)
    
    def test_identical_molecules(self):
        """Test with identical molecules."""
        identical = [TEST_SMILES[0]] * 3
        
        results, scores = prioritise_molecules(
            identical,
            return_scores=True
        )
        
        # All scores should be identical
        assert np.allclose(scores, scores[0])
    
    def test_zero_weights(self):
        """Test with zero weights."""
        results, scores = prioritise_molecules(
            TEST_SMILES,
            a=0.0,
            b=0.0,
            return_scores=True
        )
        
        # All scores should be zero
        assert np.allclose(scores, 0.0)


def test_example_from_docstring():
    """Test the example from the docstring."""
    inventory = ["FC(F)(F)C(F)(F)C(=O)O", "FC(F)(F)C(F)(F)C(F)(F)C(=O)O"]
    reference = ["FC(F)(F)C(F)(F)C(=O)O"]
    
    results, scores = prioritise_molecules(inventory, reference=reference)
    
    assert len(results) == 2
    assert len(scores) == 2
    assert results[0]['smiles'] in inventory


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
