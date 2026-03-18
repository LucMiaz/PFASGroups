"""
PFAS Molecule Prioritization - Examples
========================================

This script demonstrates how to use the prioritise_molecules function
to rank PFAS compounds based on different criteria.
"""

from PFASGroups import prioritise_molecules, get_priority_statistics
import numpy as np


def example_reference_based():
    """Example 1: Prioritize by similarity to reference compounds."""
    print("=" * 80)
    print("EXAMPLE 1: Reference-Based Prioritization")
    print("=" * 80)
    
    # Chemical inventory to prioritize
    inventory = [
        "FC(F)(F)C(F)(F)C(=O)O",  # PFPA (C3)
        "FC(F)(F)C(F)(F)C(F)(F)C(=O)O",  # PFBA (C4)
        "FC(F)(F)C(F)(F)C(F)(F)C(F)(F)C(=O)O",  # PFPeA (C5)
        "FC(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(=O)O",  # PFHxA (C6)
        "FC(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(=O)O",  # PFHpA (C7)
        "FC(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(=O)O",  # PFOA (C8)
        "FC(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(=O)O",  # PFNA (C9)
        "OCC(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)F",  # 6:2 FTOH
        "FC(F)(F)C(F)(F)S(=O)(=O)O",  # PFES (C2 sulfonic acid)
        "FC(F)(F)C(F)(F)C(F)(F)C(F)(F)S(=O)(=O)O",  # PFBS (C4 sulfonic acid)
    ]
    
    # Known compounds of concern (e.g., regulated long-chain PFCAs)
    reference = [
        "FC(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(=O)O",  # PFHpA (C7)
        "FC(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(=O)O",  # PFOA (C8)
    ]
    
    print("\nPrioritizing by similarity to PFOA and PFHpA...")
    results, scores = prioritise_molecules(
        inventory,
        reference=reference,
        group_selection='all',
        count_mode='binary',
        return_scores=True
    )
    
    print(f"\nTop 5 Most Similar Compounds:")
    print(f"{'Rank':<6} {'Score':<8} {'SMILES':<60}")
    print("-" * 74)
    for i in range(min(5, len(results))):
        smiles = results[i].smiles[:57] + "..." if len(results[i].smiles) > 60 else results[i].smiles
        print(f"{i+1:<6} {scores[i]:<8.4f} {smiles}")
    
    # Get statistics
    stats = get_priority_statistics(results, scores, top_n=3)
    print(f"\nStatistics:")
    print(f"  Mean similarity score: {stats['score_mean']:.4f}")
    print(f"  Score range: {stats['score_min']:.4f} - {stats['score_max']:.4f}")
    print(f"\n  Most common groups in top 3:")
    for group_name, count in stats['top_n_groups'][:5]:
        print(f"    {group_name}: {count}")


def example_environmental_persistence():
    """Example 2: Prioritize for environmental persistence (long chains)."""
    print("\n" + "=" * 80)
    print("EXAMPLE 2: Environmental Persistence Priority")
    print("=" * 80)
    
    inventory = [
        "FC(F)(F)C(F)(F)C(=O)O",  # Short chain (C3)
        "FC(F)(F)C(F)(F)C(F)(F)C(F)(F)C(=O)O",  # Medium chain (C5)
        "FC(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(=O)O",  # Long chain (C8)
        "FC(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(=O)O",  # Very long (C10)
        "CC(F)(F)C(F)(F)F",  # Short, partially fluorinated
        "OCC(F)(F)C(F)(F)C(F)(F)C(F)(F)F",  # Fluorotelomer alcohol
    ]
    
    print("\nPrioritizing for long-chain compounds (environmental persistence)...")
    print("Using: a=0.5, b=2.0, percentile=95")
    
    results, scores = prioritise_molecules(
        inventory,
        reference=None,
        a=0.5,  # Lower weight on total fluorination
        b=2.0,  # Higher weight on longest components
        percentile=95,  # Focus on largest components
        return_scores=True
    )
    
    print(f"\nTop 3 Highest Priority (Longest Chains):")
    print(f"{'Rank':<6} {'Score':<10} {'SMILES':<60}")
    print("-" * 76)
    for i in range(min(3, len(results))):
        smiles = results[i].smiles[:57] + "..." if len(results[i].smiles) > 60 else results[i].smiles
        print(f"{i+1:<6} {scores[i]:<10.2f} {smiles}")


def example_bioaccumulation():
    """Example 3: Prioritize for bioaccumulation potential."""
    print("\n" + "=" * 80)
    print("EXAMPLE 3: Bioaccumulation Potential")
    print("=" * 80)
    
    inventory = [
        "FC(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(=O)O",  # PFHxA
        "FC(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(=O)O",  # PFHpA
        "FC(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(=O)O",  # PFOA
        "FC(F)(F)C(F)(F)C(F)(F)C(F)(F)S(=O)(=O)O",  # PFBS
        "FC(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)S(=O)(=O)O",  # PFHxS
        "OCC(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)F",  # 8:2 FTOH
    ]
    
    print("\nPrioritizing for bioaccumulation (balanced total F and chain length)...")
    print("Using: a=1.0, b=1.0, percentile=75")
    
    results, scores = prioritise_molecules(
        inventory,
        reference=None,
        a=1.0,  # Equal weight on total fluorination
        b=1.0,  # Equal weight on chain length
        percentile=75,  # Upper quartile
        return_scores=True
    )
    
    print(f"\nPriority Ranking:")
    print(f"{'Rank':<6} {'Score':<10} {'SMILES':<60}")
    print("-" * 76)
    for i in range(len(results)):
        smiles = results[i].smiles[:57] + "..." if len(results[i].smiles) > 60 else results[i].smiles
        print(f"{i+1:<6} {scores[i]:<10.2f} {smiles}")
    
    stats = get_priority_statistics(results, scores)
    print(f"\nScore statistics: μ={stats['score_mean']:.2f}, σ={stats['score_std']:.2f}")


def example_total_fluorination():
    """Example 4: Prioritize by total fluorine content."""
    print("\n" + "=" * 80)
    print("EXAMPLE 4: Total Fluorination Priority")
    print("=" * 80)
    
    inventory = [
        "FC(F)(F)C(F)(F)C(=O)O",  # Few F atoms
        "FC(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(=O)O",  # More F atoms
        "FC(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(=O)O",  # Many F atoms
        "CC(F)(F)C(F)(F)F",  # Partially fluorinated
    ]
    
    print("\nPrioritizing for total fluorine content...")
    print("Using: a=2.0, b=0.5, percentile=50")
    
    results, scores = prioritise_molecules(
        inventory,
        reference=None,
        a=2.0,  # High weight on total fluorination
        b=0.5,  # Low weight on individual chain length
        percentile=50,  # Median
        return_scores=True
    )
    
    print(f"\nRanking by Total Fluorination:")
    print(f"{'Rank':<6} {'Score':<10} {'SMILES':<60}")
    print("-" * 76)
    for i in range(len(results)):
        smiles = results[i].smiles
        print(f"{i+1:<6} {scores[i]:<10.2f} {smiles}")


def example_comparison():
    """Example 5: Compare different prioritization strategies."""
    print("\n" + "=" * 80)
    print("EXAMPLE 5: Comparing Prioritization Strategies")
    print("=" * 80)
    
    inventory = [
        "FC(F)(F)C(F)(F)C(=O)O",  # C3
        "FC(F)(F)C(F)(F)C(F)(F)C(F)(F)C(=O)O",  # C5
        "FC(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(=O)O",  # C7
        "OCC(F)(F)C(F)(F)C(F)(F)C(F)(F)F",  # FTOH
    ]
    
    strategies = [
        ("Total Fluorination", {"a": 2.0, "b": 0.5, "percentile": 50}),
        ("Long Chains", {"a": 0.5, "b": 2.0, "percentile": 95}),
        ("Balanced", {"a": 1.0, "b": 1.0, "percentile": 75}),
    ]
    
    print(f"\n{'Molecule':<50} ", end="")
    for strategy_name, _ in strategies:
        print(f"{strategy_name:<20}", end="")
    print()
    print("-" * 110)
    
    # Store results for each strategy
    all_results = {}
    for strategy_name, params in strategies:
        results, scores = prioritise_molecules(
            inventory,
            reference=None,
            **params,
            return_scores=True
        )
        all_results[strategy_name] = (results, scores)
    
    # Display comparison
    for i, smiles in enumerate(inventory):
        short_smiles = smiles[:47] + "..." if len(smiles) > 50 else smiles
        print(f"{short_smiles:<50} ", end="")
        
        for strategy_name, _ in strategies:
            results, scores = all_results[strategy_name]
            # Find this molecule's score
            for j, result in enumerate(results):
                if result.smiles == smiles:
                    rank = j + 1
                    score = scores[j]
                    print(f"#{rank} ({score:.1f}){' ':<10}", end="")
                    break
        print()


def main():
    """Run all examples."""
    print("\n" + "=" * 80)
    print("PFAS MOLECULE PRIORITIZATION EXAMPLES")
    print("PFASGroups v3.2.0")
    print("=" * 80)
    
    example_reference_based()
    example_environmental_persistence()
    example_bioaccumulation()
    example_total_fluorination()
    example_comparison()
    
    print("\n" + "=" * 80)
    print("EXAMPLES COMPLETE")
    print("=" * 80)
    print("\nFor more information, see:")
    print("  - docs/prioritization.rst")
    print("  - help(prioritise_molecules)")
    print()


if __name__ == '__main__':
    main()
