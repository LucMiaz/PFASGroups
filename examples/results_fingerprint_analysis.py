"""
Example: ResultsFingerprint Analysis with Dimensionality Reduction

This script demonstrates how to:
1. Convert ResultsModel to ResultsFingerprint
2. Perform dimensionality reduction (PCA, kernel-PCA, t-SNE, UMAP)
3. Compare fingerprints using KL divergence
4. Save/load fingerprints to/from SQL
"""

import numpy as np
from HalogenGroups import parse_smiles

# Example PFAS molecules
smiles_list_1 = [
    "C(C(C(C(C(C(C(C(F)(F)F)(F)F)(F)F)(F)F)(F)F)(F)F)(F)F)(C(=O)O)F",  # PFOA
    "C(C(C(C(C(C(C(F)(F)F)(F)F)(F)F)(F)F)(F)F)(F)F)(C(=O)O)F",  # PFHpA
    "C(C(C(C(C(F)(F)F)(F)F)(F)F)(F)F)(C(=O)O)F",  # PFBA
    "C(C(C(C(C(C(F)(F)F)(F)F)(F)F)(F)F)(F)F)(S(=O)(=O)O)F",  # PFHxS
    "FC(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)S(=O)(=O)O",  # PFOS
]

smiles_list_2 = [
    "C(C(C(C(F)(F)F)(F)F)(F)F)(S(=O)(=O)O)F",  # PFBS
    "FC(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(=O)O",  # PFNA
    "FC(C(F)(F)F)(F)c1ccccc1",  # Fluorinated aromatic
    "C(C(C(F)(F)F)(F)F)(O)F",  # Fluorotelomer alcohol
    "FC(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(=O)O",  # PFDA
]


def example_basic_conversion():
    """Example 1: Convert results to fingerprints with different group selections."""
    print("=" * 80)
    print("Example 1: Converting ResultsModel to ResultsFingerprint")
    print("=" * 80)
    
    # Parse SMILES
    results = parse_smiles(smiles_list_1)
    print(f"Parsed {len(results)} molecules")
    
    # Convert to fingerprints with different group selections
    print("\n1. All groups:")
    fp_all = results.to_fingerprint(group_selection='all', count_mode='binary')
    print(fp_all.summary())
    
    print("\n2. OECD groups only:")
    fp_oecd = results.to_fingerprint(group_selection='oecd', count_mode='binary')
    print(fp_oecd)
    
    print("\n3. Generic groups only:")
    fp_generic = results.to_fingerprint(group_selection='generic', count_mode='binary')
    print(fp_generic)
    
    print("\n4. Count mode (instead of binary):")
    fp_count = results.to_fingerprint(group_selection='all', count_mode='count')
    print(fp_count)
    print(f"Sample fingerprint values: {fp_count.fingerprints[0][:10]}")
    
    return fp_all, fp_oecd


def example_pca_analysis(fp):
    """Example 2: Perform PCA analysis."""
    print("\n" + "=" * 80)
    print("Example 2: PCA Analysis")
    print("=" * 80)
    
    # Perform PCA with 5 components
    pca_results = fp.perform_pca(
        n_components=5,
        plot=True,
        output_file="pfas_pca_analysis.png"
    )
    
    print(f"Explained variance ratio: {pca_results['explained_variance']}")
    print(f"Cumulative variance: {np.cumsum(pca_results['explained_variance'])}")
    print(f"Transformed data shape: {pca_results['transformed'].shape}")
    
    return pca_results


def example_kernel_pca(fp):
    """Example 3: Perform kernel PCA analysis."""
    print("\n" + "=" * 80)
    print("Example 3: Kernel PCA Analysis")
    print("=" * 80)
    
    # Try different kernels
    for kernel in ['rbf', 'poly', 'sigmoid']:
        print(f"\n{kernel.upper()} kernel:")
        kpca_results = fp.perform_kernel_pca(
            n_components=2,
            kernel=kernel,
            plot=True,
            output_file=f"pfas_kpca_{kernel}.png"
        )
        print(f"Transformed data shape: {kpca_results['transformed'].shape}")


def example_tsne(fp):
    """Example 4: Perform t-SNE analysis."""
    print("\n" + "=" * 80)
    print("Example 4: t-SNE Analysis")
    print("=" * 80)
    
    # t-SNE with default parameters
    tsne_results = fp.perform_tsne(
        n_components=2,
        perplexity=30.0,
        plot=True,
        output_file="pfas_tsne_analysis.png"
    )
    
    print(f"Transformed data shape: {tsne_results['transformed'].shape}")
    print(f"Perplexity used: {tsne_results['perplexity']}")
    
    return tsne_results


def example_umap(fp):
    """Example 5: Perform UMAP analysis."""
    print("\n" + "=" * 80)
    print("Example 5: UMAP Analysis")
    print("=" * 80)
    
    try:
        # UMAP with default parameters
        umap_results = fp.perform_umap(
            n_components=2,
            n_neighbors=15,
            plot=True,
            output_file="pfas_umap_analysis.png"
        )
        
        print(f"Transformed data shape: {umap_results['transformed'].shape}")
        print(f"n_neighbors: {umap_results['n_neighbors']}")
        print(f"min_dist: {umap_results['min_dist']}")
        
        return umap_results
    except ImportError:
        print("UMAP not installed. Install with: pip install umap-learn")
        return None


def example_kl_divergence():
    """Example 6: Compare fingerprints using KL divergence."""
    print("\n" + "=" * 80)
    print("Example 6: Comparing Fingerprints with KL Divergence")
    print("=" * 80)
    
    # Create two fingerprint sets
    results_1 = parse_smiles(smiles_list_1)
    results_2 = parse_smiles(smiles_list_2)
    
    fp_1 = results_1.to_fingerprint(group_selection='all', count_mode='binary')
    fp_2 = results_2.to_fingerprint(group_selection='all', count_mode='binary')
    
    print(f"Fingerprint 1: {len(fp_1)} molecules")
    print(f"Fingerprint 2: {len(fp_2)} molecules")
    
    # Compare using different methods
    methods = ['minmax', 'forward', 'reverse', 'symmetric']
    print("\nKL Divergence comparisons:")
    for method in methods:
        kld = fp_1.compare_kld(fp_2, method=method)
        print(f"  {method:10s}: {kld:.6f}")
    
    # Interpretation
    kld_minmax = fp_1.compare_kld(fp_2, method='minmax')
    if kld_minmax < 0.1:
        print(f"\nInterpretation: Very similar distributions (KL={kld_minmax:.4f})")
    elif kld_minmax < 0.3:
        print(f"\nInterpretation: Moderately similar distributions (KL={kld_minmax:.4f})")
    else:
        print(f"\nInterpretation: Different distributions (KL={kld_minmax:.4f})")


def example_sql_save_load(fp):
    """Example 7: Save and load fingerprints to/from SQL."""
    print("\n" + "=" * 80)
    print("Example 7: SQL Save/Load")
    print("=" * 80)
    
    # Save to SQLite database
    db_file = "pfas_fingerprints.db"
    print(f"\nSaving fingerprints to {db_file}...")
    fp.to_sql(filename=db_file, if_exists='replace')
    
    # Load back from database
    print(f"\nLoading fingerprints from {db_file}...")
    fp_loaded = fp.__class__.from_sql(filename=db_file)
    
    print(f"\nOriginal fingerprints: {fp}")
    print(f"Loaded fingerprints: {fp_loaded}")
    
    # Verify they match
    if np.allclose(fp.fingerprints, fp_loaded.fingerprints):
        print("✓ Fingerprints match after save/load!")
    else:
        print("✗ Warning: Fingerprints don't match after save/load")
    
    return fp_loaded


def example_results_model_sql():
    """Example 8: Save and load ResultsModel to/from SQL."""
    print("\n" + "=" * 80)
    print("Example 8: ResultsModel SQL Save/Load")
    print("=" * 80)
    
    # Parse and save
    results = parse_smiles(smiles_list_1[:3])  # Use subset for demo
    
    db_file = "pfas_results.db"
    print(f"\nSaving results to {db_file}...")
    results.to_sql(filename=db_file, if_exists='replace')
    
    # Load back
    print(f"\nLoading results from {db_file}...")
    from PFASGroups.PFASEmbeddings import ResultsModel
    results_loaded = ResultsModel.from_sql(filename=db_file)
    
    print(f"\nOriginal results: {len(results)} molecules")
    print(f"Loaded results: {len(results_loaded)} molecules")
    
    for i, (orig, loaded) in enumerate(zip(results, results_loaded)):
        print(f"Molecule {i+1}:")
        print(f"  Original SMILES: {orig.smiles}")
        print(f"  Loaded SMILES: {loaded.smiles}")
        print(f"  Original matches: {len(orig.matches)}")
        print(f"  Loaded matches: {len(loaded.matches)}")


def main():
    """Run all examples."""
    print("\n" + "=" * 80)
    print("PFAS ResultsFingerprint Analysis Examples")
    print("=" * 80)
    
    # Example 1: Basic conversion
    fp_all, fp_oecd = example_basic_conversion()
    
    # Example 2: PCA
    pca_results = example_pca_analysis(fp_all)
    
    # Example 3: Kernel PCA
    example_kernel_pca(fp_all)
    
    # Example 4: t-SNE
    tsne_results = example_tsne(fp_all)
    
    # Example 5: UMAP
    umap_results = example_umap(fp_all)
    
    # Example 6: KL divergence comparison
    example_kl_divergence()
    
    # Example 7: SQL save/load for fingerprints
    fp_loaded = example_sql_save_load(fp_all)
    
    # Example 8: SQL save/load for ResultsModel
    example_results_model_sql()
    
    print("\n" + "=" * 80)
    print("All examples completed successfully!")
    print("=" * 80)
    print("\nGenerated files:")
    print("  - pfas_pca_analysis.png")
    print("  - pfas_kpca_*.png")
    print("  - pfas_tsne_analysis.png")
    print("  - pfas_umap_analysis.png (if umap-learn is installed)")
    print("  - pfas_fingerprints.db")
    print("  - pfas_results.db")


if __name__ == "__main__":
    main()
