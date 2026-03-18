"""
Example: PFASEmbeddingSet Analysis with Dimensionality Reduction

This script demonstrates how to:
1. Parse SMILES strings with parse_smiles() -> PFASEmbeddingSet
2. Generate embedding arrays with to_array() (replaces deprecated to_fingerprint())
3. Perform dimensionality reduction (PCA, kernel-PCA, t-SNE, UMAP)
4. Compare embedding sets using KL divergence
5. Save/load results to/from SQL
"""

import numpy as np
from PFASGroups import parse_smiles, PFASEmbeddingSet

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


def example_basic_embedding():
    """Example 1: Parse SMILES and generate embeddings with different configurations."""
    print("=" * 80)
    print("Example 1: Parsing SMILES and Generating Embeddings")
    print("=" * 80)

    # Parse SMILES -> PFASEmbeddingSet
    results = parse_smiles(smiles_list_1)
    print(f"Parsed {len(results)} molecules")

    # Inspect the first molecule (PFASEmbedding has a .summary() method)
    print("\n1. Single-molecule summary:")
    results[0].summary()

    # Generate embedding arrays. to_array() is the primary method;
    # the old to_fingerprint() is deprecated.
    print("\n2. Binary embedding (all groups):")
    arr_all = results.to_array(component_metrics=['binary'])
    print(f"  Shape: {arr_all.shape}  (n_molecules × n_groups)")

    print("\n3. OECD groups only:")
    arr_oecd = results.to_array(component_metrics=['binary'], group_selection='oecd')
    print(f"  Shape: {arr_oecd.shape}")

    print("\n4. Generic groups only:")
    arr_generic = results.to_array(component_metrics=['binary'], group_selection='generic')
    print(f"  Shape: {arr_generic.shape}")

    print("\n5. Count mode (instead of binary):")
    arr_count = results.to_array(component_metrics=['count'])
    print(f"  Shape: {arr_count.shape}")
    print(f"  Sample row (first 10 values): {arr_count[0][:10]}")

    return results


def example_pca_analysis(results):
    """Example 2: Perform PCA analysis."""
    print("\n" + "=" * 80)
    print("Example 2: PCA Analysis")
    print("=" * 80)

    # perform_pca() is called directly on PFASEmbeddingSet; it calls to_array() internally.
    pca_results = results.perform_pca(
        n_components=5,
        plot=True,
        output_file="pfas_pca_analysis.png"
    )

    print(f"Explained variance ratio: {pca_results['explained_variance']}")
    print(f"Cumulative variance: {np.cumsum(pca_results['explained_variance'])}")
    print(f"Transformed data shape: {pca_results['transformed'].shape}")

    return pca_results


def example_kernel_pca(results):
    """Example 3: Perform kernel PCA analysis."""
    print("\n" + "=" * 80)
    print("Example 3: Kernel PCA Analysis")
    print("=" * 80)

    for kernel in ['rbf', 'poly', 'sigmoid']:
        print(f"\n{kernel.upper()} kernel:")
        kpca_results = results.perform_kernel_pca(
            n_components=2,
            kernel=kernel,
            plot=True,
            output_file=f"pfas_kpca_{kernel}.png"
        )
        print(f"Transformed data shape: {kpca_results['transformed'].shape}")


def example_tsne(results):
    """Example 4: Perform t-SNE analysis."""
    print("\n" + "=" * 80)
    print("Example 4: t-SNE Analysis")
    print("=" * 80)

    tsne_results = results.perform_tsne(
        n_components=2,
        perplexity=30.0,
        plot=True,
        output_file="pfas_tsne_analysis.png"
    )

    print(f"Transformed data shape: {tsne_results['transformed'].shape}")
    print(f"Perplexity used: {tsne_results['perplexity']}")

    return tsne_results


def example_umap(results):
    """Example 5: Perform UMAP analysis."""
    print("\n" + "=" * 80)
    print("Example 5: UMAP Analysis")
    print("=" * 80)

    try:
        umap_results = results.perform_umap(
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
    """Example 6: Compare two PFASEmbeddingSets using KL divergence."""
    print("\n" + "=" * 80)
    print("Example 6: Comparing Embedding Sets with KL Divergence")
    print("=" * 80)

    # compare_kld() is called directly on PFASEmbeddingSet
    results_1 = parse_smiles(smiles_list_1)
    results_2 = parse_smiles(smiles_list_2)

    print(f"Set 1: {len(results_1)} molecules")
    print(f"Set 2: {len(results_2)} molecules")

    methods = ['minmax', 'forward', 'reverse', 'symmetric']
    print("\nKL Divergence comparisons:")
    for method in methods:
        kld = results_1.compare_kld(results_2, method=method)
        print(f"  {method:10s}: {kld:.6f}")

    kld_minmax = results_1.compare_kld(results_2, method='minmax')
    if kld_minmax < 0.1:
        print(f"\nInterpretation: Very similar distributions (KL={kld_minmax:.4f})")
    elif kld_minmax < 0.3:
        print(f"\nInterpretation: Moderately similar distributions (KL={kld_minmax:.4f})")
    else:
        print(f"\nInterpretation: Different distributions (KL={kld_minmax:.4f})")


def example_sql_save_load(results):
    """Example 7: Save and load a PFASEmbeddingSet to/from SQLite."""
    print("\n" + "=" * 80)
    print("Example 7: SQL Save/Load")
    print("=" * 80)

    db_file = "pfas_results.db"
    print(f"\nSaving {len(results)} molecules to {db_file}...")
    results.to_sql(filename=db_file, if_exists='replace')

    print(f"\nLoading from {db_file}...")
    loaded = PFASEmbeddingSet.from_sql(filename=db_file)
    print(f"Loaded {len(loaded)} molecules")

    for i, mol in enumerate(loaded):
        print(f"  Molecule {i+1}: {mol.smiles}  (matches: {len(mol.matches)})")

    return loaded


def example_per_molecule_access():
    """Example 8: Per-molecule access via PFASEmbedding."""
    print("\n" + "=" * 80)
    print("Example 8: Per-Molecule Access (PFASEmbedding)")
    print("=" * 80)

    results = parse_smiles(smiles_list_1[:3])

    for mol in results:
        print(f"\nSMILES : {mol.smiles}")
        category, total_size = mol.classify()
        print(f"  Category        : {category}")
        print(f"  Total chain size: {total_size}")
        print(f"  Group matches   : {sum(1 for m in mol.matches if m.is_group)}")
        vec = mol.to_array(component_metrics=['binary'])
        print(f"  Active groups   : {int(vec.sum())} / {len(vec)}")


def main():
    """Run all examples."""
    print("\n" + "=" * 80)
    print("PFASGroups PFASEmbeddingSet Analysis Examples")
    print("=" * 80)

    # Example 1: Parse and embed
    results_1 = example_basic_embedding()

    # Example 2: PCA
    example_pca_analysis(results_1)

    # Example 3: Kernel PCA
    example_kernel_pca(results_1)

    # Example 4: t-SNE
    example_tsne(results_1)

    # Example 5: UMAP
    example_umap(results_1)

    # Example 6: KL divergence
    example_kl_divergence()

    # Example 7: SQL save/load
    example_sql_save_load(results_1)

    # Example 8: Per-molecule API
    example_per_molecule_access()

    print("\n" + "=" * 80)
    print("All examples completed successfully!")
    print("=" * 80)
    print("\nGenerated files:")
    print("  - pfas_pca_analysis.png")
    print("  - pfas_kpca_*.png")
    print("  - pfas_tsne_analysis.png")
    print("  - pfas_umap_analysis.png  (if umap-learn is installed)")
    print("  - pfas_results.db")


if __name__ == "__main__":
    main()
