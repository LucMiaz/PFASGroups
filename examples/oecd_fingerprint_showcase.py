"""
OECD PFAS Dataset - ResultsFingerprint Showcase
================================================

This script demonstrates the new ResultsFingerprint capabilities on the OECD PFAS dataset,
showcasing dimensionality reduction, statistical comparison, and persistence features.

Author: HalogenGroups v2.2.4
Date: February 2026
"""

import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

# Import HalogenGroups
from HalogenGroups import parse_smiles

def load_oecd_dataset(filepath='tests/OECDPFAS_list_22012019.csv', limit=None):
    """Load OECD PFAS dataset and extract valid SMILES."""
    print(f"\n{'='*80}")
    print("LOADING OECD PFAS DATASET")
    print('='*80)
    
    df = pd.read_csv(filepath, encoding='latin-1')
    print(f"Total entries in dataset: {len(df)}")
    
    # Filter for valid SMILES
    valid_smiles = df[df['SMILES'].notna() & (df['SMILES'] != '-')]['SMILES'].tolist()
    
    if limit:
        valid_smiles = valid_smiles[:limit]
    
    print(f"Valid SMILES entries: {len(valid_smiles)}")
    print(f"\nExample SMILES:")
    for i, smiles in enumerate(valid_smiles[:5], 1):
        print(f"  {i}. {smiles[:60]}{'...' if len(smiles) > 60 else ''}")
    
    return valid_smiles, df

def parse_and_convert(smiles_list, group_selection='all', count_mode='binary'):
    """Parse SMILES and convert to fingerprints."""
    print(f"\n{'='*80}")
    print("PARSING MOLECULES AND GENERATING FINGERPRINTS")
    print('='*80)
    
    print(f"Group selection: {group_selection}")
    print(f"Count mode: {count_mode}")
    print(f"Processing {len(smiles_list)} molecules...")
    
    # Pre-filter valid SMILES
    from rdkit import Chem
    valid_smiles = []
    for smi in smiles_list:
        try:
            mol = Chem.MolFromSmiles(smi)
            if mol is not None:
                valid_smiles.append(smi)
        except:
            pass
    
    print(f"Valid SMILES after filtering: {len(valid_smiles)}/{len(smiles_list)}")
    
    # Parse SMILES
    results = parse_smiles(valid_smiles)
    
    # Convert to fingerprints
    fp = results.to_fingerprint(
        group_selection=group_selection,
        count_mode=count_mode
    )
    
    print("\n" + fp.summary())
    
    return results, fp

def demonstrate_pca(fp, output_dir='oecd_analysis'):
    """Demonstrate PCA analysis."""
    print(f"\n{'='*80}")
    print("PRINCIPAL COMPONENT ANALYSIS (PCA)")
    print('='*80)
    
    # Perform PCA with 10 components
    pca_result = fp.perform_pca(
        n_components=10,
        plot=True,
        output_file=os.path.join(output_dir, 'oecd_pca.png')
    )
    
    # Analyze variance
    cumulative_var = np.cumsum(pca_result['explained_variance'])
    
    print(f"\nExplained variance by component:")
    for i, (var, cum_var) in enumerate(zip(pca_result['explained_variance'], cumulative_var), 1):
        print(f"  PC{i}: {var:.2%} (cumulative: {cum_var:.2%})")
    
    print(f"\nFirst 2 PCs explain: {cumulative_var[1]:.1%} of variance")
    print(f"First 5 PCs explain: {cumulative_var[4]:.1%} of variance")
    print(f"First 10 PCs explain: {cumulative_var[9]:.1%} of variance")
    
    # Create scree plot
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    
    # Variance explained
    axes[0].bar(range(1, 11), pca_result['explained_variance'], color='steelblue', alpha=0.8)
    axes[0].set_xlabel('Principal Component', fontsize=11)
    axes[0].set_ylabel('Explained Variance Ratio', fontsize=11)
    axes[0].set_title('Scree Plot - Variance per Component', fontsize=12, fontweight='bold')
    axes[0].grid(axis='y', alpha=0.3)
    
    # Cumulative variance
    axes[1].plot(range(1, 11), cumulative_var, marker='o', linewidth=2, 
                 markersize=8, color='darkgreen')
    axes[1].axhline(y=0.8, color='red', linestyle='--', label='80% threshold')
    axes[1].axhline(y=0.9, color='orange', linestyle='--', label='90% threshold')
    axes[1].set_xlabel('Number of Components', fontsize=11)
    axes[1].set_ylabel('Cumulative Explained Variance', fontsize=11)
    axes[1].set_title('Cumulative Variance Explained', fontsize=12, fontweight='bold')
    axes[1].legend()
    axes[1].grid(alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'oecd_pca_variance.png'), dpi=300, bbox_inches='tight')
    print(f"Variance plots saved to: {output_dir}/oecd_pca_variance.png")
    plt.close()
    
    return pca_result

def demonstrate_kernel_pca(fp, output_dir='oecd_analysis'):
    """Demonstrate kernel PCA with multiple kernels."""
    print(f"\n{'='*80}")
    print("KERNEL PCA ANALYSIS")
    print('='*80)
    
    kernels = ['rbf', 'poly', 'sigmoid', 'cosine']
    
    for kernel in kernels:
        print(f"\nKernel: {kernel}")
        kpca_result = fp.perform_kernel_pca(
            n_components=2,
            kernel=kernel,
            gamma=0.1 if kernel in ['rbf', 'poly', 'sigmoid'] else None,
            plot=True,
            output_file=os.path.join(output_dir, f'oecd_kpca_{kernel}.png')
        )
        print(f"  Transformed shape: {kpca_result['transformed'].shape}")
        print(f"  Plot saved to: {output_dir}/oecd_kpca_{kernel}.png")
    
    return kpca_result

def demonstrate_tsne(fp, output_dir='oecd_analysis'):
    """Demonstrate t-SNE with different perplexities."""
    print(f"\n{'='*80}")
    print("t-SNE ANALYSIS")
    print('='*80)
    
    perplexities = [10, 30, 50]
    
    for perp in perplexities:
        print(f"\nPerplexity: {perp}")
        tsne_result = fp.perform_tsne(
            n_components=2,
            perplexity=perp,
            learning_rate=200.0,
            max_iter=1000,
            plot=True,
            output_file=os.path.join(output_dir, f'oecd_tsne_perp{perp}.png')
        )
        print(f"  Transformed shape: {tsne_result['transformed'].shape}")
        print(f"  Plot saved to: {output_dir}/oecd_tsne_perp{perp}.png")
    
    return tsne_result

def demonstrate_umap(fp, output_dir='oecd_analysis'):
    """Demonstrate UMAP analysis."""
    print(f"\n{'='*80}")
    print("UMAP ANALYSIS")
    print('='*80)
    
    try:
        # Standard UMAP
        print("\nStandard UMAP (n_neighbors=15, min_dist=0.1)")
        umap_result = fp.perform_umap(
            n_components=2,
            n_neighbors=15,
            min_dist=0.1,
            plot=True,
            output_file=os.path.join(output_dir, 'oecd_umap_standard.png')
        )
        print(f"  Transformed shape: {umap_result['transformed'].shape}")
        print(f"  Plot saved to: {output_dir}/oecd_umap_standard.png")
        
        # Tight clusters
        print("\nTight clusters (n_neighbors=5, min_dist=0.01)")
        umap_tight = fp.perform_umap(
            n_components=2,
            n_neighbors=5,
            min_dist=0.01,
            plot=True,
            output_file=os.path.join(output_dir, 'oecd_umap_tight.png')
        )
        print(f"  Plot saved to: {output_dir}/oecd_umap_tight.png")
        
        # Global structure
        print("\nGlobal structure (n_neighbors=50, min_dist=0.5)")
        umap_global = fp.perform_umap(
            n_components=2,
            n_neighbors=50,
            min_dist=0.5,
            plot=True,
            output_file=os.path.join(output_dir, 'oecd_umap_global.png')
        )
        print(f"  Plot saved to: {output_dir}/oecd_umap_global.png")
        
        return umap_result
        
    except ImportError:
        print("\n⚠️  UMAP not available. Install with: pip install umap-learn")
        return None

def demonstrate_kl_divergence(smiles_list, output_dir='oecd_analysis'):
    """Demonstrate KL divergence comparison."""
    print(f"\n{'='*80}")
    print("KL DIVERGENCE COMPARISON")
    print('='*80)
    
    # Split dataset into two halves
    mid = len(smiles_list) // 2
    smiles_1 = smiles_list[:mid]
    smiles_2 = smiles_list[mid:]
    
    print(f"\nDataset 1: {len(smiles_1)} molecules")
    print(f"Dataset 2: {len(smiles_2)} molecules")
    
    # Pre-filter valid SMILES
    from rdkit import Chem
    valid_smiles_1 = [smi for smi in smiles_1 if Chem.MolFromSmiles(smi) is not None]
    valid_smiles_2 = [smi for smi in smiles_2 if Chem.MolFromSmiles(smi) is not None]
    
    print(f"Valid SMILES - Dataset 1: {len(valid_smiles_1)}/{len(smiles_1)}")
    print(f"Valid SMILES - Dataset 2: {len(valid_smiles_2)}/{len(smiles_2)}")
    
    # Parse and convert
    results_1 = parse_smiles(valid_smiles_1)
    results_2 = parse_smiles(valid_smiles_2)
    
    fp_1 = results_1.to_fingerprint(group_selection='all', count_mode='binary')
    fp_2 = results_2.to_fingerprint(group_selection='all', count_mode='binary')
    
    # Compare using different methods
    methods = ['minmax', 'forward', 'reverse', 'symmetric']
    
    print(f"\nKL Divergence Results:")
    print(f"{'Method':<15} {'KL Divergence':<15} {'Similarity':<15}")
    print('-' * 45)
    
    for method in methods:
        kl_div = fp_1.compare_kld(fp_2, method=method)
        similarity = 1 - kl_div if method == 'minmax' else None
        
        if similarity is not None:
            print(f"{method:<15} {kl_div:<15.4f} {similarity:<15.1%}")
        else:
            print(f"{method:<15} {kl_div:<15.4f} {'N/A':<15}")
    
    # Compare different group selections
    print(f"\n\nComparison across group selections:")
    print(f"{'Selection 1':<20} {'Selection 2':<20} {'KL (minmax)':<15} {'Similarity':<15}")
    print('-' * 70)
    
    selections = ['all', 'oecd', 'generic', 'generic+telomers']
    
    for i, sel1 in enumerate(selections):
        fp_s1 = results_1.to_fingerprint(group_selection=sel1)
        
        for sel2 in selections[i:]:
            fp_s2 = results_2.to_fingerprint(group_selection=sel2)
            
            try:
                kl_div = fp_s1.compare_kld(fp_s2, method='minmax')
                similarity = 1 - kl_div
                print(f"{sel1:<20} {sel2:<20} {kl_div:<15.4f} {similarity:<15.1%}")
            except ValueError as e:
                print(f"{sel1:<20} {sel2:<20} {'Error':<15} {str(e):<15}")
    
    return fp_1, fp_2

def demonstrate_clustering(fp, pca_result, output_dir='oecd_analysis'):
    """Demonstrate clustering on PCA space."""
    print(f"\n{'='*80}")
    print("CLUSTERING ON PCA SPACE")
    print('='*80)
    
    from sklearn.cluster import KMeans, DBSCAN
    from sklearn.metrics import silhouette_score
    
    # Use first 5 PCs
    X = pca_result['transformed'][:, :5]
    
    # K-means clustering
    print("\nK-Means Clustering:")
    best_k = None
    best_score = -1
    k_scores = []
    
    for k in range(2, 11):
        kmeans = KMeans(n_clusters=k, random_state=42, n_init=10)
        clusters = kmeans.fit_predict(X)
        score = silhouette_score(X, clusters)
        k_scores.append(score)
        
        print(f"  k={k}: Silhouette score = {score:.3f}")
        
        if score > best_score:
            best_score = score
            best_k = k
    
    print(f"\nBest k: {best_k} (silhouette score: {best_score:.3f})")
    
    # Perform final clustering with best k
    kmeans = KMeans(n_clusters=best_k, random_state=42, n_init=10)
    clusters = kmeans.fit_predict(X)
    
    # Plot clustering on first 2 PCs
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    
    # K-means results
    scatter = axes[0].scatter(
        pca_result['transformed'][:, 0],
        pca_result['transformed'][:, 1],
        c=clusters,
        cmap='tab10',
        alpha=0.6,
        s=30
    )
    axes[0].set_xlabel('PC1', fontsize=11)
    axes[0].set_ylabel('PC2', fontsize=11)
    axes[0].set_title(f'K-Means Clustering (k={best_k})', fontsize=12, fontweight='bold')
    axes[0].grid(alpha=0.3)
    plt.colorbar(scatter, ax=axes[0], label='Cluster')
    
    # Silhouette scores
    axes[1].plot(range(2, 11), k_scores, marker='o', linewidth=2, markersize=8, color='darkblue')
    axes[1].axvline(x=best_k, color='red', linestyle='--', label=f'Best k={best_k}')
    axes[1].set_xlabel('Number of Clusters (k)', fontsize=11)
    axes[1].set_ylabel('Silhouette Score', fontsize=11)
    axes[1].set_title('Clustering Quality', fontsize=12, fontweight='bold')
    axes[1].legend()
    axes[1].grid(alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'oecd_clustering.png'), dpi=300, bbox_inches='tight')
    print(f"\nClustering plot saved to: {output_dir}/oecd_clustering.png")
    plt.close()
    
    # Print cluster statistics
    print(f"\nCluster sizes:")
    unique, counts = np.unique(clusters, return_counts=True)
    for cluster_id, count in zip(unique, counts):
        print(f"  Cluster {cluster_id}: {count} molecules ({count/len(clusters)*100:.1f}%)")
    
    return clusters

def demonstrate_sql_persistence(results, fp, output_dir='oecd_analysis'):
    """Demonstrate SQL save/load functionality."""
    print(f"\n{'='*80}")
    print("SQL PERSISTENCE")
    print('='*80)
    
    # Create database files
    results_db = os.path.join(output_dir, 'oecd_results.db')
    fingerprints_db = os.path.join(output_dir, 'oecd_fingerprints.db')
    
    # Save results
    print(f"\nSaving ResultsModel to: {results_db}")
    results.to_sql(filename=results_db, if_exists='replace')
    print(f"  ✓ Results saved")
    
    # Save fingerprints
    print(f"\nSaving ResultsFingerprint to: {fingerprints_db}")
    fp.to_sql(filename=fingerprints_db, if_exists='replace')
    print(f"  ✓ Fingerprints saved")
    
    # Load back and verify
    print(f"\nLoading data back from databases...")
    from HalogenGroups.results_model import ResultsModel, ResultsFingerprint
    
    loaded_results = ResultsModel.from_sql(filename=results_db)
    loaded_fp = ResultsFingerprint.from_sql(filename=fingerprints_db)
    
    print(f"\nVerification:")
    print(f"  Original results: {len(results)} molecules")
    print(f"  Loaded results: {len(loaded_results)} molecules")
    print(f"  Original fingerprints: {fp.fingerprint_matrix.shape}")
    print(f"  Loaded fingerprints: {loaded_fp.fingerprint_matrix.shape}")
    print(f"  Match: {np.allclose(fp.fingerprint_matrix.toarray(), loaded_fp.fingerprint_matrix.toarray())}")
    
    # Get file sizes
    results_size = os.path.getsize(results_db) / 1024  # KB
    fp_size = os.path.getsize(fingerprints_db) / 1024  # KB
    
    print(f"\nDatabase sizes:")
    print(f"  Results: {results_size:.1f} KB")
    print(f"  Fingerprints: {fp_size:.1f} KB")
    
    return loaded_results, loaded_fp

def generate_summary_report(smiles_list, fp, pca_result, clusters, output_dir='oecd_analysis'):
    """Generate a comprehensive summary report."""
    print(f"\n{'='*80}")
    print("GENERATING SUMMARY REPORT")
    print('='*80)
    
    report_file = os.path.join(output_dir, 'oecd_analysis_report.txt')
    
    with open(report_file, 'w') as f:
        f.write("="*80 + "\n")
        f.write("OECD PFAS DATASET - ResultsFingerprint Analysis Report\n")
        f.write("="*80 + "\n\n")
        
        f.write("Dataset Information:\n")
        f.write("-" * 40 + "\n")
        f.write(f"Total molecules analyzed: {len(smiles_list)}\n")
        f.write(f"Group selection: {fp.group_selection}\n")
        f.write(f"Count mode: {fp.count_mode}\n")
        f.write(f"Number of groups: {fp.fingerprint_matrix.shape[1]}\n")
        f.write(f"Matrix sparsity: {1 - fp.fingerprint_matrix.nnz / (fp.fingerprint_matrix.shape[0] * fp.fingerprint_matrix.shape[1]):.1%}\n\n")
        
        f.write("PCA Analysis:\n")
        f.write("-" * 40 + "\n")
        cumulative_var = np.cumsum(pca_result['explained_variance'])
        f.write(f"Variance explained by PC1-PC2: {cumulative_var[1]:.1%}\n")
        f.write(f"Variance explained by PC1-PC5: {cumulative_var[4]:.1%}\n")
        f.write(f"Variance explained by PC1-PC10: {cumulative_var[9]:.1%}\n\n")
        
        f.write("Top Principal Components:\n")
        for i, var in enumerate(pca_result['explained_variance'][:5], 1):
            f.write(f"  PC{i}: {var:.2%}\n")
        f.write("\n")
        
        f.write("Clustering Results:\n")
        f.write("-" * 40 + "\n")
        unique, counts = np.unique(clusters, return_counts=True)
        f.write(f"Number of clusters: {len(unique)}\n")
        for cluster_id, count in zip(unique, counts):
            f.write(f"  Cluster {cluster_id}: {count} molecules ({count/len(clusters)*100:.1f}%)\n")
        f.write("\n")
        
        f.write("Files Generated:\n")
        f.write("-" * 40 + "\n")
        for file in os.listdir(output_dir):
            if file.endswith(('.png', '.db')):
                f.write(f"  {file}\n")
    
    print(f"\nReport saved to: {report_file}")
    
    # Display report
    with open(report_file, 'r') as f:
        print("\n" + f.read())

def main():
    """Main execution function."""
    print("\n" + "="*80)
    print(" "*20 + "OECD PFAS DATASET SHOWCASE")
    print(" "*15 + "ResultsFingerprint Capabilities Demo")
    print(" "*25 + "Version 2.2.4")
    print("="*80)
    
    # Create output directory
    output_dir = 'oecd_analysis'
    os.makedirs(output_dir, exist_ok=True)
    print(f"\nOutput directory: {output_dir}/")
    
    # 1. Load dataset (limit to 1000 for demonstration speed)
    smiles_list, df = load_oecd_dataset(limit=1000)
    
    # 2. Parse and convert to fingerprints
    results, fp = parse_and_convert(smiles_list, group_selection='all', count_mode='binary')
    
    # 3. Dimensionality reduction analyses
    pca_result = demonstrate_pca(fp, output_dir)
    demonstrate_kernel_pca(fp, output_dir)
    demonstrate_tsne(fp, output_dir)
    demonstrate_umap(fp, output_dir)
    
    # 4. KL divergence comparison
    fp_1, fp_2 = demonstrate_kl_divergence(smiles_list, output_dir)
    
    # 5. Clustering
    clusters = demonstrate_clustering(fp, pca_result, output_dir)
    
    # 6. SQL persistence
    demonstrate_sql_persistence(results, fp, output_dir)
    
    # 7. Generate report
    generate_summary_report(smiles_list, fp, pca_result, clusters, output_dir)
    
    print(f"\n{'='*80}")
    print("SHOWCASE COMPLETE!")
    print('='*80)
    print(f"\nAll results saved to: {output_dir}/")
    print("\nGenerated files:")
    for file in sorted(os.listdir(output_dir)):
        print(f"  • {file}")
    print("\n" + "="*80 + "\n")

if __name__ == '__main__':
    main()
