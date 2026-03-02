# ResultsFingerprint: Dimensionality Reduction and Analysis

## Overview

The `ResultsFingerprint` class provides powerful tools for analyzing PFAS group fingerprints through dimensionality reduction techniques and statistical comparisons. This extends the `ResultsModel` class with methods for:

- Converting results to fingerprint representations
- Performing PCA, kernel-PCA, t-SNE, and UMAP analyses
- Comparing fingerprint distributions using KL divergence
- Saving/loading fingerprints to/from SQL databases

## Quick Start

```python
from HalogenGroups import parse_smiles

# Parse some PFAS molecules
smiles_list = [
    "C(C(C(C(C(C(C(C(F)(F)F)(F)F)(F)F)(F)F)(F)F)(F)F)(F)F)(C(=O)O)F",  # PFOA
    "FC(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)S(=O)(=O)O",  # PFOS
]

results = parse_smiles(smiles_list)

# Default fingerprint: all halogens, per-saturation, 116 groups × 4 → shape (2, 464)
fp = results.to_fingerprint()

# Fluorine only → shape (2, 116)
fp_f = results.to_fingerprint(halogens='F')

# Stacked F + Cl fingerprint → shape (2, 232)
fp_fcl = results.to_fingerprint(halogens=['F', 'Cl'])

# PCA analysis
pca_results = fp.perform_pca(n_components=2, plot=True, output_file='pca.png')

# Compare with another dataset
results2 = parse_smiles(other_smiles_list)
fp2 = results2.to_fingerprint()
kl_div = fp.compare_kld(fp2, method='minmax')
print(f"KL divergence: {kl_div:.4f}")
```

## API Reference

### ResultsModel.to_fingerprint()

Convert a `ResultsModel` to a `ResultsFingerprint` for analysis.

**Parameters:**

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `group_selection` | `str` | `'all'` | Which groups to include (see table below) |
| `count_mode` | `str` | `'binary'` | Encoding method (`'binary'`, `'count'`, `'max_component'`) |
| `selected_group_ids` | `list` | `None` | Explicit group IDs — overrides `group_selection` |
| `halogens` | `str` or `list` | `['F','Cl','Br','I']` | Halogen(s) for component SMARTS matching (see below) |
| `saturation` | `str` or `None` | `'per'` | Saturation filter: `'per'`, `'poly'`, or `None` |

**Group selection options:**

| Value | Groups included | Count |
|-------|----------------|-------|
| `'all'` | All computable groups | 116 |
| `'oecd'` | OECD-defined groups (IDs 1–28) | 28 |
| `'generic'` | Generic functional groups (IDs 29–55) | 27 |
| `'telomers'` | Telomer-related groups (IDs 74–116) | 42 |
| `'generic+telomers'` | Generic + telomers combined | 69 |

**Halogen filtering:**

- **Single halogen** (e.g. `halogens='F'`): standard fingerprint of length *n_groups*.
- **Multiple halogens** (e.g. `halogens=['F', 'Cl']`): one fingerprint vector per halogen,
  **concatenated** into a single vector of length *n_groups × n_halogens*. Group names
  are automatically suffixed with `[F]`, `[Cl]`, etc.
- Available values: `'F'`, `'Cl'`, `'Br'`, `'I'`

**Saturation filtering** (applies only to groups with component SMARTS — OECD groups 1–28):

- `'per'`: perfluorinated / perhalogenated components only (default)
- `'poly'`: polyfluorinated / polyhalogenated components only
- `None`: no saturation filter

Groups without a component SMARTS (generic, telomer groups) are unaffected by the saturation filter.

**Returns:** `ResultsFingerprint`

**Examples:**

```python
# Default: all halogens, per-saturation, 116 groups × 4 → shape (n, 464)
fp = results.to_fingerprint()

# Fluorine only → shape (n, 116)
fp = results.to_fingerprint(halogens='F')

# Stacked F + Cl fingerprint → shape (n, 232), names suffixed [F] / [Cl]
fp = results.to_fingerprint(halogens=['F', 'Cl'])

# F + Cl + Br → shape (n, 348)
fp = results.to_fingerprint(halogens=['F', 'Cl', 'Br'])

# Polyfluorinated components only
fp = results.to_fingerprint(halogens='F', saturation='poly')

# No saturation filter
fp = results.to_fingerprint(halogens='F', saturation=None)

# Count encoding with OECD groups only
fp_oecd = results.to_fingerprint(group_selection='oecd', count_mode='count')

# Custom group selection
fp_custom = results.to_fingerprint(selected_group_ids=[1, 2, 5, 10, 15])
```

### ResultsFingerprint.perform_pca()

Perform Principal Component Analysis (PCA) on fingerprints.

**Parameters:**
- `n_components` (int): Number of principal components (default: 2)
- `plot` (bool): Whether to create visualization (default: True)
- `output_file` (str): Path to save plot (None for interactive display)

**Returns:** Dictionary containing:
- `'transformed'`: PCA-transformed data
- `'explained_variance'`: Explained variance ratio for each component
- `'components'`: Principal component vectors
- `'pca_model'`: Fitted PCA model
- `'scaler'`: StandardScaler used for preprocessing

**Example:**
```python
# Basic PCA with 2 components
pca_results = fp.perform_pca(n_components=2, plot=True)

# PCA with more components
pca_results = fp.perform_pca(n_components=10, plot=False)
print(f"Variance explained: {np.cumsum(pca_results['explained_variance'])}")

# Save plot to file
pca_results = fp.perform_pca(output_file='pfas_pca.png')
```

**Scientific Background:**

PCA is a linear dimensionality reduction technique that projects data onto orthogonal axes (principal components) that maximize variance. For PFAS fingerprints:
- PC1 typically captures the most dominant structural patterns
- Explained variance indicates how much information each component retains
- Useful for identifying major sources of structural variation

### ResultsFingerprint.perform_kernel_pca()

Perform kernel PCA for non-linear dimensionality reduction.

**Parameters:**
- `n_components` (int): Number of components (default: 2)
- `kernel` (str): Kernel type (default: 'rbf')
  - `'linear'`: Linear kernel (equivalent to PCA)
  - `'poly'`: Polynomial kernel
  - `'rbf'`: Radial basis function (Gaussian) kernel
  - `'sigmoid'`: Sigmoid kernel
  - `'cosine'`: Cosine kernel
- `gamma` (float): Kernel coefficient (default: 1/n_features)
- `plot` (bool): Whether to create visualization (default: True)
- `output_file` (str): Path to save plot

**Returns:** Dictionary containing transformed data and model

**Example:**
```python
# RBF kernel (default)
kpca_rbf = fp.perform_kernel_pca(kernel='rbf', gamma=0.1)

# Polynomial kernel
kpca_poly = fp.perform_kernel_pca(kernel='poly', gamma=0.05)

# Compare different kernels
for kernel in ['rbf', 'poly', 'sigmoid']:
    kpca = fp.perform_kernel_pca(kernel=kernel, output_file=f'kpca_{kernel}.png')
```

**Scientific Background:**

Kernel PCA extends PCA to capture non-linear relationships by implicitly mapping data to a higher-dimensional space. For PFAS fingerprints:
- RBF kernel is effective for complex non-linear patterns
- Can reveal structural relationships not visible in linear PCA
- Gamma parameter controls the "reach" of each data point

### ResultsFingerprint.perform_tsne()

Perform t-Distributed Stochastic Neighbor Embedding (t-SNE).

**Parameters:**
- `n_components` (int): Number of dimensions (default: 2)
- `perplexity` (float): Perplexity parameter (default: 30.0)
  - Controls balance between local and global structure
  - Typical range: 5-50
  - Lower values focus on local structure
  - Higher values preserve global structure
- `learning_rate` (float): Learning rate (default: 200.0)
- `max_iter` (int): Maximum number of iterations (default: 1000)
- `plot` (bool): Whether to create visualization (default: True)
- `output_file` (str): Path to save plot

**Returns:** Dictionary containing transformed data and model

**Example:**
```python
# Basic t-SNE
tsne = fp.perform_tsne(perplexity=30, max_iter=1000)

# Try different perplexities
for perp in [5, 15, 30, 50]:
    tsne = fp.perform_tsne(perplexity=perp, output_file=f'tsne_perp{perp}.png')

# Higher iterations for convergence
tsne = fp.perform_tsne(max_iter=2000, learning_rate=300)
```

**Scientific Background:**

t-SNE is a non-linear technique that excels at visualizing high-dimensional data by preserving local neighborhood structure. For PFAS fingerprints:
- Excellent for revealing clusters and patterns
- Perplexity is the most important parameter
- Not deterministic (use random_state for reproducibility)
- Distances in t-SNE plots should not be over-interpreted

### ResultsFingerprint.perform_umap()

Perform Uniform Manifold Approximation and Projection (UMAP).

**Parameters:**
- `n_components` (int): Number of dimensions (default: 2)
- `n_neighbors` (int): Number of neighbors (default: 15)
  - Controls local vs global structure
  - Lower values focus on local structure
  - Higher values preserve global structure
- `min_dist` (float): Minimum distance between points (default: 0.1)
  - Controls clustering tightness
  - Lower values create tighter clusters
  - Higher values create more spread out embeddings
- `metric` (str): Distance metric (default: 'euclidean')
- `plot` (bool): Whether to create visualization (default: True)
- `output_file` (str): Path to save plot

**Returns:** Dictionary containing transformed data and model

**Example:**
```python
# Basic UMAP
umap_results = fp.perform_umap(n_neighbors=15, min_dist=0.1)

# Tighter clusters
umap_tight = fp.perform_umap(n_neighbors=5, min_dist=0.01)

# Global structure preservation
umap_global = fp.perform_umap(n_neighbors=50, min_dist=0.5)

# Custom distance metric
umap_cosine = fp.perform_umap(metric='cosine')
```

**Scientific Background:**

UMAP is a modern dimensionality reduction technique that often performs better than t-SNE. For PFAS fingerprints:
- Faster than t-SNE
- Better preserves global structure
- More deterministic (given same parameters)
- n_neighbors is analogous to t-SNE's perplexity

**Requires:** `pip install umap-learn`

### ResultsFingerprint.compare_kld()

Compare two fingerprint sets using Kullback-Leibler (KL) divergence.

**Parameters:**
- `other` (ResultsFingerprint): Other fingerprint set to compare
- `method` (str): Comparison method (default: 'minmax')
  - `'minmax'`: Min-max normalized symmetric KL divergence (recommended)
  - `'forward'`: KL(self || other)
  - `'reverse'`: KL(other || self)
  - `'symmetric'`: Average of forward and reverse

**Returns:** float - KL divergence value (lower = more similar)

**Example:**
```python
# Parse two datasets
results1 = parse_smiles(smiles_list_1)
results2 = parse_smiles(smiles_list_2)

fp1 = results1.to_fingerprint(group_selection='all')
fp2 = results2.to_fingerprint(group_selection='all')

# Compare using minmax method
kl_div = fp1.compare_kld(fp2, method='minmax')
print(f"KL divergence (minmax): {kl_div:.4f}")

# Compare using all methods
for method in ['minmax', 'forward', 'reverse', 'symmetric']:
    kl = fp1.compare_kld(fp2, method=method)
    print(f"{method:10s}: {kl:.6f}")

# Interpretation
if kl_div < 0.1:
    print("Very similar distributions")
elif kl_div < 0.3:
    print("Moderately similar distributions")
else:
    print("Different distributions")
```

**Scientific Background:**

KL divergence quantifies how one probability distribution differs from another. For PFAS fingerprints:
- Measures compositional similarity between datasets
- 0 means identical distributions
- Higher values indicate greater difference
- 'minmax' method normalizes to [0,1] range for interpretability
- Useful for comparing chemical inventories, databases, or subsets

### ResultsFingerprint.to_sql() / from_sql()

Save and load fingerprints to/from SQL databases.

**to_sql() Parameters:**
- `conn` (str or Engine): Database connection string or SQLAlchemy engine
- `filename` (str): SQLite database filename (alternative to conn)
- `table_name` (str): Name for fingerprint table (default: 'fingerprints')
- `metadata_table` (str): Name for metadata table (default: 'fingerprint_metadata')
- `if_exists` (str): Behavior if table exists: 'fail', 'replace', 'append' (default: 'append')

**from_sql() Parameters:**
- `conn` (str or Engine): Database connection string or SQLAlchemy engine
- `filename` (str): SQLite database filename
- `table_name` (str): Name of fingerprint table (default: 'fingerprints')
- `metadata_table` (str): Name of metadata table (default: 'fingerprint_metadata')
- `limit` (int): Limit number of molecules to load

**Example:**
```python
# Save to SQLite
fp.to_sql(filename='pfas_fingerprints.db', if_exists='replace')

# Save to PostgreSQL
fp.to_sql(conn='postgresql://user:pass@localhost/dbname')

# Load from database
from HalogenGroups import ResultsFingerprint
fp_loaded = ResultsFingerprint.from_sql(filename='pfas_fingerprints.db')

# Load subset
fp_subset = ResultsFingerprint.from_sql(filename='pfas_fingerprints.db', limit=100)
```

### ResultsModel.from_sql()

Load results from SQL database.

**Parameters:**
- `conn` (str or Engine): Database connection string or SQLAlchemy engine
- `filename` (str): SQLite database filename
- `components_table` (str): Name of components table (default: 'components')
- `groups_table` (str): Name of groups table (default: 'pfas_groups_in_compound')
- `limit` (int): Limit number of molecules to load

**Example:**
```python
from HalogenGroups import ResultsModel

# Save results
results.to_sql(filename='results.db', if_exists='replace')

# Load results
results_loaded = ResultsModel.from_sql(filename='results.db')

# Load subset
results_subset = ResultsModel.from_sql(filename='results.db', limit=50)
```

## Complete Workflow Example

```python
from HalogenGroups import parse_smiles
import numpy as np
import matplotlib.pyplot as plt

# 1. Parse molecules
smiles_list = [...]  # Your SMILES
results = parse_smiles(smiles_list)

# 2. Convert to fingerprints
fp = results.to_fingerprint(group_selection='all', count_mode='binary')
print(fp.summary())

# 3. Perform multiple dimensionality reductions
pca = fp.perform_pca(n_components=5, plot=True, output_file='pca.png')
tsne = fp.perform_tsne(perplexity=30, plot=True, output_file='tsne.png')
umap = fp.perform_umap(n_neighbors=15, plot=True, output_file='umap.png')

# 4. Analyze explained variance
print(f"PCA variance explained: {np.cumsum(pca['explained_variance'])}")

# 5. Compare with reference dataset
ref_results = parse_smiles(reference_smiles)
ref_fp = ref_results.to_fingerprint(group_selection='all', count_mode='binary')
kl_div = fp.compare_kld(ref_fp, method='minmax')
print(f"KL divergence vs reference: {kl_div:.4f}")

# 6. Save for later use
fp.to_sql(filename='fingerprints.db', if_exists='replace')
results.to_sql(filename='results.db', if_exists='replace')

# 7. Cluster analysis on PCA space
from sklearn.cluster import KMeans
kmeans = KMeans(n_clusters=5, random_state=42)
clusters = kmeans.fit_predict(pca['transformed'])

plt.figure(figsize=(8, 6))
scatter = plt.scatter(pca['transformed'][:, 0], pca['transformed'][:, 1], 
                     c=clusters, cmap='viridis', alpha=0.6)
plt.xlabel(f"PC1 ({pca['explained_variance'][0]:.1%})")
plt.ylabel(f"PC2 ({pca['explained_variance'][1]:.1%})")
plt.title('PFAS Clusters in PCA Space')
plt.colorbar(scatter, label='Cluster')
plt.tight_layout()
plt.savefig('pca_clusters.png', dpi=300)
plt.show()
```

## Dependencies

Required:
- numpy
- pandas
- scipy
- scikit-learn
- matplotlib
- sqlalchemy (for SQL operations)

Optional:
- umap-learn (for UMAP analysis)

Install all dependencies:
```bash
pip install numpy pandas scipy scikit-learn matplotlib sqlalchemy umap-learn
```

## Tips and Best Practices

### Choosing Group Selection
- **'all'**: Use for comprehensive analysis of all PFAS features
- **'oecd'**: Use for regulatory-focused analysis
- **'generic'**: Use for general functional group patterns
- **'generic+telomers'**: Use for transformation product analysis

### Choosing Count Mode
- **'binary'**: Best for presence/absence analysis and most ML applications
- **'count'**: Use when multiple matches matter (e.g., polymers)
- **'max_component'**: Use when chain length is important

### Dimensionality Reduction Selection
- **PCA**: Fast, interpretable, linear relationships
- **Kernel PCA**: Non-linear patterns, still relatively fast
- **t-SNE**: Best for visualization, slow for large datasets
- **UMAP**: Best balance of speed and quality, preserves global structure

### Parameter Tuning

**t-SNE perplexity:**
- Small datasets (< 100): perplexity = 5-15
- Medium datasets (100-1000): perplexity = 20-50
- Large datasets (> 1000): perplexity = 30-100

**UMAP n_neighbors:**
- Local structure: n_neighbors = 5-15
- Balanced: n_neighbors = 15-30
- Global structure: n_neighbors = 30-100

**KL Divergence interpretation:**
- < 0.1: Very similar distributions
- 0.1-0.3: Moderately similar
- 0.3-0.5: Different but related
- > 0.5: Very different distributions

## References

1. van der Maaten & Hinton (2008). "Visualizing Data using t-SNE". JMLR 9:2579-2605.
2. McInnes et al. (2018). "UMAP: Uniform Manifold Approximation and Projection". arXiv:1802.03426.
3. Schölkopf et al. (1998). "Nonlinear Component Analysis as a Kernel Eigenvalue Problem". Neural Computation.
4. Kullback & Leibler (1951). "On Information and Sufficiency". Annals of Mathematical Statistics.

## Troubleshooting

**Issue:** t-SNE takes too long
- Solution: Reduce dataset size or use UMAP instead

**Issue:** UMAP not available
- Solution: Install with `pip install umap-learn`

**Issue:** Memory error with large datasets
- Solution: Use PCA first to reduce dimensions, then apply t-SNE/UMAP

**Issue:** KL divergence is NaN
- Solution: Ensure both fingerprints have same group selection and non-empty data

## License

This functionality is part of PFASgroups and follows the same license.
