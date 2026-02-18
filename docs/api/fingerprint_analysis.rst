ResultsFingerprint Analysis
===========================

The ``ResultsFingerprint`` class provides advanced analysis capabilities for PFAS group fingerprints, including dimensionality reduction, statistical comparison, and database persistence.

Overview
--------

``ResultsFingerprint`` enables sophisticated analysis of PFAS classification results through:

- **Dimensionality Reduction**: Visualize and explore high-dimensional fingerprint data
- **Statistical Comparison**: Quantify similarity between chemical inventories
- **Flexible Encoding**: Choose appropriate fingerprint representation for your use case
- **Database Integration**: Efficient storage and retrieval of analysis results

Quick Start
-----------

.. code-block:: python

   from PFASgroups import parse_smiles
   
   # Parse molecules
   results = parse_smiles(smiles_list)
   
   # Convert to fingerprints
   fp = results.to_fingerprint(group_selection='oecd', count_mode='binary')
   
   # Perform analysis
   pca = fp.perform_pca(n_components=5, plot=True)
   tsne = fp.perform_tsne(perplexity=30, plot=True)
   
   # Compare datasets
   other_fp = other_results.to_fingerprint(group_selection='oecd')
   similarity = fp.compare_kld(other_fp, method='minmax')

ResultsModel Methods
--------------------

to_fingerprint()
^^^^^^^^^^^^^^^^

.. py:method:: ResultsModel.to_fingerprint(group_selection='all', count_mode='binary', selected_group_ids=None, halogens='F', saturation='per')

   Convert ResultsModel to ResultsFingerprint for analysis.
   
   :param group_selection: Group selection mode
   :type group_selection: str
   :param count_mode: Encoding mode (``'binary'``, ``'count'``, ``'max_component'``)
   :type count_mode: str
   :param selected_group_ids: Custom list of group IDs (overrides group_selection)
   :type selected_group_ids: list, optional
   :param halogens: Halogen(s) to use for component SMARTS matching.
       Single value → 116-column fingerprint.
       List of values → vectors stacked per halogen (length 116 × n_halogens),
       with group names suffixed ``[F]``, ``[Cl]``, etc.
   :type halogens: str or list of str
   :param saturation: Saturation filter for component SMARTS groups.
       ``'per'`` (perfluorinated/perhalogenated), ``'poly'`` (polyfluorinated/
       polyhalogenated), or ``None`` (no filter). Groups without component SMARTS
       are unaffected.
   :type saturation: str or None
   :returns: ResultsFingerprint object
   :rtype: ResultsFingerprint
   
   **Group Selection Options:**
   
   - ``'all'`` (default): All 116 computable groups
   - ``'oecd'``: OECD-defined groups (IDs 1–28, 28 groups)
   - ``'generic'``: Generic functional groups (IDs 29–55, 27 groups)
   - ``'telomers'``: Telomer-related groups (IDs 74–116, 42 groups)
   - ``'generic+telomers'``: Combined selection (69 groups)
   
   **Count Modes:**
   
   - ``'binary'``: 1 if present, 0 if absent (recommended for ML)
   - ``'count'``: Number of matches (for multiple occurrences)
   - ``'max_component'``: Maximum component size (chain length)
   
   **Example:**
   
   .. code-block:: python
   
      # Default: F only, per-saturation, 116 groups → shape (n, 116)
      fp = results.to_fingerprint()
      
      # Stacked F + Cl fingerprint → shape (n, 232)
      fp = results.to_fingerprint(halogens=['F', 'Cl'])
      
      # Polyfluorinated components only
      fp = results.to_fingerprint(halogens='F', saturation='poly')
      
      # OECD groups with count encoding
      fp = results.to_fingerprint(group_selection='oecd', count_mode='count')
      
      # Custom groups
      fp = results.to_fingerprint(selected_group_ids=[1, 2, 5, 10, 15])

from_sql()
^^^^^^^^^^

.. py:method:: ResultsModel.from_sql(conn=None, filename=None, components_table='components', groups_table='pfas_groups_in_compound', limit=None)

   Load results from SQL database.
   
   :param conn: Database connection string or SQLAlchemy engine
   :type conn: str or Engine, optional
   :param filename: SQLite database filename
   :type filename: str, optional
   :param components_table: Name of components table
   :type components_table: str
   :param groups_table: Name of groups table
   :type groups_table: str
   :param limit: Maximum number of molecules to load
   :type limit: int, optional
   :returns: Loaded results
   :rtype: ResultsModel
   
   **Example:**
   
   .. code-block:: python
   
      # Load from SQLite
      results = ResultsModel.from_sql(filename='results.db')
      
      # Load from PostgreSQL
      results = ResultsModel.from_sql(conn='postgresql://user:pass@host/db')

ResultsFingerprint Methods
--------------------------

Dimensionality Reduction
^^^^^^^^^^^^^^^^^^^^^^^^^

perform_pca()
"""""""""""""

.. py:method:: ResultsFingerprint.perform_pca(n_components=2, plot=True, output_file=None)

   Perform Principal Component Analysis.
   
   :param n_components: Number of principal components
   :type n_components: int
   :param plot: Whether to create visualization
   :type plot: bool
   :param output_file: Path to save plot (None for interactive)
   :type output_file: str, optional
   :returns: Dictionary with transformed data, explained variance, components, PCA model, and scaler
   :rtype: dict
   
   **Returns:**
   
   - ``transformed``: PCA-transformed data array
   - ``explained_variance``: Explained variance ratio per component
   - ``components``: Principal component vectors
   - ``pca_model``: Fitted scikit-learn PCA model
   - ``scaler``: StandardScaler used for preprocessing
   
   **Scientific Background:**
   
   PCA is a linear dimensionality reduction technique that projects data onto orthogonal axes 
   maximizing variance. For PFAS fingerprints:
   
   - PC1 captures the most dominant structural patterns
   - Explained variance indicates information retention
   - Useful for identifying major sources of variation
   
   **Example:**
   
   .. code-block:: python
   
      # Basic PCA
      pca = fp.perform_pca(n_components=5, plot=True)
      print(f"Variance explained: {np.cumsum(pca['explained_variance'])}")
      
      # Save plot to file
      pca = fp.perform_pca(output_file='pca_analysis.png')

perform_kernel_pca()
""""""""""""""""""""

.. py:method:: ResultsFingerprint.perform_kernel_pca(n_components=2, kernel='rbf', gamma=None, plot=True, output_file=None)

   Perform kernel PCA for non-linear dimensionality reduction.
   
   :param n_components: Number of components
   :type n_components: int
   :param kernel: Kernel type ('linear', 'poly', 'rbf', 'sigmoid', 'cosine')
   :type kernel: str
   :param gamma: Kernel coefficient (default: 1/n_features)
   :type gamma: float, optional
   :param plot: Whether to create visualization
   :type plot: bool
   :param output_file: Path to save plot
   :type output_file: str, optional
   :returns: Dictionary with transformed data and model
   :rtype: dict
   
   **Kernel Options:**
   
   - ``'rbf'``: Radial basis function (Gaussian) - default, effective for non-linear patterns
   - ``'poly'``: Polynomial kernel
   - ``'sigmoid'``: Sigmoid kernel
   - ``'cosine'``: Cosine kernel
   - ``'linear'``: Linear kernel (equivalent to PCA)
   
   **Example:**
   
   .. code-block:: python
   
      # RBF kernel (default)
      kpca = fp.perform_kernel_pca(kernel='rbf', gamma=0.1)
      
      # Compare different kernels
      for kernel in ['rbf', 'poly', 'sigmoid']:
          kpca = fp.perform_kernel_pca(
              kernel=kernel, 
              output_file=f'kpca_{kernel}.png'
          )

perform_tsne()
""""""""""""""

.. py:method:: ResultsFingerprint.perform_tsne(n_components=2, perplexity=30.0, learning_rate=200.0, n_iter=1000, plot=True, output_file=None)

   Perform t-Distributed Stochastic Neighbor Embedding.
   
   :param n_components: Number of dimensions
   :type n_components: int
   :param perplexity: Perplexity parameter (5-50 typical)
   :type perplexity: float
   :param learning_rate: Learning rate for optimization
   :type learning_rate: float
   :param n_iter: Number of iterations
   :type n_iter: int
   :param plot: Whether to create visualization
   :type plot: bool
   :param output_file: Path to save plot
   :type output_file: str, optional
   :returns: Dictionary with transformed data and model
   :rtype: dict
   
   **Perplexity Guidelines:**
   
   - Small datasets (<100): 5-15
   - Medium datasets (100-1000): 20-50
   - Large datasets (>1000): 30-100
   - Lower values emphasize local structure
   - Higher values preserve global structure
   
   **Scientific Background:**
   
   t-SNE excels at visualizing high-dimensional data by preserving local neighborhood structure.
   Excellent for revealing clusters and patterns, though distances should not be over-interpreted.
   
   **Example:**
   
   .. code-block:: python
   
      # Standard visualization
      tsne = fp.perform_tsne(perplexity=30, n_iter=1000)
      
      # Try different perplexities
      for perp in [5, 15, 30, 50]:
          tsne = fp.perform_tsne(
              perplexity=perp, 
              output_file=f'tsne_perp{perp}.png'
          )

perform_umap()
""""""""""""""

.. py:method:: ResultsFingerprint.perform_umap(n_components=2, n_neighbors=15, min_dist=0.1, metric='euclidean', plot=True, output_file=None)

   Perform Uniform Manifold Approximation and Projection.
   
   :param n_components: Number of dimensions
   :type n_components: int
   :param n_neighbors: Number of neighbors (5-100)
   :type n_neighbors: int
   :param min_dist: Minimum distance between points (0.0-1.0)
   :type min_dist: float
   :param metric: Distance metric ('euclidean', 'cosine', etc.)
   :type metric: str
   :param plot: Whether to create visualization
   :type plot: bool
   :param output_file: Path to save plot
   :type output_file: str, optional
   :returns: Dictionary with transformed data and model
   :rtype: dict
   
   **Requires:** ``pip install umap-learn``
   
   **Parameter Guidelines:**
   
   - **n_neighbors**: Controls local vs global structure
     
     - 5-15: Emphasize local structure
     - 15-30: Balanced (default)
     - 30-100: Preserve global structure
   
   - **min_dist**: Controls clustering tightness
     
     - 0.0-0.1: Tight clusters
     - 0.1-0.5: Moderate spread
     - 0.5-1.0: Loose, spread out
   
   **Advantages:**
   
   - Faster than t-SNE
   - Better preserves global structure
   - More deterministic (given same parameters)
   - Scales to larger datasets
   
   **Example:**
   
   .. code-block:: python
   
      # Standard UMAP
      umap = fp.perform_umap(n_neighbors=15, min_dist=0.1)
      
      # Tight clusters
      umap = fp.perform_umap(n_neighbors=5, min_dist=0.01)
      
      # Global structure
      umap = fp.perform_umap(n_neighbors=50, min_dist=0.5)

Statistical Comparison
^^^^^^^^^^^^^^^^^^^^^^

compare_kld()
"""""""""""""

.. py:method:: ResultsFingerprint.compare_kld(other, method='minmax')

   Compare two fingerprint sets using Kullback-Leibler divergence.
   
   :param other: Other fingerprint set to compare
   :type other: ResultsFingerprint
   :param method: Comparison method ('minmax', 'forward', 'reverse', 'symmetric')
   :type method: str
   :returns: KL divergence value (lower = more similar)
   :rtype: float
   
   **Methods:**
   
   - ``'minmax'`` (recommended): Normalized symmetric KL divergence (0-1 range)
   - ``'forward'``: KL(self || other)
   - ``'reverse'``: KL(other || self)
   - ``'symmetric'``: Average of forward and reverse
   
   **Interpretation:**
   
   ========  =========================
   Value     Interpretation
   ========  =========================
   < 0.1     Very similar distributions
   0.1-0.3   Moderately similar
   0.3-0.5   Different but related
   > 0.5     Very different
   ========  =========================
   
   **Scientific Background:**
   
   KL divergence quantifies how one probability distribution differs from another.
   For PFAS fingerprints, it measures compositional similarity between datasets,
   useful for comparing chemical inventories or databases.
   
   **Example:**
   
   .. code-block:: python
   
      # Compare two datasets
      fp1 = results1.to_fingerprint(group_selection='all')
      fp2 = results2.to_fingerprint(group_selection='all')
      
      kl_div = fp1.compare_kld(fp2, method='minmax')
      print(f"Similarity: {1 - kl_div:.1%}")
      
      # Try all methods
      for method in ['minmax', 'forward', 'reverse', 'symmetric']:
          kl = fp1.compare_kld(fp2, method=method)
          print(f"{method}: {kl:.4f}")

Database Operations
^^^^^^^^^^^^^^^^^^^

to_sql()
""""""""

.. py:method:: ResultsFingerprint.to_sql(conn=None, filename=None, table_name='fingerprints', metadata_table='fingerprint_metadata', if_exists='append')

   Save fingerprints to SQL database.
   
   :param conn: Database connection string or engine
   :type conn: str or Engine, optional
   :param filename: SQLite database filename
   :type filename: str, optional
   :param table_name: Name for fingerprint table
   :type table_name: str
   :param metadata_table: Name for metadata table
   :type metadata_table: str
   :param if_exists: Behavior if table exists ('fail', 'replace', 'append')
   :type if_exists: str
   
   **Example:**
   
   .. code-block:: python
   
      # Save to SQLite
      fp.to_sql(filename='fingerprints.db', if_exists='replace')
      
      # Save to PostgreSQL
      fp.to_sql(conn='postgresql://user:pass@localhost/dbname')

from_sql()
""""""""""

.. py:classmethod:: ResultsFingerprint.from_sql(conn=None, filename=None, table_name='fingerprints', metadata_table='fingerprint_metadata', limit=None)

   Load fingerprints from SQL database.
   
   :param conn: Database connection string or engine
   :type conn: str or Engine, optional
   :param filename: SQLite database filename
   :type filename: str, optional
   :param table_name: Name of fingerprints table
   :type table_name: str
   :param metadata_table: Name of metadata table
   :type metadata_table: str
   :param limit: Maximum number of molecules to load
   :type limit: int, optional
   :returns: Loaded fingerprints
   :rtype: ResultsFingerprint
   
   **Example:**
   
   .. code-block:: python
   
      # Load from SQLite
      fp = ResultsFingerprint.from_sql(filename='fingerprints.db')
      
      # Load subset
      fp = ResultsFingerprint.from_sql(filename='fingerprints.db', limit=100)

Utility Methods
^^^^^^^^^^^^^^^

summary()
"""""""""

.. py:method:: ResultsFingerprint.summary()

   Return text summary of fingerprint statistics.
   
   :returns: Summary string with statistics
   :rtype: str
   
   Includes:
   
   - Number of molecules and groups
   - Group selection and count mode
   - Fingerprint dimensions and sparsity
   - Most common PFAS groups

Complete Workflow Example
--------------------------

.. code-block:: python

   from PFASgroups import parse_smiles
   import numpy as np
   
   # 1. Parse molecules
   smiles_list = [...]  # Your SMILES
   results = parse_smiles(smiles_list)
   
   # 2. Convert to fingerprints
   fp = results.to_fingerprint(group_selection='all', count_mode='binary')
   print(fp.summary())
   
   # 3. Dimensionality reduction
   pca = fp.perform_pca(n_components=10, plot=True, output_file='pca.png')
   tsne = fp.perform_tsne(perplexity=30, plot=True, output_file='tsne.png')
   umap = fp.perform_umap(n_neighbors=15, plot=True, output_file='umap.png')
   
   # 4. Analyze variance
   cumulative_var = np.cumsum(pca['explained_variance'])
   print(f"First 5 PCs explain: {cumulative_var[4]:.1%} of variance")
   
   # 5. Compare with reference
   ref_results = parse_smiles(reference_smiles)
   ref_fp = ref_results.to_fingerprint(group_selection='all')
   kl_div = fp.compare_kld(ref_fp, method='minmax')
   print(f"Similarity to reference: {1 - kl_div:.1%}")
   
   # 6. Clustering on PCA space
   from sklearn.cluster import KMeans
   kmeans = KMeans(n_clusters=5, random_state=42)
   clusters = kmeans.fit_predict(pca['transformed'])
   
   # 7. Save for later
   fp.to_sql(filename='fingerprints.db', if_exists='replace')
   results.to_sql(filename='results.db', if_exists='replace')

Dependencies
------------

**Required:**

- numpy
- pandas  
- scipy
- scikit-learn
- matplotlib
- sqlalchemy

**Optional:**

- umap-learn (for UMAP analysis)

Install all dependencies:

.. code-block:: bash

   pip install numpy pandas scipy scikit-learn matplotlib sqlalchemy umap-learn

See Also
--------

- :doc:`quickstart` - Getting started guide
- :doc:`tutorial` - Detailed tutorials
- :doc:`api/core` - Core API documentation
- ``examples/results_fingerprint_analysis.py`` - Working examples
- ``docs/ResultsFingerprint_Guide.md`` - Complete guide

References
----------

1. van der Maaten, L., & Hinton, G. (2008). Visualizing Data using t-SNE. *Journal of Machine Learning Research*, 9:2579-2605.

2. McInnes, L., Healy, J., & Melville, J. (2018). UMAP: Uniform Manifold Approximation and Projection for Dimension Reduction. arXiv:1802.03426.

3. Schölkopf, B., Smola, A., & Müller, K. R. (1998). Nonlinear Component Analysis as a Kernel Eigenvalue Problem. *Neural Computation*, 10(5):1299-1319.

4. Kullback, S., & Leibler, R. A. (1951). On Information and Sufficiency. *The Annals of Mathematical Statistics*, 22(1):79-86.
