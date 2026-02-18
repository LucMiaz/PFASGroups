Changelog
=========

Version 3.1.0 (February 2026)
------------------------------

**Released:** February 18, 2026

**New Features:**

- **Multi-halogen fingerprinting**: ``generate_fingerprint`` and ``ResultsModel.to_fingerprint()``
  now accept a ``halogens`` parameter (``'F'``, ``'Cl'``, ``'Br'``, ``'I'``, or a list thereof).

  - Single halogen (default ``'F'``): standard 116-column binary/count vector.
  - Multiple halogens (e.g. ``['F', 'Cl']``): vectors are **stacked** horizontally,
    yielding a fingerprint of length 116 × n_halogens. Group names are automatically
    suffixed with ``[F]``, ``[Cl]``, etc.

- **Saturation filter**: a ``saturation`` parameter (``'per'`` | ``'poly'`` | ``None``,
  default ``'per'``) controls which component SMARTS are used for groups that have
  halogenated-chain components (OECD groups 1–28). Groups without a component SMARTS
  (generic/telomer groups) are unaffected by this filter.

- **Corrected group-selection indexing**: group selections (``'oecd'``, ``'generic'``,
  ``'telomers'``, ``'generic+telomers'``) now resolve using canonical **group IDs**
  rather than raw list positions, making them robust to future changes in the data file.

- **``ResultsFingerprint`` metadata**: the class now stores and displays ``halogens``
  and ``saturation`` in ``__repr__`` and ``summary()``.

**API Changes:**

.. code-block:: python

   # Default: F only, perfluorinated, all 116 groups → shape (n_mols, 116)
   fp = results.to_fingerprint()

   # Fluorine + chlorine stacked → shape (n_mols, 232)
   fp = results.to_fingerprint(halogens=['F', 'Cl'])

   # Polyfluorinated components only
   fp = results.to_fingerprint(halogens='F', saturation='poly')

   # No saturation filter (all component SMARTS)
   fp = results.to_fingerprint(halogens='F', saturation=None)

   # generate_fingerprint directly, same new params
   from HalogenGroups.fingerprints import generate_fingerprint
   vectors, info = generate_fingerprint(smiles_list, halogens=['F', 'Cl'], saturation='per')
   # info['halogens'] == ['F', 'Cl']
   # info['saturation'] == 'per'

**Background — group count (116 vs 117):**

The data file contains 117 group entries, but group ID 116 (*Telomers*) has
``compute=False`` — it is an *aggregate* group matched by regex against other group
names rather than parsed directly. The ``@load_HalogenGroups`` decorator excludes it,
so every fingerprint always has **116** computable columns.

Version 2.2.4 (February 2026)
------------------------------

**Released:** February 13, 2026

**New Features:**

- **ResultsFingerprint Class**: New comprehensive fingerprint analysis class with:
  
  - Flexible group selection: 'all', 'oecd', 'generic', 'telomers', or custom groups
  - Multiple encoding modes: 'binary', 'count', 'max_component'
  - Conversion method: ``ResultsModel.to_fingerprint()``

- **Dimensionality Reduction Methods**:
  
  - **PCA** (``perform_pca()``): Linear dimensionality reduction with explained variance analysis and scree plots
  - **Kernel PCA** (``perform_kernel_pca()``): Non-linear dimensionality reduction with multiple kernel options (RBF, polynomial, sigmoid, cosine)
  - **t-SNE** (``perform_tsne()``): Excellent for visualization and cluster identification with configurable perplexity
  - **UMAP** (``perform_umap()``): Fast, scalable non-linear reduction that preserves global structure (requires umap-learn)
  - All methods include automatic plot generation with customizable output

- **Statistical Comparison**:
  
  - **KL Divergence** (``compare_kld()``): Compare fingerprint distributions between datasets
  - Multiple comparison modes: 'minmax' (normalized 0-1), 'forward', 'reverse', 'symmetric'
  - Quantitative assessment of compositional similarity between chemical inventories

- **Database Persistence**:
  
  - ``ResultsFingerprint.to_sql()`` / ``from_sql()``: Efficient save/load with sparse storage
  - ``ResultsModel.from_sql()``: Load previously saved parsing results
  - Support for SQLite and PostgreSQL databases
  - Metadata preservation across save/load cycles

- **Comprehensive Documentation**:
  
  - Complete API reference in ``docs/ResultsFingerprint_Guide.md``
  - Scientific background for each dimensionality reduction method
  - Parameter tuning guidelines and best practices
  - Working examples in ``examples/results_fingerprint_analysis.py``
  - Quick reference guide in ``RESULTS_FINGERPRINT_QUICKREF.md``

**New API Methods:**

.. code-block:: python

   # Convert results to fingerprints
   fp = results.to_fingerprint(group_selection='oecd', count_mode='binary')
   
   # Dimensionality reduction
   pca_results = fp.perform_pca(n_components=5, plot=True)
   kpca_results = fp.perform_kernel_pca(kernel='rbf', plot=True)
   tsne_results = fp.perform_tsne(perplexity=30, plot=True)
   umap_results = fp.perform_umap(n_neighbors=15, plot=True)
   
   # Compare datasets
   kl_divergence = fp1.compare_kld(fp2, method='minmax')
   
   # Save/load
   fp.to_sql(filename='fingerprints.db')
   results.to_sql(filename='results.db')
   fp_loaded = ResultsFingerprint.from_sql(filename='fingerprints.db')
   results_loaded = ResultsModel.from_sql(filename='results.db')

**Use Cases:**

- Exploratory data analysis of large PFAS inventories
- Database comparison and compositional analysis  
- Cluster identification and structural pattern recognition
- Machine learning feature extraction and preprocessing
- Visualization of chemical space
- Quantitative similarity assessment between datasets

**Testing:**

- 100+ comprehensive unit tests covering all new functionality
- Integration tests for complete workflows
- Performance tests for dimensionality reduction methods
- SQL roundtrip tests for data integrity
- Edge case handling and numerical stability tests

**Dependencies:**

- Required: numpy, pandas, scipy, scikit-learn, matplotlib, sqlalchemy
- Optional: umap-learn (for UMAP analysis)

**Documentation Files:**

- ``docs/ResultsFingerprint_Guide.md`` - Complete API documentation
- ``examples/results_fingerprint_analysis.py`` - Working examples
- ``tests/test_results_fingerprint.py`` - Test suite
- ``RESULTS_FINGERPRINT_QUICKREF.md`` - Quick reference

Version 2.2.3 (February 2026)
------------------------------

**Released:** February 2026

**New Features:**

- Added a ``ResultsModel`` container around the default ``parse_mols`` / ``parse_smiles``
   output, providing convenient helpers for navigating PFAS group matches and
   components (e.g. ``show()``, ``summarise()``, ``table()``, and plotting
   functions) while remaining fully compatible with existing list-of-dicts
   consumers and JSON export.
- Extended visualisation utilities for PFAS components:
   - ``plot_pfasgroups`` now supports filtering by component path type
      (e.g. perfluoroalkyl vs polyfluoroalkyl), optional SMARTS filtering,
      panel labels, and more robust handling of invalid SMILES.
   - New figure-generation helpers in the codebase reproduce the main-text and
      supplementary figures used in the PFASgroups manuscript (overview example
      and component/path-type illustrations).
- Improved documentation and manuscript alignment:
   - Updated the main article and Additional File~1 descriptions of PFASgroups
      to match the current implementation (114 PFAS groups, component-based
      graph metrics, universal component merging, fluorotelomer linker
      validation, and regulatory definition support).
   - Clarified JSON-based configuration of PFAS groups and regulatory
      definitions and how component SMARTSs control perfluoroalkyl and
      polyfluoroalkyl detection.



License
-------

PFASgroups is licensed under CC BY-NC 4.0. See :doc:`license` for details.

Support
-------

- **Issues:** https://github.com/yourusername/PFASGroups/issues
- **Discussions:** https://github.com/yourusername/PFASGroups/discussions
- **Email:** luc.miaz@aces.su.se

Citation
--------

If you use PFASgroups in your research, please cite:

   Miaz, L.T., Cousins, I.T. (2026). Automatic Determination and Classification of Per- and Polyfluoroalkyl Substances. *Journal of Cheminformatics* (in preparation).
