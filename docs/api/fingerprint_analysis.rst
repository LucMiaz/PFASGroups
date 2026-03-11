Fingerprint Analysis
====================

.. currentmodule:: PFASGroups

.. code-block:: python

   from PFASGroups import parse_smiles

   results = parse_smiles(["CCCC(F)(F)F", "FC(F)(F)C(=O)O"])
   fp = results.to_fingerprint()

   print(fp.fingerprints.shape)        # (2, 116)  — 116 groups, F only, binary
   print(fp.group_selection)           # 'all'
   print(fp.component_metrics)         # ['binary']

   # PCA in two lines
   coords = fp.perform_pca(n_components=2)
   print(coords['transformed'].shape)  # (2, 2)

The :class:`ResultsFingerprint` class wraps a numpy fingerprint matrix and
provides group selection, three count encoding modes, dimensionality
reduction, distribution comparison, and database I/O.

.. note::

   **Attribute name**: the fingerprint matrix is stored as ``fp.fingerprints``
   (a ``numpy.ndarray``), **not** ``fp.matrix``.

.. contents:: Contents
   :local:
   :depth: 2

ResultsFingerprint
------------------

.. autoclass:: ResultsFingerprint
   :members:
   :undoc-members:
   :show-inheritance:

Creating a ResultsFingerprint
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Via :meth:`~PFASGroups.ResultsModel.to_fingerprint`:

.. code-block:: python

   from PFASGroups import parse_smiles

   smiles = ["CCCC(F)(F)F", "FC(F)(F)C(=O)O", "OCCOCCO"]
   results = parse_smiles(smiles)

   # Default: all 116 groups, fluorine only, binary encoding
   fp = results.to_fingerprint()
   print(fp.fingerprints.shape)     # (3, 116)
   print(fp.halogens)               # ['F']
   print(fp.group_names[:2])        # e.g. ['perfluoromethyl', 'perfluoroalkyl']

   # OECD groups only, count mode
   fp_oecd = results.to_fingerprint(group_selection='oecd', component_metrics=['count'])
   print(fp_oecd.fingerprints.shape)   # (3, 28)

Or use :func:`~PFASGroups.generate_fingerprint` directly for a numpy array
without the :class:`ResultsFingerprint` wrapper:

.. code-block:: python

   from PFASGroups import generate_fingerprint

   fps, info = generate_fingerprint(smiles)
   print(fps.shape)                  # (3, 116)
   print(info['group_names'][:2])    # list of group name strings

Attributes
~~~~~~~~~~~

.. list-table::
   :header-rows: 1
   :widths: 25 75

   * - Attribute
     - Description
   * - ``fingerprints``
     - ``numpy.ndarray`` of shape ``(n_mols, n_groups)``
   * - ``smiles``
     - List of SMILES (one per row)
   * - ``group_names``
     - List of group-name strings (one per column)
   * - ``group_selection``
     - Selection used: ``'all'``, ``'oecd'``, ``'generic'``, ``'telomers'``, or ``'generic+telomers'``
   * - ``halogens``
     - List of halogen symbols used, e.g. ``['F']`` or ``['F', 'Cl', 'Br', 'I']``
   * - ``saturation``
     - Saturation filter used (``'per'``, ``'poly'``, or ``None``)
   * - ``component_metrics``
     - List of metrics used: e.g. ``['binary']`` or ``['binary', 'effective_graph_resistance']``

Dimensionality reduction
-------------------------

perform_pca
~~~~~~~~~~~

.. automethod:: ResultsFingerprint.perform_pca
   :no-index:

.. code-block:: python

   result = fp.perform_pca(n_components=2, plot=True)
   coords = result['transformed']   # numpy array (n_mols, 2)

**Parameters:**

- ``n_components`` (int, default 2): Number of PCA components
- ``plot`` (bool, default True): Generate and display a scatter plot
- ``output_file`` (str, optional): File path to save the plot

*Requires*: scikit-learn

perform_tsne
~~~~~~~~~~~~

.. automethod:: ResultsFingerprint.perform_tsne
   :no-index:

.. code-block:: python

   result = fp.perform_tsne(n_components=2, perplexity=5, plot=True)
   coords = result['transformed']

**Parameters:**

- ``n_components`` (int, default 2)
- ``perplexity`` (float, default 30.0): t-SNE perplexity (typically 5–50;
  must be less than the number of molecules)
- ``n_iter`` (int, default 1000)
- ``plot`` (bool, default True)

*Requires*: scikit-learn

perform_umap
~~~~~~~~~~~~

.. automethod:: ResultsFingerprint.perform_umap
   :no-index:

.. code-block:: python

   result = fp.perform_umap(n_components=2, n_neighbors=15, plot=True)
   coords = result['transformed']

**Parameters:**

- ``n_components`` (int, default 2)
- ``n_neighbors`` (int, default 15): UMAP local neighborhood size
- ``min_dist`` (float, default 0.1): Minimum distance between embedded points
- ``plot`` (bool, default True)

*Requires*: umap-learn (``pip install umap-learn``)

Statistical comparison
-----------------------

compare_kld
~~~~~~~~~~~

.. automethod:: ResultsFingerprint.compare_kld
   :no-index:

Compute per-group KL divergence between two fingerprint distributions:

.. code-block:: python

   fp_a = parse_smiles(set_a).to_fingerprint()
   fp_b = parse_smiles(set_b).to_fingerprint()

   result = fp_a.compare_kld(fp_b, method='minmax')
   # result: dict with 'kl_divergences', 'sorted_groups', 'total_kl', ...

**Parameters:**

- ``other`` (:class:`ResultsFingerprint`): The comparison fingerprint
- ``method`` (str, default ``'minmax'``):
  - ``'minmax'``: normalise each column to [0, 1] before comparison
  - ``'forward'``: KL(self ‖ other)
  - ``'reverse'``: KL(other ‖ self)
  - ``'symmetric'``: average of forward + reverse

Database I/O
------------

to_sql / from_sql
~~~~~~~~~~~~~~~~~

.. automethod:: ResultsFingerprint.to_sql
   :no-index:
.. automethod:: ResultsFingerprint.from_sql
   :no-index:

.. code-block:: python

   # Save to SQLite
   fp.to_sql("fingerprints.db")

   # Load back
   from PFASGroups.results_model import ResultsFingerprint
   fp2 = ResultsFingerprint.from_sql("fingerprints.db")
   print(fp2.fingerprints.shape)

Matrices are stored in sparse format.  Metadata (``halogens``,
``saturation``, ``component_metrics``, ``group_names``) is preserved.

Summary
~~~~~~~

.. automethod:: ResultsFingerprint.summary
   :no-index:

.. code-block:: python

   print(fp.summary())
   # ResultsFingerprint Summary
   # ==========================
   # Molecules: 3
   # Groups: 116
   # Group selection: all
   # Halogens: F
   # Comp metr. : ['binary']
   # Fingerprint shape: (3, 116)
   # Non-zero entries: 2