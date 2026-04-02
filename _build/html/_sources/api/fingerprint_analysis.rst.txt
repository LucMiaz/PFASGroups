Embedding Analysis
==================

.. currentmodule:: PFASGroups

.. code-block:: python

   from PFASGroups import parse_smiles

   results = parse_smiles(["CCCC(F)(F)F", "FC(F)(F)C(=O)O"])
   arr = results.to_array()

   print(arr.shape)                        # (2, 112)  — 112 groups, F only, binary
   cols = results.column_names()           # list of 112 column label strings

   # PCA in two lines
   pca_result = results.perform_pca(n_components=2)
   print(pca_result['transformed'].shape)  # (2, 2)

:func:`parse_smiles` returns a :class:`PFASEmbeddingSet` — a list of
:class:`PFASEmbedding` objects, one per molecule.  Call :meth:`~PFASEmbeddingSet.to_array`
on the set to get a ``(n_mols, n_cols)`` numpy matrix, or call it on a single
:class:`PFASEmbedding` for a 1-D vector.

.. seealso::

   :doc:`/ResultsFingerprint_Guide` — complete metric reference with formulas,
   preset benchmarks, and group selection tables.

.. contents:: Contents
   :local:
   :depth: 2

PFASEmbeddingSet
----------------

.. autoclass:: PFASEmbeddingSet
   :members:
   :undoc-members:
   :show-inheritance:

Generating an embedding array
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Via :meth:`~PFASGroups.PFASEmbeddingSet.to_array`:

.. code-block:: python

   from PFASGroups import parse_smiles

   smiles = ["CCCC(F)(F)F", "FC(F)(F)C(=O)O", "OCCOCCO"]
   results = parse_smiles(smiles)

   # Default: all 112 groups, fluorine only, binary encoding
   arr = results.to_array()
   print(arr.shape)                      # (3, 112)
   cols = results.column_names()         # list of 112 column label strings

   # OECD groups only, count mode
   arr_oecd = results.to_array(group_selection='oecd', component_metrics=['count'])
   print(arr_oecd.shape)                 # (3, 28)

   # Best-performing preset (binary + effective_graph_resistance)
   arr_best = results.to_array(preset='best')
   cols_best = results.column_names(preset='best')
   print(arr_best.shape)                 # (3, 224)  — 112 groups × 2 metrics

   # Single-molecule embedding
   vec = results[0].to_array(preset='best')   # 1-D array, length 224



Key attributes
~~~~~~~~~~~~~~

:class:`PFASEmbeddingSet` is a plain list of :class:`PFASEmbedding` dicts.
The embedding matrix is computed on demand; no matrix is stored on the object.

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - Attribute / method
     - Description
   * - ``results.to_array(...)``
     - ``numpy.ndarray`` of shape ``(n_mols, n_cols)``
   * - ``results.column_names(...)``
     - List of column label strings; same arguments as ``to_array()``
   * - ``results[i]``
     - :class:`PFASEmbedding` for molecule *i* (dict subclass)
   * - ``results[i].smiles``
     - SMILES string of molecule *i*
   * - ``results[i].to_array(...)``
     - 1-D embedding vector for one molecule
   * - ``results[i].summarise()``
     - Returns a formatted string summary of matched PFAS groups
   * - ``results[i].summary()``
     - Prints the formatted summary to stdout

Dimensionality reduction
-------------------------

perform_pca
~~~~~~~~~~~

.. automethod:: PFASEmbeddingSet.perform_pca
   :no-index:

.. code-block:: python

   result = results.perform_pca(n_components=2, plot=True)
   coords = result['transformed']    # numpy array (n_mols, 2)
   evr   = result['explained_variance']   # variance ratios per component

**Parameters:**

- ``n_components`` (int, default 2): Number of PCA components
- ``plot`` (bool, default True): Generate and display a scatter plot + scree plot
- ``output_file`` (str, optional): File path to save the plot

**Returns:** dict with keys ``'transformed'``, ``'explained_variance'``,
``'components'``, ``'pca_model'``, ``'scaler'``.

*Requires*: scikit-learn, matplotlib

perform_tsne
~~~~~~~~~~~~

.. automethod:: PFASEmbeddingSet.perform_tsne
   :no-index:

.. code-block:: python

   result = results.perform_tsne(n_components=2, perplexity=5, plot=True)
   coords = result['transformed']

**Parameters:**

- ``n_components`` (int, default 2)
- ``perplexity`` (float, default 30.0): t-SNE perplexity (typically 5–50;
  must be less than the number of molecules)
- ``learning_rate`` (float, default 200.0)
- ``max_iter`` (int, default 1000)
- ``plot`` (bool, default True)

**Returns:** dict with keys ``'transformed'``, ``'tsne_model'``, ``'scaler'``,
``'perplexity'``.

*Requires*: scikit-learn, matplotlib

perform_umap
~~~~~~~~~~~~

.. automethod:: PFASEmbeddingSet.perform_umap
   :no-index:

.. code-block:: python

   result = results.perform_umap(n_components=2, n_neighbors=15, plot=True)
   coords = result['transformed']

**Parameters:**

- ``n_components`` (int, default 2)
- ``n_neighbors`` (int, default 15): UMAP local neighborhood size
- ``min_dist`` (float, default 0.1): Minimum distance between embedded points
- ``metric`` (str, default ``'euclidean'``)
- ``plot`` (bool, default True)

**Returns:** dict with keys ``'transformed'``, ``'umap_model'``, ``'scaler'``,
``'n_neighbors'``, ``'min_dist'``.

*Requires*: umap-learn (``pip install umap-learn``), matplotlib

Statistical comparison
-----------------------

compare_kld
~~~~~~~~~~~

.. automethod:: PFASEmbeddingSet.compare_kld
   :no-index:

Compute KL divergence between the group-occurrence frequencies of two sets:

.. code-block:: python

   results_a = parse_smiles(set_a)
   results_b = parse_smiles(set_b)

   kld = results_a.compare_kld(results_b, method='minmax')
   # kld: float — normalised symmetric KLD (lower = more similar distributions)

**Parameters:**

- ``other`` (:class:`PFASEmbeddingSet`): The comparison set
- ``method`` (str, default ``'minmax'``):

  - ``'minmax'``: normalised symmetric KLD ∈ [0, 1]
  - ``'forward'``: KL(self ‖ other)
  - ``'reverse'``: KL(other ‖ self)
  - ``'symmetric'``: average of forward + reverse

**Returns:** ``float``

Database I/O
------------

to_sql / from_sql
~~~~~~~~~~~~~~~~~

.. automethod:: PFASEmbeddingSet.to_sql
   :no-index:
.. automethod:: PFASEmbeddingSet.from_sql
   :no-index:

.. code-block:: python

   # Save to SQLite
   results.to_sql(filename="pfas_results.db")

   # Load back
   from PFASGroups import PFASEmbeddingSet
   results2 = PFASEmbeddingSet.from_sql(filename="pfas_results.db")

   # PostgreSQL
   results.to_sql(dbname="mydb", user="alice", password="secret", host="localhost")

Component and group-level data are stored in two tables
(``components`` and ``pfas_groups_in_compound`` by default).
Pass ``if_exists='replace'`` to overwrite existing tables.