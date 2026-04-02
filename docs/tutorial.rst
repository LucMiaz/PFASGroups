Tutorial
========

.. code-block:: python

   from PFASGroups import parse_smiles, prioritise_molecules, plot_HalogenGroups

   results  = parse_smiles(["CCCC(F)(F)F", "FC(F)(F)C(=O)O"])
   arr      = results.to_array()                     # (2, 114) binary embedding
   pca      = results.perform_pca(color_by='top_group')
   ranked   = prioritise_molecules(["CCCC(F)(F)F", "FC(F)(F)C(=O)O"])

.. contents:: Table of contents
   :local:
   :depth: 2

Sample data
-----------

All examples on this page share the following molecule list.  Copy it once at
the top of your script or notebook:

.. code-block:: python

   smiles_list = [
       "C[Si](C)(Cl)CCC(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)F",
       "FC(=C(C(F)(F)F)C(F)(F)F)C(F)(F)C(F)(F)F",
       "FC(F)(F)C(F)(C(F)(F)F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(C(F)(F)F)C(F)(F)F",
       "C[Si](Cl)(Cl)CCC(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)F",
       "FC(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)F",
       "CC(=O)C(F)(F)C(F)(F)C(F)(F)F",
       "FC(F)=C(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)F",
       "FC(F)(F)N1C(F)(F)C(F)(F)C(F)(F)C(F)(F)C1(F)F",
       "FC(F)(F)C(F)(F)C(F)(F)CCl",
       "OCC(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)F",
       "OS(=O)(=O)CCC(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)F",
       "FC(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)COC(=O)C=C",
       "FC(F)(F)C(F)(F)C(F)(F)C(F)(F)C(=O)OC=C",
       "FC(F)(F)C(I)(C(F)(F)F)C(F)(F)F",
       "OC(C(F)(F)C(F)(F)F)(C(F)(F)C(F)(F)F)C(F)(F)C(F)(F)F",
       "C[Si](Cl)(Cl)CCC(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)F",
       "FC(F)=C(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)F",
       "CCN(CC)CC.NS(=O)(=O)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)F",
       "[OH-].C[N+](C)(C)CCCNS(=O)(=O)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)F",
       "[K+].[O-]S(=O)(=O)C(F)(F)C(F)(C(F)(F)F)C(F)(F)F",
       "OCCN(CCO)S(=O)(=O)C(F)(F)C(F)(F)C(F)(F)C(F)F",
       "FC(F)(F)C(F)(F)C(F)(F)C(F)(F)N(C(F)(F)C(F)(F)C(F)(F)C(F)(F)F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)F",
       "COC(=O)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)F",
       "FC(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(=O)OCC=C",
       "CC(C)CCOC(=O)CC(NC(=O)C(F)(F)C(F)(F)C(F)(F)F)C(=O)OCCC(C)C",
       "CN(C(=O)C(F)(F)C(F)(F)C(F)(F)F)C(=O)C(F)(F)C(F)(F)C(F)(F)F",
       "FC(F)(F)C(F)(F)C(F)(F)C=CI",
       "OC(=O)C1=CC=CC=C1NC(=O)C(F)(F)C(F)(F)C(F)(F)F",
       "FC(F)(F)C1(F)C(F)(F)C(F)(F)C(F)(F)C(F)(C(F)(F)F)C1(F)F",
   ]

Parsing SMILES
--------------

:func:`~PFASGroups.parse_smiles` accepts a single SMILES string or a list and
returns a :class:`~PFASGroups.PFASEmbeddingSet`:

.. code-block:: python

   from PFASGroups import parse_smiles

   results = parse_smiles(smiles_list)
   print(results)        # PFASEmbeddingSet summary: molecule count, top groups, …
   print(results[0])     # per-molecule summary for the first entry

   results.summary()     # tabular overview of matched groups across all molecules

Iterating over results
~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   for mol in results:
       flag = " (PFAS)" if any(m.is_PFAS for m in mol.matches) else ""
       print(f"{mol.smiles}: {len(mol.matches)} match(es){flag}")

Drilling into a single match
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   mol   = results[0]
   match = mol.matches[0]

   print("Group   :", match.group_name)
   print("Group ID:", match.group_id)
   print("Is PFAS :", match.is_PFAS)
   print("Category:", match.group_category)

   for i, comp in enumerate(match.components):
       print(f"  Component {i}: atoms {comp.atoms}")
       print(f"  EGR: {comp.effective_graph_resistance}")

Exporting to DataFrame
~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   df = results.to_dataframe()
   # columns: smiles, inchi, group_name, group_id, is_PFAS,
   #          group_category, n_components, component_atoms, …
   df.to_csv("results.csv", index=False)

Parsing RDKit Mol objects
--------------------------

:func:`~PFASGroups.parse_mols` accepts a list of RDKit ``Mol`` objects directly —
useful when molecules are already in memory (e.g. read from an SD file):

.. code-block:: python

   from rdkit import Chem
   from PFASGroups import parse_mols

   mols = [Chem.MolFromSmiles(s) for s in smiles_list]
   results = parse_mols(mols)

   print(results)
   results.summary()

The returned :class:`~PFASGroups.PFASEmbeddingSet` has the same interface as
the one returned by :func:`~PFASGroups.parse_smiles`.

Plotting results
-----------------

:func:`~PFASGroups.plot_HalogenGroups` draws 2-D molecular structures with
each detected PFAS group highlighted:

.. code-block:: python

   from PFASGroups import plot_HalogenGroups

   # First three molecules, group highlights, 2-column grid
   fig, w, h = plot_HalogenGroups(smiles_list[:3], ncols=2, display=True)

   # Save to a file without displaying
   fig, w, h = plot_HalogenGroups(smiles_list[:3], path="groups.png", display=False)

   # Vector SVG output
   fig, w, h = plot_HalogenGroups(smiles_list[:3], svg=True, display=True)

   # One panel per group match (useful when a molecule has multiple groups)
   fig, w, h = plot_HalogenGroups(smiles_list[:3], split_matches=True, ncols=3)

   # Custom panel labels
   labels = [f"mol {i}" for i in range(3)]
   fig, w, h = plot_HalogenGroups(smiles_list[:3], panel_labels=labels, ncols=3)

For plain structure drawings without group highlights, use
:func:`~PFASGroups.plot_mol` (single) or :func:`~PFASGroups.plot_mols` (grid):

.. code-block:: python

   from rdkit import Chem
   from PFASGroups import plot_mol, plot_mols

   mol  = Chem.MolFromSmiles(smiles_list[0])
   mols = [Chem.MolFromSmiles(s) for s in smiles_list[:4]]

   fig, w, h = plot_mol(mol, addAtomIndices=False)
   fig, w, h = plot_mols(mols, ncols=2, addAtomIndices=False)

Component filters
-----------------

Filter parsed results by halogen, saturation level, or structural form:

.. code-block:: python

   # Fluorine only (default)
   r_f    = parse_smiles(smiles_list, halogens='F')

   # Perfluorinated alkyl chains only
   r_pfa  = parse_smiles(smiles_list, halogens='F',
                          saturation='per', form='alkyl')

   # Cyclic fluorinated structures only
   r_cyc  = parse_smiles(smiles_list, form='cyclic')

   # Polyfluorinated components only
   r_poly = parse_smiles(smiles_list, saturation='poly')

   # All four halogens
   r_all  = parse_smiles(smiles_list, halogens=['F', 'Cl', 'Br', 'I'])

Embeddings
----------

:meth:`~PFASGroups.PFASEmbeddingSet.to_array` converts a
:class:`~PFASGroups.PFASEmbeddingSet` into a 2-D
:class:`~PFASGroups.EmbeddingArray` (a numpy array subclass that carries
molecule identity metadata).  All examples reuse the pre-parsed ``results``
object — no re-parsing needed.

.. code-block:: python

   # Default: binary encoding, all groups → shape (n_mols, 114)
   arr = results.to_array()
   print(arr.shape)      # (29, 114)
   print(arr.smiles)     # list of SMILES strings, one per row

   # OECD groups only
   arr_oecd = results.to_array(group_selection='oecd')          # (n, 27)

   # Count or max-component size encoding
   arr_cnt  = results.to_array(component_metrics=['count'])
   arr_max  = results.to_array(component_metrics=['max_component'])

   # Preset 'best': binary + effective graph resistance
   # (best mean Tanimoto in benchmarks — 2 × 114 columns)
   arr_best = results.to_array(preset='best')                   # (n, 228)

   # Stacked multi-metric embedding
   arr_multi = results.to_array(
       component_metrics=['binary', 'effective_graph_resistance',
                          'n_spacer', 'ring_size'],
       molecule_metrics=['n_components', 'max_size',
                         'mean_branching', 'mean_component_fraction'],
   )

   # Column labels matching the array columns
   cols = results.column_names(
       component_metrics=['binary', 'effective_graph_resistance']
   )
   print(cols[:4])

.. note::

   ``n_spacer`` encodes the telomer CH\ :sub:`2` spacer length (the *m* in
   *m:n* telomer notation).  It is zero for non-telomers.

   ``ring_size`` is the smallest ring overlapping each matched component
   (zero for acyclic groups; 5 for azoles/furans; 6 for benzene/cyclohexane).

Saving to a database
~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   results.to_sql(filename='results.db')   # SQLite by default

Dimensionality reduction
------------------------

PCA, t-SNE, and UMAP operate directly on a
:class:`~PFASGroups.PFASEmbeddingSet` — no manual call to ``to_array()``
needed.  Each method returns a ``dict`` whose ``'transformed'`` key holds the
projected coordinates as a numpy array.

.. note::

   PCA and t-SNE require ``scikit-learn``; UMAP requires ``umap-learn``::

      pip install scikit-learn umap-learn matplotlib

PCA
~~~

.. code-block:: python

   # Basic scatter plot — displayed automatically
   pca = results.perform_pca(n_components=2, plot=True)
   print(pca['transformed'].shape)   # (29, 2)

   # Colour each point by the PFAS group with the highest match count
   pca = results.perform_pca(n_components=2, color_by='top_group')

   # Supply your own per-molecule labels (e.g. from external metadata)
   labels = ['industrial'] * 15 + ['environmental'] * (len(smiles_list) - 15)
   pca = results.perform_pca(n_components=2, color_by=labels)

   # Save figure to a file without displaying it
   pca = results.perform_pca(n_components=2, color_by='top_group',
                              output_file='pca.png', plot=False)

t-SNE
~~~~~

.. code-block:: python

   # Basic t-SNE scatter plot
   tsne = results.perform_tsne(n_components=2, perplexity=8, plot=True)
   print(tsne['transformed'].shape)  # (29, 2)

   # Colour by top PFAS group
   tsne = results.perform_tsne(perplexity=8, color_by='top_group')

   # Custom per-molecule labels
   tsne = results.perform_tsne(perplexity=8, color_by=labels)

   # Save figure
   tsne = results.perform_tsne(perplexity=8, color_by='top_group',
                                output_file='tsne.png', plot=False)

UMAP
~~~~

.. code-block:: python

   # Basic UMAP scatter plot
   umap_res = results.perform_umap(n_components=2, n_neighbors=10, plot=True)
   print(umap_res['transformed'].shape)  # (29, 2)

   # Colour by top PFAS group
   umap_res = results.perform_umap(n_neighbors=10, color_by='top_group')

   # Custom per-molecule labels
   umap_res = results.perform_umap(n_neighbors=10, color_by=labels)

   # Save figure
   umap_res = results.perform_umap(n_neighbors=10, color_by='top_group',
                                   output_file='umap.png', plot=False)

Prioritisation
--------------

:func:`~PFASGroups.prioritise_molecules` ranks molecules by structural
similarity to a reference set, or by intrinsic fluorination properties when no
reference is provided.  It accepts SMILES lists, RDKit ``Mol`` lists, or a
pre-parsed :class:`~PFASGroups.PFASEmbeddingSet`.

Ranking against a reference set
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Molecules are ranked by KL-divergence similarity to the reference distribution
(lower divergence = more similar = higher priority):

.. code-block:: python

   from PFASGroups import prioritise_molecules

   reference = [
       "FC(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(=O)O",
   ]

   ranked, scores = prioritise_molecules(
       molecules=smiles_list,
       reference=reference,
       return_scores=True,
   )

   for mol, score in zip(ranked, scores):
       print(f"{score:.4f}  {mol.smiles}")

Ranking by intrinsic fluorination
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

When no reference is given, molecules are scored by total fluorination (weight
``a``) and longest chain length (weight ``b``, at the given ``percentile``):

.. code-block:: python

   ranked, scores = prioritise_molecules(
       molecules=smiles_list,
       a=1.0,          # weight for total fluorinated content
       b=5.0,          # weight for longest chain (percentile-based)
       percentile=90,  # focus on the 90th-percentile component size
       return_scores=True,
   )

   for mol, score in zip(ranked, scores):
       print(f"{score:.4f}  {mol.smiles}")

Passing a pre-parsed set avoids double-parsing large inventories:

.. code-block:: python

   ranked, scores = prioritise_molecules(
       molecules=results,       # PFASEmbeddingSet already computed above
       a=1.0, b=5.0, percentile=90,
       return_scores=True,
   )

Comparing two datasets
-----------------------

:meth:`~PFASGroups.PFASEmbeddingSet.compare_kld` returns a single normalised
similarity score (0 = identical distributions, 1 = maximally different):

.. code-block:: python

   other_smiles = smiles_list[:15] + [
       "OC(=O)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(O)=O",
       "OCC(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)F",
   ]
   other_results = parse_smiles(other_smiles)

   sim = results.compare_kld(other_results, method='minmax')
   print(f"Similarity: {sim:.3f}")   # lower = more similar

Multi-halogen analysis
-----------------------

By default ``PFASGroups`` detects fluorine only.  Pass ``halogens`` explicitly
to include other halogens, or import from ``HalogenGroups`` where all four are
the default:

.. code-block:: python

   # Explicit multi-halogen parse via PFASGroups
   r_mh = parse_smiles(smiles_list, halogens=['F', 'Cl', 'Br', 'I'])

   # Multi-halogen embedding: 114 groups × 4 halogens = 456 columns
   arr_mh = r_mh.to_array(halogens=['F', 'Cl', 'Br', 'I'])
   print(arr_mh.shape)   # (29, 456)

See :doc:`halogengroups` for the full multi-halogen reference.

PFAS definition screening
--------------------------

.. code-block:: python

   results_defs = parse_smiles(smiles_list, include_PFAS_definitions=True)

   for mol in results_defs:
       if mol.pfas_definition_matches:
           names = [d.definition_name for d in mol.pfas_definition_matches]
           print(mol.smiles, "→", names)

Available definitions: OECD 2021, EU REACH, OPPT 2023, UK Environment Agency,
PFASSTRUCTv5.  See :doc:`pfas_definitions` for details.

Loading and inspecting groups
------------------------------

.. code-block:: python

   from PFASGroups import get_compiled_HalogenGroups

   groups = get_compiled_HalogenGroups()
   print(len(groups))        # 114 compiled groups (plus aggregate groups)

   for g in groups[:5]:
       print(g.name, g.category, g.id)

See :doc:`customization` for how to add custom group definitions.