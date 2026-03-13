Tutorial
========

.. code-block:: python

   from PFASGroups import parse_smiles, generate_fingerprint, prioritise_molecules

   results  = parse_smiles(["CCCC(F)(F)F", "FC(F)(F)C(=O)O"])
   fps, _   = generate_fingerprint(["CCCC(F)(F)F", "FC(F)(F)C(=O)O"])
   ranking  = prioritise_molecules(["CCCC(F)(F)F", "FC(F)(F)C(=O)O"])

.. contents:: Table of contents
   :local:
   :depth: 2

Setting up
----------

.. code-block:: python

   from PFASGroups import (
       parse_smiles,
       generate_fingerprint,
       get_compiled_HalogenGroups,
       prioritise_molecules,
   )
   import pandas as pd

Parsing a list of molecules
----------------------------

.. code-block:: python

   smiles_list = [
       "CCCC(F)(F)F",
       "FC(F)(F)C(=O)O",
       "OCCOCCO",      # no fluorine
   ]

   results = parse_smiles(smiles_list)
   print(len(results))          # 3 — one per input SMILES

Iterating over results
----------------------

.. code-block:: python

   for mol in results:
       n = len(mol.matches)
       flag = " (PFAS)" if any(m.is_PFAS for m in mol.matches) else ""
       print(f"{mol.smiles}: {n} match(es){flag}")

Drill into a specific match
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   mol = results[0]   # CCCC(F)(F)F

   match = mol.matches[0]
   print("Group   :", match.group_name)
   print("Group ID:", match.group_id)
   print("Is PFAS :", match.is_PFAS)
   print("Category:", match.group_category)

   # atom indices of each matched component
   for i, comp in enumerate(match.components):
       print(f"  Component {i}: atoms {comp.atoms}")
       print(f"  EGR: {comp.effective_graph_resistance}")

Exporting to DataFrame
-----------------------

.. code-block:: python

   df = results.to_dataframe()
   # columns: smiles, inchi, group_name, group_id, is_PFAS,
   #          group_category, n_components, component_atoms, ...
   df.to_csv("results.csv", index=False)

Filtering by group category
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   pfas_df = df[df["is_PFAS"] == True]
   oecd_df  = df[df["group_category"] == "OECD"]

Saturation filter
-----------------

.. code-block:: python

   # Perfluorinated components only (uses 'per' SMARTS)
   results_per  = parse_smiles(smiles_list, saturation='per')

   # Polyfluorinated components only
   results_poly = parse_smiles(smiles_list, saturation='poly')

   # No saturation filter (use all components)
   results_all  = parse_smiles(smiles_list, saturation=None)

Fingerprinting for machine learning
-------------------------------------

:func:`~PFASGroups.generate_fingerprint` converts a list of SMILES to a
2-D numpy array.  By default the output has **116 columns** (one per group,
fluorine only, binary encoding):

.. code-block:: python

   fps, info = generate_fingerprint(smiles_list)
   print(fps.shape)                  # (3, 116)
   print(info['group_names'][:3])    # e.g. ['perfluoromethyl', ...]

**Selecting groups** — pass a list or range of 0-based indices:

.. code-block:: python

   # OECD groups only (indices 0–27 map to IDs 1–28)
   fps_oecd, _ = generate_fingerprint(smiles_list, selected_groups=range(0, 28))

**Count modes:**

.. code-block:: python

   fps_bin, _ = generate_fingerprint(smiles_list, component_metrics=['binary'])       # default
   fps_cnt, _ = generate_fingerprint(smiles_list, component_metrics=['count'])
   fps_max, _ = generate_fingerprint(smiles_list, component_metrics=['max_component'])

   # Best benchmark config: binary + effective graph resistance (2 × n_groups cols)
   fps_best, _ = generate_fingerprint(smiles_list, preset='best')

Using ResultsModel.to_fingerprint
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

:meth:`~PFASGroups.ResultsModel.to_fingerprint` produces a
:class:`~PFASGroups.ResultsFingerprint` with the same matrix plus metadata:

.. code-block:: python

   results = parse_smiles(smiles_list)
   fp = results.to_fingerprint()

   print(fp.fingerprints.shape)     # (3, 116)
   print(fp.group_names[:3])        # list of group-name strings
   print(fp.halogens)               # ['F']
   print(fp.component_metrics)      # ['binary']

Group selection with to_fingerprint:

.. code-block:: python

   fp_oecd = results.to_fingerprint(group_selection='oecd')       # 28 cols
   fp_gen  = results.to_fingerprint(group_selection='generic')    # 45 cols
   fp_all  = results.to_fingerprint(group_selection='all')        # 116 cols

Dimensionality reduction
--------------------------

.. code-block:: python

   fp = results.to_fingerprint()

   coords_pca  = fp.perform_pca(n_components=2)
   coords_tsne = fp.perform_tsne(n_components=2, perplexity=5)
   coords_umap = fp.perform_umap(n_components=2)

   # each returns a dict with key 'transformed' (numpy array, shape n_mols × n_components)

Comparing distributions
~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   fp_a = parse_smiles(smiles_list[:2]).to_fingerprint()
   fp_b = parse_smiles(smiles_list[2:]).to_fingerprint()

   kld = fp_a.compare_kld(fp_b)   # dict of per-group KL divergences

Multi-halogen parsing (advanced)
----------------------------------

By default ``PFASGroups`` detects fluorine only.  To detect other halogens
pass ``halogens`` explicitly:

.. code-block:: python

   results = parse_smiles(["ClCCCl", "BrCCBr"], halogens=['F', 'Cl', 'Br', 'I'])

   # Multi-halogen fingerprint: 116 * 4 = 464 columns
   fps, info = generate_fingerprint(
       ["ClCCCl", "BrCCBr"], halogens=['F', 'Cl', 'Br', 'I'])
   print(fps.shape)     # (2, 464)

See :doc:`halogengroups` for full documentation of multi-halogen mode.

PFAS definition screening
--------------------------

.. code-block:: python

   results = parse_smiles(smiles_list, include_PFAS_definitions=True)

   for mol in results:
       if mol.pfas_definition_matches:
           names = [d.definition_name for d in mol.pfas_definition_matches]
           print(mol.smiles, "->", names)

Available definitions: OECD 2021, EU REACH, OPPT 2023, UK Environment Agency,
PFASTRUCTv5.  See :doc:`pfas_definitions` for details.

Molecule prioritization
-------------------------

.. code-block:: python

   from PFASGroups import prioritise_molecules

   molecules = ["CCCC(F)(F)F", "FC(F)(F)C(=O)O", "FCCC(F)(F)F", "ClCCCl"]
   ranking = prioritise_molecules(molecules)

   for smiles, score in ranking:
       print(smiles, score)

See :doc:`prioritization` for parameter details.

Loading and inspecting groups
------------------------------

.. code-block:: python

   from PFASGroups import get_compiled_HalogenGroups

   groups = get_compiled_HalogenGroups()
   print(len(groups))       # 116

   for g in groups[:5]:
       print(g.name, g.category, g.smarts)

Custom group collections
-------------------------

See :doc:`customization` for how to create and use custom halogen group
definitions.