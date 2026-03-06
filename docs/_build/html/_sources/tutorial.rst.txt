Tutorial
========

This tutorial walks through typical research workflows: parsing a chemical
database, generating and comparing fingerprints, applying PFAS definitions,
and exporting data.

.. contents:: Table of contents
   :local:
   :depth: 2

Setting up
----------

.. code-block:: python

   from HalogenGroups import (
       parse_smiles,
       generate_fingerprint,
       get_compiled_HalogenGroups,
   )
   import pandas as pd

Parsing a list of molecules
----------------------------

.. code-block:: python

   smiles_list = [
       "CCCC(F)(F)F",
       "FC(F)(F)C(=O)O",
       "ClC(Cl)Cl",
       "BrCCBr",
       "IC(I)(I)I",
       "OCCOCCO",      # no halogen
   ]

   results = parse_smiles(smiles_list)
   print(len(results))          # 6 — one per input SMILES

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
       print(f"  Effective graph resistance: {comp.effective_graph_resistance}")

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

Working with multiple halogens
-------------------------------

.. code-block:: python

   # All halogens (HalogenGroups default)
   results_all = parse_smiles(smiles_list, halogens=['F', 'Cl', 'Br', 'I'])

   # Just chlorine
   results_cl = parse_smiles(smiles_list, halogens='Cl')

   # Fluorine + chlorine
   results_fcl = parse_smiles(smiles_list, halogens=['F', 'Cl'])

Saturation filter
~~~~~~~~~~~~~~~~~

.. code-block:: python

   # Only detect fully saturated halogenated groups
   results_sat = parse_smiles(smiles_list, saturation='saturated')

   # Only unsaturated (e.g. fluorotelomer olefins)
   results_unsat = parse_smiles(smiles_list, saturation='unsaturated')

Fingerprinting for machine learning
-------------------------------------

.. code-block:: python

   fps, group_names = generate_fingerprint(smiles_list)
   print(fps.shape)    # (6, 464)  — 116 groups × 4 halogens

   # fluorine-only fingerprint (116 columns)
   fps_f, _ = generate_fingerprint(smiles_list, halogens='F')
   print(fps_f.shape)  # (6, 116)

count_mode controls how overlapping matches within a molecule are counted:

.. code-block:: python

   # default: max component count per group
   fps_max, _ = generate_fingerprint(smiles_list, count_mode='max_component')

   # total occurrences
   fps_all, _ = generate_fingerprint(smiles_list, count_mode='all')

   # binary (0/1 presence)
   fps_bin, _ = generate_fingerprint(smiles_list, count_mode='binary')

Using ResultsModel.to_fingerprint
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   results = parse_smiles(smiles_list)
   fp_obj = results.to_fingerprint()

   # fp_obj is a ResultsFingerprint
   print(fp_obj.matrix.shape)   # (6, 464)
   print(fp_obj.group_names)    # dict: col_index -> (group_name, halogen)

Dimensionality reduction
--------------------------

.. code-block:: python

   fp_obj = results.to_fingerprint()

   coords_pca   = fp_obj.perform_pca(n_components=2)
   coords_tsne  = fp_obj.perform_tsne(n_components=2, perplexity=5)
   coords_umap  = fp_obj.perform_umap(n_components=2)

   # each returns a numpy array of shape (n_mols, n_components)

Comparing distributions
~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   # KL divergence between two fingerprint sets
   fp_a = parse_smiles(smiles_list[:3]).to_fingerprint()
   fp_b = parse_smiles(smiles_list[3:]).to_fingerprint()

   kld = fp_a.compare_kld(fp_b)   # returns dict of per-group KL divergences

PFAS definition screening
--------------------------

.. code-block:: python

   results = parse_smiles(smiles_list, include_PFAS_definitions=True)

   for mol in results:
       if mol.pfas_definition_matches:
           names = [d.definition_name for d in mol.pfas_definition_matches]
           print(mol.smiles, "->", names)

Available definitions:

- ``OECD 2021`` — Organisation for Economic Co-operation and Development
- ``EU REACH`` — European Union REACH regulation
- ``OPPT 2023`` — US EPA Office of Pollution Prevention and Toxics
- ``UK Environment Agency`` — UK Fluorinated Polymer definition
- ``PFASTRUCTv5`` — Annotated PFAS chemical structure database

Molecule prioritization
-------------------------

Rank a set of molecules by their structural novelty relative to each other
or relative to a reference set:

.. code-block:: python

   from HalogenGroups import prioritise_molecules

   molecules = ["CCCC(F)(F)F", "FC(F)(F)C(=O)O", "FCCC(F)(F)F", "ClCCCl"]

   ranking = prioritise_molecules(molecules)
   for smiles, score in ranking:
       print(smiles, score)

See :doc:`prioritization` for parameter details.

Loading and inspecting groups
------------------------------

.. code-block:: python

   from HalogenGroups import get_compiled_HalogenGroups

   groups = get_compiled_HalogenGroups()
   print(len(groups))       # 116

   for g in groups[:5]:
       print(g.name, g.category, g.smarts)

Custom group collections
-------------------------

See :doc:`customization` for how to create and use custom halogen group
definitions.