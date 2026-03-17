Quickstart
==========

Five-minute overview of the most common PFASGroups workflows.

.. code-block:: python

   from PFASGroups import parse_smiles

   results = parse_smiles(["CCCC(F)(F)F", "FC(F)(F)C(=O)O"])
   arr, cols = results.to_array(), results.column_names()
   print(arr.shape)   # (2, n_groups)

.. contents:: Contents
   :local:
   :depth: 1

Parsing SMILES
--------------

.. code-block:: python

   from PFASGroups import parse_smiles

   smiles = [
       "CCCC(F)(F)F",       # perfluoroalkyl chain
       "FC(F)(F)C(=O)O",    # trifluoroacetic acid (TFA)
       "OCCOCCO",           # no halogen — returns no matches
   ]

   results = parse_smiles(smiles)

``results`` is a :class:`~PFASGroups.PFASEmbeddingSet` — a list-like container
of :class:`~PFASGroups.PFASEmbedding` objects (dict subclass), one per input SMILES.

Accessing matches
-----------------

.. code-block:: python

   mol = results[0]                             # first molecule (PFASEmbedding)
   print(mol.smiles)                            # canonical SMILES
   print(bool(mol.matches))                     # True if any group matched

   for match in mol.matches:                    # iterate over MatchView objects
       if match.is_group:
           print(match.group_name)              # e.g. "Perfluoroalkyl"
           print(match.group_id)               # integer group ID
           for comp in match.components:
               print(comp.atoms)               # list of atom indices

Loop over only molecules that have at least one match:

.. code-block:: python

   for mol in results:
       if mol.matches:
           print(mol.smiles, "—", len(mol.matches), "match(es)")

Converting to a DataFrame
-------------------------

.. code-block:: python

   df = results.to_dataframe()
   print(df.columns.tolist())
   # ['smiles', 'inchi', 'group_name', 'group_id', ...]

Generating embeddings
---------------------

Embeddings encode group matches as a fixed-length numeric vector suitable for
machine learning.  By default PFASGroups produces a **binary vector** with one
column per group (fluorine only):

.. code-block:: python

   from PFASGroups import parse_smiles

   smiles = ["CCCC(F)(F)F", "FC(F)(F)C(=O)O", "OCCOCCO"]

   # Convenience function — parses and returns (array, column_names)
   results = parse_smiles(smiles)
   arr, cols = results.to_array(), results.column_names()
   print(arr.shape)   # (3, n_groups) — one row per molecule
   print(type(arr))   # numpy.ndarray
   print(cols[:2])    # ['Perfluoromethyl [binary]', 'Perfluoroalkyl [binary]', ...]

   # From a pre-parsed set
   results = parse_smiles(smiles)
   arr  = results.to_array()    # (3, n_groups) matrix
   cols = results.column_names()  # matching column labels

**Group selection** — restrict to a subset of groups:

.. code-block:: python

   # OECD groups only
   arr_oecd, cols_oecd = results.to_array(group_selection='oecd'), results.column_names(group_selection='oecd')

   # From a pre-parsed set
   arr_oecd = results.to_array(group_selection='oecd')

**component_metrics** — control how matches are encoded:

.. code-block:: python

   # binary (default): 1 = present, 0 = absent
   arr_bin, _ = results.to_array(component_metrics=['binary']), results.column_names(component_metrics=['binary'])

   # count: number of independent matches
   arr_cnt, _ = results.to_array(component_metrics=['count']), results.column_names(component_metrics=['count'])

   # max_component: size (atom count) of the largest matching component
   arr_max, _ = results.to_array(component_metrics=['max_component']), results.column_names(component_metrics=['max_component'])

   # Preset combining binary + effective graph resistance ('best')
   arr_best, _ = results.to_array(preset='best'), results.column_names(preset='best')

.. note::

   For multi-halogen embeddings covering F, Cl, Br and I, see
   :ref:`multi-halogen fingerprinting <multi_halogen_fingerprint>` in
   :doc:`halogengroups`.

PFAS definition screening
--------------------------

.. code-block:: python

   from PFASGroups import parse_smiles

   results = parse_smiles(
       ["CCCC(F)(F)F", "OCCOCCO"],
       include_PFAS_definitions=True,
   )

   for mol in results:
       for match in mol.matches:
           if match.is_definition:
               print(mol.smiles, "matches", match.get("definition_name"))

Saturation filter
-----------------

.. code-block:: python

   # Only perfluorinated (fully saturated C–F) groups
   results = parse_smiles(smiles, saturation='per')

   # Polyfluorinated groups (partially substituted)
   results = parse_smiles(smiles, saturation='poly')

   # No filter — all groups (default: saturation=None for parse_smiles)
   results = parse_smiles(smiles, saturation=None)

Multi-halogen parsing (advanced)
---------------------------------

To detect Cl, Br and I groups in addition to fluorine, use the
``halogens`` argument or import from ``HalogenGroups``:

.. code-block:: python

   from PFASGroups import parse_smiles

   results = parse_smiles(["ClCCCl", "BrCCBr"], halogens=['F', 'Cl', 'Br', 'I'])

See :doc:`halogengroups` for full multi-halogen documentation.

Command-line usage
------------------

.. code-block:: bash

   # Parse a CSV of SMILES
   pfasgroups parse input.csv --output results.json

   # Generate fingerprints
   pfasgroups fingerprint input.csv --output fps.csv

   # List all 116 group names
   pfasgroups list-groups

See :doc:`cli` for the full CLI reference.