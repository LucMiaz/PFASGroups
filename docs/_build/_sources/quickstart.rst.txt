Quickstart
==========

Five-minute overview of the most common PFASGroups workflows.

.. code-block:: python

   from PFASGroups import parse_smiles, generate_fingerprint

   results = parse_smiles(["CCCC(F)(F)F", "FC(F)(F)C(=O)O"])
   fps, info = generate_fingerprint(["CCCC(F)(F)F", "FC(F)(F)C(=O)O"])
   print(fps.shape)   # (2, 116)

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

``results`` is a :class:`~PFASGroups.ResultsModel` — a list-like container
of :class:`~PFASGroups.MoleculeResult` objects, one per input SMILES.

Accessing matches
-----------------

.. code-block:: python

   mol = results[0]                        # first molecule
   print(mol.smiles)                       # canonical SMILES
   print(mol.is_PFAS)                      # True / False

   for match in mol.matches:
       print(match.group_name)             # e.g. "perfluoroalkyl"
       print(match.group_id)              # integer group ID
       print(match.is_PFAS)
       for comp in match.components:
           print(comp.atoms)              # list of atom indices

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
   # ['smiles', 'inchi', 'group_name', 'group_id', 'is_PFAS', ...]

Generating fingerprints
-----------------------

Fingerprints encode group matches as a fixed-length vector suitable for
machine learning.  By default PFASGroups produces a **116-column binary
vector** (one column per group, fluorine only):

.. code-block:: python

   from PFASGroups import generate_fingerprint

   smiles = ["CCCC(F)(F)F", "FC(F)(F)C(=O)O", "OCCOCCO"]
   fps, info = generate_fingerprint(smiles)

   print(fps.shape)                  # (3, 116) — 116 groups, F only
   print(type(fps))                  # numpy.ndarray
   print(info['group_names'][:3])    # ['perfluoromethyl', 'perfluoroalkyl', ...]

**Group selection** — restrict to a subset of groups:

.. code-block:: python

   # OECD groups only (IDs 1–28)
   fps_oecd, info = generate_fingerprint(smiles, selected_groups=range(0, 28))

**component_metrics** — control how matches are encoded:

.. code-block:: python

   # binary (default): 1 = present, 0 = absent
   fps_bin, _ = generate_fingerprint(smiles, component_metrics=['binary'])

   # count: number of independent matches
   fps_cnt, _ = generate_fingerprint(smiles, component_metrics=['count'])

   # max_component: size (atom count) of the largest matching component
   fps_max, _ = generate_fingerprint(smiles, component_metrics=['max_component'])

   # Add graph metric columns (binary + EGR, i.e. preset='best')
   fps_best, _ = generate_fingerprint(smiles, component_metrics=['binary', 'effective_graph_resistance'])

.. note::

   For multi-halogen fingerprints covering F, Cl, Br and I (4 × 116 = 464
   columns), see :ref:`multi-halogen fingerprinting <multi_halogen_fingerprint>` in
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
       for defn in mol.pfas_definition_matches:
           print(mol.smiles, "matches", defn.definition_name)

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