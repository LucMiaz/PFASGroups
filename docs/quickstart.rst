Quickstart
==========

This page shows the most common use cases in under five minutes.

Parsing SMILES
--------------

.. code-block:: python

   from HalogenGroups import parse_smiles

   smiles = [
       "CCCC(F)(F)F",               # perfluoroalkyl chain
       "ClCCCl",                     # organochlorine
       "FC(F)(F)C(=O)O",            # trifluoroacetic acid (TFA)
       "OCCOCCO",                    # no halogen — returns no matches
   ]

   results = parse_smiles(smiles)

``results`` is a :class:`~HalogenGroups.ResultsModel` (a list-like container
of :class:`~HalogenGroups.MoleculeResult` objects, one per input SMILES).

Accessing matches
-----------------

.. code-block:: python

   mol = results[0]                         # first molecule
   print(mol.smiles)                        # canonical SMILES
   print(mol.inchi)

   for match in mol.matches:
       print(match.group_name)              # e.g. "perfluoroalkyl"
       print(match.group_id)               # integer group ID
       print(match.is_PFAS)
       for comp in match.components:
           print(comp.atoms)               # list of atom indices

Skipping molecules with no match:

.. code-block:: python

   for mol in results:
       if mol.matches:
           print(mol.smiles, "—", len(mol.matches), "match(es)")

Converting to a DataFrame
-------------------------

.. code-block:: python

   df = results.to_dataframe()
   print(df.columns.tolist())
   # ['smiles', 'inchi', 'group_name', 'group_id', 'is_PFAS', 'n_components', ...]

Fluorine-only (PFASGroups mode)
--------------------------------

Use the ``PFASGroups`` import instead of ``HalogenGroups`` to restrict
detection to fluorine only (this is the default for backward compatibility):

.. code-block:: python

   from PFASGroups import parse_smiles

   results = parse_smiles(["ClCCCl", "CCCC(F)(F)F"])
   # ClCCCl has no fluorine → no matches
   # CCCC(F)(F)F → perfluoroalkyl match

Or pass ``halogens`` explicitly with either import:

.. code-block:: python

   from HalogenGroups import parse_smiles

   # fluorine only
   results_f = parse_smiles(smiles, halogens='F')

   # chlorine and bromine only
   results_clbr = parse_smiles(smiles, halogens=['Cl', 'Br'])

Generating fingerprints
-----------------------

Fingerprints encode group counts as a fixed-length vector suitable for
machine learning:

.. code-block:: python

   from HalogenGroups import generate_fingerprint

   smiles = ["CCCC(F)(F)F", "ClCCCl", "FC(F)(F)C(=O)O"]
   fps, group_names = generate_fingerprint(smiles)

   print(fps.shape)            # (3, 464) — 116 groups × 4 halogens
   print(type(fps))            # numpy.ndarray

Each column corresponds to (group, halogen) pair.  ``group_names`` is a
dict mapping column index to ``(group_name, halogen)`` tuples.

PFAS definition screening
--------------------------

Classify molecules against regulatory PFAS definitions:

.. code-block:: python

   from HalogenGroups import parse_smiles

   results = parse_smiles(
       ["CCCC(F)(F)F", "OCCOCCO"],
       include_PFAS_definitions=True,
   )

   for mol in results:
       for defn in mol.pfas_definition_matches:
           print(mol.smiles, "matches", defn.definition_name)

Saturated vs unsaturated groups
---------------------------------

Filter for fully saturated groups only:

.. code-block:: python

   results = parse_smiles(smiles, saturation='saturated')
   # or 'unsaturated' / None (default — all)

Command-line usage
------------------

.. code-block:: bash

   # Parse a CSV of SMILES and save results to JSON
   halogengroups parse input.csv --output results.json

   # Generate fingerprints
   halogengroups fingerprint input.csv --output fps.csv

   # List all 116 group names
   halogengroups list-groups

See :doc:`cli` for the full CLI reference.