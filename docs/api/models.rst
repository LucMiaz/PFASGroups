Data Models
===========

.. currentmodule:: PFASGroups

.. code-block:: python

   from PFASGroups import parse_smiles

   results = parse_smiles(["CCCC(F)(F)F", "OCCOCCO"])
   mol = results[0]                      # MoleculeResult

   print(mol.smiles, mol.is_PFAS)
   for match in mol.matches:             # GroupMatch objects
       print(match.group_name, match.group_id)
       for comp in match.components:     # MatchComponent objects
           print("  atoms:", comp.atoms)

   arr = results.to_array()             # EmbeddingArray
   print(arr.shape)                      # (2, 116)

.. contents:: Contents
   :local:
   :depth: 2

PFASEmbeddingSet
------------

.. autoclass:: PFASEmbeddingSet
   :members:
   :undoc-members:
   :show-inheritance:

``PFASEmbeddingSet`` is a list-like container of :class:`MoleculeResult` objects.
Its length equals the number of input SMILES.

**Key methods:**

.. list-table::
   :header-rows: 1
   :widths: 40 60

   * - Method
     - Description
   * - ``results[i]``
     - Access the i-th :class:`MoleculeResult`
   * - ``results.to_dataframe()``
     - Flatten all matches to a ``pandas.DataFrame``
   * - ``results.to_array(group_selection='all', component_metrics=['binary'], halogens='F', saturation='per')``
     - Convert to an :class:`EmbeddingArray` (**PFASGroups** default: ``halogens='F'``, 116 cols;
       **HalogenGroups** subclass default: ``halogens=['F','Cl','Br','I']``, 464 cols)
   * - ``results.to_sql(filename)``
     - Persist to a SQLite or PostgreSQL database
   * - ``PFASEmbeddingSet.from_sql(filename)``
     - Load from a previously saved database

to_array options:

.. note::

   The default value of ``halogens`` depends on which module the
   :class:`PFASEmbeddingSet` came from:

   * ``from PFASGroups import parse_smiles`` → ``halogens='F'`` (116 columns)
   * ``from HalogenGroups import parse_smiles`` → ``halogens=['F','Cl','Br','I']`` (464 columns)

   Always pass ``halogens`` explicitly in reusable helpers or notebook functions
   to avoid silent fingerprint-width changes when the import source changes.

.. code-block:: python

   # PFASGroups default: all 116 groups, F only, binary → (n, 116)
   arr = results.to_array()                                        # (n, 116)

   # Always-explicit (safe in any import context)
   arr = results.to_array(halogens='F')                            # (n, 116)

   # OECD groups only
   arr = results.to_array(group_selection='oecd', halogens='F')   # (n, 28)

   # Count or max-component encoding
   arr = results.to_array(component_metrics=['count'], halogens='F')
   arr = results.to_array(component_metrics=['max_component'], halogens='F')

   # Multi-halogen (advanced) — see halogengroups page
   arr = results.to_array(halogens=['F', 'Cl', 'Br', 'I'])       # (n, 464)

MoleculeResult
--------------

.. autoclass:: MoleculeResult
   :members:
   :undoc-members:
   :show-inheritance:

Represents parsing results for a single molecule.

**Attributes:**

.. list-table::
   :header-rows: 1
   :widths: 25 75

   * - Attribute
     - Description
   * - ``smiles``
     - Canonical SMILES string
   * - ``inchi``
     - InChI string
   * - ``inchikey``
     - InChIKey
   * - ``matches``
     - List of :class:`GroupMatch` objects
   * - ``pfas_definition_matches``
     - List of definition matches (populated when
       ``include_PFAS_definitions=True``)
   * - ``n_matches``
     - Number of group matches
   * - ``is_PFAS``
     - ``True`` if any match has ``is_PFAS=True``

GroupMatch
----------

Represents a single group detected in a molecule.

**Attributes:**

.. list-table::
   :header-rows: 1
   :widths: 25 75

   * - Attribute
     - Description
   * - ``group_name``
     - Human-readable group name
   * - ``group_id``
     - Integer group ID
   * - ``group_category``
     - ``'OECD'``, ``'Generic'``, or ``'Fluorotelomer'``
   * - ``is_PFAS``
     - Whether this group qualifies as PFAS
   * - ``halogen``
     - Halogen symbol matched (``'F'``, ``'Cl'``, etc.)
   * - ``components``
     - List of :class:`MatchComponent` objects
   * - ``n_components``
     - Number of components

MatchComponent
--------------

A single structural component of a group match.

**Attributes:**

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - Attribute
     - Description
   * - ``atoms``
     - List of atom indices (0-based) in the RDKit molecule
   * - ``n_atoms``
     - Number of atoms in the component
   * - ``n_halogens``
     - Number of halogen atoms
   * - ``halogen_fraction``
     - Ratio of halogen atoms to total heavy atoms
   * - ``effective_graph_resistance``
     - Kirchhoff index of the component graph (``None`` if not computed)

EmbeddingArray
------------------

See :doc:`fingerprint_analysis` for full documentation.

:class:`EmbeddingArray` is a numpy array subclass returned by
:meth:`PFASEmbeddingSet.to_array`. It carries molecule identity metadata:

.. code-block:: python

   arr = results.to_array()
   print(arr.shape)        # (n_mols, n_groups)
   print(arr.smiles)       # list of input SMILES strings

HalogenGroup
------------

.. autoclass:: HalogenGroup
   :members:
   :undoc-members:
   :show-inheritance:

Defines a single halogen structural group (SMARTS pattern + metadata).
Used to build custom group libraries.  See :doc:`../customization`.

PFASDefinition
--------------

.. autoclass:: PFASDefinition
   :members:
   :undoc-members:
   :show-inheritance:

Encapsulates a regulatory PFAS definition and its matching logic.
See :doc:`../pfas_definitions` for descriptions of the five built-in definitions.