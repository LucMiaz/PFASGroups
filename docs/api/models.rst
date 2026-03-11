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

   fp = results.to_fingerprint()         # ResultsFingerprint
   print(fp.fingerprints.shape)          # (2, 116)

.. contents:: Contents
   :local:
   :depth: 2

ResultsModel
------------

.. autoclass:: ResultsModel
   :members:
   :undoc-members:
   :show-inheritance:

``ResultsModel`` is a list-like container of :class:`MoleculeResult` objects.
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
   * - ``results.to_fingerprint(group_selection='all', component_metrics=['binary'], halogens='F', saturation='per')``
     - Convert to a :class:`ResultsFingerprint` (**PFASGroups** default: ``halogens='F'``, 116 cols;
       **HalogenGroups** subclass default: ``halogens=['F','Cl','Br','I']``, 464 cols)
   * - ``results.to_sql(filename)``
     - Persist to a SQLite or PostgreSQL database
   * - ``ResultsModel.from_sql(filename)``
     - Load from a previously saved database

to_fingerprint options:

.. note::

   The default value of ``halogens`` depends on which module the
   :class:`ResultsModel` came from:

   * ``from PFASGroups import parse_smiles`` → ``halogens='F'`` (116 columns)
   * ``from HalogenGroups import parse_smiles`` → ``halogens=['F','Cl','Br','I']`` (464 columns)

   Always pass ``halogens`` explicitly in reusable helpers or notebook functions
   to avoid silent fingerprint-width changes when the import source changes.

.. code-block:: python

   # PFASGroups default: all 116 groups, F only, binary → (n, 116)
   fp = results.to_fingerprint()                                   # (n, 116)

   # Always-explicit (safe in any import context)
   fp = results.to_fingerprint(halogens='F')                       # (n, 116)

   # OECD groups only
   fp = results.to_fingerprint(group_selection='oecd', halogens='F')  # (n, 28)

   # Count or max-component encoding
   fp = results.to_fingerprint(component_metrics=['count'], halogens='F')
   fp = results.to_fingerprint(component_metrics=['max_component'], halogens='F')

   # Multi-halogen (advanced) — see halogengroups page
   fp = results.to_fingerprint(halogens=['F', 'Cl', 'Br', 'I'])  # (n, 464)

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

ResultsFingerprint
------------------

See :doc:`fingerprint_analysis` for full documentation.

The fingerprint matrix is stored as ``fp.fingerprints`` (``numpy.ndarray``):

.. code-block:: python

   fp = results.to_fingerprint()
   print(fp.fingerprints.shape)   # (n_mols, n_groups)
   print(fp.group_names)          # list of group-name strings

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