Data Models
===========

.. currentmodule:: HalogenGroups

This page documents the data model classes returned by
:func:`~HalogenGroups.parse_smiles` and related functions.

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
   :widths: 35 65

   * - Method
     - Description
   * - ``results[i]``
     - Access the i-th :class:`MoleculeResult`
   * - ``results.to_dataframe()``
     - Flatten all matches to a ``pandas.DataFrame``
   * - ``results.to_fingerprint(...)``
     - Convert to a :class:`ResultsFingerprint`
   * - ``results.to_sql(filename)``
     - Persist to a SQLite or PostgreSQL database
   * - ``ResultsModel.from_sql(filename)``
     - Load from a previously saved database

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
     - List of :class:`PFASDefinitionMatch` objects (populated when
       ``include_PFAS_definitions=True``)
   * - ``n_matches``
     - Number of group matches
   * - ``is_PFAS``
     - ``True`` if any match has ``is_PFAS=True``

**Example:**

.. code-block:: python

   mol = results[0]
   print(mol.smiles)
   print(mol.is_PFAS)
   for match in mol.matches:
       print(match.group_name)

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

A single structural component of a group match (multiple matches of the same
group in one molecule appear as separate components).

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

HalogenGroup
------------

.. autoclass:: HalogenGroup
   :members:
   :undoc-members:
   :show-inheritance:

Defines a single halogen structural group (SMARTS pattern + metadata).
Used to build custom group libraries.

**Constructor parameters:**

.. list-table::
   :header-rows: 1
   :widths: 25 75

   * - Parameter
     - Description
   * - ``group_id``
     - Unique integer ID
   * - ``name``
     - Group name
   * - ``category``
     - ``'OECD'``, ``'Generic'``, ``'Fluorotelomer'``, or custom string
   * - ``smarts``
     - SMARTS string or list of SMARTS strings
   * - ``is_PFAS``
     - Boolean
   * - ``compute``
     - If ``False``, the group is listed in fingerprint headers but not matched
   * - ``description``
     - Optional free-text description

PFASDefinition
--------------

.. autoclass:: PFASDefinition
   :members:
   :undoc-members:
   :show-inheritance:

Encapsulates a regulatory PFAS definition and its matching logic.

See :doc:`../pfas_definitions` for descriptions of the five built-in
definitions.