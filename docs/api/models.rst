Data Models
===========

.. currentmodule:: PFASGroups

.. code-block:: python

   from PFASGroups import parse_smiles

   results = parse_smiles(["CCCC(F)(F)F", "OCCOCCO"])
   mol = results[0]                      # MoleculeResult (= PFASEmbedding)

   print(mol.smiles, mol["inchikey"])
   for match in mol.matches:             # MatchView objects
       if match.is_definition:           # PFAS definition matched
           print("definition:", match["definition_name"])
           continue
       print(match.group_name, match.group_id)
       for comp in match.components:     # ComponentView objects
           print("  atoms:", comp.atoms)

   arr = results.to_array()             # EmbeddingArray
   print(arr.shape)                      # (2, 118)

.. contents:: Contents
   :local:
   :depth: 2

PFASEmbeddingSet
----------------

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
   * - ``results.classify()``
     - ``pandas.DataFrame`` with one row per molecule (SMILES, category,
       total component size)
   * - ``results.table()`` / ``results.summary()``
     - Text table (one row per molecule) / printed coloured summary
   * - ``results.to_array(group_selection='all', component_metrics=['binary'], halogens=...)``
     - Convert to an :class:`EmbeddingArray`: 118 columns by default; with an
       explicit ``halogens`` list, one block of 117 columns per halogen
   * - ``results.to_sql(filename)``
     - Persist to a SQLite or PostgreSQL database
   * - ``PFASEmbeddingSet.from_sql(filename)``
     - Load from a previously saved database

to_array options:

.. note::

   Without ``halogens``, ``to_array()`` returns one column per group (118),
   whichever module the :class:`PFASEmbeddingSet` came from.  With
   ``halogens``, it returns one block of 117 columns per halogen, named
   ``"<group> [<metric>] (<halogen>)"``.  Pass ``halogens`` explicitly in
   reusable helpers so the width does not depend on the defaults.

.. code-block:: python

   # Default: all groups, binary → (n, 118)
   arr = results.to_array()                                        # (n, 118)

   # Explicit fluorine block
   arr = results.to_array(halogens='F')                            # (n, 117)

   # OECD groups only
   arr = results.to_array(group_selection='oecd', halogens='F')   # (n, 28)

   # Count or max-component encoding
   arr = results.to_array(component_metrics=['count'], halogens='F')
   arr = results.to_array(component_metrics=['max_component'], halogens='F')

   # Multi-halogen (advanced) — see halogengroups page
   arr = results.to_array(halogens=['F', 'Cl', 'Br', 'I'])       # (n, 468)

MoleculeResult
--------------

.. autoclass:: MoleculeResult
   :members:
   :undoc-members:
   :show-inheritance:

Represents parsing results for a single molecule.  ``MoleculeResult`` is an
alias of :class:`PFASEmbedding`, a ``dict`` subclass: fields are read as
``mol["inchikey"]``; ``smiles`` and ``matches`` are also attributes.

**Keys:**

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
   * - ``formula``
     - Molecular formula
   * - ``matches``
     - List of match objects: group matches (``type`` ``'HalogenGroup'`` or
       ``'WildcardGroup'``) and PFAS definition matches (``'PFASdefinition'``,
       see :doc:`../pfas_definitions`)

GroupMatch
----------

Represents a single match in ``mol.matches`` (class ``MatchView``, a ``dict``
subclass).  ``is_group`` / ``is_definition`` tell group and definition
matches apart; ``group_id``, ``group_name`` and ``components`` are attributes,
the other fields are keys (``m["halogen"]``).

**Fields:**

.. list-table::
   :header-rows: 1
   :widths: 25 75

   * - Attribute
     - Description
   * - ``group_name``
     - Human-readable group name
   * - ``group_id``
     - Integer group ID
   * - ``match_id``
     - ``'G<id>'`` (group), ``'W<id>'`` (wildcard group) or ``'D<id>'``
       (definition)
   * - ``type``
     - ``'HalogenGroup'``, ``'WildcardGroup'`` or ``'PFASdefinition'``
   * - ``halogen``
     - Halogen symbol matched (``'F'``, ``'Cl'``, etc.)
   * - ``components``
     - List of :class:`MatchComponent` objects
   * - ``num_components``
     - Number of components

MatchComponent
--------------

A single structural component of a group match (class ``ComponentView``).
The properties below are read as attributes; every other per-component metric
is in ``comp.data`` (e.g. ``comp.data["component_f_count"]``).

**Attributes:**

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - Attribute
     - Description
   * - ``atoms``
     - List of atom indices (0-based) in the RDKit molecule
   * - ``size``
     - Number of carbon atoms in the component
   * - ``smarts_label``
     - Component type (e.g. ``'Perfluoroalkyl'``)
   * - ``component_fraction``
     - Fraction of the molecule's atoms in the component
   * - ``branching``
     - 1.0 = linear, 0.0 = fully branched
   * - ``mean_eccentricity`` / ``diameter`` / ``radius``
     - Graph metrics of the component
   * - ``effective_graph_resistance`` / ``effective_graph_resistance_BDE``
     - Kirchhoff index of the component graph, unweighted / BDE-weighted
       (``None`` if not computed)
   * - ``min_dist_to_barycentre`` / ``min_dist_to_centre`` / ``max_dist_to_periphery``
     - Graph distances from the functional group to structural landmarks

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
See :doc:`../pfas_definitions` for descriptions of the six built-in definitions.