Customization
=============

.. code-block:: python

   from PFASGroups import parse_smiles, get_compiled_HalogenGroups, HalogenGroup

   # Filter to OECD groups only
   oecd_groups = [g for g in get_compiled_HalogenGroups() if g.category == 'OECD']
   results = parse_smiles(["CCCC(F)(F)F"], halogen_groups=oecd_groups)

The library is built on a composable set of :class:`~PFASGroups.HalogenGroup`
objects.  You can extend, filter, or replace the default group library without
modifying the source code.

.. contents:: Contents
   :local:
   :depth: 2

Inspecting the built-in groups
--------------------------------

.. code-block:: python

   from PFASGroups import get_compiled_HalogenGroups

   groups = get_compiled_HalogenGroups()
   print(len(groups))           # 116

   for g in groups:
       print(g.group_id, g.name, g.category)

Filtering groups
-----------------

Pass a custom subset to ``parse_smiles`` via its ``halogen_groups`` parameter:

.. code-block:: python

   from PFASGroups import parse_smiles, get_compiled_HalogenGroups

   all_groups = get_compiled_HalogenGroups()

   # Keep only OECD groups
   oecd_only = [g for g in all_groups if g.category == 'OECD']

   results = parse_smiles(["CCCC(F)(F)F"], halogen_groups=oecd_only)

Defining a custom group
------------------------

A :class:`~PFASGroups.HalogenGroup` object requires at minimum a name,
category, SMARTS pattern, and group ID:

.. code-block:: python

   from PFASGroups import HalogenGroup, parse_smiles

   my_group = HalogenGroup(
       group_id=200,
       name="my_custom_group",
       category="Custom",
       smarts="[CX4](F)(F)(F)[CX4](F)(F)",   # two adjacent CF3/CF2 units
       is_PFAS=True,
       compute=True,
   )

   results = parse_smiles(["CCCC(F)(F)F"], halogen_groups=[my_group])
   print(results[0].matches)

Combining built-in and custom groups
--------------------------------------

.. code-block:: python

   from PFASGroups import get_compiled_HalogenGroups, HalogenGroup, parse_smiles

   base_groups = get_compiled_HalogenGroups()
   my_group = HalogenGroup(
       group_id=201,
       name="chlorinated_alkyl",
       category="Custom",
       smarts="[CX4]([Cl])[CX4][Cl]",
       is_PFAS=False,
       compute=True,
   )
   extended = list(base_groups) + [my_group]
   results = parse_smiles(["ClCCCl"], halogen_groups=extended)

HalogenGroup attributes
------------------------

.. list-table::
   :header-rows: 1
   :widths: 25 75

   * - Attribute
     - Description
   * - ``group_id``
     - Integer identifier (must be unique in the collection you pass)
   * - ``name``
     - Human-readable group name
   * - ``category``
     - Category label, e.g. ``'OECD'``, ``'Generic'``, ``'Fluorotelomer'``
   * - ``smarts``
     - SMARTS pattern string (or list of SMARTS strings)
   * - ``is_PFAS``
     - Whether the group contributes to a PFAS classification
   * - ``compute``
     - If ``False`` the group appears in fingerprint headers but is not matched
   * - ``description``
     - Optional free-text description

Using raw (uncompiled) groups
-------------------------------

:func:`~PFASGroups.get_HalogenGroups` returns the raw (uncompiled)
group definitions as plain data objects:

.. code-block:: python

   from PFASGroups import get_HalogenGroups
   raw = get_HalogenGroups()

   for g in raw:
       print(g["name"], g["smarts"])

To compile them into :class:`~PFASGroups.HalogenGroup` instances call
:func:`~PFASGroups.get_compiled_HalogenGroups` (which internally calls
:func:`~PFASGroups.load_HalogenGroups`):

.. code-block:: python

   from PFASGroups import load_HalogenGroups, get_HalogenGroups
   raw = get_HalogenGroups()
   compiled = load_HalogenGroups(raw)

Custom SMARTS tips
-------------------

* SMARTS must be valid RDKit SMARTS (validated by ``Chem.MolFromSmarts``).
* Multiple SMARTS for one group can be supplied as a list:
  ``smarts=["[CX4](F)(F)F", "[CX4](F)(F)"]``
* Recursive SMARTS are supported.
* Test your SMARTS with ``rdkit.Chem.rdchem.Mol.GetSubstructMatches`` before
  adding them to the library.

Component-type constraints
---------------------------

Each entry in ``data/component_smarts_halogens.json`` may carry an optional
``constraints`` key that is evaluated at match time against the *full* atom
count of the component (backbone carbons **plus** directly attached halogens).

Supported constraint keys:

``gte`` (dict)
   Minimum element count.  The component must contain at least *n* atoms of the
   given element.  Example: ``{"F": 2}`` requires at least two fluorine atoms.

``exclude`` (list)
   Forbidden elements.  The component must contain none of the listed elements.
   Example: ``["Cl"]`` rejects any component that carries a chlorine atom.

Example JSON entry::

   "Polyfluoroalkyl": {
       "smarts": "...",
       "constraints": {"gte": {"F": 2}}
   }

These constraints are attached to the component type itself, so they apply
uniformly to every group that requires that component type — without needing to
repeat the check in each group definition in ``Halogen_groups_smarts.json``.

The built-in poly-halogenated alkyl component types (*Polyfluoroalkyl*,
*Polychloroalkyl*, *Polybromoalkyl*, *Polyiodoalkyl*) each carry
``{"gte": {"X": 2}}`` (where *X* is the respective halogen) to ensure that
a component labelled *poly*halogenated bears at least two halogen substituents,
distinguishing it from mono-halogenated components.