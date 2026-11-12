Customization
=============

.. code-block:: python

   from PFASGroups import parse_smiles, get_compiled_HalogenGroups, get_HalogenGroups

   # Filter to OECD groups only
   # The category of each group is in the raw JSON, under "sorting"
   oecd_ids = {g["id"] for g in get_HalogenGroups()
               if (g.get("sorting") or {}).get("category") == "OECD"}
   oecd_groups = [g for g in get_compiled_HalogenGroups() if g.id in oecd_ids]
   results = parse_smiles(["CCCC(F)(F)F"], pfas_groups=oecd_groups)

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
   print(len(groups))           # 118

   for g in groups:
       print(g.id, g.name)

Filtering groups
-----------------

Pass a custom subset to ``parse_smiles`` via its ``halogen_groups`` parameter:

.. code-block:: python

   from PFASGroups import parse_smiles, get_compiled_HalogenGroups, get_HalogenGroups

   all_groups = get_compiled_HalogenGroups()

   # Keep only OECD groups
   # The category of each group is in the raw JSON, under "sorting"
   oecd_ids = {g["id"] for g in get_HalogenGroups()
               if (g.get("sorting") or {}).get("category") == "OECD"}
   oecd_only = [g for g in all_groups if g.id in oecd_ids]

   # Pass custom groups as pfas_groups (unknown keywords are ignored silently)
   results = parse_smiles(["CCCC(F)(F)F"], pfas_groups=oecd_only)

Defining a custom group
------------------------

A :class:`~PFASGroups.HalogenGroup` object requires at minimum a name,
category, SMARTS pattern, and group ID:

.. code-block:: python

   from PFASGroups import HalogenGroup, parse_smiles

   my_group = HalogenGroup(
       id=200,
       name="my_custom_group",
       smarts={"[CX4](F)(F)(F)[CX4](F)(F)": 1},   # CF3-CF2, at least once
       componentSmarts="Perfluoroalkyl",
       componentHalogens="F",
       constraints={},
   )

   results = parse_smiles(["FC(F)(F)C(F)(F)CC"], pfas_groups=[my_group])
   print([m.group_name for m in results[0].matches if m.is_group])

Combining built-in and custom groups
--------------------------------------

.. code-block:: python

   from PFASGroups import get_compiled_HalogenGroups, HalogenGroup, parse_smiles

   base_groups = get_compiled_HalogenGroups()
   my_group = HalogenGroup(
       id=201,
       name="vicinal dichloride",
       smarts={"[CX4]([Cl])[CX4][Cl]": 1},
       constraints={},
   )
   extended = list(base_groups) + [my_group]
   results = parse_smiles(["ClCCCl"], pfas_groups=extended, halogens="Cl")

HalogenGroup attributes
------------------------

.. list-table::
   :header-rows: 1
   :widths: 25 75

   * - Attribute
     - Description
   * - ``id``
     - Integer identifier (must be unique in the collection you pass)
   * - ``name``
     - Human-readable group name
   * - ``smarts``
     - Dict ``{SMARTS: minimum count}`` of functional-group patterns
   * - ``componentSmarts``
     - Halogenated component type the group is attached to (e.g.
       ``'Perfluoroalkyl'``), or ``None``
   * - ``componentHalogens`` / ``componentSaturation`` / ``componentForm``
     - Restrict the component to halogens (``'F'``, ``['F', 'Cl']``, ...),
       saturation (``'per'``, ``'poly'``) and form (``'alkyl'``, ``'cyclic'``)
   * - ``constraints``
     - Molecular-formula constraints (``'only'``, ``'gte'``, ``'lte'``,
       ``'eq'``, ``'rel'``); see :class:`~PFASGroups.HalogenGroup`
   * - ``max_dist_from_comp`` / ``linker_smarts``
     - Allowed distance between component and functional group, and SMARTS
       the linker atoms must match
   * - ``compute``
     - If ``False`` the group is an aggregate (``re_search`` regex over group
       names, e.g. ``telomers``) and is not matched directly

Using raw (uncompiled) groups
-------------------------------

:func:`~PFASGroups.get_HalogenGroups` returns the raw (uncompiled)
group definitions as plain data objects:

.. code-block:: python

   from PFASGroups import get_HalogenGroups
   raw = get_HalogenGroups()

   for g in raw:
       print(g["id"], g["name"], g.get("smarts"))   # aggregate groups have no SMARTS

:func:`~PFASGroups.get_compiled_HalogenGroups` returns them compiled.  To
compile an edited copy of the raw definitions, build the
:class:`~PFASGroups.HalogenGroup` instances yourself:

.. code-block:: python

   from PFASGroups import HalogenGroup, get_HalogenGroups, parse_smiles
   raw = get_HalogenGroups()
   compiled = [HalogenGroup(**g) for g in raw if g.get("compute", True)]
   results = parse_smiles(["FC(F)(F)C(F)(F)C(=O)O"], pfas_groups=compiled)

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
a component labelled *polyhalogenated* bears at least two halogen substituents,
distinguishing it from mono-halogenated components.