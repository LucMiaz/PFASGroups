Multi-Halogen Analysis (Advanced)
==================================

.. note::

   This is an **advanced topic**.  For standard PFAS / fluorine-only work,
   use :doc:`quickstart` which focuses on the ``PFASGroups`` import.
   This page describes how to extend detection to **Cl, Br and I** groups.

**PFASGroups** focuses on fluorinated substances (``halogens='F'`` by default).
The same codebase can also detect chlorinated, brominated and iodinated
structural groups when you pass ``halogens`` explicitly, or by importing from
``HalogenGroups`` which sets all four halogens as the default.

.. code-block:: python

   # Option A — explicit halogens argument with PFASGroups (recommended)
   from PFASGroups import parse_smiles
   results = parse_smiles(["ClC(Cl)(Cl)Cl", "BrCCBr"], halogens=['F','Cl','Br','I'])

   # Option B — import HalogenGroups, all-halogens is the default
   from HalogenGroups import parse_smiles
   results = parse_smiles(["ClC(Cl)(Cl)Cl", "BrCCBr"])


.. list-table::
   :header-rows: 1
   :widths: 25 30 45

   * - Import
     - Default ``halogens``
     - Typical use
   * - ``HalogenGroups``
     - ``['F', 'Cl', 'Br', 'I']``
     - Broad multi-halogen screening
   * - ``PFASGroups``
     - ``'F'``
     - PFAS / fluorine-focused analysis

Either import can be forced to any halogen set by passing the ``halogens``
argument explicitly — the results will be identical.


Quick Comparison
----------------

.. code-block:: python

   from HalogenGroups import parse_smiles as hal_parse
   from PFASGroups   import parse_smiles as pfas_parse

   smiles = ["ClC(Cl)(Cl)C(Cl)(Cl)Cl",   # perchlorinated
             "BrC(Br)(Br)CBr",             # brominated
             "FC(F)(F)C(F)(F)C(=O)O"]      # PFBA

   # HalogenGroups: all halogens included by default
   results_all = hal_parse(smiles)

   # PFASGroups: fluorine only by default → Cl/Br compounds show 0 matches
   results_f   = pfas_parse(smiles)

   # Equivalent — pass halogens explicitly to PFASGroups
   results_eq  = pfas_parse(smiles, halogens=['F', 'Cl', 'Br', 'I'])


.. warning::

   **The default embedding width depends on which module you import from.**

   ``parse_smiles`` from ``HalogenGroups`` matches all four halogens at parse
   time, but ``results.to_array()`` **without** ``halogens`` still returns one
   column per group (118), as in ``PFASGroups``: a group matched through Cl
   and one matched through F set the same column.  For one block of columns
   per halogen, pass the ``halogens`` list explicitly::

      arr = results.to_array(halogens='F')                     # 117 columns
      arr = results.to_array(halogens=['F', 'Cl', 'Br', 'I'])  # 468 columns

   This applies to ``compute_config``, ``_txp_fingerprint``, and any similar
   helpers in notebooks or scripts that pass ``**kwargs`` to
   ``to_array()`` — the ``halogens`` key must be present in those kwargs.

Functions with Altered Defaults
---------------------------------

The following functions have their ``halogens`` default overridden to
``['F', 'Cl', 'Br', 'I']`` when imported from ``HalogenGroups``:

.. list-table::
   :header-rows: 1
   :widths: 40 60

   * - Function
     - Default change
   * - :func:`parse_smiles`
     - ``halogens=['F','Cl','Br','I']``
   * - :func:`parse_mols`
     - ``halogens=['F','Cl','Br','I']``
   * - :func:`generate_fingerprint`
     - ``halogens=['F','Cl','Br','I']``
   * - :class:`PFASEmbeddingSet`\ ``.to_array()``
     - ``halogens=['F','Cl','Br','I']``, ``saturation='per'``

All other functions (``parse_mol``, ``parse_groups_in_mol``, ``get_HalogenGroups``,
``prioritise_molecules``, … ) are re-exported unchanged.


Parsing Multi-Halogen Molecules
---------------------------------

.. code-block:: python

   from HalogenGroups import parse_smiles

   smiles = [
       "FC(F)(F)C(F)(F)C(=O)O",           # PFBA — fluorinated
       "ClC(Cl)(Cl)C(Cl)(Cl)C(=O)O",      # perchlorinated carboxylic acid
       "BrC(Br)(Br)CBr",                   # brominated
   ]

   # All four halogens considered — no extra argument needed
   results = parse_smiles(smiles)

   for mol in results:
       print(mol.smiles)
       for match in mol.matches:
           comps = match.components
           print(f"  {match.group_name}: {len(comps)} component(s)")

To restrict to a subset of halogens with ``HalogenGroups``:

.. code-block:: python

   # Override the default to fluorine + chlorine only
   results = parse_smiles(smiles, halogens=['F', 'Cl'])


.. _multi_halogen_fingerprint:

Multi-Halogen Fingerprints
---------------------------

:meth:`PFASEmbeddingSet.to_array` with a list of halogens **stacks** one block
of 117 columns per halogen horizontally: 4 halogens give **117 × 4 = 468
columns**.  Column names end with ``(F)``, ``(Cl)``, ``(Br)``, ``(I)``.
:func:`generate_fingerprint` currently returns one column per group (118)
whatever ``halogens`` is, so use ``to_array`` for per-halogen blocks.

.. code-block:: python

   from HalogenGroups import parse_smiles

   smiles = ["FC(F)(F)C(F)(F)C(=O)O",
             "ClC(Cl)(Cl)C(Cl)(Cl)C(=O)O"]

   results = parse_smiles(smiles)
   fps = results.to_array(halogens=['F', 'Cl', 'Br', 'I'])
   print(fps.shape)                                          # (2, 468)
   print(results.column_names(halogens=['F', 'Cl'])[117])    # 'OECD perfluoroalkyl alcohols [binary] (Cl)'



Via PFASEmbeddingSet
~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   from HalogenGroups import parse_smiles

   results = parse_smiles(smiles)

   arr_all  = results.to_array()                                   # shape (n, 118)
   arr_f    = results.to_array(halogens='F')                      # shape (n, 117)
   arr_fc   = results.to_array(halogens=['F', 'Cl'])              # shape (n, 234)
   arr_oecd = results.to_array(
       group_selection='oecd', halogens=['F', 'Cl', 'Br', 'I'])   # shape (n, 112)

   # Best preset (binary + effective_graph_resistance) with F only
   arr_best = results.to_array(preset='best', halogens='F')       # shape (n, 234)

   # Explicit component_metrics — count mode + graph metric with all halogens
   arr_cm = results.to_array(
       component_metrics=['binary', 'effective_graph_resistance'],
       halogens=['F', 'Cl', 'Br', 'I'])  # shape (n, 936)  — 2 × 117 × 4

Each entry in ``component_metrics`` adds one block of *n_groups × n_halogens*
columns. The column naming follows the pattern ``"GroupName [metric] (halogen)"``
(e.g. ``"OECD perfluoroalkyl alcohols [binary] (F)"``).


Combining with Saturation Filters
-----------------------------------

The ``saturation`` parameter (``'per'``, ``'poly'``, or ``None``) applies to all
halogens simultaneously and controls which component SMARTS are used for groups
that have halogenated-chain components (OECD groups 1-28).

.. code-block:: python

   from HalogenGroups import parse_smiles

   # Perhalogenated components only (default when using HalogenGroups fingerprinting)
   r_per  = parse_smiles(smiles, saturation='per')

   # Polyhalogenated components only
   r_poly = parse_smiles(smiles, saturation='poly')

   # All components (per + poly)
   r_all  = parse_smiles(smiles, saturation=None)


CLI with Multiple Halogens
---------------------------

The CLI always requires an explicit ``--halogens`` flag:

.. code-block:: bash

   # All four halogens
   halogengroups parse --halogens F Cl Br I "ClC(Cl)(Cl)C(Cl)(Cl)C(=O)O"

   # Per-saturation, F + Cl, OECD groups
   halogengroups parse --halogens F Cl --saturation per "FC(F)(F)C(F)(F)C(=O)O"

   # Stacked fingerprint with F + Cl
   halogengroups fingerprint --halogens F Cl "FC(F)(F)C(F)(F)C(=O)O"

See :doc:`cli` for the full CLI reference.


Implementation Details
-----------------------

``HalogenGroups/__init__.py`` uses :func:`functools.wraps`-compatible wrappers
that inject ``halogens=['F', 'Cl', 'Br', 'I']`` as a default keyword argument.
Because Python keyword defaults can always be overridden at call time, every
explicit ``halogens=...`` argument takes precedence.

The ``PFASEmbeddingSet`` subclass in ``HalogenGroups`` overrides only
``to_array()``; all other methods (``show()``, ``summary()``, ``to_sql()``,
etc.) are inherited unchanged from ``PFASGroups.PFASEmbeddings.PFASEmbeddingSet``.


When to Use Which Import
-------------------------

Use ``HalogenGroups`` when:

* Your compounds include Cl-, Br- or I-containing structures
* You want to compare fluorination vs. chlorination patterns side by side
* You are building a generic halogenated-substance screening workflow

Use ``PFASGroups`` (``halogens='F'``) when:

* You are working exclusively with PFAS (fluorine-only)
* You need strict compatibility with published PFAS fingerprint benchmarks
* You want smaller fingerprints (117-column vs. 468-column)
* You mix fluorine-only and multi-halogen calls in the same script


See Also
--------

* :doc:`quickstart` — first steps with the package
* :doc:`api/core` — full ``parse_smiles`` / ``generate_fingerprint`` reference
* :doc:`api/models` — ``PFASEmbeddingSet`` and ``EmbeddingArray`` details
* :doc:`customization` — adding custom halogen groups
