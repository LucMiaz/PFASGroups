Algorithm
=========

.. code-block:: python

   from PFASGroups import parse_smiles

   results = parse_smiles(["CCCC(F)(F)F", "FC(F)(F)C(=O)O"])
   # results[0].matches contains the detected halogen groups
   for mol in results:
       print(mol.smiles, [m.group_name for m in mol.matches])

This page describes how PFASGroups detects and classifies halogenated
structural groups in a molecule.

.. contents:: Contents
   :local:
   :depth: 2

Overview
--------

The algorithm processes a SMILES string in four stages:

1. **Parse** — convert SMILES to an RDKit molecule
2. **Match** — apply each HalogenGroup SMARTS pattern to the molecule
3. **Deduplicate** — resolve overlapping and redundant matches
4. **Score** — compute per-component structural metrics

Group library
-------------

The library contains **119 halogen groups** organised into four categories:

.. list-table::
   :header-rows: 1
   :widths: 15 15 70

   * - Category
     - Count
     - Description
   * - OECD
     - 27
     - Groups adopted from the 2021 OECD PFAS definition framework
   * - Generic
     - 48
     - Broader halogenated structural motifs (alkyl, aryl, acyl, sulfonyl, …)
   * - Fluorotelomer
     - 43
     - Fluorotelomer groups including telomer alcohols, sulfonates, amides
   * - Aggregate
     - 3
     - Pattern-matching groups (e.g. ``Telomers``) with ``compute=False``;
       included in fingerprint headers but not matched directly

When fluorine-only mode is used (the default), **114 groups** are compiled
(the 3 aggregate groups are always included in fingerprint headers but not
directly matched).

SMARTS-based matching
----------------------

Each group is defined by one or more SMARTS patterns.  The core matcher calls
:func:`~rdkit.Chem.rdchem.Mol.GetSubstructMatches` for each pattern.  If a
group has multiple SMARTS, the union of all matched atom-sets is collected.

Halogen filtering
~~~~~~~~~~~~~~~~~

After matching, each hit is filtered to keep only components that actually
contain the requested halogen(s).  For fluorine-only mode only C–F bonds are
retained; for multi-halogen mode any C–X bond where X ∈ target set is kept.

This filtering happens at the component level, so a group match can be
partially retained (some components kept, others discarded).

Saturation filtering
~~~~~~~~~~~~~~~~~~~~~

If ``saturation='saturated'``, only components whose carbon scaffold is fully
saturated (no sp2/sp3d carbons) are kept.  If ``saturation='unsaturated'``,
the complement is kept.  The default ``None`` retains all.

Overlap deduplication
-----------------------

The algorithm uses a priority-ordered merge to resolve overlapping SMARTS
matches:

1. Sort candidate groups by *specificity* (more specific first).
2. For each candidate group, remove any already-claimed atom from its
   matched components.
3. If a component shrinks below the minimum size threshold, discard it.
4. Record surviving components as confirmed matches.

This ensures that a carbon chain is attributed to the most specific group
it belongs to, and prevents double-counting.

Component graph metrics
------------------------

For each matched component a molecular graph is constructed and the following
metrics are computed (unless ``compute_component_metrics=False``):

**Effective graph resistance** (BDE-weighted Kirchhoff index)

Each bond :math:`(u,v)` with bond order :math:`b` is assigned a *conductance*
proportional to the bond's dissociation energy:

.. math::

   c_{uv} = \frac{\text{BDE}(Z_u,\, Z_v,\, b)}{\text{BDE}_\text{ref}}

where:

* :math:`\text{BDE}(Z_u, Z_v, 1)` — single-bond dissociation energy (kcal/mol)
  for the element pair :math:`(Z_u, Z_v)`, looked up from
  ``PFASGroups/data/diatomic_bonds_dict.json`` (the same reference table used
  by ``molecular_quantum_graph``).
* :math:`\text{BDE}(Z_u, Z_v, b) = \text{BDE}(Z_u, Z_v, 1) \cdot f(b)` —
  scaled by the **bond-order model** :math:`f(b)` described below.
* :math:`\text{BDE}_\text{ref}` — the C–C single-bond BDE (~83 kcal/mol),
  used as normalisation so that :math:`c_{CC,\text{single}} = 1`.

Higher BDE → higher conductance → *shorter* effective resistance path.
This means C–F bonds (BDE ≈ 130 kcal/mol, :math:`c \approx 1.56`) contribute
less resistance than C–C bonds (BDE ≈ 83 kcal/mol, :math:`c = 1.0`).

Bond-order model
~~~~~~~~~~~~~~~~

For bonds with order :math:`b > 1` (double, triple, aromatic), the single-bond
BDE is scaled by :math:`f(b)` where :math:`f(1) = 1` by construction.
PFASGroups attempts to load the calibrated model produced by
``molecular_quantum_graph``'s ``bond_order_calibration.py`` from
``molecular_quantum_graph/data/bond_order_model.json``.  Five functional forms
are supported:

.. list-table::
   :header-rows: 1
   :widths: 12 35 18

   * - Model
     - Formula :math:`f(b)`
     - Default params
   * - ``linear``  *(fallback)*
     - :math:`1 + \alpha\,(b-1)`
     - :math:`\alpha = 0.3`
   * - ``power``
     - :math:`b^{\,\beta}`
     - :math:`\beta = 0.6`
   * - ``log``
     - :math:`1 + a\,\ln b`
     - :math:`a = 1.0`
   * - ``poly2``
     - :math:`1 + a\,(b-1) + c\,(b-1)^2`
     - —
   * - ``poly3``
     - :math:`1 + a\,(b-1) + b_2\,(b-1)^2 + c\,(b-1)^3`
     - —

When ``bond_order_model.json`` is absent (e.g. ``molecular_quantum_graph`` is
not installed), the **linear model with** :math:`\alpha = 0.3` is used
automatically.

Resistance distance and Kirchhoff index
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The **weighted Laplacian** :math:`L` is assembled from the per-bond conductances.
The exact resistance distance between atoms :math:`u` and :math:`v` is then:

.. math::

   R(u,v) = L^+_{uu} + L^+_{vv} - 2\,L^+_{uv}

where :math:`L^+` is the **Moore–Penrose pseudoinverse** of :math:`L` (computed
via ``numpy.linalg.pinv``).

The **Kirchhoff index** (reported as ``effective_graph_resistance``) satisfies

.. math::

   K_f = \sum_{u < v} R(u, v) = n \sum_{i=2}^{n} \frac{1}{\lambda_i(L)}

Physical properties of :math:`R(u,v)`:

* **Symmetry**: :math:`R(u,v) = R(v,u)`
* **Positivity**: :math:`R(u,v) > 0` for :math:`u \neq v` in a connected graph
* **Triangle inequality**: :math:`R(i,k) \leq R(i,j) + R(j,k)`
* **Monotone with chain length**: for a linear homologous PFCA series,
  :math:`K_f` increases strictly with :math:`n` (fluorinated carbons)
* **Branching reduces** :math:`K_f`: branched isomers have lower :math:`K_f`
  than the linear chain of the same carbon count, because branching shortens
  the maximum pairwise path

Accessing resistance metrics
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

All resistance values are available per-component in the parse results:

.. code-block:: python

   from PFASGroups import parse_smiles

   results = parse_smiles('OC(=O)C(F)(F)C(F)(F)C(F)(F)F', halogens='F')
   comp = results[0]['matches'][0]['components'][0]

   # BDE-weighted Kirchhoff index
   print(comp['effective_graph_resistance'])

   # BDE-weighted resistance from functional group to structural landmarks
   print(comp['min_resistance_dist_to_barycentre'])
   print(comp['min_resistance_dist_to_centre'])
   print(comp['max_resistance_dist_to_periphery'])

To limit computation time on very large components:

.. code-block:: python

   # Only compute resistance for components with < 200 atoms
   results = parse_smiles(smiles, limit_effective_graph_resistance=200)

   # Skip resistance entirely
   results = parse_smiles(smiles, limit_effective_graph_resistance=0)

**Atom count** — number of heavy atoms in the component.

**Halogen count** — number of halogen atoms in the component.

**Halogen fraction** — ratio of halogen atoms to total heavy atoms.

Fingerprint generation
-----------------------

:func:`~PFASGroups.generate_fingerprint` converts a list of SMILES to a
2-D numpy array.  The default mode is fluorine-only (``halogens='F'``),
producing **116 columns** — one per group.  The column layout for multi-halogen
mode is:

``[group_0_F, group_0_Cl, group_0_Br, group_0_I,``
``  group_1_F, …,``
``  group_115_F, group_115_Cl, group_115_Br, group_115_I]``

Default F-only column count: 116 × 1 = **116**.
All-halogen column count: 116 × 4 = **464** (see :ref:`multi_halogen_fingerprint`).

Four count encoding values are available as items in ``component_metrics``:

.. list-table::
   :header-rows: 1

   * - component_metrics value
     - Cell value
   * - ``'binary'`` (default)
     - 1 if any match exists, 0 otherwise
   * - ``'count'``
     - Number of matching components
   * - ``'max_component'``
     - Size of the largest matching component (atom count)
   * - ``'total_component'``
     - Sum of all matching component sizes (atom count)

PFAS definition classification
--------------------------------

When ``include_PFAS_definitions=True``, each molecule is additionally
evaluated against five regulatory PFAS definitions.  Each definition is
encoded as a set of logical rules operating on the group matches already
computed.  No additional SMARTS matching is performed.

See :doc:`pfas_definitions` for the rule logic of each definition.