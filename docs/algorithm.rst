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

The library contains **116 halogen groups** organised into three categories:

.. list-table::
   :header-rows: 1
   :widths: 15 15 70

   * - Category
     - Count
     - Description
   * - OECD
     - 28
     - Groups adopted from the 2021 OECD PFAS definition framework
   * - Generic
     - 45
     - Broader halogenated structural motifs (alkyl, aryl, acyl, sulfonyl, …)
   * - Fluorotelomer
     - 43
     - Fluorotelomer groups including telomer alcohols, sulfonates, amides

In addition one *aggregate* group (``Telomers``) is defined with
``compute=False``; it is included in fingerprint headers but not matched
directly.

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

**Effective graph resistance** (Kirchhoff index)

.. math::

   R_e = n \sum_{i=2}^{n} \frac{1}{\lambda_i}

where :math:`\lambda_i` are the non-zero eigenvalues of the Laplacian and
:math:`n` is the number of atoms.

To limit computation time on very large components you can pass:

.. code-block:: python

   results = parse_smiles(smiles, limit_effective_graph_resistance=500)

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

Three ``count_mode`` options are available:

.. list-table::
   :header-rows: 1

   * - count_mode
     - Cell value
   * - ``'binary'`` (default)
     - 1 if any match exists, 0 otherwise
   * - ``'count'``
     - Number of matching components
   * - ``'max_component'``
     - Size of the largest matching component (atom count)

PFAS definition classification
--------------------------------

When ``include_PFAS_definitions=True``, each molecule is additionally
evaluated against five regulatory PFAS definitions.  Each definition is
encoded as a set of logical rules operating on the group matches already
computed.  No additional SMARTS matching is performed.

See :doc:`pfas_definitions` for the rule logic of each definition.