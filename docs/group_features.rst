Group Feature Extraction
========================

:func:`~PFASGroups.extract_group_features` provides a single entry point for
extracting structured numerical features from halogen-group detection.  It
is designed for machine-learning pipelines that require a **fixed-length,
named feature vector** from the six component groups and the generic
functional groups of PFASGroups.

.. contents:: Contents
   :local:
   :depth: 2

Overview
--------

Two :func:`~PFASGroups.parser.parse_mol` calls are made per molecule:

1. ``halogens=['H', 'F', 'Cl', 'Br', 'I']`` restricted to the **six component
   groups** (ids 34, 35, 37, 38, 44, 45) — captures structured halogen-chain
   information including alkyl chains detected via the ``'H'`` pseudo-halogen.
2. ``halogens=['*']`` (wildcard) — captures **generic functional group** matches
   (ids 29-76).  The six component group ids are automatically excluded from
   wildcard matching.

The result is a :class:`~PFASGroups.GroupFeatureResult` dataclass with four
named dictionaries and a :meth:`~PFASGroups.GroupFeatureResult.to_array`
method that returns a fixed-length ``float32`` array of shape ``(66,)``.

The six component groups
------------------------

These groups describe the **halogenation pattern** of carbon-chain or ring
components.  They are parameterised per halogen type and are the basis for
the :attr:`~PFASGroups.GroupFeatureResult.poly_counts` and
:attr:`~PFASGroups.GroupFeatureResult.per_halogen_sizes` feature sets.

.. list-table::
   :header-rows: 1
   :widths: 8 10 25 55

   * - ID
     - Category
     - Name
     - Description
   * - 34
     - Perhalogenated
     - perhalogenated alkyl
     - All halogen-bearing carbons in the alkyl chain carry *one single* halogen
       type (F, Cl, Br, I, or the H pseudo-halogen for un-substituted chains).
       Match is attributed to that halogen.
   * - 35
     - Polyhalogenated
     - polyhalogenated alkyl
     - The alkyl chain has halogen-bearing carbons but carries a *mix* of
       halogens or only partial halogenation.  Match is counted regardless of
       which halogen(s) are present.
   * - 37
     - Perhalogenated
     - perhalogenated aryl compounds
     - Like group 34, but the halogenated component is aromatic.
   * - 38
     - Polyhalogenated
     - polyhalogenated aryl compounds
     - Like group 35, but the halogenated component is aromatic.
   * - 44
     - Perhalogenated
     - perhalogenated cyclic compounds
     - Like group 34, but the halogenated component is cyclic (non-aromatic).
   * - 45
     - Polyhalogenated
     - polyhalogenated cyclic compounds
     - Like group 35, but the halogenated component is cyclic (non-aromatic).

.. note::

   The six component group ids (34, 35, 37, 38, 44, 45) are **excluded from
   wildcard matching** by PFASGroups design.  Their contribution to the
   generic-group feature vector is therefore always zero.

The H pseudo-halogen
--------------------

Group 34 (perhalogenated alkyl) is compiled with ``componentHalogens``
inferred from the available component SMARTS, which includes the ``'Alkyl'``
component type.  This allows PFASGroups to detect *un-substituted* (all-H)
alkyl chains and attribute them to the ``'H'`` pseudo-halogen.

The ``h_chain_sizes`` dictionary in :class:`~PFASGroups.GroupFeatureResult`
is a convenience view of the ``'H'`` column of ``per_halogen_sizes``:

.. code-block:: python

   r.h_chain_sizes['alkyl_H']   # == r.per_halogen_sizes['g34_H']
   r.h_chain_sizes['aryl_H']    # == r.per_halogen_sizes['g37_H']
   r.h_chain_sizes['cyclic_H']  # == r.per_halogen_sizes['g44_H']

Note that ``h_chain_sizes`` is **not** included in
:meth:`~PFASGroups.GroupFeatureResult.to_array`; it is purely for convenient
named access.

Generic functional groups (ids 29-76)
--------------------------------------

The wildcard call captures all 48 functional group ids between 29 and 76.
Group names and their IDs:

.. list-table::
   :header-rows: 1
   :widths: 8 30 8 30

   * - ID
     - Name
     - ID
     - Name
   * - 29
     - acrylate
     - 30
     - acyl halide
   * - 31
     - alcohol
     - 32
     - aldehyde
   * - 33
     - alkene
     - 34
     - perhalogenated alkyl *(always 0)*
   * - 35
     - polyhalogenated alkyl *(always 0)*
     - 36
     - alkyne
   * - 37
     - perhalogenated aryl compounds *(always 0)*
     - 38
     - polyhalogenated aryl compounds *(always 0)*
   * - 39
     - benzodioxole
     - 40
     - benzoyl peroxydes
   * - 41
     - bromide
     - 42
     - carboxylic acid
   * - 43
     - chloride
     - 44
     - perhalogenated cyclic compounds *(always 0)*
   * - 45
     - polyhalogenated cyclic compounds *(always 0)*
     - 46
     - ester
   * - 47
     - ether
     - 48
     - fluoride
   * - 49
     - glucuronate
     - 50
     - iodide
   * - 51
     - ketone
     - 52
     - methacrylate
   * - 53
     - peroxydes
     - 54
     - side-chain aromatics
   * - 55
     - sulfenic acid
     - 56
     - sulfenyl halide
   * - 57
     - sulfinic acid
     - 58
     - sulfinyl amido sulfonic acid
   * - 59
     - sulfonamide
     - 60
     - sulfonamidoethanol
   * - 61
     - sulfonic acid
     - 62
     - sulfonyl halide
   * - 63
     - sulfonyl propanoic acid
     - 64
     - sulfuric acid
   * - 65
     - thioester keto dicarboxylic acid
     - 66
     - thiocyanic acid
   * - 67
     - phosphinic acid
     - 68
     - phosphonic acid
   * - 69
     - amide
     - 70
     - amine
   * - 71
     - heterocyclic azine
     - 72
     - heterocyclic azole
   * - 73
     - betaine
     - 74
     - glycine
   * - 75
     - trichlorosilane
     - 76
     - silane

Array layout
------------

:meth:`~PFASGroups.GroupFeatureResult.to_array` returns 66 ``float32`` values:

.. code-block:: text

    [0]    poly_alkyl           group 35 match count
    [1]    poly_aryl            group 38 match count
    [2]    poly_cyclic          group 45 match count
    [3]    g34_F                max perfluoroalkyl component size
    [4]    g34_Cl               max perchloroalkyl component size
    [5]    g34_Br               max perbromoalkyl component size
    [6]    g34_I                max periodoalkyl component size
    [7]    g34_H                max alkyl (H chain) component size
    [8]    g37_F                max perfluorinated aryl component size
    [9]    g37_Cl               max perchlorinated aryl component size
    [10]   g37_Br               max perbrominated aryl component size
    [11]   g37_I                max periodinated aryl component size
    [12]   g37_H                max aryl (H) component size
    [13]   g44_F                max perfluorinated cyclic component size
    [14]   g44_Cl               max perchlorinated cyclic component size
    [15]   g44_Br               max perbrominated cyclic component size
    [16]   g44_I                max periodinated cyclic component size
    [17]   g44_H                max cyclic (H) component size
    [18]   g29                  acrylate count
    [19]   g30                  acyl halide count
    ...
    [65]   g76                  silane count

Feature name labels are available from
:meth:`~PFASGroups.GroupFeatureResult.feature_names`.

Quick start
-----------

Basic usage with an RDKit molecule or a SMILES string:

.. code-block:: python

   from rdkit import Chem
   from PFASGroups import extract_group_features

   # PFOA-like molecule
   mol = Chem.MolFromSmiles("FC(F)(F)C(F)(F)C(F)(F)C(=O)O")
   r = extract_group_features(mol)

   print(r)
   # GroupFeatureResult(poly_counts=1 nonzero, per_halogen_sizes=1 nonzero, generic_groups=2 nonzero)

   print(r.poly_counts)
   # {'poly_alkyl': 1.0, 'poly_aryl': 0.0, 'poly_cyclic': 0.0}

   print(r.per_halogen_sizes['g34_F'])  # largest perfluoroalkyl chain
   # 4.0

   arr = r.to_array()
   print(arr.shape, arr.dtype)
   # (66,) float32

SMILES strings are also accepted directly:

.. code-block:: python

   r = extract_group_features("CCCCCCCC")  # octane
   print(r.h_chain_sizes)
   # {'alkyl_H': 6.0, 'aryl_H': 0.0, 'cyclic_H': 0.0}

Worked examples
---------------

.. code-block:: python

   from rdkit import Chem
   from PFASGroups import extract_group_features, GENERIC_GROUP_NAMES

   molecules = {
       "PFOA":            "FC(F)(F)C(F)(F)C(F)(F)C(=O)O",
       "PFOS":            "FC(F)(F)C(F)(F)C(F)(F)S(=O)(=O)O",
       "mixed":           "CCCCC(F)(F)C(F)(F)C(=O)O",  # CF2 segment + alkyl tail
       "octane":          "CCCCCCCC",
       "hexafluorobenzene": "Fc1c(F)c(F)c(F)c(F)c1F",
       "perchloroalkyl":  "ClC(Cl)(Cl)C(Cl)(Cl)C(Cl)(Cl)Cl",
   }

   for name, smi in molecules.items():
       r = extract_group_features(smi)
       print(f"{name}")
       print(f"  poly:       {r.poly_counts}")
       print(f"  g34(alkyl): F={r.per_halogen_sizes['g34_F']}, H={r.per_halogen_sizes['g34_H']}")
       print(f"  g37(aryl):  F={r.per_halogen_sizes['g37_F']}")
       nz_gen = {GENERIC_GROUP_NAMES[int(k[1:])]: v for k, v in r.generic_groups.items() if v}
       if nz_gen:
           print(f"  generic:    {nz_gen}")
       print()

Expected output::

   PFOA
     poly:       {'poly_alkyl': 1.0, 'poly_aryl': 0.0, 'poly_cyclic': 0.0}
     g34(alkyl): F=4.0, H=0.0
     g37(aryl):  F=0.0
     generic:    {'carboxylic acid': 1.0, 'fluoride': 1.0}

   PFOS
     poly:       {'poly_alkyl': 1.0, 'poly_aryl': 0.0, 'poly_cyclic': 0.0}
     g34(alkyl): F=3.0, H=0.0
     g37(aryl):  F=0.0
     generic:    {'fluoride': 1.0, 'sulfonic acid': 1.0}

   mixed
     poly:       {'poly_alkyl': 1.0, 'poly_aryl': 0.0, 'poly_cyclic': 0.0}
     g34(alkyl): F=2.0, H=3.0
     g37(aryl):  F=0.0
     generic:    {'carboxylic acid': 1.0, 'fluoride': 1.0}

   octane
     poly:       {'poly_alkyl': 0.0, 'poly_aryl': 0.0, 'poly_cyclic': 0.0}
     g34(alkyl): F=0.0, H=6.0
     g37(aryl):  F=0.0

   hexafluorobenzene
     poly:       {'poly_alkyl': 0.0, 'poly_aryl': 1.0, 'poly_cyclic': 0.0}
     g34(alkyl): F=0.0, H=0.0
     g37(aryl):  F=6.0
     generic:    {'fluoride': 1.0}

   perchloroalkyl
     poly:       {'poly_alkyl': 1.0, 'poly_aryl': 0.0, 'poly_cyclic': 0.0}
     g34(alkyl): F=0.0, H=0.0
     g37(aryl):  F=0.0
     generic:    {'chloride': 1.0}

Using feature names
-------------------

.. code-block:: python

   from PFASGroups import extract_group_features
   from rdkit import Chem

   mol = Chem.MolFromSmiles("FC(F)(F)C(F)(F)C(=O)O")
   r = extract_group_features(mol)
   names = r.feature_names()   # list of 66 labels
   arr = r.to_array()

   # Non-zero features with their labels
   for name, val in zip(names, arr):
       if val:
           print(f"{name}: {val}")
   # poly_alkyl: 1.0
   # g34_F: 3.0
   # g42: 1.0   (carboxylic acid)
   # g48: 1.0   (fluoride)

API Reference
-------------

.. autofunction:: PFASGroups.extract_group_features

.. autoclass:: PFASGroups.GroupFeatureResult
   :members:
   :undoc-members:

.. autodata:: PFASGroups.PER_GROUP_IDS
.. autodata:: PFASGroups.POLY_GROUP_IDS
.. autodata:: PFASGroups.HALOGENS_ORDER
.. autodata:: PFASGroups.GENERIC_GROUP_VOCAB
.. autodata:: PFASGroups.GENERIC_GROUP_NAMES
