Generic Classification Groups
==============================

This section documents the 27 generic PFAS classification groups (IDs 29-57)
that extend beyond the OECD-defined groups to provide comprehensive classification.

Overview
--------

The generic groups capture broader structural patterns that may not fit into
specific OECD categories but are still important for PFAS classification. These
groups are organized by:

- **Perfluoroalkyl structures** (complete fluorination)
- **Polyfluoroalkyl structures** (partial fluorination)
- **Branched and cyclic patterns**
- **Other halogenated compounds**

Perfluoroalkyl Generic Groups
-----------------------------

Group 29: Perfluoroalkyl
^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: text

   ID: 29
   Name: Perfluoroalkyl
   Description: Generic perfluoroalkyl chain
   Pathway: Perfluoroalkyl
   SMARTS: None (matches any perfluoroalkyl pathway)
   Constraints: None

Matches any compound with a perfluorinated alkyl chain (CₙF₂ₙ₊₁-).

**Example:**

.. code-block:: python

   from PFASgroups import parse_smiles
   
   # Generic perfluoroalkyl compound
   results = parse_smiles(["FC(F)(F)C(F)(F)C(F)(F)N"])
   # Matches: Perfluoroalkyl (Group 29) if no specific OECD group applies

Group 30: Perfluoroalkyl with Acid
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: text

   ID: 30
   Name: Perfluoroalkyl-acid
   Description: Perfluoroalkyl with generic acid group
   Pathway: Perfluoroalkyl
   SMARTS: [C,S,P](=O)[O;H1]
   Constraints: None

Matches perfluoroalkyl compounds with any acid functional group.

Group 31: Perfluoroalkyl with Ester
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: text

   ID: 31
   Name: Perfluoroalkyl-ester
   Description: Perfluoroalkyl with ester linkage
   Pathway: Perfluoroalkyl
   SMARTS: C(=O)O[C,c]
   Constraints: None

Group 32: Perfluoroalkyl with Alcohol
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: text

   ID: 32
   Name: Perfluoroalkyl-OH
   Description: Perfluoroalkyl with alcohol group
   Pathway: Perfluoroalkyl
   SMARTS: [O;H1][C;X4]
   Constraints: eq: {O: 1}

Group 33: Perfluoroalkyl with Ether
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: text

   ID: 33
   Name: Perfluoroalkyl-ether
   Description: Perfluoroalkyl with ether linkage
   Pathway: Perfluoroalkyl
   SMARTS: [C;X4](F)(F)[O;X2][C;X4]
   Constraints: None

Polyfluoroalkyl Generic Groups
------------------------------

Group 34: Polyfluoroalkyl
^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: text

   ID: 34
   Name: Polyfluoroalkyl
   Description: Generic polyfluoroalkyl chain
   Pathway: Polyfluoroalkyl
   SMARTS: None (matches any polyfluoroalkyl pathway)
   Constraints: None

Matches any compound with a partially fluorinated alkyl chain containing both
C-F and C-H bonds.

**Example:**

.. code-block:: python

   from PFASgroups import parse_smiles
   
   # Polyfluoroalkyl compound
   results = parse_smiles(["FC(F)(F)C(F)(F)C(F)HCC"])
   # Matches: Polyfluoroalkyl (Group 34)

Group 35: Polyfluoroalkyl with Acid
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: text

   ID: 35
   Name: Polyfluoroalkyl-acid
   Description: Polyfluoroalkyl with generic acid group
   Pathway: Polyfluoroalkyl
   SMARTS: [C,S,P](=O)[O;H1]
   Constraints: None

Group 36: Polyfluoroalkyl with Ester
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: text

   ID: 36
   Name: Polyfluoroalkyl-ester
   Description: Polyfluoroalkyl with ester linkage
   Pathway: Polyfluoroalkyl
   SMARTS: C(=O)O[C,c]
   Constraints: None

Group 37: Polyfluoroalkyl with Alcohol
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: text

   ID: 37
   Name: Polyfluoroalkyl-OH
   Description: Polyfluoroalkyl with alcohol group
   Pathway: Polyfluoroalkyl
   SMARTS: [O;H1][C;X4]
   Constraints: eq: {O: 1}

Branched and Polyfluoro Groups
------------------------------

Group 38: Polyfluoro (Branched)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: text

   ID: 38
   Name: Polyfluoro
   Description: Generic polyfluorinated compound (branched)
   Pathway: Polyfluoro
   SMARTS: None
   Constraints: None

Matches compounds with branched fluorinated structures.

Group 39-45: Polyfluoro Functional Groups
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

These groups combine the Polyfluoro pathway with various functional groups:

.. list-table::
   :header-rows: 1
   :widths: 10 25 65

   * - ID
     - Name
     - Description
   * - 39
     - Polyfluoro-acid
     - Polyfluoro with acid group
   * - 40
     - Polyfluoro-ester
     - Polyfluoro with ester linkage
   * - 41
     - Polyfluoro-OH
     - Polyfluoro with alcohol group
   * - 42
     - Polyfluoro-ether
     - Polyfluoro with ether linkage
   * - 43
     - Polyfluoro-amine
     - Polyfluoro with amine group
   * - 44
     - Polyfluoro-halide
     - Polyfluoro with halide substituent
   * - 45
     - Polyfluoro-silane
     - Polyfluoro with silane group

Brominated Groups
-----------------

Group 46: Polyfluorobr
^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: text

   ID: 46
   Name: Polyfluorobr
   Description: Brominated polyfluorinated compound
   Pathway: Polyfluorobr
   SMARTS: None
   Constraints: None

Matches brominated polyfluorinated compounds.

Group 47-52: Polyfluorobr Functional Groups
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. list-table::
   :header-rows: 1
   :widths: 10 25 65

   * - ID
     - Name
     - Description
   * - 47
     - Polyfluorobr-acid
     - Brominated polyfluoro with acid
   * - 48
     - Polyfluorobr-ester
     - Brominated polyfluoro with ester
   * - 49
     - Polyfluorobr-OH
     - Brominated polyfluoro with alcohol
   * - 50
     - Polyfluorobr-ether
     - Brominated polyfluoro with ether
   * - 51
     - Polyfluorobr-amine
     - Brominated polyfluoro with amine
   * - 52
     - Polyfluorobr-halide
     - Brominated polyfluoro with additional halide

Cyclic and Special Structures
-----------------------------

Group 53: Perfluorocyclic
^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: text

   ID: 53
   Name: Perfluorocyclic
   Description: Perfluorinated cyclic structure
   Pathway: cyclic
   SMARTS: [C;R](F)(F)[C;R](F)(F)
   Constraints: None

Matches perfluorinated cyclic structures.

Group 54: Fluoroaromatic
^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: text

   ID: 54
   Name: Fluoroaromatic
   Description: Fluorinated aromatic ring
   Pathway: None
   SMARTS: [c]F
   Constraints: None

Matches compounds with fluorinated aromatic rings.

Group 55-57: Additional Special Groups
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. list-table::
   :header-rows: 1
   :widths: 10 25 65

   * - ID
     - Name
     - Description
   * - 55
     - Semifluorinated
     - Semifluorinated alkyl chains
   * - 56
     - Side-chain fluorinated polymer
     - Fluorinated polymer side chains
   * - 57
     - Other fluorinated
     - Other fluorinated compounds not matching above categories

Usage in Code
-------------

.. code-block:: python

   from PFASgroups import parse_smiles, get_PFASGroups, generate_fingerprint

   # Get all generic groups (IDs 29-57)
   generic_groups = [g for g in get_PFASGroups() if g.id >= 29]
   
   # Parse with all groups (OECD + generic)
   results = parse_smiles(smiles_list)
   
   # Separate OECD and generic matches
   for mol_results in results:
       oecd = [(g.name, count) for g, count, _, _ in mol_results if g.id <= 28]
       generic = [(g.name, count) for g, count, _, _ in mol_results if g.id > 28]
       print(f"OECD: {oecd}, Generic: {generic}")

   # Generate fingerprint with only generic groups
   fps, info = generate_fingerprint(smiles_list, selected_groups=range(29, 58))

Hierarchical Classification
---------------------------

Generic groups often provide broader classification when specific OECD groups
don't apply:

.. code-block:: text

   Most specific → Least specific
   
   PFOA → PFCA (OECD) → Perfluoroalkyl-acid → Perfluoroalkyl
   
   6:2 FTOH → FTOH (OECD) → Polyfluoroalkyl-OH → Polyfluoroalkyl

**Example:**

.. code-block:: python

   from PFASgroups import parse_smiles

   # This matches PFCA (specific) and Perfluoroalkyl-acid (generic)
   results = parse_smiles(["FC(F)(F)C(F)(F)C(=O)O"])
   groups = [g.name for g, _, _, _ in results[0]]
   # groups = ['PFCA', 'Perfluoroalkyl-acid', 'Perfluoroalkyl']

See Also
--------

- :doc:`oecd_groups` - OECD-defined groups
- :doc:`pathway_types` - Pathway type definitions
