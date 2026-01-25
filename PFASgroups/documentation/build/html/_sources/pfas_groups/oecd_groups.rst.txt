OECD-Defined Groups
===================

This section documents the 28 OECD-defined PFAS groups (IDs 1-28) that form the
core classification system in PFASgroups.

Overview
--------

The OECD (Organisation for Economic Co-operation and Development) has established
standardized definitions for PFAS groups based on functional groups and structural
features. These definitions are used in regulatory contexts worldwide.

The groups are organized by functional group type:

- **Acids and derivatives** (Groups 1-14)
- **Sulfonates and sulfonamides** (Groups 15-20)
- **Alcohols and ethers** (Groups 21-24)
- **Other functional groups** (Groups 25-28)

Acids and Derivatives
---------------------

Group 1: PFCA - Perfluoroalkyl Carboxylic Acids
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: text

   ID: 1
   Name: PFCA
   Full Name: Perfluoroalkyl carboxylic acids
   Pathway: Perfluoroalkyl
   SMARTS: C(=O)O
   Constraints: eq: {O: 2}, only: [C, F, H, O]

**Examples:**

- PFOA (C8): ``FC(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(=O)O``
- PFBA (C4): ``FC(F)(F)C(F)(F)C(F)(F)C(=O)O``
- PFHxA (C6): ``FC(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(=O)O``

**Code example:**

.. code-block:: python

   from PFASgroups import parse_smiles

   pfoa = "FC(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(=O)O"
   results = parse_smiles([pfoa])
   # Matches: PFCA (Group 1)

Group 2: n:2 FTCA - n:2 Fluorotelomer Carboxylic Acids
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: text

   ID: 2
   Name: n:2 FTCA
   Full Name: n:2 fluorotelomer carboxylic acids
   Pathway: Perfluoroalkyl
   SMARTS: [C;X4;!$(C(F)(F)F)](F)(F)[C;X4]([!F;!Cl;!Br;!I])([!F;!Cl;!Br;!I])C(=O)O
   Constraints: eq: {O: 2}, only: [C, F, H, O]

**Characteristic:** Contains a -CF₂-CH₂- linkage followed by carboxylic acid.

Group 3: n:3 FTCA - n:3 Fluorotelomer Carboxylic Acids
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: text

   ID: 3
   Name: n:3 FTCA
   Full Name: n:3 fluorotelomer carboxylic acids
   Pathway: Perfluoroalkyl
   SMARTS: [C;X4;!$(C(F)(F)F)](F)(F)[C;X4]([!F;!Cl;!Br;!I])([!F;!Cl;!Br;!I])[C;X4]([!F;!Cl;!Br;!I])([!F;!Cl;!Br;!I])C(=O)O
   Constraints: eq: {O: 2}, only: [C, F, H, O]

**Characteristic:** Contains a -CF₂-CH₂-CH₂- linkage followed by carboxylic acid.

Group 4: PFECA - Perfluoroalkyl Ether Carboxylic Acids
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: text

   ID: 4
   Name: PFECA
   Full Name: Perfluoroalkyl ether carboxylic acids
   Pathway: Perfluoroalkyl
   SMARTS: [O;X2]([C;X4](F)(F))C(=O)O, FC(F)(F)[O;X2]
   Constraints: gte: {O: 3}, only: [C, F, H, O]

**Examples:**

- GenX (HFPO-DA): ``OC(=O)C(F)(OC(F)(F)C(F)(F)F)C(F)(F)F``

Group 5: PFSA - Perfluoroalkyl Sulfonic Acids
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: text

   ID: 5
   Name: PFSA
   Full Name: Perfluoroalkyl sulfonic acids
   Pathway: Perfluoroalkyl
   SMARTS: S(=O)(=O)O
   Constraints: eq: {S: 1, O: 3}, only: [C, F, H, O, S]

**Examples:**

- PFOS (C8): ``FC(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)S(=O)(=O)O``
- PFBS (C4): ``FC(F)(F)C(F)(F)C(F)(F)C(F)(F)S(=O)(=O)O``
- PFHxS (C6): ``FC(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)S(=O)(=O)O``

Group 6: n:2 FTSA - n:2 Fluorotelomer Sulfonic Acids
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: text

   ID: 6
   Name: n:2 FTSA
   Full Name: n:2 fluorotelomer sulfonic acids
   Pathway: Perfluoroalkyl
   SMARTS: [C;X4;!$(C(F)(F)F)](F)(F)[C;X4]([!F;!Cl;!Br;!I])([!F;!Cl;!Br;!I])S(=O)(=O)O
   Constraints: eq: {S: 1, O: 3}, only: [C, F, H, O, S]

Group 7: PFESA - Perfluoroalkyl Ether Sulfonic Acids
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: text

   ID: 7
   Name: PFESA
   Full Name: Perfluoroalkyl ether sulfonic acids
   Pathway: Perfluoroalkyl
   SMARTS: [O;X2]([C;X4](F)(F))S(=O)(=O)[O;H1]
   Constraints: gte: {O: 4}, eq: {S: 1}, only: [C, F, H, O, S]

Sulfonamides
------------

Group 8: FASA - Perfluoroalkane Sulfonamides
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: text

   ID: 8
   Name: FASA
   Full Name: Perfluoroalkane sulfonamides
   Pathway: Perfluoroalkyl
   SMARTS: [N;!R;X3;!$(NC=O)]S(=O)(=O)
   Constraints: eq: {S: 1, O: 2, N: 1}, only: [C, F, H, O, N, S]

Group 9: FASAA - Perfluoroalkane Sulfonamido Acetic Acids
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: text

   ID: 9
   Name: FASAA
   Full Name: Perfluoroalkane sulfonamido acetic acids
   Pathway: Perfluoroalkyl
   SMARTS: [N;!R;X3]([C;!R](=O)[O;H1])S(=O)(=O)
   Constraints: eq: {S: 1, O: 4, N: 1}, only: [C, F, H, O, N, S]

Alcohols and Ethers
-------------------

Group 10: FTOH - Fluorotelomer Alcohols
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: text

   ID: 10
   Name: FTOH
   Full Name: Fluorotelomer alcohols
   Pathway: Perfluoroalkyl
   SMARTS: [C;X4;!$(C(F)(F)F)](F)(F)[C;X4]([!F;!Cl;!Br;!I])([!F;!Cl;!Br;!I])[O;H1]
   Constraints: eq: {O: 1}, only: [C, F, H, O]

**Examples:**

- 6:2 FTOH: ``FC(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)CCO``

Group 11: PFAL - Perfluoroalkyl Aldehydes
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: text

   ID: 11
   Name: PFAL
   Full Name: Perfluoroalkyl aldehydes
   Pathway: Perfluoroalkyl
   SMARTS: [C;X3;H1](=O)
   Constraints: eq: {O: 1}, only: [C, F, O]

Phosphorus-Containing Groups
----------------------------

Group 12: PAP - Polyfluoroalkyl Phosphate Esters
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: text

   ID: 12
   Name: PAP
   Full Name: Polyfluoroalkyl phosphate esters
   Pathway: Perfluoroalkyl
   SMARTS: [O;X2]([P;X4](=O)([O,OH])([O,OH]))[C;X4]
   Constraints: eq: {P: 1}, only: [C, F, H, O, P]

Group 13: diPAP - Polyfluoroalkyl Phosphate Diesters
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: text

   ID: 13
   Name: diPAP
   Full Name: Polyfluoroalkyl phosphate diesters
   Pathway: Perfluoroalkyl
   SMARTS: [O;X2]([P;X4](=O)([O;X2][C;X4])([O,OH]))[C;X4]
   Constraints: eq: {P: 1}, only: [C, F, H, O, P]

Cyclic Structures
-----------------

Group 14: PFECHS - Perfluoro-4-ethylcyclohexane Sulfonate
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: text

   ID: 14
   Name: PFECHS
   Full Name: Perfluoro-4-ethylcyclohexane sulfonate
   Pathway: cyclic
   SMARTS: [C;R;X4](F)([C;R;X4](F)[C;R;X4](F)[C;R;X4](F)[C;R;X4](F)[C;R;X4](F)S(=O)(=O)[O])[C;X4](F)(F)
   Constraints: eq: {S: 1, O: 3}

Other OECD Groups
-----------------

The remaining OECD groups (15-28) include:

.. list-table::
   :header-rows: 1
   :widths: 10 30 60

   * - ID
     - Name
     - Description
   * - 15
     - PFPIA
     - Perfluoroalkyl phosphinic acids
   * - 16
     - PFPA
     - Perfluoroalkyl phosphonic acids
   * - 17
     - PFSF
     - Perfluoroalkane sulfonyl fluorides
   * - 18
     - PASF
     - Perfluoroalkane sulfonamidoethyl acrylates
   * - 19
     - MeFASE
     - N-Methyl perfluoroalkane sulfonamidoethanols
   * - 20
     - EtFASE
     - N-Ethyl perfluoroalkane sulfonamidoethanols
   * - 21
     - MeFASA
     - N-Methyl perfluoroalkane sulfonamides
   * - 22
     - EtFASA
     - N-Ethyl perfluoroalkane sulfonamides
   * - 23
     - MeFASAA
     - N-Methyl perfluoroalkane sulfonamido acetic acids
   * - 24
     - EtFASAA
     - N-Ethyl perfluoroalkane sulfonamido acetic acids
   * - 25
     - FTI
     - Fluorotelomer iodides
   * - 26
     - FTA
     - Fluorotelomer acrylates
   * - 27
     - FTMA
     - Fluorotelomer methacrylates
   * - 28
     - FTAC
     - Fluorotelomer acrylate copolymers

Usage in Code
-------------

.. code-block:: python

   from PFASgroups import parse_smiles, get_PFASGroups

   # Get all OECD groups (IDs 1-28)
   oecd_groups = [g for g in get_PFASGroups() if g.id <= 28]
   
   # Filter results to only OECD groups
   results = parse_smiles(smiles_list)
   for mol_results in results:
       oecd_matches = [(g, count, chains, matched) 
                       for g, count, chains, matched in mol_results 
                       if g.id <= 28]
       print(oecd_matches)

   # Generate fingerprint with only OECD groups
   from PFASgroups import generate_fingerprint
   fps, info = generate_fingerprint(smiles_list, selected_groups=range(1, 29))

See Also
--------

- :doc:`generic_groups` - Generic classification groups
- :doc:`pathway_types` - Pathway type definitions
