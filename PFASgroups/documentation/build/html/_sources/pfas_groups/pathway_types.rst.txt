Pathway Types
=============

This section documents the pathway types used to classify fluorinated chains
in PFASgroups.

Overview
--------

Pathway types define the repeating structural units that form fluorinated chains.
Each pathway has:

- **Chain SMARTS**: Pattern for the repeating unit
- **End SMARTS**: Pattern for the terminal group

The module supports four main pathway types:

.. list-table::
   :header-rows: 1
   :widths: 20 40 40

   * - Pathway
     - Chain Unit
     - Description
   * - Perfluoroalkyl
     - -CF₂-CF₂-
     - Fully fluorinated chain
   * - Polyfluoroalkyl
     - -CF₂-CHF- or -CF₂-CH₂-
     - Partially fluorinated chain
   * - Polyfluoro
     - Various
     - Branched/irregular fluorination
   * - Polyfluorobr
     - Contains Br
     - Brominated polyfluorinated

Perfluoroalkyl
--------------

.. code-block:: text

   Name: Perfluoroalkyl
   Chain SMARTS: [C;X4](F)(F)!@!=!#[C;X4](F)(F)
   End SMARTS: [C;X4](F)(F)F
   Description: Fully fluorinated alkyl chain

**Structural Formula:**

.. math::

   C_nF_{2n+1}-

**Characteristics:**

- All carbon atoms are bonded to fluorine
- Linear or branched chain
- Highly stable and persistent
- Most common in legacy PFAS compounds

**Examples:**

.. code-block:: python

   # PFOA - C8 perfluoroalkyl carboxylic acid
   "FC(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(=O)O"
   
   # PFOS - C8 perfluoroalkyl sulfonic acid
   "FC(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)S(=O)(=O)O"

**Chain Length Calculation:**

The chain length counts the number of fully fluorinated carbon atoms:

.. code-block:: python

   from PFASgroups import parse_smiles

   results = parse_smiles(["FC(F)(F)C(F)(F)C(F)(F)C(=O)O"])
   for group, count, chain_lengths, matched in results[0]:
       print(f"{group.name}: chain length = {chain_lengths}")
   # PFCA: chain length = [3]

Polyfluoroalkyl
---------------

.. code-block:: text

   Name: Polyfluoroalkyl
   Chain SMARTS: [C;X4]([F,H,I,Br,Cl])([F,H,I,Br,Cl])!@!=!#[C;X4]([F,H,I,Br,Cl])([F,H,I,Br,Cl])
   End SMARTS: [C;X4]([F,H,I,Br,Cl])([F,H,I,Br,Cl])[F,H,I,Br,Cl]
   Description: Partially fluorinated alkyl chain

**Structural Formula:**

.. math::

   C_nF_xH_y-  \text{ where } x + y = 2n + 1

**Characteristics:**

- Mix of C-F and C-H bonds
- Includes fluorotelomer structures
- May contain halogens (Cl, Br, I)
- Common in newer PFAS alternatives

**Examples:**

.. code-block:: python

   # 6:2 FTOH - Fluorotelomer alcohol
   "FC(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)CCO"
   
   # 4:2 FTCA - Fluorotelomer carboxylic acid
   "FC(F)(F)C(F)(F)C(F)(F)C(F)(F)CC(=O)O"

**Difference from Perfluoroalkyl:**

.. code-block:: text

   Perfluoroalkyl: F-C(F)(F)-C(F)(F)-C(F)(F)-COOH
                   ↑ all C-F bonds
   
   Polyfluoroalkyl: F-C(F)(F)-C(F)(F)-C(F)(F)-C(H)(H)-COOH
                                              ↑ C-H bonds

Polyfluoro (Branched)
---------------------

.. code-block:: text

   Name: Polyfluoro
   Chain SMARTS: [C;X4]([F,H])([F,H])!@!=!#[C;X4;H0,H1]
   End SMARTS: [C;X4;H0](F)(F)F
   Description: Branched polyfluorinated structures

**Characteristics:**

- Branched carbon backbone
- Variable fluorination pattern
- May include tertiary carbons
- Common in complex industrial PFAS

**Example:**

.. code-block:: python

   # Branched perfluorinated compound
   "FC(F)(F)C(F)(C(F)(F)F)C(F)(F)C(=O)O"

Polyfluorobr (Brominated)
-------------------------

.. code-block:: text

   Name: Polyfluorobr
   Chain SMARTS: [C;X4]([F,H,Br])([F,H,Br])!@!=!#[C;X4]([F,H,Br])([F,H,Br])
   End SMARTS: [C;X4]([F,H,Br])([F,H,Br])[F,H,Br]
   Description: Brominated polyfluorinated structures

**Characteristics:**

- Contains bromine atoms
- Often used in flame retardants
- May have mixed Br/F substitution

Custom Pathway Definition
-------------------------

You can define custom pathway types:

.. code-block:: python

   from PFASgroups import get_smartsPaths, compile_smartsPath, parse_smiles

   # Get default paths
   paths = get_smartsPaths()

   # Add a custom path (e.g., perchlorinated)
   paths['Perchlorinated'] = compile_smartsPath(
       chain_smarts="[C;X4](Cl)(Cl)!@!=!#[C;X4](Cl)(Cl)",
       end_smarts="[C;X4](Cl)(Cl)Cl"
   )

   # Use custom paths
   results = parse_smiles(smiles_list, smartsPaths=paths)

Pathway Selection in Groups
---------------------------

Each PFAS group is associated with a specific pathway type:

.. code-block:: python

   from PFASgroups import get_PFASGroups

   groups = get_PFASGroups()
   
   # Show pathway for each group
   for group in groups[:10]:
       print(f"{group.name}: {group.smartsPath}")

   # Output:
   # PFCA: Perfluoroalkyl
   # n:2 FTCA: Perfluoroalkyl
   # PFSA: Perfluoroalkyl
   # FTOH: Perfluoroalkyl
   # ...

Algorithm Details
-----------------

**Chain Detection Algorithm:**

1. Find all matches of the chain SMARTS pattern
2. Find all matches of the end SMARTS pattern
3. Build a molecular graph
4. Find shortest path between end groups through chain atoms
5. Count chain length as number of carbon atoms in path

**Handling Cyclic Structures:**

For cyclic structures, the algorithm uses a different approach:

.. code-block:: python

   # Cyclic pathway
   smartsPath = "cyclic"
   
   # Uses connected component analysis instead of path finding
   # Matches cyclic fluorinated structures as a whole

**Multiple Chains:**

If a molecule has multiple fluorinated chains:

.. code-block:: python

   from PFASgroups import parse_smiles

   # Molecule with two fluorinated chains
   smiles = "FC(F)(F)C(F)(F)CC(F)(F)C(F)(F)F"
   results = parse_smiles([smiles])
   
   for group, count, chain_lengths, matched_chains in results[0]:
       print(f"{group.name}: count={count}, chains={chain_lengths}")
   # May match multiple times with different chain lengths

See Also
--------

- :doc:`oecd_groups` - OECD-defined groups
- :doc:`generic_groups` - Generic classification groups
- :doc:`../api/core` - Core module API
