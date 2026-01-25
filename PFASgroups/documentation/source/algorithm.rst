Algorithm Overview
==================

This section provides a detailed explanation of the algorithms used in PFASgroups
for PFAS detection and classification.

Core Algorithm
--------------

The PFASgroups algorithm consists of four main phases:

1. **SMARTS Matching**: Find functional groups using SMARTS patterns
2. **Path Finding**: Identify fluorinated chains between functional groups
3. **Constraint Validation**: Apply molecular formula constraints
4. **Classification**: Assign PFAS group labels

.. figure:: _static/algorithm_flowchart.png
   :alt: Algorithm Flowchart
   :align: center
   :width: 80%

   *High-level flowchart of the PFASgroups algorithm*

Phase 1: SMARTS Matching
------------------------

The algorithm first identifies functional groups using SMARTS (SMiles ARbitrary 
Target Specification) patterns.

**Input Processing:**

.. code-block:: python

   from rdkit import Chem

   # Convert SMILES to RDKit molecule
   mol = Chem.MolFromSmiles("FC(F)(F)C(F)(F)C(=O)O")
   
   # Add hydrogens for accurate matching
   mol = Chem.AddHs(mol)

**SMARTS Pattern Matching:**

Each PFAS group defines one or two SMARTS patterns:

- ``smarts1``: Primary functional group (e.g., carboxylic acid)
- ``smarts2``: Secondary group (e.g., trifluoromethyl terminal)

.. code-block:: python

   # Example: PFCA (Perfluoroalkyl carboxylic acid)
   smarts1 = "C(=O)O"   # Carboxylic acid
   smarts2 = None       # Uses default from pathway
   
   # Find matches
   pattern = Chem.MolFromSmarts(smarts1)
   matches = mol.GetSubstructMatches(pattern)

Phase 2: Path Finding
---------------------

The algorithm identifies fluorinated chains connecting functional groups using
graph-based pathfinding.

**Molecular Graph Construction:**

.. code-block:: python

   import networkx as nx
   from rdkit import Chem

   def mol_to_nx(mol):
       """Convert RDKit molecule to NetworkX graph."""
       G = nx.Graph()
       
       for atom in mol.GetAtoms():
           G.add_node(atom.GetIdx(), 
                     element=atom.GetAtomicNum(),
                     symbol=atom.GetSymbol())
       
       for bond in mol.GetBonds():
           G.add_edge(bond.GetBeginAtomIdx(),
                     bond.GetEndAtomIdx(),
                     order=bond.GetBondTypeAsDouble())
       
       return G

**Shortest Path Algorithm:**

The algorithm uses NetworkX's shortest path algorithm to find the minimum
path between functional groups:

.. code-block:: python

   # Find shortest path between two atom indices
   path = nx.shortest_path(G, source=smarts1_match, target=smarts2_match)
   
   # Filter to include only fluorinated carbons
   chain_length = sum(1 for idx in path 
                      if G.nodes[idx]['symbol'] == 'C' 
                      and is_fluorinated(mol, idx))

**Chain SMARTS Validation:**

Paths are validated against pathway SMARTS patterns:

.. code-block:: text

   Perfluoroalkyl chain: [C;X4](F)(F)!@!=!#[C;X4](F)(F)
   
   This matches: -CF₂-CF₂- bonds
   But not: -CF₂-CH₂- bonds (wrong fluorination)
   And not: aromatic C-C bonds (ring bonds)

Phase 3: Constraint Validation
------------------------------

Molecular formula constraints filter out false positives:

**Constraint Types:**

.. list-table::
   :header-rows: 1
   :widths: 15 35 50

   * - Type
     - Syntax
     - Description
   * - eq
     - ``{"O": 2}``
     - Exactly 2 oxygen atoms
   * - gte
     - ``{"F": 3}``
     - At least 3 fluorine atoms
   * - lte
     - ``{"Cl": 1}``
     - At most 1 chlorine atom
   * - only
     - ``["C","F","H","O"]``
     - Only these elements allowed
   * - rel
     - see below
     - Relative element ratios

**Relative Constraints:**

.. code-block:: python

   # Carbon count must equal (F + H) / 2 - 1
   constraints = {
       "rel": {
           "C": {
               "atoms": ["F", "H"],
               "div": 2,
               "add": -1
           }
       }
   }

**Constraint Validation Code:**

.. code-block:: python

   def formula_dict_satisfies_constraints(formula_dict, constraints):
       """Check if formula satisfies all constraints."""
       
       # Check exact constraints
       if 'eq' in constraints:
           for elem, count in constraints['eq'].items():
               if formula_dict.get(elem, 0) != count:
                   return False
       
       # Check minimum constraints
       if 'gte' in constraints:
           for elem, min_count in constraints['gte'].items():
               if formula_dict.get(elem, 0) < min_count:
                   return False
       
       # Check maximum constraints
       if 'lte' in constraints:
           for elem, max_count in constraints['lte'].items():
               if formula_dict.get(elem, 0) > max_count:
                   return False
       
       # Check element restrictions
       if 'only' in constraints:
           allowed = set(constraints['only'])
           for elem in formula_dict:
               if elem not in allowed and formula_dict[elem] > 0:
                   return False
       
       return True

Phase 4: Classification
-----------------------

The final phase assigns PFAS group labels and computes comprehensive metrics.

**Group Priority:**

Groups are processed in order of specificity:

1. OECD-defined groups (IDs 1-28) - most specific
2. Generic functional groups (IDs 29-57) - broader categories

**Comprehensive Metrics:**

For each matched PFAS group, the algorithm computes:

- **Structural metrics**: branching, diameter, radius, center, periphery
- **Coverage metrics**: component_fraction (per component), total_components_fraction (union)
- **Graph-theoretic metrics**: mean/median eccentricity, effective resistance
- **Positional metrics**: distances from functional groups to structural features

.. code-block:: python

   # Component fraction calculation
   full_component = set(component)  # Carbon backbone
   
   # Add all attached atoms (H, F, Cl, Br, I)
   for atom_idx in component:
       for neighbor in mol.GetAtomWithIdx(atom_idx).GetNeighbors():
           full_component.add(neighbor.GetIdx())
   
   # Add SMARTS-matched atoms
   if smarts_matches:
       full_component = full_component.union(smarts_matches)
   
   component_fraction = len(full_component) / mol.GetNumAtoms()

**Output Format:**

.. code-block:: python

   results = [
       (PFASGroup, match_count, chain_lengths, matched_components),
       ...
   ]
   
   # Each matched_component includes:
   # - component: atom indices
   # - component_fraction: 0.0-1.0 (molecule coverage)
   # - branching: 1.0=linear, 0.0=branched
   # - mean_eccentricity, median_eccentricity
   # - diameter, radius, center, periphery
   # - and more...

bycomponent Mode
----------------

The ``bycomponent`` parameter enables alternative analysis for complex structures:

**Standard Mode (bycomponent=False):**

- Analyzes the entire molecule as one unit
- Finds paths between any pair of functional groups
- Best for linear chains

**Component Mode (bycomponent=True):**

- Fragments the molecule at functional groups
- Analyzes each fragment independently
- Better for cyclic structures and multi-chain molecules

.. code-block:: python

   from PFASgroups import parse_smiles

   # Standard mode
   results_std = parse_smiles(smiles, bycomponent=False)
   
   # Component mode
   results_comp = parse_smiles(smiles, bycomponent=True)

max_dist_from_CF Parameter
--------------------------

For groups without explicit formula constraints, the ``max_dist_from_CF``
parameter limits how far a functional group can be from a fluorinated carbon:

.. code-block:: json

   {
       "id": 50,
       "name": "Custom Group",
       "smarts1": "[N;!R]",
       "smartsPath": "Perfluoroalkyl",
       "constraints": {},
       "max_dist_from_CF": 2
   }

This ensures the functional group is within 2 bonds of a CF₂ unit.

Computational Complexity
------------------------

**Time Complexity:**

.. list-table::
   :header-rows: 1
   :widths: 30 30 40

   * - Phase
     - Complexity
     - Notes
   * - SMARTS Matching
     - O(N × M × S)
     - N=atoms, M=groups, S=SMARTS size
   * - Path Finding
     - O(V + E)
     - Dijkstra's algorithm
   * - Constraint Validation
     - O(E)
     - E=elements in formula
   * - Overall
     - O(N × M × S)
     - Dominated by SMARTS matching

**Memory Usage:**

- Molecular graph: O(V + E) per molecule
- SMARTS patterns: Pre-compiled, constant
- Results: O(G) per molecule, G=matched groups

Performance Optimization
------------------------

**Pre-compilation:**

SMARTS patterns are pre-compiled on module load:

.. code-block:: python

   # Patterns compiled once
   smarts_mol = Chem.MolFromSmarts(smarts_string)
   
   # Reused for all molecules
   matches = mol.GetSubstructMatches(smarts_mol)

**Parallel Processing:**

For batch processing, use multiprocessing:

.. code-block:: python

   from multiprocessing import Pool
   from PFASgroups import parse_smiles

   def process_batch(smiles_batch):
       return parse_smiles(smiles_batch)

   with Pool(processes=4) as pool:
       results = pool.map(process_batch, smiles_batches)

See Also
--------

- :doc:`api/core` - Core module API documentation
- :doc:`api/models` - Data model documentation
- :doc:`testing` - Testing and benchmarking
