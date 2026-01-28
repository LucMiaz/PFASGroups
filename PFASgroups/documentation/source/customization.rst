Customization
=============

PFASgroups is designed to be easily extensible with custom PFAS groups and pathway types.

Creating Custom PFAS Groups
---------------------------

Python API
^^^^^^^^^^

Add custom groups programmatically:

.. code-block:: python

   from PFASgroups import get_PFASGroups, PFASGroup, parse_smiles

   # Get default groups
   groups = get_PFASGroups()

   # Create a custom group
   custom_group = PFASGroup(
       id=999,
       name="My Custom PFAS Group",
       smarts1="[C](F)(F)F",              # Primary SMARTS
       smarts2="[N+](=O)[O-]",             # Secondary SMARTS (optional)
       smartsPath="Perfluoroalkyl",        # Pathway type
       constraints={
           "gte": {"F": 3},                # At least 3 fluorines
           "eq": {"N": 1},                 # Exactly 1 nitrogen
       },
       max_dist_from_CF=2,                 # Max distance from fluorinated carbon
       linker_smarts="[CH2X4]"             # Only allow CH2 linker atoms (optional)
   )

   # Add to groups list
   groups.append(custom_group)

   # Use custom groups
   results = parse_smiles(["FC(F)(F)C(F)(F)[N+](=O)[O-]"], pfas_groups=groups)

PFASGroup Parameters
^^^^^^^^^^^^^^^^^^^^

.. list-table::
   :header-rows: 1
   :widths: 20 15 65

   * - Parameter
     - Type
     - Description
   * - ``id``
     - int
     - Unique identifier for the group
   * - ``name``
     - str
     - Descriptive name
   * - ``smarts1``
     - str
     - Primary SMARTS pattern for functional group detection
   * - ``smarts2``
     - str or None
     - Secondary SMARTS pattern (optional, for path finding)
   * - ``smartsPath``
     - str or None
     - Pathway type: 'Perfluoroalkyl', 'Polyfluoroalkyl', 'cyclic', or None
   * - ``constraints``
     - dict
     - Molecular formula constraints
   * - ``max_dist_from_CF``
     - int
     - Maximum bond distance from fluorinated carbon (default: 2)
   * - ``linker_smarts``
     - str or None
     - SMARTS pattern for validating linker atoms between fluorinated component and functional group (optional)

Constraint Types
^^^^^^^^^^^^^^^^

The ``constraints`` dictionary supports five types:

**Exact Count (eq)**

.. code-block:: python

   {"eq": {"O": 2, "S": 1}}  # Exactly 2 oxygens and 1 sulfur

**Greater Than or Equal (gte)**

.. code-block:: python

   {"gte": {"F": 3, "O": 1}}  # At least 3 fluorines and 1 oxygen

**Less Than or Equal (lte)**

.. code-block:: python

   {"lte": {"Cl": 2}}  # At most 2 chlorines

**Only These Elements (only)**

.. code-block:: python

   {"only": ["C", "F", "H", "O"]}  # Only these elements allowed

**Relative Ratios (rel)**

.. code-block:: python

   # C = (F + H) / 2 - 1
   {"rel": {"C": {"atoms": ["F", "H"], "add": -1, "div": 2}}}

   # With additional atoms to add
   {"rel": {"C": {"atoms": ["F"], "add": 0.5, "div": 2, "add_atoms": ["O"]}}}

Linker SMARTS Validation
^^^^^^^^^^^^^^^^^^^^^^^^

The ``linker_smarts`` parameter allows you to specify which atoms are permitted in the path between the fluorinated component and the functional group. This is particularly useful for fluorotelomer compounds where only specific linker atoms (e.g., CH2) should connect the perfluorinated chain to the functional group.

**Example: Fluorotelomer Alcohols**

Fluorotelomer alcohols have the structure CF3(CF2)n-(CH2)m-OH, where only CH2 units should link the perfluorinated chain to the alcohol group:

.. code-block:: python

   fluorotelomer_alcohol = PFASGroup(
       id=15,
       name="Fluorotelomer alcohols",
       smarts1="[#6X4$([C;!$([C]=O)][OH1])]",
       smartsPath="Perfluoroalkyl",
       linker_smarts="[#6H2X4]",  # Only CH2 atoms allowed in linker
       max_dist_from_CF=12,
       constraints={}
   )

**How It Works:**

1. The algorithm finds the shortest path between the fluorinated component and the functional group
2. If ``linker_smarts`` is specified and the path has intermediate atoms (length > 2)
3. Each intermediate atom (excluding endpoints) is validated against the linker pattern
4. If any atom doesn't match, the path is rejected

**Common Patterns:**

.. code-block:: python

   # Only CH2 units
   linker_smarts="[CH2X4]"
   
   # CH2 or oxygen (for ethoxylates)
   linker_smarts="[CH2X4,O$(O([CH2])[CH2])]"
   
   # Any sp3 carbon
   linker_smarts="[CX4]"
   
   # Sp3 carbon or nitrogen
   linker_smarts="[CX4,NX3]"

**Use Cases:**

- Fluorotelomer compounds with specific linker requirements
- Ensuring no heteroatoms in certain positions
- Distinguishing between direct attachment vs. linker-mediated attachment
- Validating chain composition in degradation products

JSON Configuration
^^^^^^^^^^^^^^^^^^

Create a JSON file for custom groups:

.. code-block:: json

   [
     {
       "id": 999,
       "name": "My Custom PFAS",
       "alias": "CustomPFAS",
       "base_functional_groups": ["custom"],
       "main_group": "Custom Groups",
       "smarts1": "[C](F)(F)F",
       "smarts2": null,
       "smartsPath": "Perfluoroalkyl",
       "linker_smarts": "[CH2X4]",
       "constraints": {
         "gte": {"F": 3},
         "only": ["C", "F", "H", "O"]
       },
       "max_dist_from_CF": 2
     }
   ]

Load and use:

.. code-block:: python

   from PFASgroups import get_PFASGroups, parse_smiles
   
   # Load from custom file
   groups = get_PFASGroups(filename='custom_groups.json')
   
   results = parse_smiles(smiles_list, pfas_groups=groups)

Creating Custom Pathway Types
-----------------------------

Python API
^^^^^^^^^^

.. code-block:: python

   from PFASgroups import compile_smartsPath, get_smartsPaths, parse_smiles

   # Get default pathways
   paths = get_smartsPaths()

   # Add a custom pathway (e.g., perchlorinated)
   paths['Perchlorinated'] = compile_smartsPath(
       chain_smarts="[C;X4](Cl)(Cl)!@!=!#[C;X4](Cl)(Cl)",  # Chain pattern
       end_smarts="[C;X4](Cl)(Cl)Cl"                        # End pattern
   )

   # Use with custom pathway
   results = parse_smiles(
       ["ClC(Cl)(Cl)C(Cl)(Cl)C(=O)O"],
       smartsPaths=paths
   )

Pathway SMARTS Patterns
^^^^^^^^^^^^^^^^^^^^^^^

A pathway is defined by two SMARTS patterns:

1. **Chain pattern**: Matches the repeating unit in the fluorinated chain
2. **End pattern**: Matches the terminal group of the chain

Default Perfluoroalkyl patterns:

.. code-block:: text

   Chain: [C;X4;H0](F)(F)!@!=!#[C;X4;H0](F)(F)
   End:   [C;X4;H0](F)(F)F

Default Polyfluoroalkyl patterns:

.. code-block:: text

   Chain: [C;X4;H1](F)!@!=!#[C;X4](F)
   End:   [C;X4;H1](F)F

Multiple Pathways
^^^^^^^^^^^^^^^^^

.. code-block:: python

   from PFASgroups import compile_smartsPaths

   custom_paths = {
       'Perchlorinated': {
           'chain': '[C;X4](Cl)(Cl)!@!=!#[C;X4](Cl)(Cl)',
           'end': '[C;X4](Cl)(Cl)Cl'
       },
       'MixedHalo': {
           'chain': '[C;X4]([F,Cl])([F,Cl])!@!=!#[C;X4]([F,Cl])',
           'end': '[C;X4]([F,Cl])([F,Cl])[F,Cl]'
       }
   }

   compiled = compile_smartsPaths(custom_paths)

JSON Configuration
^^^^^^^^^^^^^^^^^^

Create a custom ``fpaths.json``:

.. code-block:: json

   {
     "Perfluoroalkyl": {
       "chain": "[C;X4;H0](F)(F)!@!=!#[C;X4;H0](F)(F)",
       "end": "[C;X4;H0](F)(F)F"
     },
     "Polyfluoroalkyl": {
       "chain": "[C;X4;H1](F)!@!=!#[C;X4](F)",
       "end": "[C;X4;H1](F)F"
     },
     "Perchlorinated": {
       "chain": "[C;X4](Cl)(Cl)!@!=!#[C;X4](Cl)(Cl)",
       "end": "[C;X4](Cl)(Cl)Cl"
     }
   }

Load from file:

.. code-block:: python

   from PFASgroups import get_smartsPaths, parse_smiles

   paths = get_smartsPaths(filename='custom_fpaths.json')
   results = parse_smiles(smiles_list, smartsPaths=paths)

Extending vs Replacing
----------------------

Extending Defaults
^^^^^^^^^^^^^^^^^^

Add your custom definitions to the defaults:

.. code-block:: python

   from PFASgroups import get_PFASGroups, get_smartsPaths, PFASGroup, compile_smartsPath

   # Extend groups
   groups = get_PFASGroups()
   groups.append(PFASGroup(...))

   # Extend paths
   paths = get_smartsPaths()
   paths['Custom'] = compile_smartsPath(...)

Replacing Defaults
^^^^^^^^^^^^^^^^^^

Use only your custom definitions:

.. code-block:: python

   from PFASgroups import PFASGroup, compile_smartsPaths, parse_smiles
   import json

   # Load custom groups only
   with open('my_groups.json') as f:
       groups_data = json.load(f)
   groups = [PFASGroup(**g) for g in groups_data]

   # Load custom paths only
   with open('my_fpaths.json') as f:
       paths_data = json.load(f)
   paths = compile_smartsPaths(paths_data)

   # Use only custom definitions
   results = parse_smiles(smiles_list, pfas_groups=groups, smartsPaths=paths)

Filtering Defaults
^^^^^^^^^^^^^^^^^^

Use a subset of defaults:

.. code-block:: python

   from PFASgroups import get_PFASGroups

   # Get all groups
   all_groups = get_PFASGroups()

   # Filter to OECD groups only (IDs 1-28)
   oecd_groups = [g for g in all_groups if g.id <= 28]

   # Filter by pathway type
   perfluoro_groups = [g for g in all_groups if g.smartsPath == 'Perfluoroalkyl']

   # Filter by name pattern
   acid_groups = [g for g in all_groups if 'acid' in g.name.lower()]

Best Practices
--------------

1. **Use unique IDs**: Avoid ID conflicts with default groups (use IDs > 100)
2. **Test SMARTS patterns**: Validate patterns with RDKit before deployment
3. **Document constraints**: Include comments in JSON files explaining constraint logic
4. **Version control**: Keep configuration files under version control
5. **Validate changes**: Test custom groups against known molecules

Example: Complete Custom Configuration
--------------------------------------

.. code-block:: python

   from PFASgroups import (
       get_PFASGroups, get_smartsPaths, 
       PFASGroup, compile_smartsPath, 
       parse_smiles
   )

   # 1. Get and extend pathways
   paths = get_smartsPaths()
   paths['Perchlorinated'] = compile_smartsPath(
       "[C;X4](Cl)(Cl)!@!=!#[C;X4](Cl)(Cl)",
       "[C;X4](Cl)(Cl)Cl"
   )

   # 2. Get and extend groups
   groups = get_PFASGroups()
   
   # Add perchlorinated carboxylic acid
   groups.append(PFASGroup(
       id=101,
       name="Perchloroalkyl carboxylic acids",
       smarts1="[#6$([#6][#6](=O)([OH1,O-]))]",
       smarts2=None,
       smartsPath="Perchlorinated",
       constraints={
           "eq": {"O": 2},
           "only": ["C", "Cl", "O", "H"],
           "gte": {"Cl": 1}
       }
   ))

   # 3. Analyze molecules
   smiles_list = [
       "FC(F)(F)C(F)(F)C(=O)O",      # PFAS - matches PFCA
       "ClC(Cl)(Cl)C(Cl)(Cl)C(=O)O"  # Perchloro - matches custom group
   ]

   results = parse_smiles(
       smiles_list,
       pfas_groups=groups,
       smartsPaths=paths
   )

   for i, smiles in enumerate(smiles_list):
       print(f"\n{smiles}:")
       for group, count, lengths, chains in results[i]:
           print(f"  - {group.name}: {count} matches")
