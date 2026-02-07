Data Models API Reference
=========================

This module defines the data structures used to represent PFAS groups and definitions.

PFASGroup Class
---------------

.. autoclass:: PFASgroups.PFASGroupModel.PFASGroup
   :members:
   :special-members: __init__

**Class Description:**

Represents a PFAS functional group definition with SMARTS patterns, pathway requirements, and molecular formula constraints.

**Attributes:**

- **id** (*int*): Unique identifier for the group (1-55+ for OECD groups, higher for custom)
- **name** (*str*): Descriptive name (e.g., "Perfluoroalkyl carboxylic acids")
- **alias** (*str*, optional): Short name or abbreviation (e.g., "PFCAs")
- **smarts1** (*rdkit.Chem.Mol*): Primary SMARTS pattern (compiled RDKit Mol object)
- **smarts2** (*rdkit.Chem.Mol*, optional): Secondary SMARTS pattern for additional specificity
- **componentSmarts** (*str*): Pathway type - 'Perfluoroalkyl', 'Polyfluoroalkyl', 'cyclic', or None
- **constraints** (*dict*): Molecular formula constraints
- **max_dist_from_CF** (*int*): Maximum bond distance from fluorinated carbon terminal (default: 2)
- **base_functional_groups** (*list*, optional): Base functional group categories
- **main_group** (*str*, optional): Main classification category

**Constraint Dictionary Structure:**

.. code-block:: python

   constraints = {
       "eq": {"O": 2},  # Exactly 2 oxygen atoms
       "gte": {"F": 1},  # At least 1 fluorine atom
       "lte": {"N": 1},  # At most 1 nitrogen atom
       "only": ["C", "F", "O", "H"],  # Only these elements allowed
       "rel": {  # Relative element ratios
           "C": {
               "atoms": ["F", "H"],  # Reference atoms
               "add": 0.5,  # Addition term
               "div": 2  # Division factor
           }
       }
   }

**Relative Constraint Formula:**

.. math::

   n_{target} = \frac{\sum_{i} n_{atoms_i} + add}{div}

For example, perfluorinated compounds: C = (F + H) / 2 - 1

**Methods:**

formula_dict_satisfies_constraints
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   def formula_dict_satisfies_constraints(self, formula_dict: dict) -> bool:
       """
       Check if a molecular formula satisfies all constraints.
       
       Args:
           formula_dict: Dictionary mapping element symbols to counts
                        e.g., {'C': 8, 'F': 15, 'O': 2, 'H': 1}
       
       Returns:
           True if all constraints satisfied, False otherwise
       """

**Examples:**

Creating a custom PFAS group:

.. code-block:: python

   from PFASgroups import PFASGroup
   from rdkit import Chem
   
   # Define a custom PFAS group
   my_group = PFASGroup(
       id=999,
       name="Perfluoroalkyl nitrates",
       alias="PFANs",
       smarts1=Chem.MolFromSmarts("[C]-[O]-[N+](=O)[O-]"),
       smarts2=None,
       componentSmarts="Perfluoroalkyl",
       constraints={
           "eq": {"N": 1, "O": 3},  # One N and three O (nitrate)
           "gte": {"F": 1},  # At least one F
           "only": ["C", "F", "O", "N", "H"]
       },
       max_dist_from_CF=2
   )
   
   # Check if a formula satisfies constraints
   formula = {'C': 4, 'F': 7, 'O': 3, 'N': 1, 'H': 0}
   satisfies = my_group.formula_dict_satisfies_constraints(formula)
   print(f"Formula satisfies: {satisfies}")

Using in custom analysis:

.. code-block:: python

   from PFASgroups import get_PFASGroups, parse_smiles
   
   # Get default groups
   groups = get_PFASGroups()
   
   # Add custom group
   groups.append(my_group)
   
   # Use in parsing
   results = parse_smiles(
       "FC(F)(F)C(F)(F)C(F)(F)ON(=O)=O",
       pfas_groups=groups
   )

PFASDefinition Class
--------------------

.. autoclass:: PFASgroups.PFASDefinitionModel.PFASDefinition
   :members:
   :special-members: __init__

**Class Description:**

Represents broader PFAS definitions (e.g., OECD, EPA, EU) based on SMARTS patterns and fluorine ratios.

**Attributes:**

- **id** (*int*): Unique identifier
- **name** (*str*): Definition name (e.g., "OECD PFAS Definition")
- **smarts** (*list*): List of SMARTS pattern strings
- **fluorineRatio** (*float*): Minimum F/(C+F+H) ratio threshold (0.0-1.0)
- **description** (*str*): Human-readable description
- **includeHydrogen** (*bool*): Include H in ratio calculation
- **requireBoth** (*bool*): Require both SMARTS match AND fluorine ratio

**Methods:**

applies_to_molecule
~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   def applies_to_molecule(self, mol_or_smiles, formula=None, **kwargs) -> bool:
       """
       Check if this definition applies to a molecule.
       
       Args:
           mol_or_smiles: RDKit Mol object or SMILES string
           formula: Pre-computed molecular formula (optional)
           **kwargs: Additional parameters
       
       Returns:
           True if definition applies, False otherwise
       """

**Examples:**

Creating a custom PFAS definition:

.. code-block:: python

   from PFASgroups import PFASDefinition
   
   # OECD definition example
   oecd_def = PFASDefinition(
       id=1,
       name="OECD PFAS Definition",
       smarts=[
           "[C](F)(F)F",  # CF3 group
           "[C](F)(F)[C](F)(F)"  # CF2-CF2 segment
       ],
       fluorineRatio=0.0,  # No minimum ratio (SMARTS match sufficient)
       description="OECD: Contains at least one perfluoroalkyl moiety",
       includeHydrogen=True,
       requireBoth=False  # SMARTS match OR ratio
   )
   
   # Check if molecule matches definition
   from rdkit import Chem
   mol = Chem.MolFromSmiles("FC(F)(F)C(F)(F)C(=O)O")
   matches = oecd_def.applies_to_molecule(mol)
   print(f"Matches OECD definition: {matches}")

Using multiple definitions:

.. code-block:: python

   # Define multiple PFAS definitions
   definitions = [
       PFASDefinition(
           id=1, name="OECD",
           smarts=["[C](F)(F)F"],
           fluorineRatio=0.0,
           requireBoth=False
       ),
       PFASDefinition(
           id=2, name="High Fluorination",
           smarts=["[F]"],
           fluorineRatio=0.5,  # 50% F content
           includeHydrogen=True,
           requireBoth=True  # Must have F AND meet ratio
       )
   ]
   
   # Test molecule against all definitions
   smiles = "FC(F)(F)C(F)(F)C(F)(F)C(=O)O"
   mol = Chem.MolFromSmiles(smiles)
   
   for definition in definitions:
       matches = definition.applies_to_molecule(mol)
       print(f"{definition.name}: {matches}")

Data Model Integration
----------------------

Both PFASGroup and PFASDefinition integrate with the core parsing functions:

.. code-block:: python

   from PFASgroups import (
       parse_smiles, 
       get_PFASGroups,
       PFASGroup,
       PFASDefinition
   )
   from rdkit import Chem
   
   # Custom groups
   custom_groups = get_PFASGroups()
   custom_groups.append(PFASGroup(
       id=100,
       name="My Custom Group",
       smarts1=Chem.MolFromSmarts("[custom_pattern]"),
       componentSmarts="Perfluoroalkyl",
       constraints={"gte": {"F": 1}}
   ))
   
   # Use in parsing
   results = parse_smiles(
       "FC(F)(F)C(F)(F)C(=O)O",
       pfas_groups=custom_groups
   )
   
   # Custom definitions for broad classification
   oecd_def = PFASDefinition(
       id=1, name="OECD",
       smarts=["[C](F)(F)F"],
       fluorineRatio=0.0
   )
   
   # Check definition
   mol = Chem.MolFromSmiles("FC(F)(F)C(F)(F)C(=O)O")
   is_pfas = oecd_def.applies_to_molecule(mol)

JSON Configuration
------------------

Both models support JSON serialization for configuration files:

**PFAS Groups JSON (PFAS_groups_smarts.json):**

.. code-block:: json

   [
     {
       "id": 1,
       "name": "Perfluoroalkyl carboxylic acids",
       "alias": "PFCAs",
       "smarts1": "[#6$([#6][#6](=O)([OH1,O-]))]",
       "smarts2": null,
       "componentSmarts": "Perfluoroalkyl",
       "constraints": {
         "eq": {"O": 2},
         "gte": {"F": 1},
         "only": ["C", "F", "O", "H"]
       },
       "max_dist_from_CF": 2
     }
   ]

**PFAS Definitions JSON:**

.. code-block:: json

   [
     {
       "id": 1,
       "name": "OECD PFAS Definition",
       "smarts": ["[C](F)(F)F", "[C](F)(F)[C](F)(F)"],
       "fluorineRatio": 0.0,
       "description": "OECD definition",
       "includeHydrogen": true,
       "requireBoth": false
     }
   ]

See Also
--------

- :doc:`core`: Core parsing functions
- :doc:`../customization`: Creating custom groups and definitions
- :doc:`../algorithm`: How these models are used in detection
