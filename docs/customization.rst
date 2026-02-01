Customization
=============

PFASgroups provides extensive customization options for defining custom PFAS groups, pathways, and definitions.

Custom PFAS Groups
------------------

Creating Groups from Scratch
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Define a completely new PFAS group:

.. code-block:: python

   from PFASgroups import PFASGroup, parse_smiles
   from rdkit import Chem
   
   # Define a custom group for perfluoroalkyl nitrates
   pfan_group = PFASGroup(
       id=100,
       name="Perfluoroalkyl nitrates",
       alias="PFANs",
       smarts1=Chem.MolFromSmarts("[C]-[O]-[N+](=O)[O-]"),
       smarts2=None,  # Optional secondary pattern
       smartsPath="Perfluoroalkyl",  # Require perfluorinated chain
       constraints={
           "eq": {"N": 1, "O": 3},  # Exactly 1 N and 3 O
           "gte": {"F": 1},  # At least 1 fluorine
           "only": ["C", "F", "O", "N", "H"]  # Only these elements
       },
       max_dist_from_CF=2  # Max distance from fluorinated carbon
   )
   
   # Test the custom group
   test_smiles = "FC(F)(F)C(F)(F)C(F)(F)ON(=O)=O"
   results = parse_smiles(test_smiles, pfas_groups=[pfan_group])
   
   if results[0]:
       for group, count, lengths, _ in results[0]:
           print(f"✓ Detected: {group.name}")
           print(f"  Chain length: {lengths}")

Modifying Existing Groups
~~~~~~~~~~~~~~~~~~~~~~~~~~

Customize existing PFAS groups:

.. code-block:: python

   from PFASgroups import get_PFASGroups, parse_smiles
   from rdkit import Chem
   
   # Get default groups
   groups = get_PFASGroups()
   
   # Find and modify PFCA group (ID 1)
   pfca = next(g for g in groups if g.id == 1)
   
   # Make it less restrictive (allow more elements)
   pfca.constraints["only"] = ["C", "F", "O", "H", "N", "S"]
   
   # Or change max distance
   pfca.max_dist_from_CF = 3  # Allow groups farther from chain
   
   # Use modified groups
   results = parse_smiles("FC(F)(F)C(F)(F)C(=O)O", pfas_groups=groups)

Loading Groups from JSON
~~~~~~~~~~~~~~~~~~~~~~~~~

Store custom groups in a JSON file:

**my_groups.json:**

.. code-block:: json

   [
     {
       "id": 100,
       "name": "Perfluoroalkyl nitrates",
       "alias": "PFANs",
       "smarts1": "[C]-[O]-[N+](=O)[O-]",
       "smarts2": null,
       "smartsPath": "Perfluoroalkyl",
       "constraints": {
         "eq": {"N": 1, "O": 3},
         "gte": {"F": 1},
         "only": ["C", "F", "O", "N", "H"]
       },
       "max_dist_from_CF": 2
     },
     {
       "id": 101,
       "name": "Perfluoroalkyl thiols",
       "alias": "PFATs",
       "smarts1": "[C]-[SH1]",
       "smartsPath": "Perfluoroalkyl",
       "constraints": {
         "eq": {"S": 1},
         "gte": {"F": 1}
       }
     }
   ]

**Load and use:**

.. code-block:: python

   from PFASgroups import get_PFASGroups, parse_smiles
   
   # Load custom groups
   custom_groups = get_PFASGroups(custom_groups_file='my_groups.json')
   
   # Use in parsing
   results = parse_smiles(
       "FC(F)(F)C(F)(F)C(F)(F)ON(=O)=O",
       pfas_groups=custom_groups
   )

Constraint Examples
~~~~~~~~~~~~~~~~~~~

Different types of molecular constraints:

**Exact count (eq):**

.. code-block:: python

   constraints = {
       "eq": {"O": 2}  # Exactly 2 oxygens (e.g., carboxylic acid)
   }

**Minimum count (gte):**

.. code-block:: python

   constraints = {
       "gte": {"F": 3}  # At least 3 fluorines
   }

**Maximum count (lte):**

.. code-block:: python

   constraints = {
       "lte": {"N": 1}  # At most 1 nitrogen
   }

**Element restrictions (only):**

.. code-block:: python

   constraints = {
       "only": ["C", "F", "O", "H"]  # No other elements allowed
   }

**Relative ratios (rel):**

.. code-block:: python

   # Perfluorinated: C = (F + H) / 2 - 1
   constraints = {
       "rel": {
           "C": {
               "atoms": ["F", "H"],  # Sum F + H
               "add": -1,  # Subtract 1
               "div": 2  # Divide by 2
           }
       }
   }

Custom Pathways
---------------

Defining Pathway Types
~~~~~~~~~~~~~~~~~~~~~~~

Create custom fluorinated pathways:

.. code-block:: python

   from PFASgroups import compile_smartsPath, get_smartsPaths, parse_smiles
   
   # Get default pathways
   paths = get_smartsPaths()
   
   # Add perchlorinated pathway (Cl instead of F)
   paths['Perchlorinated'] = compile_smartsPath(
       "[C;X4](Cl)(Cl)!@!=!#[C;X4](Cl)(Cl)",  # Chain segment
       "[C;X4](Cl)(Cl)Cl"  # Terminal CCl3
   )
   
   # Add perbrominated pathway
   paths['Perbrominated'] = compile_smartsPath(
       "[C;X4](Br)(Br)!@!=!#[C;X4](Br)(Br)",
       "[C;X4](Br)(Br)Br"
   )
   
   # Use with custom groups
   from PFASgroups import PFASGroup
   from rdkit import Chem
   
   chlorinated_acid = PFASGroup(
       id=200,
       name="Perchlorinated carboxylic acids",
       smarts1=Chem.MolFromSmarts("[#6$([#6][#6](=O)([OH1,O-]))]"),
       smartsPath="Perchlorinated",
       constraints={"eq": {"O": 2}, "gte": {"Cl": 1}}
   )
   
   # Test on chlorinated compound
   smiles = "ClC(Cl)(Cl)C(Cl)(Cl)C(=O)O"
   results = parse_smiles(
       smiles,
       pfas_groups=[chlorinated_acid],
       smartsPaths=paths
   )

Loading Pathways from JSON
~~~~~~~~~~~~~~~~~~~~~~~~~~~

**fpaths.json:**

.. code-block:: json

   {
     "Perfluoroalkyl": {
       "chain": "[C;X4;H0](F)(F)!@!=!#[C;X4;H0](F)(F)",
       "end": "[C;X4;H0](F)(F)F"
     },
     "Perchlorinated": {
       "chain": "[C;X4](Cl)(Cl)!@!=!#[C;X4](Cl)(Cl)",
       "end": "[C;X4](Cl)(Cl)Cl"
     }
   }

**Load and use:**

.. code-block:: python

   from PFASgroups import compile_smartsPaths, parse_smiles
   import json
   
   # Load custom pathways
   with open('fpaths.json', 'r') as f:
       paths_dict = json.load(f)
   
   paths = compile_smartsPaths(paths_dict)
   
   # Use in parsing
   results = parse_smiles(smiles, smartsPaths=paths)

SMARTS Pattern Guide
~~~~~~~~~~~~~~~~~~~~

Key SMARTS syntax for pathway definitions:

.. code-block:: text

   [C;X4;H0]   - Carbon with 4 connections, no hydrogens
   (F)(F)      - Two fluorine substituents
   !@          - Non-ring bond
   !=          - Not double bond
   !#          - Not triple bond

**Examples:**

- Perfluoroalkyl chain: ``[C;X4;H0](F)(F)``
- Polyfluoroalkyl chain: ``[C;X4;H1](F)`` (allows one H)
- Terminal CF₃: ``[C;X4;H0](F)(F)F``
- Terminal CHF₂: ``[C;X4;H1](F)F``

Custom PFAS Definitions
------------------------

Creating Definitions
~~~~~~~~~~~~~~~~~~~~

Define broad PFAS classifications:

.. code-block:: python

   from PFASgroups import PFASDefinition
   from rdkit import Chem
   
   # OECD-style definition
   oecd_pfas = PFASDefinition(
       id=1,
       name="OECD PFAS",
       smarts=["[C](F)(F)F", "[C](F)(F)[C](F)(F)"],
       fluorineRatio=0.0,  # No minimum ratio
       description="Contains perfluoroalkyl moiety",
       includeHydrogen=True,
       requireBoth=False  # SMARTS OR ratio
   )
   
   # High fluorination definition
   high_f = PFASDefinition(
       id=2,
       name="High Fluorination",
       smarts=["[F]"],  # Must have F
       fluorineRatio=0.5,  # 50% F content
       description="Highly fluorinated compounds",
       includeHydrogen=True,
       requireBoth=True  # SMARTS AND ratio
   )
   
   # Test molecules
   test_compounds = [
       "FC(F)(F)C(F)(F)C(=O)O",  # PFBA
       "FC(F)(F)C(=O)O",  # TFA
       "CCO"  # Ethanol
   ]
   
   for smiles in test_compounds:
       mol = Chem.MolFromSmiles(smiles)
       print(f"\n{smiles}:")
       print(f"  OECD: {oecd_pfas.applies_to_molecule(mol)}")
       print(f"  High-F: {high_f.applies_to_molecule(mol)}")

Fluorine Ratio Calculation
~~~~~~~~~~~~~~~~~~~~~~~~~~~

Understanding fluorine ratio thresholds:

.. code-block:: python

   from rdkit import Chem
   
   def calculate_f_ratio(smiles, include_h=True):
       """Calculate fluorine ratio for a molecule."""
       mol = Chem.MolFromSmiles(smiles)
       if mol is None:
           return None
       
       # Count atoms
       n_c = sum(1 for atom in mol.GetAtoms() if atom.GetSymbol() == 'C')
       n_f = sum(1 for atom in mol.GetAtoms() if atom.GetSymbol() == 'F')
       n_h = sum(atom.GetTotalNumHs() for atom in mol.GetAtoms())
       
       # Calculate ratio
       if include_h:
           denominator = n_c + n_f + n_h
       else:
           denominator = n_c + n_f
       
       return n_f / denominator if denominator > 0 else 0
   
   # Test on various compounds
   compounds = [
       ("PFOA", "FC(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(=O)O"),
       ("TFA", "FC(F)(F)C(=O)O"),
       ("Benzene", "c1ccccc1")
   ]
   
   for name, smiles in compounds:
       ratio = calculate_f_ratio(smiles)
       print(f"{name}: {ratio:.2%}")

Advanced Customization
-----------------------

Combining Custom Components
~~~~~~~~~~~~~~~~~~~~~~~~~~~

Use custom groups, pathways, and definitions together:

.. code-block:: python

   from PFASgroups import (
       PFASGroup, PFASDefinition,
       compile_smartsPath, get_PFASGroups,
       parse_smiles
   )
   from rdkit import Chem
   
   # 1. Custom pathway
   paths = {
       'PartiallyFluorinated': compile_smartsPath(
           "[C;X4](F)!@!=!#[C;X4]",  # One F per carbon
           "[C;X4](F)F"  # Terminal
       )
   }
   
   # 2. Custom group using custom pathway
   partial_f_acid = PFASGroup(
       id=300,
       name="Partially fluorinated carboxylic acids",
       smarts1=Chem.MolFromSmarts("[C]-[C](=O)O"),
       smartsPath="PartiallyFluorinated",
       constraints={"eq": {"O": 2}, "gte": {"F": 1}}
   )
   
   # 3. Custom definition
   partial_f_def = PFASDefinition(
       id=10,
       name="Partial Fluorination",
       smarts=["[C][F]"],
       fluorineRatio=0.2,  # 20% minimum
       requireBoth=True
   )
   
   # Use all together
   smiles = "FC(C(F)C(F)C(=O)O)"
   
   results = parse_smiles(
       smiles,
       pfas_groups=[partial_f_acid],
       smartsPaths=paths
   )
   
   mol = Chem.MolFromSmiles(smiles)
   matches_def = partial_f_def.applies_to_molecule(mol)
   
   print(f"Groups detected: {len(results[0]) if results[0] else 0}")
   print(f"Matches definition: {matches_def}")

Configuration Management
~~~~~~~~~~~~~~~~~~~~~~~~

Manage multiple configuration sets:

.. code-block:: python

   import json
   from PFASgroups import compile_smartsPaths, get_PFASGroups
   
   class PFASConfig:
       """Manage PFAS analysis configuration."""
       
       def __init__(self, name):
           self.name = name
           self.groups = []
           self.paths = {}
           self.definitions = []
       
       def load_groups(self, json_file):
           """Load groups from JSON file."""
           self.groups = get_PFASGroups(custom_groups_file=json_file)
           return self
       
       def load_paths(self, json_file):
           """Load pathways from JSON file."""
           with open(json_file, 'r') as f:
               paths_dict = json.load(f)
           self.paths = compile_smartsPaths(paths_dict)
           return self
       
       def save_config(self, output_file):
           """Save configuration to JSON."""
           config = {
               'name': self.name,
               'groups_count': len(self.groups),
               'paths': list(self.paths.keys())
           }
           with open(output_file, 'w') as f:
               json.dump(config, f, indent=2)
       
       def parse(self, smiles):
           """Parse SMILES with this configuration."""
           from PFASgroups import parse_smiles
           return parse_smiles(
               smiles,
               pfas_groups=self.groups,
               smartsPaths=self.paths
           )
   
   # Use it
   config = PFASConfig('regulatory')
   config.load_groups('regulatory_groups.json')
   config.load_paths('regulatory_paths.json')
   config.save_config('regulatory_config.json')
   
   results = config.parse("FC(F)(F)C(F)(F)C(=O)O")

Testing Custom Groups
~~~~~~~~~~~~~~~~~~~~~~

Validate custom PFAS groups:

.. code-block:: python

   from PFASgroups import PFASGroup, parse_smiles
   from rdkit import Chem
   
   def test_pfas_group(group, test_cases):
       """Test a PFAS group on multiple test cases."""
       print(f"\nTesting: {group.name}")
       print("=" * 60)
       
       results = []
       for smiles, expected in test_cases:
           result = parse_smiles(smiles, pfas_groups=[group])
           detected = len(result[0]) > 0 if result[0] else False
           
           status = "✓" if detected == expected else "✗"
           results.append((smiles, expected, detected, status))
           
           print(f"{status} {smiles}")
           print(f"   Expected: {expected}, Got: {detected}")
       
       # Summary
       passed = sum(1 for _, exp, det, _ in results if exp == det)
       print(f"\nPassed: {passed}/{len(results)}")
       return results
   
   # Example usage
   my_group = PFASGroup(
       id=100,
       name="Test Group",
       smarts1=Chem.MolFromSmarts("[C](F)(F)F"),
       smartsPath="Perfluoroalkyl",
       constraints={"gte": {"F": 3}}
   )
   
   test_cases = [
       ("FC(F)(F)C(F)(F)C(=O)O", True),  # Should match
       ("CCO", False),  # Should not match
       ("FC(F)(F)CF", True)  # Should match
   ]
   
   test_pfas_group(my_group, test_cases)

Best Practices
--------------

SMARTS Pattern Design
~~~~~~~~~~~~~~~~~~~~~

Tips for effective SMARTS patterns:

1. **Be specific but not too restrictive:**

   .. code-block:: python
   
      # Too general
      smarts = "[C][F]"  # Matches any C-F bond
      
      # Too specific
      smarts = "[CH2][CF2][CF3]"  # Only exact pattern
      
      # Good balance
      smarts = "[C;X4](F)(F)F"  # CF3 with flexibility

2. **Use atom environment specifications:**

   .. code-block:: python
   
      "[C;X4]"  # Carbon with 4 connections
      "[C;H0]"  # Carbon with no hydrogens
      "[C;!R]"  # Non-ring carbon

3. **Test thoroughly:**

   .. code-block:: python
   
      from rdkit import Chem
      
      pattern = Chem.MolFromSmarts("[your_pattern]")
      test_mol = Chem.MolFromSmiles("test_smiles")
      matches = test_mol.GetSubstructMatches(pattern)
      print(f"Matches: {len(matches)}")

Constraint Design
~~~~~~~~~~~~~~~~~

Guidelines for molecular constraints:

1. **Start with essential constraints:**

   .. code-block:: python
   
      # Minimum viable constraints
      constraints = {
           "gte": {"F": 1}  # Must have fluorine
      }

2. **Add specificity gradually:**

   .. code-block:: python
   
      # More specific
      constraints = {
           "gte": {"F": 1},
           "eq": {"O": 2},  # Carboxylic acid
           "only": ["C", "F", "O", "H"]
      }

3. **Use relative constraints for chain validation:**

   .. code-block:: python
   
      # Perfluorinated validation
      constraints = {
           "rel": {
               "C": {"atoms": ["F"], "add": 0.5, "div": 2}
           }
      }

Documentation
~~~~~~~~~~~~~

Document your custom groups:

.. code-block:: python

   my_group = PFASGroup(
       id=100,
       name="My Custom Group",
       alias="MCG",
       smarts1=pattern,
       smartsPath="Perfluoroalkyl",
       constraints=constraints,
       # Add description in comments
       # Purpose: Detect specific fluorinated compounds
       # Created: 2026-02-01
       # Author: Your Name
       # Test cases: See tests/test_custom.py
   )

See Also
--------

- :doc:`api/core`: Core API reference
- :doc:`api/models`: Data model details
- :doc:`tutorial`: Practical examples
- :doc:`algorithm`: Algorithm explanation
