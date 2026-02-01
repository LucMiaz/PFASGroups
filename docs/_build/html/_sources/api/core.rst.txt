Core API Reference
==================

This module contains the main functions for PFAS group detection and analysis.

Main Functions
--------------

parse_smiles
~~~~~~~~~~~~

.. autofunction:: PFASgroups.parser.parse_smiles

**Detailed Description:**

The primary function for detecting PFAS groups in molecules from SMILES strings.

**Parameters:**

- **smiles** (*str* or *list*): Single SMILES string or list of SMILES strings
- **bycomponent** (*bool*, optional): Use component-based analysis for cyclic structures (default: False)
- **output_format** (*str*, optional): Output format - 'list', 'dataframe', or 'csv' (default: 'list')
- **pfas_groups** (*list*, optional): Custom PFAS group definitions (default: uses built-in groups)
- **smartsPaths** (*dict*, optional): Custom pathway definitions (default: uses built-in paths)
- **formula** (*str*, optional): Pre-computed molecular formula (for optimization)
- **include_PFAS_definitions** (*bool*, optional): Include PFAS definition checks (default: False)

**Returns:**

- If output_format='list': List of tuples ``[(PFASGroup, match_count, chain_lengths, matched_chains), ...]``
- If output_format='dataframe': pandas DataFrame with columns: smiles, group_id, group_name, match_count, chain_lengths
- If output_format='csv': CSV string representation

**Examples:**

.. code-block:: python

   from PFASgroups import parse_smiles
   
   # Single molecule
   results = parse_smiles("FC(F)(F)C(F)(F)C(=O)O")
   for group, count, lengths, chains in results[0]:
       print(f"{group.name}: {count} matches, chains: {lengths}")
   
   # Multiple molecules with DataFrame output
   smiles_list = ["FC(F)(F)C(F)(F)C(=O)O", "FC(F)(F)C(F)(F)S(=O)(=O)O"]
   df = parse_smiles(smiles_list, output_format='dataframe')
   print(df)
   
   # Component-based analysis
   results = parse_smiles("C1(F)(F)C(F)(F)C(F)(F)C1(F)(F)", bycomponent=True)

parse_groups_in_mol
~~~~~~~~~~~~~~~~~~~

.. autofunction:: PFASgroups.parser.parse_groups_in_mol

**Detailed Description:**

Parse an RDKit molecule object directly, useful when you already have a molecule object or need to set specific properties.

**Parameters:**

- **mol** (*rdkit.Chem.Mol*): RDKit molecule object
- **bycomponent** (*bool*, optional): Use component-based analysis (default: False)
- **formula** (*str*, optional): Pre-computed molecular formula
- **pfas_groups** (*list*, optional): Custom PFAS group definitions
- **smartsPaths** (*dict*, optional): Custom pathway definitions

**Returns:**

List of tuples: ``[(PFASGroup, match_count, chain_lengths, matched_chains), ...]``

**Examples:**

.. code-block:: python

   from rdkit import Chem
   from PFASgroups import parse_groups_in_mol
   
   # Create molecule from SMILES
   mol = Chem.MolFromSmiles("FC(F)(F)C(F)(F)C(=O)O")
   
   # Parse groups
   results = parse_groups_in_mol(mol)
   
   # With pre-computed formula (optimization)
   from rdkit.Chem import Descriptors
   formula = Descriptors.rdMolDescriptors.CalcMolFormula(mol)
   results = parse_groups_in_mol(mol, formula=formula)

generate_fingerprint
~~~~~~~~~~~~~~~~~~~~

.. autofunction:: PFASgroups.fingerprints.generate_fingerprint

**Detailed Description:**

Generate PFAS group fingerprints for machine learning applications. Creates binary, count-based, or chain-length-based feature vectors.

**Parameters:**

- **smiles** (*str* or *list*): SMILES string(s) to fingerprint
- **selected_groups** (*list* or *range*, optional): Specific group IDs to include (default: all 55 groups)
- **representation** (*str*, optional): Output format - 'vector' (numpy array), 'dict' (dictionary), or 'sparse' (scipy sparse matrix) (default: 'vector')
- **count_mode** (*str*, optional): Fingerprint type:
  
  - 'binary': 1 if group present, 0 otherwise
  - 'count': Number of matches for each group
  - 'max_chain': Maximum chain length for each group
  
  (default: 'binary')

**Returns:**

Tuple: ``(fingerprints, group_info)``

- **fingerprints**: numpy array, list of dicts, or scipy sparse matrix
- **group_info**: Dictionary mapping group IDs to group names

**Examples:**

.. code-block:: python

   from PFASgroups import generate_fingerprint
   import numpy as np
   
   smiles_list = [
       "FC(F)(F)C(F)(F)C(=O)O",
       "FC(F)(F)C(F)(F)S(=O)(=O)O"
   ]
   
   # Binary fingerprints (default)
   fps, info = generate_fingerprint(smiles_list)
   print(f"Shape: {fps.shape}")  # (2, 55)
   print(f"Groups: {list(info.keys())}")
   
   # Dictionary representation (sparse)
   fps_dict, info = generate_fingerprint(smiles_list, representation='dict')
   print(fps_dict[0])  # {1: 1, 33: 1, ...}
   
   # Count-based
   fps_count, info = generate_fingerprint(smiles_list, count_mode='count')
   
   # Max chain length
   fps_chain, info = generate_fingerprint(smiles_list, count_mode='max_chain')
   
   # Select specific groups (OECD groups 1-28)
   fps_oecd, info = generate_fingerprint(smiles_list, selected_groups=range(1, 29))
   print(f"OECD fingerprint shape: {fps_oecd.shape}")  # (2, 28)

plot_pfasgroups
~~~~~~~~~~~~~~~

.. autofunction:: PFASgroups.draw_mols.plot_pfasgroups

**Detailed Description:**

Visualize PFAS group assignments on molecular structures with color-coded highlights.

**Parameters:**

- **smiles** (*str* or *list*): SMILES string(s) to visualize
- **display** (*bool*, optional): Whether to display the plot (default: True in Jupyter, False otherwise)
- **path** (*str*, optional): Path to save the image file
- **svg** (*bool*, optional): Generate SVG output (default: False, generates PNG)
- **subwidth** (*int*, optional): Width of each sub-image in pixels (default: 300)
- **subheight** (*int*, optional): Height of each sub-image in pixels (default: 300)
- **ncols** (*int*, optional): Number of columns in grid layout (default: 2)
- **addAtomIndices** (*bool*, optional): Add atom index numbers to visualization (default: False)
- **pfas_groups** (*list*, optional): Custom PFAS group definitions
- **smartsPaths** (*dict*, optional): Custom pathway definitions

**Returns:**

Image object (in Jupyter) or None

**Examples:**

.. code-block:: python

   from PFASgroups import plot_pfasgroups
   
   # Single molecule
   plot_pfasgroups("FC(F)(F)C(F)(F)C(=O)O")
   
   # Save to file
   plot_pfasgroups("FC(F)(F)C(F)(F)C(=O)O", path="pfba.png")
   
   # Generate SVG (vector graphics)
   plot_pfasgroups("FC(F)(F)C(F)(F)C(=O)O", svg=True, path="pfba.svg")
   
   # Multiple molecules in grid
   smiles_list = [
       "FC(F)(F)C(F)(F)C(=O)O",
       "FC(F)(F)C(F)(F)S(=O)(=O)O",
       "C(C(F)(F)C(F)(F)F)O"
   ]
   plot_pfasgroups(smiles_list, ncols=3, subwidth=400)
   
   # With atom indices (for debugging)
   plot_pfasgroups("FC(F)(F)C(F)(F)C(=O)O", addAtomIndices=True)

Helper Functions
----------------

get_PFASGroups
~~~~~~~~~~~~~~

.. autofunction:: PFASgroups.getter.get_PFASGroups

**Description:**

Returns the list of all 55 default PFAS group definitions.

**Parameters:**

- **custom_groups_file** (*str*, optional): Path to custom JSON file with group definitions

**Returns:**

List of PFASGroup objects

**Examples:**

.. code-block:: python

   from PFASgroups import get_PFASGroups
   
   # Get default groups
   groups = get_PFASGroups()
   print(f"Total groups: {len(groups)}")
   
   # Print OECD groups (1-28)
   for group in groups[:28]:
       print(f"{group.id}: {group.name}")
   
   # Load custom groups
   custom_groups = get_PFASGroups(custom_groups_file='my_groups.json')

get_smartsPaths
~~~~~~~~~~~~~~~

.. autofunction:: PFASgroups.getter.get_smartsPaths

**Description:**

Returns the default pathway SMARTS definitions for perfluorinated and polyfluorinated chains.

**Parameters:**

- **custom_paths_file** (*str*, optional): Path to custom JSON file with pathway definitions

**Returns:**

Dictionary with keys: 'Perfluoroalkyl', 'Polyfluoroalkyl', 'Polyfluoro', 'Polyfluorobr'

**Examples:**

.. code-block:: python

   from PFASgroups import get_smartsPaths
   
   # Get default paths
   paths = get_smartsPaths()
   print(paths.keys())  # dict_keys(['Perfluoroalkyl', 'Polyfluoroalkyl', ...])
   
   # Inspect perfluoroalkyl pathway
   print(paths['Perfluoroalkyl'])
   # Output: {'chain': <rdkit.Chem.rdchem.Mol>, 'end': <rdkit.Chem.rdchem.Mol>}

compile_smartsPath
~~~~~~~~~~~~~~~~~~

.. autofunction:: PFASgroups.parser.compile_smartsPath

**Description:**

Compile SMARTS patterns for custom pathway definitions.

**Parameters:**

- **chain_smarts** (*str*): SMARTS pattern for chain segments
- **end_smarts** (*str*): SMARTS pattern for chain terminals

**Returns:**

Dictionary: ``{'chain': <Mol>, 'end': <Mol>}``

**Examples:**

.. code-block:: python

   from PFASgroups import compile_smartsPath
   
   # Create perchlorinated pathway (analogous to perfluorinated)
   perchlorinated = compile_smartsPath(
       "[C;X4](Cl)(Cl)!@!=!#[C;X4](Cl)(Cl)",  # Chain segment
       "[C;X4](Cl)(Cl)Cl"  # Terminal
   )
   
   print(perchlorinated)
   # {'chain': <rdkit.Chem.rdchem.Mol>, 'end': <rdkit.Chem.rdchem.Mol>}

compile_smartsPaths
~~~~~~~~~~~~~~~~~~~

.. autofunction:: PFASgroups.parser.compile_smartsPaths

**Description:**

Compile multiple pathway definitions from a dictionary.

**Parameters:**

- **paths_dict** (*dict*): Dictionary with pathway names as keys and {'chain': str, 'end': str} as values

**Returns:**

Dictionary of compiled pathways

**Examples:**

.. code-block:: python

   from PFASgroups import compile_smartsPaths
   
   # Define custom pathways
   custom_paths = {
       'Perchlorinated': {
           'chain': '[C;X4](Cl)(Cl)!@!=!#[C;X4](Cl)(Cl)',
           'end': '[C;X4](Cl)(Cl)Cl'
       },
       'Perbrominated': {
           'chain': '[C;X4](Br)(Br)!@!=!#[C;X4](Br)(Br)',
           'end': '[C;X4](Br)(Br)Br'
       }
   }
   
   # Compile all pathways
   compiled_paths = compile_smartsPaths(custom_paths)
   print(compiled_paths.keys())

Internal Functions
------------------

These functions are used internally but can be accessed for advanced use cases:

- ``get_default_PFAS_groups()``: Load default group definitions from JSON
- ``get_default_smarts_paths()``: Load default pathway definitions from JSON
- ``check_formula_constraints(formula_dict, constraints)``: Validate molecular formula against constraints
- ``find_fluorinated_chains(mol, source_atoms, pathway)``: Graph-based chain finding algorithm

See Also
--------

- :doc:`models`: PFAS group and definition data models
- :doc:`../algorithm`: Detailed algorithm description
- :doc:`../customization`: Creating custom groups and pathways
