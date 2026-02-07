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
- **componentSmartss** (*dict*, optional): Custom pathway definitions (default: uses built-in paths)
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
- **componentSmartss** (*dict*, optional): Custom pathway definitions

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

Generate PFAS group fingerprints for machine learning and chemoinformatics applications.
Fingerprints encode the presence, abundance, or structural characteristics of PFAS groups
as numerical feature vectors suitable for computational analysis.

**Scientific Background:**

Molecular fingerprints are binary or count-based vectors where each dimension represents
a specific structural feature. PFAS fingerprints differ from traditional molecular fingerprints
(e.g., Morgan, MACCS keys) by encoding domain-specific PFAS structural motifs rather than
generic chemical features. This enables:

- **Classification models**: Predict PFAS properties (toxicity, bioaccumulation, persistence)
- **Similarity analysis**: Quantify structural similarity using Tanimoto coefficients
- **Database clustering**: Group chemicals by structural archetypes
- **Inventory comparison**: Assess compositional differences using Kullback-Leibler divergence

**Algorithm:**

1. Parse each SMILES string to RDKit molecule object
2. Apply SMARTS pattern matching for each of 55 PFAS groups
3. Count matches or quantify chain lengths based on ``count_mode``
4. Encode results as binary (0/1) or count-based integer vector
5. Return array with dimensions (n_molecules, n_groups)

**Parameters:**

- **smiles** (*str* or *list*): SMILES string(s) to fingerprint
- **selected_groups** (*list* or *range*, optional): Specific group IDs to include (default: all 55 groups)
  
  - Groups 1-28: OECD-defined PFAS classes
  - Groups 29-55: Generic functional group classifications
  - Example: ``range(1, 29)`` for OECD groups only

- **representation** (*str*, optional): Output format:
  
  - ``'vector'``: NumPy array (n_molecules × n_groups) - **recommended for ML**
  - ``'dict'``: List of dictionaries mapping group names to counts - human-readable
  - ``'sparse'``: List of dictionaries with only non-zero entries - memory-efficient
  - ``'detailed'``: Full match information including component details
  - ``'int'``: Convert binary vector to integer representation
  
  (default: 'vector')

- **count_mode** (*str*, optional): Encoding scheme:
  
  - ``'binary'``: 1 if group present, 0 otherwise - **standard for classification**
  - ``'count'``: Integer count of matches - captures multiplicity
  - ``'max_component'``: Maximum fluorinated chain size - structural complexity metric
  
  (default: 'binary')

- **include_PFAS_definitions** (*bool*, optional): Include regulatory definition matches (IDs 1-5)
- **pfas_groups** (*list*, optional): Custom group definitions (overrides defaults)

**Returns:**

Tuple: ``(fingerprints, group_info)``

- **fingerprints**: 
  
  - If representation='vector': NumPy array with shape (n_molecules, n_groups)
  - If representation='dict': List of dictionaries
  - If representation='sparse': List of dictionaries with non-zero entries only

- **group_info**: Dictionary with keys:
  
  - ``'group_names'``: List of group names in order
  - ``'group_ids'``: List of group IDs
  - ``'selected_indices'``: Indices of selected groups
  - ``'total_groups'``: Total number of groups

**Performance:**

- Typical runtime: 10-50 ms per molecule
- Memory: ~440 bytes per fingerprint (55 groups × 8 bytes)
- Scales linearly with number of molecules
- Consider multiprocessing for >10,000 molecules

**Examples:**

.. code-block:: python

   from PFASgroups import generate_fingerprint
   import numpy as np
   from sklearn.ensemble import RandomForestClassifier
   from scipy.spatial.distance import jaccard
   
   smiles_list = [
       "FC(F)(F)C(F)(F)C(=O)O",      # PFBA
       "FC(F)(F)C(F)(F)S(=O)(=O)O",  # PFBS  
       "C(C(F)(F)C(F)(F)F)O"         # FTOH
   ]
   
   # Binary fingerprints (default) - for classification
   fps, info = generate_fingerprint(smiles_list)
   print(f"Shape: {fps.shape}")  # (3, 55)
   print(f"Data type: {fps.dtype}")  # int64
   print(f"Groups detected: {np.sum(fps, axis=1)}")  # [2 2 1]
   
   # Dictionary representation - human-readable
   fps_dict, info = generate_fingerprint(smiles_list, representation='dict')
   print("PFBA groups:", {k: v for k, v in fps_dict[0].items() if v > 0})
   # Output: {'PFCAs': 1, 'carboxylic acid': 1}
   
   # Count-based fingerprints - capture multiplicity
   fps_count, info = generate_fingerprint(smiles_list, count_mode='count')
   print(f"PFBA match counts: {fps_count[0]}")
   
   # Maximum chain length encoding
   fps_chain, info = generate_fingerprint(smiles_list, count_mode='max_component')
   print(f"Max chain lengths: {fps_chain}")
   
   # Select OECD groups only (groups 1-28, but 0-indexed)
   fps_oecd, info = generate_fingerprint(
       smiles_list, 
       selected_groups=range(0, 28)
   )
   print(f"OECD fingerprint shape: {fps_oecd.shape}")  # (3, 28)
   print(f"OECD group names: {info['group_names'][:5]}")
   
   # Machine learning application
   # Assuming you have labels (e.g., toxic=1, non-toxic=0)
   X_train, y_train = fps[:100], labels[:100]
   X_test, y_test = fps[100:], labels[100:]
   
   clf = RandomForestClassifier(n_estimators=100, random_state=42)
   clf.fit(X_train, y_train)
   accuracy = clf.score(X_test, y_test)
   print(f"Classification accuracy: {accuracy:.3f}")
   
   # Feature importance
   importances = clf.feature_importances_
   top_groups = np.argsort(importances)[-5:]
   print("Most predictive groups:")
   for idx in top_groups:
       print(f"  {info['group_names'][idx]}: {importances[idx]:.3f}")
   
   # Similarity analysis
   # Jaccard similarity between PFBA and PFBS
   sim = 1 - jaccard(fps[0], fps[1])
   print(f"PFBA-PFBS similarity: {sim:.3f}")
   
   # Tanimoto coefficient (for binary fingerprints)
   def tanimoto(a, b):
       c = np.sum(a * b)
       return c / (np.sum(a) + np.sum(b) - c)
   
   tanimoto_sim = tanimoto(fps[0], fps[1])
   print(f"PFBA-PFBS Tanimoto: {tanimoto_sim:.3f}")

**Use Cases:**

1. **QSAR Modeling**: Train regression or classification models to predict PFAS properties

   .. code-block:: python
   
      from sklearn.linear_model import Ridge
      from sklearn.model_selection import cross_val_score
      
      # Predict log Kow (octanol-water partition coefficient)
      fps, _ = generate_fingerprint(smiles_list)
      scores = cross_val_score(Ridge(), fps, log_kow, cv=5)
      print(f"R² = {np.mean(scores):.3f}")

2. **Database Screening**: Identify PFAS subclasses in large inventories

   .. code-block:: python
   
      # Screen 100,000 compounds
      fps, info = generate_fingerprint(large_smiles_list, representation='sparse')
      
      # Find molecules with perfluoroalkyl carboxylic acids
      pfca_idx = info['group_names'].index('PFCAs')
      pfcas = [i for i, fp in enumerate(fps) if pfca_idx in fp]

3. **Similarity-Based Grouping**: Cluster structurally similar PFAS

   .. code-block:: python
   
      from sklearn.cluster import KMeans
      from sklearn.decomposition import PCA
      
      fps, _ = generate_fingerprint(smiles_list)
      
      # Dimensionality reduction
      pca = PCA(n_components=10)
      fps_reduced = pca.fit_transform(fps)
      
      # Clustering
      kmeans = KMeans(n_clusters=5)
      labels = kmeans.fit_predict(fps_reduced)

4. **Regulatory Classification**: Check multiple PFAS definitions

   .. code-block:: python
   
      fps, info = generate_fingerprint(
          smiles_list,
          include_PFAS_definitions=True,
          selected_groups=range(0, 5)  # Definition IDs 1-5
      )
      
      definition_names = ['OECD', 'EU', 'OPPT', 'UK', 'PFASTRUCTv5']
      for i, smiles in enumerate(smiles_list):
          matches = [definition_names[j] for j in range(5) if fps[i, j] == 1]
          print(f"{smiles}: {', '.join(matches)}")

**See Also:**

- :doc:`../pfas_definitions` - Regulatory PFAS definitions
- :doc:`../tutorial` - Extended tutorial with examples
- :doc:`models` - PFASGroup and PFASDefinition models

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
- **componentSmartss** (*dict*, optional): Custom pathway definitions

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

get_componentSmartss
~~~~~~~~~~~~~~~

.. autofunction:: PFASgroups.getter.get_componentSmartss

**Description:**

Returns the default pathway SMARTS definitions for perfluorinated and polyfluorinated chains.

**Parameters:**

- **custom_paths_file** (*str*, optional): Path to custom JSON file with pathway definitions

**Returns:**

Dictionary with keys: 'Perfluoroalkyl', 'Polyfluoroalkyl', 'Polyfluoro', 'Polyfluorobr'

**Examples:**

.. code-block:: python

   from PFASgroups import get_componentSmartss
   
   # Get default paths
   paths = get_componentSmartss()
   print(paths.keys())  # dict_keys(['Perfluoroalkyl', 'Polyfluoroalkyl', ...])
   
   # Inspect perfluoroalkyl pathway
   print(paths['Perfluoroalkyl'])
   # Output: {'component': <rdkit.Chem.rdchem.Mol>, 'end': <rdkit.Chem.rdchem.Mol>}

compile_componentSmarts
~~~~~~~~~~~~~~~~~~

.. autofunction:: PFASgroups.parser.compile_componentSmarts

**Description:**

Compile SMARTS patterns for custom pathway definitions.

**Parameters:**

- **chain_smarts** (*str*): SMARTS pattern for chain segments
- **end_smarts** (*str*): SMARTS pattern for chain terminals

**Returns:**

Dictionary: ``{'component': <Mol>, 'end': <Mol>}``

**Examples:**

.. code-block:: python

   from PFASgroups import compile_componentSmarts
   
   # Create perchlorinated pathway (analogous to perfluorinated)
   perchlorinated = compile_componentSmarts(
       "[C;X4](Cl)(Cl)!@!=!#[C;X4](Cl)(Cl)",  # Chain segment
       "[C;X4](Cl)(Cl)Cl"  # Terminal
   )
   
   print(perchlorinated)
   # {'component': <rdkit.Chem.rdchem.Mol>, 'end': <rdkit.Chem.rdchem.Mol>}

compile_componentSmartss
~~~~~~~~~~~~~~~~~~~

.. autofunction:: PFASgroups.parser.compile_componentSmartss

**Description:**

Compile multiple pathway definitions from a dictionary.

**Parameters:**

- **paths_dict** (*dict*): Dictionary with pathway names as keys and {'component': str, 'end': str} as values

**Returns:**

Dictionary of compiled pathways

**Examples:**

.. code-block:: python

   from PFASgroups import compile_componentSmartss
   
   # Define custom pathways
   custom_paths = {
       'Perchlorinated': {
           'component': '[C;X4](Cl)(Cl)!@!=!#[C;X4](Cl)(Cl)',
           'end': '[C;X4](Cl)(Cl)Cl'
       },
       'Perbrominated': {
           'component': '[C;X4](Br)(Br)!@!=!#[C;X4](Br)(Br)',
           'end': '[C;X4](Br)(Br)Br'
       }
   }
   
   # Compile all pathways
   compiled_paths = compile_componentSmartss(custom_paths)
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
