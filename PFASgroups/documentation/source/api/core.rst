Core Module (core.py)
=====================

The core module provides the main functions for PFAS group identification and analysis.

.. module:: PFASgroups.core
   :synopsis: Core PFAS parsing and analysis functions

Main Parsing Functions
----------------------

parse_smiles
^^^^^^^^^^^^

.. py:function:: parse_smiles(smiles, bycomponent=False, output_format='list', **kwargs)

   Parse SMILES string(s) and return PFAS group information.

   :param smiles: Single SMILES string or list of SMILES strings
   :type smiles: str or list[str]
   :param bycomponent: Whether to use component-based analysis (for cyclic structures)
   :type bycomponent: bool, optional
   :param output_format: Output format - 'list', 'dataframe', or 'csv'
   :type output_format: str, optional
   :param kwargs: Additional parameters (pfas_groups, componentSmartss, etc.)
   :returns: Parsed results in specified format
   :rtype: list, pandas.DataFrame, or str

   **Example:**

   .. code-block:: python

      from PFASgroups import parse_smiles

      # Parse single SMILES
      results = parse_smiles("FC(F)(F)C(F)(F)C(=O)O")

      # Parse multiple SMILES
      results = parse_smiles(["SMILES1", "SMILES2"])

      # Get DataFrame output
      df = parse_smiles(smiles_list, output_format='dataframe')

parse_groups_in_mol
^^^^^^^^^^^^^^^^^^^

.. py:function:: parse_groups_in_mol(mol, bycomponent=False, **kwargs)

   Iterates over PFAS groups and finds the ones that match the molecule.

   :param mol: RDKit molecule object
   :type mol: rdkit.Chem.Mol
   :param bycomponent: Whether to use component-based analysis
   :type bycomponent: bool, optional
   :param kwargs: Additional parameters (formula, pfas_groups, componentSmartss, etc.)
   :returns: List of (PFASGroup, match_count, chain_lengths, matched_chains) tuples
   :rtype: list[tuple]

   **Processing Logic:**

   1. **Cyclic groups** (componentSmarts='cyclic'): Search for connected component matching smarts1
   2. **Both SMARTS defined**: Find shortest path between pairs of matches
   3. **Only smarts1, no constraints**: Use default smarts2 for each pathway type
   4. **Only smarts1 with constraints**: Substructure match with component analysis

parse_mols
^^^^^^^^^^

.. py:function:: parse_mols(mols, bycomponent=False, output_format='list', include_PFAS_definitions=True, **kwargs)

   Parse RDKit molecule(s) and return PFAS group information.

   :param mols: List of RDKit molecule objects
   :type mols: list[rdkit.Chem.Mol]
   :param bycomponent: Whether to use component-based analysis
   :type bycomponent: bool, optional
   :param output_format: Output format - 'list', 'dataframe', or 'csv'
   :type output_format: str, optional
   :param include_PFAS_definitions: Whether to include PFAS definition matches
   :type include_PFAS_definitions: bool, optional
   :param kwargs: Additional parameters
   :returns: Parsed results in specified format with comprehensive metrics including:
             - ``mean_branching``: Average branching (1.0=linear, 0.0=branched)
             - ``mean_component_fraction``: Average fraction per component
             - ``total_components_fraction``: Total union coverage
             - ``mean_eccentricity``, ``median_eccentricity``: Graph metrics
             - Graph structure metrics (diameter, radius, resistance)
             - Positional metrics (distance to barycenter, center, periphery)
             - Individual component details with all metrics
   :rtype: list, pandas.DataFrame, or str

Fingerprint Generation
----------------------

generate_fingerprint
^^^^^^^^^^^^^^^^^^^^

.. py:function:: generate_fingerprint(smiles, selected_groups=None, representation='vector', count_mode='binary', include_PFAS_definitions=False, **kwargs)

   Generate PFAS group fingerprints for machine learning and similarity analysis.
   
   Fingerprints are binary or count-based feature vectors where each dimension represents
   the presence/absence or abundance of a specific PFAS structural motif. These fingerprints
   enable quantitative structure-activity relationship (QSAR) modeling, chemical similarity
   assessment, and database clustering applications.

   :param smiles: SMILES string(s) to fingerprint
   :type smiles: str or list[str]
   :param selected_groups: Specific group IDs to include (default: all 55 groups)
   :type selected_groups: list[int] or range, optional
   :param representation: Output format - 'vector', 'dict', 'sparse', 'detailed', or 'int'
   :type representation: str, optional
   :param count_mode: Counting mode - 'binary', 'count', or 'max_component'
   :type count_mode: str, optional
   :param include_PFAS_definitions: Include PFAS definition matches (IDs 1-5)
   :type include_PFAS_definitions: bool, optional
   :returns: Tuple of (fingerprints, group_info)
   :rtype: tuple

   **Algorithm Overview:**

   1. **Structure Parsing**: Convert SMILES to RDKit molecule object
   2. **Group Detection**: Identify all matching PFAS groups via SMARTS pattern matching
   3. **Chain Analysis**: Quantify fluorinated chain components using graph algorithms
   4. **Vector Construction**: Encode detections as binary/count vectors

   **Representation Options:**

   - ``'vector'``: NumPy array with shape (n_molecules, n_groups) - suitable for ML
   - ``'dict'``: List of dictionaries mapping group names to values - human-readable
   - ``'sparse'``: List of dictionaries with only non-zero entries - memory-efficient
   - ``'detailed'``: Full match information including chain details and components
   - ``'int'``: Convert binary vector to integer representation

   **Count Mode Options:**

   - ``'binary'``: 1 if group present, 0 otherwise (default for classification)
   - ``'count'``: Number of matches for each group (captures multiplicity)
   - ``'max_component'``: Maximum chain length for each group (structural complexity)

   **Scientific Applications:**

   **Machine Learning**: Binary fingerprints serve as feature vectors for supervised learning
   models (random forests, neural networks, gradient boosting) to predict PFAS properties
   such as bioaccumulation potential, toxicity, or environmental persistence.

   **Similarity Analysis**: Pairwise Tanimoto coefficients computed from binary fingerprints
   quantify structural similarity between compounds. For two fingerprints **a** and **b**:

   .. math::

      T(a,b) = \\frac{|a \\cap b|}{|a \\cup b|} = \\frac{\\sum_i (a_i \\cdot b_i)}{\\sum_i a_i + \\sum_i b_i - \\sum_i (a_i \\cdot b_i)}

   **Database Clustering**: Fingerprints enable dimensionality reduction via PCA or t-SNE,
   followed by clustering (k-means, DBSCAN) to identify structural archetypes in PFAS
   inventories.

   **Distribution Analysis**: Frequency distributions of fingerprint patterns within chemical
   inventories can be compared using Kullback-Leibler divergence to assess compositional
   similarity between databases.

   **Example:**

   .. code-block:: python

      from PFASgroups import generate_fingerprint
      import numpy as np
      from sklearn.ensemble import RandomForestClassifier
      from sklearn.metrics.pairwise import cosine_similarity

      # Generate binary fingerprints for ML
      fps, info = generate_fingerprint(smiles_list)
      print(f"Fingerprint shape: {fps.shape}")  # (n_molecules, 55)

      # Train a classifier (example)
      X_train, X_test = fps[:800], fps[800:]
      y_train, y_test = labels[:800], labels[800:]
      clf = RandomForestClassifier(n_estimators=100)
      clf.fit(X_train, y_train)
      accuracy = clf.score(X_test, y_test)

      # Compute pairwise similarity matrix
      similarity_matrix = cosine_similarity(fps)
      
      # Count-based fingerprints for abundance analysis
      fps_count, info = generate_fingerprint(smiles_list, count_mode='count')

      # Select OECD groups only (groups 1-28)
      fps_oecd, info = generate_fingerprint(
          smiles_list, 
          selected_groups=range(0, 28)  # 0-indexed
      )
      
      # Dictionary format for inspection
      fps_dict, info = generate_fingerprint(
          ["FC(F)(F)C(F)(F)C(=O)O"],
          representation='dict'
      )
      print(fps_dict[0])  # {'PFCAs': 1, 'carboxylic acid': 1, ...}

   **Performance Considerations:**

   - Typical runtime: 10-50 ms per molecule (depends on size and complexity)
   - Memory: ~400 bytes per fingerprint (55 groups × 8 bytes/int64)
   - Parallelization: Use multiprocessing for large datasets (>10,000 molecules)

   **See Also:**

   - :doc:`../pfas_definitions` - PFAS definition matching
   - :doc:`../pfas_groups/oecd_groups` - OECD group classifications
   - :doc:`../user_guide` - Comprehensive usage examples

Visualization
-------------

plot_pfasgroups
^^^^^^^^^^^^^^^

.. py:function:: plot_pfasgroups(smiles, display=True, path=None, svg=False, ipython=False, subwidth=300, subheight=300, ncols=2, addAtomIndices=True, addBondIndices=False, paths=[0, 1, 2, 3], split_matches=False, SMARTS=None, **kwargs)

   Plot PFAS group assignments for a list of SMILES strings.

   :param smiles: SMILES string(s) to visualize
   :type smiles: str or list[str]
   :param display: Whether to display the plot
   :type display: bool, optional
   :param path: Path to save the plot image
   :type path: str, optional
   :param svg: Whether to generate SVG output
   :type svg: bool, optional
   :param ipython: Whether to display in IPython/Jupyter
   :type ipython: bool, optional
   :param subwidth: Width of each sub-image in pixels
   :type subwidth: int, optional
   :param subheight: Height of each sub-image in pixels
   :type subheight: int, optional
   :param ncols: Number of columns in grid layout
   :type ncols: int, optional
   :param addAtomIndices: Whether to add atom indices to the plot
   :type addAtomIndices: bool, optional
   :param addBondIndices: Whether to add bond indices to the plot
   :type addBondIndices: bool, optional
   :param paths: List of pathway indices to include
   :type paths: list[int], optional
   :param split_matches: Create separate images for each match
   :type split_matches: bool, optional
   :param SMARTS: Optional SMARTS pattern to highlight
   :type SMARTS: str, optional
   :returns: Image object (PIL Image or SVG)

Configuration Functions
-----------------------

get_componentSmartss
^^^^^^^^^^^^^^^

.. py:function:: get_componentSmartss(**kwargs)

   Get the default pathway SMARTS definitions.

   :param filename: Custom JSON file path (optional)
   :type filename: str, optional
   :returns: Dictionary mapping pathway names to [chain_mol, end_mol] pairs
   :rtype: dict

   **Default Pathways:**

   - ``Perfluoroalkyl``: Fully fluorinated chains (CF₂ units)
   - ``Polyfluoroalkyl``: Partially fluorinated chains (CHF units)
   - ``Polyfluoro``: Branched polyfluorinated chains
   - ``Polyfluorobr``: Brominated polyfluorinated chains

get_PFASGroups
^^^^^^^^^^^^^^

.. py:function:: get_PFASGroups(**kwargs)

   Get the default PFAS group definitions.

   :param filename: Custom JSON file path (optional)
   :type filename: str, optional
   :returns: List of PFASGroup objects
   :rtype: list[PFASGroup]

get_PFASDefinitions
^^^^^^^^^^^^^^^^^^^

.. py:function:: get_PFASDefinitions(**kwargs)

   Get the default PFAS definition objects.

   :param filename: Custom JSON file path (optional)
   :type filename: str, optional
   :returns: List of PFASDefinition objects
   :rtype: list[PFASDefinition]

compile_componentSmarts
^^^^^^^^^^^^^^^^^^

.. py:function:: compile_componentSmarts(chain_smarts, end_smarts)

   Compile a pair of SMARTS patterns into a ready-to-use path definition.

   :param chain_smarts: SMARTS pattern for the repeating chain unit
   :type chain_smarts: str
   :param end_smarts: SMARTS pattern for the terminal group
   :type end_smarts: str
   :returns: List containing [chain_mol, end_mol] preprocessed RDKit Mol objects
   :rtype: list

   **Example:**

   .. code-block:: python

      from PFASgroups import compile_componentSmarts, get_componentSmartss

      paths = get_componentSmartss()
      paths['Perchlorinated'] = compile_componentSmarts(
          "[C;X4](Cl)(Cl)!@!=!#[C;X4](Cl)(Cl)",
          "[C;X4](Cl)(Cl)Cl"
      )

compile_componentSmartss
^^^^^^^^^^^^^^^^^^^

.. py:function:: compile_componentSmartss(paths_dict)

   Compile multiple SMARTS path definitions from a dictionary.

   :param paths_dict: Dictionary with 'chain' and 'end' keys for each pathway
   :type paths_dict: dict
   :returns: Dictionary mapping path names to [chain_mol, end_mol] pairs
   :rtype: dict

Helper Functions
----------------

n_from_formula
^^^^^^^^^^^^^^

.. py:function:: n_from_formula(formula, element=None)

   Compute the number of elements from a molecular formula string.

   :param formula: Molecular formula string (e.g., "C8HF15O2")
   :type formula: str
   :param element: Specific element to count (None for all elements)
   :type element: str, optional
   :returns: Element count (int) or dictionary of all elements
   :rtype: int or dict

mol_to_nx
^^^^^^^^^

.. py:function:: mol_to_nx(mol)

   Construct a NetworkX graph from an RDKit molecule.

   :param mol: RDKit molecule object
   :type mol: rdkit.Chem.Mol
   :returns: NetworkX graph with atoms as nodes and bonds as edges
   :rtype: networkx.Graph

   Node attributes include ``element`` (atomic number) and ``symbol``.
   Edge attributes include ``order`` (bond order).

get_substruct
^^^^^^^^^^^^^

.. py:function:: get_substruct(mol, struct)

   Get indices of atoms matching a substructure.

   :param mol: RDKit molecule object
   :type mol: rdkit.Chem.Mol
   :param struct: SMARTS pattern (RDKit Mol object)
   :type struct: rdkit.Chem.Mol
   :returns: Set of matching atom indices
   :rtype: set[int]

remove_atoms
^^^^^^^^^^^^

.. py:function:: remove_atoms(mol, idxs, removable=['H', 'F', 'Cl', 'Br', 'I'], show_on_error=False)

   Remove atoms by indices and maintain molecular connectivity.

   :param mol: RDKit molecule object
   :type mol: rdkit.Chem.Mol
   :param idxs: Atom indices to remove
   :type idxs: list[int]
   :param removable: Elements that can be removed along with specified atoms
   :type removable: list[str], optional
   :param show_on_error: Display molecule on error
   :type show_on_error: bool, optional
   :returns: Modified molecule with atoms removed
   :rtype: rdkit.Chem.Mol

fragment_until_valence_is_correct
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. py:function:: fragment_until_valence_is_correct(mol, frags)

   Iteratively fragment a molecule until valence errors are resolved.

   :param mol: RDKit molecule object with potential valence issues
   :type mol: rdkit.Chem.Mol
   :param frags: Accumulator list for valid fragments
   :type frags: list
   :returns: List of valid molecular fragments
   :rtype: list[rdkit.Chem.Mol]
