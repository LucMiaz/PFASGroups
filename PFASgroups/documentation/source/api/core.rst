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
   :param kwargs: Additional parameters (pfas_groups, smartsPaths, etc.)
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
   :param kwargs: Additional parameters (formula, pfas_groups, smartsPaths, etc.)
   :returns: List of (PFASGroup, match_count, chain_lengths, matched_chains) tuples
   :rtype: list[tuple]

   **Processing Logic:**

   1. **Cyclic groups** (smartsPath='cyclic'): Search for connected component matching smarts1
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

   Generate PFAS group fingerprints for machine learning applications.

   :param smiles: SMILES string(s) to fingerprint
   :type smiles: str or list[str]
   :param selected_groups: Specific group IDs to include (default: all)
   :type selected_groups: list[int] or range, optional
   :param representation: Output format - 'vector', 'dict', or 'sparse'
   :type representation: str, optional
   :param count_mode: Counting mode - 'binary', 'count', or 'max_chain'
   :type count_mode: str, optional
   :param include_PFAS_definitions: Include PFAS definition matches
   :type include_PFAS_definitions: bool, optional
   :returns: Tuple of (fingerprints, group_info)
   :rtype: tuple

   **Representation Options:**

   - ``'vector'``: NumPy array with shape (n_molecules, n_groups)
   - ``'dict'``: List of dictionaries mapping group names to values
   - ``'sparse'``: List of dictionaries with only non-zero entries

   **Count Mode Options:**

   - ``'binary'``: 1 if group present, 0 otherwise
   - ``'count'``: Number of matches for each group
   - ``'max_chain'``: Maximum chain length for each group

   **Example:**

   .. code-block:: python

      from PFASgroups import generate_fingerprint

      # Binary fingerprints
      fps, info = generate_fingerprint(smiles_list)

      # Count-based fingerprints
      fps, info = generate_fingerprint(smiles_list, count_mode='count')

      # Select specific groups
      fps, info = generate_fingerprint(smiles_list, selected_groups=range(1, 29))

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

get_smartsPaths
^^^^^^^^^^^^^^^

.. py:function:: get_smartsPaths(**kwargs)

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

compile_smartsPath
^^^^^^^^^^^^^^^^^^

.. py:function:: compile_smartsPath(chain_smarts, end_smarts)

   Compile a pair of SMARTS patterns into a ready-to-use path definition.

   :param chain_smarts: SMARTS pattern for the repeating chain unit
   :type chain_smarts: str
   :param end_smarts: SMARTS pattern for the terminal group
   :type end_smarts: str
   :returns: List containing [chain_mol, end_mol] preprocessed RDKit Mol objects
   :rtype: list

   **Example:**

   .. code-block:: python

      from PFASgroups import compile_smartsPath, get_smartsPaths

      paths = get_smartsPaths()
      paths['Perchlorinated'] = compile_smartsPath(
          "[C;X4](Cl)(Cl)!@!=!#[C;X4](Cl)(Cl)",
          "[C;X4](Cl)(Cl)Cl"
      )

compile_smartsPaths
^^^^^^^^^^^^^^^^^^^

.. py:function:: compile_smartsPaths(paths_dict)

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
