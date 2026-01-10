# PFASgroups Function Documentation

This document provides comprehensive documentation for all functions in the PFASgroups package.

## Core Module (`core.py`)

### Component Detection and Analysis

#### `ComponentsSolver` class
**Purpose**: Manages fluorinated component detection and metrics computation in molecules.

**Key Methods**:
- `__init__(mol, **kwargs)`: Initializes with molecule, computes all component types and precomputes full component sizes
- `get(pathType, max_dist=0)`: Retrieves components of specified type, optionally extended by distance
- `extend_components(max_dist)`: Extends components to include nearby atoms within graph distance
- `get_augmented_component(pathType, max_dist, component_index, smarts_matches)`: Returns original component augmented with shortest path atoms connecting SMARTS matches
- `get_full_component_atoms(component)`: Gets all atoms in a component including attached H, F, halogens (not just carbon backbone)
- `get_total_components_fraction(matched_components_list)`: Calculates fraction of molecule covered by the union of all components
- `compute_component_metrics(component)`: Calculates graph metrics (diameter, radius, center, periphery, barycenter)
- `compute_smarts_component_metrics(component, smarts_matches)`: Calculates SMARTS-specific metrics including centrality and positional metrics
- `get_matched_component_dict(component, smarts_matches, smarts_type)`: Returns comprehensive dictionary with all metrics for a matched component

**Component Fraction Metrics**:
- `component_fraction`: Fraction of molecule covered by each component (includes all atoms: C backbone + attached H, F, Cl, Br, I + SMARTS matches)
- `total_components_fraction`: Fraction of molecule covered by union of all components (accounts for overlaps)

**Usage**:
```python
with ComponentsSolver(mol) as solver:
    perfluoro_components = solver.get('Perfluoroalkyl', max_dist=0)
    extended = solver.get('Perfluoroalkyl', max_dist=2)
    
    # Get full component with all attached atoms
    full_atoms = solver.get_full_component_atoms(component)
    
    # Calculate total coverage across all matched components
    total_fraction = solver.get_total_components_fraction(matched_components)
```

#### `find_alkyl_components(mol, smarts, components, pathType, max_dist, component_solver)`
**Purpose**: Find fluorinated alkyl components containing SMARTS pattern matches.

**Parameters**:
- `mol`: RDKit molecule
- `smarts`: SMARTS pattern for functional group
- `components`: List of component sets to search
- `pathType`: Type of component path ('Perfluoroalkyl', etc.)
- `max_dist`: Extension distance (0 = no extension)
- `component_solver`: ComponentsSolver instance

**Returns**: `(max_length, component_sizes, num_matches, matched_components_list)`

**Behavior with max_dist > 0**:
- Matches SMARTS against extended components
- Returns original component + shortest path atoms connecting SMARTS to original
- Enables detection of functional groups near (but not in) fluorinated chains

#### `find_aryl_components(mol, aryl_smarts, component_solver)`
**Purpose**: Find cyclic/aromatic fluorinated components.

**Returns**: Same format as find_alkyl_components

#### `find_components_between_smarts(mol, smarts1, smarts2, component_solver)`
**Purpose**: Find components containing both smarts1 and smarts2 matches (bridging groups).

**Use Cases**:
- Detecting molecules with functional groups at both ends
- Finding linkers between two specific groups
- When smarts2=None, behaves like single SMARTS matching

### Main Classification Functions

#### `parse_groups_in_mol(mol, **kwargs)`
**Purpose**: Main entry point for PFAS group classification.

**Process**:
1. Adds explicit hydrogens
2. Validates molecular formula
3. Fragments if valence errors exist
4. Iterates through all PFAS groups checking:
   - Formula constraints
   - SMARTS pattern matches
   - Component-based matching
5. Returns list of (PFASGroup, match_count, component_sizes, matched_components)

**Component Matching Logic by Group Type**:
- `smartsPath='cyclic'`: Find aromatic/cyclic components with smarts1
- `smarts1 + smarts2`: Find components containing both patterns
- `smarts1 only`: Find components containing smarts1 matches
- `smartsPath only`: Find all components of that type

#### `parse_smiles(smiles, bycomponent=False, output_format='list')`
**Purpose**: Parse SMILES string(s) to identify PFAS groups.

**Parameters**:
- `smiles`: Single SMILES or list
- `bycomponent`: Use component-based analysis
- `output_format`: 'list', 'dataframe', or 'csv'

**Returns**: Depends on output_format

#### `parse_mols(mols, output_format='list', include_PFAS_definitions=True)`
**Purpose**: Parse RDKit molecules for PFAS groups and definitions.

**Returns**: List of dictionaries with:
- Molecule identifiers (SMILES, InChI, InChIKey)
- Matched PFAS groups with comprehensive metrics including:
  - `mean_branching`: Average branching across all components
  - `mean_component_fraction`: Average fraction of molecule covered per component
  - `total_components_fraction`: Total fraction covered by union of all components
  - `mean_eccentricity`: Average mean eccentricity across components
  - `median_eccentricity`: Average median eccentricity across components
  - Graph metrics (diameter, radius, resistance)
  - Positional metrics (distance to barycenter, center, periphery)
  - Individual component details with all metrics
- Matched PFAS definitions (if enabled)

### Helper Functions

#### `mol_to_nx(mol)`
**Purpose**: Convert RDKit molecule to NetworkX graph for graph-theoretic analysis.

#### `remove_atoms(mol, idxs, removable=['H','F','Cl','Br','I'])`
**Purpose**: Remove specified atoms while maintaining molecular connectivity.

**Behavior**:
- Removes target atoms and their removable neighbors
- Reconnects remaining structure by finding shortest paths
- Preserves formal charges
- Handles branching and multiple removal points

#### `n_from_formula(formula, element=None)`
**Purpose**: Parse molecular formula string into element count dictionary.

**Examples**:
```python
n_from_formula("C8F17O2")  # {'C': 8, 'F': 17, 'O': 2}
n_from_formula("C8F17O2", element='F')  # 17
```

#### `fragment_on_bond(mol, a1, a2)`
**Purpose**: Fragment molecule by breaking bond between atoms a1 and a2.

#### `fragment_until_valence_is_correct(mol, frags)`
**Purpose**: Iteratively fragment molecule to resolve valence errors.

### Fingerprint Generation

#### `generate_fingerprint(smiles, selected_groups=None, representation='vector', count_mode='binary')`
**Purpose**: Generate PFAS group-based fingerprints for molecules.

**Parameters**:
- `selected_groups`: Indices of groups to include (None = all)
- `representation`: 'vector', 'dict', 'sparse', 'detailed', 'int'
- `count_mode`: 'binary', 'count', 'max_chain'

**Returns**: `(fingerprints, group_info)`

**Use Cases**:
- Machine learning features
- Similarity calculations
- Database searching

### Visualization

#### `plot_pfasgroups(smiles, **kwargs)`
**Purpose**: Visualize molecules with highlighted PFAS group matches.

**Parameters**:
- `svg`: Generate SVG (True) or raster (False)
- `subwidth`, `subheight`: Image dimensions
- `addAtomIndices`: Show atom indices
- `paths`: Which group types to highlight
- `split_matches`: Separate image per match

---

## Generate Molecule Module (`generate_mol.py`)

### Chain Generation

#### `generate_random_carbon_chain(n, cycle=False, alkene=False, alkyne=False)`
**Purpose**: Create random carbon backbone with n carbons.

**Parameters**:
- `cycle`: Start with benzene ring
- `alkene`: Allow C=C bonds (50% probability)
- `alkyne`: Allow C≡C bonds (50% probability)

**Algorithm**:
1. Start with single C or benzene
2. Add carbons one at a time
3. Randomly select attachment point (respecting valence)
4. Randomly assign bond type based on flags
5. Sanitize final structure

### Fluorination

#### `fluorinate_mol(mol, perfluorinated=True, p=0.3, phigh=1)`
**Purpose**: Replace hydrogen atoms with fluorine.

**Parameters**:
- `perfluorinated`: Full fluorination if True
- `p`: Minimum replacement probability (polyfluoro mode)
- `phigh`: Maximum replacement probability

**Probability Adaptation**:
- Formula: `max(p, min(phigh, 1 - nF/nH))`
- Adapts as fluorination progresses
- Ensures gradual tapering in polyfluoro mode

### Functional Group Attachment

#### `get_attachment(mol, m, atom_symbols=['C'], neighbors_symbols={'C':['F','H']})`
**Purpose**: Find m suitable attachment points for functional groups.

**Returns**: List of (atom_idx, neighbor_idx) tuples

#### `append_functional_group(mol, group_smiles, insertion='attach', m=1, atom_indices=None)`
**Purpose**: Add functional groups to molecule.

**Insertion Modes**:
- `'attach'`: Replace F/H neighbor with functional group
- `'insert'`: Break bond and insert functional group between atoms

#### `attach_mol(mol, submol, atom_index)`
**Purpose**: Attach submolecule by replacing F/H neighbor.

#### `insert_mol(mol, group_mol, atom_index, neighbor_index)`
**Purpose**: Insert submolecule between two bonded atoms.

#### `append_functional_groups(mol, functional_groups, **kwargs)`
**Purpose**: Add multiple functional groups with complex specifications.

**Functional Group Specification**:
```python
{
    'group_smiles': 'C(=O)O',  # SMILES string
    'n': 1 or '[1,3]',  # Count or range
    'mode': 'attach' or 'insert',
    'neighbours': ['C'],
    'items': [  # Optional sub-components
        {'smiles': 'O', 'n': '[1,5]'}
    ]
}
```

#### `generate_random_mol(n, functional_groups, perfluorinated=True, **kwargs)`
**Purpose**: Complete pipeline for PFAS molecule generation.

**Pipeline**:
1. Generate carbon chain
2. Fluorinate
3. Add functional groups

---

## Generate Homologues Module (`generate_homologues.py`)

#### `find_chain(mol, pathsmarts, endsmarts, repeating='C(F)(F)')`
**Purpose**: Find fluorinated chains with repeating units between end groups.

**Returns**: List of chain dictionaries with:
- `chain`: All atom indices in chain
- `start`, `end`: Terminal atom indices
- `belly`: Repeating unit atom indices
- `consecutive_parts`: Partitioned belly sections

**Algorithm**:
1. Find all pairs of end group matches
2. Calculate shortest paths between pairs
3. Filter paths that are fully fluorinated
4. Identify repeating unit sections
5. Remove chains that are subsets of longer chains

#### `generate_homologues(mol, smartsPathName='Perfluoroalkyl', repeating='C(F)F', base_repeating=['C'])`
**Purpose**: Generate homologous series by removing repeating units.

**Process**:
1. Find fluorinated chains with repeating units
2. Extract consecutive segments of repeating units
3. Generate all combinations of partial removals
4. Remove atoms while preserving connectivity
5. Return dictionary of {InChIKey: {formula: molecule}}

**Use Cases**:
- Degradation product prediction
- Homologous series enumeration
- Systematic structure exploration

**Example**:
```python
# From C8 PFOA, generate C7, C6, C5, C4, C3, C2 homologues
homologues = generate_homologues(pfoa_mol, repeating='C(F)F')
```

---

## Fragmentation Module (`fragmentation.py`)

### Bond Dissociation Energy (BDE)

#### `yield_scheme(name='BDE')`
**Purpose**: Decorator that provides bond dissociation energy data.

**Provides**: Dictionary of {atom1_Z: {atom2_Z: BDE_kJ/mol}}

#### `generate_fragments(mol, bondWeights=None, nb_breakingpoints=5)`
**Purpose**: Fragment molecule by breaking weakest bonds.

**Algorithm**:
1. Calculate BDE for each bond
2. Invert to get break probabilities (weaker = higher)
3. Randomly select bonds weighted by break probability
4. Fragment on selected bonds

**Returns**: List of fragment molecules

#### `fragment(mol, bde_dict, max_depth=4, nb_breakingpoints=5)`
**Purpose**: Iteratively fragment molecule to generate product tree.

**Process**:
1. Fragment parent molecule
2. Fragment each resulting fragment
3. Continue for max_depth iterations
4. Collect unique fragments by formula

**Returns**: `{InChIKey: {formula: {'mol': mol, 'count': int, 'mz': float}}}`

#### `fragment_to_distribution(mol, bde_dict, max_depth=4, nb_breakingpoints=5)`
**Purpose**: Generate mass spectrum-like distribution from fragmentation.

**Returns**: `(mz_array, intensity_array)` sorted by m/z

### PFAS-Specific Fragmentation

#### `find_fluorinated_chains(mol, smartsPathName='Perfluoroalkyl')`
**Purpose**: Identify atom indices in fluorinated chains.

**Returns**: Set of atom indices that are part of fluorinated moieties

#### `get_non_fluorinated_bonds(mol, fluorinated_atoms=None)`
**Purpose**: Get bonds not part of fluorinated chain.

**Use**: Identify where non-PFAS portions can be cleaved

#### `get_bonds_connected_to_fluorinated_path(mol, fluorinated_atoms=None)`
**Purpose**: Get bonds connecting to (but not in) fluorinated chain.

**Use**: Identify fragmentation sites for degradation products

#### `generate_degradation_products(mol, smartsPathName='Perfluoroalkyl', max_breaks=3)`
**Purpose**: Generate realistic PFAS degradation products.

**Algorithm**:
1. Identify fluorinated chain
2. Find bonds connected to chain
3. Generate all combinations up to max_breaks
4. Fragment and collect products

**Preserves**: Fluorinated chain intact (realistic degradation)

#### `generate_systematic_degradation_products(mol, preserve_fluorinated_chain=True, max_combinations=1000)`
**Purpose**: Comprehensive degradation product enumeration.

**Returns**: `{n_breaks: {InChIKey: {formula: {...}}}}`

#### `analyze_degradation_pathways(mol, max_breaks=3)`
**Purpose**: Analyze and report on degradation possibilities.

**Returns**: Dictionary with:
- `original_molecule`: Parent molecule info
- `degradation_summary`: Statistics by break count
- `products`: Detailed list of all products

---

## Drawing Module (`draw_mols.py`)

#### `merge_raster(imgs, buffer, ncols)`
**Purpose**: Arrange PIL Images in grid layout.

**Returns**: `(merged_image, total_width, total_height)`

#### `merge_svg(imgs, buffer, ncols)`
**Purpose**: Arrange SVG figures in grid layout.

**Returns**: `(merged_svg, total_width, total_height)`

#### `draw_images(imgs, buffer=1, ncols=2, svg=False)`
**Purpose**: Dispatcher for merge_raster or merge_svg.

#### `plot_mol(mol, **kwargs)`
**Purpose**: Draw single molecule with customization.

**Options**:
- `svg`: Output format
- `subwidth`, `subheight`: Dimensions
- `fixedBondLength`: Bond length in pixels
- `addAtomIndices`, `addBondIndices`: Index labeling
- `bondLineWidth`: Line thickness
- `maxFontSize`, `minFontSize`: Font size range

#### `plot_mols(mols, **kwargs)`
**Purpose**: Draw multiple molecules in grid.

**Same options as plot_mol**, plus:
- `ncols`: Grid columns
- `buffer`: Spacing between molecules

---

## Model Classes

### PFASGroup (`PFASGroupModel.py`)

**Purpose**: Represents a specific PFAS functional group.

**Attributes**:
- `id`, `name`: Identifiers
- `smarts1`, `smarts2`: Functional group patterns
- `smartsPath`: Component type to search
- `max_dist_from_CF`: Extension distance
- `constraints`: Formula constraints

**Constraint Methods**:
- `constraint_gte(formula_dict)`: Minimum element counts
- `constraint_lte(formula_dict)`: Maximum element counts
- `constraint_eq(formula_dict)`: Exact element counts
- `constraint_only(formula_dict)`: Allowed elements only
- `constraint_rel(formula_dict)`: Relational constraints (e.g., O = C/2)
- `formula_dict_satisfies_constraints(formula_dict)`: Check all constraints

**Example Constraint**:
```python
constraints = {
    'only': ['C', 'F', 'O', 'H'],
    'gte': {'C': 2, 'F': 3},
    'eq': {'O': 2},
    'rel': {'C': {'atoms': ['F'], 'div': 2, 'add': 0}}
}
```

### PFASDefinition (`PFASDefinitionModel.py`)

**Purpose**: Represents broad PFAS definition criteria.

**Attributes**:
- `id`, `name`, `description`: Identifiers
- `smarts_strings`, `smarts_patterns`: Structural patterns
- `fluorineRatio`: Minimum F ratio threshold
- `includeHydrogen`: Include H in ratio calculation
- `requireBoth`: AND vs OR logic for SMARTS + ratio

**Methods**:
- `applies_to_molecule(mol_or_smiles)`: Check if definition applies
- `_compute_formula(mol)`: Calculate formula dictionary
- `_check_fluorine_ratio(formula)`: Validate F ratio

**Logic Modes**:
- `requireBoth=False`: SMARTS OR fluorine ratio (permissive)
- `requireBoth=True`: SMARTS AND fluorine ratio (strict)

---

## Usage Examples

### Classify a Molecule
```python
from PFASgroups import parse_smiles

results = parse_smiles("FC(F)(F)C(F)(F)C(F)(F)C(F)(F)C(=O)O")
for result in results:
    for match in result['matches']:
        print(f"{match['group_name']}: {match['match_count']} matches")
```

### Generate Synthetic PFAS
```python
from PFASgroups.generate_mol import generate_random_mol

mol = generate_random_mol(
    n=8,
    functional_groups={'group_smiles': 'C(=O)O', 'n': 1, 'mode': 'attach'},
    perfluorinated=True
)
```

### Find Degradation Products
```python
from PFASgroups.fragmentation import generate_degradation_products

products = generate_degradation_products(mol, max_breaks=2)
print(f"Generated {len(products)} degradation products")
```

### Generate Homologous Series
```python
from PFASgroups.generate_homologues import generate_homologues

homologues = generate_homologues(mol, smartsPathName='Perfluoroalkyl')
for inchikey, formulas in homologues.items():
    for formula, mol in formulas.items():
        print(f"{formula}: {Chem.MolToSmiles(mol)}")
```

---

## Key Concepts

### Component-Based Matching
- **Components**: Connected fluorinated subgraphs in molecules
- **Component Types**: Perfluoroalkyl, Polyfluoroalkyl, Polyfluoro, Polyfluorobr, cyclic
- **Extension**: Components can be extended by graph distance to find nearby functional groups
- **Augmentation**: When extended components match, original component + connecting path is returned

### Graph Metrics
- **Diameter**: Maximum eccentricity (longest shortest path)
- **Radius**: Minimum eccentricity
- **Center**: Nodes with minimum eccentricity
- **Periphery**: Nodes with maximum eccentricity
- **Barycenter**: Nodes minimizing total distance to all others
- **Branching**: Measure of linearity (1.0 = linear, 0.0 = branched) - renamed from "eccentricity" for clarity
- **Mean Eccentricity**: Average graph-theoretic eccentricity across all nodes in component
- **Median Eccentricity**: Median graph-theoretic eccentricity across all nodes in component
- **SMARTS Centrality**: How central functional group is (1.0 = center, 0.0 = peripheral)
- **Component Fraction**: Fraction of total molecule atoms covered by this component (includes all attached H, F, Cl, Br, I)
- **Total Components Fraction**: Fraction of total molecule atoms covered by union of all components for a PFAS group (summary metric)

### Constraint System
Molecular formula constraints support:
- **Only**: Exclusive element list
- **GTE/LTE**: Min/max element counts
- **EQ**: Exact element counts
- **REL**: Relational formulas between elements

### Molecular Identifiers
- **SMILES**: Human-readable structure notation
- **InChI**: Canonical structure identifier
- **InChIKey**: Hashed InChI for database keys
- **Formula**: Molecular formula (e.g., C8F17O2)

---

## Performance Tips

1. **Component Caching**: ComponentsSolver caches metrics - reuse same instance
2. **Parallel Processing**: Use `parse_mols()` for batch processing
3. **max_dist_from_CF**: Keep low (0-2) unless specifically needed
4. **Fingerprints**: Use 'vector' or 'sparse' for large-scale work
5. **Fragment Limits**: Set `max_depth` and `max_combinations` conservatively

## Error Handling

Common issues and solutions:

1. **Valence Errors**: Molecule sanitization fails
   - Solution: Code auto-fragments to resolve
   
2. **No Components Found**: Empty component list
   - Check: Is molecule actually fluorinated?
   - Check: Is smartsPath correct for fluorination type?

3. **Fragment Disconnection**: Removal breaks molecule
   - Solution: Use `remove_atoms()` which maintains connectivity
   
4. **Invalid SMARTS**: Group fails to match
   - Check: SMARTS pattern is valid RDKit syntax
   - Check: Pattern actually exists in molecule

---

This documentation covers the complete PFASgroups API. For more details on specific functions, refer to the inline docstrings in each module.
