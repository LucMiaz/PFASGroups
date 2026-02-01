from rdkit import Chem
import numpy as np
import re
import networkx as nx
from .core import mol_to_nx
import logging




class PFASGroup():
    """Model class representing a specific PFAS functional group with structural patterns.
    
    A PFASGroup defines a specific fluorinated functional group using SMARTS patterns,
    component path types, and molecular formula constraints. Groups are used to classify
    molecules into specific categories (e.g., "Perfluoroalkyl carboxylic acid").
    
    Attributes
    ----------
    id : int
        Unique identifier for this PFAS group
    name : str
        Human-readable group name (e.g., "Perfluoroalkyl carboxylic acid")
    smarts : Chem.Mol or None
        SMARTS patterns (compiled RDKit molecule) for functional group detection.
        None if group is defined by smartsPath alone.
    smartsPath : str or None
        Type of fluorinated component to search:
        - 'Perfluoroalkyl': Fully fluorinated carbon chains
        - 'Polyfluoroalkyl': Partially fluorinated carbon chains  
        - 'Polyfluoro': General polyfluorinated structures
        - 'Polyfluorobr': Polyfluorobromo structures
        - 'cyclic': Aromatic/cyclic fluorinated structures
        - None: Check all component types
    max_dist_from_CF : int
        Maximum graph distance (number of bonds) from fluorinated component to functional group.
        When > 0, extends component search radius to find nearby functional groups.
    linker_smarts : Chem.Mol or None
        Compiled SMARTS pattern for validating linker atoms between fluorinated component
        and functional group. When None (default), no restriction is applied to linker atoms.
        Only used when max_dist_from_CF > 0.
    constraints : dict
        Molecular formula constraints with keys:
        - 'only': Elements that must be present exclusively (e.g., ['C', 'F', 'O'])
        - 'gte': Minimum element counts (e.g., {'C': 2})
        - 'lte': Maximum element counts (e.g., {'O': 2})
        - 'eq': Exact element counts (e.g., {'N': 1})
        - 'rel': Relational constraints (e.g., {'O': {'atoms': ['C'], 'div': 2, 'add': 0}})
    
    Examples
    --------
    >>> # Perfluoroalkyl carboxylic acid: R_F-COOH
    >>> pfaa = PFASGroup(
    ...     id=1,
    ...     name="Perfluoroalkyl carboxylic acid",
    ...     smarts={"C(=O)O":1},  # Carboxylic acid group
    ...     smartsPath="Perfluoroalkyl",
    ...     constraints={"only": ["C", "F", "O", "H"]},
    ...     max_dist_from_CF=0,
    ...     linker_smarts=None
    ... )
    
    Notes
    -----
    - SMARTS patterns are compiled on initialization for efficient matching
    - Constraints are validated when checking if a molecule belongs to this group
    - max_dist_from_CF allows finding functional groups connected via non-fluorinated linkers
    - linker_smarts restricts which atoms can be in the path between component and functional group
    """
    def __init__(self, id, name,**kwargs):
        self.id = id
        self.name = name
        smarts = kwargs.get('smarts',{})
        # Save original SMARTS strings for atom counting
        if smarts and len(smarts) > 0:
            self.smarts_str, self.smarts_count = zip(*smarts.items())
        else:
            self.smarts_str = None
            self.smarts_count = None
        self.smarts = [] if self.smarts_str else None
        self.smartsPath = kwargs.get('smartsPath',None)
        self.max_dist_from_CF = kwargs.get('max_dist_from_CF', 0)
        # Compile linker_smarts pattern if provided
        linker_smarts_str = kwargs.get('linker_smarts', None)
        self.linker_smarts = None
        if linker_smarts_str is not None:
            try:
                self.linker_smarts = Chem.MolFromSmarts(linker_smarts_str)
                self.linker_smarts.UpdatePropertyCache()
                Chem.GetSymmSSSR(self.linker_smarts)
                self.linker_smarts.GetRingInfo().NumRings()
            except:
                raise ValueError(f"Invalid linker_smarts pattern '{linker_smarts_str}' for PFASGroup '{self.name}' (ID: {self.id})")
        if self.smarts_str is not None:
            for smarts_pattern in self.smarts_str:
                if smarts_pattern and smarts_pattern != "":
                    try:
                        smarts_mol = Chem.MolFromSmarts(smarts_pattern)
                        smarts_mol.UpdatePropertyCache()
                        Chem.GetSymmSSSR(smarts_mol)
                        smarts_mol.GetRingInfo().NumRings()
                        self.smarts.append(smarts_mol)
                    except:
                        raise ValueError(f"Invalid SMARTS pattern(s) for PFASGroup '{self.name}' (ID: {self.id})")
        self.constraints = kwargs.get('constraints',{})
        # Precompute number of extra atoms in SMARTS patterns (beyond matched atom and H/F/Cl/Br/I)
        self.smarts_extra_atoms = self._count_smarts_extra_atoms(self.smarts_str)
        self.component_specific_extra_atoms = []
        self.all_matches = []
        self.compute = kwargs.get('compute',True)# whether the pfasgroups needs to be parsed, or is an aggregate group, e.g. telomers
        self.re_search = kwargs.get('re_search',None)# regex for aggregate groups, e.g. telomers
        if self.re_search is not None:
            try:
                self.re_search = re.compile(self.re_search)
            except Exception as e:
                raise Exception(f"Error for agg Group {self.id}: {self.name}\n {e}")
            
    
    def _count_smarts_extra_atoms(self, smarts_str):
        """Count number of extra carbon atoms in functional group beyond what's captured by component.
        
        Parameters
        ----------
        smarts_str : str or None
            Original SMARTS string before compilation
        Returns
        -------
        int
            Number of extra carbon atoms beyond the matched atom
        
        Notes
        -----
        The component fraction calculation is now based on carbon atoms only:
        1. Carbon atoms in component
        2. Carbon atoms in SMARTS matches
        3. Additional carbon atoms from SMARTS (this return value)
        
        
        For automatic counting (when manual_size is None), this returns 0 since we now focus only
        on carbons and they are already counted in the component and SMARTS matches.
        """
        PAT_c = re.compile(r'((C(?![adeflmnorsu]))|((?<![TAS])c)|(\#6))')  # Match 'C' not followed by a letter, or c not preceded by T,A,S or #6
        return [max(0,len(PAT_c.findall(s))-1) for s in smarts_str] if smarts_str is not None else None

    
    def __str__(self):
        return self.name
    def constraint_gte(self, formula_dict):
        """Check 'greater than or equal' constraints on element counts.
        
        Parameters
        ----------
        formula_dict : dict
            Molecular formula as {element: count} dictionary
        
        Returns
        -------
        bool
            True if all 'gte' constraints are satisfied, False otherwise
        
        Examples
        --------
        >>> # Requires at least 2 carbons and 3 fluorines
        >>> group.constraints = {'gte': {'C': 2, 'F': 3}}
        >>> group.constraint_gte({'C': 3, 'F': 5, 'O': 1})  # True
        >>> group.constraint_gte({'C': 1, 'F': 5, 'O': 1})  # False (C < 2)
        """
        success = True
        for e,n in self.constraints.get('gte',{}).items():
            success = success and formula_dict.get(e,0)>=n
        return success
    def constraint_lte(self, formula_dict):
        """Check 'less than or equal' constraints on element counts.
        
        Parameters
        ----------
        formula_dict : dict
            Molecular formula as {element: count} dictionary
        
        Returns
        -------
        bool
            True if all 'lte' constraints are satisfied, False otherwise
        
        Examples
        --------
        >>> # Requires at most 2 oxygens
        >>> group.constraints = {'lte': {'O': 2}}
        >>> group.constraint_lte({'C': 8, 'F': 15, 'O': 2})  # True
        >>> group.constraint_lte({'C': 8, 'F': 15, 'O': 3})  # False (O > 2)
        """
        success = True
        for e,n in self.constraints.get('lte',{}).items():
            success = success and formula_dict.get(e,0)<=n
        return success
    def constraint_eq(self, formula_dict):
        """Check 'equal to' constraints on element counts.
        
        Parameters
        ----------
        formula_dict : dict
            Molecular formula as {element: count} dictionary
        
        Returns
        -------
        bool
            True if all 'eq' constraints are satisfied, False otherwise
        
        Examples
        --------
        >>> # Requires exactly 1 nitrogen
        >>> group.constraints = {'eq': {'N': 1}}
        >>> group.constraint_eq({'C': 8, 'F': 15, 'N': 1})  # True
        >>> group.constraint_eq({'C': 8, 'F': 15, 'N': 2})  # False (N != 1)
        """
        success = True
        for e,n in self.constraints.get('eq',{}).items():
            success = success and formula_dict.get(e,0)==n
        return success
    def constraint_only(self, formula_dict):
        """Check 'only' constraint - molecule must contain only specified elements.
        
        Parameters
        ----------
        formula_dict : dict
            Molecular formula as {element: count} dictionary
        
        Returns
        -------
        bool
            True if molecule contains only the allowed elements, False otherwise
        
        Examples
        --------
        >>> # Molecule must contain only C, F, O, H
        >>> group.constraints = {'only': ['C', 'F', 'O', 'H']}
        >>> group.constraint_only({'C': 8, 'F': 15, 'O': 2, 'H': 1})  # True
        >>> group.constraint_only({'C': 8, 'F': 15, 'O': 2, 'S': 1})  # False (S not allowed)
        
        Notes
        -----
        Checks that sum of allowed elements equals total atoms in molecule.
        """
        success = True
        if 'only' in self.constraints.keys():
            tot = sum(formula_dict.values())
            nn = 0
            for e in self.constraints['only']:
                nn += formula_dict.get(e,0)
            success = success and tot == nn
        return success
    def constraint_rel(self, formula_dict):
        """Check relational constraints between element counts.
        
        Validates relationships of the form: count(element) = f(other_elements)
        where f can include division, addition, and summing other element counts.
        
        Parameters
        ----------
        formula_dict : dict
            Molecular formula as {element: count} dictionary
        
        Returns
        -------
        bool
            True if all relational constraints are satisfied, False otherwise
        
        Notes
        -----
        Constraint Format::
        
            'rel': {
                'ElementA': {
                    'atoms': ['ElementB', 'ElementC'],  # Elements to sum
                    'div': int,  # Divisor (default 1)
                    'add': int,  # Additive constant (default 0)
                    'add_atoms': ['ElementD']  # Additional elements to add
                }
            }
        
        Formula: count(ElementA) = (sum(atoms) / div) + add + sum(add_atoms)
        
        Examples
        --------
        >>> # Carbon count must equal half the fluorine count
        >>> group.constraints = {'rel': {'C': {'atoms': ['F'], 'div': 2, 'add': 0}}}
        >>> group.constraint_rel({'C': 4, 'F': 8, 'O': 2})  # True (4 == 8/2)
        >>> group.constraint_rel({'C': 3, 'F': 8, 'O': 2})  # False (3 != 8/2)
        
        >>> # Oxygen count must equal carbon count plus 1
        >>> group.constraints = {'rel': {'O': {'atoms': ['C'], 'div': 1, 'add': 1}}}
        >>> group.constraint_rel({'C': 3, 'F': 7, 'O': 4})  # True (4 == 3 + 1)
        """
        success = True
        for e,v in self.constraints.get('rel',{}).items():
            n = sum([formula_dict.get(x,0) for x in v.get('atoms',[])])
            success = success and formula_dict.get(e,0)==n/v.get('div',1)+v.get('add',0)+sum([formula_dict.get(x,0) for x in v.get('add_atoms',[])])
        return success
    def formula_dict_satisfies_constraints(self,formula_dict):
        """Check if a molecular formula satisfies all constraints for this PFAS group.
        
        Evaluates all constraint types in order: relational → only → equal → lte → gte.
        Stops evaluation at first failure for efficiency.
        
        Parameters
        ----------
        formula_dict : dict
            Molecular formula as {element: count} dictionary (e.g., {'C': 8, 'F': 17, 'O': 2})
        
        Returns
        -------
        bool
            True if all constraints are satisfied, False if any constraint fails
        
        Constraint Evaluation Order
        ---------------------------
        1. Relational constraints ('rel') - element count relationships
        2. 'Only' constraints - allowed elements
        3. Equality constraints ('eq') - exact element counts
        4. Upper bound constraints ('lte') - maximum element counts
        5. Lower bound constraints ('gte') - minimum element counts
        
        Examples
        --------
        >>> # Perfluoroalkyl carboxylic acid constraints
        >>> group.constraints = {
        ...     'only': ['C', 'F', 'O', 'H'],  # No other elements
        ...     'gte': {'C': 2, 'F': 3},  # At least 2 carbons, 3 fluorines
        ...     'eq': {'O': 2}  # Exactly 2 oxygens
        ... }
        >>> group.formula_dict_satisfies_constraints({'C': 8, 'F': 15, 'O': 2, 'H': 1})
        True
        >>> group.formula_dict_satisfies_constraints({'C': 8, 'F': 15, 'O': 3, 'H': 1})
        False  # Fails 'eq': {'O': 2}
        
        Notes
        -----
        - Returns True immediately if no constraints are defined
        - Short-circuits on first constraint failure for performance
        - Constraint evaluation order is fixed for consistency
        """
        if len(self.constraints.keys())==0:
            return True
        success = True
        process = [None,self.constraint_rel,self.constraint_only,self.constraint_eq, self.constraint_lte, self.constraint_gte]
        k = process.pop()
        while success and k is not None:
            success = k(formula_dict)
            k = process.pop()
        return success
    def find_matched_atoms(self, mol):
        """Find all substructure matches of this PFAS group's SMARTS patterns in a molecule.
        
        Parameters
        ----------
        mol : Chem.Mol
            RDKit molecule object to search for matches
        
        Returns
        -------
        List[List[int]]
            List of matches, where each match is a list of atom indices in the molecule
        
        Notes
        -----
        - If no SMARTS patterns are defined, returns an empty list.
        - Each SMARTS pattern is searched independently; matches from all patterns are combined.
        """
        self.all_matches = []
        self.subset = set()
        if self.smarts is not None:
            for smarts_mol,min_count in zip(self.smarts, self.smarts_count):
                matches = mol.GetSubstructMatches(smarts_mol)
                if len(matches) < min_count:
                    return False
                if len(matches)>0:
                    self.all_matches.append(set(matches))
                    self.subset.update({y for x in matches for y in x if len(x)>0})
        return True
    def component_satisfies_all_smarts(self, component):
        """Check if a fluorinated component matches all SMARTS patterns of this PFAS group.
        
        Parameters
        ----------
        component : PFASComponent
            PFASComponent object representing a fluorinated component in the molecule
        
        Returns
        -------
        bool
            True if the component matches all SMARTS patterns, False otherwise
        
        Notes
        -----
        - If no SMARTS patterns are defined for this group, returns True.
        - Each SMARTS pattern must have at least one match that includes the component's atom.
        """
        atom_count = 0
        for i, (matches, min_count) in enumerate(zip(self.all_matches,self.smarts_count)):
            # Count how many SMARTS matches have overlap with this component
            # matches is a set of tuples, each tuple represents one SMARTS match
            component_set = set(component)
            found = sum(1 for match_tuple in matches if any(atom_idx in component_set for atom_idx in match_tuple))
            
            if found < min_count:
                self.component_specific_extra_atoms.append(0)
                return False
            atom_count += found * self.smarts_extra_atoms[i]
        self.component_specific_extra_atoms.append(atom_count)
        return True
    
    def find_alkyl_components(self, mol, component_solver, **kwargs):
        """Find fluorinated components in a molecule that match this PFAS group's criteria.
        
        Parameters
        ----------
        mol : Chem.Mol
            RDKit molecule object to search
        components : List[PFASComponent]
            List of PFASComponent objects representing fluorinated components in the molecule
        
        Returns
        -------
        List[PFASComponent]
            List of PFASComponent objects that match this PFAS group's criteria
        
        Notes
        -----
        - Matches are determined based on smartsPath and max_dist_from_CF attributes.
        - If smartsPath is None, all components are considered.
        - max_dist_from_CF allows extending the search radius for functional groups.
        """
        if not self.find_matched_atoms(mol):
            return 0, [], 0, []
        
        # Clear component-specific extra atoms list for this matching attempt
        self.component_specific_extra_atoms = []
        
        if self.smartsPath is None:
            # If no smartsPath specified, only check alkyl components (not cyclic)
            # This ensures functional groups like carboxylic acid (group 33) are only
            # detected when attached to perfluoroalkyl or polyfluoroalkyl chains,
            # not when attached directly to cyclic structures
            all_components = []
            for path_type in ['Perfluoroalkyl', 'Polyfluoroalkyl']:
                path_comps = component_solver.get(path_type, max_dist = self.max_dist_from_CF, default=[])
                all_components.extend(path_comps)
            components = all_components
        else:
            # Get components for the specified path type, fallback to Polyfluoroalkyl if not found
            components = component_solver.get(self.smartsPath, max_dist = self.max_dist_from_CF, default = component_solver.get("Polyfluoroalkyl", max_dist = self.max_dist_from_CF, default=[]))
        
        if self.smartsPath is not None:
            smartsPaths = [self.smartsPath]
        else:
            smartsPaths = component_solver.smartsPaths.keys()
        
        # Filter components connected to the smarts and get augmented versions
        augmented_matched_components = []
        for _smartsPath in smartsPaths:
            extended_components = component_solver.get(_smartsPath, self.max_dist_from_CF, [])
            for i, comp in enumerate(extended_components):
                # Check if this component is connected to SMARTS matches
                if self.component_satisfies_all_smarts(comp):
                    augmented = component_solver.get_augmented_component(
                        _smartsPath, self.max_dist_from_CF, i, self.subset, self.linker_smarts
                    )
                    # Accept augmented component if valid (linker validation already done in get_augmented_component)
                    if augmented is not None and len(augmented) > 0:
                        augmented_matched_components.append(
                        component_solver.get_matched_component_dict(augmented, self.subset, _smartsPath, self, comp_id = i)
                        )
        
        if len(augmented_matched_components) == 0:
            return 0, [], 0, []
        
        # Get all component sizes from all path types
        all_components = list(set([comp for comps in augmented_matched_components for comp in comps]))
        component_sizes = [len(x) for x in all_components]

        self.all_matches = []  # Clear matches after use
        self.component_specific_extra_atoms = []
        return max([0] + component_sizes), component_sizes, len(all_components), augmented_matched_components
    
    def find_aryl_components(self,mol, component_solver=None, **kwargs):
        """Find aryl components in a molecule with comprehensive metrics."""
        matches = mol.GetSubstructMatches(self.smarts[0])
        subset = [y for x in matches for y in x]
        if len(subset)==0:
            return 0, [], 0, []
        
        components = component_solver._connected_components(subset)
        component_sizes = [len(x) for x in components]
        
        # Get molecular graph for metrics calculation
        subset_set = set(subset)
        
        # Convert components to the same format as other functions return with comprehensive metrics
        matched_components = []
        for comp in components:
            matched_components.append(
                component_solver.get_matched_component_dict(comp, subset_set, 'cyclic', self)
            )
        
        return max([0]+[len(x) for x in components]), component_sizes, len(components), matched_components
        
    def find_components(self, mol, fd, component_solver, **kwargs):
        """Find fluorinated components in a molecule that match this PFAS group's criteria."""
        group_matches = []
        if self.formula_dict_satisfies_constraints(fd) is True:
            component_sizes = []
            # Create molecular graph once for this fragment
            G = mol_to_nx(mol)
            kwargs['G'] = G
            if self.smartsPath =='cyclic':
                # treat cyclic groups separately
                # Use first SMARTS pattern for cyclic
                match_count, component_sizes, matched1_len, matched_components = self.find_aryl_components(mol, component_solver=component_solver, **kwargs)
            elif self.smarts is not None and len(self.smarts) > 0:
                # Handle groups with SMARTS patterns
                match_count, component_sizes, matched1_len, matched_components = self.find_alkyl_components(mol, component_solver, **kwargs)
            elif self.smartsPath is not None:
                # treat cases with only smartsPath defined (no SMARTS patterns), find all components of that path type
                path_components = component_solver.get(self.smartsPath, max_dist = self.max_dist_from_CF, default = component_solver.get("Polyfluoroalkyl", max_dist = self.max_dist_from_CF, default=[]))
                match_count = len(path_components)
                component_sizes = [len(x) for x in path_components]
                matched1_len = len(path_components)  # Set matched1_len to enable group matching
                matched_components = []
                for comp in path_components:
                    # Use get_matched_component_dict with no SMARTS matches
                    matched_components.append(
                        component_solver.get_matched_component_dict(comp, None, self.smartsPath, self)
                    )
            else:
                # treat cases with no SMARTS patterns, just formula constraints
                match_count = 1
                matched_components = []
                matched1_len = 1
            if match_count > 0 and matched1_len > 0:
                # add to matches if functional group was found
                group_matches.append((self, match_count, component_sizes, matched_components))
            else:
                return None
            return group_matches
        return None
        
