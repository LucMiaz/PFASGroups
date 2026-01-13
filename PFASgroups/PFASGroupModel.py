from rdkit import Chem
import numpy as np

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
    smarts1 : Chem.Mol or None
        Primary SMARTS pattern (compiled RDKit molecule) for functional group detection.
        None if group is defined by smartsPath alone.
    smarts2 : Chem.Mol or None
        Secondary SMARTS pattern for groups requiring two functional groups.
        None if only one functional group is needed.
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
    ...     smarts1="C(=O)O",  # Carboxylic acid group
    ...     smarts2=None,
    ...     smartsPath="Perfluoroalkyl",
    ...     constraints={"only": ["C", "F", "O", "H"]},
    ...     max_dist_from_CF=0
    ... )
    
    Notes
    -----
    - SMARTS patterns are compiled on initialization for efficient matching
    - Constraints are validated when checking if a molecule belongs to this group
    - max_dist_from_CF allows finding functional groups connected via non-fluorinated linkers
    """
    def __init__(self, id, name, smarts1, smarts2, smartsPath, constraints,**kwargs):
        self.id = id
        self.name = name
        # Save original SMARTS strings for atom counting
        self.smarts1_str = smarts1 if smarts1 and smarts1 != "" else None
        self.smarts2_str = smarts2 if smarts2 and smarts2 != "" else None
        # Get manual SMARTS sizes if provided in JSON
        self.smarts1_size = kwargs.get('smarts1_size', None)
        self.smarts2_size = kwargs.get('smarts2_size', None)
        self.smarts1 = smarts1
        self.smarts2 = smarts2
        self.smartsPath = smartsPath
        self.max_dist_from_CF = kwargs.get('max_dist_from_CF', 0)
        if self.smarts1 !="" and self.smarts1 is not None:
            try:
                self.smarts1 = Chem.MolFromSmarts(self.smarts1)
                self.smarts1.UpdatePropertyCache()
                Chem.GetSymmSSSR(self.smarts1)
                self.smarts1.GetRingInfo().NumRings()
                if self.smarts2 !="" and self.smarts2 is not None:
                    self.smarts2 = Chem.MolFromSmarts(self.smarts2)
                    self.smarts2.UpdatePropertyCache()
                    Chem.GetSymmSSSR(self.smarts2)
                    self.smarts2.GetRingInfo().NumRings()
                else:
                    self.smarts2 = None
            except:
                raise ValueError(f"Invalid SMARTS pattern(s) for PFASGroup '{self.name}' (ID: {self.id})")
        else:
            self.smarts1 = None
            self.smarts2 = None
        self.constraints = constraints
        # Precompute number of extra atoms in SMARTS patterns (beyond matched atom and H/F/Cl/Br/I)
        self.smarts1_extra_atoms = self._count_smarts_extra_atoms(self.smarts1, self.smarts1_str, self.smarts1_size)
        self.smarts2_extra_atoms = self._count_smarts_extra_atoms(self.smarts2, self.smarts2_str, self.smarts2_size)
    
    def _count_smarts_extra_atoms(self, smarts_mol, smarts_str, manual_size=None):
        """Count number of extra carbon atoms in functional group beyond what's captured by component.
        
        Parameters
        ----------
        smarts_mol : Chem.Mol or None
            Compiled SMARTS pattern
        smarts_str : str or None
            Original SMARTS string before compilation
        manual_size : int or None
            Manually specified number of carbon atoms in the SMARTS pattern (from JSON).
            If provided, returns manual_size - 1 (to exclude the matched carbon atom).
            
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
        
        If manual_size is provided in the JSON (smarts1_size or smarts2_size), it represents
        the total number of CARBON atoms in the functional group, so we subtract 1 for the matched carbon.
        
        For automatic counting (when manual_size is None), this returns 0 since we now focus only
        on carbons and they are already counted in the component and SMARTS matches.
        """
        # Use manual size if provided (subtract 1 for the matched carbon atom itself)
        if manual_size is not None:
            return max(0, manual_size - 1)
        
        # Since we now focus on carbon atoms only, and carbons are already counted
        # in the component and SMARTS matches, we don't need to add extra atoms
        # when manual_size is not provided
        return 0
    
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
        
        Constraint Format
        -----------------
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