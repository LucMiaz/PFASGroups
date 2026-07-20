from rdkit import Chem
from typing import List, Dict, Optional, Union

class PFASDefinition:
    """Model class representing a PFAS definition based on structural criteria.
    
    A PFAS definition identifies molecules using SMARTS patterns and/or fluorine ratio thresholds.
    Unlike HalogenGroup which focuses on specific functional groups, PFASDefinition provides
    broader chemical definitions (e.g., "contains at least one perfluoroalkyl moiety").
    
    Attributes
    ----------
    id : int
        Unique identifier for this PFAS definition
    name : str
        Human-readable name (e.g., "Per- and polyfluoroalkyl substances")
    description : str
        Detailed description of what this definition represents
    fluorineRatio : Optional[float]
        Minimum ratio of fluorine atoms required (None if not applicable)
    smarts_strings : List[str]
        Original SMARTS pattern strings for structural matching
    smarts_patterns : List[Chem.Mol]
        Compiled SMARTS molecule objects for efficient matching
    includeHydrogen : bool
        Whether to include hydrogen atoms in fluorine ratio calculations
    requireBoth : bool
        If True, requires both SMARTS match AND fluorine ratio.
        If False, requires SMARTS match OR fluorine ratio.
    
    Examples
    --------
    >>> # Definition requiring perfluoroalkyl chain OR high fluorine ratio
    >>> pfas_def = PFASDefinition(
    ...     id=1,
    ...     name="PFAS (OECD definition)",
    ...     smarts=["[CX4][CX4]([F])([F])[F]"],
    ...     fluorineRatio=0.4,
    ...     description="Contains perfluoroalkyl moiety with ≥2 carbons",
    ...     requireBoth=False
    ... )
    """
    def __init__(self, id: int, name: str, smarts: List[str], fluorineRatio: Optional[float], description: str, **kwargs):
        self.id = id
        self.name = name
        self.description = description
        self.fluorineRatio = fluorineRatio
        self.smarts_strings = smarts
        self.includeHydrogen = kwargs.get('includeHydrogen', True)
        self.requireBoth = kwargs.get('requireBoth', False)
        # Preload SMARTS patterns
        self.smarts_patterns = []
        for smarts_str in smarts:
            try:
                mol = Chem.MolFromSmarts(smarts_str)
                if mol is not None:
                    mol.UpdatePropertyCache()
                    Chem.GetSymmSSSR(mol)
                    mol.GetRingInfo().NumRings()
                    self.smarts_patterns.append(mol)
            except:
                pass
    
    def __str__(self):
        return self.name
    
    def applies_to_molecule(self,
                  mol_or_smiles: Union[Chem.Mol, str],
                  formula: Optional[Dict[str, int]] = None,
                  **kwargs) -> bool:
        """Check if this PFAS definition applies to a given molecule.

        This method evaluates whether a molecule meets the structural and/or compositional
        criteria defined by this PFASDefinition. The evaluation logic depends on the
        requireBoth flag:

        - If requireBoth=False (default): Returns True if EITHER SMARTS matches OR
          fluorine ratio is met (logical OR)
        - If requireBoth=True: Returns True only if BOTH SMARTS matches AND fluorine
          ratio are met (logical AND)

        Parameters
        ----------
        mol_or_smiles : Union[Chem.Mol, str]
            Input molecule as RDKit Mol object or SMILES string
        formula : Optional[Dict[str, int]], default=None
            Pre-computed molecular formula as {element: count} dictionary.
            If None, will be computed from the molecule.
        **kwargs : dict
            Additional parameters:

            - include_hydrogen (bool): Whether to include H in fluorine ratio calculation.
              Defaults to self.includeHydrogen
            - require_both (bool): Override the instance's requireBoth setting

        Returns
        -------
        bool
            True if the molecule meets the definition criteria, False otherwise

        Examples
        --------
        >>> pfas_def = PFASDefinition(
        ...     id=1, name="Test", smarts=["[CX4]F"],
        ...     fluorineRatio=0.3, description="Test"
        ... )
        >>> pfas_def.applies_to_molecule("FC(F)(F)C(F)(F)F")  # PFOA-like
        True
        >>> pfas_def.applies_to_molecule("CCCCCC")  # No fluorine
        False

        Notes
        -----
        - SMARTS patterns are checked using substructure matching (HasSubstructMatch)
        - Fluorine ratio is mass-weighted: F mass / total mass (see `match_details`
          for the breakdown of which criterion -- SMARTS, ratio, or both -- was met,
          and the actual computed ratio value)
        - Invalid SMILES strings return False
        """
        return self.match_details(mol_or_smiles, formula=formula, **kwargs)['applies']

    def match_details(self,
                  mol_or_smiles: Union[Chem.Mol, str],
                  formula: Optional[Dict[str, int]] = None,
                  **kwargs) -> Dict:
        """Evaluate a molecule against this definition and report each criterion separately.

        Unlike `applies_to_molecule`, which collapses everything to a single bool,
        this returns the SMARTS-match result, the fluorine-ratio result, and the
        actual computed ratio value individually -- useful for explaining *why* a
        molecule matched (e.g. distinguishing a genuine structural PFAS from one
        that only clears a fluorine-content threshold).

        Parameters
        ----------
        mol_or_smiles : Union[Chem.Mol, str]
            Input molecule as RDKit Mol object or SMILES string
        formula : Optional[Dict[str, int]], default=None
            Pre-computed molecular formula as {element: count} dictionary.
            If None, will be computed from the molecule.
        **kwargs : dict
            Same as `applies_to_molecule` (`include_hydrogen`, `require_both`).

        Returns
        -------
        dict
            - 'valid' (bool): False if the input SMILES could not be parsed
            - 'smarts_match' (bool): whether any SMARTS pattern matched
            - 'fluorine_ratio' (float or None): the computed ratio, or None if this
              definition has no fluorineRatio threshold
            - 'fluorine_ratio_threshold' (float or None): this definition's threshold
            - 'ratio_match' (bool): whether the ratio met the threshold (True if no
              threshold is defined, matching `applies_to_molecule`'s OR semantics)
            - 'applies' (bool): the same overall result `applies_to_molecule` returns
        """
        # Convert SMILES to Mol if needed
        if isinstance(mol_or_smiles, str):
            mol = Chem.MolFromSmiles(mol_or_smiles)
            if mol is None:
                return {
                    'valid': False,
                    'smarts_match': False,
                    'fluorine_ratio': None,
                    'fluorine_ratio_threshold': self.fluorineRatio,
                    'ratio_match': False,
                    'applies': False,
                }
        else:
            mol = mol_or_smiles

        # Check SMARTS matches
        smarts_match = False
        for pattern in self.smarts_patterns:
            if mol.HasSubstructMatch(pattern):
                smarts_match = True
                break

        # Check fluorine ratio if defined
        ratio_match = True  # Default to True if no ratio requirement
        fluorine_ratio_value = None
        if self.fluorineRatio is not None:
            if formula is None:
                formula = self._compute_formula(mol, kwargs.get("include_hydrogen", self.includeHydrogen))

            fluorine_ratio_value = self._compute_fluorine_ratio(formula, kwargs.get("include_hydrogen", self.includeHydrogen))
            ratio_match = fluorine_ratio_value >= self.fluorineRatio

        # Apply logic based on require_both
        if kwargs.get("require_both", self.requireBoth):
            applies = smarts_match and ratio_match
        elif self.fluorineRatio is None:
            # If no fluorine ratio is defined, only check SMARTS
            applies = smarts_match
        else:
            # Otherwise, SMARTS OR ratio
            applies = smarts_match or ratio_match

        return {
            'valid': True,
            'smarts_match': smarts_match,
            'fluorine_ratio': round(fluorine_ratio_value, 4) if fluorine_ratio_value is not None else None,
            'fluorine_ratio_threshold': self.fluorineRatio,
            'ratio_match': ratio_match,
            'applies': applies,
        }
    
    def _compute_formula(self, mol: Chem.Mol, include_hydrogen: bool) -> Dict[str, int]:
        """Compute molecular formula as element count dictionary.
        
        Parameters
        ----------
        mol : Chem.Mol
            RDKit molecule object
        include_hydrogen : bool
            If True, add explicit hydrogens before counting
        
        Returns
        -------
        Dict[str, int]
            Dictionary mapping element symbols to their counts, e.g. {'C': 8, 'F': 17, 'O': 2}
        """
        formula = {}
        
        # Add hydrogens if needed
        if include_hydrogen:
            mol = Chem.AddHs(mol)
        
        for atom in mol.GetAtoms():
            symbol = atom.GetSymbol()
            formula[symbol] = formula.get(symbol, 0) + 1
        
        return formula
    
    def _check_fluorine_ratio(self, formula: Dict[str, int], include_hydrogen: bool) -> bool:
        """Check if the fluorine ratio in a molecular formula meets the threshold.

        Mass-weighted (F mass / total mass):
        `fluorine_ratio = F_mass / (molecular_weight - H_mass)`.

        Parameters
        ----------
        formula : Dict[str, int]
            Molecular formula as {element: count} dictionary
        include_hydrogen : bool
            If True, includes hydrogen mass in the total mass.
            If False, excludes hydrogen mass from total (heavy atoms only) --
            this is the "heavy atoms only" mode PFASSTRUCTv5 uses.

        Returns
        -------
        bool
            True if (F_mass / total_mass) >= self.fluorineRatio, False otherwise

        Notes
        -----
        - Returns False if total_mass is 0
        - Formula with no fluorine (F_count=0) will fail unless fluorineRatio=0
        """
        from rdkit.Chem import GetPeriodicTable
        pt = GetPeriodicTable()

        f_count = formula.get('F', 0)
        f_mass = f_count * pt.GetAtomicWeight('F')

        if include_hydrogen:
            total_mass = sum(count * pt.GetAtomicWeight(sym) for sym, count in formula.items())
        else:
            total_mass = sum(count * pt.GetAtomicWeight(sym) for sym, count in formula.items() if sym != 'H')

        if total_mass == 0:
            return False

        ratio = f_mass / total_mass
        return ratio >= self.fluorineRatio
    
    def test(self, test_data=None):
        """Test this PFAS definition against test molecules from metadata.
        
        Validates that the definition correctly classifies true positives, true negatives,
        false positives, and false negatives based on test metadata in 
        PFAS_definitions_smarts.json.
        
        Parameters
        ----------
        test_data : dict, optional
            Test metadata dictionary. If None, will be loaded from the definition's
            entry in PFAS_definitions_smarts.json. Expected keys: ``category``,
            ``examples`` (dict with keys ``true_positives``, ``true_negatives``,
            ``false_positives``, ``false_negatives`` each a list of dicts).
        
        Returns
        -------
        dict
            Test results with keys: ``passed`` (bool), ``total_tests`` (int),
            ``failures`` (list), ``category`` (str), ``stats`` (dict with
            counts for true/false positives/negatives).
        
        Notes
        -----
        - Tests against benchmark test compounds with known PFAS/non-PFAS labels
        - Validates both SMARTS patterns and fluorine ratio criteria
        - Returns detailed failure information for debugging
        """
        from rdkit import Chem
        
        # Load test data if not provided
        if test_data is None:
            import json
            from pathlib import Path
            definitions_file = Path(__file__).parent / 'data' / 'PFAS_definitions_smarts.json'
            with open(definitions_file, 'r') as f:
                all_definitions = json.load(f)
            
            # Find this definition's test data
            test_data = None
            for def_data in all_definitions:
                if def_data['id'] == self.id:
                    test_data = def_data.get('test', {})
                    break
            
            if test_data is None or not test_data:
                return {
                    'passed': None,
                    'total_tests': 0,
                    'failures': [],
                    'category': 'unknown',
                    'error': f'No test data found for definition {self.id}'
                }
        
        results = {
            'passed': True,
            'total_tests': 0,
            'failures': [],
            'category': test_data.get('category', 'definition'),
            'stats': {
                'true_positives': 0,
                'true_negatives': 0,
                'false_positives': 0,
                'false_negatives': 0
            }
        }
        
        examples = test_data.get('examples', {})
        
        # Test positives (should match)
        for smiles in examples.get('positives', []):
            results['total_tests'] += 1
            try:
                mol = Chem.MolFromSmiles(smiles)
                if mol is None:
                    results['passed'] = False
                    results['failures'].append({
                        'smiles': smiles,
                        'expected': True,
                        'got': None,
                        'type': 'true_positive',
                        'error': 'Invalid SMILES'
                    })
                    continue
                
                # Check if definition applies
                applies = self.applies_to_molecule(mol)
                if applies:
                    results['stats']['true_positives'] += 1
                else:
                    results['passed'] = False
                    results['stats']['false_negatives'] += 1
                    results['failures'].append({
                        'smiles': smiles,
                        'expected': True,
                        'got': False,
                        'type': 'true_positive',
                        'error': 'Definition should match but did not'
                    })
            except Exception as e:
                results['passed'] = False
                results['failures'].append({
                    'smiles': smiles,
                    'expected': True,
                    'got': None,
                    'type': 'positive',
                    'error': f'Exception: {str(e)}'
                })
        
        # Test negatives (should NOT match)
        for smiles in examples.get('negatives', []):
            results['total_tests'] += 1
            try:
                mol = Chem.MolFromSmiles(smiles)
                if mol is None:
                    results['passed'] = False
                    results['failures'].append({
                        'smiles': smiles,
                        'expected': False,
                        'got': None,
                        'type': 'true_negative',
                        'error': 'Invalid SMILES'
                    })
                    continue
                
                # Check if definition applies
                applies = self.applies_to_molecule(mol)
                if not applies:
                    results['stats']['true_negatives'] += 1
                else:
                    results['passed'] = False
                    results['stats']['false_positives'] += 1
                    results['failures'].append({
                        'smiles': smiles,
                        'expected': False,
                        'got': True,
                        'type': 'true_negative',
                        'error': 'Definition should not match but did'
                    })
            except Exception as e:
                results['passed'] = False
                results['failures'].append({
                    'smiles': smiles,
                    'expected': False,
                    'got': None,
                    'type': 'negative',
                    'error': f'Exception: {str(e)}'
                })
        
        # Test false positives (known to incorrectly match - document these)
        for item in examples.get('false_positives', []):
            smiles = item if isinstance(item, str) else item.get('smiles', '')
            results['total_tests'] += 1
            try:
                mol = Chem.MolFromSmiles(smiles)
                if mol is None:
                    continue
                
                # These are expected to match (false positives)
                applies = self.applies_to_molecule(mol)
                if applies:
                    results['stats']['false_positives'] += 1
            except Exception:
                pass
        
        # Test false negatives (known to incorrectly NOT match - document these)
        for item in examples.get('false_negatives', []):
            smiles = item if isinstance(item, str) else item.get('smiles', '')
            results['total_tests'] += 1
            try:
                mol = Chem.MolFromSmiles(smiles)
                if mol is None:
                    continue
                
                # These are expected to NOT match (false negatives)
                applies = self.applies_to_molecule(mol)
                if not applies:
                    results['stats']['false_negatives'] += 1
            except Exception:
                pass
        
        return results
