from rdkit import Chem
from typing import List, Dict, Optional, Union
class PFASDefinition:
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
        """
        Check if this PFAS definition applies to a molecule.
        
        Args:
            mol_or_smiles: RDKit molecule or SMILES string
            formula: Optional pre-computed formula dict
            include_hydrogen: Whether to include H atoms in fluorine ratio calculation
            require_both: If True, requires both SMARTS match AND fluorine ratio.
                         If False (default), requires SMARTS match OR fluorine ratio.
        
        Returns:
            True if the definition applies, False otherwise
        """
        # Convert SMILES to Mol if needed
        if isinstance(mol_or_smiles, str):
            mol = Chem.MolFromSmiles(mol_or_smiles)
            if mol is None:
                return False
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
        if self.fluorineRatio is not None:
            if formula is None:
                formula = self._compute_formula(mol, kwargs.get("include_hydrogen", self.includeHydrogen))
            
            ratio_match = self._check_fluorine_ratio(formula, kwargs.get("include_hydrogen", self.includeHydrogen))
        
        # Apply logic based on require_both
        if kwargs.get("require_both", self.requireBoth):
            return smarts_match and ratio_match
        else:
            # If no fluorine ratio is defined, only check SMARTS
            if self.fluorineRatio is None:
                return smarts_match
            # Otherwise, SMARTS OR ratio
            return smarts_match or ratio_match
    
    def _compute_formula(self, mol: Chem.Mol, include_hydrogen: bool) -> Dict[str, int]:
        """Compute molecular formula as a dictionary."""
        formula = {}
        
        # Add hydrogens if needed
        if include_hydrogen:
            mol = Chem.AddHs(mol)
        
        for atom in mol.GetAtoms():
            symbol = atom.GetSymbol()
            formula[symbol] = formula.get(symbol, 0) + 1
        
        return formula
    
    def _check_fluorine_ratio(self, formula: Dict[str, int], include_hydrogen: bool) -> bool:
        """Check if fluorine ratio meets the threshold."""
        f_count = formula.get('F', 0)
        
        if include_hydrogen:
            total_atoms = sum(formula.values())
        else:
            total_atoms = sum(v for k, v in formula.items() if k != 'H')
        
        if total_atoms == 0:
            return False
        
        ratio = f_count / total_atoms
        return ratio >= self.fluorineRatio

