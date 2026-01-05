# PFASgroups package
from .PFASGroupModel import PFASGroup
from .PFASDefinitionModel import PFASDefinition
from .core import parse_smiles, parse_mols, parse_mol,parse_groups_in_mol, plot_pfasgroups, generate_fingerprint, get_smartsPaths, get_PFASDefinitions, compile_smartsPath, compile_smartsPaths
from .draw_mols import plot_mol, plot_mols
from .generate_homologues import generate_homologues
from .fragmentation import generate_degradation_products
# Optional imports for testing
try:
    from . import test_examples
    from . import generate_mol
except ImportError:
    # These modules might not be available in all environments
    pass

__version__ = "0.1.0"
__all__ = ['PFASGroup', 'PFASDefinition', 'parse_smiles', 'parse_mols','parse_mol', 'parse_groups_in_mol', 'plot_pfasgroups', 'plot_mol','plot_mols',"generate_fingerprint", 'get_smartsPaths', 'get_PFASGroups', 'get_PFASDefinitions' ,'compile_smartsPath', 'compile_smartsPaths']
