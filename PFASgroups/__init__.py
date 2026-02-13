# PFASgroups package
from .PFASGroupModel import PFASGroup
from .PFASDefinitionModel import PFASDefinition
from .ComponentsSolverModel import ComponentsSolver
from .core import rdkit_disable_log, PFAS_GROUPS_FILE
from .parser import parse_smiles, parse_mols, parse_mol, parse_groups_in_mol, parse_from_database, setup_pfas_groups_database, compile_componentSmarts, compile_componentSmartss, load_PFASGroups
from .draw_mols import plot_mol, plot_mols, plot_pfasgroups
from .getter import get_componentSmartss, get_PFASGroups, get_PFASDefinitions
from .fingerprints import generate_fingerprint
from .generate_homologues import generate_homologues
from .fragmentation import generate_degradation_products
from .results_model import ResultsModel
# Optional imports for testing
try:
    from . import test_examples
    from . import generate_mol
except ImportError:
    # These modules might not be available in all environments
    pass

__version__ = "2.2.1"
__all__ = ['PFASGroup', 'PFASDefinition', 'parse_smiles', 'parse_mols','parse_mol', 'parse_groups_in_mol', 'parse_from_database', 'setup_pfas_groups_database', 'plot_pfasgroups', 'plot_mol','plot_mols',"generate_fingerprint", 'get_componentSmartss', 'get_PFASGroups', 'get_PFASDefinitions' ,'compile_componentSmarts', 'compile_componentSmartss','ComponentsSolver', 'generate_homologues', 'generate_degradation_products',"rdkit_disable_log","load_PFASGroups", "PFAS_GROUPS_FILE"]
__all__.append('ResultsModel')
