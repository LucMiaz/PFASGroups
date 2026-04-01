# PFASGroups package
from .HalogenGroupModel import HalogenGroup
from .PFASDefinitionModel import PFASDefinition
from .ComponentsSolverModel import ComponentsSolver
from .core import rdkit_disable_log, HALOGEN_GROUPS_FILE
from .parser import parse_smiles, parse_mols, parse_mol, parse_groups_in_mol, parse_from_database, setup_halogen_groups_database, load_HalogenGroups
# PFASFingerprint: convenience alias for parse_smiles — returns a PFASEmbeddingSet
PFASFingerprint = parse_smiles
from .draw_mols import plot_mol, plot_mols, plot_HalogenGroups
from .getter import get_componentSMARTSs, get_HalogenGroups, get_compiled_HalogenGroups, get_compiled_PFASGroups, get_PFASDefinitions, get_compiled_componentSMARTSs
from .embeddings import FINGERPRINT_PRESETS, EMBEDDING_PRESETS
from .generate_homologues import generate_homologues
from .homologue_series import HomologueSeries, HomologueEntry
from .fragmentation import generate_degradation_products
# PFASEmbedding (dict subclass, primary) must be imported after embeddings to take precedence
from .PFASEmbeddings import PFASEmbedding, PFASEmbeddingSet, EmbeddingArray, ResultsModel, MoleculeResult, generate_fingerprint
from .prioritise import prioritise_molecules, prioritize_molecules, get_priority_statistics
__version__ = "3.2.0"
__all__ = ['HalogenGroup', 'PFASDefinition', 'parse_smiles', 'parse_mols','parse_mol', 'parse_groups_in_mol', 'parse_from_database', 'setup_halogen_groups_database', 'plot_HalogenGroups', 'plot_mol','plot_mols', 'FINGERPRINT_PRESETS', 'PFASFingerprint', 'generate_fingerprint', 'get_compiled_componentSMARTSs', 'get_componentSMARTSs', 'get_HalogenGroups', 'get_compiled_HalogenGroups', 'get_compiled_PFASGroups', 'get_PFASDefinitions' ,'ComponentsSolver', 'generate_homologues', 'generate_degradation_products',"rdkit_disable_log","load_HalogenGroups", "HALOGEN_GROUPS_FILE"]
__all__.extend(['PFASEmbedding', 'PFASEmbeddingSet', 'EmbeddingArray', 'ResultsModel', 'MoleculeResult', 'prioritise_molecules', 'prioritize_molecules', 'get_priority_statistics'])
__all__.extend(['HomologueSeries', 'HomologueEntry'])
