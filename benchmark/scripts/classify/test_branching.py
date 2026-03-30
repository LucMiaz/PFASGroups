
from PFASGroups.generate_mol import generate_random_mol, get_branching_index
from PFASGroups.core import mol_to_nx
from PFASGroups import parse_mol, get_compiled_PFASGroups
import numpy as np
from rdkit import Chem
from rdkit.Chem.rdMolDescriptors import CalcMolFormula
import os, sys

FILE_DIR = os.path.dirname(os.path.abspath(__file__))
SCRIPT_DIR = os.path.dirname(FILE_DIR)
REPO_DIR = os.path.dirname(os.path.dirname(SCRIPT_DIR))
DATA_DIR = os.path.join(SCRIPT_DIR, 'data')


# Try to import PFAS-Atlas
atlas_dir = os.path.join(REPO_DIR, 'PFAS-atlas')
try:
    sys.path.append(atlas_dir)  
    from classification_helper.classify_pfas import classify_pfas_molecule
    ATLAS_AVAILABLE = True
    print("✅ PFAS-Atlas available")
except ImportError:
    try:
        sys.path.append(os.path.join(atlas_dir, 'classification_helper'))
        from classify_pfas import classify_pfas_molecule
        ATLAS_AVAILABLE = True
        print("✅ PFAS-Atlas available (fallback import)")
    except ImportError:
        print("❌ PFAS-Atlas not available")
        ATLAS_AVAILABLE = False

def generate_test_cases(self,n = 300, branching_range=None, min_carbons=5, max_carbons=15, seed=2025):
    """
    Generate test cases with varying branching indices for classification testing.
    """
    if branching_range is None:
        branching_range = (0.0, 0.4)  # Default range for branching index
    np.random.seed(seed)
    PFGs = {x.id: x for x in get_compiled_PFASGroups() if x.test_dict.get('generate',{}).get('mode') == 'attach' and x.id in range(29,70)}
    test_cases = []
    duplicated = 0
    while len(test_cases) < n and duplicated < n*2:
        m = np.random.randint(min_carbons, max_carbons)
        gid = np.random.choice(list(PFGs.keys()),1)[0]
        mol = generate_random_mol(m,PFGs[gid].test_dict.get('generate').get('smiles'), mode = PFGs[gid].test_dict.get('generate').get('mode'),perfluorinated=True, branching_range=branching_range)
        smiles = Chem.MolToSmiles(mol)
        formula = CalcMolFormula(mol)
        branching_index = get_branching_index(mol)
        # check if duplicated        
        if any(tc['smiles'] == smiles for tc in test_cases):
            duplicated += 1
            continue
        test_cases.append({"group_id": gid, "PFASGroup": PFGs[gid], "smiles": smiles, "formula": formula, "branching_index": branching_index})
    return test_cases

class TestBranching:
    def __init__(self, n=300, **kwargs):
        self.n = n
        self.test_cases = generate_test_cases(n=n, **kwargs)

    def test_branched_compounds(self):
        for tc in self.test_cases:
            mol = Chem.MolFromSmiles(tc['smiles'])
            result = parse_mol(mol, include_PFAS_definitions=False)
            groups = [m for m in (result.get("matches", []) if isinstance(result, dict) else []) if m.get("type") == "HalogenGroup"]
        group_ids = [g.get("id") for g in groups]
        assert tc['group_id'] in group_ids, f"Failed to classify {tc['smiles']} with branching index {tc['branching_index']:.2f} as group {tc['group_id']}"
    
    def test_with_PFASGroups(self, smiles, include_PFAS_definitions=True,
                             limit_effective_graph_resistance=None,
                             compute_component_metrics=True):
        """Test molecule with PFASGroups detection
        
        Args:
            smiles: SMILES string
            include_PFAS_definitions: Whether to include PFAS definitions (True for accuracy, False for specificity)
            limit_effective_graph_resistance: Limit or disable graph resistance computation (None=all, False/0=skip)
            compute_component_metrics: Whether to compute component graph metrics
        """
        
        start_time = time.perf_counter()
        PFASGroups_result = {
            'detected_groups': [],
            'detected_definitions': [],
            'matches': [],
            'success': False,
            'error': None,
            'execution_time': 0.0,
            'include_definitions': include_PFAS_definitions
        }
        
        if not PFASGroups_AVAILABLE:
            PFASGroups_result['error'] = 'PFASGroups not available'
            return PFASGroups_result
        
        try:
            mol = Chem.MolFromSmiles(smiles)
            if mol is not None:
                # Use parse_mol which returns dict with new format
                results = parse_mol(
                    mol,
                    include_PFAS_definitions=include_PFAS_definitions,
                    limit_effective_graph_resistance=limit_effective_graph_resistance,
                    compute_component_metrics=compute_component_metrics
                )
                
                # Extract groups and definitions from the new dictionary format
                # results is a dict with 'matches' key containing list of match dicts
                group_ids = []
                definition_ids = []
                all_matches = []
                
                if isinstance(results, dict) and 'matches' in results:
                    for match in results['matches']:
                        if match.get('type') == 'HalogenGroup':
                            group_ids.append(match['id'])
                            all_matches.append({
                                'type': 'group',
                                'id': match['id'],
                                'match_id': match.get('match_id'),
                                'name': match.get('group_name'),
                                'match_count': match.get('match_count'),
                                'components_sizes': match.get('components_sizes', []),
                                'num_components': match.get('num_components', 0),
                                'components_types': match.get('components_types', []),
                                # Summary metrics
                                'mean_eccentricity': match.get('mean_eccentricity', 0.0),
                                'mean_diameter': match.get('mean_diameter', float('nan')),
                                'mean_radius': match.get('mean_radius', float('nan'))
                            })
                        elif match.get('type') == 'PFASdefinition':
                            definition_ids.append(match['id'])
                            all_matches.append({
                                'type': 'definition',
                                'id': match['id'],
                                'match_id': match.get('match_id'),
                                'name': match.get('definition_name')
                            })
                
                PFASGroups_result['detected_groups'] = group_ids
                PFASGroups_result['detected_definitions'] = definition_ids
                PFASGroups_result['matches'] = all_matches
                PFASGroups_result['success'] = len(group_ids) > 0 or len(definition_ids) > 0
                
        except Exception as e:
            PFASGroups_result['error'] = str(e)
        finally:
            PFASGroups_result['execution_time'] = time.perf_counter() - start_time
        
        return PFASGroups_result
    
    def test_with_atlas(self, smiles):
        """Test molecule with PFAS-Atlas classification"""
        
        start_time = time.perf_counter()
        atlas_result = {
            'first_class': None,
            'second_class': None, 
            'success': False,
            'error': None,
            'execution_time': 0.0
        }
        
        if not ATLAS_AVAILABLE:
            atlas_result['error'] = 'PFAS-Atlas not available'
            return atlas_result
        
        try:
            # Use the correct PFAS-Atlas function classify_pfas_molecule
            predictions = classify_pfas_molecule(smiles)
            
            if predictions and len(predictions) >= 2:
                atlas_result['first_class'] = predictions[0]
                atlas_result['second_class'] = predictions[1]
                atlas_result['success'] = True
            elif predictions and len(predictions) >= 1:
                atlas_result['first_class'] = predictions[0]
                atlas_result['success'] = True
                    
        except Exception as e:
            atlas_result['error'] = str(e)
        finally:
            atlas_result['execution_time'] = time.perf_counter() - start_time
        
        return atlas_result
    def run_tests(self):
        for tc in self.test_cases:
            print(f"Testing SMILES: {tc['smiles']} (Branching Index: {tc['branching_index']:.2f})")
            pf_result = self.test_with_PFASGroups(tc['smiles'])
            tc['PFASGroups_result'] = pf_result
            if ATLAS_AVAILABLE:
                atlas_result = self.test_with_atlas(tc['smiles'])
                tc['PFASAtlas_result'] = atlas_result
            print("-" * 80)
    def summarize_results(self):
        self.summary = {
            'total_cases': len(self.test_cases),
            'PFASGroups': {
                'success_count': sum(1 for tc in self.test_cases if tc.get('PFASGroups_result', {}).get('success')),
                'failure_count': sum(1 for tc in self.test_cases if not tc.get('PFASGroups_result', {}).get('success')),
                'errors': [tc['PFASGroups_result']['error'] for tc in self.test_cases if tc.get('PFASGroups_result', {}).get('error')]
            },
            'PFASAtlas': {
                'success_count': sum(1 for tc in self.test_cases if tc.get('PFASAtlas_result', {}).get('success')),
                'failure_count': sum(1 for tc in self.test_cases if not tc.get('PFASAtlas_result', {}).get('success')),
                'errors': [tc['PFASAtlas_result']['error'] for tc in self.test_cases if tc.get('PFASAtlas_result', {}).get('error')]
            } if ATLAS_AVAILABLE else None
        }
    def save(self, filename):
        """Save test cases and results to a JSON file"""
        import json
        self.summarize_results()
        results = [self.summary] + self.test_cases
        with open(filename, 'w') as f:
            json.dump(self.test_cases, f, indent=2)
        
if __name__ == "__main__":
    seed = 2025
    tester = TestBranching(n=300, branching_range=(0.0, 0.4), min_carbons=5, max_carbons=15, seed=seed)
    tester.run_tests()
    tester.save(os.path.join(DATA_DIR, f"branching_test_results_{seed}.json"))
