
from PFASGroups.generate_mol import generate_random_mol, get_branching_index
from PFASGroups.core import mol_to_nx
from PFASGroups import parse_mol, get_compiled_PFASGroups
import numpy as np
from rdkit import Chem
from rdkit.Chem.rdMolDescriptors import CalcMolFormula
from rdkit import rdBase
import os, sys
import time
rdBase.DisableLog('rdApp.warning')  # Suppress RDKit warnings


CLASSIFY_DIR = os.path.dirname(os.path.abspath(__file__))
SCRIPT_DIR = os.path.dirname(CLASSIFY_DIR)
BENCH_DIR = os.path.dirname(SCRIPT_DIR)
REPO_DIR = os.path.dirname(BENCH_DIR)
GIT_DIR = os.path.dirname(REPO_DIR)
DATA_DIR = os.path.join(BENCH_DIR, 'data')

# Try to import PFAS-Atlas
atlas_dir = os.path.join(GIT_DIR, 'PFAS-atlas')
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
        print(f"❌ PFAS-Atlas not available as {atlas_dir} or {os.path.join(atlas_dir, 'classification_helper')}")
        ATLAS_AVAILABLE = False

def generate_highly_branched_failing_cases(n, min_carbons=8, max_carbons=16, seed=None,
                                           existing_smiles=None, p_fluorination=0.3,
                                           phigh_fluorination=0.65, branching_range=(0.0, 0.2)):
    """Generate very highly branched molecules that fail at least one PFAS definition.

    Uses partial fluorination (not perfluorinated) and tight branching constraints so
    that some molecules lack the CF₂–CF₂ or CF₃ stretches required by stricter
    definitions (e.g. UK Definition 4).

    Parameters
    ----------
    n : int
        Number of molecules to attempt to generate.
    min_carbons, max_carbons : int
        Carbon-chain size range (larger chains give more branching possibility).
    seed : int or None
        Seed appended to the global RNG state (does not reset it).
    existing_smiles : set or None
        SMILES already in the dataset; duplicates are skipped.
    p_fluorination : float
        Minimum H→F probability passed to ``fluorinate_mol``.
    phigh_fluorination : float
        Maximum H→F probability passed to ``fluorinate_mol``.
    branching_range : tuple of float
        (min_BI, max_BI) – very low values (close to 0) give highly branched molecules.

    Returns
    -------
    list of dict
        Each entry has the same keys as regular test cases plus
        ``'highly_branched_failing': True`` and
        ``'detected_definitions_at_generation': [int, ...]``.
    """
    if existing_smiles is None:
        existing_smiles = set()
    if seed is not None:
        np.random.seed(seed)

    PFGs = {x.id: x for x in get_compiled_PFASGroups()
            if x.test_dict.get('generate', {}).get('mode') == 'attach'
            and x.id in [42, 61]}# keep the groups simple to be sure PFAS-Atlas should find them

    cases = []
    max_attempts = n * 40

    for _ in range(max_attempts):
        if len(cases) >= n:
            break
        m = np.random.randint(min_carbons, max_carbons)
        gid = np.random.choice(list(PFGs.keys()), 1)[0]
        fg = PFGs[gid].test_dict.get('generate', {})

        try:
            mol = generate_random_mol(
                m,
                fg.get('smiles'),
                mode=fg.get('mode', 'attach'),
                perfluorinated=False,
                p=p_fluorination,
                phigh=phigh_fluorination,
                branching_range=branching_range,
                max_matched_definitions=4,
            )
        except Exception:
            continue

        smiles = Chem.MolToSmiles(mol)
        if smiles in existing_smiles:
            continue

        # Record which definitions are matched (generate_random_mol guarantees ≤4)
        try:
            result = parse_mol(mol, include_PFAS_definitions=True)
            detected = sorted({match['id'] for match in result.get('matches', [])
                               if match.get('type') == 'PFASdefinition'})
        except Exception:
            detected = []

        branching_index = get_branching_index(mol)
        formula = CalcMolFormula(mol)
        existing_smiles.add(smiles)
        cases.append({
            'group_id': gid,
            'PFASGroup': PFGs[gid],
            'smiles': smiles,
            'formula': formula,
            'branching_index': branching_index,
            'detected_definitions_at_generation': sorted(detected),
            'highly_branched_failing': True,
        })

    if len(cases) < n:
        print(f"⚠  Could only generate {len(cases)}/{n} highly branched failing cases "
              f"(increased max_attempts or lower p_fluorination may help)")
    return cases


def generate_test_cases(n = 300, branching_range=None, min_carbons=5, max_carbons=15, seed=2025,
                        n_highly_branched_failing_proportion=0.25):
    """
    Generate test cases with varying branching indices for classification testing.
    """
    if branching_range is None:
        branching_range = (0.0, 0.4)  # Default range for branching index
    np.random.seed(seed)
    n_highly_branched_failing = int(n * n_highly_branched_failing_proportion)
    n = n - n_highly_branched_failing  # Adjust n to account for highly branched failing cases
    PFGs = {x.id: x for x in get_compiled_PFASGroups() if x.test_dict.get('generate',{}).get('mode') == 'attach' and x.id in [42, 61]}# keep the groups simple to be sure PFAS-Atlas should find them
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

    if n_highly_branched_failing > 0:
        existing = {tc['smiles'] for tc in test_cases}
        failing = generate_highly_branched_failing_cases(
            n=n_highly_branched_failing,
            min_carbons=max(min_carbons, 8),
            max_carbons=max(max_carbons, 16),
            existing_smiles=existing,
        )
        test_cases.extend(failing)
        print(f"Added {len(failing)} highly branched failing cases "
              f"(target: {n_highly_branched_failing})")

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
    
    def test_with_PFASGroups(self, test_case, include_PFAS_definitions=True,
                             limit_effective_graph_resistance=None,
                             compute_component_metrics=True):
        """Test molecule with PFASGroups detection
        
        Args:
            smiles: SMILES string
            include_PFAS_definitions: Whether to include PFAS definitions (True for accuracy, False for specificity)
            limit_effective_graph_resistance: Limit or disable graph resistance computation (None=all, False/0=skip)
            compute_component_metrics: Whether to compute component graph metrics
        """
        smiles = test_case['smiles']
        expected_group_id = test_case['group_id']
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
                            # Forward all summary metrics from parse_mol
                            _SUMMARY_KEYS = (
                                'mean_branching', 'total_branching', 'sum_component_branching_ratio',
                                'mean_smarts_centrality', 'mean_component_fraction', 'total_components_fraction',
                                'mean_eccentricity', 'median_eccentricity',
                                'mean_diameter', 'mean_radius',
                                'mean_effective_graph_resistance', 'mean_effective_graph_resistance_BDE',
                                'mean_dist_to_barycentre', 'mean_dist_to_centre', 'mean_dist_to_periphery',
                            )
                            all_matches.append({
                                'type': 'group',
                                'id': match['id'],
                                'match_id': match.get('match_id'),
                                'name': match.get('group_name'),
                                'match_count': match.get('match_count'),
                                'components_sizes': match.get('components_sizes', []),
                                'num_components': match.get('num_components', 0),
                                'components_types': match.get('components_types', []),
                                **{k: match.get(k) for k in _SUMMARY_KEYS},
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
                all_defs = len(PFASGroups_result['detected_definitions']) ==5
                PFASGroups_result['matches'] = all_matches
                PFASGroups_result['success'] = expected_group_id in group_ids
                
        except Exception as e:
            PFASGroups_result['error'] = str(e)
        finally:
            PFASGroups_result['execution_time'] = time.perf_counter() - start_time
        
        return PFASGroups_result, all_defs
    
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
                atlas_result['success'] = predictions[0] not in  ['Other PFASs', "Unknown"]
            elif predictions and len(predictions) >= 1:
                atlas_result['first_class'] = predictions[0]
                atlas_result['success'] = predictions[0] not in  ['Other PFASs', "Unknown"]  # Consider it a success if it is classified as any PFAS class
                    
        except Exception as e:
            atlas_result['error'] = str(e)
        finally:
            atlas_result['execution_time'] = time.perf_counter() - start_time
        
        return atlas_result
    def run_tests(self):
        rdBase.DisableLog('rdApp.warning')  # Suppress RDKit warnings during tests
        for tc in self.test_cases:
            pf_result, all_defs = self.test_with_PFASGroups(tc)
            tc['PFASGroups_result'] = pf_result
            tc['all_definitions_detected'] = all_defs
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
        serialised_results = []
        for item in results:
            if isinstance(item, dict):
                serialised_item = {}
                for k, v in item.items():
                    if isinstance(v, (str, int, float, bool)) or v is None:
                        serialised_item[k] = v
                    else:
                        serialised_item[k] = str(v)  # Convert non-serializable items to string
                serialised_results.append(serialised_item)
            else:
                serialised_results.append(str(item))
        with open(filename, 'w') as f:
            json.dump(serialised_results, f, indent=2)
        
if __name__ == "__main__":
    seed = 2025
    tester = TestBranching(n=300, branching_range=(0.0, 0.6), min_carbons=5, max_carbons=15, seed=seed,
                           n_highly_branched_failing_proportion=0.25)
    tester.run_tests()
    tester.save(os.path.join(DATA_DIR, f"branching_test_results_{seed}.json"))
