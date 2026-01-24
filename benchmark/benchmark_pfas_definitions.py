"""
Comprehensive PFAS Definitions Benchmarking System

This module provides thorough benchmarking for the 5 PFAS definitions:
1. OECD Definition
2. EU PFAS Restriction  
3. OPPT 2023
4. UK PFAS Definition
5. PFASTRUCTv5

The benchmark evaluates:
- Accuracy: Correct identification of known PFAS
- Specificity: Correct rejection of non-PFAS
- Sensitivity: Detection of edge cases and borderline compounds
- Inter-definition agreement: Overlap and divergence between definitions
- Performance: Speed and scalability
"""

import json
import csv
import math
import time
import sys
import os
from datetime import datetime
from collections import defaultdict, Counter
from typing import List, Dict, Tuple, Optional
# Add parent directory to path
script_dir = os.path.dirname(os.path.abspath(__file__))
parent_dir = os.path.dirname(script_dir)
sys.path.append(parent_dir)


class PFASDefinitionBenchmark:
    """
    Comprehensive benchmarking system for PFAS definitions.
    
    Tests each of the 5 PFAS definitions across multiple dimensions:
    - True positives (known PFAS compounds)
    - True negatives (non-PFAS compounds)
    - Edge cases (borderline compounds)
    - Inter-definition concordance
    - Performance metrics
    """
    
    def __init__(self, load_definitions: bool = True):
        """Initialize benchmark with test molecule datasets"""
        
        self.load_definitions = load_definitions

        # Load PFAS definitions (lazy import to allow export without RDKit/Numpy)
        if load_definitions:
            from PFASgroups.PFASDefinitionModel import PFASDefinition
            definitions_path = os.path.join(parent_dir, 'PFASgroups', 'data', 'PFAS_definitions_smarts.json')
            with open(definitions_path, 'r') as f:
                definitions_data = json.load(f)
            
            self.definitions = {d['id']: PFASDefinition(**d) for d in definitions_data}
            self.definition_names = {d['id']: d['name'] for d in definitions_data}
        else:
            self.definitions = {}
            self.definition_names = {}
        
        # Initialize test molecule sets (CSV source of truth with defaults fallback)
        self.test_compounds_csv_path = os.path.join(script_dir, 'benchmark_test_compounds.csv')
        self._initialize_test_datasets()
        
    def _initialize_default_test_datasets(self):
        """Initialize curated default (hard-coded) test molecule datasets"""
        
        # 1. TRUE POSITIVES: Known PFAS compounds (should be detected by most definitions)
        self.true_positives = {
            'perfluoroalkyl_carboxylic_acids': [
                ('OC(=O)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)F', 'PFOA (C8)'),
                ('OC(=O)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)F', 'PFHxA (C6)'),
                ('OC(=O)C(F)(F)C(F)(F)C(F)(F)F', 'PFBA (C4)'),
                ('OC(=O)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)F', 'PFDA (C10)'),
                ('OC(=O)C(F)(F)F', 'TFA (C2)'),
            ],
            'perfluoroalkyl_sulfonic_acids': [
                ('FC(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)S(=O)(=O)O', 'PFOS (C8)'),
                ('FC(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)S(=O)(=O)O', 'PFHxS (C6)'),
                ('FC(F)(F)C(F)(F)C(F)(F)C(F)(F)S(=O)(=O)O', 'PFBS (C4)'),
                ('FC(F)(F)C(F)(F)S(=O)(=O)O', 'PFEtS (C2)'),
            ],
            'fluorotelomer_alcohols': [
                ('OCCC(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)F', '6:2 FTOH'),
                ('OCCC(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)F', '8:2 FTOH'),
                ('OCCC(F)(F)C(F)(F)C(F)(F)F', '4:2 FTOH'),
                ('OCCC(F)(F)C(F)(F)C(F)(F)C(F)(F)F', '5:2 FTOH'),
            ],
            'perfluoroalkyl_sulfonamides': [
                ('FC(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)S(=O)(=O)N', 'PFOSA'),
                ('FC(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)S(=O)(=O)NC', 'N-MeFOSA'),
                ('FC(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)S(=O)(=O)NCC', 'N-EtFOSA'),
            ],
            'perfluoroalkyl_iodides': [
                ('FC(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)I', 'C8F17I'),
                ('FC(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)I', 'C6F13I'),
            ],
            'side_chain_fluorinated': [
                ('FC(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)CCOC(=O)C', '8:2 fluorotelomer acrylate'),
                ('FC(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)CCOP(=O)(O)O', '6:2 fluorotelomer phosphate'),
            ],
            'cyclic_perfluorinated': [
                ('FC1(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C1(F)F', 'Perfluorocyclohexane'),
                ('FC1(F)C(F)(F)C(F)(F)C(F)(F)C1(F)F', 'Perfluorocyclopentane'),
            ],
            'perfluoroethers': [
                ('FC(F)(F)C(F)(F)OC(F)(F)C(F)(F)F', 'Perfluorodiethyl ether'),
                ('FC(F)(F)OC(F)(F)C(F)(F)OC(F)(F)F', 'PFPE'),
            ],
        }
        
        # 2. TRUE NEGATIVES: Non-PFAS compounds (should NOT be detected)
        self.true_negatives = {
            'non_fluorinated_analogs': [
                ('OC(=O)CCCCCCC', 'Octanoic acid (non-fluorinated analog of PFOA)'),
                ('CCCCCCCCS(=O)(=O)O', 'Octane sulfonic acid (non-fluorinated analog of PFOS)'),
                ('OCCCCCCCC', 'Octanol (non-fluorinated analog of FTOH)'),
            ],
            'monofluorinated': [
                ('OC(=O)CF', 'Fluoroacetic acid'),
                ('FCC', 'Fluoroethane'),
            ],
            'difluorinated_non_geminal': [
                ('FCCCF', '1,4-difluorobutane'),
                ('FC(C)CF', '1,3-difluoropropane'),
            ],
            'pharmaceutical_fluorinated': [
                ('c1ccc(F)cc1', 'Fluorobenzene'),
            ],
            'EU_specific_exclusions': [
                # Compounds that EU definition should specifically reject due to exclusion criteria
                # These have CF2/CF3 but are bonded to excluded groups (-OR, -NRR', H, halogens)
                ('C(F)(F)(F)Cl', 'Trifluoromethyl chloride (CF3-Cl) - halogen exclusion'),
                ('C(F)(F)(F)Br', 'Trifluoromethyl bromide (CF3-Br) - halogen exclusion'),
                ('C(F)(F)(F)I', 'Trifluoromethyl iodide (CF3-I) - halogen exclusion'),
                ('C(F)(F)(F)[H]', 'Fluoroform (CHF3) - H exclusion'),
                ('C(F)(F)(F)O', 'Trifluoromethanol (CF3-OH) - alcohol exclusion'),
                ('C(F)(F)(F)OC', 'Trifluoromethyl methyl ether (CF3-OCH3) - ether exclusion'),
                ('C(F)(F)(F)N', 'Trifluoromethylamine (CF3-NH2) - amine exclusion'),
                ('C(F)(F)(Cl)C', '1-Chloro-1,1-difluoroethane (CF2Cl-CH3) - halogen exclusion'),
                ('C(F)(F)([H])C', '1,1-Difluoroethane (CF2H-CH3) - H exclusion'),
                ('C(F)(F)(OC)C', '1,1-Difluoro-1-methoxyethane (CF2-OCH3-CH3) - ether exclusion'),
                ('C(F)(F)(N)C', '1,1-Difluoro-1-aminoethane (CF2-NH2-CH3) - amine exclusion'),
                ('FC1(F)Oc2ccccc2O1', 'Benzodioxole with CF2 - aromatic ether exclusion'),
            ],
        }
        
        # 3. EDGE CASES: Borderline compounds where definitions may disagree
        self.edge_cases = {
            'short_chain_perfluorinated': [
                ('FC(F)(F)C(F)(F)F', 'Perfluoropropane (C3)'),
                ('OC(=O)C(F)(F)F', 'TFA - C2 PFCA'),
                ('FC(F)(F)S(=O)(=O)O', 'C1 perfluorosulfonic acid'),
            ],
            'polyfluorinated_not_perfluorinated': [
                ('OC(=O)CC(F)(F)C(F)(F)F', 'H-PFBA (one non-fluorinated carbon)'),
                ('FC(F)C(F)(F)C(F)(F)F', 'Polyfluorobutane (missing one F)'),
                ('OC(=O)C(F)(F)C(F)(F)CC', 'Semi-fluorinated carboxylic acid'),
            ],
            'fluorine_rich_aromatics': [
                ('FC(F)(F)c1c(F)c(F)c(F)c(F)c1F', 'Perfluorotoluene'),
                ('c1c(F)c(F)c(F)c(F)c1F', 'Pentafluorobenzene'),
            ],
            'heteroatom_connected': [
                ('FC(F)(F)C(F)(F)OC', 'Perfluoroethyl methyl ether - OECD may differ from EU'),
                ('FC(F)(F)C(F)(F)NC', 'Perfluoroethylmethylamine - OECD may differ from EU'),
            ],
            'halogen_substituted': [
                ('ClC(F)(F)C(F)(F)F', 'Chloroperfluoropropane - test halogen exclusions'),
                ('FC(F)(F)C(F)(F)C(Cl)(F)F', 'Chlorine on perfluoroalkyl - distance matters'),
            ],
            'ultra_short_chain': [
                ('FC(F)(F)F', 'CF4 - Carbon tetrafluoride'),
                ('FC(F)(F)C(F)(F)F', 'C2F6 - Hexafluoroethane'),
            ],
            'EU_specific_edge_cases': [
                # Edge cases specific to EU definition exclusions
                ('FC(F)(F)C', 'Trifluoromethyl-methane - short chain with CF3'),
                ('FC(F)(OC(=O)C)C', 'CF2 with acetoxy - EU should reject, others may accept'),
                ('FC(F)(SC)C', 'CF2 with methylthio - EU exclusion for -SR'),
                ('FC(F)(F)C(F)(F)SC', 'CF3-CF2-SCH3 - sulfur at distance, should differ'),
                ('FC(F)=C(F)OC(F)(F)F', 'Fluorovinyl ether with CF3 - unsaturated edge case'),
                ('FC(F)(N(C)C)C(F)(F)F', 'CF2 with tertiary amine - EU should reject'),
                ('FC(F)(F)C(N(C)C)C(F)(F)F', 'Amine between CF3 groups - EU may accept'),
            ],
        }
        
        # 4. DEFINITION-SPECIFIC TEST CASES
        self.definition_specific_tests = {
            'OECD_specific': [
                # OECD requires CF3 or CF2 not bonded to H, Cl, Br, I
                ('FC(F)(F)Cl', 'CF3Cl - should NOT match OECD (CF3 bonded to Cl)', False),
                ('FC(F)(F)C(F)(F)F', 'C2F6 - should match OECD', True),
                ('FC(F)C(F)(F)F', 'CF2H-CF3 - CF2 bonded to H, may not match', False),
            ],
            'EU_specific': [
                # EU has more exclusion criteria than OECD - comprehensive tests from EU SMARTS construction
                # Based on actual EU definition: excludes CF2/CF3 bonded to -OR, -NRR', H, or halogens
                
                # CF2 exclusions - should NOT match EU definition
                ('FC(F)(Cl)C', 'CF2 bonded to Cl - should NOT match', False),
                ('FC(F)(Br)C', 'CF2 bonded to Br - should NOT match', False),
                ('FC(F)(I)C', 'CF2 bonded to I - should NOT match', False),
                ('FC(F)([H])C', 'CF2 bonded to H - should NOT match', False),
                ('FC(F)(O)C', 'CF2 bonded to -OH - should NOT match', False),
                ('FC(F)(OC)C', 'CF2 bonded to -OCH3 - should NOT match', False),
                ('FC(F)(OCC)C', 'CF2 bonded to -OEt - should NOT match', False),
                ('FC(F)(Oc1ccccc1)C', 'CF2 bonded to -OPh - should NOT match', False),
                ('FC(F)(OC(=O)C)C', 'CF2 bonded to -OCOCH3 - should NOT match', False),
                ('FC(F)(N)C', 'CF2 bonded to -NH2 - should NOT match', False),
                ('FC(F)(NC)C', 'CF2 bonded to -NHCH3 - should NOT match', False),
                ('FC(F)(N(C)C)C', 'CF2 bonded to -N(CH3)2 - should NOT match', False),
                
                # CF2 with acceptable substituents - should match EU definition
                ('FC(F)(C(F)(F)F)C', 'CF2 with perfluoro chain - should match', True),
                ('FC(F)(P)C', 'CF2 bonded to P - should match', True),
                ('FC(F)(F)C(F)(F)C(F)(P)OC', 'CF2 chain with P and -OR at distance - should match', True),
                
                # CF3 exclusions - should NOT match EU definition
                ('FC(F)(F)Cl', 'CF3 bonded to Cl - should NOT match', False),
                ('FC(F)(F)Br', 'CF3 bonded to Br - should NOT match', False),
                ('FC(F)(F)I', 'CF3 bonded to I - should NOT match', False),
                ('FC(F)(F)[H]', 'CF3 bonded to H (CHF3) - should NOT match', False),
                ('FC(F)(F)O', 'CF3 bonded to -OH - should NOT match', False),
                ('FC(F)(F)OC', 'CF3 bonded to -OCH3 - should NOT match', False),
                ('FC(F)(F)OCC', 'CF3 bonded to -OEt - should NOT match', False),
                ('FC(F)(F)Oc1ccccc1', 'CF3 bonded to -OPh - should NOT match', False),
                ('FC(F)(F)OC(=O)C', 'CF3 bonded to -OCOCH3 - should NOT match', False),
                ('FC(F)(F)N', 'CF3 bonded to -NH2 - should NOT match', False),
                ('FC(F)(F)NC', 'CF3 bonded to -NHCH3 - should NOT match', False),
                ('FC(F)(F)N(C)C', 'CF3 bonded to -N(CH3)2 - should NOT match', False),
                
                # CF3 with acceptable substituents - should match EU definition
                ('FC(F)(F)F', 'CF4 - should match (no exclusions)', True),
                ('FC(F)(F)C(F)(F)F', 'Perfluoroethane - should match', True),
                ('FC(F)(F)C(F)(F)O', 'CF3-CF2-OH (ether at distance) - should match', True),
                ('FC(F)(F)C(F)(F)N', 'CF3-CF2-NH2 (amine at distance) - should match', True),
                
                # Complex molecules with multiple functional groups
                ('FC(F)=C(F)OC(F)(F)OC(F)=C(F)F', 'PFPE with ether - should match (ether not directly on CF2/CF3)', True),
                ('FC1(F)Oc2cccc(C=O)c2O1', 'Benzodioxole with CF2 bonded to aromatic O - should NOT match', False),
                
                # Real pharmaceutical example (from test file)
                ('[H]OC([H])([H])C([H])(O[H])C([H])([H])N1C2=C([H])C(F)=C(N([H])C(=O)C3(C4=C([H])C5=C(OC(F)(F)O5)C([H])=C4[H])C([H])([H])C3([H])[H])C([H])=C2C([H])=C1C(C([H])([H])[H])(C([H])([H])[H])C([H])([H])O[H]', 
                 'Complex pharma with CF2 in benzodioxole - should NOT match (CF2 bonded to aromatic O)', False),
            ],
            'OPPT_specific': [
                # OPPT has 3 inclusion criteria
                ('FC(F)(F)C(F)(F)OC(F)(F)C(F)(F)F', 'PFPE - matches criterion 2', True),
                ('FC(F)(F)C(C(F)(F)F)(C(F)(F)F)C(F)(F)F', 'Branched perfluoro - criterion 3', True),
            ],
            'UK_specific': [
                # UK is similar to OECD but different exclusions
                ('FC(F)(F)C(F)(F)F', 'Should match UK', True),
            ],
            'PFASTRUCTv5_specific': [
                # Has fluorine ratio criterion (>= 0.3) OR structural patterns
                # Basic threshold tests
                ('FC(F)(F)C(F)(F)CCCCCCCC', 'Lower F ratio - test threshold', True),
                ('FC(F)(F)CCCCCCCCCCCC', 'Very low F ratio (~0.15)', False),
                
                # Validation test set from PFASTRUCTv5 reference (23 compounds)
                # TRUE POSITIVES (should be detected by PFASTRUCTv5)
                ('C(=O)(C(F)(F)Cl)C(F)(F)Cl', '1,3-Dichloro-1,1,3,3-tetrafluoropropan-2-one (DTXSID1073159)', True),
                ('COC(C(OC(C(=O)C(F)(F)F)(F)F)F)(F)F', '1,1,1,3,3-Pentafluoro-3-(1,2,2-trifluoro-2-methoxyethoxy)propan-2-one (DTXSID40748237)', True),
                ('C(C(=O)CC(C(=O)C(F)(F)F)(F)F)C(C(=O)C(F)(F)F)(F)F', '2,5,8-Nonanetrione, decafluoro- (DTXSID60435490)', True),
                ('C(C(CC(F)(F)F)(F)F)C(CC(F)(F)F)(F)F', '1,1,1,3,3,5,5,7,7,7-Decafluoroheptane (DTXSID20776331)', True),
                ('C(=C(F)F)(OC(OC(=C(F)F)F)(F)F)F', '1-(Difluoro[(trifluoroethenyl)oxy]methoxy)-1,2,2-trifluoroethene (DTXSID00896607)', True),
                ('C(=C(F)F)(C(C(=C(F)F)F)(F)F)F', '1,1,2,3,3,4,5,5-Octafluoropenta-1,4-diene (DTXSID30529779)', True),
                ('C1(=C(C(=C(C(=C1F)F)F)F)F)C2=C(C(C(=C(C2(F)F)F)F)(F)F)F', 'Dodecafluoro-2,5-dihydro-1,1-biphenyl (DTXSID20525346)', True),
                ('C(C(=O)F)(F)F', '2,2-Difluoroacetyl fluoride (DTXSID70381191)', True),
                ('C(=C(C(F)(F)Cl)F)(C(F)(F)Cl)F', '(2E)-1,4-Dichloro-hexafluoro-2-butene (DTXSID601020902)', True),
                ('C(C(C(C(F)F)F)F)F', '1,1,2,3,4-Pentafluorobutane (DTXSID50931113)', True),
                ('CCN(CC)C1(OCCCO1)C(C(F)(F)F)F', 'N,N-Diethyl-2-(1,2,2,2-tetrafluoroethyl)-1,3-dioxan-2-amine (DTXSID40519722)', True),
                ('COC(CC(=C(F)F)C(F)(F)F)SC1=CC=CC=C1', '[4,4-Difluoro-1-methoxy-3-(trifluoromethyl)but-3-en-1-yl]sulfanyl)benzene (DTXSID20809173)', True),
                ('C1=CC(=CC(=C1)Cl)NOC(C(C(=O)C(F)(F)F)F)(F)F', '4-[(3-Chloroanilino)oxy]-1,1,1,3,4,4-hexafluorobutan-2-one (DTXSID201023166)', True),
                ('CC(C)(CS(=O)(=O)O)NC(=O)CCSCCC(C(C(C(F)(F)F)(F)F)(F)F)(F)F', '4:2 Fluorotelomer thioether amido sulfonic acid (DTXSID00892528)', True),
                ('CN(C)CCCN(CCC(=O)O)S(=O)(=O)C(C(C(C(F)(F)F)(F)F)(F)F)(F)F', 'N-(Perfluorobutanesulfonyl)-N-(3-dimethylaminopropyl)-3-aminopropanoic acid (DTXSID20882022)', True),
                ('C1=CC(=C(C=C1N)Cl)OC(C(OC(F)(F)F)F)(F)F', '3-Chloro-4-[1,1,2-trifluoro-2-(trifluoromethoxy)ethoxy]aniline (DTXSID10576660)', True),
                ('C1=CC=C2C(=C1)C(OC(O2)(C(F)(F)F)F)(F)F', '2,4,4-Trifluoro-2-(trifluoromethyl)-2H,4H-1,3-benzodioxine (DTXSID80661330)', True),
                ('CCCOC(=O)C(=C(F)F)C(F)(F)F', 'Propyl 3,3-difluoro-2-(trifluoromethyl)prop-2-enoate (DTXSID80804225)', True),
                ('C[N+](C)(CCCNS(=O)(=O)CCC(C(C(C(C(C(F)(F)F)(F)F)(F)F)(F)F)(F)F)(F)F)CC(=O)[O-]', '6:2 Fluorotelomer sulfonamide betaine (DTXSID4041284)', True),
                
                # TRUE NEGATIVES (should NOT be detected by PFASTRUCTv5)
                ('CCCCC(CC(C(F)(F)Br)(F)Cl)Br', '1,4-dibromo-2-chloro-1,1,2-trifluorooctane - low F ratio (DTXSID10382122)', False),
                ('C(C(OC(Cl)(Cl)Cl)(F)F)(F)Cl', '2-Chloro-1,1,2-trifluoro-1-(trichloromethoxy)ethane - low F ratio (DTXSID30963477)', False),
                ('C1=COC(=C1)C=C(S(=O)(=O)C(F)(F)F)S(=O)(=O)C(F)(F)F', '2-[2,2-Bis(trifluoromethanesulfonyl)ethenyl]furan - borderline (DTXSID10711579)', False),
                ('COC(=C(F)Cl)F', '1-Chloro-1,2-difluoro-2-methoxyethene - too few fluorines (DTXSID50777992)', False),
            ],
        }
        
        # 5. COMPARISON PAIRS: Known reference compounds for inter-definition comparison
        self.reference_compounds = {
            'legacy_pfas': [
                ('OC(=O)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)F', 'PFOA'),
                ('FC(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)S(=O)(=O)O', 'PFOS'),
            ],
            'short_chain_replacement': [
                ('OC(=O)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)F', 'PFHxA (C6)'),
                ('FC(F)(F)C(F)(F)C(F)(F)C(F)(F)S(=O)(=O)O', 'PFBS (C4)'),
            ],
            'telomer_based': [
                ('OCCC(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)F', '6:2 FTOH'),
                ('OCCC(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)F', '8:2 FTOH'),
            ],
        }

    def _initialize_test_datasets(self):
        """Load test molecules from CSV; fall back to defaults if missing."""
        # Build defaults first (used for fallback/export)
        self._initialize_default_test_datasets()
        self._default_test_compounds_rows = self._flatten_test_compounds_for_export(source='current')

        # Load from CSV if present
        if os.path.exists(self.test_compounds_csv_path):
            self._load_test_compounds_from_csv(self.test_compounds_csv_path)
        else:
            # Keep defaults and export a fresh CSV for user edits
            self.export_test_compounds_csv(self.test_compounds_csv_path, source='default')

    def _reset_test_datasets(self):
        """Reset all test datasets to empty containers."""
        self.true_positives = defaultdict(list)
        self.true_negatives = defaultdict(list)
        self.edge_cases = defaultdict(list)
        self.definition_specific_tests = defaultdict(list)
        self.reference_compounds = defaultdict(list)

    def _parse_expected_result(self, value) -> Optional[bool]:
        if value is None or (isinstance(value, float) and math.isnan(value)):
            return None
        text = str(value).strip().lower()
        if text in {'true', 'yes', '1'}:
            return True
        if text in {'false', 'no', '0'}:
            return False
        if text in {'variable', 'na', 'n/a', 'unknown', ''}:
            return None
        return None

    def _load_test_compounds_from_csv(self, csv_path: str):
        """Load test compounds from CSV and overwrite in-memory datasets."""
        with open(csv_path, newline='', encoding='utf-8') as f:
            reader = csv.DictReader(f)
            if reader.fieldnames is None:
                raise ValueError(f"No header row found in {csv_path}")

            normalized_cols = {c.strip().lower(): c for c in reader.fieldnames}
            required = ['smiles', 'test for definition', 'reason', 'expected_result']
            missing = [c for c in required if c not in normalized_cols]
            if missing:
                raise ValueError(f"Missing required columns in {csv_path}: {', '.join(missing)}")

            smiles_col = normalized_cols['smiles']
            test_col = normalized_cols['test for definition']
            reason_col = normalized_cols['reason']
            expected_col = normalized_cols['expected_result']

            self._reset_test_datasets()

            def_map = {
                'oecd': 'OECD_specific',
                'eu': 'EU_specific',
                'oppt': 'OPPT_specific',
                'uk': 'UK_specific',
                'pfastructv5': 'PFASTRUCTv5_specific'
            }

            for row in reader:
                smiles = str(row.get(smiles_col, '')).strip()
                if not smiles:
                    continue
                test_for = str(row.get(test_col, '')).strip()
                reason = str(row.get(reason_col, '')).strip() or 'Unnamed test'
                expected = self._parse_expected_result(row.get(expected_col, ''))

                test_key = test_for.strip().lower()

                if test_key in {'all definitions (general)', 'all definitions', 'general'}:
                    if expected is True:
                        self.true_positives['general'].append((smiles, reason))
                        self.reference_compounds['general'].append((smiles, reason))
                    elif expected is False:
                        self.true_negatives['general'].append((smiles, reason))
                    else:
                        self.edge_cases['general'].append((smiles, reason))
                elif test_key in {'edge case', 'edge cases'}:
                    self.edge_cases['edge_cases'].append((smiles, reason))
                elif test_key in def_map:
                    if expected is None:
                        continue
                    self.definition_specific_tests[def_map[test_key]].append((smiles, reason, expected))
                else:
                    self.edge_cases['unclassified'].append((smiles, reason))

    def _flatten_test_compounds_for_export(self, source: str = 'current') -> List[Dict[str, str]]:
        """Flatten test datasets into CSV rows."""
        rows: List[Dict[str, str]] = []

        def add_row(smiles: str, test_for: str, reason: str, expected: Optional[bool]):
            if expected is True:
                expected_text = 'True'
            elif expected is False:
                expected_text = 'False'
            else:
                expected_text = 'Variable'
            rows.append({
                'smiles': smiles,
                'test for definition': test_for,
                'reason': reason,
                'expected_result': expected_text
            })

        for _, molecules in self.true_positives.items():
            for smiles, name in molecules:
                add_row(smiles, 'All definitions (general)', name, True)

        for _, molecules in self.true_negatives.items():
            for smiles, name in molecules:
                add_row(smiles, 'All definitions (general)', name, False)

        for _, molecules in self.edge_cases.items():
            for smiles, name in molecules:
                add_row(smiles, 'Edge case', name, None)

        def_label_map = {
            'OECD_specific': 'OECD',
            'EU_specific': 'EU',
            'OPPT_specific': 'OPPT',
            'UK_specific': 'UK',
            'PFASTRUCTv5_specific': 'PFASTRUCTv5'
        }
        for def_key, tests in self.definition_specific_tests.items():
            label = def_label_map.get(def_key, def_key)
            for smiles, name, expected in tests:
                add_row(smiles, label, name, expected)

        return rows

    def export_test_compounds_csv(self, output_path: Optional[str] = None, source: str = 'current') -> str:
        """Export test compounds to CSV for manual editing."""
        output_path = output_path or self.test_compounds_csv_path

        if source == 'default' and hasattr(self, '_default_test_compounds_rows'):
            rows = self._default_test_compounds_rows
        else:
            rows = self._flatten_test_compounds_for_export(source=source)

        with open(output_path, 'w', newline='', encoding='utf-8') as f:
            writer = csv.DictWriter(f, fieldnames=[
                'smiles', 'test for definition', 'reason', 'expected_result'
            ])
            writer.writeheader()
            writer.writerows(rows)

        return output_path
        
    def test_single_molecule(self, smiles: str, expected_definitions: List[int] = None) -> Dict:
        """
        Test a single molecule against all PFAS definitions.
        
        Args:
            smiles: SMILES string of molecule
            expected_definitions: Optional list of definition IDs expected to match
            
        Returns:
            Dictionary with test results
        """
        from rdkit import Chem
        from rdkit.Chem import Descriptors, rdMolDescriptors
        from PFASgroups.core import parse_mol

        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return {
                'smiles': smiles,
                'valid': False,
                'error': 'Invalid SMILES'
            }
        
        # Get molecular properties
        mol_with_h = Chem.AddHs(mol)
        formula = rdMolDescriptors.CalcMolFormula(mol_with_h)
        mol_weight = Descriptors.MolWt(mol)
        num_atoms = mol.GetNumAtoms()
        num_fluorines = sum(1 for atom in mol.GetAtoms() if atom.GetSymbol() == 'F')
        fluorine_ratio = num_fluorines / max(num_atoms, 1)
        
        # Test with parse_mol (includes all definitions)
        start_time = time.perf_counter()
        result = parse_mol(mol, include_PFAS_definitions=True)
        execution_time = time.perf_counter() - start_time
        
        # Extract detected definitions
        detected_definitions = []
        if isinstance(result, dict) and 'matches' in result:
            for match in result['matches']:
                if match.get('type') == 'PFASdefinition':
                    detected_definitions.append({
                        'id': match['id'],
                        'name': self.definition_names.get(match['id'], 'Unknown'),
                        'match_id': match.get('match_id')
                    })
        
        detected_ids = [d['id'] for d in detected_definitions]
        
        # Evaluate correctness if expected_definitions provided
        correctness = None
        if expected_definitions is not None:
            true_positives = [d for d in expected_definitions if d in detected_ids]
            false_positives = [d for d in detected_ids if d not in expected_definitions]
            false_negatives = [d for d in expected_definitions if d not in detected_ids]
            
            correctness = {
                'true_positives': true_positives,
                'false_positives': false_positives,
                'false_negatives': false_negatives,
                'precision': len(true_positives) / max(len(detected_ids), 1),
                'recall': len(true_positives) / max(len(expected_definitions), 1),
            }
            if correctness['precision'] + correctness['recall'] > 0:
                correctness['f1_score'] = 2 * (correctness['precision'] * correctness['recall']) / \
                                         (correctness['precision'] + correctness['recall'])
            else:
                correctness['f1_score'] = 0.0
        
        return {
            'smiles': smiles,
            'valid': True,
            'molecular_properties': {
                'formula': formula,
                'molecular_weight': mol_weight,
                'num_atoms': num_atoms,
                'num_fluorines': num_fluorines,
                'fluorine_ratio': fluorine_ratio,
            },
            'detected_definitions': detected_definitions,
            'detected_ids': detected_ids,
            'expected_ids': expected_definitions,
            'correctness': correctness,
            'execution_time': execution_time,
        }
    
    def run_true_positives_benchmark(self) -> Dict:
        """Test all true positive cases (known PFAS)"""
        print("\n" + "="*80)
        print("TRUE POSITIVES BENCHMARK: Known PFAS Compounds")
        print("="*80)
        
        results = {
            'category': 'true_positives',
            'subcategories': {},
            'summary': {
                'total_molecules': 0,
                'total_detections': 0,
                'by_definition': {i: 0 for i in range(1, 6)},
            },
            'false_negative_details': []  # Store detailed false negative info
        }
        
        for category, molecules in self.true_positives.items():
            print(f"\n[CATEGORY] {category.replace('_', ' ').title()}")
            print(f"   Testing {len(molecules)} molecules...")
            
            category_results = []
            for smiles, name in molecules:
                # For true positives, we expect at least one definition to match
                test_result = self.test_single_molecule(smiles, expected_definitions=None)
                test_result['name'] = name
                test_result['category'] = category
                
                category_results.append(test_result)
                
                if test_result['valid']:
                    results['summary']['total_molecules'] += 1
                    num_detected = len(test_result['detected_ids'])
                    if num_detected > 0:
                        results['summary']['total_detections'] += 1
                    
                    for def_id in test_result['detected_ids']:
                        results['summary']['by_definition'][def_id] += 1
                    
                    # Print summary and track false negatives
                    detected_names = [d['name'] for d in test_result['detected_definitions']]
                    if detected_names:
                        print(f"   [OK] {name}: {', '.join(detected_names)}")
                    else:
                        print(f"   [FAIL] {name}: No definitions matched")
                        # Store false negative details
                        results['false_negative_details'].append({
                            'smiles': smiles,
                            'molecule_name': name,
                            'category': category,
                            'expected_to_detect': 'All definitions (true positive)',
                            'actually_detected': [],
                            'type': 'false_negative'
                        })
                else:
                    print(f"   [WARN] {name}: Invalid molecule")
            
            results['subcategories'][category] = category_results
        
        # Calculate summary statistics
        total_mols = results['summary']['total_molecules']
        detection_rate = (results['summary']['total_detections'] / max(total_mols, 1)) * 100
        
        print(f"\n{'='*80}")
        print(f"TRUE POSITIVES SUMMARY:")
        print(f"  - Total molecules tested: {total_mols}")
        print(f"  - Overall detection rate: {detection_rate:.1f}%")
        print(f"  - Detection by definition:")
        for def_id in range(1, 6):
            count = results['summary']['by_definition'][def_id]
            rate = (count / max(total_mols, 1)) * 100
            print(f"    - {self.definition_names[def_id]}: {count}/{total_mols} ({rate:.1f}%)")
        
        return results
    
    def run_true_negatives_benchmark(self) -> Dict:
        """Test all true negative cases (non-PFAS)"""
        print("\n" + "="*80)
        print("TRUE NEGATIVES BENCHMARK: Non-PFAS Compounds")
        print("="*80)
        
        results = {
            'category': 'true_negatives',
            'subcategories': {},
            'summary': {
                'total_molecules': 0,
                'correct_rejections': 0,
                'false_positives': 0,
                'by_definition': {i: 0 for i in range(1, 6)},
            },
            'false_positive_details': []  # Store detailed false positive info
        }
        
        for category, molecules in self.true_negatives.items():
            print(f"\n[CATEGORY] {category.replace('_', ' ').title()}")
            print(f"   Testing {len(molecules)} molecules...")
            
            category_results = []
            for smiles, name in molecules:
                # For true negatives, we expect NO definitions to match
                test_result = self.test_single_molecule(smiles, expected_definitions=[])
                test_result['name'] = name
                test_result['category'] = category
                
                category_results.append(test_result)
                
                if test_result['valid']:
                    results['summary']['total_molecules'] += 1
                    num_detected = len(test_result['detected_ids'])
                    
                    if num_detected == 0:
                        results['summary']['correct_rejections'] += 1
                        print(f"   [OK] {name}: Correctly rejected")
                    else:
                        results['summary']['false_positives'] += 1
                        detected_names = [d['name'] for d in test_result['detected_definitions']]
                        print(f"   [FAIL] {name}: False positive - {', '.join(detected_names)}")
                        
                        # Store detailed false positive info
                        results['false_positive_details'].append({
                            'smiles': smiles,
                            'molecule_name': name,
                            'category': category,
                            'detected_by': detected_names,
                            'detected_ids': test_result['detected_ids'],
                            'type': 'false_positive'
                        })
                        
                        for def_id in test_result['detected_ids']:
                            results['summary']['by_definition'][def_id] += 1
                else:
                    print(f"   [WARN] {name}: Invalid molecule")
            
            results['subcategories'][category] = category_results
        
        # Calculate summary statistics
        total_mols = results['summary']['total_molecules']
        specificity = (results['summary']['correct_rejections'] / max(total_mols, 1)) * 100
        
        print(f"\n{'='*80}")
        print(f"TRUE NEGATIVES SUMMARY:")
        print(f"  - Total molecules tested: {total_mols}")
        print(f"  - Specificity (correct rejections): {specificity:.1f}%")
        print(f"  - False positives: {results['summary']['false_positives']}")
        print(f"  - False positives by definition:")
        for def_id in range(1, 6):
            count = results['summary']['by_definition'][def_id]
            if count > 0:
                print(f"    - {self.definition_names[def_id]}: {count}")
        
        return results
    
    def run_edge_cases_benchmark(self) -> Dict:
        """Test edge cases where definitions may disagree"""
        print("\n" + "="*80)
        print("EDGE CASES BENCHMARK: Borderline Compounds")
        print("="*80)
        
        results = {
            'category': 'edge_cases',
            'subcategories': {},
            'summary': {
                'total_molecules': 0,
                'unanimous_decisions': 0,
                'split_decisions': 0,
                'disagreement_matrix': defaultdict(int),
            }
        }
        
        for category, molecules in self.edge_cases.items():
            print(f"\n[CATEGORY] {category.replace('_', ' ').title()}")
            print(f"   Testing {len(molecules)} molecules...")
            
            category_results = []
            for smiles, name in molecules:
                test_result = self.test_single_molecule(smiles, expected_definitions=None)
                test_result['name'] = name
                test_result['category'] = category
                
                category_results.append(test_result)
                
                if test_result['valid']:
                    results['summary']['total_molecules'] += 1
                    num_detected = len(test_result['detected_ids'])
                    
                    # Check for unanimous vs split decision
                    if num_detected == 0 or num_detected == 5:
                        results['summary']['unanimous_decisions'] += 1
                        decision_type = "All reject" if num_detected == 0 else "All accept"
                        print(f"     {name}: {decision_type}")
                    else:
                        results['summary']['split_decisions'] += 1
                        detected_names = [d['name'] for d in test_result['detected_definitions']]
                        print(f"    {name}: Split decision - {', '.join(detected_names)}")
                        
                        # Track disagreement patterns
                        detected_set = frozenset(test_result['detected_ids'])
                        results['summary']['disagreement_matrix'][detected_set] += 1
                else:
                    print(f"     {name}: Invalid molecule")
            
            results['subcategories'][category] = category_results
        
        # Calculate summary statistics
        total_mols = results['summary']['total_molecules']
        split_rate = (results['summary']['split_decisions'] / max(total_mols, 1)) * 100
        
        print(f"\n{'='*80}")
        print(f"EDGE CASES SUMMARY:")
        print(f"  - Total molecules tested: {total_mols}")
        print(f"  - Unanimous decisions: {results['summary']['unanimous_decisions']}")
        print(f"  - Split decisions: {results['summary']['split_decisions']} ({split_rate:.1f}%)")
        print(f"  - Common disagreement patterns:")
        
        # Show top disagreement patterns
        sorted_patterns = sorted(results['summary']['disagreement_matrix'].items(), 
                                key=lambda x: x[1], reverse=True)
        for pattern, count in sorted_patterns[:5]:
            if len(pattern) > 0:
                pattern_names = [self.definition_names[i] for i in sorted(pattern)]
                print(f"    - {count}x: {', '.join(pattern_names)}")
        
        return results
    
    def run_definition_specific_tests(self) -> Dict:
        """Test definition-specific edge cases"""
        print("\n" + "="*80)
        print("DEFINITION-SPECIFIC TESTS")
        print("="*80)
        
        results = {
            'category': 'definition_specific',
            'by_definition': {},
            'summary': {
                'total_tests': 0,
                'passed': 0,
                'failed': 0,
            }
        }
        
        for def_category, test_cases in self.definition_specific_tests.items():
            print(f"\n[CATEGORY] {def_category.replace('_', ' ').title()}")
            
            def_results = []
            for test_case in test_cases:
                if len(test_case) == 3:
                    smiles, name, expected_match = test_case
                    
                    test_result = self.test_single_molecule(smiles, expected_definitions=None)
                    test_result['name'] = name
                    test_result['expected_match'] = expected_match
                    
                    # Extract definition ID from category name
                    if 'OECD' in def_category:
                        target_def_id = 1
                    elif 'EU' in def_category:
                        target_def_id = 2
                    elif 'OPPT' in def_category:
                        target_def_id = 3
                    elif 'UK' in def_category:
                        target_def_id = 4
                    elif 'PFASTRUCTv5' in def_category:
                        target_def_id = 5
                    else:
                        target_def_id = None
                    
                    if test_result['valid'] and target_def_id:
                        actual_match = target_def_id in test_result['detected_ids']
                        test_result['actual_match'] = actual_match
                        test_result['test_passed'] = (actual_match == expected_match)
                        
                        results['summary']['total_tests'] += 1
                        if test_result['test_passed']:
                            results['summary']['passed'] += 1
                            print(f"   [OK] {name}: Expected={expected_match}, Actual={actual_match}")
                        else:
                            results['summary']['failed'] += 1
                            print(f"   [FAIL] {name}: Expected={expected_match}, Actual={actual_match}")
                    else:
                        test_result['test_passed'] = None
                        print(f"   [WARN] {name}: Invalid test")
                    
                    def_results.append(test_result)
            
            results['by_definition'][def_category] = def_results
        
        # Summary
        total_tests = results['summary']['total_tests']
        pass_rate = (results['summary']['passed'] / max(total_tests, 1)) * 100
        
        print(f"\n{'='*80}")
        print(f"DEFINITION-SPECIFIC SUMMARY:")
        print(f"  - Total tests: {total_tests}")
        print(f"  - Passed: {results['summary']['passed']} ({pass_rate:.1f}%)")
        print(f"  - Failed: {results['summary']['failed']}")
        
        return results
    
    def run_concordance_analysis(self) -> Dict:
        """Analyze inter-definition concordance"""
        import numpy as np
        print("\n" + "="*80)
        print("INTER-DEFINITION CONCORDANCE ANALYSIS")
        print("="*80)
        
        results = {
            'category': 'concordance',
            'compounds': [],
            'summary': {
                'total_compounds': 0,
                'agreement_matrix': np.zeros((5, 5), dtype=int),
                'jaccard_similarity': np.zeros((5, 5)),
            }
        }
        
        # Collect all reference compounds
        all_compounds = []
        for category, molecules in self.reference_compounds.items():
            for smiles, name in molecules:
                all_compounds.append((smiles, name, category))
        
        print(f"\n[ANALYSIS] Analyzing concordance across {len(all_compounds)} reference compounds...")
        
        # Test all compounds
        definition_detections = {i: [] for i in range(1, 6)}
        
        for smiles, name, category in all_compounds:
            test_result = self.test_single_molecule(smiles)
            if test_result['valid']:
                results['summary']['total_compounds'] += 1
                results['compounds'].append({
                    'smiles': smiles,
                    'name': name,
                    'category': category,
                    'detected_ids': test_result['detected_ids'],
                })
                
                # Track which definitions detected this compound
                for def_id in range(1, 6):
                    definition_detections[def_id].append(1 if def_id in test_result['detected_ids'] else 0)
        
        # Calculate pairwise agreement
        for i in range(1, 6):
            for j in range(1, 6):
                if i == j:
                    results['summary']['agreement_matrix'][i-1, j-1] = len(definition_detections[i])
                else:
                    # Count agreements
                    agreements = sum(1 for k in range(len(definition_detections[i])) 
                                   if definition_detections[i][k] == definition_detections[j][k])
                    results['summary']['agreement_matrix'][i-1, j-1] = agreements
                    
                    # Calculate Jaccard similarity
                    set_i = set(k for k, v in enumerate(definition_detections[i]) if v == 1)
                    set_j = set(k for k, v in enumerate(definition_detections[j]) if v == 1)
                    
                    if len(set_i.union(set_j)) > 0:
                        jaccard = len(set_i.intersection(set_j)) / len(set_i.union(set_j))
                    else:
                        jaccard = 1.0
                    
                    results['summary']['jaccard_similarity'][i-1, j-1] = jaccard
        
        # Print agreement matrix
        print(f"\n{'='*80}")
        print("AGREEMENT MATRIX (number of agreements on reference compounds):")
        print(f"{'':20}", end='')
        for i in range(1, 6):
            print(f"{self.definition_names[i][:12]:>15}", end='')
        print()
        
        for i in range(1, 6):
            print(f"{self.definition_names[i][:18]:20}", end='')
            for j in range(1, 6):
                print(f"{results['summary']['agreement_matrix'][i-1, j-1]:>15}", end='')
            print()
        
        # Print Jaccard similarity
        print(f"\nJACCARD SIMILARITY (0=no overlap, 1=perfect overlap):")
        print(f"{'':20}", end='')
        for i in range(1, 6):
            print(f"{self.definition_names[i][:12]:>15}", end='')
        print()
        
        for i in range(1, 6):
            print(f"{self.definition_names[i][:18]:20}", end='')
            for j in range(1, 6):
                print(f"{results['summary']['jaccard_similarity'][i-1, j-1]:>15.3f}", end='')
            print()
        
        return results
    
    def run_performance_benchmark(self, num_molecules: int = 100) -> Dict:
        """Benchmark performance across molecule sizes"""
        import numpy as np
        print("\n" + "="*80)
        print(f"PERFORMANCE BENCHMARK ({num_molecules} molecules)")
        print("="*80)
        
        from PFASgroups.generate_mol import generate_random_mol, generate_random_carbon_chain, fluorinate_mol
        
        results = {
            'category': 'performance',
            'measurements': [],
            'summary': {
                'total_molecules': 0,
                'mean_time': 0.0,
                'median_time': 0.0,
                'std_time': 0.0,
                'by_size': defaultdict(list),
            }
        }
        
        print(f"\n[PERFORMANCE] Generating and testing {num_molecules} random PFAS molecules...")
        
        execution_times = []
        
        for i in range(num_molecules):
            # Generate random perfluorinated molecule
            try:
                chain_length = int(np.random.randint(4, 12))
                mol = generate_random_carbon_chain(chain_length)
                mol = fluorinate_mol(mol, perfluorinated=True)
                smiles = Chem.MolToSmiles(mol)
                
                # Test performance
                test_result = self.test_single_molecule(smiles)
                
                if test_result['valid']:
                    results['summary']['total_molecules'] += 1
                    exec_time = test_result['execution_time']
                    execution_times.append(exec_time)
                    
                    num_atoms = test_result['molecular_properties']['num_atoms']
                    results['summary']['by_size'][num_atoms].append(exec_time)
                    
                    results['measurements'].append({
                        'smiles': smiles,
                        'num_atoms': num_atoms,
                        'execution_time': exec_time,
                        'num_definitions_matched': len(test_result['detected_ids']),
                    })
                    
                    if (i + 1) % 20 == 0:
                        print(f"   Processed {i+1}/{num_molecules} molecules...")
            except Exception as e:
                if i == 0:
                    print(f"   [WARN] Error generating molecules: {type(e).__name__}: {e}")
                    print(f"   [INFO] Skipping performance benchmark - generate_mol functions may not be available")
                    break
                continue
        
        # Calculate statistics
        if execution_times:
            results['summary']['mean_time'] = np.mean(execution_times)
            results['summary']['median_time'] = np.median(execution_times)
            results['summary']['std_time'] = np.std(execution_times)
            results['summary']['min_time'] = np.min(execution_times)
            results['summary']['max_time'] = np.max(execution_times)
            
            print(f"\n{'='*80}")
            print(f"PERFORMANCE SUMMARY:")
            print(f"  - Molecules tested: {results['summary']['total_molecules']}")
            print(f"  - Mean execution time: {results['summary']['mean_time']*1000:.2f} ms")
            print(f"  - Median execution time: {results['summary']['median_time']*1000:.2f} ms")
            print(f"  - Std deviation: {results['summary']['std_time']*1000:.2f} ms")
            print(f"  - Min/Max: {results['summary']['min_time']*1000:.2f} / {results['summary']['max_time']*1000:.2f} ms")
            
            # Show performance by size
            print(f"\n  - Performance by molecule size:")
            for size in sorted(results['summary']['by_size'].keys()):
                times = results['summary']['by_size'][size]
                mean_time = np.mean(times)
                print(f"    - {size} atoms: {mean_time*1000:.2f} ms (n={len(times)})")
        else:
            print(f"\n{'='*80}")
            print(f"PERFORMANCE SUMMARY:")
            print(f"  [WARN] No valid molecules generated for performance testing")
            print(f"  - Molecules tested: {results['summary']['total_molecules']}")
        
        return results
    
    def export_false_classifications(self, benchmark_results: Dict, timestamp: str) -> str:
        """Export false positives and negatives to a CSV file for detailed analysis"""
        import csv
        from pathlib import Path
        
        # Create output directory
        output_dir = Path(__file__).parent / 'data'
        output_dir.mkdir(exist_ok=True)
        output_file = output_dir / f'false_classifications_{timestamp}.csv'
        
        false_classifications = []
        
        # Collect false negatives from true positives benchmark
        true_pos_results = benchmark_results['benchmarks'].get('true_positives', {})
        for fn_detail in true_pos_results.get('false_negative_details', []):
            false_classifications.append({
                'classification_type': 'False Negative',
                'smiles': fn_detail['smiles'],
                'molecule_name': fn_detail['molecule_name'],
                'category': fn_detail['category'],
                'expected_behavior': 'Should be detected (true PFAS)',
                'actual_behavior': 'Not detected by any definition',
                'definitions_involved': 'None',
                'issue_description': 'Known PFAS compound not detected'
            })
        
        # Collect false positives from true negatives benchmark
        true_neg_results = benchmark_results['benchmarks'].get('true_negatives', {})
        for fp_detail in true_neg_results.get('false_positive_details', []):
            false_classifications.append({
                'classification_type': 'False Positive',
                'smiles': fp_detail['smiles'],
                'molecule_name': fp_detail['molecule_name'],
                'category': fp_detail['category'],
                'expected_behavior': 'Should NOT be detected (non-PFAS)',
                'actual_behavior': 'Incorrectly detected',
                'definitions_involved': ', '.join(fp_detail['detected_by']),
                'issue_description': 'Non-PFAS compound incorrectly classified as PFAS'
            })
        
        # Collect definition-specific failures
        def_specific_results = benchmark_results['benchmarks'].get('definition_specific', {})
        for def_name, tests in def_specific_results.get('by_definition', {}).items():
            for test in tests:
                if test.get('test_passed') == False:
                    expected = test.get('expected_match', 'Unknown')
                    actual = test.get('actual_match', 'Unknown')
                    
                    classification_type = 'Definition-Specific Error'
                    if expected and not actual:
                        classification_type = 'Definition False Negative'
                        issue_desc = f'Expected {def_name} to detect but it did not'
                    elif not expected and actual:
                        classification_type = 'Definition False Positive'  
                        issue_desc = f'{def_name} detected when it should not have'
                    else:
                        issue_desc = f'{def_name}: Expected {expected}, got {actual}'
                    
                    false_classifications.append({
                        'classification_type': classification_type,
                        'smiles': test.get('smiles', 'Unknown'),
                        'molecule_name': test.get('name', 'Definition-specific test'),
                        'category': f'{def_name}_specific_test',
                        'expected_behavior': f'Should match: {expected}',
                        'actual_behavior': f'Actually matched: {actual}',
                        'definitions_involved': def_name,
                        'issue_description': issue_desc
                    })
        
        # Write to CSV
        if false_classifications:
            with open(output_file, 'w', newline='', encoding='utf-8') as f:
                writer = csv.DictWriter(f, fieldnames=[
                    'classification_type', 'smiles', 'molecule_name', 'category',
                    'expected_behavior', 'actual_behavior', 'definitions_involved', 'issue_description'
                ])
                writer.writeheader()
                writer.writerows(false_classifications)
            
            print(f"\n[EXPORT] False classifications saved to: {output_file}")
            print(f"   - Total issues found: {len(false_classifications)}")
            
            # Summary by type
            by_type = {}
            for item in false_classifications:
                cls_type = item['classification_type']
                by_type[cls_type] = by_type.get(cls_type, 0) + 1
            
            for cls_type, count in by_type.items():
                print(f"   - {cls_type}: {count}")
        else:
            print(f"\n[EXPORT] No false classifications found!")
            # Still create empty file
            with open(output_file, 'w', newline='', encoding='utf-8') as f:
                writer = csv.DictWriter(f, fieldnames=[
                    'classification_type', 'smiles', 'molecule_name', 'category',
                    'expected_behavior', 'actual_behavior', 'definitions_involved', 'issue_description'
                ])
                writer.writeheader()
        
        return str(output_file)
    
    def run_complete_benchmark(self) -> Dict:
        """Run all benchmarks"""
        print("\n" + "="*80)
        print("COMPREHENSIVE PFAS DEFINITIONS BENCHMARK")
        print("="*80)
        print(f"\nTesting {len(self.definitions)} PFAS definitions:")
        for def_id, name in self.definition_names.items():
            print(f"  {def_id}. {name}")
        
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        
        # Run all benchmark suites
        results = {
            'metadata': {
                'timestamp': timestamp,
                'definitions_tested': list(self.definition_names.values()),
            },
            'benchmarks': {}
        }
        
        # 1. True positives
        results['benchmarks']['true_positives'] = self.run_true_positives_benchmark()
        
        # 2. True negatives
        results['benchmarks']['true_negatives'] = self.run_true_negatives_benchmark()
        
        # 3. Edge cases
        results['benchmarks']['edge_cases'] = self.run_edge_cases_benchmark()
        
        # 4. Definition-specific tests
        results['benchmarks']['definition_specific'] = self.run_definition_specific_tests()
        
        # 5. Concordance analysis
        results['benchmarks']['concordance'] = self.run_concordance_analysis()
        
        # 6. Performance
        results['benchmarks']['performance'] = self.run_performance_benchmark(num_molecules=100)
        
        # Convert frozensets and defaultdicts for JSON serialization
        def convert_for_json(obj):
            """Recursively convert non-JSON-serializable objects to JSON-compatible types"""
            if isinstance(obj, defaultdict):
                # Convert defaultdict to regular dict
                return convert_for_json(dict(obj))
            elif isinstance(obj, dict):
                new_dict = {}
                for key, value in obj.items():
                    if isinstance(key, frozenset):
                        # Convert frozenset to sorted tuple then to string
                        new_key = str(tuple(sorted(key)))
                    else:
                        new_key = key
                    new_dict[new_key] = convert_for_json(value)
                return new_dict
            elif isinstance(obj, list):
                return [convert_for_json(item) for item in obj]
            elif isinstance(obj, (set, frozenset)):
                return list(obj)
            else:
                return obj
        
        results = convert_for_json(results)
        
        # Save results
        output_dir = os.path.join(script_dir, 'data')
        os.makedirs(output_dir, exist_ok=True)
        output_file = os.path.join(output_dir, f'pfas_definitions_benchmark_{timestamp}.json')
        
        with open(output_file, 'w') as f:
            json.dump(results, f, indent=2, default=str)
        
        # Export false classifications to CSV
        false_classifications_file = self.export_false_classifications(results, timestamp)
        
        print(f"\n{'='*80}")
        print(f"[SUCCESS] Benchmark complete! Results saved to:")
        print(f"   - Main results: {output_file}")
        print(f"   - False classifications: {false_classifications_file}")
        
        return results, output_file


def main():
    """Main function to run benchmark"""
    if '--export-test-compounds' in sys.argv:
        benchmark = PFASDefinitionBenchmark(load_definitions=False)
        output_path = benchmark.export_test_compounds_csv(source='default')
        print(f"[EXPORT] Test compounds CSV saved to: {output_path}")
        return None, output_path
    benchmark = PFASDefinitionBenchmark()
    results, output_file = benchmark.run_complete_benchmark()
    
    return results, output_file


if __name__ == "__main__":
    main()
