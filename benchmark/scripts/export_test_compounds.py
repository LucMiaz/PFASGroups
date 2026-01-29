"""
Export all test compounds from benchmark to CSV file.
Creates a CSV with columns: smiles, definition_name, reason, expected_result
"""

import csv
from pathlib import Path

# All test compounds from benchmark_pfas_definitions.py
test_data = []

# 1. TRUE POSITIVES - expected to be detected by MOST definitions
# These are added as "All definitions (general)" since they should generally match
true_positives = [
    # perfluoroalkyl_carboxylic_acids
    ('OC(=O)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)F', 'All definitions (general)', 'PFOA (C8) - perfluoroalkyl carboxylic acid', 'True'),
    ('OC(=O)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)F', 'All definitions (general)', 'PFHxA (C6) - perfluoroalkyl carboxylic acid', 'True'),
    ('OC(=O)C(F)(F)C(F)(F)C(F)(F)F', 'All definitions (general)', 'PFBA (C4) - perfluoroalkyl carboxylic acid', 'True'),
    ('OC(=O)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)F', 'All definitions (general)', 'PFDA (C10) - perfluoroalkyl carboxylic acid', 'True'),
    ('OC(=O)C(F)(F)F', 'All definitions (general)', 'TFA (C2) - perfluoroalkyl carboxylic acid', 'True'),
    
    # perfluoroalkyl_sulfonic_acids
    ('FC(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)S(=O)(=O)O', 'All definitions (general)', 'PFOS (C8) - perfluoroalkyl sulfonic acid', 'True'),
    ('FC(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)S(=O)(=O)O', 'All definitions (general)', 'PFHxS (C6) - perfluoroalkyl sulfonic acid', 'True'),
    ('FC(F)(F)C(F)(F)C(F)(F)C(F)(F)S(=O)(=O)O', 'All definitions (general)', 'PFBS (C4) - perfluoroalkyl sulfonic acid', 'True'),
    ('FC(F)(F)C(F)(F)S(=O)(=O)O', 'All definitions (general)', 'PFEtS (C2) - perfluoroalkyl sulfonic acid', 'True'),
    
    # fluorotelomer_alcohols
    ('OCCC(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)F', 'All definitions (general)', '6:2 FTOH - fluorotelomer alcohol', 'True'),
    ('OCCC(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)F', 'All definitions (general)', '8:2 FTOH - fluorotelomer alcohol', 'True'),
    ('OCCC(F)(F)C(F)(F)C(F)(F)F', 'All definitions (general)', '4:2 FTOH - fluorotelomer alcohol', 'True'),
    ('OCCC(F)(F)C(F)(F)C(F)(F)C(F)(F)F', 'All definitions (general)', '5:2 FTOH - fluorotelomer alcohol', 'True'),
    
    # perfluoroalkyl_sulfonamides
    ('FC(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)S(=O)(=O)N', 'All definitions (general)', 'PFOSA - perfluoroalkyl sulfonamide', 'True'),
    ('FC(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)S(=O)(=O)NC', 'All definitions (general)', 'N-MeFOSA - perfluoroalkyl sulfonamide', 'True'),
    ('FC(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)S(=O)(=O)NCC', 'All definitions (general)', 'N-EtFOSA - perfluoroalkyl sulfonamide', 'True'),
    
    # perfluoroalkyl_iodides
    ('FC(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)I', 'All definitions (general)', 'C8F17I - perfluoroalkyl iodide', 'True'),
    ('FC(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)I', 'All definitions (general)', 'C6F13I - perfluoroalkyl iodide', 'True'),
    
    # side_chain_fluorinated
    ('FC(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)CCOC(=O)C', 'All definitions (general)', '8:2 fluorotelomer acrylate - side chain fluorinated', 'True'),
    ('FC(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)CCOP(=O)(O)O', 'All definitions (general)', '6:2 fluorotelomer phosphate - side chain fluorinated', 'True'),
    
    # cyclic_perfluorinated
    ('FC1(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C1(F)F', 'All definitions (general)', 'Perfluorocyclohexane - cyclic perfluorinated', 'True'),
    ('FC1(F)C(F)(F)C(F)(F)C(F)(F)C1(F)F', 'All definitions (general)', 'Perfluorocyclopentane - cyclic perfluorinated', 'True'),
    
    # perfluoroethers
    ('FC(F)(F)C(F)(F)OC(F)(F)C(F)(F)F', 'All definitions (general)', 'Perfluorodiethyl ether - perfluoroether', 'True'),
    ('FC(F)(F)OC(F)(F)C(F)(F)OC(F)(F)F', 'All definitions (general)', 'PFPE - perfluoroether', 'True'),
]

# 2. TRUE NEGATIVES - should NOT be detected by ANY definition
true_negatives = [
    # non_fluorinated_analogs
    ('OC(=O)CCCCCCC', 'All definitions (general)', 'Octanoic acid - non-fluorinated analog of PFOA', 'False'),
    ('CCCCCCCCS(=O)(=O)O', 'All definitions (general)', 'Octane sulfonic acid - non-fluorinated analog of PFOS', 'False'),
    ('OCCCCCCCC', 'All definitions (general)', 'Octanol - non-fluorinated analog of FTOH', 'False'),
    
    # monofluorinated
    ('OC(=O)CF', 'All definitions (general)', 'Fluoroacetic acid - monofluorinated', 'False'),
    ('FC', 'All definitions (general)', 'Fluoromethane - monofluorinated', 'False'),
    ('FCC', 'All definitions (general)', 'Fluoroethane - monofluorinated', 'False'),
    
    # difluorinated_non_geminal
    ('FCCCF', 'All definitions (general)', '1,4-difluorobutane - difluorinated non-geminal', 'False'),
    ('FC(C)CF', 'All definitions (general)', '1,3-difluoropropane - difluorinated non-geminal', 'False'),
    
    # inorganic_fluorides
    ('F[Na]', 'All definitions (general)', 'Sodium fluoride - inorganic', 'False'),
    ('F[Ca]F', 'All definitions (general)', 'Calcium fluoride - inorganic', 'False'),
    
    # single_CF_group
    ('CC(F)(F)CC', 'All definitions (general)', '2,2-difluorobutane - single CF2', 'False'),
    ('CC(C)(F)C', 'All definitions (general)', '2-fluoroisobutane - single CF', 'False'),
    
    # pharmaceutical_fluorinated
    ('FC(F)(F)c1ccc(cc1)C(=O)O', 'All definitions (general)', 'Trifluoroacetic acid derivative - aromatic pharmaceutical', 'False'),
    ('c1ccc(F)cc1', 'All definitions (general)', 'Fluorobenzene - aromatic pharmaceutical', 'False'),
    
    # EU_specific_exclusions
    ('C(F)(F)(F)Cl', 'All definitions (general)', 'Trifluoromethyl chloride (CF3-Cl) - halogen bonded to CF3', 'False'),
    ('C(F)(F)(F)Br', 'All definitions (general)', 'Trifluoromethyl bromide (CF3-Br) - halogen bonded to CF3', 'False'),
    ('C(F)(F)(F)I', 'All definitions (general)', 'Trifluoromethyl iodide (CF3-I) - halogen bonded to CF3', 'False'),
    ('C(F)(F)(F)[H]', 'All definitions (general)', 'Fluoroform (CHF3) - H bonded to CF3', 'False'),
    ('C(F)(F)(F)O', 'All definitions (general)', 'Trifluoromethanol (CF3-OH) - alcohol bonded to CF3', 'False'),
    ('C(F)(F)(F)OC', 'All definitions (general)', 'Trifluoromethyl methyl ether (CF3-OCH3) - ether bonded to CF3', 'False'),
    ('C(F)(F)(F)N', 'All definitions (general)', 'Trifluoromethylamine (CF3-NH2) - amine bonded to CF3', 'False'),
    ('C(F)(F)(Cl)C', 'All definitions (general)', '1-Chloro-1,1-difluoroethane (CF2Cl-CH3) - halogen bonded to CF2', 'False'),
    ('C(F)(F)([H])C', 'All definitions (general)', '1,1-Difluoroethane (CF2H-CH3) - H bonded to CF2', 'False'),
    ('C(F)(F)(OC)C', 'All definitions (general)', '1,1-Difluoro-1-methoxyethane (CF2-OCH3-CH3) - ether bonded to CF2', 'False'),
    ('C(F)(F)(N)C', 'All definitions (general)', '1,1-Difluoro-1-aminoethane (CF2-NH2-CH3) - amine bonded to CF2', 'False'),
    ('FC1(F)Oc2ccccc2O1', 'All definitions (general)', 'Benzodioxole with CF2 - aromatic ether exclusion', 'False'),
]

# 3. EDGE CASES - borderline compounds where definitions may disagree
edge_cases = [
    # short_chain_perfluorinated
    ('FC(F)(F)C(F)(F)F', 'Edge case', 'Perfluoropropane (C3) - short chain', 'Variable'),
    ('OC(=O)C(F)(F)F', 'Edge case', 'TFA - C2 PFCA - short chain', 'Variable'),
    ('FC(F)(F)S(=O)(=O)O', 'Edge case', 'C1 perfluorosulfonic acid - ultra short', 'Variable'),
    
    # polyfluorinated_not_perfluorinated
    ('OC(=O)CC(F)(F)C(F)(F)F', 'Edge case', 'H-PFBA - one non-fluorinated carbon', 'Variable'),
    ('FC(F)C(F)(F)C(F)(F)F', 'Edge case', 'Polyfluorobutane - missing one F', 'Variable'),
    ('OC(=O)C(F)(F)C(F)(F)CC', 'Edge case', 'Semi-fluorinated carboxylic acid', 'Variable'),
    
    # fluorine_rich_aromatics
    ('FC(F)(F)c1c(F)c(F)c(F)c(F)c1F', 'Edge case', 'Perfluorotoluene - aromatic', 'Variable'),
    ('c1c(F)c(F)c(F)c(F)c1F', 'Edge case', 'Pentafluorobenzene - aromatic', 'Variable'),
    
    # heteroatom_connected
    ('FC(F)(F)C(F)(F)OC', 'Edge case', 'Perfluoroethyl methyl ether - OECD vs EU', 'Variable'),
    ('FC(F)(F)C(F)(F)NC', 'Edge case', 'Perfluoroethylmethylamine - OECD vs EU', 'Variable'),
    
    # halogen_substituted
    ('ClC(F)(F)C(F)(F)F', 'Edge case', 'Chloroperfluoropropane - halogen exclusions', 'Variable'),
    ('FC(F)(F)C(F)(F)C(Cl)(F)F', 'Edge case', 'Chlorine on perfluoroalkyl - distance matters', 'Variable'),
    
    # ultra_short_chain
    ('FC(F)(F)F', 'Edge case', 'CF4 - carbon tetrafluoride - ultra short', 'Variable'),
    ('FC(F)(F)C(F)(F)F', 'Edge case', 'C2F6 - hexafluoroethane - ultra short', 'Variable'),
    
    # EU_specific_edge_cases
    ('FC(F)(F)C', 'Edge case', 'Trifluoromethyl-methane - short chain with CF3', 'Variable'),
    ('FC(F)(OC(=O)C)C', 'Edge case', 'CF2 with acetoxy - EU should reject, others may accept', 'Variable'),
    ('FC(F)(SC)C', 'Edge case', 'CF2 with methylthio - EU exclusion for -SR', 'Variable'),
    ('FC(F)(F)C(F)(F)SC', 'Edge case', 'CF3-CF2-SCH3 - sulfur at distance', 'Variable'),
    ('FC(F)=C(F)OC(F)(F)F', 'Edge case', 'Fluorovinyl ether with CF3 - unsaturated', 'Variable'),
    ('FC(F)(N(C)C)C(F)(F)F', 'Edge case', 'CF2 with tertiary amine - EU should reject', 'Variable'),
    ('FC(F)(F)C(N(C)C)C(F)(F)F', 'Edge case', 'Amine between CF3 groups - EU may accept', 'Variable'),
]

# 4. DEFINITION-SPECIFIC TESTS
definition_specific = [
    # OECD_specific
    ('FC(F)(F)Cl', 'OECD', 'CF3Cl - CF3 bonded to Cl should NOT match OECD', 'False'),
    ('FC(F)(F)C(F)(F)F', 'OECD', 'C2F6 - should match OECD', 'True'),
    ('FC(F)C(F)(F)F', 'OECD', 'CF2H-CF3 - CF2 bonded to H may not match', 'False'),
    
    # EU_specific - CF2 exclusions
    ('FC(F)(Cl)C', 'EU', 'CF2 bonded to Cl - should NOT match EU', 'False'),
    ('FC(F)(Br)C', 'EU', 'CF2 bonded to Br - should NOT match EU', 'False'),
    ('FC(F)(I)C', 'EU', 'CF2 bonded to I - should NOT match EU', 'False'),
    ('FC(F)([H])C', 'EU', 'CF2 bonded to H - should NOT match EU', 'False'),
    ('FC(F)(O)C', 'EU', 'CF2 bonded to -OH - should NOT match EU', 'False'),
    ('FC(F)(OC)C', 'EU', 'CF2 bonded to -OCH3 - should NOT match EU', 'False'),
    ('FC(F)(OCC)C', 'EU', 'CF2 bonded to -OEt - should NOT match EU', 'False'),
    ('FC(F)(Oc1ccccc1)C', 'EU', 'CF2 bonded to -OPh - should NOT match EU', 'False'),
    ('FC(F)(OC(=O)C)C', 'EU', 'CF2 bonded to -OCOCH3 - should NOT match EU', 'False'),
    ('FC(F)(N)C', 'EU', 'CF2 bonded to -NH2 - should NOT match EU', 'False'),
    ('FC(F)(NC)C', 'EU', 'CF2 bonded to -NHCH3 - should NOT match EU', 'False'),
    ('FC(F)(N(C)C)C', 'EU', 'CF2 bonded to -N(CH3)2 - should NOT match EU', 'False'),
    
    # EU_specific - CF2 acceptable
    ('FC(F)(C(F)(F)F)C', 'EU', 'CF2 with perfluoro chain - should match EU', 'True'),
    ('FC(F)(P)C', 'EU', 'CF2 bonded to P - should match EU', 'True'),
    ('FC(F)(F)C(F)(F)C(F)(P)OC', 'EU', 'CF2 chain with P and -OR at distance - should match EU', 'True'),
    
    # EU_specific - CF3 exclusions
    ('FC(F)(F)Cl', 'EU', 'CF3 bonded to Cl - should NOT match EU', 'False'),
    ('FC(F)(F)Br', 'EU', 'CF3 bonded to Br - should NOT match EU', 'False'),
    ('FC(F)(F)I', 'EU', 'CF3 bonded to I - should NOT match EU', 'False'),
    ('FC(F)(F)[H]', 'EU', 'CF3 bonded to H (CHF3) - should NOT match EU', 'False'),
    ('FC(F)(F)O', 'EU', 'CF3 bonded to -OH - should NOT match EU', 'False'),
    ('FC(F)(F)OC', 'EU', 'CF3 bonded to -OCH3 - should NOT match EU', 'False'),
    ('FC(F)(F)OCC', 'EU', 'CF3 bonded to -OEt - should NOT match EU', 'False'),
    ('FC(F)(F)Oc1ccccc1', 'EU', 'CF3 bonded to -OPh - should NOT match EU', 'False'),
    ('FC(F)(F)OC(=O)C', 'EU', 'CF3 bonded to -OCOCH3 - should NOT match EU', 'False'),
    ('FC(F)(F)N', 'EU', 'CF3 bonded to -NH2 - should NOT match EU', 'False'),
    ('FC(F)(F)NC', 'EU', 'CF3 bonded to -NHCH3 - should NOT match EU', 'False'),
    ('FC(F)(F)N(C)C', 'EU', 'CF3 bonded to -N(CH3)2 - should NOT match EU', 'False'),
    
    # EU_specific - CF3 acceptable
    ('FC(F)(F)F', 'EU', 'CF4 - should match EU (no exclusions)', 'True'),
    ('FC(F)(F)C(F)(F)F', 'EU', 'Perfluoroethane - should match EU', 'True'),
    ('FC(F)(F)C(F)(F)O', 'EU', 'CF3-CF2-OH (ether at distance) - should match EU', 'True'),
    ('FC(F)(F)C(F)(F)N', 'EU', 'CF3-CF2-NH2 (amine at distance) - should match EU', 'True'),
    
    # EU_specific - complex
    ('FC(F)=C(F)OC(F)(F)OC(F)=C(F)F', 'EU', 'PFPE with ether - should match (ether not directly on CF2/CF3)', 'True'),
    ('FC1(F)Oc2cccc(C=O)c2O1', 'EU', 'Benzodioxole with CF2 bonded to aromatic O - should NOT match EU', 'False'),
    ('[H]OC([H])([H])C([H])(O[H])C([H])([H])N1C2=C([H])C(F)=C(N([H])C(=O)C3(C4=C([H])C5=C(OC(F)(F)O5)C([H])=C4[H])C([H])([H])C3([H])[H])C([H])=C2C([H])=C1C(C([H])([H])[H])(C([H])([H])[H])C([H])([H])O[H]', 'EU', 'Complex pharma with CF2 in benzodioxole - should NOT match', 'False'),
    
    # OPPT_specific
    ('FC(F)(F)C(F)(F)OC(F)(F)C(F)(F)F', 'OPPT', 'PFPE - matches OPPT criterion 2', 'True'),
    ('FC(F)(F)C(C(F)(F)F)(C(F)(F)F)C(F)(F)F', 'OPPT', 'Branched perfluoro - OPPT criterion 3', 'True'),
    
    # UK_specific
    ('FC(F)(F)C(F)(F)F', 'UK', 'Should match UK definition', 'True'),
    
    # PFASTRUCTv5_specific - threshold tests
    ('FC(F)(F)C(F)(F)CCCCCCCC', 'PFASTRUCTv5', 'Lower F ratio - test threshold', 'True'),
    ('FC(F)(F)CCCCCCCCCCCC', 'PFASTRUCTv5', 'Very low F ratio (~0.15)', 'False'),
    
    # PFASTRUCTv5_specific - validation set TRUE POSITIVES
    ('C(=O)(C(F)(F)Cl)C(F)(F)Cl', 'PFASTRUCTv5', '1,3-Dichloro-1,1,3,3-tetrafluoropropan-2-one (DTXSID1073159)', 'True'),
    ('COC(C(OC(C(=O)C(F)(F)F)(F)F)F)(F)F', 'PFASTRUCTv5', '1,1,1,3,3-Pentafluoro-3-(1,2,2-trifluoro-2-methoxyethoxy)propan-2-one (DTXSID40748237)', 'True'),
    ('C(C(=O)CC(C(=O)C(F)(F)F)(F)F)C(C(=O)C(F)(F)F)(F)F', 'PFASTRUCTv5', '2,5,8-Nonanetrione, decafluoro- (DTXSID60435490)', 'True'),
    ('C(C(CC(F)(F)F)(F)F)C(CC(F)(F)F)(F)F', 'PFASTRUCTv5', '1,1,1,3,3,5,5,7,7,7-Decafluoroheptane (DTXSID20776331)', 'True'),
    ('C(=C(F)F)(OC(OC(=C(F)F)F)(F)F)F', 'PFASTRUCTv5', '1-(Difluoro[(trifluoroethenyl)oxy]methoxy)-1,2,2-trifluoroethene (DTXSID00896607)', 'True'),
    ('C(=C(F)F)(C(C(=C(F)F)F)(F)F)F', 'PFASTRUCTv5', '1,1,2,3,3,4,5,5-Octafluoropenta-1,4-diene (DTXSID30529779)', 'True'),
    ('C1(=C(C(=C(C(=C1F)F)F)F)F)C2=C(C(C(=C(C2(F)F)F)F)(F)F)F', 'PFASTRUCTv5', 'Dodecafluoro-2,5-dihydro-1,1-biphenyl (DTXSID20525346)', 'True'),
    ('C(C(=O)F)(F)F', 'PFASTRUCTv5', '2,2-Difluoroacetyl fluoride (DTXSID70381191)', 'True'),
    ('C(=C(C(F)(F)Cl)F)(C(F)(F)Cl)F', 'PFASTRUCTv5', '(2E)-1,4-Dichloro-hexafluoro-2-butene (DTXSID601020902)', 'True'),
    ('C(C(C(C(F)F)F)F)F', 'PFASTRUCTv5', '1,1,2,3,4-Pentafluorobutane (DTXSID50931113)', 'True'),
    ('CCN(CC)C1(OCCCO1)C(C(F)(F)F)F', 'PFASTRUCTv5', 'N,N-Diethyl-2-(1,2,2,2-tetrafluoroethyl)-1,3-dioxan-2-amine (DTXSID40519722)', 'True'),
    ('COC(CC(=C(F)F)C(F)(F)F)SC1=CC=CC=C1', 'PFASTRUCTv5', '[4,4-Difluoro-1-methoxy-3-(trifluoromethyl)but-3-en-1-yl]sulfanyl)benzene (DTXSID20809173)', 'True'),
    ('C1=CC(=CC(=C1)Cl)NOC(C(C(=O)C(F)(F)F)F)(F)F', 'PFASTRUCTv5', '4-[(3-Chloroanilino)oxy]-1,1,1,3,4,4-hexafluorobutan-2-one (DTXSID201023166)', 'True'),
    ('CC(C)(CS(=O)(=O)O)NC(=O)CCSCCC(C(C(C(F)(F)F)(F)F)(F)F)(F)F', 'PFASTRUCTv5', '4:2 Fluorotelomer thioether amido sulfonic acid (DTXSID00892528)', 'True'),
    ('CN(C)CCCN(CCC(=O)O)S(=O)(=O)C(C(C(C(F)(F)F)(F)F)(F)F)(F)F', 'PFASTRUCTv5', 'N-(Perfluorobutanesulfonyl)-N-(3-dimethylaminopropyl)-3-aminopropanoic acid (DTXSID20882022)', 'True'),
    ('C1=CC(=C(C=C1N)Cl)OC(C(OC(F)(F)F)F)(F)F', 'PFASTRUCTv5', '3-Chloro-4-[1,1,2-trifluoro-2-(trifluoromethoxy)ethoxy]aniline (DTXSID10576660)', 'True'),
    ('C1=CC=C2C(=C1)C(OC(O2)(C(F)(F)F)F)(F)F', 'PFASTRUCTv5', '2,4,4-Trifluoro-2-(trifluoromethyl)-2H,4H-1,3-benzodioxine (DTXSID80661330)', 'True'),
    ('CCCOC(=O)C(=C(F)F)C(F)(F)F', 'PFASTRUCTv5', 'Propyl 3,3-difluoro-2-(trifluoromethyl)prop-2-enoate (DTXSID80804225)', 'True'),
    ('C[N+](C)(CCCNS(=O)(=O)CCC(C(C(C(C(C(F)(F)F)(F)F)(F)F)(F)F)(F)F)(F)F)CC(=O)[O-]', 'PFASTRUCTv5', '6:2 Fluorotelomer sulfonamide betaine (DTXSID4041284)', 'True'),
    
    # PFASTRUCTv5_specific - validation set TRUE NEGATIVES
    ('CCCCC(CC(C(F)(F)Br)(F)Cl)Br', 'PFASTRUCTv5', '1,4-dibromo-2-chloro-1,1,2-trifluorooctane - low F ratio (DTXSID10382122)', 'False'),
    ('C(C(OC(Cl)(Cl)Cl)(F)F)(F)Cl', 'PFASTRUCTv5', '2-Chloro-1,1,2-trifluoro-1-(trichloromethoxy)ethane - low F ratio (DTXSID30963477)', 'False'),
    ('C1=COC(=C1)C=C(S(=O)(=O)C(F)(F)F)S(=O)(=O)C(F)(F)F', 'PFASTRUCTv5', '2-[2,2-Bis(trifluoromethanesulfonyl)ethenyl]furan - borderline (DTXSID10711579)', 'False'),
    ('COC(=C(F)Cl)F', 'PFASTRUCTv5', '1-Chloro-1,2-difluoro-2-methoxyethene - too few fluorines (DTXSID50777992)', 'False'),
]

# Combine all data
all_data = true_positives + true_negatives + edge_cases + definition_specific

# Write to CSV
output_file = Path(__file__).parent / 'benchmark_test_compounds.csv'

with open(output_file, 'w', newline='', encoding='utf-8') as f:
    writer = csv.writer(f)
    writer.writerow(['smiles', 'definition_name', 'reason', 'expected_result'])
    writer.writerows(all_data)

print(f"Exported {len(all_data)} test compounds to: {output_file}")
print(f"  - True positives: {len(true_positives)}")
print(f"  - True negatives: {len(true_negatives)}")
print(f"  - Edge cases: {len(edge_cases)}")
print(f"  - Definition-specific: {len(definition_specific)}")
