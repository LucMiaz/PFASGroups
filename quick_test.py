import sys
import os

sys.path.insert(0, '.')

try:
    from PFASgroups.test_examples import test_pfoa_like_compounds
    print("Testing PFOA-like compounds...")
    result = test_pfoa_like_compounds()
    print(f"Test result: {result}")
    
    # Also test the specificity function on a small subset
    from PFASgroups.test_examples import create_specificity_test_molecules, test_pfas_group_specificity
    import pandas as pd
    
    print("\nTesting specificity detection...")
    test_molecules = create_specificity_test_molecules()
    
    # Filter to just test a few dual-SMARTS groups
    dual_smarts_test = test_molecules[test_molecules['group_ids'].apply(lambda x: any(g in [3, 4, 8, 9] for g in x))]
    
    if len(dual_smarts_test) > 0:
        print(f"Testing {len(dual_smarts_test)} molecules for dual-SMARTS groups...")
        results = test_pfas_group_specificity(dual_smarts_test, output_file=None, verbose=True)
        
        # Summary
        detection_rate = results['expected_group_detected'].mean()
        print(f"Dual-SMARTS detection rate: {detection_rate:.1%}")
    else:
        print("No dual-SMARTS test molecules found")
        
except Exception as e:
    print(f"Error: {e}")
    import traceback
    traceback.print_exc()
