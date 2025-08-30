#!/usr/bin/env python3
"""
Standalone script to run PFASgroups tests.

This script can be run from the package root directory to test the PFASgroups functionality
without import issues.
"""

import sys
import os

# Add the package directory to Python path
package_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'PFASgroups')
sys.path.insert(0, package_dir)

try:
    # Import the test module
    import test_examples
    
    def main():
        print("PFASgroups Test Runner")
        print("=" * 30)
        
        if len(sys.argv) > 1:
            command = sys.argv[1].lower()
            
            if command == "quick":
                print("Running quick test...")
                success = test_examples.run_quick_test()
                sys.exit(0 if success else 1)
                
            elif command == "full":
                print("Running full test suite...")
                tester = test_examples.TestPFASGroups()
                success = tester.run_all_tests()
                sys.exit(0 if success else 1)
                
            elif command == "pfoa":
                print("Testing PFOA-like compounds...")
                success = test_examples.test_pfoa_like_compounds()
                sys.exit(0 if success else 1)
                
            elif command == "generate":
                n_compounds = int(sys.argv[2]) if len(sys.argv) > 2 else 5
                print(f"Generating test datasets ({n_compounds} compounds per group)...")
                test_examples.generate_oecd_test_compounds('oecd_test_compounds.csv', n_compounds)
                test_examples.generate_generic_test_compounds('generic_test_compounds.csv', n_compounds)
                print("Test datasets generated successfully!")
                
            elif command == "validate":
                oecd_file = sys.argv[2] if len(sys.argv) > 2 else 'oecd_test_compounds.csv'
                generic_file = sys.argv[3] if len(sys.argv) > 3 else 'generic_test_compounds.csv'
                
                if os.path.exists(oecd_file):
                    print(f"Validating {oecd_file}...")
                    test_examples.validate_test_compounds(oecd_file, 'oecd_validated.csv')
                
                if os.path.exists(generic_file):
                    print(f"Validating {generic_file}...")
                    test_examples.validate_test_compounds(generic_file, 'generic_validated.csv')
                    
            else:
                print(f"Unknown command: {command}")
                print_usage()
        else:
            print("Running default quick test...")
            success = test_examples.run_quick_test()
            if success:
                print("\n✓ Quick test completed successfully!")
            else:
                print("\n✗ Quick test failed!")
            print_usage()

    def print_usage():
        print("\nUsage:")
        print("  python run_tests.py [command]")
        print("\nCommands:")
        print("  quick     - Run quick functionality test (default)")
        print("  full      - Run comprehensive test suite")
        print("  pfoa      - Test PFOA-like compound detection")
        print("  generate [n] - Generate test datasets (n compounds per group)")
        print("  validate [oecd_file] [generic_file] - Validate test compounds")

    if __name__ == "__main__":
        main()

except ImportError as e:
    print(f"Error importing test_examples: {e}")
    print("\nTroubleshooting:")
    print("1. Make sure you're running this script from the PFASgroups package root directory")
    print("2. Ensure all required dependencies (RDKit, pandas, numpy, tqdm) are installed")
    print("3. Check that the PFASgroups package structure is correct")
    sys.exit(1)
