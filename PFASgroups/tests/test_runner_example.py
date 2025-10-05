"""
Example pytest configuration and usage for the enhanced test_examples.py

This script demonstrates how to run the tests and generate summary reports.
"""

import subprocess
import sys
import os

def run_pytest_with_summary():
    """Run pytest and generate summary files."""
    
    print("Running PFAS Groups tests with pytest...")
    print("=" * 50)
    
    # Change to the tests directory
    tests_dir = os.path.dirname(__file__)
    print(f"Tests directory: {tests_dir}")
    
    try:
        # Run pytest with verbose output
        result = subprocess.run([
            sys.executable, '-m', 'pytest', 
            'test_examples.py', 
            '-v',  # verbose
            '--tb=short',  # shorter traceback format
            '-s'  # don't capture output (so we can see prints)
        ], 
        cwd=tests_dir,
        capture_output=False,  # Show output in real time
        text=True
        )
        
        print(f"\nPytest finished with exit code: {result.returncode}")
        
        # Check for generated files
        summary_file = os.path.join(tests_dir, 'test_summary_report.json')
        if os.path.exists(summary_file):
            print(f"✓ Test summary generated: {summary_file}")
        else:
            print("✗ Test summary not found")
            
        # Check for CSV files
        csv_files = ['oecd_test_results.csv', 'generic_test_results.csv', 'specificity_test_results.csv']
        for csv_file in csv_files:
            full_path = os.path.join(tests_dir, csv_file)
            if os.path.exists(full_path):
                print(f"✓ Detailed results: {csv_file}")
            else:
                print(f"- Results file not found: {csv_file}")
        
        return result.returncode == 0
        
    except FileNotFoundError:
        print("Error: pytest not found. Please install pytest:")
        print("pip install pytest")
        return False
    except Exception as e:
        print(f"Error running pytest: {e}")
        return False

def run_manual_tests():
    """Run tests manually without pytest."""
    print("Running PFAS Groups tests manually...")
    print("=" * 50)
    
    try:
        # Import and run tests manually
        sys.path.append(os.path.join(os.path.dirname(__file__), '..'))
        from tests.test_examples import TestPFASGroups
        
        tester = TestPFASGroups()
        result = tester.run_all_tests()
        
        print(f"\nManual tests finished: {'PASSED' if result else 'FAILED'}")
        return result
        
    except Exception as e:
        print(f"Error running manual tests: {e}")
        return False

def display_summary(summary_file='tests/test_summary_report.json'):
    """Display the test summary in a readable format."""
    import json
    
    try:
        with open(summary_file, 'r') as f:
            summary = json.load(f)
        
        print("\n" + "=" * 60)
        print("PFAS GROUPS TEST SUMMARY")
        print("=" * 60)
        
        metadata = summary.get('test_metadata', {})
        print(f"Test Date: {metadata.get('timestamp', 'Unknown')}")
        print(f"Duration: {metadata.get('test_duration_seconds', 0):.1f} seconds")
        print(f"Python Version: {metadata.get('python_version', 'Unknown')}")
        
        overall = summary.get('overall_summary', {})
        print(f"\nOverall Status: {overall.get('test_status', 'UNKNOWN')}")
        print(f"Overall Accuracy: {overall.get('overall_accuracy', 0):.1%}")
        print(f"Total Tests: {overall.get('total_tests_run', 0)}")
        print(f"Successful: {overall.get('total_successful_detections', 0)}")
        
        # OECD Results
        oecd = summary.get('oecd_test_results', {})
        if oecd:
            print(f"\nOECD Tests:")
            print(f"  Detection Rate: {oecd.get('overall_detection_rate', 0):.1%}")
            print(f"  Tests: {oecd.get('successful_detections', 0)}/{oecd.get('total_tests', 0)}")
        
        # Generic Results
        generic = summary.get('generic_test_results', {})
        if generic:
            print(f"\nGeneric Tests:")
            print(f"  Detection Rate: {generic.get('overall_detection_rate', 0):.1%}")
            print(f"  Tests: {generic.get('successful_detections', 0)}/{generic.get('total_tests', 0)}")
        
        # Specificity Results
        specificity = summary.get('specificity_test_results', {})
        if specificity:
            print(f"\nSpecificity Tests:")
            print(f"  Detection Rate: {specificity.get('detection_rate', 0):.1%}")
            print(f"  Specificity Rate: {specificity.get('specificity_rate', 0):.1%}")
            print(f"  Avg Groups per Test: {specificity.get('average_detected_groups_per_test', 0):.1f}")
        
        print("=" * 60)
        
    except FileNotFoundError:
        print(f"Summary file not found: {summary_file}")
    except Exception as e:
        print(f"Error reading summary: {e}")

if __name__ == "__main__":
    print("PFAS Groups Test Runner")
    print("=" * 30)
    
    if len(sys.argv) > 1:
        if sys.argv[1].lower() == 'pytest':
            success = run_pytest_with_summary()
        elif sys.argv[1].lower() == 'manual':
            success = run_manual_tests()
        elif sys.argv[1].lower() == 'summary':
            display_summary()
            sys.exit(0)
        else:
            print(f"Unknown option: {sys.argv[1]}")
            print("Usage: python run_tests.py [pytest|manual|summary]")
            sys.exit(1)
    else:
        # Try pytest first, fall back to manual
        print("Attempting to run with pytest...")
        success = run_pytest_with_summary()
        
        if not success:
            print("\nPytest failed or not available, trying manual tests...")
            success = run_manual_tests()
    
    # Display summary if available
    if os.path.exists('tests/test_summary_report.json'):
        display_summary()
    
    sys.exit(0 if success else 1)