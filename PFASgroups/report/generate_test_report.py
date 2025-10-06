"""
PFAS Algorithm Test Report Generator

This is the main entry point for generating PFAS algorithm test reports.
It automatically tries to run fresh tests, but falls back to analyzing
existing data if there are dependency issues.
"""

import sys
import os
import subprocess
from datetime import datetime

def run_with_fallback():
    """
    Run automated testing with fallback to existing data analysis.
    """
    print("🚀 PFAS Algorithm Test Report Generator")
    print("=" * 50)
    print("Attempting to run fresh tests and generate comprehensive report...")
    print()
    
    # First, try to run the full automated test runner
    try:
        print("📋 Attempting to run fresh tests...")
        result = subprocess.run([sys.executable, 'automated_test_runner.py'], 
                              capture_output=True, text=True, timeout=300)
        
        if result.returncode == 0:
            print("✅ Fresh tests completed successfully!")
            print("\n" + result.stdout)
            return True
        else:
            print(f"⚠️ Fresh tests failed with return code {result.returncode}")
            print("Error output:", result.stderr)
            print("\nFalling back to existing data analysis...")
            
    except FileNotFoundError:
        print("⚠️ automated_test_runner.py not found")
        print("Falling back to existing data analysis...")
    except subprocess.TimeoutExpired:
        print("⚠️ Fresh tests timed out (>5 minutes)")
        print("Falling back to existing data analysis...")
    except Exception as e:
        print(f"⚠️ Error running fresh tests: {e}")
        print("Falling back to existing data analysis...")
    
    # Fallback to standalone analysis
    try:
        print("\n📊 Running analysis on existing test data...")
        result = subprocess.run([sys.executable, 'standalone_report_generator.py'], 
                              capture_output=True, text=True, timeout=120)
        
        if result.returncode == 0:
            print("✅ Existing data analysis completed successfully!")
            print("\n" + result.stdout)
            return True
        else:
            print(f"❌ Existing data analysis failed with return code {result.returncode}")
            print("Error output:", result.stderr)
            return False
            
    except FileNotFoundError:
        print("❌ standalone_report_generator.py not found")
        return False
    except subprocess.TimeoutExpired:
        print("❌ Analysis timed out")
        return False
    except Exception as e:
        print(f"❌ Error running analysis: {e}")
        return False

def main():
    """Main entry point."""
    print(f"Started at: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    
    # Change to script directory to ensure relative imports work
    script_dir = os.path.dirname(os.path.abspath(__file__))
    original_cwd = os.getcwd()
    os.chdir(script_dir)
    
    try:
        success = run_with_fallback()
        
        if success:
            print("\n" + "=" * 50)
            print("🎉 REPORT GENERATION COMPLETED SUCCESSFULLY!")
            print("=" * 50)
            print("\n📄 Generated Files:")
            
            # List generated files
            report_files = []
            for filename in os.listdir('.'):
                if filename.endswith('_report.html'):
                    report_files.append(filename)
                    print(f"  📊 {filename}")
            
            if report_files:
                print(f"\n🌐 Open any of the HTML files above in your web browser to view the detailed analysis.")
                print("📈 The report includes:")
                print("  • Overall performance metrics")
                print("  • Detailed problem case analysis")
                print("  • Specific molecular examples")
                print("  • Group-level performance breakdown")
                print("  • Actionable recommendations")
            else:
                print("⚠️ No report files found - check console output for errors")
            
        else:
            print("\n" + "=" * 50)
            print("❌ REPORT GENERATION FAILED")
            print("=" * 50)
            print("\n🔧 Troubleshooting:")
            print("  1. Ensure you're in the correct directory")
            print("  2. Check that test data files exist (*.csv)")
            print("  3. Verify Python dependencies are installed")
            print("  4. Run the individual scripts manually for more details")
        
        return success
        
    finally:
        # Restore original working directory
        os.chdir(original_cwd)

if __name__ == "__main__":
    success = main()
    sys.exit(0 if success else 1)