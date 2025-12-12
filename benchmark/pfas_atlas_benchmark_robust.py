#!/usr/bin/env python3
"""
PFAS-atlas benchmark with environment handling.

This version attempts to handle RDKit environment issues and provides
a robust comparison between PFASGroups and PFAS-atlas.
"""

import sys
import os
import pandas as pd
import numpy as np
from typing import List, Dict, Tuple, Any
import json
import ast
import subprocess

def check_pfas_atlas_environment():
    """Check if PFAS-atlas can be imported and used."""
    try:
        # Try to run a simple test with mamba
        result = subprocess.run([
            'mamba', 'run', '-n', 'pfasatlas', 'python', '-c', 
            'import sys; sys.path.insert(0, "/home/luc/git/PFAS-atlas"); from classification_helper import classify_pfas_molecule; print("SUCCESS")'
        ], stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
        
        if result.returncode == 0 and "SUCCESS" in result.stdout:
            print("✓ PFAS-atlas environment check passed")
            return True
        else:
            print(f"✗ PFAS-atlas environment check failed: {result.stderr}")
            return False
    except Exception as e:
        print(f"✗ Failed to check PFAS-atlas environment: {e}")
        return False

def classify_smiles_with_pfas_atlas(smiles: str) -> Tuple[str, str]:
    """Classify a single SMILES string using PFAS-atlas via subprocess."""
    try:
        # Create a temporary Python script to run the classification
        script_content = f"""
import sys
sys.path.insert(0, "/home/luc/git/PFAS-atlas")
from classification_helper import classify_pfas_molecule
import json

smiles = "{smiles}"
try:
    result = classify_pfas_molecule(smiles)
    if result and len(result) >= 2:
        print(json.dumps({{
            "first_class": result[0] if result[0] else "Unclassified",
            "second_class": result[1] if result[1] else "Unclassified",
            "success": True
        }}))
    else:
        print(json.dumps({{
            "first_class": "Unclassified",
            "second_class": "Unclassified", 
            "success": False
        }}))
except Exception as e:
    print(json.dumps({{
        "first_class": "Error",
        "second_class": "Error",
        "success": False,
        "error": str(e)
    }}))
"""
        
        # Write script to temporary file
        script_path = '/tmp/pfas_classify_temp.py'
        with open(script_path, 'w') as f:
            f.write(script_content)
        
        # Run with mamba
        result = subprocess.run([
            'mamba', 'run', '-n', 'pfasatlas', 'python', script_path
        ], stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
        
        if result.returncode == 0:
            try:
                data = json.loads(result.stdout.strip())
                return data['first_class'], data['second_class']
            except:
                return "ParseError", "ParseError"
        else:
            return "SubprocessError", "SubprocessError"
            
    except Exception as e:
        return "Error", "Error"
    finally:
        # Clean up temporary file
        try:
            os.remove('/tmp/pfas_classify_temp.py')
        except:
            pass

def classify_batch_with_pfas_atlas(smiles_list: List[str], batch_size: int = 10) -> List[Tuple[str, str]]:
    """Classify a batch of SMILES using PFAS-atlas."""
    print(f"🤖 Classifying {len(smiles_list)} molecules with PFAS-atlas...")
    
    results = []
    failed_count = 0
    
    for i in range(0, len(smiles_list), batch_size):
        batch = smiles_list[i:i+batch_size]
        batch_results = []
        
        # Create batch classification script
        script_content = f"""
import sys
sys.path.insert(0, "/home/luc/git/PFAS-atlas")
from classification_helper import classify_pfas_molecule
import json

smiles_list = {json.dumps(batch)}
results = []

for smiles in smiles_list:
    try:
        if not smiles or smiles == '':
            results.append({{"first_class": "Invalid", "second_class": "Invalid", "success": False}})
            continue
            
        result = classify_pfas_molecule(smiles)
        if result and len(result) >= 2:
            results.append({{
                "first_class": result[0] if result[0] else "Unclassified",
                "second_class": result[1] if result[1] else "Unclassified",
                "success": True
            }})
        else:
            results.append({{
                "first_class": "Unclassified",
                "second_class": "Unclassified",
                "success": False
            }})
    except Exception as e:
        results.append({{
            "first_class": "Error",
            "second_class": "Error", 
            "success": False,
            "error": str(e)
        }})

print(json.dumps(results))
"""
        
        # Write and run batch script
        script_path = f'/tmp/pfas_batch_classify_{i}.py'
        try:
            with open(script_path, 'w') as f:
                f.write(script_content)
            
            result = subprocess.run([
                'mamba', 'run', '-n', 'pfasatlas', 'python', script_path
            ], stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
            
            if result.returncode == 0:
                try:
                    batch_data = json.loads(result.stdout.strip())
                    for item in batch_data:
                        if not item['success']:
                            failed_count += 1
                        batch_results.append((item['first_class'], item['second_class']))
                except Exception as e:
                    print(f"Warning: Failed to parse batch {i}: {e}")
                    batch_results = [("ParseError", "ParseError")] * len(batch)
                    failed_count += len(batch)
            else:
                print(f"Warning: Batch {i} failed: {result.stderr}")
                batch_results = [("SubprocessError", "SubprocessError")] * len(batch)
                failed_count += len(batch)
        
        except Exception as e:
            print(f"Warning: Exception in batch {i}: {e}")
            batch_results = [("Error", "Error")] * len(batch)
            failed_count += len(batch)
        finally:
            try:
                os.remove(script_path)
            except:
                pass
        
        results.extend(batch_results)
        
        # Progress indicator
        progress = min(int(i) + int(batch_size), len(smiles_list))
        print(f"  Progress: {progress}/{len(smiles_list)} molecules...")
    
    print(f"✓ PFAS-atlas classification complete. {failed_count}/{len(smiles_list)} failed.")
    return results

def load_and_analyze_data():
    """Load data and perform comprehensive analysis."""
    
    # Load PFASGroups results
    csv_path = "/home/luc/git/PFASGroups/PFASgroups/tests/results/specificity_test_results.csv"
    df = pd.read_csv(csv_path)
    print(f"✓ Loaded {len(df)} molecules")
    
    # Parse PFASGroups results
    df['detected_groups_parsed'] = df['detected_groups'].apply(
        lambda x: ast.literal_eval(x) if pd.notna(x) and x.strip() else []
    )
    df['pfasgroups_has_classification'] = df['detected_groups_parsed'].apply(lambda x: len(x) > 0)
    
    # Check if PFAS-atlas is available
    atlas_available = check_pfas_atlas_environment()
    
    if atlas_available:
        # Run PFAS-atlas classification
        smiles_list = df['smiles'].tolist()
        atlas_results = classify_batch_with_pfas_atlas(smiles_list)
        
        df['pfas_atlas_first_class'] = [result[0] for result in atlas_results]
        df['pfas_atlas_second_class'] = [result[1] for result in atlas_results]
        df['pfas_atlas_classified'] = df['pfas_atlas_first_class'].apply(
            lambda x: x not in ['Invalid', 'Error', 'Unclassified', 'Unknown', 'ParseError', 'SubprocessError']
        )
    else:
        print("⚠️  PFAS-atlas not available - generating analysis for PFASGroups only")
        df['pfas_atlas_first_class'] = 'Unavailable'
        df['pfas_atlas_second_class'] = 'Unavailable'
        df['pfas_atlas_classified'] = False
    
    return df, atlas_available

def generate_comparison_report(df: pd.DataFrame, atlas_available: bool) -> str:
    """Generate comprehensive comparison report."""
    
    total_molecules = len(df)
    pfasgroups_classified = int(df['pfasgroups_has_classification'].sum())
    pfasgroups_coverage = pfasgroups_classified / total_molecules * 100
    
    report = f"""
# PFAS Classification Benchmark Report
## PFASGroups vs PFAS-atlas

### Summary
- **Total molecules**: {total_molecules}
- **PFASGroups classified**: {pfasgroups_classified} ({pfasgroups_coverage:.1f}%)
"""
    
    if atlas_available:
        atlas_classified = int(df['pfas_atlas_classified'].sum())
        atlas_coverage = atlas_classified / total_molecules * 100
        
        both_classified = int(((df['pfasgroups_has_classification']) & 
                          (df['pfas_atlas_classified'])).sum())
        agreement_rate = both_classified / max(pfasgroups_classified, atlas_classified, 1) * 100
        
        report += f"""- **PFAS-atlas classified**: {atlas_classified} ({atlas_coverage:.1f}%)
- **Both systems classified**: {both_classified}
- **Agreement rate**: {agreement_rate:.1f}%

### Classification Overlap Analysis
"""
        
        # Analyze overlap patterns
        pfasgroups_only = int(((df['pfasgroups_has_classification']) & 
                          (~df['pfas_atlas_classified'])).sum())
        atlas_only = int(((~df['pfasgroups_has_classification']) & 
                     (df['pfas_atlas_classified'])).sum())
        neither = int(((~df['pfasgroups_has_classification']) & 
                  (~df['pfas_atlas_classified'])).sum())
        
        report += f"""
- **Only PFASGroups classified**: {pfasgroups_only}
- **Only PFAS-atlas classified**: {atlas_only}  
- **Neither classified**: {neither}

### Top PFAS-atlas Classifications
"""
        # Add top atlas classifications
        atlas_counts = df['pfas_atlas_first_class'].value_counts().head(10)
        for class_name, count in atlas_counts.items():
            if class_name not in ['Invalid', 'Error', 'Unclassified', 'Unavailable']:
                report += f"- **{class_name}**: {count} molecules\n"
    
    else:
        report += f"""
- **PFAS-atlas**: Not available due to environment issues

### PFASGroups Analysis Only
This analysis focuses on PFASGroups performance since PFAS-atlas could not be loaded.
"""
    
    # Add PFASGroups group analysis
    group_names = {
        18: "ethene", 19: "alkene", 20: "alkane", 21: "halide",
        23: "alkane-per", 25: "iodide-per", 27: "ketone", 28: "ketone-per",
        30: "Perfluoroalkyl ketones-per", 35: "acyl halide-per",
        41: "alkene-per", 42: "iodide-per", 48: "Perfluoroalkane-per",
        49: "Perfluoroalkene-per", 50: "alkyne-per"
    }
    
    all_groups = []
    for groups_list in df['detected_groups_parsed']:
        all_groups.extend(groups_list)
    
    group_counts = {}
    for group in all_groups:
        group_counts[group] = group_counts.get(group, 0) + 1
    
    report += f"""
### Top PFASGroups Classifications
"""
    for group_id, count in sorted(group_counts.items(), key=lambda x: x[1], reverse=True)[:10]:
        group_name = group_names.get(group_id, f"Group_{group_id}")
        report += f"- **{group_name}** (ID: {group_id}): {count} detections\n"
    
    # PFASGroups performance metrics
    if 'is_specific' in df.columns:
        specificity_rate = df['is_specific'].sum() / total_molecules * 100
        report += f"""
### PFASGroups Performance
- **Specificity rate**: {specificity_rate:.1f}%
- **Detection coverage**: {pfasgroups_coverage:.1f}%
"""
    
    return report

def main():
    """Main benchmark function."""
    print("🧪 PFAS Classification Benchmark: PFASGroups vs PFAS-atlas")
    print("=" * 60)
    
    output_dir = "/home/luc/git/PFASGroups/benchmark"
    os.makedirs(output_dir, exist_ok=True)
    
    # Load and analyze data
    df, atlas_available = load_and_analyze_data()
    
    # Generate report
    print("\n📄 Generating benchmark report...")
    report = generate_comparison_report(df, atlas_available)
    
    # Save files
    report_path = os.path.join(output_dir, "pfas_atlas_benchmark_report.md")
    with open(report_path, 'w') as f:
        f.write(report)
    print(f"✓ Report saved to {report_path}")
    
    results_path = os.path.join(output_dir, "pfas_atlas_benchmark_results.csv")
    df.to_csv(results_path, index=False)
    print(f"✓ Results saved to {results_path}")
    
    print("\n✅ Benchmark complete!")

if __name__ == "__main__":
    main()