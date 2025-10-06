"""
Standalone PFAS Algorithm Report Generator

This script generates a comprehensive test report using existing test data,
avoiding RDKit dependency issues while still providing detailed analysis.
"""

import pandas as pd
import json
import os
from datetime import datetime
from collections import Counter
import traceback

def load_existing_test_data(data_dir='.'):
    """
    Load existing test data from CSV files.
    
    Args:
        data_dir: Directory containing test result files
    
    Returns:
        dict: Test data from available files
    """
    print("📂 Loading Existing Test Data...")
    
    data = {
        'specificity_results': None,
        'oecd_results': None,
        'generic_results': None,
        'test_summary': None,
        'timestamp': datetime.now(),
        'data_sources': []
    }
    
    # Look for test result files
    test_files = {
        'specificity_test_results.csv': 'specificity_results',
        'oecd_test_results.csv': 'oecd_results', 
        'generic_test_results.csv': 'generic_results',
        'test_summary_report.json': 'test_summary'
    }
    
    # Also check parent directory and PFASgroups/tests
    search_dirs = [
        data_dir,
        os.path.join(data_dir, '..', 'PFASgroups', 'tests'),
        os.path.join(data_dir, '..', '..', 'PFASgroups', 'tests'),
        os.path.dirname(data_dir)
    ]
    
    for filename, data_key in test_files.items():
        found = False
        for search_dir in search_dirs:
            filepath = os.path.join(search_dir, filename)
            if os.path.exists(filepath):
                try:
                    if filename.endswith('.csv'):
                        df = pd.read_csv(filepath)
                        data[data_key] = df
                        print(f"✓ Loaded {filename}: {len(df)} rows")
                    elif filename.endswith('.json'):
                        with open(filepath, 'r') as f:
                            json_data = json.load(f)
                        data[data_key] = json_data
                        print(f"✓ Loaded {filename}")
                    
                    data['data_sources'].append(filepath)
                    found = True
                    break
                except Exception as e:
                    print(f"⚠️ Error loading {filepath}: {e}")
        
        if not found:
            print(f"⚠️ Could not find {filename}")
    
    return data

def analyze_existing_data(data):
    """
    Analyze existing test data and extract metrics and problem cases.
    
    Args:
        data: Test data from load_existing_test_data()
    
    Returns:
        dict: Comprehensive analysis
    """
    print("\n📊 Analyzing Existing Test Data...")
    
    analysis = {
        'overall_metrics': {},
        'problem_cases': {
            'false_negatives': [],
            'false_positives': [],
            'low_specificity': [],
            'severe_specificity': []
        },
        'group_performance': {},
        'data_sources': data['data_sources']
    }
    
    # Analyze specificity results if available
    if data['specificity_results'] is not None:
        df = data['specificity_results']
        print(f"Analyzing specificity data: {len(df)} total cases")
        
        # Handle different column name variations
        valid_col = 'valid_smiles' if 'valid_smiles' in df.columns else 'valid'
        detected_col = 'expected_group_detected' if 'expected_group_detected' in df.columns else 'detected'
        specific_col = 'is_specific' if 'is_specific' in df.columns else 'specific'
        n_detected_col = 'n_detected_groups' if 'n_detected_groups' in df.columns else 'n_detected'
        
        if valid_col in df.columns:
            valid_tests = df[df[valid_col] == True]
        else:
            valid_tests = df  # Assume all are valid if column missing
        
        if len(valid_tests) > 0:
            # Calculate metrics
            detection_rate = valid_tests[detected_col].mean() if detected_col in valid_tests.columns else 0
            specificity_rate = valid_tests[specific_col].mean() if specific_col in valid_tests.columns else 0
            avg_detected = valid_tests[n_detected_col].mean() if n_detected_col in valid_tests.columns else 0
            
            analysis['overall_metrics'] = {
                'total_tests': len(valid_tests),
                'detection_rate': detection_rate,
                'specificity_rate': specificity_rate,
                'average_detected_groups': avg_detected,
                'data_source': 'specificity_test_results.csv'
            }
            
            print(f"  Detection rate: {detection_rate:.1%}")
            print(f"  Specificity rate: {specificity_rate:.1%}")
            print(f"  Average groups per test: {avg_detected:.1f}")
            
            # Identify problem cases
            if detected_col in valid_tests.columns:
                false_negatives = valid_tests[valid_tests[detected_col] == False]
                print(f"  False negatives: {len(false_negatives)}")
                
                # Extract problem case details
                for _, row in false_negatives.head(15).iterrows():
                    # Parse group lists - handle string representation
                    expected_groups = parse_group_list(row.get('group_ids', []))
                    detected_groups = parse_group_list(row.get('detected_groups', []))
                    
                    case = {
                        'expected_groups': expected_groups,
                        'detected_groups': detected_groups,
                        'smiles': row.get('smiles', 'N/A'),
                        'n_detected': row.get(n_detected_col, 0),
                        'missing_groups': list(set(expected_groups) - set(detected_groups))
                    }
                    analysis['problem_cases']['false_negatives'].append(case)
            
            if specific_col in valid_tests.columns:
                false_positives = valid_tests[valid_tests[specific_col] == False]
                print(f"  False positives: {len(false_positives)}")
                
                for _, row in false_positives.head(15).iterrows():
                    expected_groups = parse_group_list(row.get('group_ids', []))
                    detected_groups = parse_group_list(row.get('detected_groups', []))
                    
                    case = {
                        'expected_groups': expected_groups,
                        'detected_groups': detected_groups,
                        'smiles': row.get('smiles', 'N/A'),
                        'n_detected': row.get(n_detected_col, 0),
                        'extra_groups': list(set(detected_groups) - set(expected_groups))
                    }
                    analysis['problem_cases']['false_positives'].append(case)
            
            if n_detected_col in valid_tests.columns:
                low_specificity = valid_tests[valid_tests[n_detected_col] > 3]
                severe_specificity = valid_tests[valid_tests[n_detected_col] > 5]
                
                print(f"  Low specificity (>3 groups): {len(low_specificity)}")
                print(f"  Severe specificity (>5 groups): {len(severe_specificity)}")
                
                for case_type, cases_df in [('low_specificity', low_specificity), ('severe_specificity', severe_specificity)]:
                    for _, row in cases_df.head(10).iterrows():
                        expected_groups = parse_group_list(row.get('group_ids', []))
                        detected_groups = parse_group_list(row.get('detected_groups', []))
                        
                        case = {
                            'expected_groups': expected_groups,
                            'detected_groups': detected_groups,
                            'smiles': row.get('smiles', 'N/A'),
                            'n_detected': row.get(n_detected_col, 0),
                            'extra_groups': list(set(detected_groups) - set(expected_groups))
                        }
                        analysis['problem_cases'][case_type].append(case)
            
            # Group performance analysis
            missing_counter = Counter()
            extra_counter = Counter()
            
            for case in analysis['problem_cases']['false_negatives']:
                for group in case['missing_groups']:
                    missing_counter[group] += 1
            
            for case_type in ['false_positives', 'low_specificity', 'severe_specificity']:
                for case in analysis['problem_cases'][case_type]:
                    for group in case['extra_groups']:
                        extra_counter[group] += 1
            
            analysis['group_performance'] = {
                'most_missed_groups': missing_counter.most_common(10),
                'most_over_detected_groups': extra_counter.most_common(10),
                'total_missing_detections': sum(missing_counter.values()),
                'total_over_detections': sum(extra_counter.values())
            }
    
    # Analyze OECD/Generic results if available
    for result_type in ['oecd_results', 'generic_results']:
        if data[result_type] is not None:
            df = data[result_type]
            detected_col = 'detected' if 'detected' in df.columns else 'success'
            if detected_col in df.columns:
                detection_rate = df[detected_col].mean()
                analysis['overall_metrics'][f'{result_type.split("_")[0]}_detection_rate'] = detection_rate
                print(f"  {result_type.replace('_', ' ').title()}: {detection_rate:.1%}")
    
    return analysis

def parse_group_list(group_data):
    """Parse group list from various formats (string, list, etc.)."""
    if isinstance(group_data, list):
        return [int(x) for x in group_data if isinstance(x, (int, str)) and str(x).isdigit()]
    elif isinstance(group_data, str):
        try:
            # Try to evaluate as Python literal
            import ast
            parsed = ast.literal_eval(group_data)
            if isinstance(parsed, list):
                return [int(x) for x in parsed if isinstance(x, (int, str)) and str(x).isdigit()]
        except:
            pass
        
        # Try to extract numbers from string
        import re
        numbers = re.findall(r'\d+', group_data)
        return [int(x) for x in numbers]
    elif isinstance(group_data, (int, float)):
        return [int(group_data)]
    else:
        return []

def generate_report_from_existing_data(data, analysis, output_file='pfas_algorithm_existing_data_report.html'):
    """
    Generate HTML report from existing data analysis.
    
    Args:
        data: Test data
        analysis: Analysis results
        output_file: Output file path
    
    Returns:
        str: Path to generated report
    """
    print(f"\n📄 Generating Report from Existing Data: {output_file}")
    
    timestamp = data['timestamp'].strftime('%Y-%m-%d %H:%M:%S')
    metrics = analysis['overall_metrics']
    
    # Load group names if possible
    group_names = {}
    try:
        # Try to find group names file
        search_dirs = ['.', '..', '../PFASgroups/data', '../../PFASgroups/data']
        for search_dir in search_dirs:
            group_file = os.path.join(search_dir, 'PFAS_groups_smarts.json')
            if os.path.exists(group_file):
                with open(group_file, 'r') as f:
                    group_data = json.load(f)
                group_names = {group['id']: group['name'] for group in group_data}
                break
    except:
        pass
    
    def format_groups(group_list, with_names=True):
        if not group_list:
            return "None"
        if with_names and group_names:
            return ", ".join([f"G{g} ({group_names.get(g, 'Unknown')})" for g in sorted(group_list)])
        else:
            return ", ".join([f"G{g}" for g in sorted(group_list)])
    
    def create_problem_table(cases, title, description):
        if not cases:
            return f"""
            <div class="section">
                <h3>{title}</h3>
                <div class="success">✅ No issues found - {description}</div>
            </div>
            """
        
        rows = []
        for case in cases[:15]:  # Show first 15 cases
            smiles = case.get('smiles', 'N/A')
            smiles_short = smiles[:60] + "..." if len(smiles) > 60 else smiles
            
            # Create problem-specific columns
            if 'missing_groups' in case:
                problem_info = f"<td><strong>{format_groups(case['missing_groups'], False)}</strong></td>"
                problem_header = "Missing Groups"
            else:
                problem_info = f"<td><strong>{format_groups(case.get('extra_groups', []), False)}</strong></td>"
                problem_header = "Extra Groups"
            
            rows.append(f"""
                <tr>
                    <td class="smiles">{smiles_short}</td>
                    <td>{format_groups(case.get('expected_groups', []), False)}</td>
                    <td>{format_groups(case.get('detected_groups', []), False)}</td>
                    {problem_info}
                    <td>{case.get('n_detected', 0)}</td>
                </tr>
            """)
        
        showing_text = f"Showing first 15 of {len(cases)} cases" if len(cases) > 15 else f"All {len(cases)} cases"
        
        return f"""
        <div class="section">
            <h3>{title} ({len(cases)} total)</h3>
            <p>{description}</p>
            <p><em>{showing_text}</em></p>
            
            <table class="problem-table">
                <thead>
                    <tr>
                        <th>SMILES (truncated)</th>
                        <th>Expected Groups</th>
                        <th>Detected Groups</th>
                        <th>{problem_header}</th>
                        <th># Detected</th>
                    </tr>
                </thead>
                <tbody>
                    {''.join(rows)}
                </tbody>
            </table>
        </div>
        """
    
    # Generate HTML content
    html_content = f"""
<!DOCTYPE html>
<html>
<head>
    <title>PFAS Algorithm Test Report (Existing Data)</title>
    <style>
        body {{ font-family: 'Segoe UI', Arial, sans-serif; margin: 40px; line-height: 1.6; color: #333; }}
        .header {{ text-align: center; margin-bottom: 40px; padding: 30px; background: linear-gradient(135deg, #667eea 0%, #764ba2 100%); color: white; border-radius: 10px; }}
        .header h1 {{ margin: 0; font-size: 2.5em; font-weight: 300; }}
        .header p {{ margin: 10px 0 0 0; font-size: 1.1em; opacity: 0.9; }}
        .section {{ margin-bottom: 40px; }}
        .metric-grid {{ display: grid; grid-template-columns: repeat(auto-fit, minmax(250px, 1fr)); gap: 20px; margin: 30px 0; }}
        .metric-card {{ border: 1px solid #e1e5e9; padding: 25px; border-radius: 10px; background: #f8f9fa; text-align: center; box-shadow: 0 2px 4px rgba(0,0,0,0.1); }}
        .metric-title {{ font-weight: 600; margin-bottom: 15px; color: #495057; font-size: 0.9em; text-transform: uppercase; letter-spacing: 0.5px; }}
        .metric-value {{ font-size: 2.2em; font-weight: 700; margin-bottom: 5px; }}
        .metric-subtitle {{ font-size: 0.85em; color: #6c757d; }}
        .excellent {{ color: #28a745; }}
        .good {{ color: #28a745; }}
        .warning {{ color: #ffc107; }}
        .danger {{ color: #dc3545; }}
        .alert {{ padding: 20px; margin: 25px 0; border-radius: 8px; border-left: 5px solid; }}
        .alert-success {{ background-color: #d1eddd; border-color: #28a745; }}
        .alert-warning {{ background-color: #fff3cd; border-color: #ffc107; }}
        .alert-info {{ background-color: #d1ecf1; border-color: #17a2b8; }}
        .problem-table {{ width: 100%; border-collapse: collapse; margin: 20px 0; font-size: 13px; }}
        .problem-table th, .problem-table td {{ padding: 12px; text-align: left; border: 1px solid #dee2e6; }}
        .problem-table th {{ background-color: #f8f9fa; font-weight: 600; color: #495057; }}
        .problem-table tr:nth-child(even) {{ background-color: #f8f9fa; }}
        .problem-table tr:hover {{ background-color: #e9ecef; }}
        .smiles {{ font-family: 'Courier New', monospace; font-size: 11px; max-width: 300px; word-break: break-all; }}
        .summary-box {{ background: #f8f9fa; padding: 25px; border-radius: 8px; margin: 25px 0; border: 1px solid #e9ecef; }}
        .summary-box h3 {{ margin-top: 0; color: #495057; }}
        .success {{ background-color: #d1eddd; color: #155724; padding: 15px; border-radius: 5px; margin: 15px 0; }}
        .group-list {{ columns: 3; column-gap: 20px; }}
        .status-badge {{ display: inline-block; padding: 4px 12px; border-radius: 20px; font-size: 0.8em; font-weight: 600; text-transform: uppercase; }}
        .status-excellent {{ background-color: #d1eddd; color: #155724; }}
        .status-good {{ background-color: #d4edda; color: #155724; }}
        .status-warning {{ background-color: #fff3cd; color: #856404; }}
        .status-danger {{ background-color: #f8d7da; color: #721c24; }}
    </style>
</head>
<body>
    <div class="header">
        <h1>🧪 PFAS Algorithm Test Report</h1>
        <p>Analysis of Existing Test Data</p>
        <p><strong>Generated:</strong> {timestamp}</p>
    </div>
    
    <div class="alert alert-info">
        <h3>📂 Data Sources</h3>
        <p>This report analyzes existing test data from:</p>
        <ul>
    """
    
    for source in analysis['data_sources']:
        html_content += f"<li><code>{os.path.basename(source)}</code></li>"
    
    html_content += """
        </ul>
        <p>To generate a report with fresh test data, use the automated test runner script.</p>
    </div>
    """
    
    # Performance metrics
    html_content += f"""
    <div class="section">
        <h2>📊 Performance Metrics</h2>
        <div class="metric-grid">
            <div class="metric-card">
                <div class="metric-title">Detection Rate</div>
                <div class="metric-value {'excellent' if metrics.get('detection_rate', 0) >= 0.95 else 'good' if metrics.get('detection_rate', 0) >= 0.85 else 'warning' if metrics.get('detection_rate', 0) >= 0.70 else 'danger'}">{metrics.get('detection_rate', 0):.1%}</div>
                <div class="metric-subtitle">Expected groups found</div>
            </div>
            <div class="metric-card">
                <div class="metric-title">Specificity Rate</div>
                <div class="metric-value {'excellent' if metrics.get('specificity_rate', 0) >= 0.90 else 'good' if metrics.get('specificity_rate', 0) >= 0.80 else 'warning' if metrics.get('specificity_rate', 0) >= 0.70 else 'danger'}">{metrics.get('specificity_rate', 0):.1%}</div>
                <div class="metric-subtitle">Accurate targeting</div>
            </div>
            <div class="metric-card">
                <div class="metric-title">Avg Groups per Test</div>
                <div class="metric-value {'excellent' if metrics.get('average_detected_groups', 0) <= 2 else 'good' if metrics.get('average_detected_groups', 0) <= 3 else 'warning' if metrics.get('average_detected_groups', 0) <= 5 else 'danger'}">{metrics.get('average_detected_groups', 0):.1f}</div>
                <div class="metric-subtitle">Detection specificity</div>
            </div>
            <div class="metric-card">
                <div class="metric-title">Total Tests</div>
                <div class="metric-value">{metrics.get('total_tests', 0)}</div>
                <div class="metric-subtitle">Molecules analyzed</div>
            </div>
        </div>
    </div>
    """
    
    # Algorithm status
    detection_rate = metrics.get('detection_rate', 0)
    specificity_rate = metrics.get('specificity_rate', 0)
    avg_groups = metrics.get('average_detected_groups', 0)
    
    if detection_rate >= 0.95 and specificity_rate >= 0.90 and avg_groups <= 2:
        status = "EXCELLENT"
        status_class = "status-excellent"
        status_desc = "Algorithm performs exceptionally well"
    elif detection_rate >= 0.85 and specificity_rate >= 0.80 and avg_groups <= 3:
        status = "GOOD"
        status_class = "status-good"
        status_desc = "Algorithm performs well with minor improvements possible"
    elif detection_rate >= 0.70 and specificity_rate >= 0.70:
        status = "NEEDS IMPROVEMENT"
        status_class = "status-warning"
        status_desc = "Algorithm shows promise but needs optimization"
    else:
        status = "REQUIRES ATTENTION"
        status_class = "status-danger"
        status_desc = "Algorithm has significant issues"
    
    html_content += f"""
    <div class="summary-box">
        <h3>🎯 Algorithm Status</h3>
        <span class="status-badge {status_class}">{status}</span>
        <p>{status_desc}</p>
    </div>
    """
    
    # Problem cases
    problem_cases = analysis['problem_cases']
    
    html_content += f"""
    <div class="section">
        <h2>❌ Problem Case Analysis</h2>
        
        {create_problem_table(
            problem_cases['false_negatives'],
            "False Negative Cases",
            "Molecules where expected PFAS groups were not detected"
        )}
        
        {create_problem_table(
            problem_cases['false_positives'],
            "False Positive Cases", 
            "Molecules where unexpected groups were also detected"
        )}
        
        {create_problem_table(
            problem_cases['low_specificity'],
            "Low Specificity Cases (>3 groups)",
            "Molecules with too many group detections"
        )}
        
        {create_problem_table(
            problem_cases['severe_specificity'],
            "Severe Specificity Cases (>5 groups)",
            "Molecules with severe over-detection issues"
        )}
    </div>
    """
    
    # Group performance
    group_perf = analysis['group_performance']
    
    html_content += f"""
    <div class="section">
        <h2>🔍 Group Performance Analysis</h2>
        
        <div class="summary-box">
            <h3>📉 Most Frequently Missing Groups</h3>
            <div class="group-list">
                <ul>
    """
    
    for group_id, count in group_perf.get('most_missed_groups', [])[:10]:
        group_name = group_names.get(group_id, f"Unknown")
        html_content += f"<li><strong>Group {group_id}</strong> ({group_name}): {count} cases</li>"
    
    html_content += f"""
                </ul>
            </div>
        </div>
        
        <div class="summary-box">
            <h3>📈 Most Frequently Over-Detected Groups</h3>
            <div class="group-list">
                <ul>
    """
    
    for group_id, count in group_perf.get('most_over_detected_groups', [])[:10]:
        group_name = group_names.get(group_id, f"Unknown")
        html_content += f"<li><strong>Group {group_id}</strong> ({group_name}): {count} cases</li>"
    
    html_content += """
                </ul>
            </div>
        </div>
    </div>
    
    <div class="section">
        <h2>🔄 Running Fresh Tests</h2>
        <div class="alert alert-info">
            <h3>To generate a report with fresh test data:</h3>
            <p><code>python automated_test_runner.py</code></p>
            <p>This will run new tests and generate an updated report with current algorithm performance.</p>
        </div>
    </div>
    
</body>
</html>
    """
    
    # Save report
    with open(output_file, 'w', encoding='utf-8') as f:
        f.write(html_content)
    
    print(f"✅ Report generated: {output_file}")
    return output_file

def main():
    """Main function for standalone report generation."""
    print("📊 PFAS Algorithm Standalone Report Generator")
    print("=" * 50)
    print("Analyzing existing test data and generating comprehensive report...")
    print()
    
    try:
        # Load existing data
        data = load_existing_test_data()
        
        if not any(data[key] is not None for key in ['specificity_results', 'oecd_results', 'generic_results']):
            print("❌ No test data found!")
            print("Please ensure test result CSV files are available.")
            print("Try running the automated test runner first.")
            return False
        
        # Analyze data
        analysis = analyze_existing_data(data)
        
        # Generate report
        report_file = generate_report_from_existing_data(data, analysis)
        
        # Print summary
        print("\n" + "=" * 50)
        print("📋 ANALYSIS SUMMARY")
        print("=" * 50)
        
        metrics = analysis['overall_metrics']
        if metrics:
            print(f"🎯 Detection Rate: {metrics.get('detection_rate', 0):.1%}")
            print(f"🎯 Specificity Rate: {metrics.get('specificity_rate', 0):.1%}")
            print(f"🎯 Average Groups per Test: {metrics.get('average_detected_groups', 0):.1f}")
            print(f"📊 Total Tests: {metrics.get('total_tests', 0)}")
            
            # Problem summary
            false_negs = len(analysis['problem_cases']['false_negatives'])
            false_pos = len(analysis['problem_cases']['false_positives'])
            low_spec = len(analysis['problem_cases']['low_specificity'])
            
            print(f"\n🔍 PROBLEM CASES:")
            print(f"  False Negatives: {false_negs}")
            print(f"  False Positives: {false_pos}")
            print(f"  Low Specificity: {low_spec}")
        else:
            print("⚠️ Limited analysis due to data format issues")
        
        print(f"\n📄 Report saved as: {report_file}")
        print("\n✅ Analysis completed successfully!")
        
        return True
        
    except Exception as e:
        print(f"\n❌ Analysis failed: {e}")
        traceback.print_exc()
        return False

if __name__ == "__main__":
    success = main()
    exit(0 if success else 1)