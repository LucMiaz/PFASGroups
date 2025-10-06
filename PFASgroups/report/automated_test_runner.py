"""
Automated PFAS Algorithm Test Runner and Report Generator

This script automatically runs the PFAS group detection tests and generates
a comprehensive report showing accuracy, specificity, and detailed problem analysis.
"""

import sys
import os
import pandas as pd
import numpy as np
import json
from datetime import datetime
from collections import Counter, defaultdict
import traceback

# Add the PFASgroups path for imports
current_dir = os.path.dirname(os.path.abspath(__file__))
project_root = os.path.dirname(os.path.dirname(current_dir))
tests_dir = os.path.join(project_root, 'PFASgroups', 'tests')
sys.path.insert(0, project_root)
sys.path.insert(0, tests_dir)

# Try to import the test module
try:
    from PFASgroups.tests.test_examples import (
        df_test_pfas_group_specificity,
        create_specificity_test_molecules,
        TestPFASGroups,
        analyze_group_overlap,
        are_detected_groups_acceptable,
        load_graph_from_json
    )
    print("✓ Successfully imported test modules")
except ImportError as e:
    print(f"❌ Error importing test modules: {e}")
    print("Please ensure you're running this script from the correct directory")
    sys.exit(1)

def run_fresh_tests(output_dir=current_dir, refresh = False):
    """
    Run fresh PFAS group detection tests and return results.
    
    Args:
        output_dir: Directory to save test result files
    
    Returns:
        dict: Test results including specificity data and summary metrics
    """
    print("🔬 Running Fresh PFAS Group Detection Tests")
    print("=" * 50)
    
    results = {
        'timestamp': datetime.now(),
        'specificity_results': None,
        'oecd_results': None,
        'generic_results': None,
        'test_summary': None,
        'errors': []
    }
    
    try:
        # 1. Run specificity tests (most comprehensive)
        print("\n📊 Running specificity tests...")
        specificity_file = os.path.join(output_dir, 'fresh_specificity_test_results.csv')
        
        # Generate fresh test molecules
        test_molecules = create_specificity_test_molecules()
        print(f"Generated {len(test_molecules)} test molecules")
        
        # Run specificity tests
        if refresh is True:
            specificity_df = df_test_pfas_group_specificity(
            test_molecules=test_molecules, 
            output_file=specificity_file, 
            verbose=False
            )
        else:
            specificity_df= pd.read_csv(f"{tests_dir}/specificity_test_results.csv")
        results['specificity_results'] = specificity_df
        print(f"✓ Specificity tests completed - {len(specificity_df)} results")
        
    except Exception as e:
        error_msg = f"Specificity tests failed: {e}"
        print(f"❌ {error_msg}")
        results['errors'].append(error_msg)
        traceback.print_exc()
    
    try:
        # 2. Run OECD and Generic tests
        print("\n🧪 Running OECD and Generic group tests...")
        tester = TestPFASGroups()
        
        # Run OECD tests
        try:
            tester.test_oecd_pfas_groups()
            # Access the global test results
            from PFASgroups.tests.test_examples import TEST_SUMMARY_DATA
            results['oecd_results'] = TEST_SUMMARY_DATA.get('oecd_results', [])
            print(f"✓ OECD tests completed - {len(results['oecd_results'])} results")
        except Exception as e:
            error_msg = f"OECD tests failed: {e}"
            print(f"❌ {error_msg}")
            results['errors'].append(error_msg)
        
        # Run Generic tests
        try:
            tester.test_generic_pfas_groups()
            results['generic_results'] = TEST_SUMMARY_DATA.get('generic_results', [])
            print(f"✓ Generic tests completed - {len(results['generic_results'])} results")
        except Exception as e:
            error_msg = f"Generic tests failed: {e}"
            print(f"❌ {error_msg}")
            results['errors'].append(error_msg)
            
    except Exception as e:
        error_msg = f"OECD/Generic tests setup failed: {e}"
        print(f"❌ {error_msg}")
        results['errors'].append(error_msg)
    
    # 3. Generate test summary
    try:
        from PFASgroups.tests.test_examples import generate_test_summary
        summary = generate_test_summary(os.path.join(output_dir, 'fresh_test_summary.json'))
        results['test_summary'] = summary
        print("✓ Test summary generated")
    except Exception as e:
        error_msg = f"Test summary generation failed: {e}"
        print(f"❌ {error_msg}")
        results['errors'].append(error_msg)
    
    return results

def analyze_test_results(test_results):
    """
    Analyze test results and extract detailed metrics and problem cases.
    
    Args:
        test_results: Results from run_fresh_tests()
    
    Returns:
        dict: Comprehensive analysis including problem cases
    """
    print("\n📈 Analyzing Test Results...")
    
    analysis = {
        'overall_metrics': {},
        'specificity_analysis': {},
        'problem_cases': {
            'false_negatives': [],
            'false_positives': [],
            'low_specificity': [],
            'severe_specificity': []
        },
        'group_performance': {},
        'molecular_examples': {}
    }
    
    # Analyze specificity results (most detailed)
    if test_results['specificity_results'] is not None:
        df = test_results['specificity_results']
        valid_tests = df[df['valid_smiles'] == True]
        
        if len(valid_tests) > 0:
            # Overall metrics
            detection_rate = valid_tests['expected_group_detected'].mean()
            specificity_rate = valid_tests['is_specific'].mean()
            avg_detected_groups = valid_tests['n_detected_groups'].mean()
            
            analysis['overall_metrics'] = {
                'total_tests': len(valid_tests),
                'detection_rate': detection_rate,
                'specificity_rate': specificity_rate,
                'average_detected_groups': avg_detected_groups,
                'valid_smiles_rate': len(valid_tests) / len(df) if len(df) > 0 else 0
            }
            
            # Identify problem cases
            false_negatives = valid_tests[valid_tests['expected_group_detected'] == False]
            false_positives = valid_tests[valid_tests['is_specific'] == False]
            low_specificity = valid_tests[valid_tests['n_detected_groups'] > 3]
            severe_specificity = valid_tests[valid_tests['n_detected_groups'] > 5]
            
            # Store problem case details
            for case_type, cases_df in [
                ('false_negatives', false_negatives),
                ('false_positives', false_positives), 
                ('low_specificity', low_specificity),
                ('severe_specificity', severe_specificity)
            ]:
                case_list = []
                for _, row in cases_df.head(10).iterrows():  # Top 10 examples
                    case_info = {
                        'expected_groups': row['group_ids'],
                        'detected_groups': row['detected_groups'],
                        'smiles': row['smiles'],
                        'n_detected': row['n_detected_groups'],
                        'problem_type': case_type
                    }
                    
                    # Add specific analysis for each case type
                    if case_type == 'false_negatives':
                        missing_groups = list(set(row['group_ids']) - set(row['detected_groups']))
                        case_info['missing_groups'] = missing_groups
                    elif case_type in ['false_positives', 'low_specificity', 'severe_specificity']:
                        extra_groups = list(set(row['detected_groups']) - set(row['group_ids']))
                        case_info['extra_groups'] = extra_groups
                    
                    case_list.append(case_info)
                
                analysis['problem_cases'][case_type] = case_list
            
            # Group-level analysis
            group_stats = {}
            
            # Analyze missing groups (false negatives)
            missing_counter = Counter()
            for _, row in false_negatives.iterrows():
                missing_groups = set(row['group_ids']) - set(row['detected_groups'])
                for group in missing_groups:
                    missing_counter[group] += 1
            
            # Analyze extra groups (false positives)
            extra_counter = Counter()
            for _, row in false_positives.iterrows():
                extra_groups = set(row['detected_groups']) - set(row['group_ids'])
                for group in extra_groups:
                    extra_counter[group] += 1
            
            analysis['group_performance'] = {
                'most_missed_groups': missing_counter.most_common(10),
                'most_over_detected_groups': extra_counter.most_common(10),
                'total_missing_detections': sum(missing_counter.values()),
                'total_over_detections': sum(extra_counter.values())
            }
            
            print(f"✓ Specificity analysis completed")
            print(f"  - Detection rate: {detection_rate:.1%}")
            print(f"  - Specificity rate: {specificity_rate:.1%}")
            print(f"  - False negatives: {len(false_negatives)}")
            print(f"  - Low specificity cases: {len(low_specificity)}")
    
    # Analyze OECD results
    if test_results['oecd_results']:
        oecd_df = pd.DataFrame(test_results['oecd_results'])
        oecd_detection_rate = oecd_df['detected'].mean() if len(oecd_df) > 0 else 0
        analysis['overall_metrics']['oecd_detection_rate'] = oecd_detection_rate
        print(f"✓ OECD analysis: {oecd_detection_rate:.1%} detection rate")
    
    # Analyze Generic results  
    if test_results['generic_results']:
        generic_df = pd.DataFrame(test_results['generic_results'])
        generic_detection_rate = generic_df['detected'].mean() if len(generic_df) > 0 else 0
        analysis['overall_metrics']['generic_detection_rate'] = generic_detection_rate
        print(f"✓ Generic analysis: {generic_detection_rate:.1%} detection rate")
    
    return analysis

def generate_comprehensive_report(test_results, analysis, output_file=f'{current_dir}/pfas_algorithm_test_report.html'):
    """
    Generate a comprehensive HTML report with all test results and analysis.
    
    Args:
        test_results: Results from run_fresh_tests()
        analysis: Analysis from analyze_test_results()
        output_file: Output HTML file path
    
    Returns:
        str: Path to generated report file
    """
    print(f"\n📄 Generating Comprehensive Report: {output_file}")
    
    timestamp = test_results['timestamp'].strftime('%Y-%m-%d %H:%M:%S')
    metrics = analysis['overall_metrics']
    
    # Load group names for better display
    try:
        data_folder = os.path.join(os.path.dirname(current_dir), 'PFASgroups', 'data')
        with open(os.path.join(data_folder, 'PFAS_groups_smarts.json'), 'r') as f:
            group_data = json.load(f)
        group_names = {group['id']: group['name'] for group in group_data}
    except:
        group_names = {}
    
    def format_groups(group_list, with_names=True):
        if not group_list:
            return '<span class="group-names">None</span>'
        if with_names and group_names:
            formatted_groups = []
            for g in sorted(group_list):
                group_name = group_names.get(g, 'Unknown')
                # Truncate very long group names for better display
                if len(group_name) > 40:
                    group_name = group_name[:37] + "..."
                formatted_groups.append(f'<span class="group-id">G{g}</span> <span class="group-desc">({group_name})</span>')
            return '<div class="group-names">' + '<br>'.join(formatted_groups) + '</div>'
        else:
            return '<div class="group-names">' + ", ".join([f'<span class="group-id">G{g}</span>' for g in sorted(group_list)]) + '</div>'
    
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
            smiles_short = case['smiles'][:60] + "..." if len(case['smiles']) > 60 else case['smiles']
            
            # Create problem-specific columns
            problem_info = ""
            if 'missing_groups' in case:
                problem_info = f"<td><strong>{format_groups(case['missing_groups'])}</strong></td>"
            elif 'extra_groups' in case:
                problem_info = f"<td><strong>{format_groups(case['extra_groups'])}</strong></td>"
            else:
                problem_info = "<td>-</td>"
            
            rows.append(f"""
                <tr>
                    <td class="smiles">{smiles_short}</td>
                    <td>{format_groups(case['expected_groups'])}</td>
                    <td>{format_groups(case['detected_groups'])}</td>
                    {problem_info}
                    <td>{case['n_detected']}</td>
                </tr>
            """)
        
        table_header = "Missing Groups" if 'missing_groups' in cases[0] else "Extra Groups"
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
                        <th>{table_header}</th>
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
    <title>PFAS Algorithm Automated Test Report</title>
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
        .alert-danger {{ background-color: #f8d7da; border-color: #dc3545; }}
        .alert-info {{ background-color: #d1ecf1; border-color: #17a2b8; }}
        .problem-table {{ width: 100%; border-collapse: collapse; margin: 20px 0; font-size: 12px; table-layout: fixed; }}
        .problem-table th, .problem-table td {{ padding: 8px; text-align: left; border: 1px solid #dee2e6; word-wrap: break-word; vertical-align: top; }}
        .problem-table th {{ background-color: #f8f9fa; font-weight: 600; color: #495057; }}
        .problem-table tr:nth-child(even) {{ background-color: #f8f9fa; }}
        .problem-table tr:hover {{ background-color: #e9ecef; }}
        .problem-table th:nth-child(1), .problem-table td:nth-child(1) {{ width: 25%; }} /* SMILES column */
        .problem-table th:nth-child(2), .problem-table td:nth-child(2) {{ width: 25%; }} /* Expected groups */
        .problem-table th:nth-child(3), .problem-table td:nth-child(3) {{ width: 25%; }} /* Detected groups */
        .problem-table th:nth-child(4), .problem-table td:nth-child(4) {{ width: 20%; }} /* Missing/Extra groups */
        .problem-table th:nth-child(5), .problem-table td:nth-child(5) {{ width: 5%; text-align: center; }} /* Count */
        .smiles {{ font-family: 'Courier New', monospace; font-size: 10px; max-width: 100%; word-break: break-all; }}
        .group-names {{ font-size: 11px; line-height: 1.3; }}
        .group-id {{ font-weight: bold; color: #495057; }}
        .group-desc {{ color: #6c757d; font-style: italic; }}
        .summary-box {{ background: #f8f9fa; padding: 25px; border-radius: 8px; margin: 25px 0; border: 1px solid #e9ecef; }}
        .summary-box h3 {{ margin-top: 0; color: #495057; }}
        .error-box {{ background-color: #f8d7da; color: #721c24; padding: 15px; border-radius: 5px; margin: 15px 0; }}
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
        <p>Automated Testing & Analysis</p>
        <p><strong>Generated:</strong> {timestamp}</p>
    </div>
    
    <div class="section">
        <h2>📊 Executive Summary</h2>
        <div class="metric-grid">
            <div class="metric-card">
                <div class="metric-title">Overall Detection Rate</div>
                <div class="metric-value {'excellent' if metrics.get('detection_rate', 0) >= 0.95 else 'good' if metrics.get('detection_rate', 0) >= 0.85 else 'warning' if metrics.get('detection_rate', 0) >= 0.70 else 'danger'}">{metrics.get('detection_rate', 0):.1%}</div>
                <div class="metric-subtitle">Expected groups found</div>
            </div>
            <div class="metric-card">
                <div class="metric-title">Specificity Rate</div>
                <div class="metric-value {'excellent' if metrics.get('specificity_rate', 0) >= 0.90 else 'good' if metrics.get('specificity_rate', 0) >= 0.80 else 'warning' if metrics.get('specificity_rate', 0) >= 0.70 else 'danger'}">{metrics.get('specificity_rate', 0):.1%}</div>
                <div class="metric-subtitle">Accurate group targeting</div>
            </div>
            <div class="metric-card">
                <div class="metric-title">Average Groups per Test</div>
                <div class="metric-value {'excellent' if metrics.get('average_detected_groups', 0) <= 2 else 'good' if metrics.get('average_detected_groups', 0) <= 3 else 'warning' if metrics.get('average_detected_groups', 0) <= 5 else 'danger'}">{metrics.get('average_detected_groups', 0):.1f}</div>
                <div class="metric-subtitle">Detection specificity</div>
            </div>
            <div class="metric-card">
                <div class="metric-title">Total Tests Run</div>
                <div class="metric-value">{metrics.get('total_tests', 0)}</div>
                <div class="metric-subtitle">Test molecules analyzed</div>
            </div>
        </div>
        
        <div class="summary-box">
            <h3>🎯 Algorithm Status Assessment</h3>
    """
    
    # Algorithm status assessment
    detection_rate = metrics.get('detection_rate', 0)
    specificity_rate = metrics.get('specificity_rate', 0)
    avg_groups = metrics.get('average_detected_groups', 0)
    
    if detection_rate >= 0.95 and specificity_rate >= 0.90 and avg_groups <= 2:
        status = "EXCELLENT"
        status_class = "status-excellent"
        status_desc = "Algorithm performs exceptionally well with high accuracy and specificity"
    elif detection_rate >= 0.85 and specificity_rate >= 0.80 and avg_groups <= 3:
        status = "GOOD" 
        status_class = "status-good"
        status_desc = "Algorithm performs well with minor areas for improvement"
    elif detection_rate >= 0.70 and specificity_rate >= 0.70:
        status = "NEEDS IMPROVEMENT"
        status_class = "status-warning"
        status_desc = "Algorithm shows promise but requires optimization"
    else:
        status = "REQUIRES ATTENTION"
        status_class = "status-danger"
        status_desc = "Algorithm has significant issues that need addressing"
    
    html_content += f"""
            <span class="status-badge {status_class}">{status}</span>
            <p>{status_desc}</p>
        </div>
    </div>
    """
    
    # Add error reporting if any
    if test_results['errors']:
        html_content += f"""
        <div class="section">
            <h2>⚠️ Test Execution Issues</h2>
            <div class="alert alert-warning">
                <h3>Some tests encountered errors:</h3>
                <ul>
        """
        for error in test_results['errors']:
            html_content += f"<li>{error}</li>"
        html_content += """
                </ul>
                <p>Results may be incomplete. Check the console output for more details.</p>
            </div>
        </div>
        """
    
    # Problem case analysis
    problem_cases = analysis['problem_cases']
    
    html_content += f"""
    <div class="section">
        <h2>❌ Detection Problems (False Negatives)</h2>
        {create_problem_table(
            problem_cases['false_negatives'],
            "False Negative Cases", 
            "Test molecules where expected PFAS groups were not detected by the algorithm"
        )}
    </div>
    
    <div class="section">
        <h2>🎯 Specificity Problems (False Positives)</h2>
        {create_problem_table(
            problem_cases['false_positives'],
            "False Positive Cases",
            "Test molecules where unexpected groups were detected in addition to expected ones"
        )}
    </div>
    
    <div class="section">
        <h2>📈 Low Specificity Cases</h2>
        {create_problem_table(
            problem_cases['low_specificity'],
            "Low Specificity Cases (>3 groups detected)",
            "Test molecules where the algorithm detected many groups, indicating potential over-broad patterns"
        )}
    </div>
    
    <div class="section">
        <h2>🚨 Severe Specificity Issues</h2>
        {create_problem_table(
            problem_cases['severe_specificity'],
            "Severe Specificity Cases (>5 groups detected)",
            "Test molecules with the most severe specificity problems"
        )}
    </div>
    """
    
    # Group performance analysis
    group_perf = analysis['group_performance']
    
    html_content += f"""
    <div class="section">
        <h2>🔍 Group-Level Performance Analysis</h2>
        
        <div class="summary-box">
            <h3>📉 Most Frequently Missing Groups</h3>
            <p>Groups that are expected but not detected (false negatives):</p>
            <div class="group-list">
                <ul>
    """
    
    for group_id, count in group_perf.get('most_missed_groups', [])[:10]:
        group_name = group_names.get(group_id, f"Unknown Group {group_id}")
        percentage = (count / metrics.get('total_tests', 1)) * 100
        html_content += f"<li><strong>Group {group_id}</strong> ({group_name}): {count} cases ({percentage:.1f}%)</li>"
    
    html_content += f"""
                </ul>
            </div>
        </div>
        
        <div class="summary-box">
            <h3>📈 Most Frequently Over-Detected Groups</h3>
            <p>Groups that are detected but not expected (false positives):</p>
            <div class="group-list">
                <ul>
    """
    
    for group_id, count in group_perf.get('most_over_detected_groups', [])[:10]:
        group_name = group_names.get(group_id, f"Unknown Group {group_id}")
        percentage = (count / metrics.get('total_tests', 1)) * 100
        html_content += f"<li><strong>Group {group_id}</strong> ({group_name}): {count} cases ({percentage:.1f}%)</li>"
    
    html_content += """
                </ul>
            </div>
        </div>
    </div>
    """
    
    # Add OECD and Generic results if available
    if test_results.get('oecd_results') or test_results.get('generic_results'):
        html_content += """
        <div class="section">
            <h2>🧪 Additional Test Results</h2>
        """
        
        if test_results.get('oecd_results'):
            oecd_rate = metrics.get('oecd_detection_rate', 0)
            html_content += f"""
            <div class="summary-box">
                <h3>OECD PFAS Groups Test</h3>
                <p><strong>Detection Rate:</strong> <span class="{'excellent' if oecd_rate >= 0.8 else 'warning'}">{oecd_rate:.1%}</span></p>
                <p>Tests synthetic molecules designed to match specific OECD PFAS group definitions.</p>
            </div>
            """
        
        if test_results.get('generic_results'):
            generic_rate = metrics.get('generic_detection_rate', 0)
            html_content += f"""
            <div class="summary-box">
                <h3>Generic PFAS Groups Test</h3>
                <p><strong>Detection Rate:</strong> <span class="{'excellent' if generic_rate >= 0.8 else 'warning'}">{generic_rate:.1%}</span></p>
                <p>Tests synthetic molecules with generic functional groups in fluorinated contexts.</p>
            </div>
            """
        
        html_content += "</div>"
    
    # Recommendations section
    html_content += f"""
    <div class="section">
        <h2>🛠️ Recommendations</h2>
        
        <div class="alert {'alert-success' if status == 'EXCELLENT' else 'alert-warning' if status in ['GOOD', 'NEEDS IMPROVEMENT'] else 'alert-danger'}">
            <h3>Priority Actions</h3>
    """
    
    if status == "EXCELLENT":
        html_content += """
            <p>✅ <strong>Algorithm performing excellently!</strong> Consider these optional enhancements:</p>
            <ul>
                <li>Monitor performance on new compound types</li>
                <li>Document current performance for validation</li>
                <li>Consider expanding test coverage to edge cases</li>
            </ul>
        """
    else:
        # Generate specific recommendations based on problems found
        recommendations = []
        
        missing_groups = group_perf.get('most_missed_groups', [])
        if missing_groups:
            top_missing_with_names = [f"G{g[0]} ({group_names.get(g[0], 'Unknown')})" for g in missing_groups[:3]]
            recommendations.append(f"<strong>Fix detection issues for Groups {', '.join(top_missing_with_names)}</strong> - these are frequently not detected")
        
        over_detected = group_perf.get('most_over_detected_groups', [])
        if over_detected:
            top_over_with_names = [f"G{g[0]} ({group_names.get(g[0], 'Unknown')})" for g in over_detected[:3]]
            recommendations.append(f"<strong>Improve specificity for Groups {', '.join(top_over_with_names)}</strong> - these are frequently over-detected")
        
        if metrics.get('average_detected_groups', 0) > 3:
            recommendations.append("<strong>Review SMARTS patterns for over-broad matching</strong> - many tests detect too many groups")
        
        if metrics.get('detection_rate', 0) < 0.85:
            recommendations.append("<strong>Review test molecule generation</strong> - ensure expected groups are chemically valid")
        
        if recommendations:
            html_content += "<ul>" + "".join([f"<li>{rec}</li>" for rec in recommendations]) + "</ul>"
        else:
            html_content += "<p>Review the detailed problem cases above for specific areas to address.</p>"
    
    html_content += """
        </div>
    </div>
    
    <div class="section">
        <h2>📁 Data Files Generated</h2>
        <div class="summary-box">
            <h3>Output Files</h3>
            <ul>
                <li><code>fresh_specificity_test_results.csv</code> - Detailed specificity test results</li>
                <li><code>fresh_test_summary.json</code> - Machine-readable test summary</li>
                <li><code>pfas_algorithm_test_report.html</code> - This comprehensive report</li>
            </ul>
            <p>These files can be used for further analysis or tracking performance over time.</p>
        </div>
    </div>
    
    <div class="section">
        <h2>🔄 Re-running Tests</h2>
        <div class="alert alert-info">
            <h3>To regenerate this report with fresh data:</h3>
            <p><code>python automated_test_runner.py</code></p>
            <p>This script automatically runs all tests and generates an updated report with the latest algorithm performance.</p>
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
    """Main function to run automated testing and report generation."""
    print("🚀 PFAS Algorithm Automated Test Runner")
    print("=" * 50)
    item_id = 1
    print("This script will:")
    print(f"  {item_id}. Run fresh PFAS group detection tests")
    item_id+=1
    print(f"  {item_id}. Analyze results for accuracy and specificity")
    item_id+=1
    print(f"  {item_id}. Generate a comprehensive HTML report")
    item_id+=1
    print(f"  {item_id}. Identify specific problem cases with examples")
    print()
    # Check command line arguments
    try:
        # 1. Run fresh tests
        test_results = run_fresh_tests(refresh = True)

        # 2. Analyze results
        analysis = analyze_test_results(test_results)
        
        # 3. Generate report
        report_file = generate_comprehensive_report(test_results, analysis)
        
        # 4. Print summary
        print("\n" + "=" * 50)
        print("📋 TEST SUMMARY")
        print("=" * 50)
        
        metrics = analysis['overall_metrics']
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
        
        if test_results['errors']:
            print(f"\n⚠️  ERRORS ENCOUNTERED: {len(test_results['errors'])}")
            for error in test_results['errors']:
                print(f"    {error}")
        
        print(f"\n📄 Detailed report saved as: {report_file}")
        print("\n✅ Automated testing completed successfully!")
        
        return True
        
    except Exception as e:
        print(f"\n❌ Automated testing failed: {e}")
        traceback.print_exc()
        return False

if __name__ == "__main__":
    success = main()
    sys.exit(0 if success else 1)