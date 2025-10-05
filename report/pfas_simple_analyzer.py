"""
PFAS Algorithm Analysis Report Generator (using existing test results)

This script analyzes the existing test results to generate a comprehensive report
without needing to re-run the tests.
"""

import pandas as pd
import json
import sys
import os
from datetime import datetime
from collections import defaultdict, Counter

def load_existing_results():
    """Load existing test results from CSV files."""
    results = {}
    
    # Load specificity test results
    try:
        results['specificity'] = pd.read_csv('specificity_test_results.csv')
        print(f"Loaded specificity results: {len(results['specificity'])} tests")
    except FileNotFoundError:
        print("specificity_test_results.csv not found")
        results['specificity'] = None
    
    # Load test summary
    try:
        with open('test_summary_report.json', 'r') as f:
            results['summary'] = json.load(f)
        print("Loaded test summary report")
    except FileNotFoundError:
        print("test_summary_report.json not found")
        results['summary'] = None
    
    return results

def analyze_existing_results(results):
    """Analyze the existing test results."""
    analysis = {}
    
    if results['specificity'] is not None:
        df = results['specificity']
        valid_tests = df[df['valid_smiles'] == True]
        
        analysis['specificity_analysis'] = {
            'total_tests': len(df),
            'valid_tests': len(valid_tests),
            'detection_rate': valid_tests['expected_group_detected'].mean() if len(valid_tests) > 0 else 0,
            'specificity_rate': valid_tests['is_specific'].mean() if len(valid_tests) > 0 else 0,
            'avg_groups_detected': valid_tests['n_detected_groups'].mean() if len(valid_tests) > 0 else 0,
            'false_negatives': len(valid_tests[valid_tests['expected_group_detected'] == False]),
            'low_specificity_cases': len(valid_tests[valid_tests['n_detected_groups'] > 3])
        }
        
        # Analyze misidentifications
        false_negatives = valid_tests[valid_tests['expected_group_detected'] == False]
        low_specificity = valid_tests[valid_tests['n_detected_groups'] > 3]
        
        analysis['misidentifications'] = {
            'false_negatives': false_negatives.to_dict('records'),
            'low_specificity': low_specificity.to_dict('records')
        }
        
        # Group performance analysis
        group_performance = defaultdict(lambda: {'total': 0, 'detected': 0, 'other_groups': []})
        
        for _, row in valid_tests.iterrows():
            expected_groups = eval(row['group_ids']) if isinstance(row['group_ids'], str) else row['group_ids']
            detected_groups = eval(row['detected_groups']) if isinstance(row['detected_groups'], str) else row['detected_groups']
            
            for group_id in expected_groups:
                group_performance[group_id]['total'] += 1
                if group_id in detected_groups:
                    group_performance[group_id]['detected'] += 1
                
                other_groups = [g for g in detected_groups if g != group_id]
                group_performance[group_id]['other_groups'].extend(other_groups)
        
        # Calculate group detection rates
        for group_id in group_performance:
            perf = group_performance[group_id]
            perf['detection_rate'] = perf['detected'] / perf['total'] if perf['total'] > 0 else 0
            perf['avg_other_groups'] = len(perf['other_groups']) / perf['total'] if perf['total'] > 0 else 0
            
        analysis['group_performance'] = dict(group_performance)
    
    # Add summary data if available
    if results['summary'] is not None:
        analysis['test_summary'] = results['summary']
    
    return analysis

def generate_detailed_report(analysis, output_file='pfas_detailed_analysis_report.html'):
    """Generate a detailed HTML report."""
    
    html_content = f"""
<!DOCTYPE html>
<html>
<head>
    <title>PFAS Algorithm Detailed Analysis Report</title>
    <style>
        body {{ font-family: Arial, sans-serif; margin: 40px; line-height: 1.6; }}
        .header {{ text-align: center; margin-bottom: 30px; }}
        .section {{ margin-bottom: 40px; }}
        .metric-grid {{ display: grid; grid-template-columns: repeat(auto-fit, minmax(200px, 1fr)); gap: 15px; margin: 20px 0; }}
        .metric-card {{ border: 1px solid #ddd; padding: 15px; border-radius: 8px; background: #f9f9f9; }}
        .metric-title {{ font-weight: bold; margin-bottom: 10px; color: #333; }}
        .metric-value {{ font-size: 20px; font-weight: bold; }}
        .good {{ color: #28a745; }}
        .warning {{ color: #ffc107; }}
        .bad {{ color: #dc3545; }}
        table {{ width: 100%; border-collapse: collapse; margin-top: 15px; }}
        th, td {{ border: 1px solid #ddd; padding: 8px; text-align: left; }}
        th {{ background-color: #f2f2f2; font-weight: bold; }}
        .smiles {{ font-family: monospace; font-size: 11px; max-width: 200px; word-break: break-all; }}
        .group-list {{ font-size: 12px; }}
        .summary-stats {{ background: #e9ecef; padding: 20px; border-radius: 8px; }}
        .highlight {{ background-color: #fff3cd; padding: 10px; border-left: 4px solid #ffc107; margin: 10px 0; }}
        .error {{ background-color: #f8d7da; padding: 10px; border-left: 4px solid #dc3545; margin: 10px 0; }}
        .success {{ background-color: #d1eddd; padding: 10px; border-left: 4px solid #28a745; margin: 10px 0; }}
    </style>
</head>
<body>
    <div class="header">
        <h1>PFAS Group Identification Algorithm</h1>
        <h2>Detailed Performance Analysis Report</h2>
        <p>Generated on: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}</p>
    </div>
    """
    
    # Executive Summary
    if 'specificity_analysis' in analysis:
        spec_analysis = analysis['specificity_analysis']
        
        html_content += f"""
    <div class="section">
        <h2>Executive Summary</h2>
        <div class="summary-stats">
            <div class="metric-grid">
                <div class="metric-card">
                    <div class="metric-title">Total Tests</div>
                    <div class="metric-value">{spec_analysis['total_tests']}</div>
                </div>
                <div class="metric-card">
                    <div class="metric-title">Valid Tests</div>
                    <div class="metric-value">{spec_analysis['valid_tests']}</div>
                </div>
                <div class="metric-card">
                    <div class="metric-title">Detection Rate</div>
                    <div class="metric-value {'good' if spec_analysis['detection_rate'] > 0.8 else 'warning' if spec_analysis['detection_rate'] > 0.6 else 'bad'}">
                        {spec_analysis['detection_rate']:.1%}
                    </div>
                </div>
                <div class="metric-card">
                    <div class="metric-title">Specificity Rate</div>
                    <div class="metric-value {'good' if spec_analysis['specificity_rate'] > 0.7 else 'warning' if spec_analysis['specificity_rate'] > 0.5 else 'bad'}">
                        {spec_analysis['specificity_rate']:.1%}
                    </div>
                </div>
                <div class="metric-card">
                    <div class="metric-title">Avg Groups per Test</div>
                    <div class="metric-value">{spec_analysis['avg_groups_detected']:.1f}</div>
                </div>
                <div class="metric-card">
                    <div class="metric-title">False Negatives</div>
                    <div class="metric-value {'bad' if spec_analysis['false_negatives'] > 10 else 'warning' if spec_analysis['false_negatives'] > 5 else 'good'}">
                        {spec_analysis['false_negatives']}
                    </div>
                </div>
            </div>
        </div>
    </div>
        """
    
    # Group Performance Analysis
    if 'group_performance' in analysis:
        html_content += """
    <div class="section">
        <h2>Individual Group Performance Analysis</h2>
        <p>This table shows how well each PFAS group is detected across all test cases.</p>
        """
        
        group_perf = analysis['group_performance']
        sorted_groups = sorted(group_perf.items(), key=lambda x: x[1]['detection_rate'], reverse=True)
        
        html_content += """
        <table>
            <thead>
                <tr>
                    <th>Group ID</th>
                    <th>Total Tests</th>
                    <th>Successful Detections</th>
                    <th>Detection Rate</th>
                    <th>Avg Other Groups Detected</th>
                    <th>Performance</th>
                </tr>
            </thead>
            <tbody>
        """
        
        for group_id, perf in sorted_groups:
            detection_rate = perf['detection_rate']
            performance_class = 'good' if detection_rate > 0.8 else 'warning' if detection_rate > 0.6 else 'bad'
            performance_text = 'Excellent' if detection_rate > 0.9 else 'Good' if detection_rate > 0.8 else 'Fair' if detection_rate > 0.6 else 'Poor'
            
            html_content += f"""
                <tr>
                    <td>{group_id}</td>
                    <td>{perf['total']}</td>
                    <td>{perf['detected']}</td>
                    <td class="{performance_class}">{detection_rate:.1%}</td>
                    <td>{perf['avg_other_groups']:.1f}</td>
                    <td class="{performance_class}">{performance_text}</td>
                </tr>
            """
        
        html_content += "</tbody></table></div>"
    
    # Misidentification Analysis
    if 'misidentifications' in analysis:
        misid = analysis['misidentifications']
        
        html_content += f"""
    <div class="section">
        <h2>Misidentification Analysis</h2>
        
        <h3>False Negatives ({len(misid['false_negatives'])} cases)</h3>
        <p>These are cases where expected PFAS groups were not detected by the algorithm.</p>
        """
        
        if len(misid['false_negatives']) > 0:
            html_content += """
            <table>
                <thead>
                    <tr>
                        <th>Expected Groups</th>
                        <th>Detected Groups</th>
                        <th>SMILES</th>
                        <th>Groups Detected</th>
                    </tr>
                </thead>
                <tbody>
            """
            
            for case in misid['false_negatives'][:20]:  # Show first 20
                expected = case.get('group_ids', [])
                detected = case.get('detected_groups', [])
                smiles = case.get('smiles', '')
                n_detected = case.get('n_detected_groups', 0)
                
                html_content += f"""
                    <tr>
                        <td class="group-list">{expected}</td>
                        <td class="group-list">{detected if detected else 'None'}</td>
                        <td class="smiles">{smiles[:60]}{'...' if len(smiles) > 60 else ''}</td>
                        <td>{n_detected}</td>
                    </tr>
                """
            
            html_content += "</tbody></table>"
            
            if len(misid['false_negatives']) > 20:
                html_content += f"<p><em>Showing first 20 of {len(misid['false_negatives'])} false negative cases.</em></p>"
        
        html_content += f"""
        <h3>Low Specificity Cases ({len(misid['low_specificity'])} cases)</h3>
        <p>These are cases where too many PFAS groups were detected (more than 3), indicating low specificity.</p>
        """
        
        if len(misid['low_specificity']) > 0:
            html_content += """
            <table>
                <thead>
                    <tr>
                        <th>Expected Groups</th>
                        <th>All Detected Groups</th>
                        <th>Total Detected</th>
                        <th>SMILES</th>
                        <th>Specificity Issue</th>
                    </tr>
                </thead>
                <tbody>
            """
            
            # Sort by number of detected groups (worst first)
            sorted_low_spec = sorted(misid['low_specificity'], 
                                   key=lambda x: x.get('n_detected_groups', 0), reverse=True)
            
            for case in sorted_low_spec[:15]:  # Show first 15
                expected = case.get('group_ids', [])
                detected = case.get('detected_groups', [])
                smiles = case.get('smiles', '')
                n_detected = case.get('n_detected_groups', 0)
                
                specificity_issue = 'Severe' if n_detected > 8 else 'High' if n_detected > 5 else 'Moderate'
                issue_class = 'bad' if n_detected > 8 else 'warning' if n_detected > 5 else 'good'
                
                html_content += f"""
                    <tr>
                        <td class="group-list">{expected}</td>
                        <td class="group-list">{detected}</td>
                        <td class="{issue_class}">{n_detected}</td>
                        <td class="smiles">{smiles[:60]}{'...' if len(smiles) > 60 else ''}</td>
                        <td class="{issue_class}">{specificity_issue}</td>
                    </tr>
                """
            
            html_content += "</tbody></table>"
            
            if len(sorted_low_spec) > 15:
                html_content += f"<p><em>Showing first 15 of {len(misid['low_specificity'])} low specificity cases.</em></p>"
    
    # Recommendations
    html_content += generate_recommendations_section(analysis)
    
    # Close HTML
    html_content += """
</body>
</html>
    """
    
    with open(output_file, 'w', encoding='utf-8') as f:
        f.write(html_content)
    
    return output_file

def generate_recommendations_section(analysis):
    """Generate recommendations based on the analysis."""
    recommendations = """
    <div class="section">
        <h2>Recommendations for Algorithm Improvement</h2>
    """
    
    if 'specificity_analysis' in analysis:
        spec = analysis['specificity_analysis']
        
        if spec['detection_rate'] < 0.8:
            recommendations += """
            <div class="error">
                <h3>🔴 Critical: Low Detection Rate</h3>
                <p>Detection rate is below 80%. Priority actions:</p>
                <ul>
                    <li>Review SMARTS patterns for frequently missed groups</li>
                    <li>Check if molecular preprocessing is removing essential features</li>
                    <li>Validate that test molecules are representative of real PFAS compounds</li>
                </ul>
            </div>
            """
        
        if spec['specificity_rate'] < 0.7:
            recommendations += """
            <div class="highlight">
                <h3>⚠️ Warning: Low Specificity</h3>
                <p>Specificity rate is below 70%. Recommended actions:</p>
                <ul>
                    <li>Make SMARTS patterns more specific to avoid false positives</li>
                    <li>Review group hierarchy and overlapping definitions</li>
                    <li>Consider implementing group priority rules</li>
                </ul>
            </div>
            """
        
        if spec['avg_groups_detected'] > 4:
            recommendations += """
            <div class="highlight">
                <h3>⚠️ Warning: High Average Group Detection</h3>
                <p>Algorithm detects too many groups per molecule on average. Consider:</p>
                <ul>
                    <li>Implementing mutual exclusion rules for conflicting groups</li>
                    <li>Adding confidence scoring to rank group matches</li>
                    <li>Reviewing overlapping group definitions</li>
                </ul>
            </div>
            """
    
    if 'group_performance' in analysis:
        group_perf = analysis['group_performance']
        worst_performers = sorted(group_perf.items(), key=lambda x: x[1]['detection_rate'])[:5]
        
        if worst_performers and worst_performers[0][1]['detection_rate'] < 0.5:
            recommendations += f"""
            <div class="error">
                <h3>🔴 Priority Groups for Improvement</h3>
                <p>The following groups show consistently low detection rates:</p>
                <ul>
            """
            
            for group_id, perf in worst_performers:
                if perf['detection_rate'] < 0.7:
                    recommendations += f"""
                    <li>Group {group_id}: {perf['detection_rate']:.1%} detection rate 
                        ({perf['detected']}/{perf['total']} successful)</li>
                    """
            
            recommendations += """
                </ul>
                <p>Review SMARTS patterns and test cases for these groups.</p>
            </div>
            """
    
    # General recommendations
    recommendations += """
        <div class="success">
            <h3>✅ General Improvement Strategies</h3>
            <ul>
                <li><strong>SMARTS Pattern Review:</strong> Regularly validate SMARTS patterns against known PFAS structures</li>
                <li><strong>Test Case Expansion:</strong> Add more diverse test molecules, especially for poorly performing groups</li>
                <li><strong>Hierarchical Classification:</strong> Implement parent-child relationships between groups to reduce conflicts</li>
                <li><strong>Confidence Scoring:</strong> Add confidence metrics to help users interpret results</li>
                <li><strong>Performance Monitoring:</strong> Regular automated testing with performance benchmarks</li>
            </ul>
        </div>
    </div>
    """
    
    return recommendations

def main():
    """Main function to generate the analysis report."""
    print("PFAS Algorithm Analysis Report Generator")
    print("=" * 50)
    
    # Load existing results
    results = load_existing_results()
    
    if results['specificity'] is None and results['summary'] is None:
        print("No test results found. Please run the tests first.")
        return
    
    # Analyze results
    print("\nAnalyzing test results...")
    analysis = analyze_existing_results(results)
    
    # Generate report
    print("Generating detailed HTML report...")
    report_file = generate_detailed_report(analysis)
    
    print(f"\n✅ Analysis complete!")
    print(f"📊 Detailed report saved as: {report_file}")
    
    # Print summary to console
    if 'specificity_analysis' in analysis:
        spec = analysis['specificity_analysis']
        print(f"\n📈 Key Performance Metrics:")
        print(f"   Detection Rate: {spec['detection_rate']:.1%}")
        print(f"   Specificity Rate: {spec['specificity_rate']:.1%}")
        print(f"   False Negatives: {spec['false_negatives']}")
        print(f"   Low Specificity Cases: {spec['low_specificity_cases']}")
        print(f"   Average Groups per Test: {spec['avg_groups_detected']:.1f}")
    
    return report_file

if __name__ == "__main__":
    main()