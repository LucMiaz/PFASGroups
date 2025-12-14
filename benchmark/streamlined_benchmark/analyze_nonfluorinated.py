#!/usr/bin/env python3
"""
Non-Fluorinated PFAS Benchmark Analysis
Analyzes results from the non-fluorinated functional group exclusion test
"""

import json
import sys
from datetime import datetime

def load_nonfluorinated_results(filename):
    """Load non-fluorinated benchmark results from JSON"""
    
    with open(filename, 'r') as f:
        results = json.load(f)
    
    print(f"📊 Loaded non-fluorinated benchmark with {len(results)} functional group tests")
    return results

def analyze_nonfluorinated_results(results):
    """Analyze non-fluorinated benchmark results"""
    
    total_molecules = sum(r['molecules_tested'] for r in results)
    total_pfas_false_positives = sum(r['pfasgroups_detections'] for r in results)
    total_atlas_false_positives = sum(r['atlas_detections'] for r in results)
    
    analysis = {
        'total_functional_groups': len(results),
        'total_molecules': total_molecules,
        'pfasgroups_false_positives': total_pfas_false_positives,
        'atlas_false_positives': total_atlas_false_positives,
        'pfasgroups_fpr': (total_pfas_false_positives / max(total_molecules, 1)) * 100,
        'atlas_fpr': (total_atlas_false_positives / max(total_molecules, 1)) * 100,
        'group_details': []
    }
    
    for group_result in results:
        group_analysis = {
            'group_id': group_result['group_id'],
            'group_name': group_result['group_name'],
            'molecules_tested': group_result['molecules_tested'],
            'pfasgroups_false_positives': group_result['pfasgroups_detections'],
            'atlas_false_positives': group_result['atlas_detections'],
            'pfasgroups_fpr': (group_result['pfasgroups_detections'] / max(group_result['molecules_tested'], 1)) * 100,
            'atlas_fpr': (group_result['atlas_detections'] / max(group_result['molecules_tested'], 1)) * 100,
        }
        analysis['group_details'].append(group_analysis)
    
    return analysis

def create_nonfluorinated_html_report(results, analysis, timestamp):
    """Create HTML report for non-fluorinated analysis"""
    
    html_content = f"""
<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Non-Fluorinated PFAS Exclusion Analysis</title>
    <style>
        body {{
            font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif;
            margin: 0;
            padding: 20px;
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            min-height: 100vh;
            color: #333;
        }}
        .container {{
            max-width: 1200px;
            margin: 0 auto;
            background-color: white;
            padding: 40px;
            border-radius: 15px;
            box-shadow: 0 10px 30px rgba(0,0,0,0.3);
        }}
        h1 {{
            color: #2c5aa0;
            text-align: center;
            font-size: 2.8em;
            margin-bottom: 10px;
        }}
        h2 {{
            color: #34495e;
            border-left: 4px solid #667eea;
            padding-left: 15px;
            margin-top: 50px;
        }}
        .summary-grid {{
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(200px, 1fr));
            gap: 20px;
            margin: 30px 0;
        }}
        .summary-card {{
            background: linear-gradient(45deg, #667eea 0%, #764ba2 100%);
            color: white;
            padding: 25px;
            border-radius: 12px;
            text-align: center;
            box-shadow: 0 5px 15px rgba(0,0,0,0.2);
        }}
        .summary-number {{
            font-size: 2.5em;
            font-weight: bold;
            margin-bottom: 10px;
        }}
        .summary-label {{
            font-size: 1.0em;
            opacity: 0.95;
        }}
        .results-table {{
            width: 100%;
            border-collapse: collapse;
            margin: 30px 0;
            background-color: white;
            border-radius: 8px;
            overflow: hidden;
            box-shadow: 0 5px 15px rgba(0,0,0,0.1);
        }}
        .results-table th, .results-table td {{
            padding: 15px;
            text-align: left;
            border-bottom: 1px solid #ddd;
        }}
        .results-table th {{
            background: linear-gradient(45deg, #667eea 0%, #764ba2 100%);
            color: white;
            font-weight: bold;
        }}
        .results-table tr:hover {{
            background-color: #f5f5f5;
        }}
        .excellent {{ background-color: #d4edda; color: #155724; font-weight: bold; }}
        .good {{ background-color: #d1ecf1; color: #0c5460; font-weight: bold; }}
        .warning {{ background-color: #fff3cd; color: #856404; font-weight: bold; }}
        .danger {{ background-color: #f8d7da; color: #721c24; font-weight: bold; }}
        .highlight-section {{
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            color: white;
            padding: 40px;
            border-radius: 15px;
            margin: 50px 0;
            text-align: center;
        }}
        .interpretation {{
            background: linear-gradient(135deg, #28a745 0%, #20c997 100%);
            color: white;
            padding: 30px;
            border-radius: 12px;
            margin: 40px 0;
        }}
    </style>
</head>
<body>
    <div class="container">
        <h1>🚫 Non-Fluorinated PFAS Exclusion Test</h1>
        
        <div class="highlight-section">
            <h3>🎯 Testing PFAS Detection System Specificity</h3>
            <p style="font-size: 1.2em; margin: 20px 0;">
                Evaluating whether PFASGroups and PFAS-Atlas correctly exclude non-fluorinated molecules<br>
                Analysis Date: {datetime.now().strftime('%B %d, %Y at %H:%M')}
            </p>
            <p><strong>Expected Result:</strong> 0% detection (these molecules should NOT be classified as PFAS)</p>
        </div>

        <div class="summary-grid">
            <div class="summary-card">
                <div class="summary-number">{analysis['total_molecules']}</div>
                <div class="summary-label">Non-Fluorinated<br>Molecules Tested</div>
            </div>
            <div class="summary-card">
                <div class="summary-number">{analysis['total_functional_groups']}</div>
                <div class="summary-label">Functional Group<br>Types Tested</div>
            </div>
            <div class="summary-card">
                <div class="summary-number">{analysis['pfasgroups_fpr']:.1f}%</div>
                <div class="summary-label">PFASGroups<br>False Positive Rate</div>
            </div>
            <div class="summary-card">
                <div class="summary-number">{analysis['atlas_fpr']:.1f}%</div>
                <div class="summary-label">PFAS-Atlas<br>False Positive Rate</div>
            </div>
        </div>

        <h2>📊 Detailed Results by Functional Group</h2>
        <table class="results-table">
            <thead>
                <tr>
                    <th>Functional Group</th>
                    <th>Molecules Tested</th>
                    <th>PFASGroups False Positives</th>
                    <th>PFASGroups FPR</th>
                    <th>Atlas False Positives</th>
                    <th>Atlas FPR</th>
                    <th>Overall Assessment</th>
                </tr>
            </thead>
            <tbody>"""
    
    # Add group details
    for group in analysis['group_details']:
        # Determine assessment class
        max_fpr = max(group['pfasgroups_fpr'], group['atlas_fpr'])
        if max_fpr == 0:
            assessment_class = 'excellent'
            assessment_text = '✅ Perfect'
        elif max_fpr < 5:
            assessment_class = 'good'
            assessment_text = '✅ Good'
        elif max_fpr < 15:
            assessment_class = 'warning'
            assessment_text = '⚠️ Acceptable'
        else:
            assessment_class = 'danger'
            assessment_text = '❌ Poor'
        
        html_content += f"""
                <tr>
                    <td><strong>Group {group['group_id']}</strong><br>{group['group_name']}</td>
                    <td>{group['molecules_tested']}</td>
                    <td>{group['pfasgroups_false_positives']}</td>
                    <td>{group['pfasgroups_fpr']:.1f}%</td>
                    <td>{group['atlas_false_positives']}</td>
                    <td>{group['atlas_fpr']:.1f}%</td>
                    <td class="{assessment_class}">{assessment_text}</td>
                </tr>"""
    
    # Overall interpretation
    if analysis['pfasgroups_fpr'] == 0 and analysis['atlas_fpr'] == 0:
        interpretation_text = "🎉 <strong>PERFECT PERFORMANCE!</strong> Both systems correctly identified that all non-fluorinated molecules are not PFAS. This demonstrates excellent specificity."
        interpretation_class = "excellent"
    elif analysis['pfasgroups_fpr'] < 5 and analysis['atlas_fpr'] < 5:
        interpretation_text = "✅ <strong>EXCELLENT PERFORMANCE!</strong> Both systems have very low false positive rates (&lt;5%), showing good specificity for PFAS detection."
        interpretation_class = "good"
    elif analysis['pfasgroups_fpr'] < 15 and analysis['atlas_fpr'] < 15:
        interpretation_text = "⚠️ <strong>ACCEPTABLE PERFORMANCE.</strong> False positive rates are moderate (&lt;15%). Some calibration may be beneficial."
        interpretation_class = "warning"
    else:
        interpretation_text = "❌ <strong>CONCERNING PERFORMANCE.</strong> High false positive rates indicate potential issues with PFAS specificity. Review and calibration recommended."
        interpretation_class = "danger"
    
    html_content += f"""
            </tbody>
        </table>

        <div class="interpretation">
            <h3>🔍 Performance Interpretation</h3>
            <p style="font-size: 1.1em; margin: 0;">{interpretation_text}</p>
        </div>

        <h2>📋 Test Methodology</h2>
        <div style="background: #f8f9fa; padding: 30px; border-radius: 12px; margin: 30px 0;">
            <h3>Functional Groups Tested:</h3>
            <ul>"""
    
    for group in analysis['group_details']:
        html_content += f"<li><strong>Group {group['group_id']} ({group['group_name']}):</strong> {group['molecules_tested']} non-fluorinated molecules</li>"
    
    html_content += f"""
            </ul>
            <h3>Test Design:</h3>
            <ul>
                <li>Generated molecules with target functional groups but <strong>without fluorination</strong></li>
                <li>Used same molecular generation parameters as fluorinated tests</li>
                <li>Verified all test molecules contain zero fluorine atoms</li>
                <li>Expected result: 0% PFAS detection (perfect specificity)</li>
            </ul>
        </div>

        <h2>🎯 Key Findings</h2>
        <div style="background: #e3f2fd; padding: 25px; border-radius: 12px; margin: 30px 0;">
            <ul style="margin: 0; padding-left: 20px;">
                <li><strong>Total Specificity:</strong> {100 - analysis['pfasgroups_fpr']:.1f}% (PFASGroups), {100 - analysis['atlas_fpr']:.1f}% (PFAS-Atlas)</li>
                <li><strong>False Discovery Rate:</strong> {analysis['pfasgroups_false_positives']}/{analysis['total_molecules']} (PFASGroups), {analysis['atlas_false_positives']}/{analysis['total_molecules']} (Atlas)</li>
                <li><strong>System Reliability:</strong> Both systems {'demonstrate excellent' if max(analysis['pfasgroups_fpr'], analysis['atlas_fpr']) < 5 else 'show acceptable' if max(analysis['pfasgroups_fpr'], analysis['atlas_fpr']) < 15 else 'require attention for'} specificity for PFAS detection</li>
            </ul>
        </div>
    </div>
</body>
</html>"""
    
    return html_content

def main():
    """Main function to analyze non-fluorinated benchmark results"""
    
    if len(sys.argv) != 2:
        print("Usage: python analyze_nonfluorinated.py <nonfluorinated_results.json>")
        sys.exit(1)
    
    filename = sys.argv[1]
    
    print("🚫 NON-FLUORINATED PFAS EXCLUSION ANALYSIS")
    print("=" * 50)
    
    # Load and analyze results
    results = load_nonfluorinated_results(filename)
    analysis = analyze_nonfluorinated_results(results)
    
    # Print summary
    print(f"\n📊 Analysis Summary:")
    print(f"   • Total molecules tested: {analysis['total_molecules']}")
    print(f"   • Functional group types: {analysis['total_functional_groups']}")
    print(f"   • PFASGroups false positive rate: {analysis['pfasgroups_fpr']:.1f}%")
    print(f"   • PFAS-Atlas false positive rate: {analysis['atlas_fpr']:.1f}%")
    
    # Group-by-group summary
    print(f"\n📋 Results by Functional Group:")
    for group in analysis['group_details']:
        print(f"   • Group {group['group_id']} ({group['group_name']}): {group['molecules_tested']} molecules")
        print(f"     - PFASGroups: {group['pfasgroups_false_positives']} false positives ({group['pfasgroups_fpr']:.1f}%)")
        print(f"     - Atlas: {group['atlas_false_positives']} false positives ({group['atlas_fpr']:.1f}%)")
    
    # Generate HTML report
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    html_content = create_nonfluorinated_html_report(results, analysis, timestamp)
    html_filename = f"nonfluorinated_analysis_{timestamp}.html"
    
    with open(html_filename, 'w') as f:
        f.write(html_content)
    
    print(f"\n✅ Non-Fluorinated Analysis Complete!")
    print(f"📁 Generated Report: {html_filename}")
    
    # Performance interpretation
    if analysis['pfasgroups_fpr'] == 0 and analysis['atlas_fpr'] == 0:
        print(f"🎉 PERFECT: Both systems correctly excluded all non-fluorinated molecules!")
    elif analysis['pfasgroups_fpr'] < 5 and analysis['atlas_fpr'] < 5:
        print(f"✅ EXCELLENT: Both systems show very low false positive rates (<5%)")
    elif analysis['pfasgroups_fpr'] < 15 and analysis['atlas_fpr'] < 15:
        print(f"⚠️  ACCEPTABLE: Moderate false positive rates (<15%) - consider calibration")
    else:
        print(f"❌ WARNING: High false positive rates detected - systems need review")

if __name__ == "__main__":
    main()