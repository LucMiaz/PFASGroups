#!/usr/bin/env python3
"""
Generate telomer validation report for review-app.

Reads telomer_validation_results.json and generates an HTML report
with statistics, group detection rates, and examples.
"""

import json
from pathlib import Path
from datetime import datetime

def generate_telomer_report():
    """Generate HTML report for telomer validation results."""
    
    # Load results
    results_file = Path('data/telomer_validation_results.json')
    if not results_file.exists():
        print(f"Error: {results_file} not found")
        return
    
    with open(results_file) as f:
        data = json.load(f)
    
    # Extract data
    test_date = data.get('test_date', 'Unknown')
    total = data.get('total_molecules', 0)
    detected = data.get('telomer_detected', 0)
    detection_rate = data.get('detection_rate', 0)
    group_counts = data.get('group_counts', [])
    results = data.get('results', [])
    
    # Calculate additional statistics
    missed = total - detected
    miss_rate = 100 - detection_rate
    
    # Get examples of detected and missed
    detected_examples = [r for r in results if r.get('detected')][:5]
    missed_examples = [r for r in results if not r.get('detected')][:5]
    
    # Generate HTML
    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Telomer Detection Validation Report</title>
    <style>
        body {{
            font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, Oxygen, Ubuntu, sans-serif;
            max-width: 1200px;
            margin: 0 auto;
            padding: 20px;
            background: #f5f5f5;
        }}
        .header {{
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            color: white;
            padding: 30px;
            border-radius: 10px;
            margin-bottom: 30px;
        }}
        .header h1 {{
            margin: 0 0 10px 0;
        }}
        .header p {{
            margin: 5px 0;
            opacity: 0.9;
        }}
        .stats-grid {{
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(200px, 1fr));
            gap: 20px;
            margin-bottom: 30px;
        }}
        .stat-card {{
            background: white;
            padding: 20px;
            border-radius: 8px;
            box-shadow: 0 2px 4px rgba(0,0,0,0.1);
        }}
        .stat-value {{
            font-size: 36px;
            font-weight: bold;
            color: #667eea;
            margin: 10px 0;
        }}
        .stat-label {{
            color: #666;
            font-size: 14px;
            text-transform: uppercase;
            letter-spacing: 0.5px;
        }}
        .section {{
            background: white;
            padding: 25px;
            border-radius: 8px;
            box-shadow: 0 2px 4px rgba(0,0,0,0.1);
            margin-bottom: 20px;
        }}
        h2 {{
            color: #333;
            margin-top: 0;
            border-bottom: 2px solid #667eea;
            padding-bottom: 10px;
        }}
        table {{
            width: 100%;
            border-collapse: collapse;
            margin-top: 15px;
        }}
        th, td {{
            text-align: left;
            padding: 12px;
            border-bottom: 1px solid #ddd;
        }}
        th {{
            background: #f8f9fa;
            font-weight: 600;
            color: #333;
        }}
        tr:hover {{
            background: #f8f9fa;
        }}
        .progress-bar {{
            height: 30px;
            background: #e0e0e0;
            border-radius: 15px;
            overflow: hidden;
            margin: 10px 0;
        }}
        .progress-fill {{
            height: 100%;
            background: linear-gradient(90deg, #667eea 0%, #764ba2 100%);
            display: flex;
            align-items: center;
            padding: 0 15px;
            color: white;
            font-weight: bold;
            font-size: 14px;
        }}
        .smiles {{
            font-family: 'Courier New', monospace;
            font-size: 12px;
            color: #666;
            word-break: break-all;
        }}
        .group-badge {{
            display: inline-block;
            background: #667eea;
            color: white;
            padding: 4px 8px;
            border-radius: 4px;
            font-size: 11px;
            margin: 2px;
        }}
        .success {{ color: #28a745; }}
        .warning {{ color: #ffc107; }}
        .info {{ color: #17a2b8; }}
    </style>
</head>
<body>
    <div class="header">
        <h1>\ud83e\uddea Telomer Detection Validation Report</h1>
        <p><strong>Dataset:</strong> PubChem Fluorotelomers</p>
        <p><strong>Test Date:</strong> {datetime.fromisoformat(test_date).strftime('%Y-%m-%d %H:%M:%S')}</p>
        <p><strong>Groups Tested:</strong> IDs 62-87 (including 13 new unsaturated and generic telomer groups)</p>
    </div>

    <div class="stats-grid">
        <div class="stat-card">
            <div class="stat-label">Total Molecules</div>
            <div class="stat-value">{total}</div>
        </div>
        <div class="stat-card">
            <div class="stat-label">Detected</div>
            <div class="stat-value success">{detected}</div>
        </div>
        <div class="stat-card">
            <div class="stat-label">Missed</div>
            <div class="stat-value warning">{missed}</div>
        </div>
        <div class="stat-card">
            <div class="stat-label">Detection Rate</div>
            <div class="stat-value info">{detection_rate:.1f}%</div>
        </div>
    </div>

    <div class="section">
        <h2>\ud83c\udfaf Overall Performance</h2>
        <p><strong>Detection Rate:</strong></p>
        <div class="progress-bar">
            <div class="progress-fill" style="width: {detection_rate}%">
                {detection_rate:.1f}% ({detected}/{total})
            </div>
        </div>
        <p style="margin-top: 15px;"><strong>Miss Rate:</strong></p>
        <div class="progress-bar">
            <div class="progress-fill" style="width: {miss_rate}%; background: linear-gradient(90deg, #ffc107 0%, #ff6347 100%)">
                {miss_rate:.1f}% ({missed}/{total})
            </div>
        </div>
    </div>

    <div class="section">
        <h2>\ud83d\udccb Telomer Groups Detected</h2>
        <p>Frequency of detection by PFAS group (most common first):</p>
        <table>
            <thead>
                <tr>
                    <th>Rank</th>
                    <th>ID</th>
                    <th>Group Name</th>
                    <th>Detections</th>
                    <th>Percentage</th>
                </tr>
            </thead>
            <tbody>"""
    
    for i, group in enumerate(group_counts, 1):
        percentage = (group['count'] / total) * 100
        html += f"""
                <tr>
                    <td>{i}</td>
                    <td><strong>{group['id']}</strong></td>
                    <td>{group['name']}</td>
                    <td>{group['count']}</td>
                    <td>{percentage:.1f}%</td>
                </tr>"""
    
    html += """
            </tbody>
        </table>
    </div>
    
    <div class="section">
        <h2>\u2705 Examples of Detected Telomers</h2>
        <p>Sample molecules successfully identified as fluorotelomers:</p>
        <table>
            <thead>
                <tr>
                    <th style="width: 60%">SMILES</th>
                    <th>Detected Groups</th>
                </tr>
            </thead>
            <tbody>"""
    
    for example in detected_examples:
        smiles = example.get('smiles', '')
        groups = example.get('telomer_groups', [])
        groups_html = ''.join([f'<span class="group-badge">{g["id"]}: {g["name"]}</span>' for g in groups])
        html += f"""
                <tr>
                    <td class="smiles">{smiles}</td>
                    <td>{groups_html}</td>
                </tr>"""
    
    html += """
            </tbody>
        </table>
    </div>"""
    
    if missed_examples:
        html += """
    <div class="section">
        <h2>\u26a0\ufe0f Examples of Missed Molecules</h2>
        <p>Sample molecules from PubChem fluorotelomer search that were not detected:</p>
        <table>
            <thead>
                <tr>
                    <th>SMILES</th>
                    <th>Reason</th>
                </tr>
            </thead>
            <tbody>"""
        
        for example in missed_examples:
            smiles = example.get('smiles', '')
            error = example.get('error', 'No telomer groups matched')
            html += f"""
                <tr>
                    <td class="smiles">{smiles}</td>
                    <td>{error}</td>
                </tr>"""
        
        html += """
            </tbody>
        </table>
    </div>"""
    
    html += """
    <div class="section">
        <h2>\ud83d\udcc8 Key Findings</h2>
        <ul>
            <li><strong>Generic telomer group (ID 86)</strong> provides broad coverage, catching most fluorotelomer structures</li>
            <li><strong>Specific functional groups</strong> (alcohols, ethoxylates, sulfonic acids) enable detailed classification</li>
            <li><strong>Unsaturated telomer groups (IDs 75-85)</strong> distinguish double-bond containing variants</li>
            <li><strong>Linker validation</strong> ensures only true telomers (with CH\u2082 linker chains) are identified</li>
        </ul>
    </div>

    <div class="section">
        <h2>\ud83d\udd27 Methodology</h2>
        <p><strong>Dataset Source:</strong> PubChem compounds from search term "fluorotelomer"</p>
        <p><strong>Groups Tested:</strong></p>
        <ul>
            <li><strong>Saturated telomers (62-74):</strong> Established fluorotelomer groups with specific functional groups</li>
            <li><strong>Unsaturated telomers (75-85):</strong> New groups featuring C=C double bonds (C(F)=C pattern)</li>
            <li><strong>Generic telomers (86-87):</strong> Catch-all groups for any fluorotelomer structure</li>
        </ul>
        <p><strong>Detection Method:</strong> SMARTS pattern matching with linker validation (CH\u2082 units between perfluorinated chain and functional group)</p>
    </div>

</body>
</html>"""
    
    # Save report
    output_file = Path('review-app/analysis_reports/telomer_validation.html')
    output_file.parent.mkdir(parents=True, exist_ok=True)
    
    with open(output_file, 'w') as f:
        f.write(html)
    
    print(f"\u2705 Telomer validation report generated: {output_file}")
    print(f"   Detection rate: {detection_rate:.1f}% ({detected}/{total})")
    print(f"   Groups detected: {len(group_counts)}")

if __name__ == '__main__':
    generate_telomer_report()
