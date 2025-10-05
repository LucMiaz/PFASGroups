"""
Comprehensive PFAS Group Algorithm Analysis Report Generator

This script generates detailed analysis of the PFAS group identification algorithm's
specificity and accuracy using the test framework from test_examples.py.
"""

import sys
import os
import json
import pandas as pd
import numpy as np
from datetime import datetime
import matplotlib.pyplot as plt
import seaborn as sns
from collections import defaultdict, Counter

# Add the PFASgroups directory to the path
sys.path.append(os.path.dirname(os.path.abspath(__file__)))

# Import required modules
try:
    from PFASgroups.core import parse_PFAS_groups
    from PFASgroups.tests.test_examples import (
        TestPFASGroups, create_specificity_test_molecules,
        df_test_pfas_group_specificity, OECD_PFAS_GROUPS, 
        GENERIC_PFAS_GROUPS, IGNORE_GROUPS, load_equivalent_groups_from_json
    )
    from rdkit import Chem
    from rdkit.Chem.rdMolDescriptors import CalcMolFormula
except ImportError as e:
    print(f"Import error: {e}")
    print("Please ensure PFASgroups module is properly installed and in the Python path")
    sys.exit(1)

class PFASAlgorithmAnalyzer:
    """
    Comprehensive analyzer for PFAS group identification algorithm performance.
    """
    
    def __init__(self):
        self.test_results = {
            'oecd_results': [],
            'generic_results': [],
            'specificity_results': None,
            'oecd_only_results': None,
            'generic_only_results': None
        }
        self.group_info = self._load_group_info()
        self.report_data = {}
        
    def _load_group_info(self):
        """Load PFAS group information for analysis."""
        group_info = {}
        
        # Add OECD groups
        for group_id, group_name, template, pathtype in OECD_PFAS_GROUPS:
            group_info[group_id] = {
                'name': group_name,
                'type': 'OECD',
                'pathtype': pathtype,
                'template': template
            }
        
        # Add generic groups
        for group_id, group_name, template, insertion_mode in GENERIC_PFAS_GROUPS:
            group_info[group_id] = {
                'name': group_name,
                'type': 'Generic',
                'pathtype': 'Both',
                'template': template,
                'insertion_mode': insertion_mode
            }
        
        return group_info
    
    def run_comprehensive_tests(self):
        """Run all tests and collect comprehensive data."""
        print("Running comprehensive PFAS algorithm analysis...")
        print("=" * 60)
        
        # 1. Run standard specificity tests
        print("1. Running standard specificity tests...")
        self.test_results['specificity_results'] = df_test_pfas_group_specificity(verbose=False)
        
        # 2. Run OECD-only tests
        print("2. Running OECD-only group tests...")
        self.test_results['oecd_only_results'] = self._run_oecd_only_tests()
        
        # 3. Run Generic-only tests
        print("3. Running Generic-only group tests...")
        self.test_results['generic_only_results'] = self._run_generic_only_tests()
        
        # 4. Generate additional test molecules for misidentification analysis
        print("4. Generating additional test data...")
        self._generate_additional_test_data()
        
        print("All tests completed!")
        
    def _run_oecd_only_tests(self):
        """Run specificity tests using only OECD PFAS groups."""
        print("  Creating test molecules...")
        test_molecules = create_specificity_test_molecules()
        
        # Filter to only include OECD groups
        oecd_test_molecules = test_molecules[
            test_molecules['group_ids'].apply(
                lambda x: any(gid in self.group_info and self.group_info[gid]['type'] == 'OECD' for gid in x)
            )
        ].copy()
        
        print(f"  Testing {len(oecd_test_molecules)} OECD test molecules...")
        results = []
        
        for i, (inchi, inchikey, formula, pathtype, group_ids, smiles) in oecd_test_molecules.iterrows():
            try:
                mol = Chem.MolFromSmiles(smiles)
                if mol is None:
                    continue
                
                # Get OECD groups only
                oecd_groups = [g for g in self._get_oecd_pfas_groups() if g.id not in IGNORE_GROUPS]
                
                # Test with OECD groups only
                formula = CalcMolFormula(mol)
                matches = parse_PFAS_groups(mol, formula, pfas_groups=oecd_groups)
                
                detected_groups = [match[0].id for match in matches]
                expected_oecd_groups = [gid for gid in group_ids if gid in self.group_info and self.group_info[gid]['type'] == 'OECD']
                
                expected_group_detected = len(set(expected_oecd_groups).intersection(detected_groups)) == len(expected_oecd_groups)
                
                results.append({
                    'group_ids': group_ids,
                    'expected_oecd_groups': expected_oecd_groups,
                    'smiles': smiles,
                    'inchi': inchi,
                    'valid_smiles': True,
                    'expected_group_detected': expected_group_detected,
                    'detected_groups': detected_groups,
                    'n_detected_groups': len(detected_groups),
                    'is_specific': len(detected_groups) <= 2,  # Allow some overlap
                    'error': None
                })
                
            except Exception as e:
                results.append({
                    'group_ids': group_ids,
                    'expected_oecd_groups': [],
                    'smiles': smiles,
                    'inchi': inchi,
                    'valid_smiles': False,
                    'expected_group_detected': False,
                    'detected_groups': [],
                    'n_detected_groups': 0,
                    'is_specific': False,
                    'error': str(e)
                })
        
        return pd.DataFrame(results)
    
    def _run_generic_only_tests(self):
        """Run specificity tests using only Generic PFAS groups."""
        print("  Creating test molecules...")
        test_molecules = create_specificity_test_molecules()
        
        # Filter to only include Generic groups
        generic_test_molecules = test_molecules[
            test_molecules['group_ids'].apply(
                lambda x: any(gid in self.group_info and self.group_info[gid]['type'] == 'Generic' for gid in x)
            )
        ].copy()
        
        print(f"  Testing {len(generic_test_molecules)} Generic test molecules...")
        results = []
        
        for i, (inchi, inchikey, formula, pathtype, group_ids, smiles) in generic_test_molecules.iterrows():
            try:
                mol = Chem.MolFromSmiles(smiles)
                if mol is None:
                    continue
                
                # Get Generic groups only
                generic_groups = [g for g in self._get_generic_pfas_groups() if g.id not in IGNORE_GROUPS]
                
                # Test with Generic groups only
                formula = CalcMolFormula(mol)
                matches = parse_PFAS_groups(mol, formula, pfas_groups=generic_groups)
                
                detected_groups = [match[0].id for match in matches]
                expected_generic_groups = [gid for gid in group_ids if gid in self.group_info and self.group_info[gid]['type'] == 'Generic']
                
                expected_group_detected = len(set(expected_generic_groups).intersection(detected_groups)) == len(expected_generic_groups)
                
                results.append({
                    'group_ids': group_ids,
                    'expected_generic_groups': expected_generic_groups,
                    'smiles': smiles,
                    'inchi': inchi,
                    'valid_smiles': True,
                    'expected_group_detected': expected_group_detected,
                    'detected_groups': detected_groups,
                    'n_detected_groups': len(detected_groups),
                    'is_specific': len(detected_groups) <= 2,  # Allow some overlap
                    'error': None
                })
                
            except Exception as e:
                results.append({
                    'group_ids': group_ids,
                    'expected_generic_groups': [],
                    'smiles': smiles,
                    'inchi': inchi,
                    'valid_smiles': False,
                    'expected_group_detected': False,
                    'detected_groups': [],
                    'n_detected_groups': 0,
                    'is_specific': False,
                    'error': str(e)
                })
        
        return pd.DataFrame(results)
    
    def _get_oecd_pfas_groups(self):
        """Get OECD PFAS groups for testing."""
        try:
            # Try to load from the standard location
            from PFASgroups.core import load_PFASGroups
            @load_PFASGroups()
            def get_all_groups(pfas_groups=None):
                return pfas_groups
            
            all_groups = get_all_groups()
            oecd_groups = [g for g in all_groups if g.id <= 28 and g.id not in IGNORE_GROUPS]
            return oecd_groups
        except:
            # Fallback: create minimal group objects for testing
            from PFASgroups.PFASGroupModel import PFASGroup
            groups = []
            for group_id, group_name, template, pathtype in OECD_PFAS_GROUPS:
                if group_id not in IGNORE_GROUPS:
                    groups.append(PFASGroup(id=group_id, name=group_name))
            return groups
    
    def _get_generic_pfas_groups(self):
        """Get Generic PFAS groups for testing."""
        try:
            # Try to load from the standard location
            from PFASgroups.core import load_PFASGroups
            @load_PFASGroups()
            def get_all_groups(pfas_groups=None):
                return pfas_groups
            
            all_groups = get_all_groups()
            generic_groups = [g for g in all_groups if g.id >= 29 and g.id not in IGNORE_GROUPS]
            return generic_groups
        except:
            # Fallback: create minimal group objects for testing
            from PFASgroups.PFASGroupModel import PFASGroup
            groups = []
            for group_id, group_name, template, insertion_mode in GENERIC_PFAS_GROUPS:
                if group_id not in IGNORE_GROUPS:
                    groups.append(PFASGroup(id=group_id, name=group_name))
            return groups
    
    def _generate_additional_test_data(self):
        """Generate additional test data for misidentification analysis."""
        # This could include edge cases, borderline molecules, etc.
        # For now, we'll use the existing test framework
        pass
    
    def analyze_misidentifications(self):
        """Analyze all misidentified molecules in detail."""
        print("Analyzing misidentifications...")
        
        misidentifications = {
            'all_groups': [],
            'oecd_only': [],
            'generic_only': []
        }
        
        # Analyze standard specificity results
        if self.test_results['specificity_results'] is not None:
            valid_tests = self.test_results['specificity_results'][
                self.test_results['specificity_results']['valid_smiles'] == True
            ]
            
            # False negatives (expected groups not detected)
            false_negatives = valid_tests[valid_tests['expected_group_detected'] == False]
            for _, row in false_negatives.iterrows():
                misidentifications['all_groups'].append({
                    'type': 'False Negative',
                    'expected_groups': row['group_ids'],
                    'detected_groups': row['detected_groups'],
                    'smiles': row['smiles'],
                    'inchi': row['inchi'],
                    'n_detected': row['n_detected_groups']
                })
            
            # Low specificity (too many groups detected)
            low_specificity = valid_tests[valid_tests['n_detected_groups'] > 3]
            for _, row in low_specificity.iterrows():
                misidentifications['all_groups'].append({
                    'type': 'Low Specificity',
                    'expected_groups': row['group_ids'],
                    'detected_groups': row['detected_groups'],
                    'smiles': row['smiles'],
                    'inchi': row['inchi'],
                    'n_detected': row['n_detected_groups']
                })
        
        # Analyze OECD-only results
        if self.test_results['oecd_only_results'] is not None:
            valid_oecd = self.test_results['oecd_only_results'][
                self.test_results['oecd_only_results']['valid_smiles'] == True
            ]
            
            false_negatives_oecd = valid_oecd[valid_oecd['expected_group_detected'] == False]
            for _, row in false_negatives_oecd.iterrows():
                misidentifications['oecd_only'].append({
                    'type': 'False Negative',
                    'expected_groups': row['expected_oecd_groups'],
                    'detected_groups': row['detected_groups'],
                    'smiles': row['smiles'],
                    'inchi': row['inchi'],
                    'n_detected': row['n_detected_groups']
                })
            
            low_specificity_oecd = valid_oecd[valid_oecd['n_detected_groups'] > 2]
            for _, row in low_specificity_oecd.iterrows():
                misidentifications['oecd_only'].append({
                    'type': 'Low Specificity',
                    'expected_groups': row['expected_oecd_groups'],
                    'detected_groups': row['detected_groups'],
                    'smiles': row['smiles'],
                    'inchi': row['inchi'],
                    'n_detected': row['n_detected_groups']
                })
        
        # Analyze Generic-only results
        if self.test_results['generic_only_results'] is not None:
            valid_generic = self.test_results['generic_only_results'][
                self.test_results['generic_only_results']['valid_smiles'] == True
            ]
            
            false_negatives_generic = valid_generic[valid_generic['expected_group_detected'] == False]
            for _, row in false_negatives_generic.iterrows():
                misidentifications['generic_only'].append({
                    'type': 'False Negative',
                    'expected_groups': row['expected_generic_groups'],
                    'detected_groups': row['detected_groups'],
                    'smiles': row['smiles'],
                    'inchi': row['inchi'],
                    'n_detected': row['n_detected_groups']
                })
            
            low_specificity_generic = valid_generic[valid_generic['n_detected_groups'] > 2]
            for _, row in low_specificity_generic.iterrows():
                misidentifications['generic_only'].append({
                    'type': 'Low Specificity',
                    'expected_groups': row['expected_generic_groups'],
                    'detected_groups': row['detected_groups'],
                    'smiles': row['smiles'],
                    'inchi': row['inchi'],
                    'n_detected': row['n_detected_groups']
                })
        
        self.report_data['misidentifications'] = misidentifications
        return misidentifications
    
    def calculate_performance_metrics(self):
        """Calculate comprehensive performance metrics."""
        print("Calculating performance metrics...")
        
        metrics = {
            'all_groups': {},
            'oecd_only': {},
            'generic_only': {}
        }
        
        # All groups metrics
        if self.test_results['specificity_results'] is not None:
            valid_tests = self.test_results['specificity_results'][
                self.test_results['specificity_results']['valid_smiles'] == True
            ]
            
            if len(valid_tests) > 0:
                metrics['all_groups'] = {
                    'total_tests': len(valid_tests),
                    'detection_rate': valid_tests['expected_group_detected'].mean(),
                    'specificity_rate': valid_tests['is_specific'].mean(),
                    'avg_groups_detected': valid_tests['n_detected_groups'].mean(),
                    'false_negative_rate': 1 - valid_tests['expected_group_detected'].mean(),
                    'low_specificity_rate': (valid_tests['n_detected_groups'] > 3).mean()
                }
        
        # OECD-only metrics
        if self.test_results['oecd_only_results'] is not None:
            valid_oecd = self.test_results['oecd_only_results'][
                self.test_results['oecd_only_results']['valid_smiles'] == True
            ]
            
            if len(valid_oecd) > 0:
                metrics['oecd_only'] = {
                    'total_tests': len(valid_oecd),
                    'detection_rate': valid_oecd['expected_group_detected'].mean(),
                    'specificity_rate': valid_oecd['is_specific'].mean(),
                    'avg_groups_detected': valid_oecd['n_detected_groups'].mean(),
                    'false_negative_rate': 1 - valid_oecd['expected_group_detected'].mean(),
                    'low_specificity_rate': (valid_oecd['n_detected_groups'] > 2).mean()
                }
        
        # Generic-only metrics
        if self.test_results['generic_only_results'] is not None:
            valid_generic = self.test_results['generic_only_results'][
                self.test_results['generic_only_results']['valid_smiles'] == True
            ]
            
            if len(valid_generic) > 0:
                metrics['generic_only'] = {
                    'total_tests': len(valid_generic),
                    'detection_rate': valid_generic['expected_group_detected'].mean(),
                    'specificity_rate': valid_generic['is_specific'].mean(),
                    'avg_groups_detected': valid_generic['n_detected_groups'].mean(),
                    'false_negative_rate': 1 - valid_generic['expected_group_detected'].mean(),
                    'low_specificity_rate': (valid_generic['n_detected_groups'] > 2).mean()
                }
        
        self.report_data['metrics'] = metrics
        return metrics
    
    def analyze_group_performance(self):
        """Analyze performance of individual PFAS groups."""
        print("Analyzing individual group performance...")
        
        group_performance = {}
        
        # Analyze all groups performance
        if self.test_results['specificity_results'] is not None:
            valid_tests = self.test_results['specificity_results'][
                self.test_results['specificity_results']['valid_smiles'] == True
            ]
            
            for _, row in valid_tests.iterrows():
                for group_id in row['group_ids']:
                    if group_id not in group_performance:
                        group_performance[group_id] = {
                            'total_tests': 0,
                            'successful_detections': 0,
                            'false_negatives': 0,
                            'group_name': self.group_info.get(group_id, {}).get('name', f'Unknown {group_id}'),
                            'group_type': self.group_info.get(group_id, {}).get('type', 'Unknown'),
                            'avg_other_groups_detected': []
                        }
                    
                    group_performance[group_id]['total_tests'] += 1
                    
                    if group_id in row['detected_groups']:
                        group_performance[group_id]['successful_detections'] += 1
                    else:
                        group_performance[group_id]['false_negatives'] += 1
                    
                    # Track how many other groups are typically detected
                    other_groups = len([g for g in row['detected_groups'] if g != group_id])
                    group_performance[group_id]['avg_other_groups_detected'].append(other_groups)
            
            # Calculate averages
            for group_id in group_performance:
                if group_performance[group_id]['avg_other_groups_detected']:
                    group_performance[group_id]['avg_other_groups_detected'] = np.mean(
                        group_performance[group_id]['avg_other_groups_detected']
                    )
                else:
                    group_performance[group_id]['avg_other_groups_detected'] = 0
                
                group_performance[group_id]['detection_rate'] = (
                    group_performance[group_id]['successful_detections'] / 
                    group_performance[group_id]['total_tests']
                )
        
        self.report_data['group_performance'] = group_performance
        return group_performance
    
    def generate_report(self, output_file='pfas_algorithm_analysis_report.html'):
        """Generate comprehensive HTML report."""
        print(f"Generating comprehensive report: {output_file}")
        
        # Run all analyses
        self.run_comprehensive_tests()
        misidentifications = self.analyze_misidentifications()
        metrics = self.calculate_performance_metrics()
        group_performance = self.analyze_group_performance()
        
        # Generate HTML report
        html_content = self._generate_html_report(misidentifications, metrics, group_performance)
        
        with open(output_file, 'w', encoding='utf-8') as f:
            f.write(html_content)
        
        print(f"Report saved to: {output_file}")
        
        # Also save raw data
        self._save_raw_data()
        
        return output_file
    
    def _generate_html_report(self, misidentifications, metrics, group_performance):
        """Generate the HTML content for the report."""
        
        html = f"""
<!DOCTYPE html>
<html>
<head>
    <title>PFAS Algorithm Performance Analysis Report</title>
    <style>
        body {{ font-family: Arial, sans-serif; margin: 40px; }}
        .header {{ text-align: center; margin-bottom: 30px; }}
        .section {{ margin-bottom: 40px; }}
        .metric-grid {{ display: grid; grid-template-columns: repeat(auto-fit, minmax(250px, 1fr)); gap: 20px; }}
        .metric-card {{ border: 1px solid #ddd; padding: 15px; border-radius: 8px; background: #f9f9f9; }}
        .metric-title {{ font-weight: bold; margin-bottom: 10px; color: #333; }}
        .metric-value {{ font-size: 24px; font-weight: bold; }}
        .good {{ color: #28a745; }}
        .warning {{ color: #ffc107; }}
        .bad {{ color: #dc3545; }}
        table {{ width: 100%; border-collapse: collapse; margin-top: 10px; }}
        th, td {{ border: 1px solid #ddd; padding: 8px; text-align: left; }}
        th {{ background-color: #f2f2f2; }}
        .smiles {{ font-family: monospace; font-size: 12px; }}
        .group-list {{ font-size: 12px; }}
        .summary-stats {{ background: #e9ecef; padding: 15px; border-radius: 8px; }}
    </style>
</head>
<body>
    <div class="header">
        <h1>PFAS Group Identification Algorithm</h1>
        <h2>Performance Analysis Report</h2>
        <p>Generated on: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}</p>
    </div>
    
    <div class="section">
        <h2>Executive Summary</h2>
        <div class="summary-stats">
            {self._generate_executive_summary(metrics)}
        </div>
    </div>
    
    <div class="section">
        <h2>Performance Metrics Comparison</h2>
        <div class="metric-grid">
            {self._generate_metrics_cards(metrics)}
        </div>
    </div>
    
    <div class="section">
        <h2>Individual Group Performance</h2>
        {self._generate_group_performance_table(group_performance)}
    </div>
    
    <div class="section">
        <h2>Misidentification Analysis</h2>
        {self._generate_misidentification_analysis(misidentifications)}
    </div>
    
    <div class="section">
        <h2>Low Specificity Examples</h2>
        {self._generate_low_specificity_examples(misidentifications)}
    </div>
    
    <div class="section">
        <h2>Recommendations</h2>
        {self._generate_recommendations(metrics, group_performance)}
    </div>
    
</body>
</html>
        """
        
        return html
    
    def _generate_executive_summary(self, metrics):
        """Generate executive summary section."""
        all_metrics = metrics.get('all_groups', {})
        oecd_metrics = metrics.get('oecd_only', {})
        generic_metrics = metrics.get('generic_only', {})
        
        summary = f"""
        <h3>Key Findings</h3>
        <ul>
            <li><strong>Overall Detection Rate:</strong> {all_metrics.get('detection_rate', 0):.1%} 
                ({all_metrics.get('total_tests', 0)} tests)</li>
            <li><strong>Overall Specificity Rate:</strong> {all_metrics.get('specificity_rate', 0):.1%}</li>
            <li><strong>OECD Groups Performance:</strong> {oecd_metrics.get('detection_rate', 0):.1%} detection, 
                {oecd_metrics.get('specificity_rate', 0):.1%} specificity</li>
            <li><strong>Generic Groups Performance:</strong> {generic_metrics.get('detection_rate', 0):.1%} detection, 
                {generic_metrics.get('specificity_rate', 0):.1%} specificity</li>
            <li><strong>Average Groups per Test:</strong> {all_metrics.get('avg_groups_detected', 0):.1f}</li>
        </ul>
        """
        
        return summary
    
    def _generate_metrics_cards(self, metrics):
        """Generate metrics comparison cards."""
        cards = ""
        
        for test_type, metric_data in metrics.items():
            if not metric_data:
                continue
                
            title = test_type.replace('_', ' ').title()
            detection_class = self._get_metric_class(metric_data.get('detection_rate', 0), 0.8, 0.6)
            specificity_class = self._get_metric_class(metric_data.get('specificity_rate', 0), 0.7, 0.5)
            
            cards += f"""
            <div class="metric-card">
                <div class="metric-title">{title}</div>
                <div>Detection Rate: <span class="metric-value {detection_class}">
                    {metric_data.get('detection_rate', 0):.1%}</span></div>
                <div>Specificity Rate: <span class="metric-value {specificity_class}">
                    {metric_data.get('specificity_rate', 0):.1%}</span></div>
                <div>Total Tests: {metric_data.get('total_tests', 0)}</div>
                <div>Avg Groups: {metric_data.get('avg_groups_detected', 0):.1f}</div>
            </div>
            """
        
        return cards
    
    def _get_metric_class(self, value, good_threshold, warning_threshold):
        """Get CSS class based on metric value."""
        if value >= good_threshold:
            return "good"
        elif value >= warning_threshold:
            return "warning"
        else:
            return "bad"
    
    def _generate_group_performance_table(self, group_performance):
        """Generate individual group performance table."""
        if not group_performance:
            return "<p>No group performance data available.</p>"
        
        # Sort by detection rate
        sorted_groups = sorted(group_performance.items(), 
                             key=lambda x: x[1]['detection_rate'], reverse=True)
        
        table = """
        <table>
            <thead>
                <tr>
                    <th>Group ID</th>
                    <th>Group Name</th>
                    <th>Type</th>
                    <th>Tests</th>
                    <th>Detection Rate</th>
                    <th>False Negatives</th>
                    <th>Avg Other Groups</th>
                </tr>
            </thead>
            <tbody>
        """
        
        for group_id, perf in sorted_groups:
            detection_class = self._get_metric_class(perf['detection_rate'], 0.8, 0.6)
            
            table += f"""
                <tr>
                    <td>{group_id}</td>
                    <td>{perf['group_name']}</td>
                    <td>{perf['group_type']}</td>
                    <td>{perf['total_tests']}</td>
                    <td class="{detection_class}">{perf['detection_rate']:.1%}</td>
                    <td>{perf['false_negatives']}</td>
                    <td>{perf['avg_other_groups_detected']:.1f}</td>
                </tr>
            """
        
        table += "</tbody></table>"
        return table
    
    def _generate_misidentification_analysis(self, misidentifications):
        """Generate misidentification analysis section."""
        analysis = ""
        
        for test_type, mistakes in misidentifications.items():
            if not mistakes:
                continue
                
            title = test_type.replace('_', ' ').title()
            false_negatives = [m for m in mistakes if m['type'] == 'False Negative']
            
            analysis += f"""
            <h3>{title} - False Negatives ({len(false_negatives)} cases)</h3>
            """
            
            if false_negatives:
                analysis += """
                <table>
                    <thead>
                        <tr>
                            <th>Expected Groups</th>
                            <th>Detected Groups</th>
                            <th>SMILES</th>
                            <th>InChI</th>
                        </tr>
                    </thead>
                    <tbody>
                """
                
                for mistake in false_negatives[:10]:  # Show first 10
                    expected_names = [self.group_info.get(g, {}).get('name', f'Group {g}') for g in mistake['expected_groups']]
                    detected_names = [self.group_info.get(g, {}).get('name', f'Group {g}') for g in mistake['detected_groups']]
                    
                    analysis += f"""
                        <tr>
                            <td class="group-list">{', '.join(expected_names)}</td>
                            <td class="group-list">{', '.join(detected_names) if detected_names else 'None'}</td>
                            <td class="smiles">{mistake['smiles'][:50]}...</td>
                            <td class="smiles">{mistake['inchi'][:50]}...</td>
                        </tr>
                    """
                
                analysis += "</tbody></table>"
                
                if len(false_negatives) > 10:
                    analysis += f"<p><em>Showing first 10 of {len(false_negatives)} false negatives.</em></p>"
        
        return analysis
    
    def _generate_low_specificity_examples(self, misidentifications):
        """Generate low specificity examples section."""
        analysis = ""
        
        for test_type, mistakes in misidentifications.items():
            if not mistakes:
                continue
                
            title = test_type.replace('_', ' ').title()
            low_specificity = [m for m in mistakes if m['type'] == 'Low Specificity']
            
            if low_specificity:
                analysis += f"""
                <h3>{title} - Low Specificity ({len(low_specificity)} cases)</h3>
                <table>
                    <thead>
                        <tr>
                            <th>Expected Groups</th>
                            <th>All Detected Groups</th>
                            <th>Count</th>
                            <th>SMILES</th>
                        </tr>
                    </thead>
                    <tbody>
                """
                
                # Sort by number of detected groups (worst first)
                low_specificity.sort(key=lambda x: x['n_detected'], reverse=True)
                
                for mistake in low_specificity[:10]:  # Show first 10
                    expected_names = [self.group_info.get(g, {}).get('name', f'Group {g}') for g in mistake['expected_groups']]
                    detected_names = [self.group_info.get(g, {}).get('name', f'Group {g}') for g in mistake['detected_groups']]
                    
                    analysis += f"""
                        <tr>
                            <td class="group-list">{', '.join(expected_names)}</td>
                            <td class="group-list">{', '.join(detected_names)}</td>
                            <td>{mistake['n_detected']}</td>
                            <td class="smiles">{mistake['smiles'][:50]}...</td>
                        </tr>
                    """
                
                analysis += "</tbody></table>"
                
                if len(low_specificity) > 10:
                    analysis += f"<p><em>Showing first 10 of {len(low_specificity)} low specificity cases.</em></p>"
        
        return analysis
    
    def _generate_recommendations(self, metrics, group_performance):
        """Generate recommendations based on analysis."""
        recommendations = "<h3>Algorithm Improvement Recommendations</h3><ul>"
        
        all_metrics = metrics.get('all_groups', {})
        
        # Detection rate recommendations
        if all_metrics.get('detection_rate', 0) < 0.8:
            recommendations += "<li><strong>Improve Detection Rate:</strong> Current detection rate is below 80%. Consider reviewing SMARTS patterns for frequently missed groups.</li>"
        
        # Specificity recommendations
        if all_metrics.get('specificity_rate', 0) < 0.7:
            recommendations += "<li><strong>Improve Specificity:</strong> High false positive rate detected. Consider making SMARTS patterns more specific or adjusting group hierarchy.</li>"
        
        # Group-specific recommendations
        if group_performance:
            worst_performers = sorted(group_performance.items(), 
                                    key=lambda x: x[1]['detection_rate'])[:3]
            
            if worst_performers and worst_performers[0][1]['detection_rate'] < 0.5:
                recommendations += f"<li><strong>Priority Groups for Improvement:</strong> "
                group_names = [self.group_info.get(g[0], {}).get('name', f'Group {g[0]}') for g in worst_performers]
                recommendations += f"{', '.join(group_names)} show consistently low detection rates.</li>"
        
        # OECD vs Generic comparison
        oecd_metrics = metrics.get('oecd_only', {})
        generic_metrics = metrics.get('generic_only', {})
        
        if oecd_metrics and generic_metrics:
            if oecd_metrics.get('detection_rate', 0) > generic_metrics.get('detection_rate', 0) + 0.1:
                recommendations += "<li><strong>Generic Group Enhancement:</strong> OECD groups significantly outperform generic groups. Consider improving generic group SMARTS patterns.</li>"
            elif generic_metrics.get('detection_rate', 0) > oecd_metrics.get('detection_rate', 0) + 0.1:
                recommendations += "<li><strong>OECD Group Enhancement:</strong> Generic groups outperform OECD groups. Review OECD-specific patterns for improvements.</li>"
        
        recommendations += "</ul>"
        return recommendations
    
    def _save_raw_data(self):
        """Save raw analysis data to files."""
        timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
        
        # Save test results
        if self.test_results['specificity_results'] is not None:
            self.test_results['specificity_results'].to_csv(f'all_groups_results_{timestamp}.csv', index=False)
        
        if self.test_results['oecd_only_results'] is not None:
            self.test_results['oecd_only_results'].to_csv(f'oecd_only_results_{timestamp}.csv', index=False)
        
        if self.test_results['generic_only_results'] is not None:
            self.test_results['generic_only_results'].to_csv(f'generic_only_results_{timestamp}.csv', index=False)
        
        # Save analysis data
        with open(f'analysis_report_data_{timestamp}.json', 'w') as f:
            # Convert numpy types to native Python types for JSON serialization
            def convert_numpy(obj):
                if isinstance(obj, np.integer):
                    return int(obj)
                elif isinstance(obj, np.floating):
                    return float(obj)
                elif isinstance(obj, np.ndarray):
                    return obj.tolist()
                else:
                    return obj
            
            json.dump(self.report_data, f, indent=2, default=convert_numpy)
        
        print(f"Raw data saved with timestamp: {timestamp}")

def main():
    """Main function to run the analysis."""
    print("PFAS Algorithm Performance Analysis")
    print("=" * 50)
    
    analyzer = PFASAlgorithmAnalyzer()
    
    try:
        report_file = analyzer.generate_report()
        print(f"\nAnalysis complete! Report available at: {report_file}")
        print("\nKey files generated:")
        print("- HTML report with detailed analysis")
        print("- CSV files with raw test results")
        print("- JSON file with analysis data")
        
    except Exception as e:
        print(f"Error during analysis: {e}")
        import traceback
        traceback.print_exc()

if __name__ == "__main__":
    main()