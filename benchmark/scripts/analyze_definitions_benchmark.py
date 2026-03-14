"""
PFAS Definitions Benchmark Analysis and Reporting

Creates comprehensive analysis reports with visualizations:
- Performance metrics per definition
- Confusion matrices
- Agreement/disagreement patterns
- Interactive HTML report
"""

import json
import sys
import os
from datetime import datetime
from collections import defaultdict, Counter
import numpy as np
import pandas as pd

# Plotting libraries
import matplotlib.pyplot as plt
import seaborn as sns
import plotly.graph_objects as go
import plotly.express as px
from plotly.subplots import make_subplots

# Add parent directory to path
script_dir = os.path.dirname(os.path.abspath(__file__))
parent_dir = os.path.dirname(script_dir)
sys.path.append(parent_dir)


class DefinitionBenchmarkAnalyzer:
    """Analyze and visualize PFAS definitions benchmark results"""
    
    def __init__(self, benchmark_file: str):
        """Initialize with benchmark results file"""
        with open(benchmark_file, 'r') as f:
            self.results = json.load(f)
        
        self.definition_names = {
            1: "OECD Definition",
            2: "EU PFAS Restriction",
            3: "OPPT 2023",
            4: "UK PFAS Definition",
            5: "PFASTRUCTv5"
        }
        
        self.output_dir = os.path.dirname(benchmark_file)
        self.figures_dir = os.path.join(self.output_dir, 'figures_definitions')
        os.makedirs(self.figures_dir, exist_ok=True)
        
    def generate_summary_statistics(self) -> pd.DataFrame:
        """Generate summary statistics table"""
        print("\n Generating summary statistics...")
        
        stats_data = []
        
        for def_id in range(1, 6):
            def_name = self.definition_names[def_id]
            
            # True positives
            tp_count = self.results['benchmarks']['true_positives']['summary']['by_definition'][str(def_id)]
            tp_total = self.results['benchmarks']['true_positives']['summary']['total_molecules']
            tp_rate = (tp_count / max(tp_total, 1)) * 100
            
            # False positives (from true negatives)
            fp_count = self.results['benchmarks']['true_negatives']['summary']['by_definition'][str(def_id)]
            tn_total = self.results['benchmarks']['true_negatives']['summary']['total_molecules']
            fp_rate = (fp_count / max(tn_total, 1)) * 100
            specificity = 100 - fp_rate
            
            stats_data.append({
                'Definition': def_name,
                'Sensitivity (%)': f"{tp_rate:.1f}",
                'TP Count': f"{tp_count}/{tp_total}",
                'Specificity (%)': f"{specificity:.1f}",
                'FP Count': f"{fp_count}/{tn_total}",
            })
        
        df = pd.DataFrame(stats_data)
        
        print("\n Summary Statistics:")
        print(df.to_string(index=False))
        
        return df
    
    def plot_sensitivity_specificity(self):
        """Create sensitivity vs specificity plot"""
        print("\n Creating sensitivity vs specificity plot...")
        
        sensitivities = []
        specificities = []
        names = []
        
        for def_id in range(1, 6):
            # Sensitivity (true positive rate)
            tp_count = self.results['benchmarks']['true_positives']['summary']['by_definition'][str(def_id)]
            tp_total = self.results['benchmarks']['true_positives']['summary']['total_molecules']
            sensitivity = (tp_count / max(tp_total, 1)) * 100
            
            # Specificity (true negative rate)
            fp_count = self.results['benchmarks']['true_negatives']['summary']['by_definition'][str(def_id)]
            tn_total = self.results['benchmarks']['true_negatives']['summary']['total_molecules']
            specificity = 100 - (fp_count / max(tn_total, 1)) * 100
            
            sensitivities.append(sensitivity)
            specificities.append(specificity)
            names.append(self.definition_names[def_id])
        
        # Create scatter plot
        fig = go.Figure()
        
        fig.add_trace(go.Scatter(
            x=specificities,
            y=sensitivities,
            mode='markers+text',
            marker=dict(size=20, color=list(range(1, 6)), colorscale='Viridis'),
            text=names,
            textposition='top center',
            textfont=dict(size=10),
            hovertemplate='<b>%{text}</b><br>Specificity: %{x:.1f}%<br>Sensitivity: %{y:.1f}%<extra></extra>'
        ))
        
        # Add ideal point reference
        fig.add_trace(go.Scatter(
            x=[100], y=[100],
            mode='markers',
            marker=dict(size=15, color='red', symbol='star'),
            name='Ideal (100%, 100%)',
            hovertemplate='Ideal Point<extra></extra>'
        ))
        
        fig.update_layout(
            title='PFAS Definitions: Sensitivity vs Specificity',
            xaxis_title='Specificity (%) - Correctly Reject Non-PFAS',
            yaxis_title='Sensitivity (%) - Correctly Identify PFAS',
            xaxis=dict(range=[0, 105]),
            yaxis=dict(range=[0, 105]),
            showlegend=True,
            width=800,
            height=700,
        )
        
        # Add diagonal reference line (sensitivity = specificity)
        fig.add_shape(
            type="line",
            x0=0, y0=0, x1=100, y1=100,
            line=dict(color="gray", width=1, dash="dash"),
        )
        
        output_file = os.path.join(self.figures_dir, 'sensitivity_specificity.html')
        fig.write_html(output_file)
        print(f"   Saved to {output_file}")
        
        return fig
    
    def plot_detection_heatmap(self):
        """Create heatmap showing which definitions detect which compound categories"""
        print("\n Creating detection heatmap...")
        
        # Collect detection data
        categories = []
        detection_matrix = []
        
        # True positives
        tp_data = self.results['benchmarks']['true_positives']
        for category, results_list in tp_data['subcategories'].items():
            categories.append(category.replace('_', ' ').title())
            
            category_detections = [0] * 5
            for result in results_list:
                if result['valid']:
                    for def_id in result['detected_ids']:
                        category_detections[def_id - 1] += 1
            
            # Normalize by number of molecules
            total_mols = len([r for r in results_list if r['valid']])
            if total_mols > 0:
                category_detections = [d / total_mols * 100 for d in category_detections]
            
            detection_matrix.append(category_detections)
        
        # Create heatmap
        fig = go.Figure(data=go.Heatmap(
            z=detection_matrix,
            x=[self.definition_names[i] for i in range(1, 6)],
            y=categories,
            colorscale='RdYlGn',
            text=[[f"{val:.0f}%" for val in row] for row in detection_matrix],
            texttemplate='%{text}',
            textfont={"size": 10},
            colorbar=dict(title="Detection Rate (%)"),
            hovertemplate='Category: %{y}<br>Definition: %{x}<br>Detection Rate: %{z:.1f}%<extra></extra>'
        ))
        
        fig.update_layout(
            title='Detection Rates by PFAS Category',
            xaxis_title='PFAS Definition',
            yaxis_title='PFAS Category',
            height=600,
            width=1000,
        )
        
        output_file = os.path.join(self.figures_dir, 'detection_heatmap.html')
        fig.write_html(output_file)
        print(f"   Saved to {output_file}")
        
        return fig
    
    def plot_concordance_matrix(self):
        """Visualize inter-definition concordance"""
        print("\nCreating concordance matrix...")
        
        concordance_data = self.results['benchmarks']['concordance']
        jaccard_matrix = concordance_data['summary']['jaccard_similarity']
        
        # Convert to numpy array and ensure it's numeric
        try:
            if isinstance(jaccard_matrix, list):
                # Try to convert each element to float
                jaccard_matrix = [[float(x) for x in row] for row in jaccard_matrix]
                jaccard_matrix = np.array(jaccard_matrix, dtype=float)
            elif isinstance(jaccard_matrix, str):
                # Handle string representation of matrix - try to parse it
                try:
                    # Try to evaluate the string as a Python literal
                    import ast
                    matrix_data = ast.literal_eval(jaccard_matrix)
                    jaccard_matrix = np.array(matrix_data, dtype=float)
                    print("   Successfully parsed matrix from string")
                except:
                    # If parsing fails, use identity matrix
                    print("   Warning: Could not parse matrix string, using identity matrix")
                    jaccard_matrix = np.eye(5)
            else:
                jaccard_matrix = np.array(jaccard_matrix, dtype=float)
        except (ValueError, TypeError) as e:
            # If conversion fails, create a default matrix
            print(f"   Warning: Could not convert jaccard matrix to float ({e}), using identity matrix")
            jaccard_matrix = np.eye(5)  # Identity matrix as fallback
        
        # Create heatmap using matplotlib for better control
        fig, ax = plt.subplots(figsize=(10, 8))
        
        im = ax.imshow(jaccard_matrix, cmap='YlOrRd', aspect='auto', vmin=0, vmax=1)
        
        # Set ticks and labels
        ax.set_xticks(np.arange(5))
        ax.set_yticks(np.arange(5))
        ax.set_xticklabels([self.definition_names[i] for i in range(1, 6)], rotation=45, ha='right')
        ax.set_yticklabels([self.definition_names[i] for i in range(1, 6)])
        
        # Add colorbar
        cbar = plt.colorbar(im, ax=ax)
        cbar.set_label('Jaccard Similarity (0=no overlap, 1=perfect overlap)', rotation=270, labelpad=20)
        
        # Add text annotations
        for i in range(5):
            for j in range(5):
                text = ax.text(j, i, f"{jaccard_matrix[i, j]:.2f}",
                             ha="center", va="center", color="black", fontsize=12, fontweight='bold')
        
        ax.set_title('Inter-Definition Concordance (Jaccard Similarity)', fontsize=14, fontweight='bold')
        plt.tight_layout()
        
        output_file = os.path.join(self.figures_dir, 'concordance_matrix.png')
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        print(f"   Saved to {output_file}")
        plt.close()
        
        return fig
    
    def plot_edge_case_disagreements(self):
        """Visualize edge case disagreement patterns"""
        print("\n Creating edge case disagreement analysis...")
        
        edge_data = self.results['benchmarks']['edge_cases']
        disagreement_matrix = edge_data['summary'].get('disagreement_matrix', {})
        
        # Convert frozenset keys back to regular sets
        disagreement_counts = []
        pattern_labels = []
        
        for pattern_str, count in disagreement_matrix.items():
            # Parse the pattern
            if isinstance(pattern_str, str):
                # Pattern might be stored as string representation
                continue
            pattern = set(pattern_str) if pattern_str else set()
            
            if len(pattern) > 0 and len(pattern) < 5:
                pattern_names = [self.definition_names[i] for i in sorted(pattern)]
                pattern_label = ' + '.join(pattern_names)
                disagreement_counts.append(count)
                pattern_labels.append(pattern_label)
        
        if disagreement_counts:
            # Create bar chart
            fig = go.Figure(data=[
                go.Bar(
                    x=pattern_labels,
                    y=disagreement_counts,
                    marker_color='indianred',
                    text=disagreement_counts,
                    textposition='auto',
                )
            ])
            
            fig.update_layout(
                title='Edge Case Disagreement Patterns<br><sub>Which definitions agreed on borderline compounds</sub>',
                xaxis_title='Definition Combination',
                yaxis_title='Number of Compounds',
                height=600,
                width=1000,
            )
            
            output_file = os.path.join(self.figures_dir, 'edge_case_disagreements.html')
            fig.write_html(output_file)
            print(f"   Saved to {output_file}")
        else:
            print("   No disagreement data available")
            fig = None
        
        return fig
    
    def plot_performance_by_size(self):
        """Plot execution time vs molecule size"""
        print("\n  Creating performance analysis...")
        
        perf_data = self.results['benchmarks']['performance']
        measurements = perf_data['measurements']
        
        # Check if we have performance data
        if not measurements:
            print("   No performance data available")
            return
        
        # Extract data with type conversion to handle string data from JSON
        sizes = [int(m['num_atoms']) for m in measurements]
        times_ms = [float(m['execution_time']) * 1000 for m in measurements]
        
        # Additional safety check
        if not sizes or not times_ms:
            print("   No valid performance measurements found")
            return
        
        # Create scatter plot with trend
        fig = go.Figure()
        
        fig.add_trace(go.Scatter(
            x=sizes,
            y=times_ms,
            mode='markers',
            marker=dict(size=5, color='steelblue', opacity=0.5),
            name='Measurements',
            hovertemplate='Atoms: %{x}<br>Time: %{y:.2f} ms<extra></extra>'
        ))
        
        # Add trend line
        z = np.polyfit(sizes, times_ms, 2)
        p = np.poly1d(z)
        x_trend = np.linspace(min(sizes), max(sizes), 100)
        y_trend = p(x_trend)
        
        fig.add_trace(go.Scatter(
            x=x_trend,
            y=y_trend,
            mode='lines',
            line=dict(color='red', width=2),
            name='Trend (polynomial)',
            hovertemplate='Trend<extra></extra>'
        ))
        
        fig.update_layout(
            title='Definition Checking Performance vs Molecule Size',
            xaxis_title='Number of Atoms',
            yaxis_title='Execution Time (ms)',
            height=600,
            width=900,
            showlegend=True,
        )
        
        output_file = os.path.join(self.figures_dir, 'performance_by_size.html')
        fig.write_html(output_file)
        print(f"   Saved to {output_file}")
        
        return fig
    
    def generate_html_report(self):
        """Generate comprehensive HTML report"""
        print("\n Generating HTML report...")
        
        timestamp = self.results['metadata']['timestamp']
        
        # Generate all figures first
        summary_df = self.generate_summary_statistics()
        
        html_content = f"""
<!DOCTYPE html>
<html>
<head>
    <title>PFAS Definitions Benchmark Report</title>
    <style>
        body {{
            font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif;
            margin: 40px;
            background-color: #f5f5f5;
        }}
        .container {{
            max-width: 1200px;
            margin: 0 auto;
            background-color: white;
            padding: 40px;
            box-shadow: 0 0 20px rgba(0,0,0,0.1);
            border-radius: 10px;
        }}
        h1 {{
            color: #2c3e50;
            border-bottom: 3px solid #3498db;
            padding-bottom: 10px;
        }}
        h2 {{
            color: #34495e;
            margin-top: 40px;
            border-left: 5px solid #3498db;
            padding-left: 15px;
        }}
        h3 {{
            color: #7f8c8d;
        }}
        .metric-card {{
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            color: white;
            padding: 20px;
            border-radius: 10px;
            margin: 10px 0;
            box-shadow: 0 5px 15px rgba(0,0,0,0.2);
        }}
        .metric-value {{
            font-size: 2.5em;
            font-weight: bold;
        }}
        .metric-label {{
            font-size: 1.2em;
            opacity: 0.9;
        }}
        table {{
            width: 100%;
            border-collapse: collapse;
            margin: 20px 0;
            box-shadow: 0 2px 8px rgba(0,0,0,0.1);
        }}
        th, td {{
            padding: 12px 15px;
            text-align: left;
            border-bottom: 1px solid #ddd;
        }}
        th {{
            background-color: #3498db;
            color: white;
            font-weight: bold;
        }}
        tr:hover {{
            background-color: #f5f5f5;
        }}
        .definition-box {{
            background-color: #ecf0f1;
            border-left: 4px solid #3498db;
            padding: 15px;
            margin: 10px 0;
            border-radius: 5px;
        }}
        .definition-name {{
            font-weight: bold;
            color: #2c3e50;
            font-size: 1.1em;
        }}
        .definition-desc {{
            color: #7f8c8d;
            margin-top: 5px;
        }}
        .section {{
            margin: 40px 0;
        }}
        .highlight {{
            background-color: #fff9c4;
            padding: 2px 5px;
            border-radius: 3px;
        }}
        .good {{
            color: #27ae60;
            font-weight: bold;
        }}
        .warning {{
            color: #f39c12;
            font-weight: bold;
        }}
        .bad {{
            color: #e74c3c;
            font-weight: bold;
        }}
        iframe {{
            width: 100%;
            height: 700px;
            border: 1px solid #ddd;
            border-radius: 5px;
            margin: 20px 0;
        }}
        .footer {{
            margin-top: 60px;
            padding-top: 20px;
            border-top: 2px solid #ecf0f1;
            color: #95a5a6;
            text-align: center;
        }}
    </style>
</head>
<body>
    <div class="container">
        <h1>🧪 PFAS Definitions Benchmark Report</h1>
        
        <div class="section">
            <p><strong>Generated:</strong> {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}</p>
            <p><strong>Benchmark Run:</strong> {timestamp}</p>
        </div>
        
        <h2>📋 Definitions Tested</h2>
        <div class="section">
"""
        
        # Add definition descriptions
        definitions_desc = {
            1: "At least one CF₃ or CF₂ group (not bonded to H, Cl, Br, I)",
            2: "Similar to OECD but with more exclusion criteria",
            3: "Three inclusion criteria for perfluorinated compounds",
            4: "At least one CF₃ or CF₂ group (not bonded to Cl, Br, I)",
            5: "Custom structural patterns OR fluorine ratio ≥ 0.3"
        }
        
        for def_id in range(1, 6):
            html_content += f"""
            <div class="definition-box">
                <div class="definition-name">{def_id}. {self.definition_names[def_id]}</div>
                <div class="definition-desc">{definitions_desc[def_id]}</div>
            </div>
"""
        
        html_content += """
        </div>
        
        <h2> Summary Statistics</h2>
        <div class="section">
"""
        
        # Add summary table
        html_content += summary_df.to_html(index=False, classes='summary-table')
        
        # Add key metrics
        tp_summary = self.results['benchmarks']['true_positives']['summary']
        tn_summary = self.results['benchmarks']['true_negatives']['summary']
        edge_summary = self.results['benchmarks']['edge_cases']['summary']
        perf_summary = self.results['benchmarks']['performance']['summary']
        
        html_content += f"""
        </div>
        
        <h2>🎯 Key Findings</h2>
        <div class="section">
            <div class="metric-card">
                <div class="metric-label">True Positives Tested</div>
                <div class="metric-value">{tp_summary['total_molecules']}</div>
                <div class="metric-label">Known PFAS compounds</div>
            </div>
            
            <div class="metric-card">
                <div class="metric-label">True Negatives Tested</div>
                <div class="metric-value">{tn_summary['total_molecules']}</div>
                <div class="metric-label">Non-PFAS compounds</div>
            </div>
            
            <div class="metric-card">
                <div class="metric-label">Edge Cases Analyzed</div>
                <div class="metric-value">{edge_summary['total_molecules']}</div>
                <div class="metric-label">Borderline compounds</div>
            </div>
            
            <div class="metric-card">
                <div class="metric-label">Performance Testing</div>
                <div class="metric-value">{perf_summary['total_molecules']}</div>
                <div class="metric-label">Mean: {perf_summary['mean_time']*1000:.2f} ms</div>
            </div>
        </div>
        
        <h2> Sensitivity vs Specificity</h2>
        <div class="section">
            <p>This plot shows how well each definition balances correctly identifying PFAS (sensitivity) 
            with correctly rejecting non-PFAS (specificity). The ideal position is the top-right corner (100%, 100%).</p>
            <iframe src="figures_definitions/sensitivity_specificity.html"></iframe>
        </div>
        
        <h2> Detection Heatmap</h2>
        <div class="section">
            <p>Shows which definitions successfully identify different categories of PFAS compounds.</p>
            <iframe src="figures_definitions/detection_heatmap.html"></iframe>
        </div>
        
        <h2> Inter-Definition Concordance</h2>
        <div class="section">
            <p>Jaccard similarity between definitions (1.0 = perfect agreement, 0.0 = no overlap).</p>
            <img src="figures_definitions/concordance_matrix.png" style="max-width: 100%; height: auto;">
        </div>
        
        <h2> Performance Analysis</h2>
        <div class="section">
            <p>Execution time vs molecule size for definition checking.</p>
            <iframe src="figures_definitions/performance_by_size.html"></iframe>
        </div>
        
        <h2>🔍 Detailed Results</h2>
        <div class="section">
            <h3>True Positives (Known PFAS)</h3>
            <ul>
"""
        
        # Add detailed true positive results by category
        for category, count in Counter([r['category'] for subcat in tp_summary.get('subcategories', {}).values() 
                                       for r in subcat if r['valid']]).most_common():
            html_content += f"<li><strong>{category.replace('_', ' ').title()}</strong>: {count} compounds tested</li>\n"
        
        html_content += """
            </ul>
            
            <h3>True Negatives (Non-PFAS)</h3>
            <ul>
"""
        
        # Add false positive information
        for def_id in range(1, 6):
            fp_count = tn_summary['by_definition'][str(def_id)]
            if fp_count > 0:
                html_content += f"<li class='warning'>{self.definition_names[def_id]}: {fp_count} false positives</li>\n"
            else:
                html_content += f"<li class='good'>{self.definition_names[def_id]}: No false positives ✓</li>\n"
        
        html_content += f"""
            </ul>
            
            <h3>Edge Cases</h3>
            <p>Split decisions (definitions disagree): <span class="highlight">{edge_summary['split_decisions']}/{edge_summary['total_molecules']}</span></p>
            <p>Unanimous decisions (all agree): <span class="highlight">{edge_summary['unanimous_decisions']}/{edge_summary['total_molecules']}</span></p>
        </div>
        
        <div class="footer">
            <p>Generated by PFASGroups Benchmark Suite</p>
            <p>Comprehensive PFAS Definitions Benchmarking System</p>
        </div>
    </div>
</body>
</html>
"""
        
        # Save HTML report
        output_file = os.path.join(self.output_dir, f'definitions_benchmark_report_{timestamp}.html')
        with open(output_file, 'w', encoding='utf-8') as f:
            f.write(html_content)
        
        print(f"\n HTML report saved to: {output_file}")
        
        return output_file
    
    def run_complete_analysis(self):
        """Run all analysis and generate report"""
        print("\n" + "="*80)
        print(" ANALYZING PFAS DEFINITIONS BENCHMARK RESULTS")
        print("="*80)
        
        # Generate all visualizations
        self.plot_sensitivity_specificity()
        self.plot_detection_heatmap()
        self.plot_concordance_matrix()
        self.plot_edge_case_disagreements()
        self.plot_performance_by_size()
        
        # Generate HTML report
        report_file = self.generate_html_report()
        
        print("\n" + "="*80)
        print(" Analysis complete!")
        print(f"[FIGURES] Saved to: {self.figures_dir}")
        print(f" Report saved to: {report_file}")
        print("="*80)
        
        return report_file


def main():
    """Main function"""
    if len(sys.argv) < 2:
        print("Usage: python analyze_definitions_benchmark.py <benchmark_results.json>")
        print("\nSearching for latest benchmark file...")
        
        # Find latest benchmark file
        data_dir = os.path.join(parent_dir, 'data')
        benchmark_files = [f for f in os.listdir(data_dir) if f.startswith('pfas_definitions_benchmark_') and f.endswith('.json')]
        
        if not benchmark_files:
            print("❌ No benchmark files found. Run benchmark_pfas_definitions.py first.")
            return
        
        latest_file = sorted(benchmark_files)[-1]
        benchmark_file = os.path.join(data_dir, latest_file)
        print(f" Using: {latest_file}")
    else:
        benchmark_file = sys.argv[1]
    
    if not os.path.exists(benchmark_file):
        print(f"❌ File not found: {benchmark_file}")
        return
    
    analyzer = DefinitionBenchmarkAnalyzer(benchmark_file)
    report_file = analyzer.run_complete_analysis()
    
    return report_file


if __name__ == "__main__":
    main()
