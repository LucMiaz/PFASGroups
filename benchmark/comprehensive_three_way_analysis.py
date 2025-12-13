import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from collections import defaultdict, Counter
import json
from matplotlib.patches import FancyBboxPatch, ConnectionPatch
import matplotlib.patches as mpatches

def load_group_definitions():
    """Load PFAS group definitions from JSON file"""
    with open('/home/luc/git/PFASGroups/PFASgroups/data/PFAS_groups_smarts.json', 'r') as f:
        groups_data = json.load(f)
    
    # Create lookup dictionary
    groups_lookup = {}
    for group in groups_data:
        groups_lookup[group['id']] = {
            'name': group['name'],
            'alias': group.get('alias', ''),
            'main_group': group.get('main_group', ''),
            'base_functional_groups': group.get('base_functional_groups', [])
        }
    
    return groups_lookup

def load_and_prepare_data():
    """Load and prepare the benchmark data"""
    print("🔄 Loading benchmark data...")
    
    # Load data
    df = pd.read_csv('pfas_atlas_benchmark_results.csv')
    groups_lookup = load_group_definitions()
    
    print(f"📊 Data loaded: {len(df):,} molecules")
    print(f"📚 Group definitions: {len(groups_lookup)} groups")
    
    # Parse detected groups
    def parse_groups(groups_str):
        if isinstance(groups_str, str) and groups_str.startswith('['):
            try:
                return eval(groups_str)
            except:
                return []
        return []
    
    df['parsed_groups'] = df['group_ids'].apply(parse_groups)
    
    return df, groups_lookup

def analyze_origin_to_atlas_accuracy():
    """Analyze how well PFAS-atlas predicts origin classification"""
    df, groups_lookup = load_and_prepare_data()
    
    print("\n🎯 ORIGIN vs PFAS-ATLAS ACCURACY ANALYSIS")
    print("=" * 60)
    
    # Create confusion matrix data
    origin_to_atlas = defaultdict(lambda: defaultdict(int))
    
    for _, row in df.iterrows():
        origin = row['origin']
        atlas_class = row['pfas_atlas_second_class']
        origin_to_atlas[origin][atlas_class] += 1
    
    # Calculate accuracy metrics
    accuracy_stats = []
    total_correct = 0
    total_molecules = len(df)
    
    for origin in sorted(origin_to_atlas.keys()):
        class_counts = origin_to_atlas[origin]
        total_in_origin = sum(class_counts.values())
        
        # Find most common atlas prediction for this origin
        most_common_atlas = max(class_counts.items(), key=lambda x: x[1])
        atlas_class, count = most_common_atlas
        
        accuracy = count / total_in_origin * 100
        total_correct += count
        
        accuracy_stats.append({
            'origin': origin,
            'total_molecules': total_in_origin,
            'best_atlas_match': atlas_class,
            'correct_predictions': count,
            'accuracy_percent': accuracy
        })
        
        print(f"🔸 {origin} → {atlas_class}: {count}/{total_in_origin} ({accuracy:.1f}%)")
    
    overall_accuracy = total_correct / total_molecules * 100
    print(f"\n📈 Overall PFAS-Atlas accuracy: {total_correct:,}/{total_molecules:,} ({overall_accuracy:.1f}%)")
    
    return origin_to_atlas, accuracy_stats, overall_accuracy

def create_three_way_sankey():
    """Create a three-way Sankey diagram showing all relationships"""
    df, groups_lookup = load_and_prepare_data()
    origin_to_atlas, accuracy_stats, overall_accuracy = analyze_origin_to_atlas_accuracy()
    
    print("\n🎨 Creating comprehensive three-way Sankey diagram...")
    
    # Prepare data structures
    origin_categories = sorted(df['origin'].unique())
    atlas_categories = sorted(df['pfas_atlas_second_class'].unique()) 
    
    # Get top functional groups for visualization clarity
    all_groups = []
    for _, row in df.iterrows():
        all_groups.extend(row['parsed_groups'])
    
    group_counts = Counter(all_groups)
    top_groups = dict(sorted(group_counts.items(), key=lambda x: x[1], reverse=True)[:20])
    
    print(f"📊 Selected {len(top_groups)} most frequent functional groups for visualization")
    
    # Create the visualization
    fig = plt.figure(figsize=(20, 14))
    gs = fig.add_gridspec(2, 2, height_ratios=[3, 1], width_ratios=[1, 1])
    
    # Main Sankey-style diagram
    ax_main = fig.add_subplot(gs[0, :])
    
    # Position calculations
    origin_x = 0.1
    atlas_x = 0.5
    groups_x = 0.9
    
    origin_y_positions = np.linspace(0.1, 0.9, len(origin_categories))
    atlas_y_positions = np.linspace(0.2, 0.8, len(atlas_categories))
    group_y_positions = np.linspace(0.1, 0.9, len(top_groups))
    
    # Color schemes
    origin_colors = plt.cm.Set3(np.linspace(0, 1, len(origin_categories)))
    atlas_colors = plt.cm.Set2(np.linspace(0, 1, len(atlas_categories)))
    group_colors = plt.cm.Oranges(np.linspace(0.3, 0.9, len(top_groups)))
    
    # Draw origin categories
    origin_positions = {}
    for i, (origin, y_pos) in enumerate(zip(origin_categories, origin_y_positions)):
        count = len(df[df['origin'] == origin])
        
        # Box size proportional to count
        box_height = 0.06 + (count / len(df)) * 0.1
        box_width = 0.12
        
        box = FancyBboxPatch((origin_x - box_width/2, y_pos - box_height/2),
                           box_width, box_height,
                           boxstyle="round,pad=0.01",
                           facecolor=origin_colors[i],
                           edgecolor='black',
                           linewidth=1)
        ax_main.add_patch(box)
        
        # Text label
        text = f"{origin}\n{count:,} molecules"
        ax_main.text(origin_x, y_pos, text, ha='center', va='center', 
                    fontsize=8, fontweight='bold', wrap=True)
        
        origin_positions[origin] = (origin_x, y_pos)
    
    # Draw atlas categories
    atlas_positions = {}
    for i, (atlas_cat, y_pos) in enumerate(zip(atlas_categories, atlas_y_positions)):
        count = len(df[df['pfas_atlas_second_class'] == atlas_cat])
        
        box_height = 0.06 + (count / len(df)) * 0.1
        box_width = 0.14
        
        box = FancyBboxPatch((atlas_x - box_width/2, y_pos - box_height/2),
                           box_width, box_height,
                           boxstyle="round,pad=0.01",
                           facecolor=atlas_colors[i],
                           edgecolor='black',
                           linewidth=1)
        ax_main.add_patch(box)
        
        # Truncate long names
        display_name = atlas_cat if len(atlas_cat) <= 25 else atlas_cat[:22] + "..."
        text = f"{display_name}\n{count:,} molecules"
        ax_main.text(atlas_x, y_pos, text, ha='center', va='center', 
                    fontsize=8, fontweight='bold', wrap=True)
        
        atlas_positions[atlas_cat] = (atlas_x, y_pos)
    
    # Draw functional groups
    group_positions = {}
    for i, (group_id, y_pos) in enumerate(zip(sorted(top_groups.keys()), group_y_positions)):
        count = top_groups[group_id]
        group_info = groups_lookup.get(group_id, {'name': f'Group {group_id}', 'alias': ''})
        
        box_height = 0.04 + (count / max(top_groups.values())) * 0.08
        box_width = 0.12
        
        box = FancyBboxPatch((groups_x - box_width/2, y_pos - box_height/2),
                           box_width, box_height,
                           boxstyle="round,pad=0.01",
                           facecolor=group_colors[i],
                           edgecolor='black',
                           linewidth=1)
        ax_main.add_patch(box)
        
        # Group name (truncated)
        group_name = group_info['name']
        if len(group_name) > 20:
            group_name = group_name[:17] + "..."
        
        text = f"G{group_id}: {group_name}\n{count:,} molecules"
        ax_main.text(groups_x, y_pos, text, ha='center', va='center', 
                    fontsize=7, wrap=True)
        
        group_positions[group_id] = (groups_x, y_pos)
    
    # Draw flows: Origin → Atlas (accuracy assessment)
    print("🔗 Drawing Origin → Atlas accuracy flows...")
    for origin in origin_categories:
        origin_x_pos, origin_y = origin_positions[origin]
        
        for atlas_cat in atlas_categories:
            atlas_x_pos, atlas_y = atlas_positions[atlas_cat]
            count = origin_to_atlas[origin][atlas_cat]
            
            if count > 0:
                total_in_origin = len(df[df['origin'] == origin])
                
                # Line thickness proportional to count
                thickness = max(1, (count / total_in_origin) * 15)
                
                # Color intensity based on whether this is the best match
                best_match = max(origin_to_atlas[origin].items(), key=lambda x: x[1])[0]
                alpha = 0.8 if atlas_cat == best_match else 0.3
                color = 'green' if atlas_cat == best_match else 'red'
                
                # Draw connection
                con = ConnectionPatch((origin_x_pos + 0.06, origin_y), 
                                    (atlas_x_pos - 0.07, atlas_y),
                                    "data", "data",
                                    arrowstyle="->", 
                                    shrinkA=0, shrinkB=0,
                                    color=color, alpha=alpha,
                                    linewidth=thickness)
                ax_main.add_patch(con)
    
    # Draw flows: Atlas → Groups
    print("🔗 Drawing Atlas → Groups correspondence flows...")
    atlas_to_groups = defaultdict(lambda: defaultdict(int))
    
    for _, row in df.iterrows():
        atlas_cat = row['pfas_atlas_second_class']
        for group_id in row['parsed_groups']:
            if group_id in top_groups:
                atlas_to_groups[atlas_cat][group_id] += 1
    
    for atlas_cat in atlas_categories:
        atlas_x_pos, atlas_y = atlas_positions[atlas_cat]
        total_in_atlas = len(df[df['pfas_atlas_second_class'] == atlas_cat])
        
        for group_id in sorted(top_groups.keys()):
            if group_id in group_positions:
                group_x_pos, group_y = group_positions[group_id]
                count = atlas_to_groups[atlas_cat][group_id]
                
                if count > 5:  # Only show significant connections
                    thickness = max(1, (count / total_in_atlas) * 10)
                    
                    con = ConnectionPatch((atlas_x_pos + 0.07, atlas_y), 
                                        (group_x_pos - 0.06, group_y),
                                        "data", "data",
                                        arrowstyle="->", 
                                        shrinkA=0, shrinkB=0,
                                        color='blue', alpha=0.5,
                                        linewidth=thickness)
                    ax_main.add_patch(con)
    
    # Formatting
    ax_main.set_xlim(0, 1)
    ax_main.set_ylim(0, 1)
    ax_main.set_aspect('equal')
    ax_main.axis('off')
    
    # Add title and labels
    ax_main.text(0.5, 0.98, 'Three-Way PFAS Classification Correspondence Analysis', 
               ha='center', va='top', fontsize=18, fontweight='bold')
    ax_main.text(0.5, 0.95, f'Origin vs PFAS-Atlas Accuracy • Atlas vs PFASGroups Correspondence • {len(df):,} molecules analyzed', 
               ha='center', va='top', fontsize=12)
    
    ax_main.text(origin_x, 0.02, 'Origin\nClassification', ha='center', va='bottom', 
               fontsize=14, fontweight='bold')
    ax_main.text(atlas_x, 0.02, 'PFAS-Atlas\nSecond Class', ha='center', va='bottom', 
               fontsize=14, fontweight='bold')
    ax_main.text(groups_x, 0.02, 'PFASGroups\nFunctional Groups', ha='center', va='bottom', 
               fontsize=14, fontweight='bold')
    
    # Create legend
    legend_elements = [
        plt.Line2D([0], [0], color='green', alpha=0.8, linewidth=3, label='Best Atlas Match (High Accuracy)'),
        plt.Line2D([0], [0], color='red', alpha=0.3, linewidth=3, label='Poor Atlas Match (Low Accuracy)'),
        plt.Line2D([0], [0], color='blue', alpha=0.5, linewidth=3, label='Atlas → Groups Correspondence')
    ]
    ax_main.legend(handles=legend_elements, loc='upper left', fontsize=10)
    
    # Add accuracy statistics subplot
    ax_stats = fig.add_subplot(gs[1, :])
    ax_stats.axis('off')
    
    # Create accuracy summary table
    stats_text = f"📈 PFAS-ATLAS ACCURACY SUMMARY (Overall: {overall_accuracy:.1f}%)\n"
    stats_text += "─" * 100 + "\n"
    
    for stat in accuracy_stats[:10]:  # Show top 10
        stats_text += f"{stat['origin']:<30} → {stat['best_atlas_match']:<35} : {stat['correct_predictions']:>4}/{stat['total_molecules']:<4} ({stat['accuracy_percent']:>5.1f}%)\n"
    
    ax_stats.text(0.05, 0.95, stats_text, transform=ax_stats.transAxes, 
                 fontsize=9, fontfamily='monospace', va='top')
    
    plt.tight_layout()
    
    # Save in multiple formats
    plt.savefig('three_way_sankey_analysis.png', dpi=300, bbox_inches='tight')
    plt.savefig('three_way_sankey_analysis.pdf', bbox_inches='tight')
    plt.savefig('three_way_sankey_analysis.svg', bbox_inches='tight')
    
    print("✅ Three-way Sankey diagram saved: PNG, PDF, SVG")
    
    return fig, origin_to_atlas, atlas_to_groups, accuracy_stats

def create_detailed_accuracy_heatmap():
    """Create detailed accuracy heatmap"""
    df, groups_lookup = load_and_prepare_data()
    origin_to_atlas, accuracy_stats, overall_accuracy = analyze_origin_to_atlas_accuracy()
    
    print("🔥 Creating detailed accuracy heatmap...")
    
    # Create confusion matrix
    origins = sorted(origin_to_atlas.keys())
    atlas_cats = sorted(set().union(*[set(d.keys()) for d in origin_to_atlas.values()]))
    
    matrix = np.zeros((len(origins), len(atlas_cats)))
    
    for i, origin in enumerate(origins):
        total_in_origin = sum(origin_to_atlas[origin].values())
        for j, atlas_cat in enumerate(atlas_cats):
            count = origin_to_atlas[origin].get(atlas_cat, 0)
            matrix[i, j] = count / total_in_origin * 100 if total_in_origin > 0 else 0
    
    # Create heatmap
    plt.figure(figsize=(16, 10))
    
    # Prepare labels
    origin_labels = [f"{origin}\n({sum(origin_to_atlas[origin].values())} molecules)" for origin in origins]
    atlas_labels = [f"{cat}\n({len(df[df['pfas_atlas_second_class'] == cat])} molecules)" for cat in atlas_cats]
    
    # Create heatmap
    ax = sns.heatmap(matrix, 
                    xticklabels=atlas_labels,
                    yticklabels=origin_labels,
                    annot=True, 
                    fmt='.1f',
                    cmap='RdYlGn',
                    center=50,
                    cbar_kws={'label': 'Percentage of origin molecules predicted as Atlas class'},
                    square=False)
    
    plt.title(f'PFAS-Atlas Classification Accuracy Heatmap\nOverall Accuracy: {overall_accuracy:.1f}% ({len(df):,} molecules)', 
              fontsize=14, pad=20)
    plt.xlabel('PFAS-Atlas Second Class Predictions', fontsize=12)
    plt.ylabel('Origin (True) Classifications', fontsize=12)
    
    plt.xticks(rotation=45, ha='right')
    plt.yticks(rotation=0)
    
    plt.tight_layout()
    
    # Save
    plt.savefig('accuracy_heatmap_detailed.png', dpi=300, bbox_inches='tight')
    plt.savefig('accuracy_heatmap_detailed.pdf', bbox_inches='tight')
    plt.savefig('accuracy_heatmap_detailed.svg', bbox_inches='tight')
    
    print("✅ Accuracy heatmap saved: PNG, PDF, SVG")
    
    return matrix, origins, atlas_cats

def generate_comprehensive_statistics():
    """Generate comprehensive statistics"""
    df, groups_lookup = load_and_prepare_data()
    origin_to_atlas, accuracy_stats, overall_accuracy = analyze_origin_to_atlas_accuracy()
    
    print("📊 Generating comprehensive statistics...")
    
    # Create detailed statistics
    detailed_stats = []
    
    for _, row in df.iterrows():
        origin = row['origin']
        atlas_class = row['pfas_atlas_second_class'] 
        groups = row['parsed_groups']
        
        # Check if atlas prediction matches origin expectation
        origin_atlas_pairs = origin_to_atlas[origin]
        best_atlas_match = max(origin_atlas_pairs.items(), key=lambda x: x[1])[0]
        is_accurate_prediction = (atlas_class == best_atlas_match)
        
        # Group information
        group_names = []
        for group_id in groups:
            group_info = groups_lookup.get(group_id, {'name': f'Group {group_id}'})
            group_names.append(f"G{group_id}:{group_info['name']}")
        
        detailed_stats.append({
            'origin': origin,
            'pfas_atlas_prediction': atlas_class,
            'best_atlas_match_for_origin': best_atlas_match,
            'is_accurate_prediction': is_accurate_prediction,
            'detected_groups': groups,
            'group_names': '; '.join(group_names),
            'n_groups_detected': len(groups)
        })
    
    stats_df = pd.DataFrame(detailed_stats)
    stats_df.to_csv('comprehensive_classification_analysis.csv', index=False)
    
    # Summary statistics
    accuracy_by_origin = stats_df.groupby('origin').agg({
        'is_accurate_prediction': ['count', 'sum', 'mean'],
        'n_groups_detected': ['mean', 'std']
    }).round(3)
    
    accuracy_by_origin.columns = ['total_molecules', 'correct_predictions', 'accuracy_rate', 
                                'avg_groups_detected', 'std_groups_detected']
    accuracy_by_origin.to_csv('accuracy_by_origin_summary.csv')
    
    print(f"💾 Comprehensive statistics saved:")
    print(f"  • comprehensive_classification_analysis.csv")
    print(f"  • accuracy_by_origin_summary.csv")
    
    return stats_df, accuracy_by_origin

def create_html_report():
    """Create comprehensive HTML report"""
    df, groups_lookup = load_and_prepare_data()
    origin_to_atlas, accuracy_stats, overall_accuracy = analyze_origin_to_atlas_accuracy()
    
    html_content = f"""
    <!DOCTYPE html>
    <html>
    <head>
        <title>Comprehensive PFAS Classification Analysis</title>
        <style>
            body {{ font-family: Arial, sans-serif; margin: 40px; line-height: 1.6; }}
            h1 {{ color: #2c3e50; text-align: center; }}
            h2 {{ color: #3498db; border-bottom: 2px solid #ecf0f1; padding-bottom: 10px; }}
            h3 {{ color: #e74c3c; }}
            table {{ border-collapse: collapse; width: 100%; margin: 20px 0; }}
            th, td {{ border: 1px solid #bdc3c7; padding: 12px; text-align: left; }}
            th {{ background-color: #3498db; color: white; }}
            tr:nth-child(even) {{ background-color: #f8f9fa; }}
            .summary-box {{ background-color: #ecf0f1; padding: 20px; border-radius: 8px; margin: 20px 0; }}
            .metric {{ display: inline-block; margin: 10px 20px; text-align: center; }}
            .metric-value {{ font-size: 2em; font-weight: bold; color: #e74c3c; }}
            .metric-label {{ font-size: 0.9em; color: #7f8c8d; }}
            .accuracy-good {{ color: #27ae60; font-weight: bold; }}
            .accuracy-poor {{ color: #e74c3c; font-weight: bold; }}
            img {{ max-width: 100%; height: auto; margin: 20px 0; border: 1px solid #bdc3c7; }}
        </style>
    </head>
    <body>
        <h1>Comprehensive PFAS Classification Analysis</h1>
        
        <div class="summary-box">
            <h2>Executive Summary</h2>
            <p>This analysis evaluates the accuracy of PFAS-Atlas classifications against origin (ground truth) classifications 
            and examines correspondence with PFASGroups functional group detections across <strong>{len(df):,} molecules</strong>.</p>
            
            <div class="metric">
                <div class="metric-value">{overall_accuracy:.1f}%</div>
                <div class="metric-label">Overall Atlas Accuracy</div>
            </div>
            <div class="metric">
                <div class="metric-value">{len(origin_to_atlas)}</div>
                <div class="metric-label">Origin Categories</div>
            </div>
            <div class="metric">
                <div class="metric-value">{len(set(df['pfas_atlas_second_class']))}</div>
                <div class="metric-label">Atlas Categories</div>
            </div>
            <div class="metric">
                <div class="metric-value">{len(groups_lookup)}</div>
                <div class="metric-label">Total Func Groups</div>
            </div>
        </div>
        
        <h2>PFAS-Atlas Accuracy by Origin Category</h2>
        <table>
            <tr><th>Origin Category</th><th>Best Atlas Match</th><th>Correct/Total</th><th>Accuracy</th><th>Assessment</th></tr>
    """
    
    for stat in sorted(accuracy_stats, key=lambda x: x['accuracy_percent'], reverse=True):
        accuracy_class = 'accuracy-good' if stat['accuracy_percent'] >= 70 else 'accuracy-poor'
        assessment = 'Excellent' if stat['accuracy_percent'] >= 90 else 'Good' if stat['accuracy_percent'] >= 70 else 'Poor'
        
        html_content += f"""
            <tr>
                <td><strong>{stat['origin']}</strong></td>
                <td>{stat['best_atlas_match']}</td>
                <td>{stat['correct_predictions']:,}/{stat['total_molecules']:,}</td>
                <td class="{accuracy_class}">{stat['accuracy_percent']:.1f}%</td>
                <td class="{accuracy_class}">{assessment}</td>
            </tr>
        """
    
    html_content += f"""
        </table>
        
        <h2>Key Findings</h2>
        <div class="summary-box">
            <ul>
                <li><strong>Overall Performance:</strong> PFAS-Atlas achieves {overall_accuracy:.1f}% accuracy in predicting origin classifications</li>
                <li><strong>Best Performing:</strong> {max(accuracy_stats, key=lambda x: x['accuracy_percent'])['origin']} ({max(accuracy_stats, key=lambda x: x['accuracy_percent'])['accuracy_percent']:.1f}% accuracy)</li>
                <li><strong>Challenging Categories:</strong> {min(accuracy_stats, key=lambda x: x['accuracy_percent'])['origin']} ({min(accuracy_stats, key=lambda x: x['accuracy_percent'])['accuracy_percent']:.1f}% accuracy)</li>
                <li><strong>Total Molecules Analyzed:</strong> {len(df):,}</li>
                <li><strong>Classification Systems:</strong> 3-way comparison (Origin, PFAS-Atlas, PFASGroups)</li>
            </ul>
        </div>
        
        <h2>Visualizations</h2>
        <h3>Three-Way Sankey Diagram</h3>
        <img src="three_way_sankey_analysis.png" alt="Three-way Sankey Analysis">
        
        <h3>Accuracy Heatmap</h3>
        <img src="accuracy_heatmap_detailed.png" alt="Accuracy Heatmap">
        
        <h2>Files Generated</h2>
        <ul>
            <li><strong>three_way_sankey_analysis.[png/pdf/svg]</strong> - Comprehensive three-way flow diagram</li>
            <li><strong>accuracy_heatmap_detailed.[png/pdf/svg]</strong> - Accuracy confusion matrix</li>
            <li><strong>comprehensive_classification_analysis.csv</strong> - Detailed per-molecule statistics</li>
            <li><strong>accuracy_by_origin_summary.csv</strong> - Summary statistics by origin category</li>
            <li><strong>comprehensive_analysis_report.html</strong> - This report</li>
        </ul>
        
        <div class="summary-box">
            <p><em>Analysis generated on 2025-12-13 using PFASGroups benchmark suite</em></p>
        </div>
        
    </body>
    </html>
    """
    
    with open('comprehensive_analysis_report.html', 'w') as f:
        f.write(html_content)
    
    print("✅ HTML report saved: comprehensive_analysis_report.html")

def main():
    """Main analysis function"""
    print("🚀 Starting Comprehensive PFAS Classification Analysis...")
    print("=" * 80)
    
    try:
        # Create three-way Sankey diagram
        fig, origin_to_atlas, atlas_to_groups, accuracy_stats = create_three_way_sankey()
        
        # Create detailed accuracy heatmap
        matrix, origins, atlas_cats = create_detailed_accuracy_heatmap()
        
        # Generate comprehensive statistics
        stats_df, accuracy_by_origin = generate_comprehensive_statistics()
        
        # Create HTML report
        create_html_report()
        
        print(f"\n🎉 COMPREHENSIVE ANALYSIS COMPLETE!")
        print("=" * 80)
        print(f"Generated deliverables:")
        print(f"  🎨 Three-way Sankey: three_way_sankey_analysis.[png/pdf/svg]")
        print(f"  🔥 Accuracy heatmap: accuracy_heatmap_detailed.[png/pdf/svg]")
        print(f"  📊 Statistics: comprehensive_classification_analysis.csv")
        print(f"  📈 Summary: accuracy_by_origin_summary.csv")
        print(f"  🌐 Report: comprehensive_analysis_report.html")
        
        # Print key insights
        overall_accuracy = stats_df['is_accurate_prediction'].mean() * 100
        best_category = accuracy_stats[max(range(len(accuracy_stats)), key=lambda x: accuracy_stats[x]['accuracy_percent'])]
        worst_category = accuracy_stats[min(range(len(accuracy_stats)), key=lambda x: accuracy_stats[x]['accuracy_percent'])]
        
        print(f"\n💡 Key Insights:")
        print(f"  • Overall PFAS-Atlas accuracy: {overall_accuracy:.1f}%")
        print(f"  • Best performing category: {best_category['origin']} ({best_category['accuracy_percent']:.1f}%)")
        print(f"  • Most challenging category: {worst_category['origin']} ({worst_category['accuracy_percent']:.1f}%)")
        print(f"  • Total molecules analyzed: {len(stats_df):,}")
        
    except Exception as e:
        print(f"❌ Error in analysis: {e}")
        import traceback
        traceback.print_exc()

if __name__ == "__main__":
    main()