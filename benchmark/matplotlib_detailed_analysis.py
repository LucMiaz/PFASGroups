import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.sankey import Sankey
import seaborn as sns
from collections import Counter, defaultdict
import json

def get_group_names():
    """Map group IDs to their chemical meanings"""
    group_names = {
        2: "Carboxylic acid",
        7: "Sulfonic acid", 
        10: "Polyfluoroalkylether sulfonic acid",
        14: "Alcohol-perfluoro",
        17: "Ether-perfluoro", 
        18: "Ethene-perfluoro",
        19: "Ethene-polyfluoro",
        20: "Alkane-perfluoro",
        21: "Polyfluoroalkylether",
        22: "Polyfluoroalkyl",
        23: "Alkane-perfluoro",
        25: "Perfluoroalkyl iodides",
        27: "Ketone-perfluoro",
        28: "Ketone-polyfluoro", 
        29: "Alcohol-polyfluoro",
        30: "Perfluoroalkyl ketones",
        31: "Ether",
        32: "Ester-perfluoro",
        33: "Ester",
        34: "Amide-perfluoro",
        35: "Acyl halide-perfluoro",
        36: "Sulfonic acid derivatives",
        37: "Sulfenic acid",
        38: "Sulfinic acid",
        39: "Phosphinic acid", 
        40: "Phosphonic acid",
        41: "Alkene-perfluoro",
        42: "Iodide-perfluoro",
        44: "Azine",
        45: "Azole", 
        46: "Ether-aromatic",
        47: "Aromatic",
        48: "Perfluoroalkyl",
        49: "Perfluoroalkene",
        50: "Alkyne-perfluoro",
        51: "Side-chain aromatics"
    }
    return group_names

def prepare_data():
    """Prepare data for analysis"""
    df = pd.read_csv('direct_benchmark_results.csv')
    group_names = get_group_names()
    
    print("🔄 Analyzing correspondence between PFAS-atlas and PFASGroups...")
    
    # Get atlas second class categories
    atlas_categories = df['atlas_second_class'].unique()
    print(f"Atlas categories found: {list(atlas_categories)}")
    
    # Parse PFASGroups data and count occurrences
    atlas_to_groups = defaultdict(lambda: defaultdict(int))
    
    for idx, row in df.iterrows():
        atlas_class = row['atlas_second_class']
        groups_str = row['detected_groups']
        
        if isinstance(groups_str, str) and groups_str.startswith('['):
            try:
                groups = eval(groups_str)
                for group_id in groups:
                    atlas_to_groups[atlas_class][group_id] += 1
            except:
                continue
    
    # Get top functional groups (limit for readability)
    all_group_counts = defaultdict(int)
    for class_groups in atlas_to_groups.values():
        for group_id, count in class_groups.items():
            all_group_counts[group_id] += count
    
    # Take top 15 groups for clarity
    top_groups = dict(sorted(all_group_counts.items(), key=lambda x: x[1], reverse=True)[:15])
    
    print(f"Top 15 functional groups selected for analysis:")
    for group_id, total_count in list(top_groups.items())[:10]:
        group_name = group_names.get(group_id, f"Group {group_id}")
        print(f"  G{group_id}: {group_name} ({total_count} total occurrences)")
    
    return df, atlas_to_groups, top_groups, group_names

def create_correspondence_heatmap():
    """Create a correspondence heatmap"""
    df, atlas_to_groups, top_groups, group_names = prepare_data()
    
    # Create correspondence matrix
    atlas_categories = sorted(atlas_to_groups.keys())
    group_ids = sorted(top_groups.keys())
    
    # Build matrix
    matrix = np.zeros((len(atlas_categories), len(group_ids)))
    
    for i, atlas_class in enumerate(atlas_categories):
        total_in_class = len(df[df['atlas_second_class'] == atlas_class])
        for j, group_id in enumerate(group_ids):
            count = atlas_to_groups[atlas_class].get(group_id, 0)
            # Store as percentage of molecules in this atlas class
            matrix[i, j] = count / total_in_class * 100 if total_in_class > 0 else 0
    
    # Create heatmap
    plt.figure(figsize=(16, 8))
    
    # Create custom labels
    atlas_labels = []
    for cat in atlas_categories:
        count = len(df[df['atlas_second_class'] == cat])
        atlas_labels.append(f"{cat}\n({count:,} molecules)")
    
    group_labels = []
    for group_id in group_ids:
        group_name = group_names.get(group_id, f"Group {group_id}")
        if len(group_name) > 20:
            group_name = group_name[:17] + "..."
        group_labels.append(f"G{group_id}\n{group_name}")
    
    # Create heatmap
    ax = sns.heatmap(matrix, 
                    xticklabels=group_labels,
                    yticklabels=atlas_labels,
                    annot=True, 
                    fmt='.1f',
                    cmap='YlOrRd',
                    cbar_kws={'label': 'Percentage of molecules in Atlas class'},
                    square=False)
    
    plt.title('PFAS-Atlas Second Class vs PFASGroups Functional Groups\nPercentage of molecules in each Atlas category containing each functional group', 
              fontsize=14, pad=20)
    plt.xlabel('PFASGroups Functional Groups', fontsize=12)
    plt.ylabel('PFAS-Atlas Second Class Categories', fontsize=12)
    
    # Rotate x-axis labels for better readability
    plt.xticks(rotation=45, ha='right')
    plt.yticks(rotation=0)
    
    plt.tight_layout()
    
    # Save in multiple formats
    plt.savefig('detailed_correspondence_heatmap.png', dpi=300, bbox_inches='tight')
    plt.savefig('detailed_correspondence_heatmap.pdf', bbox_inches='tight')
    plt.savefig('detailed_correspondence_heatmap.svg', bbox_inches='tight')
    
    print("✅ Correspondence heatmap saved: PNG, PDF, SVG")
    
    plt.show()
    return matrix, atlas_categories, group_ids

def create_flow_diagram():
    """Create a simplified flow diagram using matplotlib"""
    df, atlas_to_groups, top_groups, group_names = prepare_data()
    
    # Create a flow visualization
    fig, ax = plt.subplots(figsize=(16, 12))
    
    atlas_categories = sorted(atlas_to_groups.keys())
    group_ids = sorted(list(top_groups.keys())[:10])  # Top 10 for readability
    
    # Positions
    atlas_x = 0.1
    group_x = 0.8
    
    atlas_y_positions = np.linspace(0.1, 0.9, len(atlas_categories))
    group_y_positions = np.linspace(0.1, 0.9, len(group_ids))
    
    # Draw atlas categories
    atlas_colors = ['lightblue', 'lightcoral', 'lightgreen']
    for i, (cat, y_pos) in enumerate(zip(atlas_categories, atlas_y_positions)):
        count = len(df[df['atlas_second_class'] == cat])
        
        # Draw rectangle
        width = 0.15
        height = 0.08
        rect = plt.Rectangle((atlas_x - width/2, y_pos - height/2), 
                           width, height, 
                           facecolor=atlas_colors[i % len(atlas_colors)],
                           edgecolor='black',
                           linewidth=1)
        ax.add_patch(rect)
        
        # Add text
        text = f"{cat}\n{count:,} molecules"
        ax.text(atlas_x, y_pos, text, ha='center', va='center', 
                fontsize=9, fontweight='bold', wrap=True)
    
    # Draw functional groups
    group_colors = plt.cm.Oranges(np.linspace(0.3, 0.8, len(group_ids)))
    for i, (group_id, y_pos) in enumerate(zip(group_ids, group_y_positions)):
        group_name = group_names.get(group_id, f"Group {group_id}")
        count = top_groups[group_id]
        
        # Draw rectangle
        width = 0.15
        height = 0.06
        rect = plt.Rectangle((group_x - width/2, y_pos - height/2), 
                           width, height, 
                           facecolor=group_colors[i],
                           edgecolor='black',
                           linewidth=1)
        ax.add_patch(rect)
        
        # Add text
        if len(group_name) > 15:
            group_name = group_name[:12] + "..."
        text = f"G{group_id}: {group_name}\n{count:,} total"
        ax.text(group_x, y_pos, text, ha='center', va='center', 
                fontsize=8, wrap=True)
    
    # Draw connections (flows)
    for i, atlas_cat in enumerate(atlas_categories):
        atlas_y = atlas_y_positions[i]
        total_in_class = len(df[df['atlas_second_class'] == atlas_cat])
        
        for j, group_id in enumerate(group_ids):
            group_y = group_y_positions[j]
            count = atlas_to_groups[atlas_cat].get(group_id, 0)
            
            if count > 0:
                # Line thickness proportional to count
                max_count = max(atlas_to_groups[atlas_cat].values()) if atlas_to_groups[atlas_cat] else 1
                thickness = max(1, (count / max_count) * 8)
                
                # Draw connection
                ax.plot([atlas_x + 0.075, group_x - 0.075], 
                       [atlas_y, group_y], 
                       color='gray', 
                       alpha=0.6, 
                       linewidth=thickness)
                
                # Add count label for strong connections
                if count / total_in_class > 0.1:  # Show label if > 10%
                    mid_x = (atlas_x + group_x) / 2
                    mid_y = (atlas_y + group_y) / 2
                    percentage = count / total_in_class * 100
                    ax.text(mid_x, mid_y, f'{count}\n({percentage:.1f}%)', 
                           ha='center', va='center', fontsize=7,
                           bbox=dict(boxstyle='round,pad=0.2', 
                                   facecolor='white', alpha=0.8))
    
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_aspect('equal')
    ax.axis('off')
    
    # Add title and labels
    ax.text(0.5, 0.98, 'PFAS-Atlas Second Class → PFASGroups Functional Groups', 
           ha='center', va='top', fontsize=16, fontweight='bold')
    ax.text(0.5, 0.95, 'Flow diagram showing correspondence between classification systems', 
           ha='center', va='top', fontsize=12)
    
    ax.text(atlas_x, 0.02, 'PFAS-Atlas\nClassifications', ha='center', va='bottom', 
           fontsize=14, fontweight='bold')
    ax.text(group_x, 0.02, 'PFASGroups\nFunctional Groups', ha='center', va='bottom', 
           fontsize=14, fontweight='bold')
    
    plt.tight_layout()
    
    # Save in multiple formats
    plt.savefig('detailed_flow_diagram.png', dpi=300, bbox_inches='tight')
    plt.savefig('detailed_flow_diagram.pdf', bbox_inches='tight')
    plt.savefig('detailed_flow_diagram.svg', bbox_inches='tight')
    
    print("✅ Flow diagram saved: PNG, PDF, SVG")
    plt.show()

def create_statistics_table():
    """Create detailed statistics"""
    df, atlas_to_groups, top_groups, group_names = prepare_data()
    
    stats_data = []
    
    for atlas_class in sorted(atlas_to_groups.keys()):
        class_data = df[df['atlas_second_class'] == atlas_class]
        total_molecules = len(class_data)
        
        # Top 10 groups for this class
        class_group_counts = atlas_to_groups[atlas_class]
        top_groups_for_class = sorted(class_group_counts.items(), key=lambda x: x[1], reverse=True)[:10]
        
        for group_id, count in top_groups_for_class:
            group_name = group_names.get(group_id, f"Group {group_id}")
            percentage = count / total_molecules * 100
            
            stats_data.append({
                'Atlas_Class': atlas_class,
                'Total_Molecules_in_Class': total_molecules,
                'Group_ID': group_id,
                'Group_Name': group_name,
                'Molecules_with_Group': count,
                'Percentage_of_Class': round(percentage, 1),
                'Global_Rank': sorted(top_groups.keys(), key=lambda x: top_groups[x], reverse=True).index(group_id) + 1 if group_id in top_groups else 'N/A'
            })
    
    stats_df = pd.DataFrame(stats_data)
    stats_df.to_csv('detailed_correspondence_statistics.csv', index=False)
    
    print("\n📊 DETAILED CORRESPONDENCE STATISTICS")
    print("=" * 60)
    
    for atlas_class in sorted(atlas_to_groups.keys()):
        print(f"\n🎯 {atlas_class}:")
        class_group_counts = atlas_to_groups[atlas_class]
        top_5 = sorted(class_group_counts.items(), key=lambda x: x[1], reverse=True)[:5]
        
        total_in_class = len(df[df['atlas_second_class'] == atlas_class])
        
        for rank, (group_id, count) in enumerate(top_5, 1):
            group_name = group_names.get(group_id, f"Group {group_id}")
            percentage = count / total_in_class * 100
            print(f"  {rank}. G{group_id} ({group_name}): {count:,} molecules ({percentage:.1f}%)")
    
    print(f"\n💾 Detailed statistics saved: detailed_correspondence_statistics.csv")
    return stats_df

def create_html_summary():
    """Create an HTML summary page"""
    df, atlas_to_groups, top_groups, group_names = prepare_data()
    
    html_content = f"""
    <!DOCTYPE html>
    <html>
    <head>
        <title>PFAS-Atlas vs PFASGroups Correspondence Analysis</title>
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
            img {{ max-width: 100%; height: auto; margin: 20px 0; border: 1px solid #bdc3c7; }}
        </style>
    </head>
    <body>
        <h1>PFAS-Atlas vs PFASGroups Correspondence Analysis</h1>
        
        <div class="summary-box">
            <h2>Overview</h2>
            <p>This analysis examines how PFAS-Atlas second class classifications correspond to PFASGroups functional group detections across <strong>{len(df):,} molecules</strong>.</p>
            
            <div class="metric">
                <div class="metric-value">{len(df):,}</div>
                <div class="metric-label">Total Molecules</div>
            </div>
            <div class="metric">
                <div class="metric-value">{len(atlas_to_groups)}</div>
                <div class="metric-label">Atlas Categories</div>
            </div>
            <div class="metric">
                <div class="metric-value">{len(top_groups)}</div>
                <div class="metric-label">Top Functional Groups</div>
            </div>
        </div>
        
        <h2>Atlas Categories Summary</h2>
        <table>
            <tr><th>Atlas Category</th><th>Molecules</th><th>Percentage</th><th>Top Functional Groups</th></tr>
    """
    
    for atlas_class in sorted(atlas_to_groups.keys()):
        count = len(df[df['atlas_second_class'] == atlas_class])
        percentage = count / len(df) * 100
        
        # Get top 3 groups for this category
        class_group_counts = atlas_to_groups[atlas_class]
        top_3_groups = sorted(class_group_counts.items(), key=lambda x: x[1], reverse=True)[:3]
        top_groups_text = ", ".join([f"G{gid} ({group_names.get(gid, 'Unknown')})" for gid, _ in top_3_groups])
        
        html_content += f"""
            <tr>
                <td><strong>{atlas_class}</strong></td>
                <td>{count:,}</td>
                <td>{percentage:.1f}%</td>
                <td>{top_groups_text}</td>
            </tr>
        """
    
    html_content += """
        </table>
        
        <h2>Visualizations</h2>
        <h3>Correspondence Heatmap</h3>
        <img src="detailed_correspondence_heatmap.png" alt="Correspondence Heatmap">
        
        <h3>Flow Diagram</h3>
        <img src="detailed_flow_diagram.png" alt="Flow Diagram">
        
        <h2>Key Findings</h2>
        <div class="summary-box">
    """
    
    # Add key findings
    complex_structure_count = len(df[df['atlas_second_class'] == 'Complex structure'])
    aromatic_count = len(df[df['atlas_second_class'] == 'Aromatic PFASs'])
    
    html_content += f"""
            <ul>
                <li><strong>Complex Structure</strong> ({complex_structure_count:,} molecules): Shows diverse functional group patterns, dominated by perfluoroalkyl groups</li>
                <li><strong>Aromatic PFASs</strong> ({aromatic_count} molecules): Strong correspondence with polyfluoroalkyl and side-chain aromatic groups</li>
                <li><strong>Functional Group Diversity</strong>: Top 15 groups account for most molecular features across categories</li>
            </ul>
    """
    
    html_content += """
        </div>
        
        <h2>Files Generated</h2>
        <ul>
            <li><strong>detailed_correspondence_heatmap.png/pdf/svg</strong> - Correspondence heatmap</li>
            <li><strong>detailed_flow_diagram.png/pdf/svg</strong> - Flow diagram visualization</li>
            <li><strong>detailed_correspondence_statistics.csv</strong> - Detailed statistics table</li>
            <li><strong>detailed_summary.html</strong> - This summary report</li>
        </ul>
        
        <div class="summary-box">
            <p><em>Analysis generated on 2025-12-13 using PFASGroups benchmark data</em></p>
        </div>
        
    </body>
    </html>
    """
    
    with open('detailed_summary.html', 'w') as f:
        f.write(html_content)
    
    print("✅ HTML summary saved: detailed_summary.html")

def main():
    """Main function"""
    print("🚀 Creating detailed correspondence analysis...")
    print("=" * 60)
    
    try:
        # Create heatmap
        matrix, atlas_categories, group_ids = create_correspondence_heatmap()
        
        # Create flow diagram
        create_flow_diagram()
        
        # Create statistics
        stats_df = create_statistics_table()
        
        # Create HTML summary
        create_html_summary()
        
        print(f"\n🎉 CORRESPONDENCE ANALYSIS COMPLETE!")
        print("=" * 60)
        print(f"Generated files:")
        print(f"  📊 Heatmap: detailed_correspondence_heatmap.[png/pdf/svg]")
        print(f"  🔄 Flow diagram: detailed_flow_diagram.[png/pdf/svg]")
        print(f"  📈 Statistics: detailed_correspondence_statistics.csv")
        print(f"  🌐 Summary: detailed_summary.html")
        
        print(f"\n💡 All visualizations saved in PNG, PDF, and SVG formats")
        print(f"📱 Open detailed_summary.html in a browser for interactive overview")
        
    except Exception as e:
        print(f"❌ Error in analysis: {e}")
        import traceback
        traceback.print_exc()

if __name__ == "__main__":
    main()