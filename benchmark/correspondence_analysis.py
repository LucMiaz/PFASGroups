import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from collections import Counter, defaultdict
import plotly.graph_objects as go
import plotly.express as px
from plotly.subplots import make_subplots

def analyze_correspondence():
    """Analyze correspondence between PFAS-atlas and PFASGroups classifications"""
    
    # Load data
    df = pd.read_csv('direct_benchmark_results.csv')
    
    print("🗺️ CORRESPONDENCE ANALYSIS: PFAS-ATLAS VS PFASGROUPS")
    print("=" * 60)
    
    # Examine unique classifications
    print("\n📊 PFAS-Atlas Second Class Categories:")
    atlas_second = df['atlas_second_class'].value_counts()
    for cat, count in atlas_second.items():
        print(f"  {cat}: {count:,} molecules")
    
    # Parse PFASGroups functional groups
    print(f"\n🧪 PFASGroups Functional Group Analysis:")
    all_groups = []
    group_by_molecule = {}
    
    for idx, groups_str in enumerate(df['detected_groups'].dropna()):
        if isinstance(groups_str, str) and groups_str.startswith('['):
            try:
                groups = eval(groups_str)
                all_groups.extend(groups)
                group_by_molecule[idx] = groups
            except:
                pass
    
    group_counts = Counter(all_groups)
    print(f"Total unique groups: {len(group_counts)}")
    print(f"Most common groups:")
    for group_id, count in group_counts.most_common(10):
        print(f"  Group {group_id}: {count:,} molecules")
    
    return df, group_counts, group_by_molecule

def create_correspondence_matrix():
    """Create correspondence matrix between systems"""
    
    df, group_counts, group_by_molecule = analyze_correspondence()
    
    # Create mapping dictionary for group combinations
    atlas_to_groups = defaultdict(list)
    
    for idx, row in df.iterrows():
        atlas_class = row['atlas_second_class']
        if idx in group_by_molecule:
            groups = tuple(sorted(group_by_molecule[idx]))
            atlas_to_groups[atlas_class].append(groups)
    
    print(f"\n🔗 CORRESPONDENCE PATTERNS:")
    print("-" * 40)
    
    # Analyze most common group combinations for each atlas class
    correspondence_data = []
    
    for atlas_class, group_combinations in atlas_to_groups.items():
        combo_counts = Counter(group_combinations)
        print(f"\n{atlas_class} ({len(group_combinations)} molecules):")
        
        # Show top group combinations
        for combo, count in combo_counts.most_common(10):
            percentage = count / len(group_combinations) * 100
            groups_str = ', '.join([f"G{g}" for g in combo])
            print(f"  {groups_str}: {count} molecules ({percentage:.1f}%)")
            
            # Store for visualization
            correspondence_data.append({
                'atlas_class': atlas_class,
                'groups': groups_str,
                'group_tuple': combo,
                'count': count,
                'percentage': percentage,
                'total_in_class': len(group_combinations)
            })
    
    return df, correspondence_data

def create_sankey_diagram():
    """Create Sankey diagram showing correspondence"""
    
    df, correspondence_data = create_correspondence_matrix()
    
    # Prepare data for Sankey diagram
    # We'll show Atlas classes -> Most common group combinations
    
    atlas_classes = df['atlas_second_class'].unique()
    
    # Get top group combinations across all classes
    all_groups_flat = []
    for groups_str in df['detected_groups'].dropna():
        if isinstance(groups_str, str) and groups_str.startswith('['):
            try:
                groups = eval(groups_str)
                all_groups_flat.extend(groups)
            except:
                pass
    
    top_groups = Counter(all_groups_flat).most_common(15)
    
    # Create simplified mapping: Atlas class -> Individual groups (not combinations)
    atlas_to_individual_groups = defaultdict(lambda: defaultdict(int))
    
    for idx, row in df.iterrows():
        atlas_class = row['atlas_second_class']
        groups_str = row['detected_groups']
        
        if isinstance(groups_str, str) and groups_str.startswith('['):
            try:
                groups = eval(groups_str)
                for group in groups:
                    atlas_to_individual_groups[atlas_class][group] += 1
            except:
                pass
    
    # Build Sankey data
    source_nodes = []
    target_nodes = []
    values = []
    
    # Node labels
    node_labels = list(atlas_classes) + [f"Group {g}" for g, _ in top_groups]
    
    # Create mappings
    atlas_indices = {cls: i for i, cls in enumerate(atlas_classes)}
    group_indices = {g: len(atlas_classes) + i for i, (g, _) in enumerate(top_groups)}
    
    # Build flows
    for atlas_class, group_counts in atlas_to_individual_groups.items():
        atlas_idx = atlas_indices[atlas_class]
        
        for group_id, count in group_counts.items():
            if group_id in group_indices:  # Only include top groups
                group_idx = group_indices[group_id]
                source_nodes.append(atlas_idx)
                target_nodes.append(group_idx)
                values.append(count)
    
    # Create Sankey diagram
    fig = go.Figure(data=[go.Sankey(
        node=dict(
            pad=15,
            thickness=20,
            line=dict(color="black", width=0.5),
            label=node_labels,
            color=["lightblue"] * len(atlas_classes) + ["lightcoral"] * len(top_groups)
        ),
        link=dict(
            source=source_nodes,
            target=target_nodes,
            value=values,
            color="rgba(0, 100, 80, 0.3)"
        )
    )])
    
    fig.update_layout(
        title_text="PFAS-Atlas Classifications → PFASGroups Functional Groups",
        font_size=12,
        width=1200,
        height=800
    )
    
    fig.write_html("atlas_pfasgroups_sankey.html")
    print(f"\n✅ Sankey diagram saved as: atlas_pfasgroups_sankey.html")
    
    return fig

def create_heatmap_correspondence():
    """Create heatmap showing correspondence patterns"""
    
    df, correspondence_data = create_correspondence_matrix()
    
    # Create matrix: Atlas classes vs Top PFASGroups
    atlas_classes = df['atlas_second_class'].unique()
    
    # Get individual group counts per atlas class
    atlas_to_groups = defaultdict(lambda: defaultdict(int))
    
    for idx, row in df.iterrows():
        atlas_class = row['atlas_second_class']
        groups_str = row['detected_groups']
        
        if isinstance(groups_str, str) and groups_str.startswith('['):
            try:
                groups = eval(groups_str)
                for group in groups:
                    atlas_to_groups[atlas_class][group] += 1
            except:
                pass
    
    # Get top 20 groups across all classes
    all_group_counts = defaultdict(int)
    for class_groups in atlas_to_groups.values():
        for group, count in class_groups.items():
            all_group_counts[group] += count
    
    top_groups = [g for g, _ in sorted(all_group_counts.items(), 
                                     key=lambda x: x[1], reverse=True)[:20]]
    
    # Create matrix
    matrix_data = []
    for atlas_class in sorted(atlas_classes):
        row = []
        class_total = sum(atlas_to_groups[atlas_class].values())
        
        for group in top_groups:
            count = atlas_to_groups[atlas_class][group]
            percentage = (count / class_total * 100) if class_total > 0 else 0
            row.append(percentage)
        
        matrix_data.append(row)
    
    # Create heatmap
    plt.figure(figsize=(16, 10))
    
    heatmap_df = pd.DataFrame(matrix_data, 
                             index=sorted(atlas_classes),
                             columns=[f"Group {g}" for g in top_groups])
    
    sns.heatmap(heatmap_df, annot=True, fmt='.1f', cmap='YlOrRd', 
                cbar_kws={'label': 'Percentage (%)'})
    
    plt.title('Correspondence Heatmap: PFAS-Atlas Classifications vs PFASGroups\n(Percentage of molecules in each Atlas class containing each group)')
    plt.xlabel('PFASGroups Functional Groups')
    plt.ylabel('PFAS-Atlas Second Class')
    plt.xticks(rotation=45, ha='right')
    plt.yticks(rotation=0)
    plt.tight_layout()
    
    plt.savefig('atlas_pfasgroups_heatmap.png', dpi=300, bbox_inches='tight')
    plt.savefig('atlas_pfasgroups_heatmap.pdf', bbox_inches='tight')
    print(f"📊 Heatmap saved as: atlas_pfasgroups_heatmap.png/.pdf")
    
    return heatmap_df

def create_summary_table():
    """Create detailed correspondence summary"""
    
    df, correspondence_data = create_correspondence_matrix()
    
    # Create summary table
    summary_data = []
    
    atlas_classes = df['atlas_second_class'].unique()
    
    for atlas_class in sorted(atlas_classes):
        class_data = df[df['atlas_second_class'] == atlas_class]
        
        # Get all groups for this class
        all_groups_in_class = []
        for groups_str in class_data['detected_groups'].dropna():
            if isinstance(groups_str, str) and groups_str.startswith('['):
                try:
                    groups = eval(groups_str)
                    all_groups_in_class.extend(groups)
                except:
                    pass
        
        group_counts = Counter(all_groups_in_class)
        most_common = group_counts.most_common(5)
        
        avg_groups = class_data['n_detected_groups'].mean()
        total_molecules = len(class_data)
        specific_rate = (class_data['is_specific'] == True).mean() * 100
        
        summary_data.append({
            'Atlas_Class': atlas_class,
            'Molecules': total_molecules,
            'Avg_Groups': round(avg_groups, 1),
            'Specificity_Rate': round(specific_rate, 1),
            'Top_Groups': ', '.join([f"G{g}({c})" for g, c in most_common[:3]]),
            'Unique_Groups': len(group_counts)
        })
    
    summary_df = pd.DataFrame(summary_data)
    summary_df.to_csv('atlas_pfasgroups_correspondence_summary.csv', index=False)
    
    print(f"\n📋 CORRESPONDENCE SUMMARY TABLE:")
    print(summary_df.to_string(index=False))
    print(f"\n💾 Summary saved as: atlas_pfasgroups_correspondence_summary.csv")
    
    return summary_df

def main():
    """Main function to create all correspondence visualizations"""
    
    print("Creating correspondence analysis between PFAS-atlas and PFASGroups...")
    
    # Create all visualizations
    try:
        sankey_fig = create_sankey_diagram()
        print("✅ Sankey diagram created")
    except Exception as e:
        print(f"❌ Sankey diagram failed: {e}")
    
    try:
        heatmap_df = create_heatmap_correspondence()
        print("✅ Heatmap created")
    except Exception as e:
        print(f"❌ Heatmap failed: {e}")
    
    try:
        summary_df = create_summary_table()
        print("✅ Summary table created")
    except Exception as e:
        print(f"❌ Summary table failed: {e}")
    
    print(f"\n🎯 CORRESPONDENCE ANALYSIS COMPLETE!")
    print(f"Generated files:")
    print(f"  📊 atlas_pfasgroups_sankey.html - Interactive Sankey diagram")
    print(f"  🔥 atlas_pfasgroups_heatmap.png/.pdf - Correspondence heatmap")
    print(f"  📋 atlas_pfasgroups_correspondence_summary.csv - Summary statistics")

if __name__ == "__main__":
    main()