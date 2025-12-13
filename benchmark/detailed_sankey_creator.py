import pandas as pd
import numpy as np
import plotly.graph_objects as go
import plotly.express as px
from plotly.subplots import make_subplots
import matplotlib.pyplot as plt
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

def prepare_sankey_data():
    """Prepare data for detailed Sankey diagram"""
    
    df = pd.read_csv('direct_benchmark_results.csv')
    group_names = get_group_names()
    
    print("🔄 Preparing detailed Sankey diagram data...")
    
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
    
    # Take top 20 groups for clarity
    top_groups = dict(sorted(all_group_counts.items(), key=lambda x: x[1], reverse=True)[:20])
    
    print(f"Top 20 functional groups selected for Sankey:")
    for group_id, total_count in list(top_groups.items())[:10]:
        group_name = group_names.get(group_id, f"Group {group_id}")
        print(f"  G{group_id}: {group_name} ({total_count} total occurrences)")
    
    return df, atlas_to_groups, top_groups, group_names

def create_enhanced_sankey():
    """Create enhanced Sankey diagram with better styling"""
    
    df, atlas_to_groups, top_groups, group_names = prepare_sankey_data()
    
    # Prepare node labels and colors
    atlas_categories = sorted(df['atlas_second_class'].unique())
    
    # Create comprehensive node labels
    source_labels = []
    target_labels = []
    
    # Atlas categories (source nodes)
    for cat in atlas_categories:
        count = len(df[df['atlas_second_class'] == cat])
        source_labels.append(f"{cat}\n({count:,} molecules)")
    
    # Functional groups (target nodes)  
    for group_id in sorted(top_groups.keys()):
        group_name = group_names.get(group_id, f"Group {group_id}")
        total_count = top_groups[group_id]
        # Truncate long names
        if len(group_name) > 15:
            group_name = group_name[:12] + "..."
        target_labels.append(f"G{group_id}: {group_name}\n({total_count:,} molecules)")
    
    # All node labels
    all_labels = source_labels + target_labels
    
    # Create index mappings
    atlas_indices = {cat: i for i, cat in enumerate(atlas_categories)}
    group_indices = {group_id: len(atlas_categories) + i for i, group_id in enumerate(sorted(top_groups.keys()))}
    
    # Build flows (links)
    source_nodes = []
    target_nodes = []
    values = []
    link_labels = []
    
    for atlas_class, group_counts in atlas_to_groups.items():
        atlas_idx = atlas_indices[atlas_class]
        
        for group_id, count in group_counts.items():
            if group_id in top_groups:  # Only include top groups
                group_idx = group_indices[group_id]
                source_nodes.append(atlas_idx)
                target_nodes.append(group_idx)
                values.append(count)
                
                # Calculate percentage for this atlas class
                total_in_class = len(df[df['atlas_second_class'] == atlas_class])
                percentage = count / total_in_class * 100
                link_labels.append(f"{count} molecules ({percentage:.1f}%)")
    
    # Color schemes
    atlas_colors = {
        'Complex structure': 'rgba(31, 119, 180, 0.8)',
        'Aromatic PFASs': 'rgba(255, 127, 14, 0.8)', 
        'Not PFAS by current definition': 'rgba(44, 160, 44, 0.8)'
    }
    
    # Node colors
    node_colors = []
    for cat in atlas_categories:
        node_colors.append(atlas_colors.get(cat, 'rgba(128, 128, 128, 0.8)'))
    
    # Functional group colors (gradient from light to dark coral)
    group_color_base = 'rgba(255, 99, 71, {})'  # Tomato color
    for i, group_id in enumerate(sorted(top_groups.keys())):
        opacity = 0.6 + (0.3 * (top_groups[group_id] / max(top_groups.values())))
        node_colors.append(group_color_base.format(opacity))
    
    # Create Sankey diagram
    fig = go.Figure(data=[go.Sankey(
        arrangement="snap",
        node=dict(
            pad=20,
            thickness=25,
            line=dict(color="black", width=1),
            label=all_labels,
            color=node_colors,
            x=[0.1] * len(atlas_categories) + [0.9] * len(top_groups),
            y=list(np.linspace(0.1, 0.9, len(atlas_categories))) + list(np.linspace(0.05, 0.95, len(top_groups)))
        ),
        link=dict(
            source=source_nodes,
            target=target_nodes,
            value=values,
            color=['rgba(0, 100, 80, 0.4)'] * len(values),
            hovertemplate='<b>%{source.label}</b> → <b>%{target.label}</b><br>' +
                         'Molecules: %{value}<br>' +
                         '<extra></extra>'
        )
    )])
    
    # Update layout with enhanced styling
    fig.update_layout(
        title={
            'text': 'PFAS-Atlas Second Class Classifications → PFASGroups Functional Groups<br>' +
                   '<sub>Molecular correspondence analysis across 1,832 molecules</sub>',
            'x': 0.5,
            'font': {'size': 18, 'family': 'Arial, sans-serif'}
        },
        font=dict(size=12, family='Arial, sans-serif'),
        width=1400,
        height=900,
        margin=dict(l=50, r=50, t=100, b=50),
        plot_bgcolor='white',
        paper_bgcolor='white'
    )
    
    # Add annotations
    annotations = [
        dict(
            x=0.05, y=0.95,
            text='<b>PFAS-Atlas<br>Classifications</b>',
            showarrow=False,
            font=dict(size=14, color='black'),
            align='center'
        ),
        dict(
            x=0.95, y=0.95,
            text='<b>PFASGroups<br>Functional Groups</b>',
            showarrow=False,
            font=dict(size=14, color='black'),
            align='center'
        ),
        dict(
            x=0.5, y=0.02,
            text='Flow width represents number of molecules containing each functional group in each Atlas category',
            showarrow=False,
            font=dict(size=10, color='gray'),
            align='center'
        )
    ]
    
    fig.update_layout(annotations=annotations)
    
    return fig, df, atlas_to_groups, top_groups, group_names

def save_multiple_formats():
    """Save Sankey diagram in multiple formats"""
    
    print("🎨 Creating enhanced Sankey diagram...")
    fig, df, atlas_to_groups, top_groups, group_names = create_enhanced_sankey()
    
    # Save HTML (interactive) - always works
    html_file = "detailed_sankey_atlas_pfasgroups.html"
    fig.write_html(html_file, include_plotlyjs=True, config={'displayModeBar': True})
    print(f"✅ Interactive HTML saved: {html_file}")
    
    # Try to save image formats (requires kaleido)
    png_file = pdf_file = svg_file = None
    
    try:
        # Save PNG (high resolution)
        png_file = "detailed_sankey_atlas_pfasgroups.png"
        fig.write_image(png_file, format="png", width=1400, height=900, scale=2)
        print(f"✅ PNG image saved: {png_file}")
        
        # Save PDF
        pdf_file = "detailed_sankey_atlas_pfasgroups.pdf" 
        fig.write_image(pdf_file, format="pdf", width=1400, height=900)
        print(f"✅ PDF saved: {pdf_file}")
        
        # Save SVG (vector graphics)
        svg_file = "detailed_sankey_atlas_pfasgroups.svg"
        fig.write_image(svg_file, format="svg", width=1400, height=900)
        print(f"✅ SVG saved: {svg_file}")
        
    except Exception as e:
        print(f"⚠️  Image export failed (kaleido not available): {e}")
        print(f"📱 Only interactive HTML format generated")
        print(f"💡 To generate PNG/PDF/SVG formats, install kaleido: pip install kaleido")
    
    return fig, html_file, png_file, pdf_file, svg_file

def create_correspondence_statistics():
    """Generate detailed statistics for the correspondence"""
    
    df, atlas_to_groups, top_groups, group_names = prepare_sankey_data()
    
    print(f"\n📊 DETAILED CORRESPONDENCE STATISTICS")
    print("=" * 60)
    
    # Create detailed statistics table
    stats_data = []
    
    for atlas_class in sorted(atlas_to_groups.keys()):
        class_data = df[df['atlas_second_class'] == atlas_class]
        total_molecules = len(class_data)
        
        # Top 5 groups for this class
        class_group_counts = atlas_to_groups[atlas_class]
        top_5_groups = sorted(class_group_counts.items(), key=lambda x: x[1], reverse=True)[:5]
        
        for group_id, count in top_5_groups:
            if group_id in top_groups:  # Only include groups in our top list
                group_name = group_names.get(group_id, f"Group {group_id}")
                percentage = count / total_molecules * 100
                
                stats_data.append({
                    'Atlas_Class': atlas_class,
                    'Total_Molecules_in_Class': total_molecules,
                    'Group_ID': group_id,
                    'Group_Name': group_name,
                    'Molecules_with_Group': count,
                    'Percentage_of_Class': round(percentage, 1),
                    'Global_Group_Total': top_groups[group_id]
                })
    
    stats_df = pd.DataFrame(stats_data)
    stats_df.to_csv('detailed_sankey_statistics.csv', index=False)
    
    # Print summary
    print(f"Atlas Categories and their strongest group associations:")
    for atlas_class in sorted(atlas_to_groups.keys()):
        print(f"\n🎯 {atlas_class}:")
        class_group_counts = atlas_to_groups[atlas_class]
        top_3 = sorted(class_group_counts.items(), key=lambda x: x[1], reverse=True)[:3]
        
        total_in_class = len(df[df['atlas_second_class'] == atlas_class])
        
        for rank, (group_id, count) in enumerate(top_3, 1):
            group_name = group_names.get(group_id, f"Group {group_id}")
            percentage = count / total_in_class * 100
            print(f"  {rank}. G{group_id} ({group_name}): {count:,} molecules ({percentage:.1f}%)")
    
    print(f"\n💾 Detailed statistics saved: detailed_sankey_statistics.csv")
    return stats_df

def create_sankey_metadata():
    """Create metadata file describing the visualization"""
    
    df = pd.read_csv('direct_benchmark_results.csv')
    
    metadata = {
        "diagram_info": {
            "title": "PFAS-Atlas vs PFASGroups Correspondence",
            "description": "Sankey diagram showing flow from PFAS-Atlas second class classifications to PFASGroups functional groups",
            "created_date": "2025-12-13",
            "total_molecules": len(df),
            "data_source": "direct_benchmark_results.csv"
        },
        "atlas_categories": {
            "Complex structure": {
                "molecules": len(df[df['atlas_second_class'] == 'Complex structure']),
                "description": "Diverse PFAS structures with complex molecular features"
            },
            "Aromatic PFASs": {
                "molecules": len(df[df['atlas_second_class'] == 'Aromatic PFASs']),
                "description": "PFAS compounds containing aromatic ring systems"
            },
            "Not PFAS by current definition": {
                "molecules": len(df[df['atlas_second_class'] == 'Not PFAS by current definition']),
                "description": "Compounds with PFAS-like features but classified as non-PFAS"
            }
        },
        "visualization_notes": {
            "flow_interpretation": "Flow width represents number of molecules in each Atlas category that contain each PFASGroups functional group",
            "color_scheme": "Atlas categories in blue/orange/green, functional groups in coral gradient by frequency",
            "top_groups_shown": 20,
            "total_unique_groups": 51
        },
        "files_generated": [
            "detailed_sankey_atlas_pfasgroups.html",
            "detailed_sankey_atlas_pfasgroups.png", 
            "detailed_sankey_atlas_pfasgroups.pdf",
            "detailed_sankey_atlas_pfasgroups.svg",
            "detailed_sankey_statistics.csv"
        ]
    }
    
    with open('sankey_metadata.json', 'w') as f:
        json.dump(metadata, f, indent=2)
    
    print(f"📝 Metadata saved: sankey_metadata.json")
    return metadata

def main():
    """Main function to create all Sankey outputs"""
    
    print("🚀 Creating detailed Sankey diagram in multiple formats...")
    print("=" * 60)
    
    try:
        # Create and save Sankey in all formats
        fig, html_file, png_file, pdf_file, svg_file = save_multiple_formats()
        
        # Generate statistics
        stats_df = create_correspondence_statistics()
        
        # Create metadata
        metadata = create_sankey_metadata()
        
        print(f"\n🎉 SANKEY DIAGRAM CREATION COMPLETE!")
        print("=" * 60)
        print(f"Generated files:")
        print(f"  📱 Interactive: {html_file}")
        
        if png_file:
            print(f"  🖼️  High-res PNG: {png_file}")
        if pdf_file:
            print(f"  📄 PDF: {pdf_file}")
        if svg_file:
            print(f"  🎨 Vector SVG: {svg_file}")
            
        print(f"  📊 Statistics: detailed_sankey_statistics.csv")
        print(f"  📝 Metadata: sankey_metadata.json")
        
        if png_file:  # If image exports worked
            print(f"\n💡 Usage recommendations:")
            print(f"  • Use HTML file for interactive exploration")
            print(f"  • Use PNG for presentations and reports")
            print(f"  • Use PDF for publications") 
            print(f"  • Use SVG for web embedding and scaling")
        else:
            print(f"\n💡 Usage recommendations:")
            print(f"  • HTML file provides full interactive experience")
            print(f"  • For static images, consider screenshot from browser or install kaleido")
        
    except Exception as e:
        print(f"❌ Error creating Sankey diagram: {e}")
        import traceback
        traceback.print_exc()

if __name__ == "__main__":
    main()