#!/usr/bin/env python3
"""Analyze HalogenGroups results from clinventory database.

Creates statistical analysis, visualizations, and LaTeX output for groups with id > 28.
"""

import os
import psycopg2
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from collections import Counter
import json

# Set plot style
sns.set_style("whitegrid")
plt.rcParams['figure.figsize'] = (12, 8)
plt.rcParams['font.size'] = 10

def connect_db():
    """Connect to clinventory database."""
    return psycopg2.connect(
        dbname='clinventory',
        user='django',
        password=os.environ.get('DJANGO_USER', ''),
        host='localhost',
        port=5432
    )

def fetch_data(conn):
    """Fetch all groups results from new database schema."""
    # Join HalogenGroups_in_molecules with components_in_molecules to get full details
    query = """
    SELECT 
        pm.molecule_id,
        pm.smiles,
        pm.group_id,
        pm.group_name,
        pm.match_count,
        COUNT(cm.id) as num_components,
        STRING_AGG(DISTINCT cm.smarts_label, '; ') as component_smarts,
        ARRAY_AGG(LENGTH(cm.component_atoms) - LENGTH(REPLACE(cm.component_atoms, ',', '')) + 1) as component_sizes
    FROM HalogenGroups_in_molecules pm
    LEFT JOIN components_in_molecules cm ON pm.smiles = cm.smiles AND pm.group_id = cm.group_id
    WHERE pm.group_id IS NOT NULL
    GROUP BY pm.molecule_id, pm.smiles, pm.group_id, pm.group_name, pm.match_count
    ORDER BY pm.group_id, pm.molecule_id
    """
    return pd.read_sql_query(query, conn)

def classify_group(row):
    """Classify group saturation (per/poly) and halogen type (F/Cl/Br/I/Mixed).
    
    Saturation:
    - 'Per' = fully saturated (all H replaced by halogens)
    - 'Poly' = partially saturated (some H replaced by halogens)
    
    Halogen:
    - 'F' = contains fluorine
    - 'Cl' = contains chlorine
    - 'Br' = contains bromine
    - 'I' = contains iodine
    - 'Mixed' = contains multiple halogen types
    - 'Unknown' = cannot determine
    """
    group_name = str(row['group_name']).lower()
    comp_smarts = str(row['component_smarts']) if pd.notna(row['component_smarts']) else ''
    
    # Saturation: per vs poly (independent of halogen type)
    # Per = all carbons fully halogenated (perfluoro, perchloro, etc.)
    # Poly = partial halogenation (polyfluoro, polychloro, etc.)
    if 'per' in group_name and 'poly' not in group_name:
        # Could be perfluoro, perchloro, etc.
        componentSaturation = 'Per'
    elif 'poly' in group_name or 'semi' in group_name:
        componentSaturation = 'Poly'
    else:
        componentSaturation = 'Unknown'
    
    # Halogen type: which halogen(s) are in the component
    # Check component_smarts for element patterns
    halogens_found = set()
    
    # Check for fluorine
    if 'fluor' in group_name or '[F]' in comp_smarts or '#9' in comp_smarts or 'F' in group_name:
        halogens_found.add('F')
    # Check for chlorine
    if 'chlor' in group_name or '[Cl]' in comp_smarts or '#17' in comp_smarts:
        halogens_found.add('Cl')
    # Check for bromine
    if 'brom' in group_name or '[Br]' in comp_smarts or '#35' in comp_smarts:
        halogens_found.add('Br')
    # Check for iodine
    if 'iod' in group_name or '[I]' in comp_smarts or '#53' in comp_smarts:
        halogens_found.add('I')
    
    # Determine halogen classification
    if len(halogens_found) == 0:
        componentHalogen = 'Unknown'
    elif len(halogens_found) == 1:
        componentHalogen = list(halogens_found)[0]
    else:
        componentHalogen = 'Mixed'
    
    return pd.Series({'componentSaturation': componentSaturation, 'componentHalogen': componentHalogen})

def extract_component_sizes(component_sizes_array):
    """Extract component sizes from PostgreSQL array."""
    try:
        if pd.isna(component_sizes_array):
            return []
        # PostgreSQL returns arrays as strings like '{1,2,3}' or as lists
        if isinstance(component_sizes_array, str):
            # Parse PostgreSQL array format
            sizes_str = component_sizes_array.strip('{}').split(',')
            return [int(s) for s in sizes_str if s and s != 'NULL']
        elif isinstance(component_sizes_array, list):
            return [int(s) for s in component_sizes_array if s is not None]
        else:
            return []
    except:
        return []

def analyze_groups(df, focus_groups_only=True):
    """Analyze groups focusing on id > 28."""
    if focus_groups_only:
        df_focus = df[df['group_id'] > 28].copy()
    else:
        df_focus = df.copy()
    
    # Add classifications
    df_focus[['componentSaturation', 'componentHalogen']] = df_focus.apply(classify_group, axis=1)
    
    # Extract component sizes
    df_focus['component_sizes'] = df_focus['component_sizes'].apply(extract_component_sizes)
    df_focus['avg_component_size'] = df_focus['component_sizes'].apply(
        lambda x: np.mean(x) if len(x) > 0 else 0
    )
    df_focus['max_component_size'] = df_focus['component_sizes'].apply(
        lambda x: max(x) if len(x) > 0 else 0
    )
    
    return df_focus

def create_visualizations(df_focus, output_dir='/home/luc/git/classification_article/imgs'):
    """Create visualization plots."""
    os.makedirs(output_dir, exist_ok=True)
    
    # 1. Distribution of groups
    plt.figure(figsize=(14, 6))
    group_counts = df_focus.groupby('group_id')['molecule_id'].nunique().sort_values(ascending=False)
    group_names = [df_focus[df_focus['group_id'] == gid]['group_name'].iloc[0] for gid in group_counts.index]
    
    plt.subplot(1, 2, 1)
    plt.barh(range(len(group_counts)), group_counts.values, color='steelblue')
    plt.yticks(range(len(group_counts)), [f"{gid}: {name[:30]}" for gid, name in zip(group_counts.index, group_names)], fontsize=8)
    plt.xlabel('Number of Molecules')
    plt.title('Groups Distribution (id > 28)')
    plt.tight_layout()
    
    # 2. Saturation vs Halogen Cross-tabulation
    plt.subplot(1, 2, 2)
    unique_matches = df_focus.drop_duplicates(subset=['molecule_id', 'group_id'])
    sat_hal_table = pd.crosstab(unique_matches['componentSaturation'], unique_matches['componentHalogen'])
    
    # Create grouped bar chart
    sat_hal_table.plot(kind='bar', ax=plt.gca())
    plt.xlabel('Saturation (Per=Fully halogenated, Poly=Partially halogenated)')
    plt.ylabel('Number of Group Matches')
    plt.title('Saturation vs. Component Halogen Type')
    plt.legend(title='Halogen Element', bbox_to_anchor=(1.05, 1))
    plt.xticks(rotation=0)
    plt.tight_layout()
    
    plt.savefig(f'{output_dir}/HalogenGroups_distribution.pdf', dpi=300, bbox_inches='tight')
    plt.savefig(f'{output_dir}/HalogenGroups_distribution.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    # 3. Component sizes
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    
    # Component size distribution
    all_sizes = [size for sizes in df_focus['component_sizes'] for size in sizes if size > 0]
    axes[0, 0].hist(all_sizes, bins=50, edgecolor='black', alpha=0.7)
    axes[0, 0].set_xlabel('Component Size')
    axes[0, 0].set_ylabel('Frequency')
    axes[0, 0].set_title('Distribution of Component Sizes')
    axes[0, 0].set_xlim(0, min(50, max(all_sizes) if all_sizes else 50))
    
    # Average component size by group
    avg_by_group = df_focus.groupby('group_id').agg({
        'avg_component_size': 'mean',
        'group_name': 'first'
    }).sort_values('avg_component_size', ascending=False).head(15)
    
    axes[0, 1].barh(range(len(avg_by_group)), avg_by_group['avg_component_size'].values, color='coral')
    axes[0, 1].set_yticks(range(len(avg_by_group)))
    axes[0, 1].set_yticklabels([f"{idx}: {name[:25]}" for idx, name in zip(avg_by_group.index, avg_by_group['group_name'])], fontsize=8)
    axes[0, 1].set_xlabel('Average Component Size')
    axes[0, 1].set_title('Top 15 Groups by Avg Component Size')
    
    # Number of components per match
    axes[1, 0].hist(df_focus['num_components'], bins=20, edgecolor='black', alpha=0.7, color='green')
    axes[1, 0].set_xlabel('Number of Components')
    axes[1, 0].set_ylabel('Frequency')
    axes[1, 0].set_title('Distribution of Component Counts')
    
    # Match count distribution
    axes[1, 1].hist(df_focus['match_count'], bins=30, edgecolor='black', alpha=0.7, color='purple')
    axes[1, 1].set_xlabel('Match Count')
    axes[1, 1].set_ylabel('Frequency')
    axes[1, 1].set_title('Distribution of Match Counts per Group')
    axes[1, 1].set_xlim(0, min(100, df_focus['match_count'].max()))
    
    plt.tight_layout()
    plt.savefig(f'{output_dir}/HalogenGroups_components.pdf', dpi=300, bbox_inches='tight')
    plt.savefig(f'{output_dir}/HalogenGroups_components.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"Plots saved to {output_dir}")

def create_latex_tables(df_focus, output_file='/home/luc/git/classification_article/HalogenGroups_tables.tex'):
    """Create LaTeX tables."""
    latex_content = []
    
    # Table 1: Top groups by frequency
    latex_content.append("\\begin{table}[h]")
    latex_content.append("\\centering")
    latex_content.append("\\caption{Top 15 Groups (ID > 28) by Molecule Frequency}")
    latex_content.append("\\label{tab:pfas_groups_freq}")
    latex_content.append("\\begin{tabular}{llrrr}")
    latex_content.append("\\hline")
    latex_content.append("ID & Group Name & Molecules & Total Matches & Avg Comp. Size \\\\")
    latex_content.append("\\hline")
    
    top_groups = df_focus.groupby('group_id').agg({
        'molecule_id': 'nunique',
        'group_name': 'first',
        'match_count': 'sum',
        'avg_component_size': 'mean'
    }).sort_values('molecule_id', ascending=False).head(15)
    
    for idx, row in top_groups.iterrows():
        name = row['group_name'][:35] + ('...' if len(row['group_name']) > 35 else '')
        name = name.replace('_', '\\_')
        latex_content.append(f"{idx} & {name} & {int(row['molecule_id'])} & {int(row['match_count'])} & {row['avg_component_size']:.1f} \\\\")
    
    latex_content.append("\\hline")
    latex_content.append("\\end{tabular}")
    latex_content.append("\\end{table}")
    latex_content.append("")
    
    # Table 2: Saturation vs Halogen cross-tabulation
    latex_content.append("\\begin{table}[h]")
    latex_content.append("\\centering")
    latex_content.append("\\caption{Group Classification: Saturation (Per/Poly) vs. Halogen Element}")
    latex_content.append("\\label{tab:pfas_saturation_halogen}")
    latex_content.append("\\small")
    
    sat_hal = df_focus.drop_duplicates(subset=['molecule_id', 'group_id']).groupby(['componentSaturation', 'componentHalogen']).size().unstack(fill_value=0)
    
    # Build dynamic column header
    columns = sorted([col for col in sat_hal.columns if col != 'Unknown'])
    if 'Unknown' in sat_hal.columns:
        columns.append('Unknown')
    
    col_spec = 'l' + 'r' * (len(columns) + 1)
    latex_content.append(f"\\begin{{tabular}}{{{col_spec}}}")
    latex_content.append("\\hline")
    
    header = "Saturation"
    for col in columns:
        header += f" & {col}"
    header += " & Total \\\\"
    latex_content.append(header)
    latex_content.append("\\hline")
    
    # Add note about saturation
    latex_content.append("\\multicolumn{" + str(len(columns) + 2) + "}{l}{\\textit{Per = fully halogenated; Poly = partially halogenated}} \\\\")
    latex_content.append("\\hline")
    
    for sat in ['Per', 'Poly', 'Unknown']:
        if sat in sat_hal.index:
            row = sat
            total = 0
            for col in columns:
                val = int(sat_hal.loc[sat, col]) if col in sat_hal.columns else 0
                row += f" & {val}"
                total += val
            row += f" & {total} \\\\"
            latex_content.append(row)
    
    # Add totals row
    total_row = "\\hline\nTotal"
    grand_total = 0
    for col in columns:
        col_sum = int(sat_hal[col].sum()) if col in sat_hal.columns else 0
        total_row += f" & {col_sum}"
        grand_total += col_sum
    total_row += f" & {grand_total} \\\\"
    latex_content.append(total_row)
    
    latex_content.append("\\hline")
    latex_content.append("\\end{tabular}")
    latex_content.append("\\end{table}")
    
    # Write to file
    with open(output_file, 'w') as f:
        f.write('\n'.join(latex_content))
    
    print(f"LaTeX tables saved to {output_file}")
    return '\n'.join(latex_content)

def create_latex_text(df_focus, df_all, output_file='/home/luc/git/classification_article/HalogenGroups_text.tex'):
    """Create LaTeX text describing results."""
    
    # Calculate statistics
    total_molecules = df_all['molecule_id'].nunique()
    molecules_with_groups = df_all[df_all['group_id'].notna()]['molecule_id'].nunique()
    total_groups_found = df_all['group_id'].nunique()
    groups_above_28 = df_focus['group_id'].nunique()
    molecules_in_focus = df_focus['molecule_id'].nunique()
    total_matches_focus = len(df_focus)
    
    top_group = df_focus.groupby('group_id').agg({
        'molecule_id': 'nunique',
        'group_name': 'first'
    }).sort_values('molecule_id', ascending=False).iloc[0]
    
    avg_component_size = df_focus['avg_component_size'].mean()
    max_component_size = df_focus['max_component_size'].max()
    
    # Saturation analysis (Per vs Poly)
    unique_matches = df_focus.drop_duplicates(subset=['molecule_id', 'group_id'])
    per_count = len(unique_matches[unique_matches['componentSaturation'] == 'Per'])
    poly_count = len(unique_matches[unique_matches['componentSaturation'] == 'Poly'])
    
    # Halogen analysis
    f_count = len(unique_matches[unique_matches['componentHalogen'] == 'F'])
    cl_count = len(unique_matches[unique_matches['componentHalogen'] == 'Cl'])
    br_count = len(unique_matches[unique_matches['componentHalogen'] == 'Br'])
    i_count = len(unique_matches[unique_matches['componentHalogen'] == 'I'])
    mixed_count = len(unique_matches[unique_matches['componentHalogen'] == 'Mixed'])
    
    # Cross-tabulation
    sat_hal_crosstab = pd.crosstab(unique_matches['componentSaturation'], unique_matches['componentHalogen'])
    per_f = sat_hal_crosstab.loc['Per', 'F'] if 'Per' in sat_hal_crosstab.index and 'F' in sat_hal_crosstab.columns else 0
    poly_f = sat_hal_crosstab.loc['Poly', 'F'] if 'Poly' in sat_hal_crosstab.index and 'F' in sat_hal_crosstab.columns else 0
    
    text = f"""\\subsection{{HalogenGroups Classification Results}}

    We applied the HalogenGroups algorithm to {total_molecules:,} halogenated molecules from the CL inventory database. The algorithm successfully identified groups in {molecules_with_groups:,} molecules ({molecules_with_groups/total_molecules*100:.1f}\\%), detecting a total of {total_groups_found} different group types.
    
    Focusing on advanced groups (ID > 28), we identified {groups_above_28} distinct groups in {molecules_in_focus:,} molecules, resulting in {total_matches_focus:,} total matches. The most prevalent group was \\textit{{{top_group['group_name'][:50]}}} (ID {top_group.name}), found in {int(top_group['molecule_id'])} molecules.
    
    \\subsubsection{{Structural Characteristics}}
    
    Component size analysis revealed an average component size of {avg_component_size:.1f} atoms, with the largest component containing {int(max_component_size)} atoms. The distribution of component sizes (Figure~\\ref{{fig:HalogenGroups_components}}) shows that most groups consist of small to medium-sized halogenated chains, with a notable tail extending to larger structures.
    
    \\subsubsection{{Saturation and Halogen Composition}}
    
    Classification by saturation patterns revealed {per_count} matches for fully saturated (per-halogenated) groups and {poly_count} for partially saturated (poly-halogenated) groups. The ratio of {poly_count/(per_count+poly_count)*100:.1f}\\% polyfluorinated indicates that most detected PFAS contain partially halogenated carbon chains rather than fully halogenated structures.
    
    Regarding halogen composition, fluorine-containing groups dominate with {f_count} matches ({f_count/len(unique_matches)*100:.1f}\\% of total), followed by chlorine ({cl_count} matches), bromine ({br_count} matches), and iodine ({i_count} matches). Mixed-halogen compounds (containing multiple halogen types) account for {mixed_count} matches. Notably, perfluorinated groups (fully saturated with fluorine) comprise {per_f} matches, while polyfluorinated groups (partially fluorinated) comprise {poly_f} matches, reflecting the diversity of fluorination patterns in the chemical inventory.
    
    The algorithm detected groups across a range of structural complexities, enabling comprehensive characterization of halogenated substances in the inventory beyond simple fluorinated compounds.
    
    See Tables~\\ref{{tab:pfas_groups_freq}} and \\ref{{tab:pfas_saturation_halogen}} for detailed breakdowns of group frequencies and chemical classifications. Figures~\\ref{{fig:HalogenGroups_distribution}} and \\ref{{fig:HalogenGroups_components}} illustrate the distribution patterns and structural characteristics of the classification system.
    """
    
    with open(output_file, 'w') as f:
        f.write(text)
    
    print(f"LaTeX text saved to {output_file}")
    return text

def main():
    """Main analysis function."""
    print("Connecting to database...")
    conn = connect_db()
    
    print("Fetching data...")
    df_all = fetch_data(conn)
    conn.close()
    
    print(f"Total records: {len(df_all)}")
    print(f"Unique molecules: {df_all['molecule_id'].nunique()}")
    print(f"Unique groups: {df_all['group_id'].nunique()}")
    
    print("\nAnalyzing groups with ID > 28...")
    df_focus = analyze_groups(df_all, focus_groups_only=True)
    
    print(f"Focus dataset: {len(df_focus)} records")
    print(f"Groups in focus: {df_focus['group_id'].nunique()}")
    print(f"Molecules in focus: {df_focus['molecule_id'].nunique()}")
    
    print("\nCreating visualizations...")
    create_visualizations(df_focus)
    
    print("\nCreating LaTeX tables...")
    create_latex_tables(df_focus)
    
    print("\nCreating LaTeX text...")
    create_latex_text(df_focus, df_all)
    
    print("\n=== Summary Statistics ===")
    print(f"\nGroups with ID > 28:")
    print(df_focus.groupby('group_id')['group_name'].first().to_string())
    print(f"\nSaturation breakdown:")
    print(df_focus.drop_duplicates(subset=['molecule_id', 'group_id'])['componentSaturation'].value_counts())
    print(f"\nHalogen breakdown:")
    print(df_focus.drop_duplicates(subset=['molecule_id', 'group_id'])['componentHalogen'].value_counts())
    
    print("\n✓ Analysis complete!")

if __name__ == "__main__":
    main()
