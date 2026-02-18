#!/usr/bin/env python3
"""Analyze HalogenGroups results from clinventory database and generate plots/tables."""

import sys
import psycopg2
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import seaborn as sns
import numpy as np
from pathlib import Path

# Set style
sns.set_style("whitegrid")
plt.rcParams['figure.dpi'] = 300
plt.rcParams['font.size'] = 10

# Output directory
output_dir = Path("/home/luc/git/classification_article/imgs")
output_dir.mkdir(exist_ok=True, parents=True)
tables_dir = Path("/home/luc/git/classification_article/tables")
tables_dir.mkdir(exist_ok=True, parents=True)

# Connect to database
print("Connecting to database...")
conn = psycopg2.connect(
    dbname='clinventory',
    user='luc',
    password='.Lylix1990',
    host='localhost',
    port=5432
)

# Get summary statistics
print("\n" + "="*70)
print("SUMMARY STATISTICS")
print("="*70)

stats = {}

# Total molecules in database
df = pd.read_sql('SELECT COUNT(*) FROM molecules', conn)
stats['total_db_molecules'] = int(df.iloc[0, 0])
print(f"Total molecules in database: {stats['total_db_molecules']:,}")

# Total molecules with PFAS groups
df = pd.read_sql('SELECT COUNT(DISTINCT molecule_id) FROM HalogenGroups_in_molecules WHERE molecule_id IS NOT NULL', conn)
stats['total_molecules'] = int(df.iloc[0, 0])
print(f"Total molecules with PFAS groups: {stats['total_molecules']:,}")

# Total group matches
df = pd.read_sql('SELECT COUNT(*) FROM HalogenGroups_in_molecules', conn)
stats['total_groups'] = int(df.iloc[0, 0])
print(f"Total PFAS group matches: {stats['total_groups']:,}")

# Total components
df = pd.read_sql('SELECT COUNT(*) FROM components_in_molecules', conn)
stats['total_components'] = int(df.iloc[0, 0])
stats['total_components'] = int(df.iloc[0, 0])
print(f"Total fluorinated components: {stats['total_components']:,}")

# Unique group types
df = pd.read_sql('SELECT COUNT(DISTINCT group_id) FROM HalogenGroups_in_molecules', conn)
stats['unique_groups'] = int(df.iloc[0, 0])
print(f"Unique PFAS group types detected: {stats['unique_groups']}")

# Average groups per molecule
stats['avg_groups_per_mol'] = stats['total_groups'] / stats['total_molecules']
print(f"Average groups per molecule: {stats['avg_groups_per_mol']:.2f}")

# Max groups per molecule
df = pd.read_sql('''SELECT MAX(num_groups) FROM (
    SELECT molecule_id, COUNT(DISTINCT group_id) as num_groups
    FROM HalogenGroups_in_molecules
    GROUP BY molecule_id
) t''', conn)
stats['max_groups_per_mol'] = int(df.iloc[0, 0])
print(f"Max groups per molecule: {stats['max_groups_per_mol']}")

# Percentage of database that is PFAS
stats['pfas_percent'] = 100.0 * stats['total_molecules'] / stats['total_db_molecules']
print(f"PFAS fraction of database: {stats['pfas_percent']:.1f}%")

# Top 15 most common PFAS groups
print("\n" + "="*70)
print("TOP 15 MOST COMMON PFAS GROUPS")
print("="*70)

df_top_groups = pd.read_sql('''
    SELECT group_name, group_id, COUNT(*) as count
    FROM HalogenGroups_in_molecules
    GROUP BY group_name, group_id
    ORDER BY count DESC
    LIMIT 15
''', conn)
print(df_top_groups.to_string(index=False))

# Distribution of number of groups per molecule
print("\n" + "="*70)
print("DISTRIBUTION OF GROUPS PER MOLECULE")
print("="*70)

df_groups_per_mol = pd.read_sql('''
    SELECT num_groups, COUNT(*) as num_molecules
    FROM (
        SELECT molecule_id, COUNT(DISTINCT group_id) as num_groups
        FROM HalogenGroups_in_molecules
        GROUP BY molecule_id
    ) t
    GROUP BY num_groups
    ORDER BY num_groups
''', conn)
print(df_groups_per_mol.to_string(index=False))

# OECD vs Generic vs Telomer vs Other groups
print("\n" + "="*70)
print("BREAKDOWN BY GROUP CATEGORY")
print("="*70)

df_categories = pd.read_sql('''
    SELECT
        CASE
            WHEN group_id BETWEEN 1 AND 28 THEN 'OECD'
            WHEN group_id BETWEEN 29 AND 73 THEN 'Generic'
            WHEN group_id BETWEEN 74 AND 116 THEN 'Telomer'
            ELSE 'Other'
        END as category,
        COUNT(*) as count,
        COUNT(DISTINCT molecule_id) as unique_molecules
    FROM HalogenGroups_in_molecules
    GROUP BY category
    ORDER BY count DESC
''', conn)
print(df_categories.to_string(index=False))

# ===========================
# GENERATE PLOTS
# ===========================

print("\n" + "="*70)
print("GENERATING PLOTS")
print("="*70)

# Color palette consistent across plots
PALETTE = {
    'OECD': '#3498db',
    'Generic': '#e74c3c',
    'Telomer': '#2ecc71',
    'Other': '#f39c12',
}
BAR_COLOR = '#2c7bb6'

# Plot 1: Top 15 PFAS groups bar chart
fig, ax = plt.subplots(figsize=(11, 7))
colors_bar = []
for gid in df_top_groups['group_id']:
    if 1 <= gid <= 28:
        colors_bar.append(PALETTE['OECD'])
    elif 29 <= gid <= 73:
        colors_bar.append(PALETTE['Generic'])
    elif 74 <= gid <= 116:
        colors_bar.append(PALETTE['Telomer'])
    else:
        colors_bar.append(PALETTE['Other'])

bars = ax.barh(df_top_groups['group_name'][::-1], df_top_groups['count'][::-1],
               color=colors_bar[::-1], edgecolor='white', linewidth=0.5)
ax.set_xlabel('Number of Molecules', fontsize=12)
ax.set_title('Top 15 Most Common PFAS Groups in Clinventory Database', fontsize=13, fontweight='bold')
ax.grid(axis='x', alpha=0.3)
ax.set_xlim(0, df_top_groups['count'].max() * 1.18)

for bar, val in zip(bars, df_top_groups['count'][::-1]):
    ax.text(val + df_top_groups['count'].max() * 0.01, bar.get_y() + bar.get_height() / 2,
            f'{val:,}', va='center', fontsize=9)

# Legend for categories
from matplotlib.patches import Patch
legend_elements = [Patch(facecolor=PALETTE['OECD'], label='OECD'),
                   Patch(facecolor=PALETTE['Generic'], label='Generic'),
                   Patch(facecolor=PALETTE['Telomer'], label='Telomer'),
                   Patch(facecolor=PALETTE['Other'], label='Other')]
ax.legend(handles=legend_elements, loc='lower right', fontsize=9)

plt.tight_layout()
plt.savefig(output_dir / 'clinventory_top15_groups.pdf', bbox_inches='tight')
plt.savefig(output_dir / 'clinventory_top15_groups.png', bbox_inches='tight', dpi=300)
print(f"✓ Saved: {output_dir / 'clinventory_top15_groups.pdf'}")
plt.close()

# Plot 2: Distribution of groups per molecule
fig, ax = plt.subplots(figsize=(10, 6))
ax.bar(df_groups_per_mol['num_groups'], df_groups_per_mol['num_molecules'],
       color=BAR_COLOR, edgecolor='white', alpha=0.85)
ax.set_xlabel('Number of PFAS Groups per Molecule', fontsize=12)
ax.set_ylabel('Number of Molecules', fontsize=12)
ax.set_title('Distribution of PFAS Groups per Molecule', fontsize=13, fontweight='bold')
ax.grid(axis='y', alpha=0.3)
ax.xaxis.set_major_locator(mticker.MaxNLocator(integer=True))

for x, y in zip(df_groups_per_mol['num_groups'], df_groups_per_mol['num_molecules']):
    if y > 0:
        ax.text(x, y + df_groups_per_mol['num_molecules'].max() * 0.01, f'{y:,}',
                ha='center', va='bottom', fontsize=8)

plt.tight_layout()
plt.savefig(output_dir / 'clinventory_groups_distribution.pdf', bbox_inches='tight')
plt.savefig(output_dir / 'clinventory_groups_distribution.png', bbox_inches='tight', dpi=300)
print(f"✓ Saved: {output_dir / 'clinventory_groups_distribution.pdf'}")
plt.close()

# Plot 3: Category breakdown (horizontal bar chart)
cat_order = ['OECD', 'Telomer', 'Other', 'Generic']
df_cat_sorted = df_categories.set_index('category').reindex(cat_order).reset_index().dropna()

fig, ax = plt.subplots(figsize=(9, 5))
cat_colors = [PALETTE.get(c, '#95a5a6') for c in df_cat_sorted['category']]
bars = ax.barh(df_cat_sorted['category'], df_cat_sorted['count'],
               color=cat_colors, edgecolor='white', linewidth=0.8)
ax.set_xlabel('Number of Group Matches', fontsize=12)
ax.set_title('PFAS Group Category Distribution in Clinventory', fontsize=13, fontweight='bold')
ax.grid(axis='x', alpha=0.3)
ax.set_xlim(0, df_cat_sorted['count'].max() * 1.22)

total = df_categories['count'].sum()
for bar, (_, row) in zip(bars, df_cat_sorted.iterrows()):
    pct = 100 * row['count'] / total
    ax.text(row['count'] + df_cat_sorted['count'].max() * 0.01, bar.get_y() + bar.get_height() / 2,
            f"{row['count']:,}  ({pct:.1f}%)", va='center', fontsize=10)

plt.tight_layout()
plt.savefig(output_dir / 'clinventory_category_bar.pdf', bbox_inches='tight')
plt.savefig(output_dir / 'clinventory_category_bar.png', bbox_inches='tight', dpi=300)
print(f"✓ Saved: {output_dir / 'clinventory_category_bar.pdf'}")
plt.close()

# Plot 4: Log-scale distribution
fig, ax = plt.subplots(figsize=(10, 6))
ax.bar(df_groups_per_mol['num_groups'], df_groups_per_mol['num_molecules'],
       color='mediumpurple', edgecolor='white', alpha=0.85)
ax.set_xlabel('Number of PFAS Groups per Molecule', fontsize=12)
ax.set_ylabel('Number of Molecules (log scale)', fontsize=12)
ax.set_yscale('log')
ax.set_title('Distribution of PFAS Groups per Molecule (Log Scale)', fontsize=13, fontweight='bold')
ax.grid(axis='y', alpha=0.3, which='both')
ax.xaxis.set_major_locator(mticker.MaxNLocator(integer=True))

plt.tight_layout()
plt.savefig(output_dir / 'clinventory_groups_distribution_log.pdf', bbox_inches='tight')
plt.savefig(output_dir / 'clinventory_groups_distribution_log.png', bbox_inches='tight', dpi=300)
print(f"✓ Saved: {output_dir / 'clinventory_groups_distribution_log.pdf'}")
plt.close()

# ===========================
# GENERATE LATEX TABLES
# ===========================

print("\n" + "="*70)
print("GENERATING LATEX TABLES")
print("="*70)

# Table 1: Summary statistics
latex_summary = f"""\\begin{{table}}[h]
\\centering
\\caption{{Summary statistics of PFASGroups analysis on Clinventory database.}}
\\label{{tab:clinventory_summary}}
\\begin{{tabular}}{{lr}}
\\toprule
\\textbf{{Metric}} & \\textbf{{Value}} \\\\
\\midrule
Total molecules in database & {stats['total_db_molecules']:,} \\\\
Molecules with PFAS groups & {stats['total_molecules']:,} ({stats['pfas_percent']:.1f}\\%) \\\\
Total PFAS group matches & {stats['total_groups']:,} \\\\
Total fluorinated components & {stats['total_components']:,} \\\\
Unique PFAS group types detected & {stats['unique_groups']} \\\\
Average PFAS groups per molecule & {stats['avg_groups_per_mol']:.2f} \\\\
Maximum PFAS groups in one molecule & {stats['max_groups_per_mol']} \\\\
\\bottomrule
\\end{{tabular}}
\\end{{table}}
"""

with open(tables_dir / "clinventory_summary.tex", "w") as f:
    f.write(latex_summary)
print(f"✓ Saved: {tables_dir / 'clinventory_summary.tex'}")

# Table 2: Top 15 groups
def escape_latex(s):
    return str(s).replace('_', r'\_').replace('%', r'\%').replace('&', r'\&')

latex_top15 = "\\begin{table}[h]\n\\centering\n"
latex_top15 += "\\caption{Top 15 most common PFAS groups detected in the Clinventory database.}\n"
latex_top15 += "\\label{tab:clinventory_top15}\n"
latex_top15 += "\\begin{tabular}{clr}\n\\toprule\n"
latex_top15 += "\\textbf{ID} & \\textbf{PFAS Group} & \\textbf{Molecules} \\\\\n\\midrule\n"
for _, row in df_top_groups.iterrows():
    latex_top15 += f"{int(row['group_id'])} & {escape_latex(row['group_name'])} & {int(row['count']):,} \\\\\n"
latex_top15 += "\\bottomrule\n\\end{tabular}\n\\end{table}\n"

with open(tables_dir / "clinventory_top15.tex", "w") as f:
    f.write(latex_top15)
print(f"✓ Saved: {tables_dir / 'clinventory_top15.tex'}")

# Table 3: Category breakdown
total_matches = df_categories['count'].sum()
latex_categories = "\\begin{table}[h]\n\\centering\n"
latex_categories += "\\caption{PFAS group category breakdown in the Clinventory database.}\n"
latex_categories += "\\label{tab:clinventory_categories}\n"
latex_categories += "\\begin{tabular}{lrrr}\n\\toprule\n"
latex_categories += "\\textbf{Category} & \\textbf{Group Matches} & \\textbf{Share (\\%)} & \\textbf{Unique Molecules} \\\\\n\\midrule\n"
for _, row in df_categories.iterrows():
    pct = 100 * row['count'] / total_matches
    latex_categories += f"{row['category']} & {int(row['count']):,} & {pct:.1f} & {int(row['unique_molecules']):,} \\\\\n"
latex_categories += f"\\midrule\nTotal & {int(total_matches):,} & 100.0 & {stats['total_molecules']:,} \\\\\n"
latex_categories += "\\bottomrule\n\\end{tabular}\n\\end{table}\n"

with open(tables_dir / "clinventory_categories.tex", "w") as f:
    f.write(latex_categories)
print(f"✓ Saved: {tables_dir / 'clinventory_categories.tex'}")

# Save statistics to JSON for easy access
import json
with open(output_dir.parent / "clinventory_stats.json", "w") as f:
    json.dump(stats, f, indent=2)
print(f"✓ Saved: {output_dir.parent / 'clinventory_stats.json'}")

conn.close()

print("\n" + "="*70)
print("ANALYSIS COMPLETE!")
print("="*70)
print(f"\nGenerated files:")
print(f"  - Plots: {output_dir}")
print(f"  - Tables: {tables_dir}")
print(f"  - Stats: {output_dir.parent / 'clinventory_stats.json'}")


# Connect to database
print("Connecting to database...")
conn = psycopg2.connect(
    dbname='clinventory',
    user='luc',
    password='.Lylix1990',
    host='localhost',
    port=5432
)

# Get summary statistics
print("\n" + "="*70)
print("SUMMARY STATISTICS")
print("="*70)

stats = {}

# Total molecules with PFAS groups
df = pd.read_sql('SELECT COUNT(DISTINCT molecule_id) FROM HalogenGroups_in_molecules', conn)
stats['total_molecules'] = int(df.iloc[0, 0])
print(f"Total molecules with PFAS groups: {stats['total_molecules']:,}")

# Total group matches
df = pd.read_sql('SELECT COUNT(*) FROM HalogenGroups_in_molecules', conn)
stats['total_groups'] = int(df.iloc[0, 0])
print(f"Total PFAS group matches: {stats['total_groups']:,}")

# Total components
df = pd.read_sql('SELECT COUNT(*) FROM components_in_molecules', conn)
stats['total_components'] = int(df.iloc[0, 0])
print(f"Total fluorinated components: {stats['total_components']:,}")

# Unique group types
df = pd.read_sql('SELECT COUNT(DISTINCT group_id) FROM HalogenGroups_in_molecules', conn)
stats['unique_groups'] = int(df.iloc[0, 0])
print(f"Unique PFAS group types detected: {stats['unique_groups']}")

# Average groups per molecule
stats['avg_groups_per_mol'] = stats['total_groups'] / stats['total_molecules']
print(f"Average groups per molecule: {stats['avg_groups_per_mol']:.2f}")

# Top 15 most common PFAS groups
print("\n" + "="*70)
print("TOP 15 MOST COMMON PFAS GROUPS")
print("="*70)

df_top_groups = pd.read_sql('''
    SELECT group_name, COUNT(*) as count 
    FROM HalogenGroups_in_molecules 
    GROUP BY group_name 
    ORDER BY count DESC 
    LIMIT 15
''', conn)
print(df_top_groups.to_string(index=False))

# Distribution of number of groups per molecule
print("\n" + "="*70)
print("DISTRIBUTION OF GROUPS PER MOLECULE")
print("="*70)

df_groups_per_mol = pd.read_sql('''
    SELECT num_groups, COUNT(*) as num_molecules
    FROM (
        SELECT molecule_id, COUNT(DISTINCT group_id) as num_groups
        FROM HalogenGroups_in_molecules
        GROUP BY molecule_id
    ) t
    GROUP BY num_groups
    ORDER BY num_groups
''', conn)
print(df_groups_per_mol.to_string(index=False))

# OECD vs Generic vs Telomer groups
print("\n" + "="*70)
print("BREAKDOWN BY GROUP CATEGORY")
print("="*70)

df_categories = pd.read_sql('''
    SELECT 
        CASE 
            WHEN group_id BETWEEN 1 AND 28 THEN 'OECD'
            WHEN group_id BETWEEN 29 AND 73 THEN 'Generic'
            WHEN group_id BETWEEN 74 AND 116 THEN 'Telomer'
            ELSE 'Other'
        END as category,
        COUNT(*) as count,
        COUNT(DISTINCT molecule_id) as unique_molecules
    FROM HalogenGroups_in_molecules
    GROUP BY category
    ORDER BY count DESC
''', conn)
print(df_categories.to_string(index=False))

# ===========================
# GENERATE PLOTS
# ===========================

print("\n" + "="*70)
print("GENERATING PLOTS")
print("="*70)

# Plot 1: Top 15 PFAS groups bar chart
fig, ax = plt.subplots(figsize=(10, 6))
bars = ax.barh(df_top_groups['group_name'][::-1], df_top_groups['count'][::-1], color='steelblue')
ax.set_xlabel('Number of Molecules', fontsize=12)
ax.set_ylabel('PFAS Group', fontsize=12)
ax.set_title('Top 15 Most Common PFAS Groups in Clinventory Database', fontsize=14, fontweight='bold')
ax.grid(axis='x', alpha=0.3)

# Add value labels
for i, (bar, val) in enumerate(zip(bars, df_top_groups['count'][::-1])):
    ax.text(val + max(df_top_groups['count'])*0.01, bar.get_y() + bar.get_height()/2, 
            f'{val:,}', va='center', fontsize=9)

plt.tight_layout()
plt.savefig(output_dir / 'clinventory_top15_groups.pdf', bbox_inches='tight')
plt.savefig(output_dir / 'clinventory_top15_groups.png', bbox_inches='tight', dpi=300)
print(f"✓ Saved: {output_dir / 'clinventory_top15_groups.pdf'}")
plt.close()

# Plot 2: Distribution of groups per molecule
fig, ax = plt.subplots(figsize=(10, 6))
ax.bar(df_groups_per_mol['num_groups'], df_groups_per_mol['num_molecules'], 
       color='coral', edgecolor='darkred', alpha=0.7)
ax.set_xlabel('Number of PFAS Groups per Molecule', fontsize=12)
ax.set_ylabel('Number of Molecules', fontsize=12)
ax.set_title('Distribution of PFAS Groups per Molecule', fontsize=14, fontweight='bold')
ax.grid(axis='y', alpha=0.3)

# Add value labels
for x, y in zip(df_groups_per_mol['num_groups'], df_groups_per_mol['num_molecules']):
    if y > 0:
        ax.text(x, y + max(df_groups_per_mol['num_molecules'])*0.01, f'{y:,}', 
                ha='center', va='bottom', fontsize=9)

plt.tight_layout()
plt.savefig(output_dir / 'clinventory_groups_distribution.pdf', bbox_inches='tight')
plt.savefig(output_dir / 'clinventory_groups_distribution.png', bbox_inches='tight', dpi=300)
print(f"✓ Saved: {output_dir / 'clinventory_groups_distribution.pdf'}")
plt.close()

# Plot 3: Category breakdown (horizontal bar chart)
fig, ax = plt.subplots(figsize=(10, 6))
colors = ['#3498db', '#e74c3c', '#2ecc71']

# Sort by count for better visualization
df_cat_sorted = df_categories.sort_values('count', ascending=True)

bars = ax.barh(df_cat_sorted['category'], df_cat_sorted['count'], 
               color=colors[:len(df_cat_sorted)], edgecolor='black', alpha=0.8, linewidth=1.5)
ax.set_xlabel('Number of Group Matches', fontsize=12, fontweight='bold')
ax.set_ylabel('Category', fontsize=12, fontweight='bold')
ax.set_title('PFAS Group Categories Distribution', fontsize=14, fontweight='bold')
ax.grid(axis='x', alpha=0.3)

# Add value labels and percentages
total = df_categories['count'].sum()
for i, (bar, val) in enumerate(zip(bars, df_cat_sorted['count'])):
    percentage = (val / total) * 100
    ax.text(val + max(df_cat_sorted['count'])*0.01, bar.get_y() + bar.get_height()/2,
            f"{val:,} ({percentage:.1f}%)", va='center', fontsize=11, fontweight='bold')

plt.tight_layout()
plt.savefig(output_dir / 'clinventory_category_bar.pdf', bbox_inches='tight')
plt.savefig(output_dir / 'clinventory_category_bar.png', bbox_inches='tight', dpi=300)
print(f"✓ Saved: {output_dir / 'clinventory_category_bar.pdf'}")
plt.close()

# Plot 4: Log-scale distribution for better visibility
if len(df_groups_per_mol) > 0:
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.bar(df_groups_per_mol['num_groups'], df_groups_per_mol['num_molecules'], 
           color='mediumpurple', edgecolor='darkviolet', alpha=0.7)
    ax.set_xlabel('Number of PFAS Groups per Molecule', fontsize=12)
    ax.set_ylabel('Number of Molecules (log scale)', fontsize=12)
    ax.set_yscale('log')
    ax.set_title('Distribution of PFAS Groups per Molecule (Log Scale)', fontsize=14, fontweight='bold')
    ax.grid(axis='y', alpha=0.3, which='both')
    
    plt.tight_layout()
    plt.savefig(output_dir / 'clinventory_groups_distribution_log.pdf', bbox_inches='tight')
    plt.savefig(output_dir / 'clinventory_groups_distribution_log.png', bbox_inches='tight', dpi=300)
    print(f"✓ Saved: {output_dir / 'clinventory_groups_distribution_log.pdf'}")
    plt.close()

# ===========================
# GENERATE LATEX TABLES
# ===========================

print("\n" + "="*70)
print("GENERATING LATEX TABLES")
print("="*70)

# Table 1: Summary statistics
latex_summary = f"""\\begin{{table}}[h]
\\centering
\\caption{{Summary statistics of HalogenGroups analysis on Clinventory database}}
\\label{{tab:clinventory_summary}}
\\begin{{tabular}}{{lr}}
\\toprule
\\textbf{{Metric}} & \\textbf{{Value}} \\\\
\\midrule
Total molecules with PFAS groups & {stats['total_molecules']:,} \\\\
Total PFAS group matches & {stats['total_groups']:,} \\\\
Total fluorinated components & {stats['total_components']:,} \\\\
Unique PFAS group types detected & {stats['unique_groups']} \\\\
Average groups per molecule & {stats['avg_groups_per_mol']:.2f} \\\\
\\bottomrule
\\end{{tabular}}
\\end{{table}}
"""

with open(output_dir.parent / "tables" / "clinventory_summary.tex", "w") as f:
    f.write(latex_summary)
print(f"✓ Saved: {output_dir.parent / 'tables' / 'clinventory_summary.tex'}")

# Table 2: Top 15 groups
latex_top15 = "\\begin{table}[h]\n\\centering\n"
latex_top15 += "\\caption{Top 15 most common PFAS groups in Clinventory database}\n"
latex_top15 += "\\label{tab:clinventory_top15}\n"
latex_top15 += "\\begin{tabular}{lr}\n\\toprule\n"
latex_top15 += "\\textbf{PFAS Group} & \\textbf{Count} \\\\\n\\midrule\n"
for _, row in df_top_groups.iterrows():
    latex_top15 += f"{row['group_name']} & {row['count']:,} \\\\\n"
latex_top15 += "\\bottomrule\n\\end{tabular}\n\\end{table}\n"

with open(output_dir.parent / "tables" / "clinventory_top15.tex", "w") as f:
    f.write(latex_top15)
print(f"✓ Saved: {output_dir.parent / 'tables' / 'clinventory_top15.tex'}")

# Table 3: Category breakdown
latex_categories = "\\begin{table}[h]\n\\centering\n"
latex_categories += "\\caption{PFAS group category breakdown in Clinventory database}\n"
latex_categories += "\\label{tab:clinventory_categories}\n"
latex_categories += "\\begin{tabular}{lrr}\n\\toprule\n"
latex_categories += "\\textbf{Category} & \\textbf{Total Matches} & \\textbf{Unique Molecules} \\\\\n\\midrule\n"
for _, row in df_categories.iterrows():
    latex_categories += f"{row['category']} & {row['count']:,} & {row['unique_molecules']:,} \\\\\n"
latex_categories += "\\bottomrule\n\\end{tabular}\n\\end{table}\n"

with open(output_dir.parent / "tables" / "clinventory_categories.tex", "w") as f:
    f.write(latex_categories)
print(f"✓ Saved: {output_dir.parent / 'tables' / 'clinventory_categories.tex'}")

# Save statistics to JSON for easy access
import json
with open(output_dir.parent / "clinventory_stats.json", "w") as f:
    json.dump(stats, f, indent=2)
print(f"✓ Saved: {output_dir.parent / 'clinventory_stats.json'}")

conn.close()

print("\n" + "="*70)
print("ANALYSIS COMPLETE!")
print("="*70)
print(f"\nGenerated files:")
print(f"  - Plots: {output_dir}")
print(f"  - Tables: {output_dir.parent / 'tables'}")
print(f"  - Stats: {output_dir.parent / 'clinventory_stats.json'}")
