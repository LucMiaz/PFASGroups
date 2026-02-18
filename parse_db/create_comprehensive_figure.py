#!/usr/bin/env python3
"""Create a comprehensive multi-panel figure for the Clinventory analysis."""

import psycopg2
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from pathlib import Path

# Set style
sns.set_style("whitegrid")
plt.rcParams['figure.dpi'] = 300
plt.rcParams['font.size'] = 9

# Output
output_dir = Path("/home/luc/git/classification_article/imgs")

# Connect to database
conn = psycopg2.connect(
    dbname='clinventory',
    user='luc',
    password='.Lylix1990',
    host='localhost',
    port=5432
)

# Get data
df_top_groups = pd.read_sql('''
    SELECT group_name, COUNT(*) as count 
    FROM HalogenGroups_in_molecules 
    GROUP BY group_name 
    ORDER BY count DESC 
    LIMIT 10
''', conn)

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

df_categories = pd.read_sql('''
    SELECT 
        CASE 
            WHEN group_id BETWEEN 1 AND 28 THEN 'OECD'
            WHEN group_id BETWEEN 29 AND 73 THEN 'Generic'
            WHEN group_id BETWEEN 74 AND 116 THEN 'Telomer'
        END as category,
        COUNT(*) as count
    FROM HalogenGroups_in_molecules
    GROUP BY category
''', conn)

# Component size distribution
df_component_size = pd.read_sql('''
    SELECT 
        LENGTH(component_atoms) - LENGTH(REPLACE(component_atoms, ',', '')) + 1 as component_size,
        COUNT(*) as count
    FROM components_in_molecules
    GROUP BY component_size
    ORDER BY component_size
''', conn)

conn.close()

# Create comprehensive 2x2 figure
fig = plt.figure(figsize=(14, 10))
gs = fig.add_gridspec(2, 2, hspace=0.3, wspace=0.3)

# Panel A: Top 10 PFAS groups
ax1 = fig.add_subplot(gs[0, 0])
y_pos = np.arange(len(df_top_groups))
bars = ax1.barh(y_pos, df_top_groups['count'], color='steelblue', edgecolor='darkblue', alpha=0.8)
ax1.set_yticks(y_pos)
ax1.set_yticklabels(df_top_groups['group_name'], fontsize=8)
ax1.set_xlabel('Number of Molecules', fontsize=10, fontweight='bold')
ax1.set_title('(A) Top 10 PFAS Groups', fontsize=11, fontweight='bold', loc='left')
ax1.grid(axis='x', alpha=0.3)
ax1.invert_yaxis()

# Add counts
for i, (bar, val) in enumerate(zip(bars, df_top_groups['count'])):
    ax1.text(val + max(df_top_groups['count'])*0.015, bar.get_y() + bar.get_height()/2,
            f'{val:,}', va='center', fontsize=7, fontweight='bold')

# Panel B: Category horizontal bar chart
ax2 = fig.add_subplot(gs[0, 1])
colors = ['#3498db', '#e74c3c', '#2ecc71']

# Sort by count
df_cat_sorted = df_categories.sort_values('count', ascending=True)

bars = ax2.barh(df_cat_sorted['category'], df_cat_sorted['count'],
                color=colors[:len(df_cat_sorted)], edgecolor='black', alpha=0.8, linewidth=1.2)
ax2.set_xlabel('Number of Matches', fontsize=10, fontweight='bold')
ax2.set_title('(B) Category Distribution', fontsize=11, fontweight='bold', loc='left')
ax2.grid(axis='x', alpha=0.3)

# Add value labels with percentages
total = df_categories['count'].sum()
for bar, val in zip(bars, df_cat_sorted['count']):
    percentage = (val / total) * 100
    ax2.text(val + max(df_cat_sorted['count'])*0.015, bar.get_y() + bar.get_height()/2,
            f"{val:,}\n({percentage:.1f}%)", va='center', ha='left', fontsize=8, fontweight='bold')

# Panel C: Groups per molecule distribution
ax3 = fig.add_subplot(gs[1, 0])
bars = ax3.bar(df_groups_per_mol['num_groups'], 
               df_groups_per_mol['num_molecules'],
               color='coral', edgecolor='darkred', alpha=0.7, width=0.7)
ax3.set_xlabel('Number of PFAS Groups per Molecule', fontsize=10, fontweight='bold')
ax3.set_ylabel('Number of Molecules', fontsize=10, fontweight='bold')
ax3.set_title('(C) Distribution of Groups per Molecule', fontsize=11, fontweight='bold', loc='left')
ax3.grid(axis='y', alpha=0.3)
ax3.set_xticks(df_groups_per_mol['num_groups'])

# Add counts
for bar, val in zip(bars, df_groups_per_mol['num_molecules']):
    if val > 0:
        ax3.text(bar.get_x() + bar.get_width()/2, val + max(df_groups_per_mol['num_molecules'])*0.01,
                f'{val:,}', ha='center', va='bottom', fontsize=7, fontweight='bold')

# Panel D: Component size distribution
ax4 = fig.add_subplot(gs[1, 1])
# Limit to reasonable range
df_comp_plot = df_component_size[df_component_size['component_size'] <= 25]
bars = ax4.bar(df_comp_plot['component_size'], 
               df_comp_plot['count'],
               color='mediumpurple', edgecolor='darkviolet', alpha=0.7, width=0.8)
ax4.set_xlabel('Fluorinated Component Size (atoms)', fontsize=10, fontweight='bold')
ax4.set_ylabel('Number of Components', fontsize=10, fontweight='bold')
ax4.set_title('(D) Component Size Distribution', fontsize=11, fontweight='bold', loc='left')
ax4.grid(axis='y', alpha=0.3)

# Add median line
median_size = 8  # From earlier analysis
ax4.axvline(median_size, color='red', linestyle='--', linewidth=2, alpha=0.7, label=f'Median = {median_size}')
ax4.legend(fontsize=9)

# Overall title
fig.suptitle('Clinventory Database: PFAS Classification Analysis', 
             fontsize=14, fontweight='bold', y=0.98)

# Add summary text box
summary_text = (
    f"Total molecules: 1,436\n"
    f"Total matches: 2,865\n"
    f"Unique groups: 22\n"
    f"Avg groups/mol: 2.00"
)
fig.text(0.99, 0.01, summary_text, 
         ha='right', va='bottom', fontsize=8,
         bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.3))

plt.savefig(output_dir / 'clinventory_comprehensive_analysis.pdf', bbox_inches='tight')
plt.savefig(output_dir / 'clinventory_comprehensive_analysis.png', bbox_inches='tight', dpi=300)
print(f"✓ Saved comprehensive figure to {output_dir}")
plt.close()

print("\n" + "="*70)
print("COMPREHENSIVE FIGURE CREATED")
print("="*70)
print(f"File: clinventory_comprehensive_analysis.pdf")
print(f"Location: {output_dir}")
print(f"\nThis 4-panel figure includes:")
print("  (A) Top 10 PFAS groups - horizontal bar chart")
print("  (B) Category distribution - horizontal bar chart")
print("  (C) Groups per molecule - distribution")
print("  (D) Component size - size distribution")
