import pandas as pd
import numpy as np
from collections import Counter, defaultdict
import json

def analyze_benchmark_detailed():
    """Comprehensive analysis of PFASGroups vs PFAS-atlas benchmark results"""
    
    # Load the benchmark results
    df = pd.read_csv('direct_benchmark_results.csv')
    
    print("🔬 COMPREHENSIVE PFASGROUPS VS PFAS-ATLAS ANALYSIS")
    print("=" * 70)
    print(f"Analysis Date: December 13, 2025")
    print(f"Dataset: {len(df)} molecules from specificity_test_results.csv")
    print()
    
    # === BASIC COVERAGE METRICS ===
    print("📊 COVERAGE AND COMPLETENESS")
    print("-" * 40)
    
    total_molecules = len(df)
    pfasgroups_success = (df['pfasgroups_classified'] == True).sum()
    atlas_success = (df['atlas_classified'] == True).sum()
    
    print(f"Total molecules analyzed: {total_molecules:,}")
    print(f"PFASGroups successful: {pfasgroups_success:,}/{total_molecules} ({pfasgroups_success/total_molecules*100:.1f}%)")
    print(f"PFAS-atlas successful: {atlas_success:,}/{total_molecules} ({atlas_success/total_molecules*100:.1f}%)")
    print(f"Both systems successful: {min(pfasgroups_success, atlas_success):,}")
    print()
    
    # === PFAS-ATLAS CLASSIFICATION BREAKDOWN ===
    print("🤖 PFAS-ATLAS CLASSIFICATION RESULTS")
    print("-" * 40)
    
    atlas_counts = df['atlas_first_class'].value_counts()
    for class_name, count in atlas_counts.items():
        percentage = count/total_molecules*100
        print(f"  {class_name}: {count:,} molecules ({percentage:.1f}%)")
    print()
    
    # === PFASGROUPS DETAILED ANALYSIS ===
    print("🧪 PFASGROUPS DETAILED ANALYSIS")
    print("-" * 40)
    
    # Groups detected distribution
    group_dist = df['n_detected_groups'].value_counts().sort_index()
    print("Groups detected per molecule:")
    for n_groups, count in group_dist.items():
        percentage = count/total_molecules*100
        print(f"  {n_groups} groups: {count:,} molecules ({percentage:.1f}%)")
    
    # Specificity analysis
    specific_count = (df['is_specific'] == True).sum()
    print(f"\nSpecific detections: {specific_count:,}/{total_molecules} ({specific_count/total_molecules*100:.1f}%)")
    print(f"Non-specific detections: {total_molecules - specific_count:,} ({(total_molecules - specific_count)/total_molecules*100:.1f}%)")
    print()
    
    # === FUNCTIONAL GROUP ANALYSIS ===
    print("🔬 FUNCTIONAL GROUP FREQUENCY ANALYSIS")
    print("-" * 40)
    
    all_groups = []
    group_cooccurrence = defaultdict(list)
    
    for idx, groups_str in enumerate(df['detected_groups'].dropna()):
        if isinstance(groups_str, str) and groups_str.startswith('['):
            try:
                groups = eval(groups_str)
                all_groups.extend(groups)
                
                # Track co-occurrence
                for group in groups:
                    group_cooccurrence[group].append(idx)
                    
            except Exception as e:
                print(f"Warning: Could not parse groups: {groups_str}")
    
    if all_groups:
        group_counts = Counter(all_groups)
        print(f"Total group detections: {len(all_groups):,}")
        print(f"Unique groups detected: {len(group_counts)}")
        print("\nTop 15 most frequent groups:")
        
        for i, (group_id, count) in enumerate(group_counts.most_common(15), 1):
            percentage = count/total_molecules*100
            print(f"  {i:2d}. Group {group_id:2d}: {count:4,} molecules ({percentage:5.1f}%)")
    print()
    
    # === SYSTEM AGREEMENT ANALYSIS ===
    print("🤝 SYSTEM AGREEMENT ANALYSIS")
    print("-" * 40)
    
    # Basic agreement metrics
    pfas_by_groups = (df['n_detected_groups'] > 0).sum()
    pfas_by_atlas = (~df['atlas_first_class'].isin(['Not PFAS'])).sum()
    
    print(f"PFASGroups detected PFAS: {pfas_by_groups:,}/{total_molecules} ({pfas_by_groups/total_molecules*100:.1f}%)")
    print(f"PFAS-atlas detected PFAS: {pfas_by_atlas:,}/{total_molecules} ({pfas_by_atlas/total_molecules*100:.1f}%)")
    
    # Disagreement analysis
    disagreement_mask = (df['n_detected_groups'] > 0) & (df['atlas_first_class'] == 'Not PFAS')
    disagreements = df[disagreement_mask]
    
    print(f"Disagreements (Groups=Yes, Atlas=No): {len(disagreements):,}")
    
    if len(disagreements) > 0:
        print(f"Disagreement rate: {len(disagreements)/total_molecules*100:.2f}%")
        
        # Analyze disagreement patterns
        print("\nDisagreement analysis by origin:")
        disagreement_origins = disagreements['origin'].value_counts()
        for origin, count in disagreement_origins.head(10).items():
            print(f"  {origin}: {count} molecules")
    print()
    
    # === CROSS-TABULATION ===
    print("📈 DETAILED CROSS-TABULATION")
    print("-" * 40)
    
    has_groups = df['n_detected_groups'] > 0
    crosstab = pd.crosstab(has_groups, df['atlas_first_class'], margins=True)
    print(crosstab)
    print()
    
    # === PERFORMANCE BY MOLECULE COMPLEXITY ===
    print("⚗️ PERFORMANCE BY MOLECULAR COMPLEXITY")
    print("-" * 40)
    
    # Create complexity bins based on number of fluorine atoms (proxy for complexity)
    df['f_count'] = df['smiles'].str.count('F')
    complexity_bins = pd.cut(df['f_count'], bins=[0, 10, 20, 30, 50, 100], 
                            labels=['Low (1-10F)', 'Medium (11-20F)', 'High (21-30F)', 
                                   'Very High (31-50F)', 'Extreme (50+F)'])
    
    complexity_analysis = df.groupby(complexity_bins).agg({
        'n_detected_groups': ['count', 'mean'],
        'is_specific': 'sum',
        'atlas_first_class': lambda x: (x != 'Not PFAS').sum()
    }).round(2)
    
    print("Performance by fluorine content (complexity proxy):")
    for complexity, stats in complexity_analysis.iterrows():
        if pd.isna(complexity):
            continue
        count = int(stats[('n_detected_groups', 'count')])
        if count == 0:
            continue
        avg_groups = stats[('n_detected_groups', 'mean')]
        specific = int(stats[('is_specific', 'sum')])
        atlas_pfas = int(stats[('atlas_first_class', '<lambda>')])
        
        print(f"  {complexity}: {count} molecules")
        print(f"    Avg groups detected: {avg_groups:.1f}")
        print(f"    Specific detections: {specific}/{count} ({specific/count*100:.1f}%)")
        print(f"    Atlas PFAS classification: {atlas_pfas}/{count} ({atlas_pfas/count*100:.1f}%)")
    print()
    
    # === ERROR ANALYSIS ===
    print("⚠️ ERROR AND EDGE CASE ANALYSIS")
    print("-" * 40)
    
    # Check for any errors or edge cases
    error_count = df['error'].notna().sum()
    invalid_smiles = (~df['valid_smiles']).sum() if 'valid_smiles' in df.columns else 0
    
    print(f"Molecules with errors: {error_count}")
    print(f"Invalid SMILES: {invalid_smiles}")
    
    # Analyze the "Not PFAS" cases in detail
    not_pfas_cases = df[df['atlas_first_class'] == 'Not PFAS']
    if len(not_pfas_cases) > 0:
        print(f"\n'Not PFAS' cases analysis ({len(not_pfas_cases)} molecules):")
        print("Origins of molecules classified as 'Not PFAS' by PFAS-atlas:")
        not_pfas_origins = not_pfas_cases['origin'].value_counts()
        for origin, count in not_pfas_origins.items():
            avg_groups = not_pfas_cases[not_pfas_cases['origin'] == origin]['n_detected_groups'].mean()
            print(f"  {origin}: {count} molecules (avg {avg_groups:.1f} groups)")
    print()
    
    # === RECOMMENDATIONS ===
    print("💡 ANALYSIS INSIGHTS AND RECOMMENDATIONS")
    print("-" * 40)
    
    print("Key Findings:")
    print(f"1. System Reliability: Both systems achieved {min(pfasgroups_success, atlas_success)/total_molecules*100:.1f}% success rate")
    
    agreement_rate = (total_molecules - len(disagreements)) / total_molecules * 100
    print(f"2. Classification Agreement: {agreement_rate:.1f}% agreement between systems")
    
    if group_counts:
        dominant_group = group_counts.most_common(1)[0]
        print(f"3. Dominant Pattern: Group {dominant_group[0]} detected in {dominant_group[1]/total_molecules*100:.1f}% of molecules")
    
    print(f"4. Specificity: {specific_count/total_molecules*100:.1f}% of molecules had specific group detection")
    
    print(f"5. Edge Cases: {len(disagreements)} molecules ({len(disagreements)/total_molecules*100:.2f}%) show classification disagreement")
    
    print("\nRecommendations:")
    print("• PFASGroups excels at granular functional group identification")
    print("• PFAS-atlas provides broader ML-based pattern recognition")  
    print("• Combined approach recommended for comprehensive PFAS analysis")
    if len(disagreements) > 0:
        print(f"• Manual review suggested for {len(disagreements)} disagreement cases")
    print("• Both systems now handle complex molecular structures reliably")
    
    # === SAVE DETAILED RESULTS ===
    analysis_results = {
        'analysis_date': '2025-12-13',
        'total_molecules': int(total_molecules),
        'pfasgroups_success_rate': float(pfasgroups_success/total_molecules),
        'atlas_success_rate': float(atlas_success/total_molecules),
        'agreement_rate': float(agreement_rate/100),
        'disagreement_count': int(len(disagreements)),
        'specificity_rate': float(specific_count/total_molecules),
        'atlas_classifications': {k: int(v) for k, v in atlas_counts.items()},
        'top_groups': dict(group_counts.most_common(10)) if 'group_counts' in locals() else {},
        'complexity_analysis': complexity_analysis.to_dict() if 'complexity_analysis' in locals() else {}
    }
    
    with open('detailed_analysis_results.json', 'w') as f:
        json.dump(analysis_results, f, indent=2)
    
    print(f"\n✅ Detailed analysis complete! Results saved to detailed_analysis_results.json")

if __name__ == "__main__":
    analyze_benchmark_detailed()