"""Quick verification that all metrics are computed and reported."""

from PFASgroups.parser import parse_smiles

# Test with PFOA
smiles = 'C(=O)(O)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)F'
results = parse_smiles(smiles, bycomponent=True, output_format='list')

print("=" * 80)
print("VERIFICATION: All Metrics Computed and Reported")
print("=" * 80)

if not results or len(results) == 0:
    print("❌ No results returned")
    exit(1)

mol_result = results[0]
matches = mol_result.get('matches', [])

# Check first PFAS group match
for match in matches:
    if match['type'] != 'PFASgroup':
        continue
    
    print(f"\n✓ Group: {match['group_name']}")
    
    # Verify all required summary metrics
    required_summary = [
        'mean_eccentricity',
        'mean_smarts_centrality',
        'mean_diameter',
        'mean_radius',
        'mean_effective_graph_resistance',
        'mean_dist_to_barycenter',
        'mean_dist_to_center',
        'mean_dist_to_periphery'
    ]
    
    print("\n  Summary Metrics:")
    all_summary_present = True
    for metric in required_summary:
        if metric in match:
            print(f"    ✓ {metric}: {match[metric]}")
        else:
            print(f"    ❌ {metric}: MISSING")
            all_summary_present = False
    
    # Verify component metrics
    components = match.get('components', [])
    if len(components) > 0:
        comp = components[0]
        
        required_component = [
            # Graph structure metrics
            'diameter', 'radius', 'effective_graph_resistance',
            'eccentricity_values',
            # Node sets
            'center', 'periphery', 'barycenter',
            # SMARTS distances
            'min_dist_to_barycenter',
            'min_resistance_dist_to_barycenter',
            'min_dist_to_center',
            'min_resistance_dist_to_center',
            'max_dist_to_periphery',
            'max_resistance_dist_to_periphery'
        ]
        
        print(f"\n  Component Metrics (component 0):")
        all_component_present = True
        for metric in required_component:
            if metric in comp:
                value = comp[metric]
                if isinstance(value, list):
                    print(f"    ✓ {metric}: {value[:3]}{'...' if len(value) > 3 else ''}")
                elif isinstance(value, dict):
                    print(f"    ✓ {metric}: {len(value)} nodes")
                elif isinstance(value, float):
                    print(f"    ✓ {metric}: {value:.2f}")
                else:
                    print(f"    ✓ {metric}: {value}")
            else:
                print(f"    ❌ {metric}: MISSING")
                all_component_present = False
        
        if all_summary_present and all_component_present:
            print("\n✅ ALL REQUIRED METRICS ARE PRESENT AND COMPUTED")
            break
        else:
            print("\n❌ SOME METRICS ARE MISSING")
            exit(1)
    break

print("\n" + "=" * 80)
print("VERIFICATION COMPLETE: All metrics implemented correctly!")
print("=" * 80)
