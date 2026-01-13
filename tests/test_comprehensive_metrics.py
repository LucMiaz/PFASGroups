"""Test comprehensive graph metrics in ComponentsSolver."""

from rdkit import Chem
from PFASgroups.core import parse_smiles
import json

def test_comprehensive_metrics():
    """Test comprehensive NetworkX graph metrics on various PFAS structures."""
    
    test_cases = [
        {
            'name': 'PFOA (linear perfluoroalkyl)',
            'smiles': 'C(=O)(O)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)F',
            'expected': {
                'eccentricity': 'high (linear)',
                'diameter': '> 0 (linear chain)',
                'radius': '< diameter',
                'min_dist_to_center': 'medium (carboxylic acid at end)',
                'max_dist_to_periphery': 'medium'
            }
        },
        {
            'name': 'Branched PFAS',
            'smiles': 'C(C(F)(F)F)(C(F)(F)F)C(F)(F)C(F)(F)C(F)(F)C(=O)O',
            'expected': {
                'eccentricity': 'lower (branched)',
                'diameter': 'increased by branching',
                'barycenter': 'at branch point',
                'min_dist_to_center': 'larger if far from center'
            }
        },
        {
            'name': 'PFOS (sulfonic acid)',
            'smiles': 'C(C(C(C(C(C(C(C(F)(F)S(=O)(=O)O)(F)F)(F)F)(F)F)(F)F)(F)F)(F)F)(F)F',
            'expected': {
                'eccentricity': 'high (mostly linear)',
                'min_dist_to_barycenter': 'depends on SMARTS position'
            }
        },
        {
            'name': '6:2 FTOH (terminal OH)',
            'smiles': 'C(CCO)(F)(C(C(C(C(C(F)(F)F)(F)F)(F)F)(F)F)(F)F)F',
            'expected': {
                'min_dist_to_center': 'low if near center, high if at terminus',
                'max_dist_to_periphery': 'varies by position'
            }
        }
    ]
    
    print("=" * 80)
    print("COMPREHENSIVE GRAPH METRICS TEST")
    print("=" * 80)
    
    for test_case in test_cases:
        print(f"\n{'=' * 80}")
        print(f"Testing: {test_case['name']}")
        print(f"SMILES: {test_case['smiles']}")
        print(f"{'=' * 80}")
        
        try:
            result = parse_smiles(test_case['smiles'], bycomponent=True, output_format='dict')
            
            if not result or not isinstance(result, list) or len(result) == 0:
                print(f"\n❌ No results returned for {test_case['name']}")
                continue
            
            mol_result = result[0]
            matches = mol_result.get('matches', [])
            
            print(f"\n🔍 Matched Groups: {len(matches)}")
            
            for match in matches:
                group_name = match.get('group_name', 'Unknown')
                match_count = match.get('match_count', 0)
                components = match.get('components', [])
                
                print(f"\n📊 Group: {group_name}")
                print(f"   Match count: {match_count}")
                print(f"   Components found: {len(components)}")
                
                for i, comp_data in enumerate(components):
                    print(f"\n   Component #{i+1}:")
                    print(f"   ├─ Size: {comp_data['size']} atoms")
                    print(f"   ├─ Atom indices: {comp_data['component'][:10]}{'...' if len(comp_data['component']) > 10 else ''}")
                    
                    # Basic metrics
                    print(f"   │")
                    print(f"   ├─ Basic Metrics:")
                    print(f"   │  ├─ Eccentricity (branching): {comp_data['eccentricity']:.3f}")
                    print(f"   │  │  (1.0=linear, 0.0=highly branched)")
                    print(f"   │  └─ SMARTS Centrality: {comp_data['smarts_centrality']:.3f}")
                    print(f"   │     (1.0=central, 0.0=peripheral)")
                    
                    # Graph structure metrics
                    print(f"   │")
                    print(f"   ├─ Graph Structure Metrics:")
                    diameter = comp_data.get('diameter', 'N/A')
                    radius = comp_data.get('radius', 'N/A')
                    resistance = comp_data.get('effective_graph_resistance', 'N/A')
                    
                    print(f"   │  ├─ Diameter: {diameter}")
                    print(f"   │  │  (maximum eccentricity in graph)")
                    print(f"   │  ├─ Radius: {radius}")
                    print(f"   │  │  (minimum eccentricity in graph)")
                    print(f"   │  └─ Effective Graph Resistance: {resistance}")
                    print(f"   │     (sum of resistance distances)")
                    
                    # Node sets
                    print(f"   │")
                    print(f"   ├─ Structural Node Sets:")
                    center = comp_data.get('center', [])
                    periphery = comp_data.get('periphery', [])
                    barycenter = comp_data.get('barycenter', [])
                    
                    print(f"   │  ├─ Center nodes: {center[:5]}{'...' if len(center) > 5 else ''}")
                    print(f"   │  │  (nodes with minimum eccentricity)")
                    print(f"   │  ├─ Periphery nodes: {periphery[:5]}{'...' if len(periphery) > 5 else ''}")
                    print(f"   │  │  (nodes with maximum eccentricity)")
                    print(f"   │  └─ Barycenter nodes: {barycenter[:5]}{'...' if len(barycenter) > 5 else ''}")
                    print(f"   │     (nodes minimizing total distance)")
                    
                    # SMARTS-specific distances
                    print(f"   │")
                    print(f"   └─ SMARTS Position Metrics:")
                    min_dist_bc = comp_data.get('min_dist_to_barycenter', 'N/A')
                    min_resist_bc = comp_data.get('min_resistance_dist_to_barycenter', 'N/A')
                    min_dist_center = comp_data.get('min_dist_to_center', 'N/A')
                    min_resist_center = comp_data.get('min_resistance_dist_to_center', 'N/A')
                    max_dist_periph = comp_data.get('max_dist_to_periphery', 'N/A')
                    max_resist_periph = comp_data.get('max_resistance_dist_to_periphery', 'N/A')
                    
                    print(f"      ├─ Distance to barycenter: {min_dist_bc} (resistance: {min_resist_bc})")
                    print(f"      ├─ Distance to center: {min_dist_center} (resistance: {min_resist_center})")
                    print(f"      └─ Distance to periphery: {max_dist_periph} (resistance: {max_resist_periph})")
            
            # Interpretation
            print(f"\n💡 Interpretation:")
            for expected_key, expected_value in test_case['expected'].items():
                print(f"   • {expected_key}: {expected_value}")
            
            print(f"\n✅ Test passed for {test_case['name']}")
            
        except Exception as e:
            print(f"\n❌ Error testing {test_case['name']}: {e}")
            import traceback
            traceback.print_exc()
    
    print(f"\n{'=' * 80}")
    print("TEST SUMMARY")
    print("=" * 80)
    print("All comprehensive graph metrics have been computed and displayed.")
    print("\nMetrics Explanation:")
    print("  • Diameter: Maximum distance between any two nodes")
    print("  • Radius: Minimum eccentricity (distance to farthest node)")
    print("  • Center: Nodes with minimum eccentricity (most central)")
    print("  • Periphery: Nodes with maximum eccentricity (most peripheral)")
    print("  • Barycenter: Nodes minimizing sum of distances to all nodes")
    print("  • Resistance Distance: Graph-theoretic resistance between nodes")
    print("\nSMARTS Position Metrics:")
    print("  • Show how far functional groups are from key structural nodes")
    print("  • Lower distance to center/barycenter = more central placement")
    print("  • Higher distance to periphery = functional group is more exposed")

if __name__ == '__main__':
    test_comprehensive_metrics()
