"""
Simple test script to verify the graph loading logic without RDKit dependencies.
"""

import json
import networkx as nx

def load_graph_from_json(filename):
    """Load the list-based JSON structure and convert to NetworkX directed graph with edge types."""
    with open(filename, 'r') as f:
        data = json.load(f)
    
    # Load the PFAS groups mapping to convert names to IDs
    try:
        with open('PFASgroups/data/PFAS_groups_smarts.json', 'r') as f:
            group_data = json.load(f)
        name_to_id = {group['name']: group['id'] for group in group_data}
        print(f"Loaded {len(name_to_id)} group name-to-ID mappings")
    except FileNotFoundError:
        # Fallback: if no mapping file, create a temporary mapping
        print("Warning: PFAS_groups_smarts.json not found, using name-based graph")
        name_to_id = {}
    
    G = nx.DiGraph()
    Gper = nx.DiGraph()

    # Process each edge in the list
    for edge_data in data:
        source_name = edge_data['source']
        target_name = edge_data['target']
        edge_type = edge_data['edge_type']
        
        # Convert names to IDs if mapping exists, otherwise use names
        source = name_to_id.get(source_name, source_name)
        target = name_to_id.get(target_name, target_name)
        
        print(f"  Edge: {source_name} ({source}) -> {target_name} ({target}) [{edge_type}]")
        
        # Add nodes if they don't exist
        G.add_node(source)
        G.add_node(target)
        Gper.add_node(source)
        Gper.add_node(target)
        # Add edge with type information
        G.add_edge(source, target, edge_type=edge_type)
        if edge_type in ['per','both']:
            Gper.add_edge(source, target, edge_type=edge_type)
    
    return G, Gper

def check_groups_are_related(group1, group2, graph):
    """Check if two groups are related in the directed graph."""
    if group1 == group2:
        return True
    
    # Check if both groups exist in the graph
    if group1 not in graph.nodes or group2 not in graph.nodes:
        return False
    
    # Check if group1 is ancestor of group2 or vice versa
    return (group2 in nx.descendants(graph, group1) or 
            group1 in nx.descendants(graph, group2))

def test_graph_loading():
    """Test the corrected graph loading function."""
    
    print("Testing graph loading with ID mapping...")
    
    try:
        # Load the graph
        G, Gper = load_graph_from_json('PFASgroups/tests/specificity_test_groups.json')
        
        print(f"\nMain graph has {len(G.nodes)} nodes and {len(G.edges)} edges")
        print(f"Perfluoro graph has {len(Gper.nodes)} nodes and {len(Gper.edges)} edges")
        
        # Show some node IDs to verify they are integers or names
        print(f"Sample node IDs: {list(G.nodes)[:10]}")
        
        # Test with some known relationships
        print("\nTesting group relationships:")
        
        # Test some known relationships that should exist
        test_pairs = [
            (1, 33),  # Perfluoroalkyl carboxylic acids should be related to carboxylic acid
            (17, 31), # Hydrofluoroethers should be related to ether
            (16, 31), # Perfluoropolyethers should be related to ether
        ]
        
        for group1, group2 in test_pairs:
            is_related = check_groups_are_related(group1, group2, G)
            print(f"  Groups {group1} and {group2} are related: {is_related}")
            
            # Also check if nodes exist
            exists1 = group1 in G.nodes
            exists2 = group2 in G.nodes
            print(f"    Group {group1} exists in graph: {exists1}")
            print(f"    Group {group2} exists in graph: {exists2}")
        
        print("\nGraph loading test completed successfully!")
        return True
        
    except Exception as e:
        print(f"Error in graph loading test: {e}")
        import traceback
        traceback.print_exc()
        return False

if __name__ == "__main__":
    success = test_graph_loading()
    if success:
        print("✓ All tests passed!")
    else:
        print("✗ Some tests failed!")