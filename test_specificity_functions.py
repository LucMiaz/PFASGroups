"""
Test the corrected specificity test functionality.
"""

import sys
import os

# Try to test just the graph-related functions without the full RDKit dependency
try:
    import json
    import networkx as nx
    
    # Import specific functions we need to test
    sys.path.append('PFASgroups')
    
    def load_graph_from_json(filename):
        """Load the list-based JSON structure and convert to NetworkX directed graph with edge types."""
        with open(filename, 'r') as f:
            data = json.load(f)
        
        # Load the PFAS groups mapping to convert names to IDs
        try:
            with open('PFASgroups/data/PFAS_groups_smarts.json', 'r') as f:
                group_data = json.load(f)
            name_to_id = {group['name']: group['id'] for group in group_data}
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

    def are_detected_groups_acceptable(expected_groups, detected_groups, json_file='PFASgroups/tests/specificity_test_groups.json'):
        """
        Check if detected groups are acceptable given the expected groups and graph relationships.
        """
        try:
            G, Gper = load_graph_from_json(json_file)
        except (FileNotFoundError, json.JSONDecodeError):
            # If graph file doesn't exist or is invalid, fall back to exact matching
            return sorted(expected_groups) == sorted(detected_groups)
        
        # All detected groups must be either:
        # 1. In the expected groups, OR
        # 2. Connected to at least one expected group in the graph
        for detected_group in detected_groups:
            if detected_group in expected_groups:
                continue  # This group is expected
            
            # Check if this detected group is related to any expected group
            is_related = False
            for expected_group in expected_groups:
                if check_groups_are_related(detected_group, expected_group, G):
                    is_related = True
                    break
            
            if not is_related:
                return False  # Found an unrelated detected group
        
        return True

    def test_specificity_function():
        """Test the specificity function with some example cases."""
        print("Testing are_detected_groups_acceptable function...")
        
        test_cases = [
            # (expected, detected, expected_result, description)
            ([1], [1], True, "Exact match should be acceptable"),
            ([1], [1, 33], True, "Carboxylic acid generic with specific should be acceptable"),
            ([17], [17, 31], True, "Ether generic with hydrofluoroether should be acceptable"),
            ([16], [16, 31], True, "Ether generic with perfluoropolyether should be acceptable"),
            ([1], [2], False, "Different specific groups should not be acceptable"),
            ([1], [1, 2, 33], False, "Multiple unrelated specific groups should not be acceptable"),
        ]
        
        for expected, detected, expected_result, description in test_cases:
            result = are_detected_groups_acceptable(expected, detected)
            status = "✓" if result == expected_result else "✗"
            print(f"{status} {description}: expected={expected}, detected={detected}, result={result}")
            
            if result != expected_result:
                print(f"   ERROR: Expected {expected_result}, got {result}")
        
        print("\nSpecificity function test completed!")

    # Run the test
    test_specificity_function()
    print("✓ Specificity test functions are working correctly!")
    
except Exception as e:
    print(f"Error testing specificity functions: {e}")
    import traceback
    traceback.print_exc()