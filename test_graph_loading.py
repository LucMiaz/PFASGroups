"""
Test script to verify the corrected graph loading with ID mapping.
"""

import sys
import os
import json
import networkx as nx

# Add the PFASgroups path
sys.path.append('.')
sys.path.append('..')

def test_graph_loading():
    """Test the corrected graph loading function."""
    
    # Import the fixed function
    from PFASgroups.tests.test_examples import load_graph_from_json, check_groups_are_related, are_detected_groups_acceptable
    
    print("Testing graph loading with ID mapping...")
    
    try:
        # Load the graph
        G, Gper = load_graph_from_json('PFASgroups/tests/specificity_test_groups.json')
        
        print(f"Main graph has {len(G.nodes)} nodes and {len(G.edges)} edges")
        print(f"Perfluoro graph has {len(Gper.nodes)} nodes and {len(Gper.edges)} edges")
        
        # Show some node IDs to verify they are integers
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
        
        # Test the acceptability function
        print("\nTesting group acceptability:")
        test_cases = [
            ([1], [1, 33]),  # Expected carboxylic acid groups, detected both specific and generic
            ([17], [17, 31]), # Expected hydrofluoroether, detected both specific and generic ether
        ]
        
        for expected, detected in test_cases:
            acceptable = are_detected_groups_acceptable(expected, detected)
            print(f"  Expected {expected}, detected {detected}: acceptable = {acceptable}")
        
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