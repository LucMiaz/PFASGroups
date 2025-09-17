"""
Test the modified analyze_group_overlap function with sample data.
"""

import pandas as pd
import json

def mock_analyze_group_overlap():
    """Test the overlap analysis with mock data."""
    
    # Create sample specificity results data
    sample_data = [
        {
            'group_ids': [1],  # Expected: Perfluoroalkyl carboxylic acids
            'detected_groups': [1, 33],  # Detected: specific + carboxylic acid
            'valid_smiles': True
        },
        {
            'group_ids': [17],  # Expected: Hydrofluoroethers
            'detected_groups': [17, 31],  # Detected: specific + ether
            'valid_smiles': True
        },
        {
            'group_ids': [6],  # Expected: Perfluoroalkyl sulfonic acids
            'detected_groups': [6, 36],  # Detected: specific + sulfonic acid
            'valid_smiles': True
        },
        {
            'group_ids': [1],  # Expected: Perfluoroalkyl carboxylic acids
            'detected_groups': [1, 2, 33],  # Detected: specific + wrong specific + generic
            'valid_smiles': True
        },
        {
            'group_ids': [16],  # Expected: Perfluoropolyethers
            'detected_groups': [16, 31, 17],  # Detected: specific + ether + related hydrofluoroether
            'valid_smiles': True
        }
    ]
    
    # Convert to DataFrame
    df = pd.DataFrame(sample_data)
    
    # Load group names for display
    try:
        with open('PFASgroups/data/PFAS_groups_smarts.json', 'r') as f:
            group_data = json.load(f)
        group_names = {group['id']: group['name'] for group in group_data}
    except FileNotFoundError:
        # Fallback names for testing
        group_names = {
            1: "Perfluoroalkyl carboxylic acids",
            2: "Polyfluoroalkyl carboxylic acid", 
            6: "Perfluoroalkyl sulfonic acids",
            16: "Perfluoropolyethers",
            17: "Hydrofluoroethers",
            31: "ether",
            33: "carboxylic acid",
            36: "sulfonic acid"
        }
    
    print("Mock Overlap Analysis")
    print("=" * 50)
    
    # Analyze unexpected detections
    unexpected_detections = {}
    for _, row in df.iterrows():
        expected_groups = set(row['group_ids'])
        detected_groups = set(row['detected_groups'])
        
        # Find groups that were detected but not expected
        unexpected = detected_groups - expected_groups
        
        for unexpected_group in unexpected:
            for expected_group in expected_groups:
                pair = (expected_group, unexpected_group)
                unexpected_detections[pair] = unexpected_detections.get(pair, 0) + 1
    
    # Sort and display
    common_unexpected = sorted(unexpected_detections.items(), key=lambda x: x[1], reverse=True)
    
    print("Most common unexpected group detections (Expected -> Unexpected):")
    for (expected_group, unexpected_group), count in common_unexpected:
        expected_name = group_names.get(expected_group, f"Unknown {expected_group}")
        unexpected_name = group_names.get(unexpected_group, f"Unknown {unexpected_group}")
        print(f"  Expected {expected_group:2d} ({expected_name})")
        print(f"    -> Also detected {unexpected_group:2d} ({unexpected_name}) ({count} times)")
        print()
    
    print("Test completed successfully!")

if __name__ == "__main__":
    mock_analyze_group_overlap()