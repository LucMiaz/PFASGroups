#!/usr/bin/env python3
"""
Test the improved group formatting in the automated test runner
"""

import sys
import os
import json

# Add path for imports
sys.path.insert(0, os.getcwd())

def test_group_formatting():
    """Test the new group formatting function"""
    
    print("🔍 Testing Enhanced Group Formatting")
    print("=" * 50)
    
    # Load group names
    try:
        data_folder = os.path.join('PFASgroups', 'data')
        with open(os.path.join(data_folder, 'PFAS_groups_smarts.json'), 'r') as f:
            group_data = json.load(f)
        group_names = {group['id']: group['name'] for group in group_data}
        print(f"✓ Loaded {len(group_names)} group names")
    except Exception as e:
        print(f"❌ Error loading group names: {e}")
        return
    
    # Test the formatting function (copy from automated_test_runner.py)
    def format_groups(group_list, with_names=True):
        if not group_list:
            return '<span class="group-names">None</span>'
        if with_names and group_names:
            formatted_groups = []
            for g in sorted(group_list):
                group_name = group_names.get(g, 'Unknown')
                # Truncate very long group names for better display
                if len(group_name) > 40:
                    group_name = group_name[:37] + "..."
                formatted_groups.append(f'<span class="group-id">G{g}</span> <span class="group-desc">({group_name})</span>')
            return '<div class="group-names">' + '<br>'.join(formatted_groups) + '</div>'
        else:
            return '<div class="group-names">' + ", ".join([f'<span class="group-id">G{g}</span>' for g in sorted(group_list)]) + '</div>'
    
    # Test cases
    test_cases = [
        [],  # Empty list
        [8],  # Single group (the one we just fixed)
        [7, 8, 10, 31, 36],  # Multiple groups from our previous example
        [1, 2, 3],  # Some other groups
        [15, 25],  # Test with longer names
    ]
    
    print("\n🔍 Testing Different Group Combinations:")
    print("-" * 50)
    
    for i, groups in enumerate(test_cases, 1):
        print(f"\nTest Case {i}: {groups}")
        formatted = format_groups(groups)
        # Remove HTML tags for clean console display
        clean_text = formatted.replace('<div class="group-names">', '').replace('</div>', '')
        clean_text = clean_text.replace('<span class="group-id">', '').replace('</span>', '')
        clean_text = clean_text.replace('<span class="group-desc">', '').replace('<br>', ' | ')
        print(f"Formatted: {clean_text}")
        
        # Show the actual HTML that would be in the report
        print(f"HTML: {formatted}")
    
    print(f"\n🔍 Sample Group Names:")
    print("-" * 30)
    for group_id in [6, 7, 8, 10, 15, 31, 36]:
        name = group_names.get(group_id, "Unknown")
        print(f"Group {group_id:2d}: {name}")

if __name__ == "__main__":
    test_group_formatting()