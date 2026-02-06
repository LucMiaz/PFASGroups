"""
Command-line interface for PFASgroups.

Provides command-line tools for parsing PFAS structures and generating fingerprints.
"""

import argparse
import sys
import json
from pathlib import Path
from typing import Optional

from .parser import parse_smiles, generate_fingerprint, get_componentSmartss, get_PFASGroups
from rdkit import Chem


def parse_args():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description='PFASgroups - Parse and analyze PFAS structures',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Parse SMILES from command line
  pfasgroups parse "C(C(F)(F)F)F" "FC(F)(F)C(F)(F)C(=O)O"
  
  # Parse with component metrics
  pfasgroups parse --bycomponent "FC(F)(F)C(F)(F)C(=O)O" --pretty
  
  # Parse SMILES from file
  pfasgroups parse --input smiles.txt --output results.json
  
  # Use custom configuration files
  pfasgroups parse --groups-file custom_groups.json "CCF"
  
  # Generate fingerprints
  pfasgroups fingerprint "C(C(F)(F)F)F" --output fp.json
  
  # Generate fingerprints with custom groups
  pfasgroups fingerprint --input smiles.txt --groups 28-52 --format dict
  
  # List available groups
  pfasgroups list-groups
  
  # List available path types
  pfasgroups list-paths
  
Note: Use get_componentSmartss() and get_PFASGroups() in Python to extend defaults.
        """
    )
    
    # Global options
    parser.add_argument(
        '--fpaths-file',
        type=str,
        help='Path to custom fpaths.json file (default: use package default)'
    )
    parser.add_argument(
        '--groups-file',
        type=str,
        help='Path to custom PFAS_groups_smarts.json file (default: use package default)'
    )
    
    subparsers = parser.add_subparsers(dest='command', help='Available commands')
    
    # Parse command
    parse_parser = subparsers.add_parser(
        'parse',
        help='Parse SMILES strings and identify PFAS groups'
    )
    parse_parser.add_argument(
        'smiles',
        nargs='*',
        help='SMILES strings to parse (use --input for file input)'
    )
    parse_parser.add_argument(
        '-i', '--input',
        type=str,
        help='Input file containing SMILES strings (one per line)'
    )
    parse_parser.add_argument(
        '-o', '--output',
        type=str,
        help='Output file for results (JSON format, default: stdout)'
    )
    parse_parser.add_argument(
        '--bycomponent',
        action='store_true',
        help='Use component-based analysis (provides comprehensive metrics including component_fraction, branching, eccentricity)'
    )
    parse_parser.add_argument(
        '--format',
        choices=['json', 'csv'],
        default='json',
        help='Output format (default: json)'
    )
    parse_parser.add_argument(
        '--pretty',
        action='store_true',
        help='Pretty-print JSON output (only for JSON format)'
    )
    
    # Fingerprint command
    fp_parser = subparsers.add_parser(
        'fingerprint',
        help='Generate PFAS group fingerprints'
    )
    fp_parser.add_argument(
        'smiles',
        nargs='*',
        help='SMILES strings to fingerprint (use --input for file input)'
    )
    fp_parser.add_argument(
        '-i', '--input',
        type=str,
        help='Input file containing SMILES strings (one per line)'
    )
    fp_parser.add_argument(
        '-o', '--output',
        type=str,
        help='Output file for fingerprints (JSON format, default: stdout)'
    )
    fp_parser.add_argument(
        '-g', '--groups',
        type=str,
        help='Selected groups as range (e.g., "28-52") or comma-separated indices (e.g., "28,29,30")'
    )
    fp_parser.add_argument(
        '-f', '--format',
        choices=['vector', 'dict', 'sparse', 'detailed', 'int'],
        default='vector',
        help='Fingerprint representation format (default: vector)'
    )
    fp_parser.add_argument(
        '--count-mode',
        choices=['binary', 'count', 'max_chain'],
        default='binary',
        help='How to count matches (default: binary)'
    )
    fp_parser.add_argument(
        '--output-format',
        choices=['json', 'csv'],
        default='json',
        help='Output file format (default: json)'
    )
    fp_parser.add_argument(
        '--pretty',
        action='store_true',
        help='Pretty-print JSON output (only for JSON format)'
    )
    
    # List groups command
    list_parser = subparsers.add_parser(
        'list-groups',
        help='List available PFAS groups (use in Python to extend with get_PFASGroups)'
    )
    list_parser.add_argument(
        '-o', '--output',
        type=str,
        help='Output file (default: stdout)'
    )
    list_parser.add_argument(
        '--pretty',
        action='store_true',
        default=True,
        help='Pretty-print JSON output'
    )
    
    # List paths command
    list_paths_parser = subparsers.add_parser(
        'list-paths',
        help='List available path types (use in Python to extend with get_componentSmartss)'
    )
    list_paths_parser.add_argument(
        '-o', '--output',
        type=str,
        help='Output file (default: stdout)'
    )
    list_paths_parser.add_argument(
        '--pretty',
        action='store_true',
        default=True,
        help='Pretty-print JSON output'
    )
    
    # Validate config command
    validate_parser = subparsers.add_parser(
        'validate-config',
        help='Validate custom configuration files'
    )
    
    return parser.parse_args()


def read_smiles_file(filepath: str) -> list:
    """Read SMILES strings from file."""
    with open(filepath, 'r') as f:
        return [line.strip() for line in f if line.strip()]


def parse_group_selection(groups_str: str) -> list:
    """
    Parse group selection string.
    
    Examples:
        "28-52" -> range(28, 53)
        "28,29,30" -> [28, 29, 30]
    """
    if '-' in groups_str:
        start, end = groups_str.split('-')
        return list(range(int(start), int(end) + 1))
    elif ',' in groups_str:
        return [int(x.strip()) for x in groups_str.split(',')]
    else:
        return [int(groups_str)]


def cmd_parse(args):
    """Execute parse command."""
    # Load custom configuration if provided
    kwargs = {}
    if args.fpaths_file:
        kwargs['componentSmartss'] = get_componentSmartss(filename=args.fpaths_file)
    if args.groups_file:
        kwargs['pfas_groups'] = get_PFASGroups(filename=args.groups_file)
    
    # Get SMILES from command line or file
    if args.input:
        smiles_list = read_smiles_file(args.input)
    elif args.smiles:
        smiles_list = args.smiles
    else:
        print("Error: Provide SMILES as arguments or use --input", file=sys.stderr)
        sys.exit(1)
    
    # Determine output format
    if args.format == 'csv':
        output_format = 'csv'
    elif args.format == 'json':
        output_format = 'list'
    else:
        output_format = 'list'
    
    # Parse PFAS with specified output format
    if args.format == 'csv':
        # Use CSV output format from parse_smiles
        result = parse_smiles(smiles_list, bycomponent=args.bycomponent, 
                            output_format='csv', **kwargs)
        
        if args.output:
            with open(args.output, 'w', encoding='utf-8') as f:
                f.write(result)
            print(f"Results written to {args.output}")
        else:
            print(result, end='')
    else:
        # Use list output format and convert to JSON
        results = parse_smiles(smiles_list, bycomponent=args.bycomponent, 
                             output_format='list', **kwargs)
        
        # Convert to JSON format
        output_data = []
        for smiles, matches in zip(smiles_list, results):
            result_entry = {
                'smiles': smiles,
                'groups': []
            }
            
            for group, match_count, chain_lengths, matched_chains in matches:
                result_entry['groups'].append({
                    'name': group.name,
                    'id': group.id,
                    'match_count': match_count,
                    'chain_lengths': chain_lengths,
                    'num_chains': len(matched_chains)
                })
            
            output_data.append(result_entry)
        
        # Output JSON
        indent = 2 if args.pretty else None
        output_json = json.dumps(output_data, indent=indent)
        
        if args.output:
            with open(args.output, 'w') as f:
                f.write(output_json)
            print(f"Results written to {args.output}")
        else:
            print(output_json)


def cmd_fingerprint(args):
    """Execute fingerprint command."""
    # Load custom configuration if provided
    kwargs = {}
    if args.fpaths_file:
        kwargs['componentSmartss'] = get_componentSmartss(filename=args.fpaths_file)
    if args.groups_file:
        kwargs['pfas_groups'] = get_PFASGroups(filename=args.groups_file)
    
    # Get SMILES from command line or file
    if args.input:
        smiles_list = read_smiles_file(args.input)
    elif args.smiles:
        smiles_list = args.smiles
    else:
        print("Error: Provide SMILES as arguments or use --input", file=sys.stderr)
        sys.exit(1)
    
    # Parse group selection
    selected_groups = None
    if args.groups:
        selected_groups = parse_group_selection(args.groups)
    
    # Generate fingerprints
    fps, group_info = generate_fingerprint(
        smiles_list,
        selected_groups=selected_groups,
        representation=args.format,
        count_mode=args.count_mode,
        **kwargs
    )
    
    # Prepare output
    output_data = {
        'fingerprints': [],
        'group_info': group_info
    }
    
    # Handle different representations
    if args.format == 'vector':
        # Convert numpy arrays to lists
        import numpy as np
        if isinstance(fps, list):
            output_data['fingerprints'] = [fp.tolist() if isinstance(fp, np.ndarray) else fp for fp in fps]
        else:
            output_data['fingerprints'] = fps.tolist() if isinstance(fps, np.ndarray) else fps
    elif args.format in ['int']:
        output_data['fingerprints'] = fps if isinstance(fps, list) else [fps]
    else:  # dict, sparse, detailed
        output_data['fingerprints'] = fps if isinstance(fps, list) else [fps]
    
    # Add SMILES to output for reference
    output_data['smiles'] = smiles_list
    
    # Output based on requested format
    if args.output_format == 'csv':
        # CSV output for fingerprints
        csv_rows = []
        if args.format == 'vector':
            # For vector representation, create columns for each group
            import numpy as np
            fps_list = fps if isinstance(fps, list) else [fps]
            for i, (smiles, fp) in enumerate(zip(smiles_list, fps_list)):
                row = {'smiles': smiles}
                fp_array = fp.tolist() if isinstance(fp, np.ndarray) else fp
                for j, val in enumerate(fp_array):
                    group_idx = group_info['selected_indices'][j] if 'selected_indices' in group_info else j
                    group_name = group_info['group_names'][group_idx] if group_idx < len(group_info['group_names']) else f'group_{group_idx}'
                    row[f'{group_name}'] = val
                csv_rows.append(row)
        elif args.format == 'dict':
            # For dict representation
            fps_list = fps if isinstance(fps, list) else [fps]
            for smiles, fp_dict in zip(smiles_list, fps_list):
                row = {'smiles': smiles}
                row.update(fp_dict)
                csv_rows.append(row)
        else:
            print(f"Warning: CSV output not optimized for '{args.format}' representation, using JSON", file=sys.stderr)
            args.output_format = 'json'
        
        if args.output_format == 'csv':
            fieldnames = list(csv_rows[0].keys()) if csv_rows else ['smiles']
            if args.output:
                with open(args.output, 'w', newline='', encoding='utf-8') as f:
                    writer = csv.DictWriter(f, fieldnames=fieldnames)
                    writer.writeheader()
                    writer.writerows(csv_rows)
                print(f"Results written to {args.output}")
            else:
                writer = csv.DictWriter(sys.stdout, fieldnames=fieldnames)
                writer.writeheader()
                writer.writerows(csv_rows)
            return
    
    # JSON output (default or fallback)
    indent = 2 if args.pretty else None
    output_json = json.dumps(output_data, indent=indent)
    
    if args.output:
        with open(args.output, 'w') as f:
            f.write(output_json)
        print(f"Fingerprints written to {args.output}")
    else:
        print(output_json)


def cmd_list_groups(args):
    """Execute list-groups command."""
    # Load groups (custom or default)
    if args.groups_file:
        groups = get_PFASGroups(filename=args.groups_file)
    else:
        groups = get_PFASGroups()
    
    # Create output
    groups_list = []
    for i, group in enumerate(groups):
        groups_list.append({
            'index': i,
            'id': group.id,
            'name': group.name,
            'smarts1': group.smarts1,
            'smarts2': group.smarts2,
            'componentSmarts': group.componentSmarts
        })
    
    output_data = {
        'total_groups': len(groups),
        'groups': groups_list,
        'note': 'Use get_PFASGroups() in Python to load and extend these groups'
    }
    
    # Output results
    indent = 2 if args.pretty else None
    output_json = json.dumps(output_data, indent=indent)
    
    if args.output:
        with open(args.output, 'w') as f:
            f.write(output_json)
        print(f"Groups list written to {args.output}")
    else:
        print(output_json)


def cmd_list_paths(args):
    """Execute list-paths command."""
    # Load paths (custom or default)
    if args.fpaths_file:
        paths = get_componentSmartss(filename=args.fpaths_file)
    else:
        paths = get_componentSmartss()
    
    # Create output - convert RDKit mols to SMARTS strings for display
    paths_list = []
    for name, (chain_mol, end_mol) in paths.items():
        paths_list.append({
            'name': name,
            'chain_smarts': Chem.MolToSmarts(chain_mol),
            'end_smarts': Chem.MolToSmarts(end_mol)
        })
    
    output_data = {
        'total_paths': len(paths),
        'paths': paths_list,
        'note': 'Use get_componentSmartss() in Python to load and extend these paths'
    }
    
    # Output results
    indent = 2 if args.pretty else None
    output_json = json.dumps(output_data, indent=indent)
    
    if args.output:
        with open(args.output, 'w') as f:
            f.write(output_json)
        print(f"Paths list written to {args.output}")
    else:
        print(output_json)


def cmd_validate_config(args):
    """Execute validate-config command."""
    try:
        print("Validating configuration files...")
        
        # Validate fpaths if provided
        if args.fpaths_file:
            paths = get_componentSmartss(filename=args.fpaths_file)
            print(f"✓ fpaths.json loaded successfully from: {args.fpaths_file}")
            print(f"  Found {len(paths)} path types")
        
        # Validate groups if provided
        if args.groups_file:
            groups = get_PFASGroups(filename=args.groups_file)
            print(f"✓ PFAS_groups_smarts.json loaded successfully from: {args.groups_file}")
            print(f"  Found {len(groups)} PFAS groups")
        
        if not args.fpaths_file and not args.groups_file:
            # Validate defaults
            paths = get_componentSmartss()
            groups = get_PFASGroups()
            print(f"✓ Default configuration loaded successfully")
            print(f"  Path types: {len(paths)}")
            print(f"  PFAS groups: {len(groups)}")
        
        print("\nConfiguration is valid!")
        
    except FileNotFoundError as e:
        print(f"✗ Error: {e}", file=sys.stderr)
        sys.exit(1)
    except json.JSONDecodeError as e:
        print(f"✗ JSON parsing error: {e}", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"✗ Validation error: {e}", file=sys.stderr)
        sys.exit(1)


def main():
    """Main CLI entry point."""
    args = parse_args()
    
    if args.command == 'parse':
        cmd_parse(args)
    elif args.command == 'fingerprint':
        cmd_fingerprint(args)
    elif args.command == 'list-groups':
        cmd_list_groups(args)
    elif args.command == 'list-paths':
        cmd_list_paths(args)
    elif args.command == 'validate-config':
        cmd_validate_config(args)
    else:
        print("Error: No command specified. Use --help for usage information.", file=sys.stderr)
        sys.exit(1)


if __name__ == '__main__':
    main()
