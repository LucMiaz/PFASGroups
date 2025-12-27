"""
Command-line interface for PFASgroups.

Provides command-line tools for parsing PFAS structures and generating fingerprints.
"""

import argparse
import sys
import json
from pathlib import Path
from typing import Optional

from .core import parse_pfas, generate_pfas_fingerprint, get_smartsPaths, get_PFASGroups
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
  
Note: Use get_smartsPaths() and get_PFASGroups() in Python to extend defaults.
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
        help='Use component-based analysis'
    )
    parse_parser.add_argument(
        '--pretty',
        action='store_true',
        help='Pretty-print JSON output'
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
        '--pretty',
        action='store_true',
        help='Pretty-print JSON output'
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
        help='List available path types (use in Python to extend with get_smartsPaths)'
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
        kwargs['smartsPaths'] = get_smartsPaths(filename=args.fpaths_file)
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
    
    # Parse PFAS
    results = parse_pfas(smiles_list, bycomponent=args.bycomponent, **kwargs)
    
    # Convert results to JSON-serializable format
    output_data = []
    for i, (smiles, matches) in enumerate(zip(smiles_list, results)):
        result_entry = {
            'smiles': smiles,
            'groups': []
        }
        
        for group, n_matches, n_CFchain, chains in matches:
            result_entry['groups'].append({
                'name': group.name,
                'id': group.id,
                'n_matches': n_matches,
                'n_CFchain': n_CFchain,
                'chains_count': len(chains)
            })
        
        output_data.append(result_entry)
    
    # Output results
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
        kwargs['smartsPaths'] = get_smartsPaths(filename=args.fpaths_file)
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
    fps, group_info = generate_pfas_fingerprint(
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
    
    # Output results
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
            'smartsPath': group.smartsPath
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
        paths = get_smartsPaths(filename=args.fpaths_file)
    else:
        paths = get_smartsPaths()
    
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
        'note': 'Use get_smartsPaths() in Python to load and extend these paths'
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
            paths = get_smartsPaths(filename=args.fpaths_file)
            print(f"✓ fpaths.json loaded successfully from: {args.fpaths_file}")
            print(f"  Found {len(paths)} path types")
        
        # Validate groups if provided
        if args.groups_file:
            groups = get_PFASGroups(filename=args.groups_file)
            print(f"✓ PFAS_groups_smarts.json loaded successfully from: {args.groups_file}")
            print(f"  Found {len(groups)} PFAS groups")
        
        if not args.fpaths_file and not args.groups_file:
            # Validate defaults
            paths = get_smartsPaths()
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
