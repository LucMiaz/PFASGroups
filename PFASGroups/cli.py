"""
Command-line interface for PFASGroups.

Provides command-line tools for parsing PFAS structures and generating fingerprints.
"""

import argparse
import csv
import sys
import json


from .parser import parse_smiles
from .getter import get_componentSMARTSs, get_HalogenGroups
from .PFASEmbeddings import PFASEmbedding
from rdkit import Chem


def parse_args():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description='PFASGroups - Parse and analyze PFAS structures',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Parse SMILES from command line
  PFASGroups parse "C(C(F)(F)F)F" "FC(F)(F)C(F)(F)C(=O)O"

  # Parse with component metrics
  PFASGroups parse --bycomponent "FC(F)(F)C(F)(F)C(=O)O" --pretty

  # Parse SMILES from file
  PFASGroups parse --input smiles.txt --output results.json

  # Use custom configuration files
  PFASGroups parse --groups-file custom_groups.json "CCF"

  # Generate fingerprints
  PFASGroups fingerprint "C(C(F)(F)F)F" --output fp.json

  # Generate fingerprints with custom groups
  PFASGroups fingerprint --input smiles.txt --groups 28-52 --format dict

  # List available groups
  PFASGroups list-groups

  # List available path types
  PFASGroups list-paths

Note: Use get_componentSMARTSs() and get_PFASGroups() in Python to extend defaults.
        """
    )

    # Global options
    parser.add_argument(
        '--component_smarts-file',
        type=str,
        help='Path to custom component_smarts.json file (default: use package default)'
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
        '--no-component-metrics',
        action='store_true',
        help='Skip all component graph metrics (fastest)'
    )
    parse_parser.add_argument(
        '--limit-effective-graph-resistance',
        type=int,
        help='Only compute effective graph resistance for components smaller than this size (0 disables it)'
    )
    parse_parser.add_argument(
        '--halogens',
        nargs='+',
        help='Filter components by halogen element symbol(s), e.g. F or F Cl'
    )
    parse_parser.add_argument(
        '--form',
        nargs='+',
        choices=['alkyl', 'cyclic'],
        help='Filter components by form (alkyl, cyclic)'
    )
    parse_parser.add_argument(
        '--saturation',
        nargs='+',
        choices=['per', 'poly'],
        help='Filter components by saturation (per, poly)'
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
        '--halogens',
        nargs='+',
        help='Filter components by halogen element symbol(s), e.g. F or F Cl Br I'
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
        help='List available path types (use in Python to extend with get_componentSMARTSs)'
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
    validate_parser = subparsers.add_parser(  # pylint: disable=unused-variable
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
    if ',' in groups_str:
        return [int(x.strip()) for x in groups_str.split(',')]
    else:
        return [int(groups_str)]


def cmd_parse(args):
    """Execute parse command."""
    # Load custom configuration if provided
    kwargs = {}
    if args.component_smarts_file:
        kwargs['componentSmartss'] = get_componentSMARTSs(filename=args.component_smarts_file)
    if args.groups_file:
        kwargs['pfas_groups'] = get_HalogenGroups(filename=args.groups_file)
    kwargs['compute_component_metrics'] = not args.no_component_metrics
    kwargs['limit_effective_graph_resistance'] = args.limit_effective_graph_resistance
    if args.halogens:
        kwargs['halogens'] = args.halogens
    if args.form:
        kwargs['form'] = args.form
    if args.saturation:
        kwargs['saturation'] = args.saturation

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
    else:
        output_format = 'list'  # pylint: disable=unused-variable

    # Parse PFAS — always returns a PFASEmbeddingSet
    results = parse_smiles(smiles_list, bycomponent=args.bycomponent,
                           **kwargs)

    if args.format == 'csv':
        # Build CSV output from PFASEmbeddingSet
        import io
        writer_buf = io.StringIO()
        csv_writer = csv.writer(writer_buf)
        csv_writer.writerow(['smiles', 'group_id', 'group_name', 'match_count',
                             'component_idx', 'component_smarts', 'size',
                             'branching', 'mean_eccentricity', 'component_fraction',
                             'diameter', 'radius', 'effective_graph_resistance',
                             'n_spacer', 'ring_size'])
        for embedding in results:
            smiles_val = embedding['smiles']
            for match in embedding.get('matches', []):
                if match.get('type') != 'HalogenGroup':
                    continue
                components = match.get('components', [])
                if not components:
                    csv_writer.writerow([
                        smiles_val, match['id'], match['group_name'],
                        match['match_count'], '', '', '', '', '', '', '', '', '', '', ''
                    ])
                for idx, comp in enumerate(components):
                    csv_writer.writerow([
                        smiles_val,
                        match['id'],
                        match['group_name'],
                        match['match_count'],
                        idx,
                        comp.get('SMARTS', ''),
                        comp.get('size', ''),
                        comp.get('branching', ''),
                        comp.get('mean_eccentricity', ''),
                        comp.get('component_fraction', ''),
                        comp.get('diameter', ''),
                        comp.get('radius', ''),
                        comp.get('effective_graph_resistance', ''),
                        comp.get('n_spacer', ''),
                        comp.get('ring_size', ''),
                    ])
        result = writer_buf.getvalue()

        if args.output:
            with open(args.output, 'w', encoding='utf-8') as f:
                f.write(result)
            print(f"Results written to {args.output}")
        else:
            print(result, end='')
    else:
        # Convert to JSON format
        output_data = []
        for embedding in results:
            result_entry = {
                'smiles': embedding['smiles'],
                'groups': []
            }

            for match in embedding.get('matches', []):
                if match.get('type') != 'HalogenGroup':
                    continue
                components_out = []
                for comp in match.get('components', []):
                    components_out.append({
                        'smarts': comp.get('SMARTS'),
                        'size': comp.get('size'),
                        'branching': comp.get('branching'),
                        'mean_eccentricity': comp.get('mean_eccentricity'),
                        'component_fraction': comp.get('component_fraction'),
                        'diameter': comp.get('diameter'),
                        'radius': comp.get('radius'),
                        'effective_graph_resistance': comp.get('effective_graph_resistance'),
                        'n_spacer': comp.get('n_spacer'),
                        'ring_size': comp.get('ring_size'),
                    })
                result_entry['groups'].append({
                    'name': match['group_name'],
                    'id': match['id'],
                    'match_count': match['match_count'],
                    'num_components': match['num_components'],
                    'components_types': match.get('components_types', []),
                    'components': components_out,
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
    if args.component_smarts_file:
        kwargs['componentSmartss'] = get_componentSMARTSs(filename=args.component_smarts_file)
    if args.groups_file:
        kwargs['pfas_groups'] = get_HalogenGroups(filename=args.groups_file)

    # Get SMILES from command line or file
    if args.input:
        smiles_list = read_smiles_file(args.input)
    elif args.smiles:
        smiles_list = args.smiles
    else:
        print("Error: Provide SMILES as arguments or use --input", file=sys.stderr)
        sys.exit(1)

    # Determine halogens
    halogens = getattr(args, 'halogens', None) or 'F'

    # Parse SMILES → PFASEmbeddingSet
    embs = parse_smiles(smiles_list, halogens=halogens, **kwargs)

    # Build fingerprint array and column names
    array_kwargs = {}
    if args.groups:
        array_kwargs['selected_group_ids'] = parse_group_selection(args.groups)
    arr = embs.to_array(**array_kwargs)
    col_names = embs.column_names(**{k: v for k, v in array_kwargs.items()
                                    if k in ('selected_group_ids',)})

    # Format fingerprints per molecule
    fingerprints = []
    for i, smiles_val in enumerate(smiles_list):
        row = arr[i].tolist()
        if args.format == 'dict':
            fingerprints.append({
                'smiles': smiles_val,
                'fingerprint': {col_names[j]: row[j] for j in range(len(col_names))}
            })
        else:  # vector
            fingerprints.append({
                'smiles': smiles_val,
                'fingerprint': row
            })

    # Build output
    output_data = {
        'smiles': smiles_list,
        'column_names': col_names,
        'fingerprints': fingerprints
    }

    # Output based on requested format
    if args.output_format == 'csv':
        csv_rows = []
        for entry in fingerprints:
            row = {'smiles': entry['smiles']}
            fp = entry['fingerprint']
            if isinstance(fp, dict):
                row.update(fp)
            else:
                for j, val in enumerate(fp):
                    row[col_names[j] if j < len(col_names) else f'col_{j}'] = val
            csv_rows.append(row)

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

    # JSON output (default)
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
        groups = get_HalogenGroups(filename=args.groups_file)
    else:
        groups = get_HalogenGroups()

    # Create output
    groups_list = []
    for i, group in enumerate(groups):
        groups_list.append({
            'index': i,
            'id': group['id'],
            'name': group['name'],
            'smarts': list(group.get('smarts', {}).keys()),
            'componentSmarts': group.get('componentSmarts')
        })

    output_data = {
        'total_groups': len(groups),
        'groups': groups_list,
        'note': 'Use get_HalogenGroups() in Python to load and extend these groups'
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
    if args.component_smarts_file:
        paths = get_componentSMARTSs(filename=args.component_smarts_file)
    else:
        paths = get_componentSMARTSs()

    # Create output - convert RDKit mols to SMARTS strings for display
    paths_list = []
    for name, path_info in paths.items():
        component_mol = path_info.get('component')
        paths_list.append({
            'name': name,
            'smarts': Chem.MolToSmarts(component_mol) if component_mol else None,
            'halogen': path_info.get('halogen'),
            'form': path_info.get('form'),
            'saturation': path_info.get('saturation')
        })

    output_data = {
        'total_paths': len(paths),
        'paths': paths_list,
        'note': 'Use get_componentSMARTSs() in Python to load and extend these paths'
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

        # Validate component_smarts if provided
        if args.component_smarts_file:
            paths = get_componentSMARTSs(filename=args.component_smarts_file)
            print(f"✓ component_smarts.json loaded successfully from: {args.component_smarts_file}")
            print(f"  Found {len(paths)} path types")

        # Validate groups if provided
        if args.groups_file:
            groups = get_HalogenGroups(filename=args.groups_file)
            print(f"✓ PFAS_groups_smarts.json loaded successfully from: {args.groups_file}")
            print(f"  Found {len(groups)} PFAS groups")

        if not args.component_smarts_file and not args.groups_file:
            # Validate defaults
            paths = get_componentSMARTSs()
            groups = get_HalogenGroups()
            print("✓ Default configuration loaded successfully")
            print(f"  Path types: {len(paths)}")
            print(f"  PFAS groups: {len(groups)}")

        print("\nConfiguration is valid!")

    except FileNotFoundError as e:
        print("✗ Error: " + str(e), file=sys.stderr)
        sys.exit(1)


def main(default_halogens=None):
    """Main CLI entry point."""
    args = parse_args()

    # Apply default halogens when the entry point provides one and the user
    # did not explicitly pass --halogens on the command line.
    if default_halogens and args.command in ('parse', 'fingerprint'):
        if not getattr(args, 'halogens', None):
            args.halogens = default_halogens

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


def main_halogen():
    """Entry point for the ``halogengroups`` CLI.

    Identical to ``pfasgroups`` but defaults to all four halogens
    (F, Cl, Br, I) when ``--halogens`` is not specified on the command line.
    """
    main(default_halogens=['F', 'Cl', 'Br', 'I'])


if __name__ == '__main__':
    main()
