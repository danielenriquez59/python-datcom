"""
Command-line interface for PyDATCOM.

Usage:
    pydatcom parse input.inp          # Parse and display input file
    pydatcom run input.inp            # Run analysis (future)
    pydatcom convert input.inp -f yaml # Convert to YAML (future)
"""

import argparse
import sys
from pathlib import Path
from pydatcom.io import NamelistParser, StateManager
import json


def parse_command(args):
    """Parse DATCOM input file and display contents."""
    parser = NamelistParser()
    input_path = Path(args.input)
    
    if not input_path.exists():
        print(f"Error: Input file '{input_path}' not found")
        return 1
    
    print(f"Parsing: {input_path}")
    print("=" * 60)
    
    cases = parser.parse_file(input_path)
    print(f"\nFound {len(cases)} case(s)")
    
    for i, case in enumerate(cases, 1):
        print(f"\n--- Case {i} ---")
        if case.get('caseid'):
            print(f"ID: {case['caseid']}")
        
        print(f"Namelists: {', '.join(case['namelists'].keys())}")
        
        if args.verbose:
            print("\nNamelist Details:")
            for name, params in case['namelists'].items():
                print(f"\n  ${name}")
                for key, value in params.items():
                    if isinstance(value, list) and len(value) > 5:
                        print(f"    {key}: [{value[0]}, {value[1]}, ... {value[-1]}] (len={len(value)})")
                    else:
                        print(f"    {key}: {value}")
        
        if case.get('commands'):
            print(f"Commands: {', '.join(case['commands'])}")
    
    # Convert first case to state dict
    if cases and args.state:
        print("\n" + "=" * 60)
        print("State Dictionary (Case 1):")
        print("=" * 60)
        state_dict = parser.to_state_dict(cases[0])
        
        # Group by prefix
        prefixes = {}
        for key, value in state_dict.items():
            prefix = key.split('_')[0] if '_' in key else 'other'
            if prefix not in prefixes:
                prefixes[prefix] = {}
            prefixes[prefix][key] = value
        
        for prefix in sorted(prefixes.keys()):
            print(f"\n{prefix.upper()}:")
            for key, value in sorted(prefixes[prefix].items()):
                if isinstance(value, list) and len(value) > 5:
                    print(f"  {key}: [{value[0]}, ... {value[-1]}] (len={len(value)})")
                else:
                    print(f"  {key}: {value}")
    
    return 0


def run_command(args):
    """Run DATCOM analysis (placeholder)."""
    print("Error: Run command not yet implemented")
    print("Current status: Phase 1 (I/O and utilities) complete")
    print("Aerodynamic calculations coming in Phase 4")
    return 1


def convert_command(args):
    """Convert input file to different format."""
    parser = NamelistParser()
    input_path = Path(args.input)
    
    if not input_path.exists():
        print(f"Error: Input file '{input_path}' not found")
        return 1
    
    cases = parser.parse_file(input_path)
    
    if args.format == 'yaml':
        import yaml
        output = {'cases': []}
        for case in cases:
            state_dict = parser.to_state_dict(case)
            output['cases'].append(state_dict)
        
        output_path = Path(args.output) if args.output else input_path.with_suffix('.yaml')
        with open(output_path, 'w') as f:
            yaml.dump(output, f, default_flow_style=False, sort_keys=True)
        print(f"Converted to YAML: {output_path}")
    
    elif args.format == 'json':
        output = {'cases': []}
        for case in cases:
            state_dict = parser.to_state_dict(case)
            output['cases'].append(state_dict)
        
        output_path = Path(args.output) if args.output else input_path.with_suffix('.json')
        with open(output_path, 'w') as f:
            json.dump(output, f, indent=2)
        print(f"Converted to JSON: {output_path}")
    
    return 0


def main():
    """Main CLI entry point."""
    parser = argparse.ArgumentParser(
        description='PyDATCOM - Python Digital DATCOM',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  pydatcom parse input.inp -v          # Parse with verbose output
  pydatcom parse input.inp --state     # Show state dictionary
  pydatcom convert input.inp -f yaml   # Convert to YAML
  pydatcom convert input.inp -f json -o output.json
        """
    )
    
    subparsers = parser.add_subparsers(dest='command', help='Command to execute')
    
    # Parse command
    parse_parser = subparsers.add_parser('parse', help='Parse DATCOM input file')
    parse_parser.add_argument('input', help='Input .inp file')
    parse_parser.add_argument('-v', '--verbose', action='store_true', 
                             help='Verbose output')
    parse_parser.add_argument('-s', '--state', action='store_true',
                             help='Show state dictionary')
    
    # Run command (placeholder)
    run_parser = subparsers.add_parser('run', help='Run DATCOM analysis')
    run_parser.add_argument('input', help='Input .inp file')
    run_parser.add_argument('-o', '--output', help='Output file')
    
    # Convert command
    conv_parser = subparsers.add_parser('convert', help='Convert input format')
    conv_parser.add_argument('input', help='Input .inp file')
    conv_parser.add_argument('-f', '--format', choices=['yaml', 'json'],
                            required=True, help='Output format')
    conv_parser.add_argument('-o', '--output', help='Output file')
    
    args = parser.parse_args()
    
    if not args.command:
        parser.print_help()
        return 0
    
    if args.command == 'parse':
        return parse_command(args)
    elif args.command == 'run':
        return run_command(args)
    elif args.command == 'convert':
        return convert_command(args)
    
    return 0


if __name__ == '__main__':
    sys.exit(main())

