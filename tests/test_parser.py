"""
Tests for namelist parser.
"""

from pathlib import Path
import sys

# Add parent directory to path
sys.path.insert(0, str(Path(__file__).parent.parent))

from pydatcom.io.namelist_parser import NamelistParser
from pydatcom.io.state_manager import StateManager


def test_parser_basic():
    """Test basic parsing of example file."""
    parser = NamelistParser()
    
    # Parse ex1.inp
    fixture_path = Path(__file__).parent / 'fixtures' / 'ex1.inp'
    cases = parser.parse_file(fixture_path)
    
    print(f"Parsed {len(cases)} cases")
    
    # Check first case
    assert len(cases) > 0, "Should parse at least one case"
    
    first_case = cases[0]
    print(f"\nFirst case ID: {first_case['caseid']}")
    print(f"Namelists: {list(first_case['namelists'].keys())}")
    
    # Check FLTCON namelist
    if 'FLTCON' in first_case['namelists']:
        fltcon = first_case['namelists']['FLTCON']
        print(f"\nFLTCON parameters:")
        for key, value in fltcon.items():
            print(f"  {key}: {value}")
    
    # Check BODY namelist
    if 'BODY' in first_case['namelists']:
        body = first_case['namelists']['BODY']
        print(f"\nBODY parameters:")
        for key, value in body.items():
            if isinstance(value, list) and len(value) > 3:
                print(f"  {key}: [{value[0]}, {value[1]}, ... {value[-1]}] (len={len(value)})")
            else:
                print(f"  {key}: {value}")


def test_state_manager():
    """Test state manager."""
    state = StateManager()
    
    # Test get/set
    state.set('test_value', 42.0)
    assert state.get('test_value') == 42.0
    
    # Test component get
    flight_vars = state.get_component('flight')
    print(f"\nFlight variables: {list(flight_vars.keys())}")
    
    # Test update
    state.update({
        'flight_mach': [0.6, 0.9],
        'flight_nalpha': 5,
    })
    
    assert state.get('flight_mach') == [0.6, 0.9]
    print("\nState manager tests passed!")


def test_parser_to_state():
    """Test converting parsed case to state dict."""
    parser = NamelistParser()
    
    fixture_path = Path(__file__).parent / 'fixtures' / 'ex1.inp'
    cases = parser.parse_file(fixture_path)
    
    if cases:
        state_dict = parser.to_state_dict(cases[0])
        print(f"\nState dict keys: {len(state_dict)}")
        
        # Show some key values
        for key in ['flight_nmach', 'flight_mach', 'options_sref', 'body_nx']:
            if key in state_dict:
                print(f"  {key}: {state_dict[key]}")


if __name__ == '__main__':
    print("Testing Namelist Parser...")
    print("=" * 60)
    
    test_parser_basic()
    print("\n" + "=" * 60)
    
    test_state_manager()
    print("\n" + "=" * 60)
    
    test_parser_to_state()
    print("\n" + "=" * 60)
    print("\nAll tests passed!")

