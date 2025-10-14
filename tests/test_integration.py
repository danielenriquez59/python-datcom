"""
Integration tests combining parser, state manager, and geometry calculations.
"""

from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).parent.parent))

from pydatcom.io import NamelistParser, StateManager
from pydatcom.geometry import (
    calculate_body_geometry,
    calculate_wing_geometry,
    generate_naca_airfoil,
    BodyGeometry
)


def test_parse_and_calculate_body():
    """Test parsing EX1 and calculating body geometry."""
    print("\n" + "="*60)
    print("Integration Test: Parse + Body Geometry")
    print("="*60)
    
    # Parse example file
    parser = NamelistParser()
    fixture_path = Path(__file__).parent / 'fixtures' / 'ex1.inp'
    cases = parser.parse_file(fixture_path)
    
    # Get first case
    state_dict = parser.to_state_dict(cases[0])
    
    print(f"\nParsed case 1:")
    print(f"  Body stations (NX): {state_dict.get('body_nx')}")
    print(f"  Body length: {state_dict.get('body_x', [])[-1] if state_dict.get('body_x') else 'N/A'}")
    
    # Calculate body geometry
    body_props = calculate_body_geometry(state_dict)
    
    print(f"\nCalculated Body Properties:")
    print(f"  Length: {body_props['length']:.3f} ft")
    print(f"  Max area: {body_props['max_area']:.3f} ft²")
    print(f"  Volume: {body_props['volume']:.3f} ft³")
    print(f"  Centroid: {body_props['centroid']:.3f} ft from nose")
    print(f"  Fineness ratio: {body_props['fineness_ratio']:.3f}")
    
    assert body_props['length'] == 6.26
    assert body_props['max_area'] == 1.032
    
    print("  [PASS] Parse + Body integration successful")


def test_full_pipeline():
    """Test complete pipeline: Parse -> State -> Geometry."""
    print("\n" + "="*60)
    print("Integration Test: Complete Pipeline")
    print("="*60)
    
    # Step 1: Parse
    parser = NamelistParser()
    fixture_path = Path(__file__).parent / 'fixtures' / 'ex1.inp'
    cases = parser.parse_file(fixture_path)
    
    print(f"\nStep 1: Parsed {len(cases)} case(s)")
    
    # Step 2: Load into state manager
    state_mgr = StateManager()
    case_state = parser.to_state_dict(cases[0])
    state_mgr.update(case_state)
    
    print(f"Step 2: Loaded into state manager")
    print(f"  Flight Mach: {state_mgr.get('flight_mach')}")
    print(f"  Wing area (SREF): {state_mgr.get('options_sref')} ft²")
    
    # Step 3: Calculate body geometry
    body = BodyGeometry(state_mgr.get_all())
    body_props = body.calculate_properties()
    
    print(f"\nStep 3: Body Geometry Calculated")
    print(f"  Length: {body_props['length']:.2f} ft")
    print(f"  Volume: {body_props['volume']:.2f} ft³")
    
    # Step 4: Export computed values back to state
    state_mgr.update(body.to_state_dict())
    
    print(f"\nStep 4: Updated state with computed values")
    print(f"  body_volume: {state_mgr.get('body_volume'):.2f}")
    print(f"  body_fineness_ratio: {state_mgr.get('body_fineness_ratio'):.2f}")
    
    # Step 5: Export to YAML
    output_path = Path(__file__).parent / 'integration_test_output.yaml'
    state_mgr.export_to_yaml(output_path)
    
    print(f"\nStep 5: Exported to {output_path.name}")
    
    assert output_path.exists()
    
    print("\n  [PASS] Complete pipeline successful")
    
    return state_mgr


def test_multiple_cases():
    """Test handling multiple cases from input file."""
    print("\n" + "="*60)
    print("Integration Test: Multiple Cases")
    print("="*60)
    
    parser = NamelistParser()
    fixture_path = Path(__file__).parent / 'fixtures' / 'ex1.inp'
    cases = parser.parse_file(fixture_path)
    
    print(f"\nProcessing {len(cases)} cases:\n")
    
    for i, case in enumerate(cases, 1):
        case_state = parser.to_state_dict(case)
        
        print(f"Case {i}:")
        print(f"  ID: {case.get('caseid', 'Unnamed')[:50]}")
        print(f"  Namelists: {', '.join(case['namelists'].keys())}")
        
        # Check if body geometry is defined
        if case_state.get('body_nx', 0) > 0:
            body_props = calculate_body_geometry(case_state)
            print(f"  Body: {body_props['length']:.2f} ft, {body_props['volume']:.2f} ft³")
        else:
            print(f"  Body: Not defined")
        
        print()
    
    assert len(cases) == 4  # EX1.INP has 4 cases
    
    print("  [PASS] Multiple cases handled successfully")


if __name__ == '__main__':
    print("\nPyDATCOM Integration Tests")
    print("=" * 60)
    
    test_parse_and_calculate_body()
    test_full_pipeline()
    test_multiple_cases()
    
    print("\n" + "="*60)
    print("All integration tests passed! [PASS]")
    print("="*60)
    print("\nIntegration Status:")
    print("  [PASS] Parser -> Body geometry")
    print("  [PASS] State manager -> Geometry")
    print("  [PASS] Complete pipeline (Parse -> Calculate -> Export)")
    print("  [PASS] Multi-case handling")
    print("\n  Ready for Phase 3!")
    print("="*60)

