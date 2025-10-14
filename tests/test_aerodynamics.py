"""
Tests for aerodynamic calculations - lift, drag, moments.
"""

from pathlib import Path
import sys
import numpy as np
from typing import Dict

sys.path.insert(0, str(Path(__file__).parent.parent))

from pydatcom.aerodynamics import (
    AerodynamicCalculator,
    DragCalculator,
    StabilityCalculator,
    calculate_aero_coefficients,
    calculate_subsonic_coefficients,
    calculate_supersonic_coefficients,
    calculate_transonic_coefficients,
    calculate_hypersonic_coefficients,
)
from pydatcom.io import StateManager


def create_test_state() -> Dict:
    """Create test state with typical aircraft geometry."""
    return {
        'wing_aspect_ratio': 6.0,
        'wing_taper_ratio': 0.5,
        'wing_savsi': 25.0,  # 25 deg sweep
        'wing_tovc': 0.12,   # 12% thick
        'wing_alphai': 0.0,
        'wing_span': 50.0,
        'wing_area': 375.0,
        'wing_mac': 7.78,
        'options_sref': 375.0,
        'options_cbarr': 7.78,
        'options_blref': 50.0,
        'synths_xcg': 15.0,
        'synths_xw': 12.0,
        'body_length': 40.0,
        'body_max_area': 10.0,
        'body_centroid': 18.0,
        'flight_rnnub': [5e6],
    }


def test_subsonic_calculations():
    """Test subsonic aerodynamic calculations."""
    print("\n" + "="*60)
    print("Testing Subsonic Aerodynamics")
    print("="*60)
    
    state = create_test_state()
    
    # Calculate at M=0.6, alpha=5°
    result = calculate_subsonic_coefficients(state, alpha_deg=5.0, mach=0.6, reynolds=5e6)
    
    print(f"\nSubsonic Results (M=0.6, alpha=5°):")
    print(f"  CL: {result['cl']:.4f}")
    print(f"  CD: {result['cd']:.4f}")
    print(f"  Cm: {result['cm']:.4f}")
    print(f"  CL_alpha: {result['cla']:.4f} /deg")
    print(f"  CD_friction: {result['cd_friction']:.4f}")
    print(f"  CD_induced: {result['cd_induced']:.4f}")
    
    # Sanity checks
    assert 0.3 < result['cl'] < 0.8  # Reasonable CL for 5 deg
    assert 0.01 < result['cd'] < 0.08  # Reasonable CD
    assert abs(result['cm']) < 0.5  # Reasonable moment
    assert result['cla'] > 0.04  # Positive lift slope
    
    print("  [PASS] Subsonic calculations successful")


def test_supersonic_calculations():
    """Test supersonic aerodynamic calculations."""
    print("\n" + "="*60)
    print("Testing Supersonic Aerodynamics")
    print("="*60)
    
    state = create_test_state()
    state['wing_tovc'] = 0.08  # Thinner for supersonic
    
    # Calculate at M=2.0, alpha=5°
    result = calculate_supersonic_coefficients(state, alpha_deg=5.0, mach=2.0, reynolds=3e6)
    
    print(f"\nSupersonic Results (M=2.0, alpha=5°):")
    print(f"  CL: {result['cl']:.4f}")
    print(f"  CD: {result['cd']:.4f}")
    print(f"  Cm: {result['cm']:.4f}")
    print(f"  CL_alpha: {result['cla']:.4f} /deg")
    print(f"  Beta: {result['beta']:.4f}")
    print(f"  CD_wave_volume: {result['cd_wave_volume']:.4f}")
    print(f"  CD_wave_lift: {result['cd_wave_lift']:.4f}")
    
    assert 0.15 < result['cl'] < 0.6  # Lower slope in supersonic
    assert 0.01 < result['cd'] < 0.15
    assert result['beta'] > 1.0
    
    print("  [PASS] Supersonic calculations successful")


def test_transonic_calculations():
    """Test transonic aerodynamic calculations."""
    print("\n" + "="*60)
    print("Testing Transonic Aerodynamics")
    print("="*60)
    
    state = create_test_state()
    
    # Calculate at M=1.05, alpha=3°
    result = calculate_transonic_coefficients(state, alpha_deg=3.0, mach=1.05, reynolds=4e6)
    
    print(f"\nTransonic Results (M=1.05, alpha=3°):")
    print(f"  CL: {result['cl']:.4f}")
    print(f"  CD: {result['cd']:.4f}")
    print(f"  Cm: {result['cm']:.4f}")
    print(f"  CD_divergence: {result['cd_divergence']:.4f}")
    
    assert 0.1 < result['cl'] < 0.5
    assert 0.02 < result['cd'] < 0.2  # Higher drag in transonic
    
    print("  [PASS] Transonic calculations successful")


def test_hypersonic_calculations():
    """Test hypersonic aerodynamic calculations."""
    print("\n" + "="*60)
    print("Testing Hypersonic Aerodynamics")
    print("="*60)
    
    state = create_test_state()
    state['wing_tovc'] = 0.05  # Very thin for hypersonic
    
    # Calculate at M=6.0, alpha=10°
    result = calculate_hypersonic_coefficients(state, alpha_deg=10.0, mach=6.0)
    
    print(f"\nHypersonic Results (M=6.0, alpha=10°):")
    print(f"  CL: {result['cl']:.4f}")
    print(f"  CD: {result['cd']:.4f}")
    print(f"  Cm: {result['cm']:.4f}")
    print(f"  CN: {result['cn']:.4f}")
    print(f"  CP_max: {result['cp_max']:.4f}")
    
    assert result['cn'] > 0  # Positive normal force
    assert result['cd'] > 0  # Positive drag
    
    print("  [PASS] Hypersonic calculations successful")


def test_calculator_routing():
    """Test that AerodynamicCalculator routes to correct regimes."""
    print("\n" + "="*60)
    print("Testing Regime Routing")
    print("="*60)
    
    state = create_test_state()
    calc = AerodynamicCalculator(state)
    
    # Test regime identification
    test_cases = [
        (0.5, 'subsonic'),
        (0.85, 'subsonic'),
        (1.0, 'transonic'),
        (1.1, 'transonic'),
        (2.5, 'supersonic'),
        (6.0, 'hypersonic'),
    ]
    
    print("\nRegime Identification:")
    for mach, expected_regime in test_cases:
        regime = calc.identify_regime(mach)
        print(f"  M={mach:.1f} -> {regime:12s} (expected: {expected_regime})")
        assert regime == expected_regime
    
    print("\n  [PASS] Regime routing successful")


def test_alpha_sweep():
    """Test coefficient calculation across alpha range."""
    print("\n" + "="*60)
    print("Testing Alpha Sweep")
    print("="*60)
    
    state = create_test_state()
    calc = AerodynamicCalculator(state)
    
    # Subsonic alpha sweep
    alpha_range = np.array([-5, 0, 5, 10, 15])
    result = calc.calculate_alpha_sweep(alpha_range, mach=0.6)
    
    print(f"\nAlpha Sweep (M=0.6):")
    print(f"  Alpha (deg): {alpha_range}")
    print(f"  CL:          {result['cl']}")
    print(f"  CD:          {result['cd']}")
    
    # Check that CL increases with alpha
    assert np.all(np.diff(result['cl']) > 0)
    
    # Check that CD increases with |alpha| (roughly)
    assert result['cd'][2] > result['cd'][1]  # Higher at 5° than 0°
    
    print("  [PASS] Alpha sweep successful")


def test_mach_sweep():
    """Test coefficient calculation across Mach range."""
    print("\n" + "="*60)
    print("Testing Mach Sweep")
    print("="*60)
    
    state = create_test_state()
    calc = AerodynamicCalculator(state)
    
    # Mach sweep at alpha=5°
    mach_range = np.array([0.3, 0.6, 0.85, 1.5, 2.5])
    result = calc.calculate_mach_sweep(alpha_deg=5.0, mach_range=mach_range)
    
    print(f"\nMach Sweep (alpha=5°):")
    print(f"  Mach:  {result['mach']}")
    print(f"  CL:    {result['cl']}")
    print(f"  CD:    {result['cd']}")
    
    # CD should increase dramatically in transonic region
    assert result['cd'][2] > result['cd'][1]  # M=0.85 > M=0.6
    
    print("  [PASS] Mach sweep successful")


def test_drag_polar():
    """Test drag polar generation."""
    print("\n" + "="*60)
    print("Testing Drag Polar")
    print("="*60)
    
    state = create_test_state()
    calc = AerodynamicCalculator(state)
    drag_calc = DragCalculator(state)
    
    # Generate drag polar
    cl_range = np.linspace(0, 1.5, 11)
    polar = drag_calc.calculate_drag_polar(mach=0.6, reynolds=5e6, cl_range=cl_range)
    
    print(f"\nDrag Polar (M=0.6):")
    print(f"  CL range: {cl_range[0]:.1f} to {cl_range[-1]:.1f}")
    print(f"  CD range: {polar['cd'].min():.4f} to {polar['cd'].max():.4f}")
    
    # CD should increase with CL (induced drag)
    assert polar['cd'][-1] > polar['cd'][0]
    
    # Find minimum drag
    cd_min_idx = np.argmin(polar['cd'])
    print(f"  Minimum CD: {polar['cd'][cd_min_idx]:.4f} at CL={polar['cl'][cd_min_idx]:.2f}")
    
    print("  [PASS] Drag polar generation successful")


def test_stability_derivatives():
    """Test stability derivative calculations."""
    print("\n" + "="*60)
    print("Testing Stability Derivatives")
    print("="*60)
    
    state = create_test_state()
    # Add tail for stability
    state['htail_area'] = 75.0
    state['vtail_area'] = 40.0
    state['synths_xh'] = 35.0
    state['synths_xv'] = 36.0
    
    stab_calc = StabilityCalculator(state)
    derivs = stab_calc.calculate_derivatives(mach=0.6)
    
    print(f"\nStability Derivatives (M=0.6):")
    print(f"  Cmq (pitch damping): {derivs['cmq']:.4f} /rad")
    print(f"  Clp (roll damping): {derivs['clp']:.4f} /rad")
    print(f"  Cnr (yaw damping): {derivs['cnr']:.4f} /rad")
    print(f"  Cn_beta (directional): {derivs['cn_beta']:.4f} /rad")
    print(f"  Static margin: {derivs['static_margin']*100:.1f}%")
    
    # Check for stability
    assessment = stab_calc.assess_stability(mach=0.6)
    print(f"\nStability Assessment:")
    print(f"  Longitudinal: {assessment['longitudinal']}")
    print(f"  Directional: {assessment['directional']}")
    print(f"  Overall: {assessment['overall']}")
    
    assert derivs['cmq'] < 0  # Pitch damping should be negative
    assert derivs['cn_beta'] > 0  # Directionally stable
    
    print("  [PASS] Stability calculations successful")


if __name__ == '__main__':
    print("\nTesting PyDATCOM Aerodynamic Calculations")
    print("=" * 60)
    
    test_subsonic_calculations()
    test_supersonic_calculations()
    test_transonic_calculations()
    test_hypersonic_calculations()
    test_calculator_routing()
    test_alpha_sweep()
    test_mach_sweep()
    test_drag_polar()
    test_stability_derivatives()
    
    print("\n" + "="*60)
    print("All aerodynamic tests passed! [PASS]")
    print("="*60)
    print("\nPhase 4 Progress:")
    print("  [PASS] Subsonic calculations")
    print("  [PASS] Supersonic calculations")
    print("  [PASS] Transonic calculations")
    print("  [PASS] Hypersonic calculations")
    print("  [PASS] Regime routing")
    print("  [PASS] Alpha/Mach sweeps")
    print("  [PASS] Drag polar generation")
    print("  [PASS] Stability derivatives")
    print("\n  Phase 4 100% COMPLETE!")
    print("="*60)

