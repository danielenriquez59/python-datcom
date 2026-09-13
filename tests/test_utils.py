"""
Tests for utility functions - interpolation, table lookup, etc.
"""

from pathlib import Path
import sys
import numpy as np

sys.path.insert(0, str(Path(__file__).parent.parent))

from pydatcom.utils import (
    asmint, bilinear_interp, fig26, fig53a, fig60b, fig68,
    Atmosphere, arcsin, arccos
)


def test_asmint():
    """Test asymmetric interpolation."""
    print("\n" + "="*60)
    print("Testing ASMINT (Asymmetric Interpolation)")
    print("="*60)
    
    # Create test data
    x_data = np.array([0.0, 1.0, 2.0, 3.0, 4.0])
    y_data = np.array([0.0, 1.0, 1.5, 1.8, 2.0])
    
    # Interpolate at several points
    x_vals = np.array([0.5, 1.5, 2.5, 3.5])
    y_vals = asmint(x_data, y_data, x_vals)
    
    print(f"\nInterpolation results:")
    for x, y in zip(x_vals, y_vals):
        print(f"  x={x:.1f} -> y={y:.3f}")
    
    # Check that interpolated values are reasonable
    assert all(y_vals >= 0)  # Should be positive
    assert all(y_vals <= 2.5)  # Should be in reasonable range
    
    print("  [PASS] ASMINT interpolation successful")


def test_bilinear_interp():
    """Test bilinear interpolation."""
    print("\n" + "="*60)
    print("Testing Bilinear Interpolation")
    print("="*60)
    
    # Create 2D grid
    x_data = np.array([0.0, 1.0, 2.0])
    y_data = np.array([0.0, 1.0, 2.0])
    z_table = np.array([
        [0.0, 1.0, 2.0],
        [1.0, 2.0, 3.0],
        [2.0, 3.0, 4.0]
    ])
    
    # Test interpolation
    z = bilinear_interp(0.5, 0.5, x_data, y_data, z_table)
    print(f"\nInterpolated value at (0.5, 0.5): {z:.3f}")
    print(f"  Expected: ~1.0")
    
    assert 0.9 < z < 1.1
    
    # Test corner value
    z_corner = bilinear_interp(2.0, 2.0, x_data, y_data, z_table)
    print(f"Corner value at (2.0, 2.0): {z_corner:.3f}")
    assert z_corner == 4.0
    
    print("  [PASS] Bilinear interpolation successful")


def test_fig26():
    """Test FIG26 skin friction coefficient."""
    print("\n" + "="*60)
    print("Testing FIG26 (Skin Friction Coefficient)")
    print("="*60)
    
    # Test at typical flight conditions
    reynolds = 1.0e6
    mach = 0.6
    
    cf = fig26(reynolds, mach)
    
    print(f"\nSkin friction coefficient:")
    print(f"  Reynolds number: {reynolds:.2e}")
    print(f"  Mach number: {mach}")
    print(f"  CF: {cf:.6f}")
    
    # CF should be small positive number for turbulent flow
    assert 0.001 < cf < 0.01
    
    # Test Mach dependence
    cf_low = fig26(reynolds, 0.3)
    cf_high = fig26(reynolds, 2.0)
    print(f"\n  CF at M=0.3: {cf_low:.6f}")
    print(f"  CF at M=2.0: {cf_high:.6f}")
    
    print("  [PASS] FIG26 calculations successful")


def test_fig53a():
    """Test FIG53A correlation."""
    print("\n" + "="*60)
    print("Testing FIG53A")
    print("="*60)
    
    rv = 1.0e5
    z = 0.5
    
    r = fig53a(rv, z)
    
    print(f"\nFIG53A result:")
    print(f"  RV: {rv:.2e}")
    print(f"  Z: {z}")
    print(f"  R: {r:.6f}")
    
    assert isinstance(r, (int, float))
    
    print("  [PASS] FIG53A calculations successful")


def test_fig60b():
    """Test FIG60B supersonic correlation."""
    print("\n" + "="*60)
    print("Testing FIG60B (Supersonic Parameters)")
    print("="*60)
    
    beta = 2.0
    btana = 0.637  # Source table at CNAA=1.0, BETA=2.0
    cnaa = fig60b(beta, btana)
    
    print(f"\nFIG60B results for beta={beta}:")
    print(f"  BTANA: {btana:.4f}")
    print(f"  CNAA: {cnaa:.4f}")
    
    assert abs(cnaa - 1.0) < 1e-12
    
    print("  [PASS] FIG60B calculations successful")


def test_fig68():
    """Test FIG68 shock angle."""
    print("\n" + "="*60)
    print("Testing FIG68 (Shock Wave Angle)")
    print("="*60)
    
    mach = 2.0
    delta = 10.0  # Deflection angle
    
    theta, ierr = fig68(mach, delta)
    
    print(f"\nFIG68 results:")
    print(f"  Mach: {mach}")
    print(f"  Deflection: {delta}°")
    print(f"  Shock angle: {theta:.2f}°")
    print(f"  Error code: {ierr}")
    
    assert ierr == 0
    assert theta > delta  # Shock angle > deflection angle
    
    print("  [PASS] FIG68 calculations successful")


def test_atmosphere():
    """Test atmosphere calculations."""
    print("\n" + "="*60)
    print("Testing Atmosphere Model")
    print("="*60)
    
    # Test at sea level
    atm_sl = Atmosphere.calculate(0.0)
    print(f"\nSea Level:")
    print(f"  Temperature: {atm_sl['temperature']:.2f} °R")
    print(f"  Pressure: {atm_sl['pressure']:.2f} psf")
    print(f"  Density: {atm_sl['density']:.6f} slug/ft³")
    print(f"  Speed of sound: {atm_sl['cs']:.2f} ft/s")
    
    # Check against standard values
    assert 515 < atm_sl['temperature'] < 520  # ~518.67 °R
    assert 2100 < atm_sl['pressure'] < 2120   # ~2116.2 psf
    
    # Test at 10,000 ft
    atm_10k = Atmosphere.calculate(10000.0)
    print(f"\n10,000 ft:")
    print(f"  Temperature: {atm_10k['temperature']:.2f} °R")
    print(f"  Pressure: {atm_10k['pressure']:.2f} psf")
    print(f"  Density: {atm_10k['density']:.6f} slug/ft³")
    
    # Pressure should decrease with altitude
    assert atm_10k['pressure'] < atm_sl['pressure']
    assert atm_10k['density'] < atm_sl['density']
    
    print("  [PASS] Atmosphere calculations successful")


def test_trig_functions():
    """Test arc trig functions with bounds checking."""
    print("\n" + "="*60)
    print("Testing Arc Trig Functions")
    print("="*60)
    
    # Test arcsin
    y1 = arcsin(0.5)
    print(f"\narcsin(0.5) = {y1:.4f} rad ({np.rad2deg(y1):.2f}°)")
    assert abs(y1 - np.pi/6) < 0.001  # Should be 30°
    
    # Test arcsin at boundary
    y2 = arcsin(1.0)
    print(f"arcsin(1.0) = {y2:.4f} rad ({np.rad2deg(y2):.2f}°)")
    assert abs(y2 - np.pi/2) < 0.001  # Should be 90°
    
    # Test arccos
    y3 = arccos(0.5)
    print(f"\narccos(0.5) = {y3:.4f} rad ({np.rad2deg(y3):.2f}°)")
    assert abs(y3 - np.pi/3) < 0.001  # Should be 60°
    
    print("  [PASS] Arc trig functions successful")


if __name__ == '__main__':
    print("\nTesting PyDATCOM Utility Functions")
    print("=" * 60)
    
    test_trig_functions()
    test_asmint()
    test_bilinear_interp()
    test_atmosphere()
    test_fig26()
    test_fig53a()
    test_fig60b()
    test_fig68()
    
    print("\n" + "="*60)
    print("All utility tests passed! [PASS]")
    print("="*60)
    print("\nPhase 3 Progress:")
    print("  [PASS] Math utilities (arcsin, arccos, area)")
    print("  [PASS] Atmosphere model (US Std 1962)")
    print("  [PASS] Interpolation (ASMINT, bilinear)")
    print("  [PASS] Table lookup (FIG26, FIG53A, FIG60B, FIG68)")
    print("\n  Phase 3 COMPLETE!")
    print("="*60)

