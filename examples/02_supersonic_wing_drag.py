"""
Example 2: Supersonic wave and friction drag of a straight tapered wing.

Exercises the translated SUPDRG paths:

  SUPDRG   Figure 4.1.5.2-58 zero-lift wave drag, sharp and round leading
           edge, with P/RLW geometry and SREF normalization
  SUPDRG   skin friction: exposed-MAC Reynolds number, the source roughness
           cutoff, Figure 4.1.5.1-27 Mach factor and the Mach 3 cap
  CORDSP   supersonic section coordinates for the same airfoil

Run:  python examples/02_supersonic_wing_drag.py
"""

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from pydatcom.aerodynamics.supersonic import (
    calculate_supdrg_straight_wing, calculate_supdrg_skin_friction,
)
from pydatcom.geometry.wing import calculate_straight_exposed_geometry
from pydatcom.geometry.airfoil import NACAGenerator


def build_wing(sharp=True):
    """A thin straight tapered supersonic wing, dimensions in feet."""
    return {
        'wing_type': 1.0,
        'wing_chrdr': 2.90, 'wing_chrdtp': 0.64,
        'wing_sspn': 1.59, 'wing_sspne': 1.59,
        'wing_savsi': 55.0, 'wing_chstat': 0.0,
        'wing_tovc': 0.05,
        'wing_leri': 0.0 if sharp else 0.0134,
        'options_sref': 8.85, 'options_cbarr': 2.46,
        'options_rougfc': 1.6e-4,
    }


def main():
    print("=" * 68)
    print("EXAMPLE 2  Supersonic straight tapered wing drag")
    print("=" * 68)

    state = build_wing(sharp=True)
    geometry = calculate_straight_exposed_geometry(state)

    print("\n-- Exposed planform (WTGEOM) " + "-" * 39)
    print(f"  exposed area        = {geometry['area']:8.4f} ft2")
    print(f"  exposed AR          = {geometry['aspect_ratio']:8.4f}")
    print(f"  exposed MAC         = {geometry['mac']:8.4f} ft")
    print(f"  tan(LE sweep)       = {geometry['tan_le']:8.4f}")
    print(f"  SREF                = {state['options_sref']:8.4f} ft2")

    print("\n-- Section coordinates (CORDSP) " + "-" * 36)
    generator = NACAGenerator(num_points=51)
    for shape in ('double_wedge', 'biconvex', 'hexagonal'):
        kwargs = {'flat_length': 0.2} if shape == 'hexagonal' else {}
        coords = generator.supersonic_airfoil(
            thickness_ratio=0.05, shape=shape, **kwargs)
        full = 2.0 * float(np.max(coords.thickness))
        print(f"  {shape:14} max full thickness t/c = {full:.4f}")

    print("\n-- Zero-lift wave drag and friction vs Mach " + "-" * 24)
    header = (f"{'Mach':>6}{'Re_MAC':>12}{'Re_used':>12}{'Cf':>10}"
              f"{'CD_fric':>10}{'CD_wave':>10}{'CD_0':>10}")
    print(header)
    reynolds_per_ft = 2.0e6
    for mach in (1.2, 1.5, 2.0, 2.5, 3.0, 3.5):
        reynolds = reynolds_per_ft * geometry['mac'] * mach
        wave = calculate_supdrg_straight_wing(state, mach, cl_wing=0.0)
        friction = calculate_supdrg_skin_friction(state, mach, reynolds)
        cd0 = wave['cd_wave_total'] + friction['cd_friction']
        print(f"{mach:6.1f}{reynolds:12.3e}"
              f"{friction['reynolds_used']:12.3e}{friction['cf']:10.5f}"
              f"{friction['cd_friction']:10.5f}"
              f"{wave['cd_wave_total']:10.5f}{cd0:10.5f}")

    print("\n-- Leading edge effect at Mach 2.0 " + "-" * 33)
    for label, sharp in (("sharp  (LERI=0)", True),
                         ("round  (LERI=0.0134)", False)):
        wing = build_wing(sharp=sharp)
        wave = calculate_supdrg_straight_wing(wing, 2.0, cl_wing=0.0)
        print(f"  {label:24} CD_wave = {wave['cd_wave_total']:.6f}"
              f"   ({wave.get('method', 'n/a')})")

    print("\n-- Roughness-limited Reynolds number " + "-" * 31)
    print("  SUPDRG caps the Reynolds number at a surface-roughness cutoff,")
    print("  so a smoother surface admits a higher effective Reynolds number.")
    for roughness in (0.4e-4, 1.6e-4, 6.4e-4):
        wing = build_wing()
        wing['options_rougfc'] = roughness
        friction = calculate_supdrg_skin_friction(wing, 2.0, 1.0e9)
        print(f"  ROUGFC={roughness:8.1e}  cutoff Re = "
              f"{friction['roughness_cutoff_reynolds']:.4e}"
              f"   Cf = {friction['cf']:.5f}")

    print("\nNote: the Mach lookup is capped at 3.0 by the source, so the")
    print("Mach 3.5 row reuses the Mach 3.0 friction factor by design.")


if __name__ == '__main__':
    main()
