"""
Example 1: Longitudinal analysis of a conventional wing-tail aircraft.

Exercises the translated routines that make a whole-aircraft longitudinal
result possible:

  WTGEOM   exposed and theoretical planform geometry (both surfaces)
  INFTGM   downwash synthesizing dimensions: tail arm and tail height
  DWASH    DATCOM Section 4.4.1 downwash gradient at the tail
  4.3.1.2  wing/tail-body lift carryover factors
  CLWBT    wing-body-tail lift buildup
  CMALPH   wing zero-lift pitching moment

Run:  python examples/01_wing_tail_longitudinal.py
"""

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from pydatcom.geometry.wing import calculate_straight_exposed_geometry
from pydatcom.aerodynamics.downwash import (
    calculate_downwash, calculate_downwash_geometry,
)
from pydatcom.interactions.carryover import calculate_carryover_factors
from pydatcom.aerodynamics.moment import calculate_total_pitching_moment
from pydatcom.aerodynamics.lift import calculate_wing_lift_subsonic


def build_aircraft():
    """A light twin-class straight-tapered wing with a conventional tail.

    Dimensions are in feet.  The wing is a 30 ft span, 6 ft root / 3 ft tip
    trapezoid with an unswept quarter chord; the tail is a 12 ft span
    geometrically similar surface on a 25 ft arm.
    """
    return {
        # Wing planform: WGPLNF
        'wing_type': 1.0,
        'wing_chrdr': 6.0, 'wing_chrdtp': 3.0,
        'wing_sspn': 15.0, 'wing_sspne': 13.5,
        'wing_savsi': 0.0, 'wing_chstat': 0.25,
        'wing_cmo': -0.02, 'wing_twista': 0.0,
        # Horizontal tail: HTPLNF
        'htail_type': 1.0,
        'htail_chrdr': 3.0, 'htail_chrdtp': 1.5,
        'htail_sspn': 6.0, 'htail_sspne': 5.4,
        'htail_savsi': 0.0, 'htail_chstat': 0.25,
        # Synthesis: SYNTHS
        'synths_xw': 10.0, 'synths_zw': 0.0, 'synths_aliw': 0.0,
        'synths_xh': 35.0, 'synths_zh': 2.0, 'synths_alih': -1.0,
        'synths_xcg': 12.0,
        # Reference quantities: OPTINS
        'options_sref': 135.0, 'options_cbarr': 4.6667,
        'flight_mach': 0.3,
    }


def main():
    state = build_aircraft()
    mach = state['flight_mach']

    wing = calculate_straight_exposed_geometry(state, component='wing')
    tail = calculate_straight_exposed_geometry(state, component='htail')
    geometry = calculate_downwash_geometry(state)
    carryover = calculate_carryover_factors(state, component='htail')
    downwash = calculate_downwash(state, 0.0)

    print("=" * 68)
    print("EXAMPLE 1  Conventional wing-tail aircraft, longitudinal analysis")
    print("=" * 68)

    print("\n-- Planform geometry (WTGEOM) " + "-" * 38)
    print(f"{'':22}{'wing':>12}{'tail':>12}")
    for label, key in (("exposed area  [ft2]", 'area'),
                       ("exposed AR", 'aspect_ratio'),
                       ("exposed taper", 'taper_ratio'),
                       ("exposed MAC   [ft]", 'mac'),
                       ("MAC y-station [ft]", 'mac_span_location'),
                       ("theoretical AR", 'aspect_ratio_theoretical')):
        print(f"{label:22}{wing[key]:12.4f}{tail[key]:12.4f}")
    print(f"{'SREF          [ft2]':22}{state['options_sref']:12.4f}")
    print(f"{'CBARR         [ft]':22}{state['options_cbarr']:12.4f}")

    print("\n-- Downwash geometry (INFTGM) " + "-" * 38)
    print(f"  tail arm      A(24) = {geometry['tail_arm']:8.4f} ft")
    print(f"  tail height   A(12) = {geometry['tail_height']:8.4f} ft")
    print(f"  tail angle    A(11) = {np.rad2deg(geometry['tail_angle']):8.4f} deg")

    print("\n-- Downwash gradient (DWASH Section 4.4.1) " + "-" * 25)
    print(f"  K_A                 = {downwash['k_a']:8.4f}")
    print(f"  K_lambda            = {downwash['k_lambda']:8.4f}")
    print(f"  K_H                 = {downwash['k_h']:8.4f}")
    print(f"  de/da               = {downwash['deda']:8.4f}")

    print("\n-- Lift carryover (Section 4.3.1.2) " + "-" * 32)
    print(f"  r/s                 = {carryover['ratio']:8.4f}")
    print(f"  K_H(B)  KWB         = {carryover['kwb']:8.4f}")
    print(f"  K_B(H)  KBW         = {carryover['kbw']:8.4f}")
    print(f"  sum vs (1+r/s)^2    = {carryover['kwb'] + carryover['kbw']:8.4f}"
          f" vs {(1 + carryover['ratio'])**2:.4f}")

    print("\n-- Longitudinal sweep " + "-" * 46)
    header = (f"{'alpha':>7}{'eps':>9}{'alpha_t':>9}{'CL_wing':>10}"
              f"{'Cm_wing':>10}{'Cm_tail':>10}{'Cm_tot':>10}")
    print(header)
    print(f"{'[deg]':>7}{'[deg]':>9}{'[deg]':>9}{'':>10}{'':>10}{'':>10}{'':>10}")
    rows = []
    for alpha in (-4.0, -2.0, 0.0, 2.0, 4.0, 6.0, 8.0, 10.0):
        lift = calculate_wing_lift_subsonic(state, alpha, mach)
        cl = lift['cl']
        moment = calculate_total_pitching_moment(state, cl, alpha, mach)
        eps = moment.get('eps_deg', 0.0)
        alpha_t = alpha - eps + state['synths_alih']
        rows.append((alpha, cl, moment['cm_total']))
        print(f"{alpha:7.1f}{eps:9.3f}{alpha_t:9.3f}{cl:10.4f}"
              f"{moment['cm_wing']:10.4f}{moment['cm_tail']:10.4f}"
              f"{moment['cm_total']:10.4f}")

    print("\n-- Static longitudinal stability " + "-" * 35)
    alphas = np.array([r[0] for r in rows])
    cls = np.array([r[1] for r in rows])
    cms = np.array([r[2] for r in rows])
    cla = np.polyfit(alphas, cls, 1)[0]
    cma = np.polyfit(alphas, cms, 1)[0]
    dcm_dcl = np.polyfit(cls, cms, 1)[0]
    cbar = state['options_cbarr']
    print(f"  CL_alpha            = {cla:8.5f} /deg")
    print(f"  Cm_alpha            = {cma:8.5f} /deg")
    print(f"  dCm/dCL             = {dcm_dcl:8.5f}")
    print(f"  static margin       = {-dcm_dcl * 100:8.3f} % MAC")
    xnp = state['synths_xcg'] - dcm_dcl * cbar
    print(f"  neutral point       = {xnp:8.4f} ft from nose datum")
    print(f"  CG                  = {state['synths_xcg']:8.4f} ft")
    verdict = "STATICALLY STABLE" if cma < 0 else "STATICALLY UNSTABLE"
    print(f"\n  verdict: {verdict} (Cm_alpha {'<' if cma < 0 else '>='} 0)")

    print("\n-- Tail sizing check " + "-" * 47)
    print("  Cm_tail with the tail removed would be 0.0000 by construction;")
    print(f"  the translated buildup gives {rows[-1][2]:.4f} total Cm at "
          f"{rows[-1][0]:.0f} deg.")
    print("\nNote: wing/body zero-lift drag and the body load are not part of")
    print("this longitudinal summary; see PHYSICS_REVIEW.md and")
    print("TRANSLATION_STATUS.md for what remains approximate.")


if __name__ == '__main__':
    main()
