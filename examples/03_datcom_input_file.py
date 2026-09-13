"""
Example 3: Run an original Digital DATCOM input file end to end.

Parses the AFFDL-TR-79-3032 Example Problem 2 input deck shipped in
tests/fixtures/ex2.inp and runs the translated pipeline over its flight
condition schedule.  This exercises the namelist parser (including the
unindexed array continuations that the MACH and ALSCHD cards use), the
state manager, and the regime-dispatching aerodynamic calculator.

Run:  python examples/03_datcom_input_file.py
"""

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from pydatcom.io import NamelistParser, StateManager
from pydatcom.aerodynamics.calculator import AerodynamicCalculator

FIXTURE = Path(__file__).resolve().parent.parent / 'tests' / 'fixtures' / 'ex2.inp'


def main():
    print("=" * 68)
    print("EXAMPLE 3  Digital DATCOM input deck, Example Problem 2")
    print("=" * 68)
    print(f"\ninput: {FIXTURE}")

    parser = NamelistParser()
    cases = parser.parse_file(str(FIXTURE))
    print(f"parsed {len(cases)} case(s)")

    case = cases[0]
    state = StateManager()
    state.update(parser.to_state_dict(case))

    mach_list = state.get('flight_mach')
    alphas = state.get('flight_alschd')
    altitudes = state.get('flight_alt')
    mach_list = np.atleast_1d(mach_list)
    alphas = np.atleast_1d(alphas)
    altitudes = np.atleast_1d(altitudes)

    print("\n-- Flight condition schedule (FLTCON) " + "-" * 30)
    print(f"  Mach numbers        : {', '.join(f'{m:g}' for m in mach_list)}")
    print(f"  altitudes      [ft] : {', '.join(f'{a:g}' for a in altitudes)}")
    print(f"  angle schedule [deg]: {', '.join(f'{a:g}' for a in alphas)}")
    print(f"  ({len(alphas)} angles, {len(mach_list)} Mach numbers)")

    print("\n-- Reference quantities (OPTINS) " + "-" * 35)
    for label, key in (("SREF  [ft2]", 'options_sref'),
                       ("CBARR [ft]", 'options_cbarr'),
                       ("BLREF [ft]", 'options_blref')):
        print(f"  {label:14}= {state.get(key)}")

    print("\n-- Wing planform (WGPLNF) " + "-" * 42)
    for label, key in (("CHRDR  [ft]", 'wing_chrdr'),
                       ("CHRDTP [ft]", 'wing_chrdtp'),
                       ("SSPN   [ft]", 'wing_sspn'),
                       ("SSPNE  [ft]", 'wing_sspne'),
                       ("SAVSI [deg]", 'wing_savsi'),
                       ("CHSTAT", 'wing_chstat')):
        print(f"  {label:14}= {state.get(key)}")

    calculator = AerodynamicCalculator(state.to_dict()
                                       if hasattr(state, 'to_dict')
                                       else state._state)

    print("\n-- Regime dispatch across the Mach schedule " + "-" * 24)
    for mach in mach_list:
        print(f"  Mach {mach:5.2f}  ->  {calculator.identify_regime(float(mach))}")

    print("\n-- Coefficients at each Mach, alpha = 4 deg " + "-" * 24)
    print(f"{'Mach':>7}{'regime':>12}{'CL':>10}{'CD':>10}{'Cm':>10}")
    for mach in mach_list:
        mach = float(mach)
        try:
            result = calculator.calculate_at_condition(4.0, mach)
            print(f"{mach:7.2f}{calculator.identify_regime(mach):>12}"
                  f"{result.get('cl', float('nan')):10.4f}"
                  f"{result.get('cd', float('nan')):10.5f}"
                  f"{result.get('cm', float('nan')):10.4f}")
        except Exception as error:  # noqa: BLE001 - example diagnostics
            print(f"{mach:7.2f}{calculator.identify_regime(mach):>12}"
                  f"   -- not available: {type(error).__name__}: {error}")

    print("\n-- Angle sweep at the first Mach number " + "-" * 28)
    mach = float(mach_list[0])
    print(f"{'alpha':>8}{'CL':>10}{'CD':>10}{'Cm':>10}")
    for alpha in alphas:
        try:
            result = calculator.calculate_at_condition(float(alpha), mach)
            print(f"{float(alpha):8.1f}{result.get('cl', float('nan')):10.4f}"
                  f"{result.get('cd', float('nan')):10.5f}"
                  f"{result.get('cm', float('nan')):10.4f}")
        except Exception as error:  # noqa: BLE001 - example diagnostics
            print(f"{float(alpha):8.1f}   -- {type(error).__name__}")

    print("\nThis deck is an exposed-wing-alone case, so there is no tail or")
    print("body buildup here; see examples/01 for the wing-tail path.")
    print("Whole-aircraft values remain subject to the approximations")
    print("listed in PHYSICS_REVIEW.md and TRANSLATION_STATUS.md.")


if __name__ == '__main__':
    main()
