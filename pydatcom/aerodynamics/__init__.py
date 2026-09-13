"""Aerodynamic calculations for PyDATCOM - lift, drag, moments, stability."""

from pydatcom.aerodynamics.calculator import AerodynamicCalculator, calculate_aero_coefficients
from pydatcom.aerodynamics.lift import LiftCalculator, calculate_wing_lift_subsonic
from pydatcom.aerodynamics.drag import DragCalculator, calculate_total_drag
from pydatcom.aerodynamics.moment import (
    MomentCalculator, calculate_total_pitching_moment,
    calculate_cmalph_zero_lift_moment,
)
from pydatcom.aerodynamics.subsonic import calculate_subsonic_coefficients
from pydatcom.aerodynamics.transonic import calculate_transonic_coefficients
from pydatcom.aerodynamics.supersonic import (
    calculate_supersonic_coefficients, calculate_supdrg_straight_wing,
    calculate_supdrg_skin_friction,
)
from pydatcom.aerodynamics.hypersonic import calculate_hypersonic_coefficients
from pydatcom.aerodynamics.stability import StabilityCalculator, calculate_all_stability_derivatives
from pydatcom.aerodynamics.body_alone import has_wing_or_tail, calculate_body_alone_coefficients
from pydatcom.aerodynamics.bodyrt import calculate_bodyrt
from pydatcom.aerodynamics.downwash import (
    calculate_downwash, calculate_downwash_geometry,
    calculate_downwash_gradient_441, fig4417_68a, fig4417_68b,
)
from pydatcom.aerodynamics.wing_body_tail import (
    calculate_clwbt, calculate_cdwbt, calculate_tail_load,
)

__all__ = [
    'AerodynamicCalculator',
    'calculate_aero_coefficients',
    'LiftCalculator',
    'DragCalculator',
    'MomentCalculator',
    'StabilityCalculator',
    'calculate_wing_lift_subsonic',
    'calculate_total_drag',
    'calculate_total_pitching_moment',
    'calculate_cmalph_zero_lift_moment',
    'calculate_all_stability_derivatives',
    'calculate_subsonic_coefficients',
    'calculate_transonic_coefficients',
    'calculate_supersonic_coefficients',
    'calculate_supdrg_straight_wing',
    'calculate_supdrg_skin_friction',
    'calculate_hypersonic_coefficients',
    'has_wing_or_tail',
    'calculate_body_alone_coefficients',
    'calculate_bodyrt',
    'calculate_downwash',
    'calculate_downwash_geometry',
    'calculate_downwash_gradient_441',
    'fig4417_68a',
    'fig4417_68b',
    'calculate_clwbt',
    'calculate_cdwbt',
    'calculate_tail_load',
]

