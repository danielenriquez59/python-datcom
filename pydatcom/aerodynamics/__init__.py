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
from pydatcom.aerodynamics.dynbod import (
    calculate_dynbod_subsonic, calculate_dynbod_hypersonic,
)
from pydatcom.aerodynamics.dynamic_buildup import (
    calculate_dnpawb, calculate_dnpwbt,
)
from pydatcom.aerodynamics.clrder import (
    calculate_clr_wing, calculate_clr_panel_increment,
)
from pydatcom.aerodynamics.vertical_panel import (
    calculate_vertical_panel_cdo, calculate_vertical_panel_drag,
)
from pydatcom.aerodynamics.vortex_core import calculate_sddvc
from pydatcom.aerodynamics.cnca import calculate_cnca
from pydatcom.aerodynamics.maxcl import calculate_maxcl
from pydatcom.aerodynamics.calca0 import calculate_calca0
from pydatcom.aerodynamics.wbtcdo import (
    calculate_wbtcdo, calculate_drag_divergence_mach,
)
from pydatcom.aerodynamics.subpaw import calculate_subpaw
from pydatcom.aerodynamics.subwbt import calculate_subwbt
from pydatcom.aerodynamics.sssym import calculate_sssym
from pydatcom.aerodynamics.tablec import calculate_tablec
from pydatcom.aerodynamics.tbsub import calculate_tbsub
from pydatcom.aerodynamics.tbsup import calculate_tbsup
from pydatcom.aerodynamics.tbtrn import calculate_tbtrn
from pydatcom.aerodynamics.vertical_lift_figures import (
    fig4132_56a, fig4132_56g, fig4132_60a, fig4132_60b,
    fig4132_61, fig4132_62, fig4132_63,
)
from pydatcom.aerodynamics.wingcl import (
    calculate_wingcl, calculate_wingcl_clb, calculate_wingcl_cdl,
)
from pydatcom.aerodynamics.supersonic_downwash import (
    calculate_sdwa, calculate_sdwb, calculate_sdwc,
    calculate_sdwd, calculate_sdwe,
)
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
    'calculate_dynbod_subsonic',
    'calculate_dynbod_hypersonic',
    'calculate_dnpawb',
    'calculate_dnpwbt',
    'calculate_clr_wing',
    'calculate_clr_panel_increment',
    'calculate_vertical_panel_cdo',
    'calculate_vertical_panel_drag',
    'calculate_sddvc',
    'calculate_cnca',
    'calculate_maxcl',
    'calculate_calca0',
    'calculate_wbtcdo',
    'calculate_drag_divergence_mach',
    'calculate_subpaw',
    'calculate_subwbt',
    'calculate_sssym',
    'calculate_tablec',
    'calculate_tbsub',
    'calculate_tbsup',
    'calculate_tbtrn',
    'fig4132_56a', 'fig4132_56g',
    'fig4132_60a', 'fig4132_60b',
    'fig4132_61', 'fig4132_62', 'fig4132_63',
    'calculate_wingcl',
    'calculate_wingcl_clb',
    'calculate_wingcl_cdl',
    'calculate_sdwa',
    'calculate_sdwb',
    'calculate_sdwc',
    'calculate_sdwd',
    'calculate_sdwe',
    'calculate_downwash',
    'calculate_downwash_geometry',
    'calculate_downwash_gradient_441',
    'fig4417_68a',
    'fig4417_68b',
    'calculate_clwbt',
    'calculate_cdwbt',
    'calculate_tail_load',
]

