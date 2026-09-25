"""
Supersonic aerodynamic calculations for PyDATCOM.

Implements supersonic flow methods (Mach > 1.2):
- Linearized supersonic theory
- Shock-expansion method
- Prandtl-Meyer relations

Reference: datcom.f supersonic calculation sections
"""

import math
import numpy as np
from typing import Dict
import logging

from pydatcom.utils.table_lookup import fig26, fig60b, fig68
from pydatcom.aerodynamics.lift import resolve_wing_lift_inputs, wing_reference_ratio
from pydatcom.utils.constants import UNUSED
from pydatcom.utils.legacy_tables import tlin1x
from pydatcom.utils.legacy_numeric import tbfunx
from pydatcom.geometry.wing import calculate_straight_exposed_geometry

logger = logging.getLogger(__name__)

_FIG_415258_X = np.array([0.0, 0.28, 0.40, 0.50, 1.87])
_FIG_415258_SHARP = np.array([0.54, 0.54, 0.559, 0.60, 1.95])
_FIG_415258_ROUND = np.array([0.54, 0.54, 0.593, 0.68, 2.0])
_FIG_415127_MACH = np.array([0.0, 1.0, 2.0, 3.0])
_FIG_415127_CEPT = np.array([1.57780, 1.67221, 1.98509, 2.28874])


def calculate_supdrg_skin_friction(state: Dict, mach: float,
                                   reynolds_mac: float) -> Dict[str, float]:
    """Translate SUPDRG's straight-wing skin-friction branch."""
    if mach <= 1.0 or reynolds_mac <= 0.0:
        raise ValueError("SUPDRG friction requires Mach > 1 and positive Reynolds")
    geometry = calculate_straight_exposed_geometry(state)
    sref = float(state.get('options_sref', geometry['area']) or geometry['area'])
    if sref <= 0.0:
        raise ValueError("SREF must be positive")

    reynolds_used = float(reynolds_mac)
    roughness = float(state.get('options_rougfc', 1.6e-4) or 1.6e-4)
    if roughness < 0.0:
        raise ValueError("roughness height cannot be negative")
    mach_lookup = min(float(mach), 3.0)
    cept, _ = tbfunx(_FIG_415127_MACH, _FIG_415127_CEPT,
                     mach_lookup, lower=0, upper=0)
    cutoff = (12.0 * geometry['mac'] / roughness)**1.0482 * 10.0**cept
    reynolds_used = min(reynolds_used, cutoff)

    cf = fig26(reynolds_used, mach_lookup)
    cdf = 2.0 * cf * geometry['area'] / sref
    return {
        'cd_friction': cdf,
        'cf': cf,
        'reynolds_used': reynolds_used,
        'roughness_cutoff_reynolds': cutoff,
        'mach_lookup': mach_lookup,
        'method': 'legacy_supdrg_straight',
    }


def calculate_supdrg_straight_wing(state: Dict, mach: float,
                                   cl_wing: float) -> Dict[str, float]:
    """Translate SUPDRG's straight-tapered wing wave-drag branches.

    The input lift coefficient is based on exposed wing area. Returned drag
    coefficients use the aircraft SREF basis, matching SUPLNG label 1540.
    """
    if mach <= 1.0:
        raise ValueError("SUPDRG requires Mach > 1")
    if float(state.get('wing_type', 1.0) or 1.0) != 1.0:
        raise ValueError("this SUPDRG translation supports straight tapered wings")

    geometry = calculate_straight_exposed_geometry(state)
    root = geometry['root_chord']
    tip = geometry['tip_chord']
    semispan = geometry['semispan']
    wing_area = geometry['area']
    sref = float(state.get('options_sref', wing_area) or wing_area)
    if wing_area <= 0.0 or sref <= 0.0:
        raise ValueError("wing area and SREF must be positive")

    beta = math.sqrt(mach**2 - 1.0)
    tan_le = geometry['tan_le']
    # WTGEOM A(18)=tip-leading-edge offset/root chord for a single panel.
    sigma = semispan * tan_le / root
    rlw = tip + sigma * root
    if rlw <= 0.0:
        raise ValueError("SUPDRG RLW must be positive")
    aspect_ratio = geometry['aspect_ratio']
    p = wing_area / (rlw * 2.0 * semispan)
    if aspect_ratio <= 0.0 or p <= 0.0:
        raise ValueError("SUPDRG aspect ratio and P must be positive")

    ksharp = state.get('wing_ksharp')
    sharp = (ksharp is not None and np.isfinite(ksharp) and
             abs(float(ksharp)) != UNUSED)
    table = _FIG_415258_SHARP if sharp else _FIG_415258_ROUND
    figure_argument = beta * semispan / rlw
    dragc = tlin1x(_FIG_415258_X, table, figure_argument, lower=0, upper=1)
    induced_factor = dragc * (1.0 + p) / (np.pi * aspect_ratio * p)
    area_ratio = wing_area / sref
    cd_wave_lift = induced_factor * area_ratio * cl_wing**2

    # SUPDRG labels 1050-1090: zero-lift wave drag.
    tan_for_branch = tan_le if tan_le != 0.0 else 1.0e-5
    sonic_leading_edge = beta / tan_for_branch >= 1.0
    wave_drag_denominator = beta if sonic_leading_edge else tan_for_branch
    tceff = float(state.get('wing_tceff', state.get('wing_tovc', 0.0)) or 0.0)
    if sharp:
        volume_numerator = float(ksharp) * tceff ** 2 * area_ratio
        cd_wave_volume = volume_numerator / wave_drag_denominator
    else:
        cos_le = 1.0 / math.sqrt(1.0 + tan_le ** 2)
        leri = float(state.get('wing_leri', 0.0) or 0.0)
        average_chord = (root + tip) / 2.0
        lerbw = leri * average_chord
        cd_le = (1.28 * mach ** 3 * cos_le ** 6 /
                 (1.0 + mach ** 3 * cos_le ** 3) *
                 (2.0 * lerbw * (2.0 * semispan) / (sref * cos_le)))
        volume_numerator = 16.0 * tceff ** 2 * area_ratio / 3.0
        cd_wave_volume = cd_le + volume_numerator / wave_drag_denominator

    return {
        'cd_wave_volume': cd_wave_volume,
        'cd_wave_lift': cd_wave_lift,
        'cd_wave_total': cd_wave_volume + cd_wave_lift,
        'dragc': dragc,
        'figure_argument': figure_argument,
        'p': p,
        'beta': beta,
        'method': 'legacy_supdrg_straight',
    }


def calculate_supersonic_lift_slope(mach: float, aspect_ratio: float,
                                    sweep_deg: float = 0.0) -> float:
    """
    Calculate lift curve slope for supersonic flow.
    
    Uses linearized supersonic theory.
    
    Args:
        mach: Mach number (> 1.0)
        aspect_ratio: Wing aspect ratio
        sweep_deg: Leading edge sweep angle (degrees)
        
    Returns:
        Lift curve slope (per radian)
    """
    if mach <= 1.0:
        logger.warning(f"Supersonic method called with M={mach} < 1.0")
        return 2.0 * np.pi
    
    beta = math.sqrt(mach ** 2 - 1.0)
    cla_2d = 4.0 / beta

    if aspect_ratio > 0:
        finite_span_factor = aspect_ratio / (aspect_ratio + 2.0 / beta)
        cla_3d = cla_2d * finite_span_factor
    else:
        finite_span_factor = 1.0
        cla_3d = cla_2d

    if abs(sweep_deg) > 1.0:
        sweep_rad = np.deg2rad(abs(sweep_deg))
        mach_normal = mach * math.cos(sweep_rad)
        if mach_normal > 1.0:
            beta_normal = math.sqrt(mach_normal ** 2 - 1.0)
            cla_3d = 4.0 / beta_normal * finite_span_factor

    return cla_3d


def calculate_supersonic_wave_drag(mach: float, thickness_ratio: float,
                                   aspect_ratio: float, lift_coef: float) -> Dict[str, float]:
    """
    Calculate wave drag for supersonic flow.
    
    Separates volume (thickness) and lift-dependent wave drag.
    
    Args:
        mach: Mach number
        thickness_ratio: Airfoil thickness ratio (t/c)
        aspect_ratio: Wing aspect ratio
        lift_coef: Lift coefficient
        
    Returns:
        Dictionary with wave drag components
    """
    if mach <= 1.0:
        return {'cd_wave_volume': 0.0, 'cd_wave_lift': 0.0, 'cd_wave_total': 0.0}
    
    beta = math.sqrt(mach**2 - 1.0)
    
    # Volume wave drag (thickness effect)
    # Linearized theory: CD_volume ≈ k * (t/c)² / β
    k_volume = 3.5  # Shape factor (depends on airfoil profile)
    cd_wave_volume = k_volume * thickness_ratio**2 / beta
    
    # Lift-dependent wave drag.  The beta/4 term is the linearized
    # two-dimensional limit; the finite-span term tends to zero as AR grows.
    if aspect_ratio > 0:
        cd_wave_lift = lift_coef**2 * (
            beta / 4.0 + 1.0 / (np.pi * aspect_ratio)
        )
    else:
        cd_wave_lift = 0.0
    
    cd_wave_total = cd_wave_volume + cd_wave_lift
    
    return {
        'cd_wave_volume': cd_wave_volume,
        'cd_wave_lift': cd_wave_lift,
        'cd_wave_total': cd_wave_total,
        'beta': beta,
    }


def calculate_supersonic_coefficients(state: Dict, alpha_deg: float,
                                      mach: float, reynolds: float) -> Dict[str, float]:
    """
    Calculate all aerodynamic coefficients for supersonic flow.
    
    Args:
        state: State dictionary
        alpha_deg: Angle of attack (degrees)
        mach: Mach number (> 1.2)
        reynolds: Reynolds number
        
    Returns:
        Dictionary with CL, CD, Cm
    """
    # Get geometry
    aspect_ratio = state.get('wing_aspect_ratio', 4.0)
    taper_ratio = state.get('wing_taper_ratio', 0.5)
    thickness_ratio = state.get('wing_tovc', 0.10)  # Thinner for supersonic
    sweep_deg = state.get('wing_savsi', 0.0)
    
    # Calculate lift curve slope
    cla = calculate_supersonic_lift_slope(mach, aspect_ratio, sweep_deg)
    
    # SETUP1 adds root incidence to the flight alpha schedule.  ALPHAI is
    # instead the section angle corresponding to the design lift CLI.
    section = resolve_wing_lift_inputs(state, mach)
    alpha_zero = section['alpha_zero']
    
    # Calculate lift
    alpha_eff_rad = np.deg2rad(alpha_deg - alpha_zero)
    cl_wing = cla * alpha_eff_rad
    reference_ratio = wing_reference_ratio(state)
    cl = cl_wing * reference_ratio
    
    # Calculate wave drag
    straight_wing = float(state.get('wing_type', 1.0) or 1.0) == 1.0
    has_supdrg_geometry = all(
        state.get(key) is not None
        for key in ('wing_chrdr', 'wing_chrdtp', 'wing_sspn')
    )
    use_supdrg = straight_wing and has_supdrg_geometry
    if use_supdrg:
        wave_drag = calculate_supdrg_straight_wing(state, mach, cl_wing)
    else:
        wave_drag = calculate_supersonic_wave_drag(
            mach, thickness_ratio, aspect_ratio, cl_wing,
        )
        for key in ('cd_wave_volume', 'cd_wave_lift', 'cd_wave_total'):
            wave_drag[key] *= reference_ratio
        wave_drag['method'] = 'linearized_fallback'

    from pydatcom.aerodynamics.drag import calculate_skin_friction_drag
    if use_supdrg:
        friction = calculate_supdrg_skin_friction(state, mach, reynolds)
        cd_friction = friction['cd_friction']
        friction_method = friction['method']
    else:
        cd_friction = calculate_skin_friction_drag(state, mach, reynolds)
        friction_method = 'generic_fallback'
    
    # Total drag
    cd = cd_friction + wave_drag['cd_wave_total']
    
    # Moment (simplified for supersonic)
    xcg = float(state.get('synths_xcg', 0.0) or 0.0)
    cbar = float(state.get('options_cbarr', 1.0) or 1.0)
    wing_chord = float(state.get('wing_mac', cbar) or cbar)
    mac_le = (float(state.get('synths_xw', 0.0) or 0.0) +
              float(state.get('wing_mac_location', 0.0) or 0.0))
    xac = float(state.get('wing_xac_abs', mac_le + 0.5 * wing_chord))
    
    if cbar > 0:
        cm = cl * (xcg - xac) / cbar
    else:
        cm = 0.0
    
    return {
        'cl': cl,
        'cl_wing': cl_wing,
        'cd': cd,
        'cm': cm,
        'cla': cla * reference_ratio * np.deg2rad(1.0),
        'cla_wing': cla * np.deg2rad(1.0),
        'wing_reference_ratio': reference_ratio,
        'xac': xac,
        'alpha_zero': alpha_zero,
        'section_alpha_zero': section['section_alpha_zero'],
        'incidence': section['incidence'],
        'cd_friction': cd_friction,
        'friction_method': friction_method,
        'cd_wave': wave_drag['cd_wave_total'],
        'cd_wave_volume': wave_drag['cd_wave_volume'],
        'cd_wave_lift': wave_drag['cd_wave_lift'],
        'wave_drag_method': wave_drag['method'],
        'regime': 'supersonic',
        'mach': mach,
        'alpha': alpha_deg,
        'beta': wave_drag['beta'],
    }

