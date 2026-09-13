"""
Supersonic aerodynamic calculations for PyDATCOM.

Implements supersonic flow methods (Mach > 1.2):
- Linearized supersonic theory
- Shock-expansion method
- Prandtl-Meyer relations

Reference: datcom.f supersonic calculation sections
"""

import numpy as np
from typing import Dict
import logging

from pydatcom.utils.table_lookup import fig60b, fig68
from pydatcom.aerodynamics.lift import resolve_wing_lift_inputs, wing_reference_ratio

logger = logging.getLogger(__name__)


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
    
    beta = np.sqrt(mach**2 - 1.0)
    
    # Linearized supersonic theory
    # CL_α = 4 / β for infinite aspect ratio
    cla_2d = 4.0 / beta
    
    # Finite span correction
    if aspect_ratio > 0:
        # Simplified correction
        ar_correction = aspect_ratio / (aspect_ratio + 2.0 / beta)
        cla_3d = cla_2d * ar_correction
    else:
        cla_3d = cla_2d
    
    # Sweep correction
    if abs(sweep_deg) > 1.0:
        sweep_rad = np.deg2rad(abs(sweep_deg))
        # Component of Mach normal to leading edge
        mach_normal = mach * np.cos(sweep_rad)
        if mach_normal > 1.0:
            beta_normal = np.sqrt(mach_normal**2 - 1.0)
            cla_3d = 4.0 / beta_normal * ar_correction
    
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
    
    beta = np.sqrt(mach**2 - 1.0)
    
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
    wave_drag = calculate_supersonic_wave_drag(
        mach, thickness_ratio, aspect_ratio, cl_wing
    )
    for key in ('cd_wave_volume', 'cd_wave_lift', 'cd_wave_total'):
        wave_drag[key] *= reference_ratio
    
    # Skin friction (still present in supersonic)
    from pydatcom.aerodynamics.drag import calculate_skin_friction_drag
    cd_friction = calculate_skin_friction_drag(state, mach, reynolds)
    
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
        'cd_wave': wave_drag['cd_wave_total'],
        'cd_wave_volume': wave_drag['cd_wave_volume'],
        'cd_wave_lift': wave_drag['cd_wave_lift'],
        'regime': 'supersonic',
        'mach': mach,
        'alpha': alpha_deg,
        'beta': wave_drag['beta'],
    }

