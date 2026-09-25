"""
Body-alone aerodynamic calculations for PyDATCOM.

Implements axisymmetric and asymmetric body aerodynamics without wing/tail.
Uses slender body theory and DATCOM empirical correlations.

Reference: datcom.f lines 2326-2563 (BODYRT), 2248-2325 (BODYJM)
"""

import math
import numpy as np
from typing import Dict
import logging

from pydatcom.utils.table_lookup import fig26
from pydatcom.utils.legacy_numeric import tbfunx

logger = logging.getLogger(__name__)


def calculate_body_alone_subsonic(state: Dict, alpha_deg: float, 
                                  mach: float, reynolds: float = None) -> Dict[str, float]:
    """
    Calculate aerodynamic coefficients for body-alone configuration (subsonic).
    
    Uses slender body theory with empirical corrections.
    Reference: FORTRAN BODYRT subroutine, datcom.f line 2326
    
    Args:
        state: State dictionary with body geometry
        alpha_deg: Angle of attack (degrees)
        mach: Mach number
        reynolds: Reynolds number
        
    Returns:
        Dictionary with CL, CD, CM, CN, CA
    """
    # Get body geometry
    nx = int(state.get('body_nx', 0))
    if nx < 2:
        logger.error("Insufficient body geometry for body-alone calculation")
        return {'cl': 0.0, 'cd': 0.02, 'cm': 0.0, 'cn': 0.0, 'ca': 0.02}
    
    x = np.array(state.get('body_x', []))
    s = np.array(state.get('body_s', []))  # Cross-sectional area
    r = np.array(state.get('body_r', []))  # Half-width
    
    # Reference values
    sref = state.get('options_sref', 1.0) or 1.0
    cbar = state.get('options_cbarr', 1.0) or 1.0
    xcg = state.get('synths_xcg', x[-1]/2.0) or 0.0
    rougfc = state.get('options_rougfc', 1.6e-4) or 1.6e-4
    
    # Body length and max area
    length = x[-1]
    max_area = np.max(s)
    max_area_idx = np.argmax(s)
    
    # Base area (tail)
    base_area = s[-1]
    
    # Adjust base area if boat-tailed (line 2416)
    if base_area <= 0.3 * max_area:
        base_area = 0.3 * max_area
    
    # Equivalent diameters
    d_max = math.sqrt(4.0 * max_area / np.pi)  # Maximum diameter
    d_base = math.sqrt(4.0 * base_area / np.pi)  # Base diameter
    
    # Fineness ratio
    fineness = length / d_max if d_max > 0 else 0.0
    
    # Calculate body volume and centroid using integration
    try:
        volume = np.trapezoid(s, x)
        r_integral = np.trapezoid(r, x)
        rx_integral = np.trapezoid(r * x, x)
    except AttributeError:
        volume = np.trapz(s, x)
        r_integral = np.trapz(r, x)
        rx_integral = np.trapz(r * x, x)
    
    # Normal force coefficient slope (per radian)
    # From DATCOM Figure 4.2.1.1-20 and slender body theory
    # CNA ≈ 2 * Smax / Sref for slender bodies
    # Reference: line 2462, BODY(101) = 2 * BD(9) * TMP1 / (RAD * SREF)
    
    # Empirical factor from Figure 4.2.1.2-35A (depends on fineness ratio)
    # For fineness 4-20, factor is approximately 0.6-0.8
    if 4.0 <= fineness <= 28.0:
        # Approximate curve fit from Y1217A data (line 2383)
        k_fineness = 0.56 + 0.01 * (fineness - 4.0) / 24.0
        k_fineness = np.clip(k_fineness, 0.56, 0.79)
    else:
        k_fineness = 0.7  # Default
    
    # Mach effect on normal force (Figure 4.2.1.2-35B, line 2375-2377)
    # For subsonic: k_mach ≈ 1.2
    mach_abs_sin = mach * abs(math.sin(np.deg2rad(alpha_deg)))
    if mach_abs_sin <= 1.0:
        # Interpolate from X1217B, Y1217B data
        k_mach_data = np.array([0., 0.2, 0.3, 0.36, 0.4, 0.5, 0.6, 0.7, 0.77, 0.8, 0.86, 0.9, 0.98, 1.0])
        k_mach_vals = np.array([1.2, 1.2, 1.21, 1.23, 1.27, 1.36, 1.5, 1.67, 1.75, 1.77, 1.8, 1.8, 1.8, 1.79])
        k_mach = np.interp(mach_abs_sin, k_mach_data, k_mach_vals)
    else:
        k_mach = 1.2
    
    # Normal force per radian (line 2462)
    cna_per_rad = 2.0 * k_fineness * max_area / sref
    
    # Normal force at this alpha (line 2547)
    alpha_rad = np.deg2rad(alpha_deg)
    sin_alpha = math.sin(alpha_rad)
    sin_alpha_squared = sin_alpha ** 2

    # Cross-flow drag term (line 2547)
    # BD(J+194) = 2 * sin²(α) * BD(76) * BD(J+134) * BD(88) / SREF * SGN
    cn_crossflow = (
        2.0 * sin_alpha_squared * k_fineness * k_mach * r_integral / sref
    )
    if alpha_deg < 0:
        cn_crossflow *= -1.0

    cn_potential = cna_per_rad * alpha_rad
    cn_total = cn_potential + cn_crossflow
    
    # Skin friction drag (lines 2508-2521)
    # Use FIG26 for turbulent flat plate CF
    reynolds_length = reynolds if reynolds else 1e6 * length
    
    # Roughness Reynolds number (line 2398, 2504)
    re_roughness = 12.0 * length / rougfc
    re_use = min(reynolds_length, re_roughness)
    
    # Get skin friction coefficient
    mach_for_cf = mach if not state.get('flags_transn', False) else 0.6
    cf = fig26(re_use, mach_for_cf)
    
    # Calculate wetted area from perimeter integration (line 2467)
    p = np.array(state.get('body_p', []))
    if len(p) > 0:
        try:
            perimeter_integral = np.trapezoid(p, x)
        except AttributeError:
            perimeter_integral = np.trapz(p, x)
    else:
        # Estimate from radius
        perimeter_integral = 2.0 * np.pi * r_integral
    
    # Zero-lift drag calculation (lines 2517-2521)
    # BD(59) = CF * (1 + 60/F³ + 0.0025*F) * P_integral / Smax * Smax/Sref
    # where F is fineness ratio
    form_factor = (1.0 + 60.0 / (fineness**3) + 0.0025 * fineness) if fineness > 1.0 else 1.2
    
    # Friction drag normalized by max area then by ref area
    cd_friction_norm = cf * form_factor * perimeter_integral / max_area if max_area > 0 else cf
    cd_friction = cd_friction_norm * max_area / sref if sref > 0 else cd_friction_norm
    
    # Base drag (line 2518, 2520)
    # BD(60) = 0.029 * (d_base/d_max)³ / sqrt(CD_friction) * Smax/Sref
    if d_max > 0.01 and cd_friction_norm > 0.001:
        base_drag_factor = 0.029 * ((d_base / d_max)**3) / math.sqrt(cd_friction_norm)
        cd_base = base_drag_factor * max_area / sref
    else:
        cd_base = 0.0
    
    # Total zero-lift drag (line 2521)
    cd0 = cd_friction + cd_base
    
    # Drag due to normal force (line 2557)
    cd_normal = cn_total * sin_alpha
    cd = cd0 + cd_normal
    
    # Lift (line 2559: CL = CN*cos(α) + CD*sin(α) but that's CN, we want CL)
    # Actually line 2559: BODY(J+60) = BODY(J+20)*COS + BODY(J)*SINA
    # Where BODY(J+20) is CL, BODY(J) is CD
    # So: CN = CL*cos(α) + CD*sin(α), thus CL = (CN - CD*sin(α))/cos(α)
    # But simpler: CL = CN*cos(α) - CA*sin(α) where CA is axial force
    
    # BODYRT stores CN first, adds its lift-dependent drag, and then applies
    # the exact orthogonal CN/CD -> CL/CA transform at labels 1100.
    cl = (cn_total - cd * sin_alpha) / math.cos(alpha_rad)
    ca = cd * math.cos(alpha_rad) - cl * sin_alpha
    
    # Pitching moment (lines 2484-2488, 2552-2553)
    # Calculate centroid moment arm effect
    dsdx = np.array([tbfunx(x, s, station)[1] for station in x])
    try:
        moment_integral = np.trapezoid(dsdx * x, x)
    except AttributeError:
        moment_integral = np.trapz(dsdx * x, x)
    
    # CMA (moment curve slope) - line 2488
    const = 2.0 * k_fineness / (sref * cbar) if cbar > 0 else 0.0
    cma_per_rad = ((xcg / cbar) * cna_per_rad - const * moment_integral
                   if cbar > 0 else 0.0)
    
    # Pitching moment at this alpha (line 2552-2553)
    sign = -1.0 if alpha_deg < 0.0 else 1.0
    cm = (
        cma_per_rad * alpha_rad
        - 2.0 * sin_alpha_squared * k_mach * k_fineness
        * (rx_integral - xcg * r_integral) / (cbar * sref) * sign
    )
    
    # CLA per degree (line 2462 converted to per degree)
    cla_per_deg = cna_per_rad * np.deg2rad(1.0)
    
    return {
        'cl': cl,
        'cd': cd,
        'cm': cm,
        'cn': cn_total,
        'ca': ca,
        'cla_per_deg': cla_per_deg,
        'cd0': cd0,
        'cd_friction': cd_friction,
        'cd_base': cd_base,
        'cd_normal': cd_normal,
        'regime': 'body_alone_subsonic',
        'mach': mach,
        'alpha': alpha_deg,
    }


def has_wing_or_tail(state: Dict) -> bool:
    """
    Check if configuration has wing or tail surfaces.
    
    Args:
        state: State dictionary
        
    Returns:
        True if wing/tail present, False if body-only
    """
    # Check for wing namelists/parameters (from WGPLNF parsing)
    has_wing = (state.get('wing_chrdr') is not None or 
                state.get('wing_sspn') is not None or
                state.get('wing_chrdtp') is not None or
                state.get('wing_aspect_ratio') is not None or
                state.get('wing_area') is not None or
                state.get('wing_span') is not None)
    
    # Check for tail namelists (from HTPLNF, VTPLNF parsing)
    has_htail = (state.get('htail_chrdr') is not None or
                 state.get('htail_sspn') is not None or
                 state.get('htail_area') is not None)
    
    has_vtail = (state.get('vtail_chrdr') is not None or
                 state.get('vtail_sspn') is not None or
                 state.get('vtail_area') is not None)
    
    return has_wing or has_htail or has_vtail


def calculate_body_alone_coefficients(state: Dict, alpha_deg: float,
                                      mach: float, reynolds: float = None) -> Dict[str, float]:
    """
    Calculate aerodynamic coefficients for body-alone configuration.
    
    Automatically routes to appropriate method based on Mach number.
    
    Args:
        state: State dictionary
        alpha_deg: Angle of attack (degrees)
        mach: Mach number
        reynolds: Reynolds number (estimated if None)
        
    Returns:
        Dictionary with aerodynamic coefficients
    """
    # Estimate Reynolds if not provided
    if reynolds is None:
        body_length = state.get('body_length')
        if body_length is None and state.get('body_x'):
            body_length = state['body_x'][-1]
        if body_length is None:
            body_length = 10.0
        reynolds = 1e6 * body_length * mach
    
    if mach < 0.9:
        return calculate_body_alone_subsonic(state, alpha_deg, mach, reynolds)
    if mach < 1.2:
        result = calculate_body_alone_subsonic(state, alpha_deg, 0.85, reynolds)
        result['regime'] = 'body_alone_transonic'
        result['mach'] = mach
        return result
    if mach < 5.0:
        return calculate_body_alone_supersonic(state, alpha_deg, mach)
    return calculate_body_alone_hypersonic(state, alpha_deg, mach)


def calculate_body_alone_supersonic(state: Dict, alpha_deg: float, mach: float) -> Dict[str, float]:
    """
    Calculate body-alone coefficients for supersonic flow.
    
    Args:
        state: State dictionary
        alpha_deg: Angle of attack
        mach: Mach number (> 1.2)
        
    Returns:
        Aerodynamic coefficients
    """
    # Get geometry
    max_area = state.get('body_max_area', 0.0) or np.max(state.get('body_s', [0]))
    sref = state.get('options_sref', 1.0) or 1.0
    
    alpha_rad = np.deg2rad(alpha_deg)
    sin_alpha = np.sin(alpha_rad)
    cos_alpha = np.cos(alpha_rad)

    # Normal force (supersonic slender body; reduced vs subsonic shock effects)
    cn = 1.5 * (max_area / sref) * np.sin(2.0 * alpha_rad)

    length = state.get('body_length', 10.0) or 10.0
    max_diameter = np.sqrt(4.0 * max_area / np.pi)
    fineness = length / max_diameter if max_diameter > 0 else 5.0

    cd_wave = 0.15 / fineness**2 if fineness > 0 else 0.02
    cd_friction = 0.01  # Simplified
    cd_normal = cn * sin_alpha
    cd = cd_friction + cd_wave + cd_normal

    ca = cd * cos_alpha - cn * sin_alpha
    cl = cn * cos_alpha - ca * sin_alpha

    xcg = state.get('synths_xcg', length / 2.0) or 0.0
    cbar = state.get('options_cbarr', 1.0) or 1.0
    cm = -cn * (length/2.0 - xcg) / cbar if cbar > 0 else 0.0
    
    return {
        'cl': cl,
        'cd': cd,
        'cm': cm,
        'cn': cn,
        'ca': ca,
        'cla_per_deg': 1.5 * (max_area / sref) * 2.0 * np.cos(2.0 * alpha_rad) * np.deg2rad(1.0),
        'regime': 'body_alone_supersonic',
        'mach': mach,
        'alpha': alpha_deg,
    }


def calculate_body_alone_hypersonic(state: Dict, alpha_deg: float, mach: float) -> Dict[str, float]:
    """
    Calculate body-alone coefficients for hypersonic flow.
    
    Uses Newtonian impact theory.
    
    Args:
        state: State dictionary
        alpha_deg: Angle of attack
        mach: Mach number (> 5.0)
        
    Returns:
        Aerodynamic coefficients
    """
    # Use hypersonic module but with body-specific parameters
    from pydatcom.aerodynamics.hypersonic import calculate_hypersonic_coefficients
    
    result = calculate_hypersonic_coefficients(state, alpha_deg, mach)
    result['regime'] = 'body_alone_hypersonic'
    
    return result

