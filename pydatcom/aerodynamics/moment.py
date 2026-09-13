"""
Pitching moment coefficient calculations for PyDATCOM.

Implements moment calculations:
- Wing pitching moment
- Body contribution
- Tail contribution
- Center of gravity effects

Reference: datcom.f lines 4551 (CMALPH), 5105 (CMALPO), 5386 (CNCA)
"""

import numpy as np
from typing import Dict
import logging

from pydatcom.utils.legacy_numeric import tbfunx
from pydatcom.geometry.wing import calculate_straight_exposed_geometry

logger = logging.getLogger(__name__)

_CMALPH_MACH = np.array([0.0, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50,
                         0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.90])
_CMALPH_CALM = np.array([1.000, 1.000, 1.005, 1.017, 1.031, 1.050,
                         1.072, 1.101, 1.132, 1.162, 1.197, 1.237,
                         1.287, 1.355, 1.445])


def calculate_cmalph_zero_lift_moment(state: Dict, mach: float) -> float:
    """Translate CMALPH's constant-section, untwisted CMO calculation."""
    geometry = calculate_straight_exposed_geometry(state)
    area = geometry['area']
    mac = geometry['mac']
    sref = float(state.get('options_sref', area) or area)
    cbar = float(state.get('options_cbarr', mac) or mac)
    ar = geometry['aspect_ratio']
    cmo_root = float(state.get('wing_cmo', 0.0) or 0.0)
    cmo_tip = float(state.get('wing_cmot', 0.0) or 0.0)
    # CMALPH labels 269-281 use the mean only when both section values are
    # nonzero; otherwise the surface is treated as constant-section.
    cmo = (0.5 * (cmo_root + cmo_tip)
           if abs(cmo_root) >= 1.0e-10 and abs(cmo_tip) >= 1.0e-10
           else cmo_root)
    if min(area, mac, sref, cbar, ar) <= 0.0:
        raise ValueError("CMALPH CMO requires complete positive straight-wing geometry")
    if abs(float(state.get('wing_twista', 0.0) or 0.0)) >= 1.0e-20:
        raise ValueError("twisted-wing CMALPH CMO correction is not translated")

    tan_c4 = geometry['tan_c4']
    cos_c4 = 1.0 / np.sqrt(1.0 + tan_c4**2)
    calm, _ = tbfunx(_CMALPH_MACH, _CMALPH_CALM, mach, lower=0, upper=0)
    calm *= area * mac / (sref * cbar)
    return ar * cos_c4**2 / (ar + 2.0 * cos_c4) * cmo * calm


def resolve_wing_xac(state: Dict, fraction: float = 0.25) -> float:
    """Resolve a wing aerodynamic-center fraction to an absolute X."""
    absolute = state.get('wing_xac_abs')
    if absolute is not None:
        return float(absolute)
    xw = float(state.get('synths_xw', 0.0) or 0.0)
    # CMALPH's C(6) is XAC/root chord from the root leading edge.
    root_fraction = state.get('wing_xac_root_fraction')
    if root_fraction is not None:
        geometry = calculate_straight_exposed_geometry(state)
        return (geometry['exposed_root_x'] +
                float(root_fraction) * geometry['root_chord'])
    if 'wing_xac' in state:
        # Preserve the original Python API's documented CBARR fraction.
        cbar = float(state.get('options_cbarr', 0.0) or 0.0)
        return xw + float(state['wing_xac']) * cbar
    try:
        geometry = calculate_straight_exposed_geometry(state)
        return (geometry['exposed_root_x'] + geometry['mac_le_location'] +
                fraction * geometry['mac'])
    except ValueError:
        wing_mac = float(state.get('wing_mac', state.get('options_cbarr', 1.0)) or 1.0)
        mac_le = float(state.get('wing_mac_location', 0.0) or 0.0)
        return xw + mac_le + fraction * wing_mac


def calculate_wing_moment_coefficient(cl: float, xac: float, xcg: float,
                                      mac: float, cbar: float,
                                      cm_ac: float = 0.0) -> float:
    """
    Calculate wing pitching moment about CG.
    
    Cm_cg = Cm_ac + CL * (xcg - xac) / cbar

    This is the sign and coordinate convention used by ``CMALPH``:
    ``DCMDCL=(DXCG-XAC)/CBARR`` after both locations have been put in
    the same dimensional coordinate system.
    
    Args:
        cl: Lift coefficient
        xac: Aerodynamic center location (from nose)
        xcg: Center of gravity location (from nose)
        mac: Mean aerodynamic chord (retained for API compatibility)
        cbar: Reference chord
        cm_ac: Pitching moment coefficient about the aerodynamic center
        
    Returns:
        Pitching moment coefficient about CG
    """
    del mac
    if cbar <= 0:
        raise ValueError("reference chord must be positive")
    return cm_ac + cl * (xcg - xac) / cbar


def calculate_body_pitching_moment(state: Dict, alpha_deg: float,
                                   mach: float) -> float:
    """
    Calculate body contribution to pitching moment.
    
    Based on body geometry and flow conditions.
    
    Args:
        state: State dictionary with body geometry
        alpha_deg: Angle of attack (degrees)
        mach: Mach number
        
    Returns:
        Body pitching moment coefficient
    """
    # Get body geometry
    length = state.get('body_length', 0.0)
    max_area = state.get('body_max_area', 0.0)
    centroid = state.get('body_centroid', length / 2.0)
    
    # Reference values
    sref = state.get('options_sref', 1.0) or 1.0
    cbar = state.get('options_cbarr', 1.0) or 1.0
    xcg = state.get('synths_xcg', length / 2.0 if length else 0.0) or 0.0
    
    if sref <= 0 or cbar <= 0:
        return 0.0
    
    # Body normal force coefficient (simplified)
    # CN_body ≈ k * (S_max / S_ref) * sin(2α)
    alpha_rad = np.deg2rad(alpha_deg)
    k_body = 1.0  # Body shape factor
    
    if max_area > 0 and sref > 0:
        cn_body = k_body * (max_area / sref) * np.sin(2.0 * alpha_rad)
    else:
        cn_body = 0.0
    
    # Moment arm from CG to body centroid
    moment_arm = (centroid - xcg) / cbar
    
    # Body moment
    cm_body = -cn_body * moment_arm
    
    return cm_body


def calculate_tail_moment_contribution(state: Dict, cl_tail: float,
                                       tail_type: str = 'htail') -> float:
    """
    Calculate tail contribution to pitching moment.
    
    Args:
        state: State dictionary
        cl_tail: Tail lift coefficient
        tail_type: 'htail' or 'vtail'
        
    Returns:
        Tail pitching moment coefficient
    """
    # Get tail geometry
    area_tail = state.get(f'{tail_type}_area', 0.0)
    x_tail = state.get(f'synths_x{"h" if tail_type == "htail" else "v"}', 0.0)
    
    # Reference values
    sref = state.get('options_sref', 1.0)
    cbar = state.get('options_cbarr', 1.0)
    xcg = state.get('synths_xcg', 0.0)
    
    if sref <= 0 or cbar <= 0:
        raise ValueError("reference area and chord must be positive")
    
    # Tail volume coefficient
    if area_tail > 0 and x_tail > xcg:
        # XH/XV is the surface longitudinal reference location in /SYNTSS/.
        # A translated aerodynamic-center offset may be supplied in length
        # units; omitting it means the force is referenced at XH/XV.
        xac_tail = x_tail + state.get(f'{tail_type}_xac_offset', 0.0)
        moment_arm = (xac_tail - xcg) / cbar
        volume_coef = (area_tail / sref) * moment_arm
        
        # Tail moment
        cm_tail = -cl_tail * volume_coef
    else:
        cm_tail = 0.0
    
    return cm_tail


def calculate_total_pitching_moment(state: Dict, cl_wing: float, 
                                    alpha_deg: float, mach: float) -> Dict[str, float]:
    """
    Calculate total aircraft pitching moment about CG.
    
    Cm_total = Cm_wing + Cm_body + Cm_tail
    
    Args:
        state: State dictionary
        cl_wing: Wing lift coefficient
        alpha_deg: Angle of attack (degrees)
        mach: Mach number
        
    Returns:
        Dictionary with moment components
    """
    # Wing contribution
    xcg = state.get('synths_xcg', 0.0) or 0.0
    cbar = state.get('options_cbarr', 1.0) or 1.0
    
    # Callers with a translated dimensional result can provide the
    # unambiguous absolute coordinate directly.  The older public key is
    # retained as a reference-chord fraction for compatibility.
    xac_abs = resolve_wing_xac(state)

    cm_method = 'supplied_aircraft_cm_ac'
    cm_ac = state.get('wing_cm_ac')
    if cm_ac is None:
        # WGPLNF CMO is an airfoil/wing-reference coefficient.  Put it on
        # the aircraft SREF/CBARR basis before adding it to aircraft Cm.
        wing_cmo = state.get('wing_cmo', 0.0) or 0.0
        wing_area = state.get('wing_area', state.get('options_sref', 1.0)) or 0.0
        wing_mac = state.get('wing_mac', cbar) or cbar
        sref = state.get('options_sref', wing_area) or wing_area
        complete_cmalph = all(state.get(key) is not None for key in (
            'wing_chrdr', 'wing_chrdtp', 'wing_sspn')) and float(
                state.get('wing_type', 1.0) or 1.0) == 1.0
        if complete_cmalph and abs(float(state.get('wing_twista', 0.0) or 0.0)) < 1.0e-20:
            cm_ac = calculate_cmalph_zero_lift_moment(state, mach)
            cm_method = 'legacy_cmalph_constant_section'
        else:
            cm_ac = (wing_cmo * wing_area * wing_mac / (sref * cbar)
                     if sref > 0.0 else 0.0)
            cm_method = 'reference_scaled_fallback'
    
    cm_wing = calculate_wing_moment_coefficient(
        cl_wing, xac_abs, xcg, cbar, cbar, cm_ac
    )
    
    # Body contribution
    cm_body = calculate_body_pitching_moment(state, alpha_deg, mach)
    
    # W B TAIL computes the horizontal-tail load after downwash and
    # interference.  Do not synthesize that load here.  Include it only when
    # an upstream translated routine has supplied a tail lift coefficient.
    cl_tail = state.get('aero_cl_tail')
    cm_tail = (calculate_tail_moment_contribution(state, cl_tail)
               if cl_tail is not None else 0.0)
    
    # Total moment
    cm_total = cm_wing + cm_body + cm_tail
    
    return {
        'cm_total': cm_total,
        'cm_wing': cm_wing,
        'cm_body': cm_body,
        'cm_tail': cm_tail,
        'xcg': xcg,
        'xac': xac_abs,
        'tail_supported': cl_tail is not None,
        'wing_cm_method': cm_method,
    }


def calculate_normal_force_coefficient(cl: float, cd: float, 
                                       alpha_deg: float) -> float:
    """
    Calculate normal force coefficient.
    
    CN = CL*cos(α) + CD*sin(α)
    
    Reference: FORTRAN CNCA subroutine, datcom.f line 5386
    
    Args:
        cl: Lift coefficient
        cd: Drag coefficient
        alpha_deg: Angle of attack (degrees)
        
    Returns:
        Normal force coefficient
    """
    alpha_rad = np.deg2rad(alpha_deg)
    
    cn = cl * np.cos(alpha_rad) + cd * np.sin(alpha_rad)
    
    return cn


def calculate_axial_force_coefficient(cl: float, cd: float,
                                      alpha_deg: float) -> float:
    """
    Calculate axial force coefficient.
    
    CA = CD*cos(α) - CL*sin(α)
    
    Args:
        cl: Lift coefficient
        cd: Drag coefficient
        alpha_deg: Angle of attack (degrees)
        
    Returns:
        Axial force coefficient
    """
    alpha_rad = np.deg2rad(alpha_deg)
    
    ca = cd * np.cos(alpha_rad) - cl * np.sin(alpha_rad)
    
    return ca


class MomentCalculator:
    """
    Comprehensive pitching moment calculations.
    
    Reference: datcom.f CMALPH (line 4551), CMALPO (line 5105)
    """
    
    def __init__(self, state: Dict):
        """
        Initialize moment calculator.
        
        Args:
            state: Global state dictionary
        """
        self.state = state
    
    def calculate_moment(self, cl: float, alpha_deg: float, 
                        mach: float) -> Dict[str, float]:
        """
        Calculate pitching moment at given conditions.
        
        Args:
            cl: Lift coefficient
            alpha_deg: Angle of attack (degrees)
            mach: Mach number
            
        Returns:
            Dictionary with moment results
        """
        return calculate_total_pitching_moment(self.state, cl, alpha_deg, mach)
    
    def calculate_moment_curve_slope(self, mach: float) -> float:
        """
        Calculate dCm/dα (moment curve slope).
        
        Args:
            mach: Mach number
            
        Returns:
            Moment curve slope (per radian)
        """
        # Get neutral point and CG locations
        xcg = self.state.get('synths_xcg', 0.0)
        cbar = self.state.get('options_cbarr', 1.0)
        xnp = self.state.get('wing_xnp_abs')
        if xnp is None:
            if 'wing_xnp' in self.state:
                xw = self.state.get('synths_xw', 0.0) or 0.0
                xnp = xw + self.state['wing_xnp'] * cbar
            else:
                xnp = resolve_wing_xac(self.state)
        
        # Get lift curve slope
        aspect_ratio = self.state.get('wing_aspect_ratio', 6.0)
        taper_ratio = self.state.get('wing_taper_ratio', 0.5)
        
        from pydatcom.aerodynamics.lift import calculate_lift_curve_slope_compressible
        cla = calculate_lift_curve_slope_compressible(
            aspect_ratio, taper_ratio, mach
        )
        
        # Cm_α = -CL_α * (xnp - xcg) / cbar
        if cbar > 0:
            cma = -cla * (xnp - xcg) / cbar
        else:
            cma = 0.0
        
        return cma

