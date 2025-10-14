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

logger = logging.getLogger(__name__)


def calculate_wing_moment_coefficient(cl: float, xac: float, xcg: float,
                                      mac: float, cbar: float) -> float:
    """
    Calculate wing pitching moment about CG.
    
    Cm = Cm0 + CL * (xac - xcg) / cbar
    
    Args:
        cl: Lift coefficient
        xac: Aerodynamic center location (from nose)
        xcg: Center of gravity location (from nose)
        mac: Mean aerodynamic chord
        cbar: Reference chord
        
    Returns:
        Pitching moment coefficient about CG
    """
    # Moment arm
    if cbar > 0:
        moment_arm = (xac - xcg) / cbar
    else:
        moment_arm = 0.0
    
    # Zero-lift pitching moment (wing camber effect)
    cm0 = 0.0  # Will be added from state or airfoil data
    
    # Total moment
    cm = cm0 - cl * moment_arm
    
    return cm


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
        return 0.0
    
    # Tail volume coefficient
    if area_tail > 0 and x_tail > xcg:
        moment_arm = (x_tail - xcg) / cbar
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
    xac_wing = state.get('wing_xac', 0.25)  # Typical AC at 0.25c
    xcg = state.get('synths_xcg', 0.0) or 0.0
    xw = state.get('synths_xw', 0.0) or 0.0  # Wing location
    cbar = state.get('options_cbarr', 1.0) or 1.0
    
    # Wing AC location from nose
    xac_abs = xw + xac_wing * cbar
    
    cm_wing = calculate_wing_moment_coefficient(
        cl_wing, xac_abs, xcg, cbar, cbar
    )
    
    # Body contribution
    cm_body = calculate_body_pitching_moment(state, alpha_deg, mach)
    
    # Tail contribution (simplified - assumes downwash effects included)
    # Full implementation would calculate tail lift from downwash
    cm_tail = 0.0  # Placeholder - Phase 5 will add tail effects
    
    # Total moment
    cm_total = cm_wing + cm_body + cm_tail
    
    return {
        'cm_total': cm_total,
        'cm_wing': cm_wing,
        'cm_body': cm_body,
        'cm_tail': cm_tail,
        'xcg': xcg,
        'xac': xac_abs,
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
        xnp = self.state.get('wing_xnp', xcg + 0.1)  # Default slightly aft
        cbar = self.state.get('options_cbarr', 1.0)
        
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

