"""
Stability derivative calculations for PyDATCOM.

Computes dynamic stability derivatives:
- Pitching moment derivatives (Cm_q, Cm_alpha_dot)
- Damping derivatives (Cn_r, Cl_p)
- Control derivatives

Reference: datcom.f DAMP-related subroutines
"""

import numpy as np
from typing import Dict
import logging

logger = logging.getLogger(__name__)


def calculate_pitch_damping(state: Dict, mach: float) -> Dict[str, float]:
    """
    Calculate pitch damping derivative Cm_q.
    
    Cm_q is the change in pitching moment due to pitch rate.
    
    Args:
        state: State dictionary
        mach: Mach number
        
    Returns:
        Dictionary with damping derivatives
    """
    # Get geometry
    aspect_ratio = state.get('wing_aspect_ratio', 6.0)
    taper_ratio = state.get('wing_taper_ratio', 0.5)
    
    # Wing contribution to Cm_q (empirical)
    # Typical values: -3 to -20 per radian
    
    # Base value from aspect ratio
    cmq_wing = -2.0 * aspect_ratio / (aspect_ratio + 2.0)
    
    # Compressibility effect
    if mach < 0.9:
        beta = np.sqrt(1.0 - mach**2)
        cmq_wing /= beta
    
    # Tail contribution (if present)
    htail_area = state.get('htail_area', 0.0)
    if htail_area and htail_area > 0:
        sref = state.get('options_sref', 1.0)
        tail_volume = htail_area / sref
        cmq_tail = -10.0 * tail_volume  # Tail provides significant damping
    else:
        cmq_tail = 0.0
    
    cmq_total = cmq_wing + cmq_tail
    
    return {
        'cmq': cmq_total,
        'cmq_wing': cmq_wing,
        'cmq_tail': cmq_tail,
    }


def calculate_roll_damping(state: Dict, mach: float) -> float:
    """
    Calculate roll damping derivative Cl_p.
    
    Cl_p is the rolling moment due to roll rate.
    
    Args:
        state: State dictionary
        mach: Mach number
        
    Returns:
        Cl_p (per radian)
    """
    # Get geometry
    aspect_ratio = state.get('wing_aspect_ratio', 6.0)
    taper_ratio = state.get('wing_taper_ratio', 0.5)
    
    # Empirical formula (simplified)
    # Cl_p ≈ -CL_α / 12 for straight wing
    
    # Approximate lift slope
    if mach < 0.9:
        beta = np.sqrt(1.0 - mach**2)
        cla = (2.0 * np.pi * aspect_ratio) / (2.0 + np.sqrt(4.0 + aspect_ratio**2)) / beta
    else:
        beta_super = np.sqrt(mach**2 - 1.0) if mach > 1.0 else 0.1
        cla = 4.0 / beta_super
    
    clp = -cla / 12.0
    
    return clp


def calculate_yaw_damping(state: Dict, mach: float) -> float:
    """
    Calculate yaw damping derivative Cn_r.
    
    Cn_r is the yawing moment due to yaw rate.
    
    Args:
        state: State dictionary
        mach: Mach number
        
    Returns:
        Cn_r (per radian)
    """
    # Vertical tail contribution dominates
    vtail_area = state.get('vtail_area', 0.0)
    sref = state.get('options_sref', 1.0)
    
    if vtail_area and vtail_area > 0 and sref > 0:
        # Tail volume coefficient
        vtail_volume = vtail_area / sref
        
        # Simplified: Cn_r ≈ -2 * V_v
        cnr = -2.0 * vtail_volume
    else:
        cnr = -0.1  # Minimal without vertical tail
    
    return cnr


def calculate_static_stability_margin(state: Dict, mach: float) -> Dict[str, float]:
    """
    Calculate static stability margin.
    
    Static margin = (xnp - xcg) / cbar
    where xnp is neutral point location.
    
    Args:
        state: State dictionary
        mach: Mach number
        
    Returns:
        Dictionary with stability parameters
    """
    xcg = state.get('synths_xcg', 0.0) or 0.0
    cbar = state.get('options_cbarr', 1.0) or 1.0
    
    # Estimate neutral point location
    # For wing alone: xnp ≈ xac_wing ≈ 0.25c from wing LE
    xw = state.get('synths_xw', 0.0) or 0.0
    xnp_wing = xw + 0.25 * cbar
    
    # Tail contribution moves NP aft
    htail_area = state.get('htail_area', 0.0)
    sref = state.get('options_sref', 1.0)
    
    if htail_area and htail_area > 0 and sref > 0:
        xh = state.get('synths_xh', xw + 2.0 * cbar) or (xw + 2.0 * cbar)
        tail_volume = (htail_area / sref) * (xh - xcg) / cbar
        
        # NP shift due to tail
        xnp = xnp_wing + tail_volume * cbar
    else:
        xnp = xnp_wing
    
    # Static margin
    if cbar > 0:
        static_margin = (xnp - xcg) / cbar
    else:
        static_margin = 0.0
    
    return {
        'xnp': xnp,
        'xcg': xcg,
        'static_margin': static_margin,
        'stable': static_margin > 0.05,  # At least 5% margin for stability
    }


def calculate_directional_stability(state: Dict, mach: float) -> float:
    """
    Calculate directional stability derivative Cn_beta.
    
    Cn_beta is the yawing moment due to sideslip.
    
    Args:
        state: State dictionary
        mach: Mach number
        
    Returns:
        Cn_beta (per radian)
    """
    # Vertical tail contribution
    vtail_area = state.get('vtail_area', 0.0)
    sref = state.get('options_sref', 1.0)
    bref = state.get('options_blref', 10.0)
    
    if vtail_area and vtail_area > 0 and sref > 0:
        # Tail volume coefficient
        xv = state.get('synths_xv', 20.0) or 20.0
        xcg = state.get('synths_xcg', 10.0) or 10.0
        
        tail_arm = xv - xcg
        volume_v = (vtail_area / sref) * (tail_arm / bref)
        
        # Simplified: Cn_beta ≈ V_v * CL_alpha_tail
        # Assume tail lift slope ≈ 3 per radian
        cnbeta = volume_v * 3.0
    else:
        cnbeta = 0.05  # Minimal without tail
    
    return cnbeta


def calculate_all_stability_derivatives(state: Dict, mach: float) -> Dict[str, float]:
    """
    Calculate complete set of stability derivatives.
    
    Args:
        state: State dictionary
        mach: Mach number
        
    Returns:
        Dictionary with all stability derivatives
    """
    # Longitudinal derivatives
    pitch_damp = calculate_pitch_damping(state, mach)
    stability = calculate_static_stability_margin(state, mach)
    
    # Lateral-directional derivatives
    clp = calculate_roll_damping(state, mach)
    cnr = calculate_yaw_damping(state, mach)
    cnbeta = calculate_directional_stability(state, mach)
    
    return {
        # Longitudinal
        'cmq': pitch_damp['cmq'],
        'cm_alpha_dot': pitch_damp['cmq'] / 2.0,  # Approximation
        'xnp': stability['xnp'],
        'static_margin': stability['static_margin'],
        'longitudinally_stable': stability['stable'],
        
        # Lateral-directional
        'clp': clp,
        'cnr': cnr,
        'cn_beta': cnbeta,
        'directionally_stable': cnbeta > 0,
    }


class StabilityCalculator:
    """
    Comprehensive stability derivative calculations.
    
    Computes dynamic and static stability characteristics.
    """
    
    def __init__(self, state: Dict):
        """
        Initialize stability calculator.
        
        Args:
            state: Global state dictionary
        """
        self.state = state
    
    def calculate_derivatives(self, mach: float) -> Dict[str, float]:
        """
        Calculate all stability derivatives at given Mach.
        
        Args:
            mach: Mach number
            
        Returns:
            Dictionary with all derivatives
        """
        return calculate_all_stability_derivatives(self.state, mach)
    
    def assess_stability(self, mach: float) -> Dict[str, any]:
        """
        Assess overall aircraft stability.
        
        Args:
            mach: Mach number
            
        Returns:
            Dictionary with stability assessment
        """
        derivs = self.calculate_derivatives(mach)
        
        assessment = {
            'longitudinal': 'stable' if derivs['longitudinally_stable'] else 'unstable',
            'directional': 'stable' if derivs['directionally_stable'] else 'unstable',
            'static_margin_percent': derivs['static_margin'] * 100,
            'derivatives': derivs,
        }
        
        # Overall assessment
        if derivs['longitudinally_stable'] and derivs['directionally_stable']:
            assessment['overall'] = 'stable'
        elif derivs['longitudinally_stable'] or derivs['directionally_stable']:
            assessment['overall'] = 'partially stable'
        else:
            assessment['overall'] = 'unstable'
        
        return assessment

