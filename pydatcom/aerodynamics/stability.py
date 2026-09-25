"""
Stability derivative calculations for PyDATCOM.

Computes dynamic stability derivatives:
- Pitching moment derivatives (Cm_q, Cm_alpha_dot)
- Damping derivatives (Cn_r, Cl_p)
- Control derivatives

Reference: datcom.f DAMP-related subroutines
"""

import math
import numpy as np
from typing import Dict
import logging

from pydatcom.aerodynamics.moment import resolve_wing_xac

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
    # DATCOM obtains these terms from DNPWBT after the wing/body and tail
    # interference calculations.  Prefer those translated outputs whenever
    # they are available.
    translated = [state.get('wing_cmq'), state.get('body_cmq'),
                  state.get('htail_cmq')]
    if any(value is not None for value in translated):
        cmq_wing = state.get('wing_cmq', 0.0) or 0.0
        cmq_body = state.get('body_cmq', 0.0) or 0.0
        cmq_tail = state.get('htail_cmq', 0.0) or 0.0
        return {
            'cmq': cmq_wing + cmq_body + cmq_tail,
            'cmq_wing': cmq_wing,
            'cmq_body': cmq_body,
            'cmq_tail': cmq_tail,
            'datcom_supported': True,
        }

    # Compatibility estimate used until DNPWBT is fully translated.
    aspect_ratio = state.get('wing_aspect_ratio', 6.0)
    taper_ratio = state.get('wing_taper_ratio', 0.5)
    
    # Wing contribution to Cm_q (empirical)
    # Typical values: -3 to -20 per radian
    
    # Base value from aspect ratio
    cmq_wing = -2.0 * aspect_ratio / (aspect_ratio + 2.0)
    
    # Compressibility effect
    if mach < 0.9:
        beta = math.sqrt(1.0 - mach**2)
        cmq_wing /= beta
    
    # Tail contribution (if present)
    htail_area = state.get('htail_area', 0.0)
    if htail_area and htail_area > 0:
        sref = state.get('options_sref', 1.0)
        tail_area_ratio = htail_area / sref
        cmq_tail = -10.0 * tail_area_ratio
    else:
        cmq_tail = 0.0
    
    cmq_total = cmq_wing + cmq_tail
    
    return {
        'cmq': cmq_total,
        'cmq_wing': cmq_wing,
        'cmq_body': 0.0,
        'cmq_tail': cmq_tail,
        'datcom_supported': False,
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
    aspect_ratio = state.get('wing_aspect_ratio', 6.0)

    if mach < 0.9:
        beta = math.sqrt(1.0 - mach ** 2)
        lift_slope_per_rad = (
            (2.0 * np.pi * aspect_ratio) /
            (2.0 + math.sqrt(4.0 + aspect_ratio ** 2)) / beta
        )
    else:
        supersonic_beta = np.sqrt(mach ** 2 - 1.0) if mach > 1.0 else 0.1
        lift_slope_per_rad = 4.0 / supersonic_beta

    clp = -lift_slope_per_rad / 12.0
    
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
    
    # CMALPH uses CMA = CLA * (XCG-XAC)/CBARR.  If translated total
    # derivatives are present, invert that exact relation.
    cla_total = state.get('aero_cla')
    cma_total = state.get('aero_cma')
    if cla_total is not None and cma_total is not None and cla_total != 0.0:
        xnp = xcg - cbar * cma_total / cla_total
        static_margin = (xnp - xcg) / cbar
        return {
            'xnp': xnp,
            'xcg': xcg,
            'static_margin': static_margin,
            'stable': static_margin > 0.0,
            'meets_five_percent_margin': static_margin >= 0.05,
            'tail_supported': state.get('htail_cla') is not None,
            'method': 'datcom_derivatives',
        }

    # Wing-alone XAC.  ``wing_xac_abs`` is the unambiguous dimensional
    # form; the historical ``wing_xac`` public key is a CBARR fraction.
    xnp_wing = resolve_wing_xac(state)
    
    # Tail contribution moves NP aft
    htail_area = state.get('htail_area', 0.0)
    sref = state.get('options_sref', 1.0)
    
    wing_cla = state.get('wing_cla')
    tail_cla = state.get('htail_cla')
    tail_supported = bool(htail_area and htail_area > 0 and sref > 0 and
                          wing_cla is not None and wing_cla > 0 and
                          tail_cla is not None)
    if tail_supported:
        xh = state.get('synths_xh')
        if xh is None:
            raise ValueError("synths_xh is required with htail_cla")
        xac_h = xh + state.get('htail_xac_offset', 0.0)
        q_ratio = state.get('htail_dynamic_pressure_ratio', 1.0)
        downwash = state.get('htail_downwash_gradient', 0.0)
        interference = state.get('htail_lift_interference_factor', 1.0)
        effective_tail_cla = (
            tail_cla * htail_area / sref * q_ratio *
            (1.0 - downwash) * interference
        )
        total_cla = wing_cla + effective_tail_cla
        if total_cla != 0.0:
            xnp = ((wing_cla * xnp_wing + effective_tail_cla * xac_h) /
                   total_cla)
        else:
            xnp = xnp_wing
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
        'stable': static_margin > 0.0,
        'meets_five_percent_margin': static_margin >= 0.05,
        'tail_supported': tail_supported,
        'method': 'component_derivatives' if tail_supported else 'wing_xac',
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

