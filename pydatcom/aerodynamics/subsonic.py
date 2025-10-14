"""
Subsonic aerodynamic calculations for PyDATCOM.

Implements subsonic flow methods (Mach < 0.9):
- Prandtl-Glauert compressibility corrections
- Lifting line theory
- Vortex lattice approximations

Reference: datcom.f subsonic calculation sections
"""

import numpy as np
from typing import Dict
import logging

from pydatcom.aerodynamics.lift import calculate_wing_lift_subsonic
from pydatcom.aerodynamics.drag import calculate_total_drag
from pydatcom.aerodynamics.moment import calculate_total_pitching_moment

logger = logging.getLogger(__name__)


def calculate_prandtl_glauert_factor(mach: float) -> float:
    """
    Calculate Prandtl-Glauert compressibility factor.
    
    β = sqrt(1 - M²)
    
    Args:
        mach: Mach number (< 1.0)
        
    Returns:
        Beta factor
    """
    if mach >= 1.0:
        logger.warning(f"Prandtl-Glauert called with M={mach} >= 1.0")
        return 0.1  # Avoid division by zero
    
    beta = np.sqrt(1.0 - mach**2)
    return beta


def calculate_subsonic_coefficients(state: Dict, alpha_deg: float,
                                    mach: float, reynolds: float) -> Dict[str, float]:
    """
    Calculate all aerodynamic coefficients for subsonic flow.
    
    Args:
        state: State dictionary
        alpha_deg: Angle of attack (degrees)
        mach: Mach number (< 0.9)
        reynolds: Reynolds number
        
    Returns:
        Dictionary with CL, CD, Cm, and components
    """
    # Calculate lift
    lift_result = calculate_wing_lift_subsonic(state, alpha_deg, mach)
    cl = lift_result['cl']
    
    # Calculate drag
    drag_result = calculate_total_drag(state, cl, mach, reynolds)
    cd = drag_result['cd_total']
    
    # Calculate moment
    moment_result = calculate_total_pitching_moment(state, cl, alpha_deg, mach)
    cm = moment_result['cm_total']
    
    # Combine results
    return {
        'cl': cl,
        'cd': cd,
        'cm': cm,
        'cla': lift_result['cla_per_deg'],  # Per degree
        'alpha_zero': lift_result['alpha_zero'],
        'cd_friction': drag_result['cd_friction'],
        'cd_induced': drag_result['cd_induced'],
        'cd_wave': drag_result['cd_wave'],
        'cm_wing': moment_result['cm_wing'],
        'cm_body': moment_result['cm_body'],
        'regime': 'subsonic',
        'mach': mach,
        'alpha': alpha_deg,
        'reynolds': reynolds,
    }


def calculate_lift_distribution_subsonic(state: Dict, cl: float) -> np.ndarray:
    """
    Calculate spanwise lift distribution for subsonic flow.
    
    Uses elliptical approximation for first-order estimate.
    
    Args:
        state: State dictionary with wing geometry
        cl: Total lift coefficient
        
    Returns:
        Array of spanwise lift distribution
    """
    # Get wing parameters
    span = state.get('wing_span', 50.0)
    taper = state.get('wing_taper_ratio', 0.5)
    
    # Create spanwise stations
    n_stations = 20
    y = np.linspace(0, span / 2.0, n_stations)  # Half-span
    
    # Elliptical distribution (classical result)
    # cl_local / cl_avg = π/4 * sqrt(1 - (y/(b/2))²)
    cl_distribution = (np.pi / 4.0) * np.sqrt(1.0 - (y / (span / 2.0))**2)
    
    # Scale to match total CL
    cl_distribution *= cl
    
    # Taper effect (modify distribution)
    # More taper → more root loading
    taper_effect = 1.0 + 0.3 * (1.0 - taper) * (1.0 - 2.0 * y / span)
    cl_distribution *= taper_effect
    
    return cl_distribution

