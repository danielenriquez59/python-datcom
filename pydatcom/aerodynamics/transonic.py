"""
Transonic aerodynamic calculations for PyDATCOM.

Implements transonic flow methods (0.9 < Mach < 1.2):
- Interpolation between subsonic and supersonic
- Drag divergence effects
- Shock formation

Reference: datcom.f transonic calculation sections
"""

import numpy as np
from typing import Dict
import logging

from pydatcom.aerodynamics.subsonic import calculate_subsonic_coefficients
from pydatcom.aerodynamics.supersonic import calculate_supersonic_coefficients

logger = logging.getLogger(__name__)


def calculate_transonic_coefficients(state: Dict, alpha_deg: float,
                                     mach: float, reynolds: float) -> Dict[str, float]:
    """
    Calculate aerodynamic coefficients in transonic regime.
    
    Uses interpolation between subsonic and supersonic methods.
    
    Args:
        state: State dictionary
        alpha_deg: Angle of attack (degrees)
        mach: Mach number (0.9 < M < 1.2)
        reynolds: Reynolds number
        
    Returns:
        Dictionary with interpolated coefficients
    """
    if mach < 0.9:
        return calculate_subsonic_coefficients(state, alpha_deg, mach, reynolds)
    elif mach > 1.2:
        return calculate_supersonic_coefficients(state, alpha_deg, mach, reynolds)
    
    # Interpolate between subsonic (M=0.9) and supersonic (M=1.2)
    sub_result = calculate_subsonic_coefficients(state, alpha_deg, 0.9, reynolds)
    sup_result = calculate_supersonic_coefficients(state, alpha_deg, 1.2, reynolds)
    
    # Interpolation factor
    frac = (mach - 0.9) / 0.3
    
    # Interpolate coefficients
    cl = sub_result['cl'] + frac * (sup_result['cl'] - sub_result['cl'])
    cd = sub_result['cd'] + frac * (sup_result['cd'] - sub_result['cd'])
    cm = sub_result['cm'] + frac * (sup_result['cm'] - sub_result['cm'])
    
    # Drag divergence effect (additional drag in transonic)
    cd_divergence = 0.01 * np.sin(np.pi * frac)**2  # Peak at M=1.05
    cd += cd_divergence
    
    return {
        'cl': cl,
        'cd': cd,
        'cm': cm,
        'cd_divergence': cd_divergence,
        'regime': 'transonic',
        'mach': mach,
        'alpha': alpha_deg,
        'interpolation_factor': frac,
    }

