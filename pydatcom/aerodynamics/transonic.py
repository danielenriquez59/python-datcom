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
    
    subsonic_anchor = calculate_subsonic_coefficients(
        state, alpha_deg, 0.9, reynolds,
    )
    supersonic_anchor = calculate_supersonic_coefficients(
        state, alpha_deg, 1.2, reynolds,
    )

    mach_fraction = (mach - 0.9) / 0.3

    cl = (subsonic_anchor['cl'] +
          mach_fraction * (supersonic_anchor['cl'] - subsonic_anchor['cl']))
    cd = (subsonic_anchor['cd'] +
          mach_fraction * (supersonic_anchor['cd'] - subsonic_anchor['cd']))
    cm = (subsonic_anchor['cm'] +
          mach_fraction * (supersonic_anchor['cm'] - subsonic_anchor['cm']))

    cd_divergence = 0.01 * np.sin(np.pi * mach_fraction) ** 2
    cd += cd_divergence

    return {
        'cl': cl,
        'cd': cd,
        'cm': cm,
        'cd_divergence': cd_divergence,
        'regime': 'transonic',
        'mach': mach,
        'alpha': alpha_deg,
        'interpolation_factor': mach_fraction,
    }

