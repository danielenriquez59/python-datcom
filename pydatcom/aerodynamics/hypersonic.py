"""
Hypersonic aerodynamic calculations for PyDATCOM.

Implements hypersonic flow methods (Mach > 5):
- Newtonian impact theory
- Modified Newtonian method
- Hypersonic shock relations

Reference: datcom.f HYPBOD, HYPFLP, HYPROP subroutines
"""

import numpy as np
from typing import Dict
import logging

logger = logging.getLogger(__name__)


def calculate_hypersonic_coefficients(state: Dict, alpha_deg: float,
                                      mach: float) -> Dict[str, float]:
    """
    Calculate aerodynamic coefficients for hypersonic flow.
    
    Uses modified Newtonian impact theory.
    
    Args:
        state: State dictionary
        alpha_deg: Angle of attack (degrees)
        mach: Mach number (> 5.0)
        
    Returns:
        Dictionary with hypersonic coefficients
    """
    alpha_rad = np.deg2rad(alpha_deg)
    
    # Modified Newtonian theory
    # CP_max = 2 for M → ∞, approximation: CP_max ≈ 1.84 for M=5
    if mach >= 5.0:
        cp_max = 1.84 + 0.16 * (mach - 5.0) / 5.0
        cp_max = min(cp_max, 2.0)
    else:
        # Transition region (M=3 to M=5)
        cp_max = 1.5 + 0.34 * (mach - 3.0) / 2.0
    
    # Normal force coefficient (Newtonian)
    # CN = CP_max * sin²(α)
    cn = cp_max * np.sin(alpha_rad)**2
    
    # Axial force (base drag dominates)
    ca_base = 0.2  # Typical base drag coefficient
    
    # Lift and drag
    cl = cn * np.cos(alpha_rad) - ca_base * np.sin(alpha_rad)
    cd = ca_base * np.cos(alpha_rad) + cn * np.sin(alpha_rad)
    
    # Moment (center of pressure near 0.5c for hypersonic)
    xcp = 0.5  # Center of pressure location
    xcg = state.get('synths_xcg', 0.0)
    cbar = state.get('options_cbarr', 1.0)
    
    if cbar > 0:
        cm = -cn * (xcp - xcg) / cbar
    else:
        cm = 0.0
    
    return {
        'cl': cl,
        'cd': cd,
        'cm': cm,
        'cn': cn,
        'ca': ca_base,
        'cp_max': cp_max,
        'xcp': xcp,
        'regime': 'hypersonic',
        'mach': mach,
        'alpha': alpha_deg,
    }

