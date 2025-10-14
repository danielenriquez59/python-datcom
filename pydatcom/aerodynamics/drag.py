"""
Drag coefficient calculations for PyDATCOM.

Implements drag buildup method:
- Zero-lift drag (skin friction + form + wave)
- Induced drag (lift-dependent)
- Compressibility effects
- Interference drag

Reference: datcom.f lines 3239 (CDRAG), 3696 (CDWBT)
"""

import numpy as np
from typing import Dict, Optional
import logging

from pydatcom.utils.table_lookup import fig26

logger = logging.getLogger(__name__)


def calculate_skin_friction_drag(state: Dict, mach: float, reynolds: float) -> float:
    """
    Calculate skin friction drag coefficient.
    
    Uses flat plate turbulent skin friction with form factor.
    
    Args:
        state: State dictionary with geometry
        mach: Mach number
        reynolds: Reynolds number
        
    Returns:
        Skin friction drag coefficient
    """
    # Get skin friction coefficient from FIG26
    cf = fig26(reynolds, mach)
    
    # Get wetted area ratio (default 2.5 for typical aircraft)
    swet_sref = state.get('body_swet_sref', 2.5)
    
    # Form factor (accounts for 3D effects)
    # Typical values: 1.2-1.5 for fuselages, 1.1-1.3 for wings
    form_factor = state.get('body_form_factor', 1.3)
    
    # CD0_friction = CF * (Swet/Sref) * Form_Factor
    cd_friction = cf * swet_sref * form_factor
    
    return cd_friction


def calculate_induced_drag(cl: float, aspect_ratio: float, 
                           efficiency_factor: float = None) -> float:
    """
    Calculate induced drag coefficient.
    
    CD_i = CL² / (π * AR * e)
    
    Args:
        cl: Lift coefficient
        aspect_ratio: Wing aspect ratio
        efficiency_factor: Oswald efficiency (calculated if None)
        
    Returns:
        Induced drag coefficient
    """
    if aspect_ratio <= 0:
        return 0.0
    
    # Calculate Oswald efficiency if not provided
    if efficiency_factor is None:
        efficiency_factor = calculate_oswald_efficiency(aspect_ratio)
    
    cdi = cl**2 / (np.pi * aspect_ratio * efficiency_factor)
    
    return cdi


def calculate_oswald_efficiency(aspect_ratio: float, 
                                taper_ratio: float = 0.5,
                                sweep_deg: float = 0.0) -> float:
    """
    Calculate Oswald efficiency factor.
    
    Args:
        aspect_ratio: Wing aspect ratio
        taper_ratio: Wing taper ratio
        sweep_deg: Quarter-chord sweep angle (degrees)
        
    Returns:
        Oswald efficiency factor (0.7-0.98)
    """
    # Empirical correlation (Raymer)
    e = 1.78 * (1.0 - 0.045 * aspect_ratio**0.68) - 0.64
    
    # Taper correction
    e *= (1.0 - 0.05 * abs(taper_ratio - 0.4))
    
    # Sweep correction
    if abs(sweep_deg) > 1.0:
        sweep_rad = np.deg2rad(abs(sweep_deg))
        e *= (1.0 - 0.1 * sweep_rad / (np.pi / 4.0))
    
    # Clip to reasonable range
    e = np.clip(e, 0.7, 0.98)
    
    return e


def calculate_wave_drag_subsonic(mach: float, thickness_ratio: float) -> float:
    """
    Calculate wave drag for high subsonic Mach numbers.
    
    Uses critical Mach number and drag divergence.
    
    Args:
        mach: Mach number
        thickness_ratio: Airfoil thickness ratio (t/c)
        
    Returns:
        Wave drag coefficient
    """
    # Critical Mach number (Korn equation approximation)
    mcrit = 0.87 - thickness_ratio - 0.1 * 0.0  # Simplified, no CL effect
    
    if mach < mcrit:
        return 0.0
    
    # Drag divergence (Lockwood-Rubert)
    mdd = mcrit + 0.1  # Drag divergence Mach
    
    if mach < mdd:
        # Gradual rise to drag divergence
        cd_wave = 20.0 * (mach - mcrit)**4
    else:
        # Post drag divergence
        cd_wave = 0.002 + 20.0 * (mdd - mcrit)**4
        cd_wave += 0.01 * (mach - mdd)
    
    return cd_wave


def calculate_wave_drag_supersonic(mach: float, thickness_ratio: float,
                                   aspect_ratio: float) -> float:
    """
    Calculate wave drag for supersonic flow.
    
    Uses linearized theory with finite span corrections.
    
    Args:
        mach: Mach number (> 1.0)
        thickness_ratio: Airfoil thickness ratio
        aspect_ratio: Wing aspect ratio
        
    Returns:
        Wave drag coefficient
    """
    if mach <= 1.0:
        return 0.0
    
    beta = np.sqrt(mach**2 - 1.0)
    
    # Linearized wave drag (simplified)
    # CD_wave ≈ 4 * (t/c)² / beta
    cd_wave = 4.0 * thickness_ratio**2 / beta
    
    # Finite span correction
    if aspect_ratio > 0:
        ar_correction = aspect_ratio / (aspect_ratio + 4.0 / beta)
        cd_wave *= ar_correction
    
    return cd_wave


def calculate_total_drag(state: Dict, cl: float, mach: float, 
                        reynolds: float) -> Dict[str, float]:
    """
    Calculate total drag coefficient.
    
    Drag buildup: CD = CD0 + CDi + CD_wave
    
    Args:
        state: State dictionary with geometry
        cl: Lift coefficient
        mach: Mach number
        reynolds: Reynolds number
        
    Returns:
        Dictionary with drag components
    """
    # Get geometry parameters
    aspect_ratio = state.get('wing_aspect_ratio', 6.0)
    taper_ratio = state.get('wing_taper_ratio', 0.5)
    thickness_ratio = state.get('wing_tovc', 0.12)
    sweep_deg = state.get('wing_savsi', 0.0)
    
    # Calculate Oswald efficiency
    oswald_e = calculate_oswald_efficiency(aspect_ratio, taper_ratio, sweep_deg)
    
    # Component drag buildup
    cd_friction = calculate_skin_friction_drag(state, mach, reynolds)
    cd_induced = calculate_induced_drag(cl, aspect_ratio, oswald_e)
    
    # Wave drag
    if mach < 0.9:
        # Subsonic
        cd_wave = calculate_wave_drag_subsonic(mach, thickness_ratio)
    elif mach > 1.2:
        # Supersonic
        cd_wave = calculate_wave_drag_supersonic(mach, thickness_ratio, aspect_ratio)
    else:
        # Transonic (interpolate)
        cd_sub = calculate_wave_drag_subsonic(0.9, thickness_ratio)
        cd_sup = calculate_wave_drag_supersonic(1.2, thickness_ratio, aspect_ratio)
        frac = (mach - 0.9) / 0.3
        cd_wave = cd_sub + frac * (cd_sup - cd_sub)
    
    # Miscellaneous drag (interference, protuberances, etc.)
    cd_misc = state.get('drag_misc', 0.0015)  # Typical small aircraft
    
    # Total drag
    cd_total = cd_friction + cd_induced + cd_wave + cd_misc
    
    return {
        'cd_total': cd_total,
        'cd_friction': cd_friction,
        'cd_induced': cd_induced,
        'cd_wave': cd_wave,
        'cd_misc': cd_misc,
        'oswald_e': oswald_e,
    }


class DragCalculator:
    """
    Comprehensive drag calculations for aircraft.
    
    Handles drag buildup across all flight regimes.
    Reference: datcom.f CDRAG (line 3239), CDWBT (line 3696)
    """
    
    def __init__(self, state: Dict):
        """
        Initialize drag calculator.
        
        Args:
            state: Global state dictionary
        """
        self.state = state
    
    def calculate_drag(self, cl: float, mach: float, reynolds: float) -> Dict[str, float]:
        """
        Calculate total drag at given conditions.
        
        Args:
            cl: Lift coefficient
            mach: Mach number
            reynolds: Reynolds number
            
        Returns:
            Dictionary with drag breakdown
        """
        return calculate_total_drag(self.state, cl, mach, reynolds)
    
    def calculate_drag_polar(self, mach: float, reynolds: float,
                            cl_range: np.ndarray = None) -> Dict[str, np.ndarray]:
        """
        Calculate drag polar (CD vs CL).
        
        Args:
            mach: Mach number
            reynolds: Reynolds number
            cl_range: Array of CL values (default: -0.5 to 2.0)
            
        Returns:
            Dictionary with CL and CD arrays
        """
        if cl_range is None:
            cl_range = np.linspace(-0.5, 2.0, 26)
        
        cd_values = np.zeros_like(cl_range)
        
        for i, cl in enumerate(cl_range):
            drag_result = self.calculate_drag(cl, mach, reynolds)
            cd_values[i] = drag_result['cd_total']
        
        return {
            'cl': cl_range,
            'cd': cd_values,
            'mach': mach,
            'reynolds': reynolds,
        }

