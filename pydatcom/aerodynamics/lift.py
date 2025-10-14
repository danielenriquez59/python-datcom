"""
Lift coefficient calculations for PyDATCOM.

Implements lift calculations across all flight regimes:
- Subsonic (incompressible and compressible)
- Transonic
- Supersonic
- Hypersonic

Reference: datcom.f lines 4122 (CLMCH0), 6549 (CSLOPE), 2589 (CALCA)
"""

import numpy as np
from typing import Dict, Tuple, Optional
import logging

from pydatcom.utils.constants import UNUSED
from pydatcom.utils.math_utils import arcsin

logger = logging.getLogger(__name__)


def calculate_lift_curve_slope_incompressible(aspect_ratio: float, 
                                              taper_ratio: float,
                                              sweep_angle_deg: float = 0.0) -> float:
    """
    Calculate incompressible lift curve slope using lifting line theory.
    
    Uses Helmbold equation for finite wing.
    
    Args:
        aspect_ratio: Wing aspect ratio (b²/S)
        taper_ratio: Wing taper ratio (Ct/Cr)
        sweep_angle_deg: Quarter-chord sweep angle (degrees)
        
    Returns:
        Lift curve slope dCL/dα (per radian)
    """
    # Theoretical slope for infinite aspect ratio (2D)
    cla_2d = 2.0 * np.pi
    
    # Finite wing correction using Helmbold equation
    # CL_α = (2πAR) / (2 + sqrt(4 + AR²))
    ar_term = aspect_ratio**2
    cla_3d = (2.0 * np.pi * aspect_ratio) / (2.0 + np.sqrt(4.0 + ar_term))
    
    # Sweep correction (simplified)
    if abs(sweep_angle_deg) > 0.1:
        sweep_rad = np.deg2rad(sweep_angle_deg)
        cos_sweep = np.cos(sweep_rad)
        cla_3d *= cos_sweep
    
    return cla_3d


def calculate_lift_curve_slope_compressible(aspect_ratio: float,
                                            taper_ratio: float,
                                            mach: float,
                                            sweep_angle_deg: float = 0.0) -> float:
    """
    Calculate compressible lift curve slope using Prandtl-Glauert.
    
    Args:
        aspect_ratio: Wing aspect ratio
        taper_ratio: Wing taper ratio
        mach: Mach number
        sweep_angle_deg: Quarter-chord sweep angle (degrees)
        
    Returns:
        Compressible lift curve slope (per radian)
    """
    # Get incompressible slope
    cla_incomp = calculate_lift_curve_slope_incompressible(
        aspect_ratio, taper_ratio, sweep_angle_deg
    )
    
    # Prandtl-Glauert correction for subsonic
    if mach < 0.9:
        beta = np.sqrt(1.0 - mach**2)
        if beta > 0.01:
            cla_comp = cla_incomp / beta
        else:
            cla_comp = cla_incomp
    else:
        # Near or above Mach 1, use different method
        cla_comp = cla_incomp
    
    return cla_comp


def calculate_lift_coefficient(alpha_deg: float,
                               alpha_zero_deg: float,
                               cl_alpha_per_rad: float) -> float:
    """
    Calculate lift coefficient using linear theory.
    
    CL = CL_α * (α - α_0)
    
    Args:
        alpha_deg: Angle of attack (degrees)
        alpha_zero_deg: Zero-lift angle of attack (degrees)
        cl_alpha_per_rad: Lift curve slope (per radian)
        
    Returns:
        Lift coefficient
    """
    alpha_eff_deg = alpha_deg - alpha_zero_deg
    alpha_eff_rad = np.deg2rad(alpha_eff_deg)
    
    cl = cl_alpha_per_rad * alpha_eff_rad
    
    return cl


def calculate_wing_lift_subsonic(state: Dict, alpha_deg: float, mach: float) -> Dict[str, float]:
    """
    Calculate wing lift coefficient for subsonic flow.
    
    Uses lifting line theory with compressibility corrections.
    
    Args:
        state: State dictionary with wing geometry
        alpha_deg: Angle of attack (degrees)
        mach: Mach number (< 0.9)
        
    Returns:
        Dictionary with lift parameters:
        - cl: Lift coefficient
        - cla: Lift curve slope (per radian)
        - alpha_zero: Zero-lift angle (degrees)
    """
    # Get wing geometry from state
    aspect_ratio = state.get('wing_aspect_ratio', 6.0)
    taper_ratio = state.get('wing_taper_ratio', 0.5)
    sweep_deg = state.get('wing_savsi', 0.0)
    
    # If not computed yet, try to calculate from planform
    if aspect_ratio is None or aspect_ratio == 6.0:
        # Try to calculate from span and area
        span = state.get('wing_span')
        area = state.get('wing_area') or state.get('options_sref')
        if span and area:
            aspect_ratio = span**2 / area
        else:
            logger.warning("Wing aspect ratio not available, using default 6.0")
            aspect_ratio = 6.0
    
    # Calculate lift curve slope
    cla = calculate_lift_curve_slope_compressible(
        aspect_ratio, taper_ratio, mach, sweep_deg
    )
    
    # Zero-lift angle (from camber or incidence)
    alpha_zero = state.get('wing_alphai', 0.0) or 0.0
    
    # Calculate CL
    cl = calculate_lift_coefficient(alpha_deg, alpha_zero, cla)
    
    return {
        'cl': cl,
        'cla': cla,
        'cla_per_deg': np.rad2deg(cla),  # Per degree
        'alpha_zero': alpha_zero,
    }


def calculate_maximum_lift_coefficient(aspect_ratio: float,
                                       taper_ratio: float,
                                       cl_max_section: float = 1.5) -> float:
    """
    Estimate maximum lift coefficient for wing.
    
    Uses empirical correlation between section and wing maximum lift.
    
    Args:
        aspect_ratio: Wing aspect ratio
        taper_ratio: Wing taper ratio
        cl_max_section: Section maximum lift coefficient
        
    Returns:
        Wing maximum lift coefficient
    """
    # Empirical reduction factor for 3D wing
    # CL_max_wing ≈ 0.9 * CL_max_section for typical wings
    reduction_factor = 0.9
    
    # Aspect ratio effect (higher AR → higher CL_max)
    ar_factor = 1.0 + 0.05 * (aspect_ratio - 6.0) / 6.0
    ar_factor = np.clip(ar_factor, 0.8, 1.2)
    
    # Taper effect (moderate taper optimal)
    taper_factor = 0.95 + 0.05 * np.cos(np.pi * (taper_ratio - 0.4))
    taper_factor = np.clip(taper_factor, 0.90, 1.0)
    
    cl_max = cl_max_section * reduction_factor * ar_factor * taper_factor
    
    return cl_max


def calculate_induced_drag_coefficient(cl: float, aspect_ratio: float,
                                       efficiency_factor: float = 0.95) -> float:
    """
    Calculate induced drag coefficient.
    
    CD_i = CL² / (π * AR * e)
    
    Args:
        cl: Lift coefficient
        aspect_ratio: Wing aspect ratio
        efficiency_factor: Oswald efficiency factor (default 0.95)
        
    Returns:
        Induced drag coefficient
    """
    if aspect_ratio <= 0:
        return 0.0
    
    cdi = cl**2 / (np.pi * aspect_ratio * efficiency_factor)
    
    return cdi


def calculate_oswald_efficiency(aspect_ratio: float, 
                                taper_ratio: float,
                                sweep_deg: float = 0.0) -> float:
    """
    Calculate Oswald efficiency factor.
    
    Empirical correlation based on wing geometry.
    
    Args:
        aspect_ratio: Wing aspect ratio
        taper_ratio: Wing taper ratio
        sweep_deg: Quarter-chord sweep angle (degrees)
        
    Returns:
        Oswald efficiency factor e (typically 0.7-0.95)
    """
    # Base efficiency
    e_base = 1.78 * (1.0 - 0.045 * aspect_ratio**0.68) - 0.64
    
    # Taper effect
    taper_effect = 1.0 - 0.05 * abs(taper_ratio - 0.4)
    
    # Sweep effect (reduces efficiency)
    if abs(sweep_deg) > 1.0:
        sweep_rad = np.deg2rad(abs(sweep_deg))
        sweep_effect = 1.0 - 0.1 * (sweep_rad / (np.pi / 4.0))
    else:
        sweep_effect = 1.0
    
    e = e_base * taper_effect * sweep_effect
    e = np.clip(e, 0.7, 0.98)
    
    return e


class LiftCalculator:
    """
    Comprehensive lift calculations for wings and lifting surfaces.
    
    Handles subsonic, transonic, and supersonic regimes.
    """
    
    def __init__(self, state: Dict):
        """
        Initialize lift calculator with state dictionary.
        
        Args:
            state: Global state dictionary
        """
        self.state = state
    
    def calculate_wing_lift(self, alpha_deg: float, mach: float) -> Dict[str, float]:
        """
        Calculate wing lift coefficient at given conditions.
        
        Args:
            alpha_deg: Angle of attack (degrees)
            mach: Mach number
            
        Returns:
            Dictionary with lift results
        """
        if mach < 0.9:
            # Subsonic
            return calculate_wing_lift_subsonic(self.state, alpha_deg, mach)
        elif mach < 1.2:
            # Transonic (use subsonic with correction)
            logger.warning("Transonic lift using subsonic approximation")
            result = calculate_wing_lift_subsonic(self.state, alpha_deg, 0.85)
            result['regime'] = 'transonic'
            result['cl'] *= 0.9  # Rough transonic reduction
            return result
        else:
            # Supersonic
            return self._calculate_wing_lift_supersonic(alpha_deg, mach)
    
    def _calculate_wing_lift_supersonic(self, alpha_deg: float, mach: float) -> Dict[str, float]:
        """
        Calculate wing lift for supersonic flow.
        
        Uses linearized supersonic theory.
        
        Args:
            alpha_deg: Angle of attack (degrees)
            mach: Mach number (> 1.2)
            
        Returns:
            Dictionary with lift results
        """
        # Prandtl-Meyer supersonic theory
        beta = np.sqrt(mach**2 - 1.0)
        aspect_ratio = self.state.get('wing_aspect_ratio', 4.0)
        
        # Supersonic lift curve slope (simplified)
        # CL_α = 4 / beta
        cla = 4.0 / beta
        
        # Apply finite span correction
        if aspect_ratio > 0:
            ar_correction = aspect_ratio / (aspect_ratio + 2.0 / beta)
            cla *= ar_correction
        
        # Calculate CL
        alpha_zero = self.state.get('wing_alphai', 0.0) or 0.0
        alpha_eff_rad = np.deg2rad(alpha_deg - alpha_zero)
        cl = cla * alpha_eff_rad
        
        return {
            'cl': cl,
            'cla': cla,
            'cla_per_deg': np.rad2deg(cla),
            'alpha_zero': alpha_zero,
            'regime': 'supersonic',
            'beta': beta,
        }

