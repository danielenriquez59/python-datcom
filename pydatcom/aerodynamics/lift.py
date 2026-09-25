"""
Lift coefficient calculations for PyDATCOM.

Implements lift calculations across all flight regimes:
- Subsonic (incompressible and compressible)
- Transonic
- Supersonic
- Hypersonic

Reference: datcom.f lines 4122 (CLMCH0), 6549 (CSLOPE), 2589 (CALCA)
"""

import math
import numpy as np
from typing import Dict, Optional, Sequence, Tuple, Union
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
                                            sweep_angle_deg: float = 0.0,
                                            section_cla_per_deg: Optional[float] = None) -> float:
    """
    Calculate the DATCOM linear finite-wing lift-curve slope.
    
    Args:
        aspect_ratio: Wing aspect ratio
        taper_ratio: Wing taper ratio
        mach: Mach number
        sweep_angle_deg: Quarter-chord sweep angle (degrees)
        
    Returns:
        Compressible lift curve slope (per radian)

    Notes:
        This is the analytic finite-wing expression used by ``TRSONI`` to
        establish the subsonic lift-slope points.  ``sweep_angle_deg`` is the
        half-chord sweep used by that expression.  DATCOM's optional section
        ``CLALPA`` input is per degree; when supplied it defines the airfoil
        factor ``K = CLALPA * RAD / (2*pi)``.
    """
    del taper_ratio  # Taper enters other DATCOM corrections, not this equation.
    if aspect_ratio <= 0.0:
        return 0.0
    if mach < 0.0 or mach >= 1.0:
        raise ValueError("DATCOM subsonic lift slope requires 0 <= Mach < 1")

    if section_cla_per_deg is None:
        section_factor = 1.0
    else:
        if section_cla_per_deg <= 0.0:
            raise ValueError("section CLALPA must be positive")
        section_factor = section_cla_per_deg * np.rad2deg(1.0) / (2.0 * np.pi)

    beta_squared = 1.0 - mach**2
    tan_half_chord_sweep = math.tan(np.deg2rad(sweep_angle_deg))
    ar_over_k_squared = (aspect_ratio / section_factor) ** 2
    denominator = 2.0 + math.sqrt(
        4.0 + ar_over_k_squared *
        (beta_squared + tan_half_chord_sweep**2)
    )
    return 2.0 * np.pi * aspect_ratio / denominator


def _mach_indexed_value(
    value: Union[float, Sequence[Optional[float]], None],
    state: Dict,
    mach: float,
) -> Optional[float]:
    """Resolve a scalar or Mach-indexed DATCOM input from the state."""
    if value is None:
        return None
    if np.isscalar(value):
        return float(value)

    usable = [(mach_slot, item) for mach_slot, item in enumerate(value)
              if item is not None and abs(float(item)) != UNUSED]
    if not usable:
        return None

    mach_schedule = state.get('flight_mach') or []
    candidates = [(mach_slot, item) for mach_slot, item in usable
                  if (mach_slot < len(mach_schedule) and
                      mach_schedule[mach_slot] is not None)]
    if candidates:
        _, item = min(candidates,
                      key=lambda pair: abs(float(mach_schedule[pair[0]]) - mach))
        return float(item)
    return float(usable[0][1])


def _half_chord_sweep_deg(state: Dict) -> float:
    """Convert WGPLNF sweep at CHSTAT to the half-chord sweep used by DATCOM."""
    sweep_reference = float(state.get('wing_savsi', 0.0) or 0.0)
    chord_station = float(state.get('wing_chstat', 0.25) or 0.0)
    root_chord = state.get('wing_chrdr')
    tip_chord = state.get('wing_chrdtp')
    semispan = state.get('wing_sspn')
    if root_chord is None or tip_chord is None or not semispan:
        return sweep_reference

    chord_gradient = (float(tip_chord) - float(root_chord)) / float(semispan)
    tangent = (math.tan(np.deg2rad(sweep_reference)) +
               (0.5 - chord_station) * chord_gradient)
    return float(np.rad2deg(math.atan(tangent)))


def resolve_wing_lift_inputs(state: Dict, mach: float) -> Dict[str, Optional[float]]:
    """Resolve DATCOM section design-point and incidence inputs for wing lift."""
    section_cla_per_deg = _mach_indexed_value(state.get('wing_clalpa'), state, mach)
    if section_cla_per_deg is None:
        section_cla_per_deg = 2.0 * np.pi / np.rad2deg(1.0)

    design_cl = float(state.get('wing_cli', 0.0) or 0.0)
    design_alpha_deg = float(state.get('wing_alphai', 0.0) or 0.0)
    incidence_deg = float(state.get('synths_aliw', 0.0) or 0.0)
    section_alpha_zero_deg = design_alpha_deg - design_cl / section_cla_per_deg
    body_alpha_zero_deg = section_alpha_zero_deg - incidence_deg
    return {
        'section_cla_per_deg': section_cla_per_deg,
        'design_cl': design_cl,
        'design_alpha': design_alpha_deg,
        'incidence': incidence_deg,
        'section_alpha_zero': section_alpha_zero_deg,
        'alpha_zero': body_alpha_zero_deg,
    }


def wing_reference_ratio(state: Dict) -> float:
    """Return SW/SREF, the component-to-aircraft coefficient scale."""
    wing_area = float(state.get('wing_area', 0.0) or 0.0)
    sref = float(state.get('options_sref', wing_area) or wing_area)
    if sref <= 0.0:
        raise ValueError("options_sref must be positive")
    if wing_area == 0.0:
        # A standalone wing coefficient already uses its own implicit area.
        wing_area = sref
    if wing_area < 0.0:
        raise ValueError("wing_area cannot be negative")
    return wing_area / sref


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
    sweep_deg = _half_chord_sweep_deg(state)
    
    # If not computed yet, derive it from the planform.  Prefer the
    # translated WTGEOM exposed geometry over the span/area shortcut: a deck
    # that supplies only CHRDR/CHRDTP/SSPN/SSPNE otherwise silently fell back
    # to AR=6.0, which is far from the exposed AR of a typical DATCOM case.
    if aspect_ratio is None or aspect_ratio == 6.0:
        resolved = None
        if float(state.get('wing_type', 1.0) or 1.0) == 1.0:
            try:
                from pydatcom.geometry.wing import calculate_straight_exposed_geometry
                geometry = calculate_straight_exposed_geometry(state)
                resolved = geometry['aspect_ratio']
                if taper_ratio is None or taper_ratio == 0.5:
                    taper_ratio = geometry['taper_ratio']
            except ValueError:
                resolved = None
        if resolved is None:
            span = state.get('wing_span')
            area = state.get('wing_area') or state.get('options_sref')
            if span and area:
                resolved = span**2 / area
        if resolved is None:
            logger.warning("Wing aspect ratio not available, using default 6.0")
            resolved = 6.0
        aspect_ratio = resolved
    
    # Calculate lift curve slope
    section = resolve_wing_lift_inputs(state, mach)
    cla = calculate_lift_curve_slope_compressible(
        aspect_ratio, taper_ratio, mach, sweep_deg,
        section['section_cla_per_deg'],
    )

    alpha_zero = section['alpha_zero']
    
    # DATCOM surface methods first produce coefficients on the exposed wing
    # area, then multiply by SRSTAR/SR when assembling aircraft outputs.
    reference_ratio = wing_reference_ratio(state)
    cl_wing = calculate_lift_coefficient(alpha_deg, alpha_zero, cla)
    cl = cl_wing * reference_ratio
    cla_aircraft = cla * reference_ratio
    
    return {
        'cl': cl,
        'cl_wing': cl_wing,
        'cla': cla_aircraft,
        'cla_wing': cla,
        'cla_per_deg': cla_aircraft * np.deg2rad(1.0),
        'cla_wing_per_deg': cla * np.deg2rad(1.0),
        'wing_reference_ratio': reference_ratio,
        'alpha_zero': alpha_zero,
        'section_alpha_zero': section['section_alpha_zero'],
        'incidence': section['incidence'],
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
            # Match the public transonic calculator's endpoint convention so
            # both APIs are continuous at Mach 0.9 and 1.2.
            lower = calculate_wing_lift_subsonic(self.state, alpha_deg, 0.9)
            upper = self._calculate_wing_lift_supersonic(alpha_deg, 1.2)
            fraction = (mach - 0.9) / 0.3
            return {
                'cl': lower['cl'] + fraction * (upper['cl'] - lower['cl']),
                'cla': lower['cla'] + fraction * (upper['cla'] - lower['cla']),
                'cla_per_deg': (lower['cla_per_deg'] + fraction *
                                (upper['cla_per_deg'] - lower['cla_per_deg'])),
                'alpha_zero': lower['alpha_zero'],
                'regime': 'transonic',
            }
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
        
        # Calculate CL on the wing reference, then assemble it on SREF.
        section = resolve_wing_lift_inputs(self.state, mach)
        alpha_zero = section['alpha_zero']
        alpha_eff_rad = np.deg2rad(alpha_deg - alpha_zero)
        ratio = wing_reference_ratio(self.state)
        cl_wing = cla * alpha_eff_rad
        cla_aircraft = cla * ratio
        cl = cl_wing * ratio
        
        return {
            'cl': cl,
            'cl_wing': cl_wing,
            'cla': cla_aircraft,
            'cla_per_deg': cla_aircraft * np.deg2rad(1.0),
            'alpha_zero': alpha_zero,
            'regime': 'supersonic',
            'beta': beta,
        }

