"""
Main aerodynamic calculator for PyDATCOM.

Routes calculations to appropriate flow regime (subsonic, transonic, supersonic, hypersonic).
Coordinates lift, drag, and moment calculations.

Reference: datcom.f main computational flow
"""

import numpy as np
from typing import Dict, List
import logging

from pydatcom.aerodynamics.subsonic import calculate_subsonic_coefficients
from pydatcom.aerodynamics.transonic import calculate_transonic_coefficients
from pydatcom.aerodynamics.supersonic import calculate_supersonic_coefficients
from pydatcom.aerodynamics.hypersonic import calculate_hypersonic_coefficients
from pydatcom.aerodynamics.body_alone import has_wing_or_tail, calculate_body_alone_coefficients
from pydatcom.utils.atmosphere import Atmosphere
from pydatcom.geometry.wing import calculate_straight_exposed_geometry

logger = logging.getLogger(__name__)


class AerodynamicCalculator:
    """
    Main aerodynamic calculator for PyDATCOM.
    
    Automatically selects appropriate methods based on flight regime.
    Computes CL, CD, Cm across all Mach numbers and angles of attack.
    """
    
    def __init__(self, state: Dict):
        """
        Initialize calculator with state dictionary.
        
        Args:
            state: Global state dictionary
        """
        self.state = state
    
    def identify_regime(self, mach: float) -> str:
        """
        Identify flow regime from Mach number.
        
        Args:
            mach: Mach number
            
        Returns:
            Regime name: 'subsonic', 'transonic', 'supersonic', or 'hypersonic'
        """
        if mach < 0.9:
            return 'subsonic'
        elif mach < 1.2:
            return 'transonic'
        elif mach < 5.0:
            return 'supersonic'
        else:
            return 'hypersonic'
    
    def calculate_at_condition(self, alpha_deg: float, mach: float,
                               reynolds: float = None,
                               condition_index: int = None) -> Dict[str, float]:
        """
        Calculate aerodynamic coefficients at single flight condition.
        
        Automatically detects body-only vs wing configurations.
        
        Args:
            alpha_deg: Angle of attack (degrees)
            mach: Mach number
            reynolds: Dimensionless Reynolds number based on the relevant
                component/reference length. When omitted, DATCOM's RNNUB
                (Reynolds number per unit length) is resolved for this Mach
                and multiplied by that length.
            condition_index: Optional zero-based FLTCON condition index.
            
        Returns:
            Dictionary with CL, CD, Cm and components
        """
        # Estimate Reynolds number if not provided
        if reynolds is None:
            reynolds = self._estimate_reynolds(mach, condition_index)
        
        # Check if this is a body-only configuration
        if not has_wing_or_tail(self.state):
            # Use body-alone methods (like BODYRT in FORTRAN)
            logger.info("Body-only configuration detected, using body-alone methods")
            result = calculate_body_alone_coefficients(self.state, alpha_deg, mach, reynolds)
            return result
        
        # Identify regime for wing configurations
        regime = self.identify_regime(mach)
        
        # Route to appropriate calculator
        if regime == 'subsonic':
            result = calculate_subsonic_coefficients(self.state, alpha_deg, mach, reynolds)
        elif regime == 'transonic':
            result = calculate_transonic_coefficients(self.state, alpha_deg, mach, reynolds)
        elif regime == 'supersonic':
            result = calculate_supersonic_coefficients(self.state, alpha_deg, mach, reynolds)
        else:  # hypersonic
            result = calculate_hypersonic_coefficients(self.state, alpha_deg, mach)
            result['reynolds'] = reynolds
        
        # Add regime info
        result['regime'] = regime
        
        return result
    
    def calculate_alpha_sweep(self, alpha_range: np.ndarray, mach: float,
                              reynolds: float = None) -> Dict[str, np.ndarray]:
        """
        Calculate coefficients across range of angles of attack.
        
        Args:
            alpha_range: Array of angles of attack (degrees)
            mach: Mach number
            reynolds: Reynolds number
            
        Returns:
            Dictionary with arrays of CL, CD, Cm vs alpha
        """
        n_alpha = len(alpha_range)
        
        # Pre-allocate arrays
        cl_array = np.zeros(n_alpha)
        cd_array = np.zeros(n_alpha)
        cm_array = np.zeros(n_alpha)
        
        # Calculate at each alpha
        for i, alpha in enumerate(alpha_range):
            result = self.calculate_at_condition(alpha, mach, reynolds)
            cl_array[i] = result['cl']
            cd_array[i] = result['cd']
            cm_array[i] = result['cm']
        
        return {
            'alpha': alpha_range,
            'cl': cl_array,
            'cd': cd_array,
            'cm': cm_array,
            'mach': mach,
            'reynolds': reynolds if reynolds is not None else self._estimate_reynolds(mach),
            'regime': self.identify_regime(mach),
        }
    
    def calculate_mach_sweep(self, alpha_deg: float, mach_range: np.ndarray) -> Dict[str, np.ndarray]:
        """
        Calculate coefficients across range of Mach numbers.
        
        Args:
            alpha_deg: Angle of attack (degrees)
            mach_range: Array of Mach numbers
            
        Returns:
            Dictionary with arrays of CL, CD, Cm vs Mach
        """
        n_mach = len(mach_range)
        
        cl_array = np.zeros(n_mach)
        cd_array = np.zeros(n_mach)
        cm_array = np.zeros(n_mach)
        
        reynolds_array = np.zeros(n_mach)
        for i, mach in enumerate(mach_range):
            reynolds = self._estimate_reynolds(mach, i)
            result = self.calculate_at_condition(alpha_deg, mach, reynolds, i)
            cl_array[i] = result['cl']
            cd_array[i] = result['cd']
            cm_array[i] = result['cm']
            reynolds_array[i] = reynolds
        
        return {
            'mach': mach_range,
            'cl': cl_array,
            'cd': cd_array,
            'cm': cm_array,
            'alpha': alpha_deg,
            'reynolds': reynolds_array,
        }
    
    @staticmethod
    def _as_list(value) -> List[float]:
        """Return a scalar/array state entry as a plain list."""
        if value is None:
            return []
        if np.isscalar(value):
            return [value]
        return list(value)

    def _condition_index(self, mach: float, requested: int = None) -> int:
        machs = self._as_list(self.state.get('flight_mach'))
        if requested is not None and 0 <= requested < len(machs):
            if np.isclose(machs[requested], mach, rtol=1e-9, atol=1e-12):
                return requested
        for index, scheduled_mach in enumerate(machs):
            if np.isclose(scheduled_mach, mach, rtol=1e-9, atol=1e-12):
                return index
        # A direct off-schedule API call has no FORTRAN loop index. Retain
        # the first supplied condition rather than inventing interpolation.
        return 0

    def _reference_length(self) -> float:
        """Length used to dimensionalize RNNUB for the active configuration."""
        if not has_wing_or_tail(self.state):
            body_x = self._as_list(self.state.get('body_x'))
            length = self.state.get('body_length')
            if length is None and len(body_x) >= 2:
                length = body_x[-1] - body_x[0]
        else:
            # CDRAG uses A(16), the exposed mean aerodynamic chord. CBARR is
            # the best available fallback until the complete geometry COMMON
            # layout is translated.
            length = (self.state.get('wing_mac') or
                      self.state.get('options_cbarr'))
            if (not self.state.get('wing_mac') and
                    float(self.state.get('wing_type', 1.0) or 1.0) == 1.0):
                try:
                    length = calculate_straight_exposed_geometry(self.state)['mac']
                except ValueError:
                    pass
            if not length:
                root = self.state.get('wing_chrdr')
                tip = self.state.get('wing_chrdtp')
                if root and tip is not None and root > 0.0:
                    taper = float(tip) / float(root)
                    length = (2.0 / 3.0) * float(root) * (
                        (1.0 + taper + taper**2) / (1.0 + taper))
                elif root:
                    # Partial continuation cases can omit the prior panel
                    # definition. Root chord is the only dimensional surface
                    # length available until case inheritance is translated.
                    length = root
        if length is None or not np.isfinite(length) or length <= 0:
            raise ValueError("A positive characteristic length is required to convert RNNUB")
        return float(length)

    def _estimate_reynolds(self, mach: float, condition_index: int = None) -> float:
        """
        Estimate Reynolds number from Mach and state.
        
        Args:
            mach: Mach number
            
        Returns:
            Estimated Reynolds number
        """
        index = self._condition_index(mach, condition_index)
        rnnub = self._as_list(self.state.get('flight_rnnub'))
        if rnnub:
            per_length = rnnub[min(index, len(rnnub) - 1)]
        else:
            pressures = self._as_list(self.state.get('flight_pinf'))
            temperatures = self._as_list(self.state.get('flight_tinf'))
            altitudes = self._as_list(self.state.get('flight_alt'))
            loop = int(self.state.get('flight_loop', 1) or 1)
            atmosphere_index = index if loop == 1 else 0
            if pressures and temperatures:
                pressure = pressures[min(atmosphere_index, len(pressures) - 1)]
                temperature = temperatures[min(atmosphere_index, len(temperatures) - 1)]
            elif altitudes:
                altitude = altitudes[min(atmosphere_index, len(altitudes) - 1)]
                atmosphere = Atmosphere.calculate(float(altitude))
                pressure = atmosphere['pressure']
                temperature = atmosphere['temperature']
            else:
                # INPUT labels 1130-1160 use this per-unit-length fallback.
                per_length = 5.0e6
                return per_length * self._reference_length()
            # Main program labels 1030/1050/1080. PINF is psf, TINF Rankine.
            per_length = (1.2527e6 * pressure * mach *
                          (temperature + 198.6) / temperature**2)
        if not np.isfinite(per_length) or per_length <= 0:
            raise ValueError("RNNUB must be a positive Reynolds number per unit length")
        return float(per_length) * self._reference_length()


def calculate_aero_coefficients(state: Dict, alpha_deg: float, mach: float,
                                reynolds: float = None) -> Dict[str, float]:
    """
    Convenience function to calculate coefficients at any Mach number.
    
    Automatically selects appropriate regime.
    
    Args:
        state: State dictionary
        alpha_deg: Angle of attack (degrees)
        mach: Mach number
        reynolds: Reynolds number (estimated if None)
        
    Returns:
        Dictionary with CL, CD, Cm
    """
    calculator = AerodynamicCalculator(state)
    return calculator.calculate_at_condition(alpha_deg, mach, reynolds)

