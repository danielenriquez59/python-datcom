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
                               reynolds: float = None) -> Dict[str, float]:
        """
        Calculate aerodynamic coefficients at single flight condition.
        
        Automatically detects body-only vs wing configurations.
        
        Args:
            alpha_deg: Angle of attack (degrees)
            mach: Mach number
            reynolds: Reynolds number (estimated if None)
            
        Returns:
            Dictionary with CL, CD, Cm and components
        """
        # Estimate Reynolds number if not provided
        if reynolds is None:
            reynolds = self._estimate_reynolds(mach)
        
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
            'reynolds': reynolds or self._estimate_reynolds(mach),
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
        
        for i, mach in enumerate(mach_range):
            reynolds = self._estimate_reynolds(mach)
            result = self.calculate_at_condition(alpha_deg, mach, reynolds)
            cl_array[i] = result['cl']
            cd_array[i] = result['cd']
            cm_array[i] = result['cm']
        
        return {
            'mach': mach_range,
            'cl': cl_array,
            'cd': cd_array,
            'cm': cm_array,
            'alpha': alpha_deg,
        }
    
    def _estimate_reynolds(self, mach: float) -> float:
        """
        Estimate Reynolds number from Mach and state.
        
        Args:
            mach: Mach number
            
        Returns:
            Estimated Reynolds number
        """
        # Try to get from state
        rnnub_list = self.state.get('flight_rnnub', [])
        if rnnub_list and len(rnnub_list) > 0:
            return rnnub_list[0]
        
        # Estimate from altitude and Mach
        altitude = self.state.get('flight_alt', [0.0])
        if isinstance(altitude, list) and len(altitude) > 0:
            alt = altitude[0]
        else:
            alt = 0.0
        
        # Rough estimate: Re ≈ 1e6 per ft of characteristic length at sea level
        char_length = self.state.get('options_cbarr', 10.0) or 10.0
        
        # Simple atmospheric correction
        reynolds = 1e6 * char_length * mach * np.exp(-alt / 30000.0)
        
        return reynolds


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

