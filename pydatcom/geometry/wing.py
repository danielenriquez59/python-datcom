"""
Wing geometry calculations for PyDATCOM.

Handles wing planform geometry including:
- Reference areas and chords
- Sweep angles
- Taper ratios
- Aspect ratio
- Mean aerodynamic chord

Reference: datcom.f WINGI common block and related calculations
"""

import numpy as np
from typing import Dict, Tuple, Optional
import logging

from pydatcom.utils.constants import PI, DEG, RAD, UNUSED

logger = logging.getLogger(__name__)


class WingGeometry:
    """
    Wing planform geometry calculations.
    
    Handles:
    - Wing planform parameters from WGPLNF namelist
    - Geometric calculations (area, MAC, aspect ratio, taper)
    - Sweep angle conversions
    - Dihedral effects
    
    Reference: COMMON /WINGI/ in datcom.f
    """
    
    def __init__(self, state: Dict, component: str = 'wing'):
        """
        Initialize wing geometry from state dictionary.
        
        Args:
            state: State dictionary
            component: Component prefix ('wing', 'htail', 'vtail', 'vfin')
        """
        self.state = state
        self.prefix = component
        
        # Planform parameters (WGPLNF namelist)
        self.chrdtp = state.get(f'{component}_chrdtp', None)  # Tip chord
        self.sspnop = state.get(f'{component}_sspnop', None)  # Semispan outboard panel
        self.sspne = state.get(f'{component}_sspne', None)    # Semispan exposed
        self.sspn = state.get(f'{component}_sspn', None)      # Semispan theoretical
        self.chrdbp = state.get(f'{component}_chrdbp', None)  # Chord at breakpoint
        self.chrdr = state.get(f'{component}_chrdr', None)    # Root chord
        self.savsi = state.get(f'{component}_savsi', 0.0)     # Inboard sweep angle (deg)
        self.savso = state.get(f'{component}_savso', 0.0)     # Outboard sweep angle (deg)
        self.chstat = state.get(f'{component}_chstat', 0.25)  # Sweep reference station
        self.twista = state.get(f'{component}_twista', 0.0)   # Twist angle (deg)
        self.dhdadi = state.get(f'{component}_dhdadi', 0.0)   # Inboard dihedral (deg)
        self.dhdado = state.get(f'{component}_dhdado', 0.0)   # Outboard dihedral (deg)
        self.ptype = state.get(f'{component}_type', 1.0)      # Planform type
        
        # Section characteristics (WGSCHR namelist)
        self.tovc = state.get(f'{component}_tovc', None)      # Thickness ratio
        self.xovc = state.get(f'{component}_xovc', None)      # Location of max thickness
        
        # Computed properties
        self.area = None
        self.span = None
        self.aspect_ratio = None
        self.taper_ratio = None
        self.mac = None          # Mean aerodynamic chord
        self.mac_location = None  # MAC X location
    
    def calculate_planform_properties(self) -> Dict[str, float]:
        """
        Calculate wing planform geometric properties.
        
        Computes:
        - Reference area
        - Span
        - Aspect ratio
        - Taper ratio
        - Mean aerodynamic chord (MAC)
        - MAC location
        
        Returns:
            Dictionary with computed properties
        """
        # Total span (double the semispan)
        if self.sspn is not None:
            self.span = 2.0 * self.sspn
        elif self.sspne is not None:
            self.span = 2.0 * self.sspne
        else:
            logger.warning(f"No span defined for {self.prefix}")
            self.span = 0.0
        
        # Calculate reference area
        if self.chrdr is not None and self.sspn is not None:
            # Simple trapezoidal wing
            if self.chrdtp is not None:
                # Straight tapered wing
                self.area = self.sspn * (self.chrdr + self.chrdtp)
                self.taper_ratio = self.chrdtp / self.chrdr if self.chrdr > 0 else 0.0
            else:
                # Need more info for complex planforms
                logger.warning("Complex planform - using simplified area calculation")
                self.area = self.chrdr * self.sspn * 1.5  # Rough estimate
                self.taper_ratio = 0.5  # Assume moderate taper
        else:
            self.area = 0.0
            self.taper_ratio = 0.0
        
        # Aspect ratio
        if self.area > 0:
            self.aspect_ratio = self.span**2 / self.area
        else:
            self.aspect_ratio = 0.0
        
        # Mean aerodynamic chord (MAC)
        if self.chrdr is not None and self.taper_ratio is not None:
            # Standard formula for trapezoidal wing
            lambda_ratio = self.taper_ratio
            self.mac = (2.0 / 3.0) * self.chrdr * (
                (1.0 + lambda_ratio + lambda_ratio**2) / (1.0 + lambda_ratio)
            )
            
            # MAC spanwise location
            if self.sspn and self.sspn > 0:
                y_mac = (self.sspn / 3.0) * (
                    (1.0 + 2.0 * lambda_ratio) / (1.0 + lambda_ratio)
                )
            else:
                y_mac = 0.0
            
            # MAC chordwise location (depends on sweep)
            # Simplified - full calculation requires sweep angle
            self.mac_location = 0.0  # From root LE
        else:
            self.mac = 0.0
            self.mac_location = 0.0
        
        return {
            'area': self.area,
            'span': self.span,
            'aspect_ratio': self.aspect_ratio,
            'taper_ratio': self.taper_ratio,
            'mac': self.mac,
            'mac_location': self.mac_location,
        }
    
    def calculate_sweep_at_station(self, x_c: float) -> float:
        """
        Calculate sweep angle at specific chord station.
        
        Args:
            x_c: Chord fraction (0 = LE, 0.25 = quarter-chord, 1.0 = TE)
            
        Returns:
            Sweep angle in degrees
        """
        if self.chstat is None:
            return 0.0
        
        # Get sweep at reference station
        sweep_ref = self.savsi if self.savsi != 0.0 else self.savso
        
        # For simple wing, use reference sweep
        # Full implementation would convert between different sweep reference lines
        return sweep_ref
    
    def calculate_panel_areas(self) -> Dict[str, float]:
        """
        Calculate individual panel areas for multi-panel wings.
        
        For cranked wings (TYPE=3.0) or double-delta (TYPE=2.0).
        
        Returns:
            Dictionary with panel areas
        """
        areas = {}
        
        if self.ptype == 1.0:
            # Straight tapered - single panel
            areas['total'] = self.area if self.area else 0.0
            areas['inboard'] = 0.0
            areas['outboard'] = areas['total']
        
        elif self.ptype == 2.0:
            # Double delta - two panels
            if self.chrdbp and self.sspnop and self.chrdr and self.chrdtp:
                # Inboard panel
                areas['inboard'] = (self.sspnop or 0) * (self.chrdr + self.chrdbp) / 2.0
                # Outboard panel  
                span_out = (self.sspn or 0) - (self.sspnop or 0)
                areas['outboard'] = span_out * (self.chrdbp + self.chrdtp) / 2.0
                areas['total'] = areas['inboard'] + areas['outboard']
            else:
                areas['total'] = 0.0
        
        elif self.ptype == 3.0:
            # Cranked wing
            logger.warning("Cranked wing area calculation simplified")
            areas['total'] = self.area if self.area else 0.0
        
        return areas
    
    def to_state_dict(self) -> Dict[str, any]:
        """
        Export computed properties back to state dictionary format.
        
        Returns:
            Dictionary with component_* prefixed keys
        """
        props = self.calculate_planform_properties()
        
        return {
            f'{self.prefix}_area': props.get('area'),
            f'{self.prefix}_span': props.get('span'),
            f'{self.prefix}_aspect_ratio': props.get('aspect_ratio'),
            f'{self.prefix}_taper_ratio': props.get('taper_ratio'),
            f'{self.prefix}_mac': props.get('mac'),
            f'{self.prefix}_mac_location': props.get('mac_location'),
        }


class TailGeometry(WingGeometry):
    """
    Tail geometry calculations (horizontal and vertical).
    
    Inherits from WingGeometry as tails are treated as lifting surfaces
    with similar planform calculations.
    """
    
    def __init__(self, state: Dict, tail_type: str = 'htail'):
        """
        Initialize tail geometry.
        
        Args:
            state: State dictionary
            tail_type: 'htail' for horizontal tail, 'vtail' for vertical tail
        """
        super().__init__(state, component=tail_type)
        self.tail_type = tail_type


def calculate_wing_geometry(state: Dict) -> Dict[str, float]:
    """
    Convenience function to calculate wing geometry from state dict.
    
    Args:
        state: State dictionary with wing parameters
        
    Returns:
        Dictionary with computed wing properties
    """
    wing = WingGeometry(state, component='wing')
    return wing.calculate_planform_properties()


def calculate_tail_geometry(state: Dict, tail_type: str = 'htail') -> Dict[str, float]:
    """
    Convenience function to calculate tail geometry from state dict.
    
    Args:
        state: State dictionary with tail parameters
        tail_type: 'htail' or 'vtail'
        
    Returns:
        Dictionary with computed tail properties
    """
    tail = TailGeometry(state, tail_type=tail_type)
    return tail.calculate_planform_properties()

