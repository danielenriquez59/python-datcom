"""
Body geometry calculations for PyDATCOM.

This module handles axisymmetric and asymmetric body geometry,
including area calculations, volume, and centroid computations.

Reference: datcom.f lines 1336 (BDAREA), 1805 (BODOPT), 2248 (BODYJM), 2326 (BODYRT)
"""

import numpy as np
from typing import Dict, Tuple, Optional, List
import logging

from pydatcom.utils.constants import UNUSED

logger = logging.getLogger(__name__)


class BodyGeometry:
    """
    Body geometry calculations for DATCOM analysis.
    
    Handles both axisymmetric and asymmetric (cambered) bodies.
    Reference: datcom.f BODOPT (line 1805), BODYJM (line 2248), BODYRT (line 2326)
    """
    
    def __init__(self, state: Dict):
        """
        Initialize body geometry from state dictionary.
        
        Args:
            state: State dictionary with body_* parameters
        """
        self.state = state
        self.nx = int(state.get('body_nx', 0))
        
        # Body stations
        self.x = np.array(state.get('body_x', []))
        self.s = np.array(state.get('body_s', []))  # Cross-sectional area
        self.p = np.array(state.get('body_p', []))  # Perimeter
        self.r = np.array(state.get('body_r', []))  # Half-width
        self.zu = np.array(state.get('body_zu', []))  # Upper Z coordinate
        self.zl = np.array(state.get('body_zl', []))  # Lower Z coordinate
        
        # Body parameters
        self.bnose = state.get('body_bnose', 1.0)  # Nose type (1=conical, 2=ogive)
        self.btail = state.get('body_btail', 1.0)  # Tail type (1=conical, 2=ogive)
        self.bln = state.get('body_bln', None)     # Nose length
        self.bla = state.get('body_bla', None)     # Afterbody length
        self.ds = state.get('body_ds', 0.0)        # Nose bluntness diameter
        self.itype = state.get('body_itype', 2)    # Type (1=straight, 2=swept no AR, 3=swept AR)
        self.method = state.get('body_method', 1)  # Method (1=existing, 2=Jorgensen)
        
        # Computed properties
        self.volume = None
        self.centroid = None
        self.max_area = None
        self.max_area_location = None
        self.length = None
        self.fineness_ratio = None
    
    def calculate_properties(self) -> Dict[str, float]:
        """
        Calculate basic body geometric properties.
        
        Computes:
        - Body volume (integration of cross-sectional area)
        - Centroid location
        - Maximum area and its location
        - Body length
        - Fineness ratio
        
        Returns:
            Dictionary with computed properties
        """
        if self.nx < 2:
            logger.warning("Insufficient body stations (nx < 2)")
            return {}
        
        # Body length
        self.length = self.x[-1] - self.x[0]
        
        # Find maximum area and location
        max_idx = np.argmax(self.s)
        self.max_area = self.s[max_idx]
        self.max_area_location = self.x[max_idx]
        
        # Calculate equivalent radius at max area
        req_max = np.sqrt(self.max_area / np.pi) if self.max_area > 0 else 0.0
        
        # Fineness ratio (length / max diameter)
        if req_max > 0:
            self.fineness_ratio = self.length / (2.0 * req_max)
        else:
            self.fineness_ratio = 0.0
        
        # Calculate volume using trapezoidal integration
        if len(self.x) > 1 and len(self.s) > 1:
            # Use NumPy trapz directly
            try:
                self.volume = np.trapezoid(self.s, self.x)
            except AttributeError:
                self.volume = np.trapz(self.s, self.x)
        else:
            self.volume = 0.0
        
        # Calculate centroid (first moment / volume)
        if self.volume > 1e-10:
            x_times_area = self.x * self.s
            try:
                first_moment = np.trapezoid(x_times_area, self.x)
            except AttributeError:
                first_moment = np.trapz(x_times_area, self.x)
            self.centroid = first_moment / self.volume
        else:
            self.centroid = 0.0
        
        return {
            'length': self.length,
            'volume': self.volume,
            'centroid': self.centroid,
            'max_area': self.max_area,
            'max_area_location': self.max_area_location,
            'fineness_ratio': self.fineness_ratio,
        }
    
    def calculate_equivalent_body(self) -> Dict[str, np.ndarray]:
        """
        Calculate equivalent axisymmetric body parameters.
        
        For asymmetric bodies (with ZU, ZL), calculates equivalent
        axisymmetric representation.
        
        Reference: datcom.f BODYJM line 2281-2290
        
        Returns:
            Dictionary with equivalent body parameters
        """
        if self.nx < 2:
            return {}
        
        equiv_radius = np.sqrt(self.s / np.pi)
        radius_times_x = equiv_radius * self.x

        try:
            planform_area = 2.0 * np.trapezoid(equiv_radius, self.x)
            volume_integrand = np.trapezoid(equiv_radius, self.x)
            centroid_x = (2.0 * np.trapezoid(radius_times_x, self.x) /
                          planform_area if planform_area > 0 else 0.0)
        except AttributeError:
            planform_area = 2.0 * np.trapz(equiv_radius, self.x)
            volume_integrand = np.trapz(equiv_radius, self.x)
            centroid_x = (2.0 * np.trapz(radius_times_x, self.x) /
                          planform_area if planform_area > 0 else 0.0)

        return {
            'equivalent_radius': equiv_radius,
            'rx': radius_times_x,
            'planform_area': planform_area,
            'volume_integration': volume_integrand,
            'centroid': centroid_x,
        }
    
    def is_asymmetric(self) -> bool:
        """
        Check if body is asymmetric (has camber).
        
        Returns:
            True if ZU or ZL arrays are defined
        """
        return len(self.zu) > 0 or len(self.zl) > 0
    
    def calculate_cross_sectional_properties(self, station_idx: int) -> Dict[str, float]:
        """
        Calculate properties at a specific body station.
        
        Args:
            station_idx: Index of body station (0-based)
            
        Returns:
            Dictionary with station properties
        """
        if station_idx < 0 or station_idx >= self.nx:
            raise ValueError(f"Station index {station_idx} out of range [0, {self.nx-1}]")
        
        # Basic properties
        props = {
            'x': self.x[station_idx],
            'area': self.s[station_idx],
            'perimeter': self.p[station_idx] if station_idx < len(self.p) else 0.0,
            'half_width': self.r[station_idx] if station_idx < len(self.r) else 0.0,
        }
        
        # Calculate equivalent radius
        if props['area'] > 0:
            props['equivalent_radius'] = np.sqrt(props['area'] / np.pi)
        else:
            props['equivalent_radius'] = 0.0
        
        # Add asymmetric coordinates if available
        if station_idx < len(self.zu):
            props['z_upper'] = self.zu[station_idx]
        if station_idx < len(self.zl):
            props['z_lower'] = self.zl[station_idx]
        
        # Calculate centroid offset for asymmetric body
        if 'z_upper' in props and 'z_lower' in props:
            props['z_centroid'] = (props['z_upper'] + props['z_lower']) / 2.0
        
        return props
    
    def calculate_nose_properties(self) -> Dict[str, float]:
        """
        Calculate nose geometry properties.
        
        Uses BNOSE type (1=conical, 2=ogive) and BLN (nose length).
        
        Returns:
            Dictionary with nose properties
        """
        if self.bln is None or self.bln <= 0:
            return {'type': 'none'}
        
        nose_type = 'conical' if self.bnose == 1.0 else 'ogive'
        
        # Find base of nose (first station with significant area)
        base_idx = 0
        for station in range(self.nx):
            if self.s[station] > 0.01 * self.max_area:
                base_idx = station
                break
        
        base_area = self.s[base_idx] if base_idx < self.nx else 0.0
        base_radius = np.sqrt(base_area / np.pi) if base_area > 0 else 0.0
        
        # Nose fineness ratio
        nose_fineness = self.bln / (2.0 * base_radius) if base_radius > 0 else 0.0
        
        return {
            'type': nose_type,
            'length': self.bln,
            'bluntness_diameter': self.ds,
            'base_radius': base_radius,
            'base_area': base_area,
            'fineness_ratio': nose_fineness,
        }
    
    def calculate_tail_properties(self) -> Dict[str, float]:
        """
        Calculate tail/afterbody geometry properties.
        
        Uses BTAIL type (1=conical, 2=ogive) and BLA (afterbody length).
        
        Returns:
            Dictionary with tail properties
        """
        if self.bla is None or self.bla <= 0:
            return {'type': 'none'}
        
        tail_type = 'conical' if self.btail == 1.0 else 'ogive'
        
        # Base area at start of afterbody
        # Typically at or near maximum area
        base_area = self.max_area
        base_radius = np.sqrt(base_area / np.pi) if base_area > 0 else 0.0
        
        # Tail fineness ratio
        tail_fineness = self.bla / (2.0 * base_radius) if base_radius > 0 else 0.0
        
        return {
            'type': tail_type,
            'length': self.bla,
            'base_radius': base_radius,
            'base_area': base_area,
            'fineness_ratio': tail_fineness,
        }
    
    def to_state_dict(self) -> Dict[str, any]:
        """
        Export computed properties back to state dictionary format.
        
        Returns:
            Dictionary with body_* prefixed keys
        """
        props = self.calculate_properties()
        
        return {
            'body_length': props.get('length'),
            'body_volume': props.get('volume'),
            'body_centroid': props.get('centroid'),
            'body_max_area': props.get('max_area'),
            'body_max_area_x': props.get('max_area_location'),
            'body_fineness_ratio': props.get('fineness_ratio'),
        }


def calculate_body_geometry(state: Dict) -> Dict[str, float]:
    """
    Convenience function to calculate body geometry from state dict.
    
    Args:
        state: State dictionary with body parameters
        
    Returns:
        Dictionary with computed body properties
    """
    body = BodyGeometry(state)
    return body.calculate_properties()


def get_body_cross_section(state: Dict, x_location: float) -> Dict[str, float]:
    """
    Get body cross-sectional properties at specific X location.
    
    Interpolates between defined stations if necessary.
    
    Args:
        state: State dictionary with body parameters
        x_location: X coordinate for cross-section
        
    Returns:
        Dictionary with interpolated cross-section properties
    """
    body = BodyGeometry(state)
    
    if body.nx < 2:
        return {}
    
    # Find surrounding stations
    if x_location <= body.x[0]:
        return body.calculate_cross_sectional_properties(0)
    elif x_location >= body.x[-1]:
        return body.calculate_cross_sectional_properties(body.nx - 1)
    
    station_index = np.searchsorted(body.x, x_location)
    if station_index >= body.nx:
        station_index = body.nx - 1

    # For now, return nearest station
    # Full implementation would interpolate
    if station_index > 0:
        if (abs(body.x[station_index] - x_location) <
                abs(body.x[station_index - 1] - x_location)):
            return body.calculate_cross_sectional_properties(station_index)
        return body.calculate_cross_sectional_properties(station_index - 1)

    return body.calculate_cross_sectional_properties(station_index)

