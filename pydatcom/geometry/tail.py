"""
Tail geometry calculations for PyDATCOM.

Horizontal and vertical tail geometry is handled as a specialized
lifting surface, using the wing geometry framework.

Reference: datcom.f COMMON /HTI/, /IHT/, /VTI/, /IVT/ blocks
"""

from typing import Dict
from pydatcom.geometry.wing import TailGeometry, calculate_tail_geometry

# Re-export for convenience
__all__ = ['TailGeometry', 'calculate_tail_geometry', 
           'calculate_horizontal_tail', 'calculate_vertical_tail']


def calculate_horizontal_tail(state: Dict) -> Dict[str, float]:
    """
    Calculate horizontal tail geometry.
    
    Args:
        state: State dictionary with htail_* parameters
        
    Returns:
        Dictionary with computed H-tail properties
    """
    return calculate_tail_geometry(state, tail_type='htail')


def calculate_vertical_tail(state: Dict) -> Dict[str, float]:
    """
    Calculate vertical tail geometry.
    
    Args:
        state: State dictionary with vtail_* parameters
        
    Returns:
        Dictionary with computed V-tail properties
    """
    return calculate_tail_geometry(state, tail_type='vtail')

