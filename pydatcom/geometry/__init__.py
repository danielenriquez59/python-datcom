"""Geometry calculations for PyDATCOM - airfoils, bodies, wings, tails."""

from pydatcom.geometry.airfoil import NACAGenerator, generate_naca_airfoil, AirfoilCoordinates
from pydatcom.geometry.body import BodyGeometry, calculate_body_geometry
from pydatcom.geometry.wing import (
    WingGeometry, calculate_wing_geometry, calculate_straight_exposed_geometry,
)
from pydatcom.geometry.tail import (
    TailGeometry, calculate_horizontal_tail, calculate_vertical_tail
)

__all__ = [
    'NACAGenerator', 'generate_naca_airfoil', 'AirfoilCoordinates',
    'BodyGeometry', 'calculate_body_geometry',
    'WingGeometry', 'calculate_wing_geometry', 'calculate_straight_exposed_geometry',
    'TailGeometry', 'calculate_horizontal_tail', 'calculate_vertical_tail',
]

from pydatcom.geometry.cslope import calculate_cslope

__all__ = list(__all__) + ['calculate_cslope']
