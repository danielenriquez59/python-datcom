"""Geometry for PyDATCOM - airfoils, sections, bodies, wings and tails."""

from pydatcom.geometry.airfoil import (
    NACAGenerator, generate_naca_airfoil, AirfoilCoordinates,
)
from pydatcom.geometry.body import BodyGeometry, calculate_body_geometry
from pydatcom.geometry.wing import (
    WingGeometry, calculate_wing_geometry, calculate_straight_exposed_geometry,
)
from pydatcom.geometry.tail import (
    TailGeometry, calculate_horizontal_tail, calculate_vertical_tail,
)
from pydatcom.geometry.cslope import calculate_cslope
from pydatcom.geometry.section_params import calculate_dely, calculate_arclss
from pydatcom.geometry.ideal import calculate_ideal
from pydatcom.geometry.slope import calculate_slope

__all__ = [
    # Airfoil coordinate generation
    'NACAGenerator', 'generate_naca_airfoil', 'AirfoilCoordinates',
    # Body, wing and tail planforms
    'BodyGeometry', 'calculate_body_geometry',
    'WingGeometry', 'calculate_wing_geometry',
    'calculate_straight_exposed_geometry',
    'TailGeometry', 'calculate_horizontal_tail', 'calculate_vertical_tail',
    # Section parameters and the Weber section chain
    'calculate_cslope', 'calculate_dely', 'calculate_arclss',
    'calculate_ideal', 'calculate_slope',
]
