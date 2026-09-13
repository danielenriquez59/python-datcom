"""Numeric CORDSP branches; NACA cards and COMMON input updates are not ported."""

import numpy as np
import pytest

from pydatcom.geometry.airfoil import NACAGenerator


def section_at(stations, **parameters):
    generator = NACAGenerator()
    generator.x_stations = np.array(stations, dtype=float)
    return generator.supersonic_airfoil(**parameters)


def test_double_wedge_source_vertices_and_straight_segments():
    coords = section_at([0, .15, .3, .65, 1], thickness_ratio=.08,
                        max_thickness_location=.3)
    np.testing.assert_allclose(coords.yu, [0, .02, .04, .02, 0], atol=1e-16)
    np.testing.assert_allclose(coords.yl, -coords.yu)
    np.testing.assert_array_equal(coords.xu, coords.x)
    np.testing.assert_array_equal(coords.xl, coords.x)
    np.testing.assert_array_equal(coords.thickness, coords.yu)
    assert coords.thickness_ratio == .08
    assert coords.max_thickness_location == .3
    assert coords.sharpness_parameter == pytest.approx(1 / (.3 * .7))


def test_default_diamond_has_requested_full_thickness():
    coords = NACAGenerator(num_points=61).supersonic_airfoil()
    assert np.max(coords.yu - coords.yl) == pytest.approx(.05)
    assert coords.sharpness_parameter == 4


def test_biconvex_points_lie_on_source_circle():
    coords = section_at([0, .1, .25, .5, .75, .9, 1], shape='biconvex',
                        thickness_ratio=.2, max_thickness_location=.3)
    # Circle through (0,0), (.5,.1), (1,0): center (.5,-1.2), radius 1.3.
    np.testing.assert_allclose((coords.x - .5) ** 2 + (coords.yu + 1.2) ** 2,
                               1.3 ** 2, atol=1e-14)
    np.testing.assert_allclose(coords.yu[[0, 3, -1]], [0, .1, 0], atol=1e-15)
    np.testing.assert_allclose(coords.yl, -coords.yu)
    assert coords.max_thickness_location == .5
    assert coords.sharpness_parameter == pytest.approx(16 / 3)


def test_hexagonal_source_vertices_and_straight_segments():
    coords = section_at([0, .1, .2, .35, .5, .75, 1], shape='hexagonal',
                        thickness_ratio=.08, max_thickness_location=.2,
                        flat_length=.3)
    np.testing.assert_allclose(coords.yu, [0, .02, .04, .04, .04, .02, 0],
                               atol=1e-16)
    np.testing.assert_allclose(coords.yl, -coords.yu)
    assert coords.sharpness_parameter == pytest.approx(7)


def test_zero_flat_length_hexagon_reduces_to_double_wedge():
    generator = NACAGenerator()
    wedge = generator.supersonic_airfoil(.06, max_thickness_location=.4)
    hexagon = generator.supersonic_airfoil(.06, shape='hexagonal',
                                         max_thickness_location=.4, flat_length=0)
    np.testing.assert_allclose(hexagon.yu, wedge.yu)
    assert hexagon.sharpness_parameter == pytest.approx(wedge.sharpness_parameter)


@pytest.mark.parametrize('parameters', [
    {'shape': 'unknown'}, {'thickness_ratio': 0}, {'thickness_ratio': -.05},
    {'thickness_ratio': np.nan}, {'max_thickness_location': 0},
    {'max_thickness_location': 1}, {'max_thickness_location': np.inf},
    {'shape': 'hexagonal', 'flat_length': -.1},
    {'shape': 'hexagonal', 'flat_length': .5},
    {'shape': 'hexagonal', 'flat_length': np.nan},
])
def test_invalid_section_parameters_fail_explicitly(parameters):
    with pytest.raises(ValueError):
        NACAGenerator().supersonic_airfoil(**parameters)
