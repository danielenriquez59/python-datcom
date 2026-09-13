"""Source-level checks for the legacy TRANF interpolation routine."""

import numpy as np
import pytest

from pydatcom.utils.legacy_numeric import tranf


def _hermite(x0, x1, y0, y1, slope0, slope1, query):
    width = x1 - x0
    t = (query - x0) / width
    return ((2*t**3 - 3*t**2 + 1)*y0 +
            (t**3 - 2*t**2 + t)*width*slope0 +
            (-2*t**3 + 3*t**2)*y1 +
            (t**3 - t**2)*width*slope1)


def test_tranf_uses_supplied_endpoint_slopes_and_extrapolates():
    x = np.array([0.0, 1.0, 2.0])
    y = np.array([1.0, 2.0, 4.0])
    interior = np.tan((np.arctan(1.0) + np.arctan(2.0)) / 2.0)
    assert tranf(x, y, 0.5, 3.0, -0.25) == pytest.approx(
        _hermite(0.0, 1.0, 1.0, 2.0, 0.5, interior, -0.25))
    assert tranf(x, y, 0.5, 3.0, 2.25) == pytest.approx(
        _hermite(1.0, 2.0, 2.0, 4.0, interior, 3.0, 2.25))


def test_tranf_preserves_knots_and_flattens_a_slope_reversal():
    x = np.array([0.0, 1.0, 2.0])
    y = np.array([1.0, 3.0, 2.0])
    assert tranf(x, y, 2.0, -1.0, 1.0) == pytest.approx(3.0)
    assert tranf(x, y, 2.0, -1.0, 1.5) == pytest.approx(
        _hermite(1.0, 2.0, 3.0, 2.0, 0.0, -1.0, 1.5))


def test_tranf_flattens_nonpositive_points_and_clamps_output():
    x = np.array([0.0, 1.0, 2.0])
    y = np.array([1.0, 0.0, 1.0])
    assert tranf(x, y, -10.0, 10.0, 1.0) == 0.0
    assert tranf([0.0, 1.0], [0.0, 0.0], -1.0, -1.0, 0.5) == 0.0


def test_tranf_long_table_zero_secant_rule():
    x = np.arange(11.0)
    y = np.array([1.0, 2.0, 2.0, 3.0, 4.0, 5.0,
                  6.0, 7.0, 8.0, 9.0, 10.0])
    expected = _hermite(0.0, 1.0, 1.0, 2.0, 1.0, 0.0, 0.5)
    assert tranf(x, y, 1.0, 1.0, 0.5) == pytest.approx(expected)


@pytest.mark.parametrize('x,y', [
    ([0.0], [1.0]),
    ([0.0, 0.0], [1.0, 2.0]),
    ([1.0, 0.0], [1.0, 2.0]),
])
def test_tranf_rejects_invalid_grids(x, y):
    with pytest.raises(ValueError):
        tranf(x, y, 0.0, 0.0, 0.5)
