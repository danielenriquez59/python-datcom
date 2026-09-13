"""AREA2 parity with its body-shadow FORTRAN algorithm and analytic moments."""

import math

import numpy as np
import pytest

from pydatcom.utils.math_utils import area2


def _legacy_area2(x, y, inum):
    """Literal scalar source equations; fixtures avoid singular median slopes."""
    a = math.hypot(x[1] - x[0], y[1] - y[0])
    b = math.hypot(x[2] - x[1], y[2] - y[1])
    c = math.hypot(x[0] - x[2], y[0] - y[2])
    s = (a + b + c) / 2
    area = math.sqrt(s * (s - a) * (s - b) * (s - c))
    x1, y1 = (x[0] + x[1]) / 2, (y[0] + y[1]) / 2
    x2, y2 = (x[0] + x[2]) / 2, (y[0] + y[2]) / 2
    dnum1 = (y[2] - y1) / (x[2] - x1)
    dnum2 = (y2 - y[1]) / (x2 - x[1])
    anum1, anum2 = dnum1 * x1 - y1, dnum2 * x[1] - y[1]
    xbar = (anum1 - anum2) / (dnum1 - dnum2)
    ybar = (xbar - x1) * dnum1 + y1
    ax, ay = area * xbar, area * ybar
    if inum != 3:
        a2 = math.hypot(x[3] - x[2], y[3] - y[2])
        b2 = math.hypot(x[0] - x[3], y[0] - y[3])
        c2 = math.hypot(x[2] - x[0], y[2] - y[0])
        s2 = (a2 + b2 + c2) / 2
        area22 = math.sqrt(s2 * (s2 - a2) * (s2 - b2) * (s2 - c2))
        ybar2 = y[2] / 3
        xbar2 = ((ybar2 - y[0]) * (x[2] - (x[0] + x[3]) / 2)
                 / (y[2] - y[0]) + (x[0] + x[3]) / 2)
        area += area22
        ax += area22 * xbar2
        ay += area22 * ybar2
    return area, ax, ay


@pytest.mark.parametrize("x,y,inum", [
    ([1., 2., 6.], [2., 7., 1.], 3),
    ([0., 1., 5., 4.], [0., 2., 3., 0.], 2),
    ([0., 1., 5., 4.], [0., -2., -3., 0.], 2),
    ([3., 5., 8., 7.], [0., 1., 4., 0.], 2),
])
def test_area2_matches_legacy_equations(x, y, inum):
    assert area2(x, y, inum) == pytest.approx(_legacy_area2(x, y, inum))


def test_triangle_returns_first_moments_and_ignores_fourth_vertex():
    # Triangle area=6, centroid=(7/3, 3); fourth entry is unused FORTRAN storage.
    assert area2([1, 5, 1, np.nan], [2, 2, 5, np.nan], 3) == pytest.approx((6, 14, 18))


@pytest.mark.parametrize("height", [3., -3.])
def test_body_rectangle_and_lower_profile(height):
    assert area2([0, 0, 4, 4], [0, height, height, 0], 2) == pytest.approx(
        (12, 24, 6 * height))


def test_vertical_median_and_collapsed_triangle_are_well_defined():
    assert area2([0, 2, 1], [0, 0, 3], 3) == pytest.approx((3, 3, 3))
    assert area2([1, 2, 3], [2, 4, 6], 3) == (0, 0, 0)
    assert area2([0, 1, 4, 3], [0, 0, 0, 0], 2) == (0, 0, 0)


def test_body_shadow_domain_is_explicit():
    with pytest.raises(ValueError, match="inum"):
        area2([0, 0, 1, 1], [0, 1, 1, 0], 4)
    with pytest.raises(ValueError, match="quadrilateral requires"):
        area2([0, 0, 1, 1], [2, 3, 3, 2], 2)
