"""Analytic checks of the distinct QUAD/TRAPZ/TBFUNX legacy contracts."""

import numpy as np
import pytest

from pydatcom.utils.legacy_numeric import quad, trapz, tbfunx


def test_quad_values_and_derivatives_with_shifted_origin():
    x = np.array([0., 1., 3.])
    y = 2*x*x - 3*x + 4
    for q in [-2., .5, 4.]:
        for origin in [0., 1e9]:
            np.testing.assert_allclose(quad(x+origin, y, q+origin), 2*q*q-3*q+4)
            np.testing.assert_allclose(quad(x+origin, y, q+origin, True), 4*q-3)


def test_trapz_integral_modes_and_frustum_volume():
    np.testing.assert_allclose(trapz([0, 1, 3], [0, 1, 3], 0), [0, .5, 4.5])
    np.testing.assert_allclose(trapz([0, 1, 3], [0, 1, 3], 1), [4.5])
    np.testing.assert_allclose(trapz([0, 1, 3], [0, 1, 3], 2), [0, .5, 4.5])
    # A cone is integrated exactly from its radii, not by trapezoids of pi*r^2.
    np.testing.assert_allclose(trapz([0, 2], [0, 6], -1), [8*np.pi])
    np.testing.assert_allclose(trapz([2, 2], [0, 6], -2), [24*np.pi])
    np.testing.assert_allclose(trapz([2, 0], [6, 0], -1), [-8*np.pi])


def test_tbfunx_values_are_linear_but_derivatives_quadratic():
    x = np.array([0., 1., 2., 4.])
    y = x*x
    np.testing.assert_allclose(tbfunx(x, y, 1.5), [2.5, 3.])
    for q, expected in [(-1, [0, -2]), (5, [16, 10])]:
        np.testing.assert_allclose(tbfunx(x, y, q), expected)
    np.testing.assert_allclose(tbfunx(x, y, -1, 1, 1), [-1, 1])
    np.testing.assert_allclose(tbfunx(x, y, 5, 1, 1), [22, 6])
    np.testing.assert_allclose(tbfunx(x, y, -1, 2, 2), [1, -2])
    np.testing.assert_allclose(tbfunx(x, y, 5, 2, 2), [25, 10])
    np.testing.assert_allclose(tbfunx([0, 2], [1, 5], -2), [-3, 2])
    assert tbfunx([0], [7], 10) == (7., 0.)
    # Source's LE==NP3 condition is true at BOTH ends when NP=3.
    np.testing.assert_allclose(tbfunx([0, 1, 2], [0, 1, 4], -1, 1, 1), [-5, 3])


def test_tbfunx_unordered_follows_the_source_search():
    """With ordered=False the interior search takes the last point at or
    below the query, and below XA(2) interpolates on XA(L-1), XA(L).

    The table turns over, as a lift curve past the stall does.  Only the
    end points 0.0 and 0.9 decide that a query is interior.
    """
    x = [0.0, 0.8, 1.1, 0.4, 0.9]
    y = [0.0, 8.0, 12.0, 16.0, 20.0]
    # 0.85: the last point at or below it is XA(4)=0.4 (L=4), and it is not
    # below XA(2), so the source uses XX(2), XX(3) = XA(4), XA(5).
    value, _ = tbfunx(x, y, 0.85, 1, 1, ordered=False)
    assert value == pytest.approx(16.0 + 4.0 * (0.85 - 0.4) / (0.9 - 0.4))
    # 0.5: still L=4, but below XA(2)=0.8, so the override interpolates on
    # XX(1), XX(2) = XA(3), XA(4), not on the table's first pair.
    value, _ = tbfunx(x, y, 0.5, 1, 1, ordered=False)
    assert value == pytest.approx(12.0 + 4.0 * (0.5 - 1.1) / (0.4 - 1.1))
    with pytest.raises(ValueError):
        tbfunx(x, y, 0.5, 1, 1)


def test_tbfunx_unordered_agrees_on_ordered_tables():
    x = np.linspace(-4.0, 20.0, 9)
    y = np.sin(x / 7.0)
    for q in np.linspace(-6.0, 22.0, 57):
        assert tbfunx(x, y, q, 1, 2, ordered=False) == tbfunx(x, y, q, 1, 2)
