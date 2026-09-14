"""
Regression tests for the TRANAC translation.

Source of truth: datcom-legacy/datcom_2000/tranac.f.

TRANAC is a sibling of TRANF: same angular-average interior slopes and
per-interval cubic, but with an explicit coefficient solve and the DELY
offset mechanism. These tests cover the properties the method guarantees
and the DELY branch that distinguishes it.
"""

import numpy as np
import pytest

from pydatcom.utils.tranac import tranac, _DELY_POINTS, _DELY_INDEX


def _cubic(t):
    return 2.0 * t**3 - 3.0 * t**2 + t + 5.0


def _cubic_slope(t):
    return 6.0 * t**2 - 6.0 * t + 1.0


def _grid(points=7):
    x = np.linspace(0.0, 2.0, points)
    return x, _cubic(x)


# --------------------------------------------------------------------------
# Interpolation properties
# --------------------------------------------------------------------------

@pytest.mark.parametrize("index", range(7))
def test_every_node_is_reproduced_exactly(index):
    """The cubic is fitted through the interval endpoints, so nodes are exact."""
    x, y = _grid()
    result = tranac(x, y, _cubic_slope(x[0]), _cubic_slope(x[-1]),
                    float(x[index]))
    assert result['value'] == pytest.approx(y[index], abs=1e-12)


def test_the_curve_is_continuous_across_interior_nodes():
    x = np.linspace(0.0, 2.0, 7)
    y = np.array([1., 1.4, 1.8, 2.0, 1.9, 1.5, 0.8])
    for node in (2, 3, 4):
        low = tranac(x, y, 0.5, -1.0, float(x[node]) - 1e-9)['value']
        high = tranac(x, y, 0.5, -1.0, float(x[node]) + 1e-9)['value']
        assert low == pytest.approx(high, abs=1e-7)


def test_interior_points_do_not_reproduce_an_arbitrary_cubic():
    """Interior slopes are angular averages, not the true derivatives.

    TRANAC is shape-preserving rather than an exact cubic fit, so it matches
    at the nodes but departs between them. Asserting exactness there would
    be asserting against the method.
    """
    x, y = _grid()
    result = tranac(x, y, _cubic_slope(x[0]), _cubic_slope(x[-1]), 0.15)
    assert result['value'] != pytest.approx(_cubic(0.15), abs=1e-6)
    assert abs(result['value'] - _cubic(0.15)) < 0.05


def test_end_slopes_are_taken_from_the_caller():
    """The first and last points use DYL and DYR rather than an average."""
    x, y = _grid()
    first = tranac(x, y, 3.5, -1.0, float(x[0]) + 1e-6)
    assert first['slopes'][0] == pytest.approx(3.5)
    last = tranac(x, y, 1.0, -2.5, float(x[-1]) - 1e-6)
    assert last['slopes'][1] == pytest.approx(-2.5)


def test_slope_reversal_flattens_the_interior_slope():
    """A sign change between adjacent secants zeroes the slope."""
    x = np.array([0., 1., 2., 3.])
    y = np.array([0., 1., 0., 1.])
    result = tranac(x, y, 1.0, 1.0, 1.5)
    assert result['slopes'] == (0.0, 0.0)


# --------------------------------------------------------------------------
# The DELY mechanism
# --------------------------------------------------------------------------

def _eight():
    x = np.linspace(0.0, 7.0, _DELY_POINTS)
    y = np.array([1., 1.3, 1.6, 1.9, 2.1, 2.0, 1.7, 1.2])
    return x, y


def test_dely_shifts_the_fifth_point_of_an_eight_point_table():
    x, y = _eight()
    node = float(x[_DELY_INDEX - 1])
    base = tranac(x, y, 0.3, -0.5, node)['value']
    for offset in (0.25, -0.25):
        shifted = tranac(x, y, 0.3, -0.5, node, dely=offset)['value']
        assert shifted == pytest.approx(base + offset)


def test_dely_is_ignored_on_other_table_sizes():
    """The source gates the offset on NPT == 8."""
    x = np.linspace(0.0, 6.0, 7)
    y = np.array([1., 1.3, 1.6, 1.9, 2.1, 2.0, 1.7])
    plain = tranac(x, y, 0.3, -0.5, 4.0)
    offset = tranac(x, y, 0.3, -0.5, 4.0, dely=0.5)
    assert offset['value'] == pytest.approx(plain['value'])
    assert not offset['dely_applied']


def test_negative_dely_flattens_the_fifth_point_slope():
    """The source zeroes YP there independently of the offset itself.

    A monotonic table is used so the slope-reversal rule cannot flatten it
    instead, isolating the DELY branch.
    """
    x = np.linspace(0.0, 7.0, _DELY_POINTS)
    y = np.array([1.0, 1.3, 1.6, 1.9, 2.2, 2.5, 2.8, 3.1])
    node = _DELY_INDEX - 1
    rising = tranac(x, y, 0.3, 0.3, float(x[node]) + 0.5, dely=0.0)
    flattened = tranac(x, y, 0.3, 0.3, float(x[node]) + 0.5, dely=-0.1)
    assert rising['slopes'][0] > 0.0
    assert flattened['slopes'][0] == 0.0


def test_positive_dely_does_not_flatten_the_slope():
    x = np.linspace(0.0, 7.0, _DELY_POINTS)
    y = np.array([1.0, 1.3, 1.6, 1.9, 2.2, 2.5, 2.8, 3.1])
    node = _DELY_INDEX - 1
    result = tranac(x, y, 0.3, 0.3, float(x[node]) + 0.5, dely=0.1)
    assert result['slopes'][0] > 0.0


def test_the_callers_array_is_not_mutated():
    """The source adds DELY to Y in place and subtracts it again."""
    x, y = _eight()
    original = y.copy()
    tranac(x, y, 0.3, -0.5, 4.0, dely=0.5)
    assert np.array_equal(y, original)


# --------------------------------------------------------------------------
# Validation
# --------------------------------------------------------------------------

def test_rejects_bad_input():
    x, y = _grid()
    with pytest.raises(ValueError):
        tranac(x, y[:3], 1.0, 1.0, 0.5)
    with pytest.raises(ValueError):
        tranac([0.0], [1.0], 1.0, 1.0, 0.5)
    with pytest.raises(ValueError, match="increasing"):
        tranac([0.0, 2.0, 1.0], [1.0, 2.0, 3.0], 1.0, 1.0, 0.5)
