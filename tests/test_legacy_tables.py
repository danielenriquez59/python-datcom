"""Source-contract checks for the legacy interpolation family."""

import numpy as np
import pytest

from pydatcom.utils.legacy_tables import (glook, switch, tlin1x, tlinex,
                                          tlinex_flat)


@pytest.mark.parametrize("x", [[1., 2., 4.], [4., 2., 1.]])
def test_glook_interpolation_snapping_and_clamping(x):
    index, fraction = glook(x, 3.)
    assert x[index-1] + fraction*(x[index]-x[index-1]) == 3.
    assert glook(x, 2.001) == (x.index(2.), None)
    assert glook(x, -10.) == (x.index(min(x)), None)
    assert glook(x, 10.) == (x.index(max(x)), None)


def test_glook_small_denominator_rule_and_first_match():
    # DMG becomes 1 for a tiny query: tolerance becomes absolute 0.001.
    assert glook([0., .01], .00005) == (0, None)
    # Zero query instead takes DMG from each candidate coordinate.
    assert glook([-.01, .01], 0.) == (1, .5)
    assert glook([.00005, .01], 0.) == (0, None)
    assert glook([1., 1.0005, 2.], 1.0005) == (0, None)
    assert glook([3.], 4.) == (0, None)


@pytest.mark.parametrize("x", [[1., 2., 4.], [4., 2., 1.]])
@pytest.mark.parametrize("mode", [-1, 0, 1, 2])
def test_switch_table_end_modes(x, mode):
    direction = np.sign(x[-1]-x[0])
    for first in [True, False]:
        query = x[0]-direction if first else x[-1]+direction
        flags = switch(x, query, mode if first else -1, -1 if first else mode)
        assert flags.before_first == first
        assert flags.after_last == (not first)
        assert flags.message_requested == (mode >= 0)
        assert flags.extrapolate == (mode > 0)
        assert flags.use_extrapolation == flags.extrapolate
        assert flags.ascending == (direction > 0)
        assert not flags.no_interpolation
    assert not any(switch(x, 2., mode, mode)[:5])


@pytest.mark.parametrize("reverse", [False, True])
def test_tlin1x_interpolates_linearly_and_extrapolates_by_mode(reverse):
    x = np.array([1., 2., 4., 7.])
    y = x*x + 2*x - 3
    if reverse:
        x, y = x[::-1], y[::-1]
    assert tlin1x(x, y, 3.) == 13.
    assert tlin1x(x, y, 2.001) == 5.
    for query in [-1., 9.]:
        index = np.argmin(abs(x-query))
        assert tlin1x(x, y, query, -1, -1) == y[index]
        assert tlin1x(x, y, query, 0, 0) == y[index]
        assert tlin1x(x, y, query, 2, 2) == pytest.approx(query*query+2*query-3)
        neighbor = 1 if index == 0 else len(x)-2
        expected = y[index] + (query-x[index])*(y[neighbor]-y[index])/(x[neighbor]-x[index])
        assert tlin1x(x, y, query, 1, 1) == pytest.approx(expected)


def test_tlin1x_small_tables_and_extrapolation_precedes_snapping():
    assert tlin1x([1.], [5.], 1., 2, 2) == 5.
    assert tlin1x([1.], [5.], 10.) == 5.
    with pytest.raises(ValueError, match="at least two"):
        tlin1x([1.], [5.], 10., upper=1)
    assert tlin1x([1., 2.], [2., 5.], 3., upper=2) == 8.
    assert tlin1x([1., 2.], [2., 5.], 2.0001, upper=1) == pytest.approx(5.0003)
    assert tlin1x([1., 2.], [2., 5.], 2.0001, upper=0) == 5.


@pytest.mark.parametrize("reverse1", [False, True])
@pytest.mark.parametrize("reverse2", [False, True])
def test_tlinex_orientation_and_bilinear_interpolation(reverse1, reverse2):
    x1 = np.array([1., 3., 5.])
    x2 = np.array([-2., 0., 4., 8.])
    if reverse1:
        x1 = x1[::-1]
    if reverse2:
        x2 = x2[::-1]
    surface = lambda a, b: 7. + 2.*a - 3.*b + 4.*a*b
    y = surface(x1[None, :], x2[:, None])
    assert tlinex(x1, x2, y, 2., 1.) == pytest.approx(surface(2., 1.))
    # Linear extrapolation on both axes retains the bilinear cross term.
    assert tlinex(x1, x2, y, -1., 10., 1, 1, 1, 1) == pytest.approx(surface(-1., 10.))
    # Clamp each axis independently.
    assert tlinex(x1, x2, y, -1., 10.) == pytest.approx(surface(1., 8.))
    # GLOOK applies its snapping separately on both axes.
    assert tlinex(x1, x2, y, 3.001, 4.001) == pytest.approx(surface(3., 4.))


@pytest.mark.parametrize("reverse", [False, True])
def test_tlinex_quadratic_extrapolation_on_both_ends(reverse):
    x1 = np.array([1., 2., 4.])
    x2 = np.array([-2., 0., 3.])
    if reverse:
        x1, x2 = x1[::-1], x2[::-1]
    surface = lambda a, b: a*a + a*b + 2*b*b + a*a*b*b
    y = surface(x1[None, :], x2[:, None])
    for a in [-3., 7.]:
        for b in [-4., 6.]:
            assert tlinex(x1, x2, y, a, b, 2, 2, 2, 2) == pytest.approx(surface(a, b))


def test_tlinex_mixed_modes_and_singletons():
    x1, x2 = [1., 2., 4.], [-2., 0., 3.]
    y = np.array(x1)[None, :]**2 + np.array(x2)[:, None]**2
    # First end X1 quadratic; last end X2 clamp.
    assert tlinex(x1, x2, y, -1., 5., lower1=2) == 10.
    # Last end X1 clamp; first end X2 linear.
    assert tlinex(x1, x2, y, 8., -3., lower2=1) == 22.
    assert tlinex([1.], [2.], [[7.]], 4., 5.) == 7.


@pytest.mark.parametrize("x", [[], [1., 1.], [1., 3., 2.], [1., np.nan]])
def test_invalid_grids_are_explicit_errors(x):
    with pytest.raises(ValueError):
        glook(x, 0.)


def test_tlinex_flat_reads_the_source_data_statement_order():
    # Each run of len(x2) values is one X1 column: Y(j,i) = flat[i*NX2+j].
    x1, x2 = [1., 2., 4.], [-2., 0.]
    flat = [10., 11., 20., 21., 40., 41.]
    y = np.array(flat).reshape(3, 2).T
    for q1, q2 in [(1., -2.), (3., -1.), (4., 0.), (6., 1.)]:
        assert (tlinex_flat(x1, x2, flat, q1, q2, 1, 1, 1, 1) ==
                tlinex(x1, x2, y, q1, q2, 1, 1, 1, 1))
    assert tlinex_flat(x1, x2, flat, 2., 0.) == 21.
    with pytest.raises(ValueError, match="needs 6 values"):
        tlinex_flat(x1, x2, flat[:5], 2., 0.)
