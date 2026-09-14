"""
Regression tests for the TLIN3X and TLIN4X translations.

Source of truth: datcom-legacy/datcom_2000/tlin3x.f, tlin4x.f.

TLIN3X's effective Y layout was established by compiling and running the
legacy routine, since its declaration is dead code. tlin4x.f corroborates
that finding in writing: its header states the structure is
Y(NX2,NX1,NX3,NX4).
"""

import numpy as np
import pytest

from pydatcom.utils.legacy_tables import tlin1x, tlinex, tlin3x, tlin4x

_X1 = [0., 1., 2.]
_X2 = [0., 10., 20., 30.]
_X3 = [0., 100.]
_X4 = [0., 1000., 2000.]


def _separable(a, b, c, d=0.0):
    """A multilinear function, which multilinear interpolation is exact on."""
    return 1.0 + a + b + c + d


def _cube():
    y = np.zeros((len(_X2), len(_X1), len(_X3)))
    for i, b in enumerate(_X2):
        for j, a in enumerate(_X1):
            for k, c in enumerate(_X3):
                y[i, j, k] = _separable(a, b, c)
    return y


def _hypercube():
    y = np.zeros((len(_X2), len(_X1), len(_X3), len(_X4)))
    for i, b in enumerate(_X2):
        for j, a in enumerate(_X1):
            for k, c in enumerate(_X3):
                for l, d in enumerate(_X4):
                    y[i, j, k, l] = _separable(a, b, c, d)
    return y


# --------------------------------------------------------------------------
# Layout
# --------------------------------------------------------------------------

def test_tlin3x_expects_x2_fastest():
    """Y is (len(x2), len(x1), len(x3)), not the declared (NX1,NX2,NX3)."""
    y = _cube()
    assert y.shape == (len(_X2), len(_X1), len(_X3))
    assert tlin3x(_X1, _X2, _X3, y, 1.0, 10.0, 0.0) == pytest.approx(
        _separable(1.0, 10.0, 0.0))


def test_tlin4x_expects_the_header_layout():
    """tlin4x.f states outright: Y(NX2,NX1,NX3,NX4)."""
    y = _hypercube()
    assert y.shape == (len(_X2), len(_X1), len(_X3), len(_X4))


def test_wrong_shapes_are_rejected():
    y = _hypercube()
    with pytest.raises(ValueError, match="shape"):
        tlin4x(_X1, _X2, _X3, _X4, y[:, :, :, :2], 0., 0., 0., 0.)
    with pytest.raises(ValueError, match="shape"):
        tlin3x(_X1, _X2, _X3, _cube()[:, :, :1], 0., 0., 0.)


# --------------------------------------------------------------------------
# Exactness on a multilinear function
# --------------------------------------------------------------------------

@pytest.mark.parametrize("a,b,c", [
    (0.0, 0.0, 0.0), (2.0, 30.0, 100.0), (0.5, 5.0, 50.0), (1.5, 25.0, 0.0),
])
def test_tlin3x_is_exact_on_a_separable_function(a, b, c):
    assert tlin3x(_X1, _X2, _X3, _cube(), a, b, c) == pytest.approx(
        _separable(a, b, c))


@pytest.mark.parametrize("a,b,c,d", [
    (0.0, 0.0, 0.0, 0.0), (2.0, 30.0, 100.0, 2000.0),
    (0.5, 5.0, 50.0, 500.0), (1.5, 25.0, 0.0, 1500.0),
])
def test_tlin4x_is_exact_on_a_separable_function(a, b, c, d):
    assert tlin4x(_X1, _X2, _X3, _X4, _hypercube(), a, b, c, d) == \
        pytest.approx(_separable(a, b, c, d))


# --------------------------------------------------------------------------
# Consistency down the family
# --------------------------------------------------------------------------

def test_tlin4x_at_an_x4_node_matches_tlin3x_on_that_slice():
    y = _hypercube()
    for index, node in enumerate(_X4):
        assert tlin4x(_X1, _X2, _X3, _X4, y, 0.5, 5.0, 50.0, node) == \
            pytest.approx(tlin3x(_X1, _X2, _X3, y[:, :, :, index],
                                 0.5, 5.0, 50.0))


def test_tlin3x_at_an_x3_node_matches_tlinex_on_that_slice():
    y = _cube()
    for index, node in enumerate(_X3):
        assert tlin3x(_X1, _X2, _X3, y, 0.5, 5.0, node) == pytest.approx(
            tlinex(_X1, _X2, y[:, :, index], 0.5, 5.0))


def test_a_degenerate_third_axis_reduces_to_the_two_variable_case():
    """With one X3 station TLIN3X must return the TLINEX result."""
    y = _cube()[:, :, :1]
    assert tlin3x(_X1, _X2, [0.0], y, 1.2, 12.0, 0.0) == pytest.approx(
        tlinex(_X1, _X2, y[:, :, 0], 1.2, 12.0))


def test_interx_three_variable_branch_still_matches_the_fortran_probe():
    """The refactor routed INTERX through tlin3x; the probe data must hold."""
    from pydatcom.utils.legacy_interp import interx
    table = [[10.0, 20.0, 30.0], [1.0, 2.0], [100.0, 200.0]]
    dep = [float(i + 1) for i in range(12)]
    # Values produced by the compiled legacy TLIN3X; see
    # tests/test_interx_translation.py for provenance.
    assert interx(3, table, [10.0, 2.0, 100.0], [3, 2, 2], dep) == \
        pytest.approx(4.0)
    assert interx(3, table, [30.0, 2.0, 200.0], [3, 2, 2], dep) == \
        pytest.approx(12.0)
