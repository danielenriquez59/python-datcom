"""Checks of source contracts and numerical invariants for translated routines."""

import numpy as np

from pydatcom.utils.interpolation import asmint
from pydatcom.utils.table_lookup import fig53a, fig60b, get_table_manager


def test_asmint_end_parabolas_and_join_slopes():
    x = np.array([0., 1., 2., 3., 4.])
    y = np.array([0., 1., 1.5, 1.8, 2.])
    # ASMINT labels 1020/1070: the left parabola uses the slope at X(2).
    slope = np.tan((np.arctan(1.) + np.arctan(.5)) / 2)
    query = np.array([-1., 0., .5, 1.])
    expected = (slope - 1) * query**2 + (2 - slope) * query
    np.testing.assert_allclose(asmint(x, y, query), expected, atol=1e-14)
    np.testing.assert_array_equal(asmint(x, y, x), y)
    for center in x[1:-1]:
        h = 1e-5
        values = asmint(x, y, np.array([center-h, center, center+h]))
        assert abs((values[1]-values[0])/h - (values[2]-values[1])/h) < 2e-5


def test_asmint_small_intervals_preserve_curve():
    x = np.array([0., 1., 2., 3., 4.])
    y = np.array([0., 1., 1.5, 1.8, 2.])
    query = np.array([-.5, .5, 1.5, 2.5, 3.5, 4.5])
    # Equal scaling of both axes preserves the secant angles.
    np.testing.assert_allclose(
        asmint(x * 1e-4, y * 1e-4, query * 1e-4) / 1e-4,
        asmint(x, y, query), rtol=1e-12, atol=1e-12,
    )


def test_fig53a_source_table_and_return_branches():
    # FIG53A DATA polynomials evaluated at LOG10(RV)=5. These values
    # distinguish the original coefficients and Z axis from the old fit.
    z = [0., 1., 2., 4., 10.]
    expected = [.86601125, .90500875, .92399, .9451375, .9597125]
    np.testing.assert_allclose([fig53a(1e5, zi) for zi in z], expected, atol=1e-12)
    for i in range(4):
        np.testing.assert_allclose(
            fig53a(1e5, (z[i] + z[i+1]) / 2),
            (expected[i] + expected[i+1]) / 2, atol=1e-12,
        )
    assert fig53a(0., 2.) == 0.
    assert fig53a(1., 2.) == 0.
    assert fig53a(10., 1.) == 0.  # Negative polynomial is clipped.
    assert fig53a(1e5, .0005) == fig53a(1e5, 0.)
    assert fig53a(1e5, 12.) == fig53a(1e5, 10.)
    # Legacy search also falls through to the last polynomial below Z=0.
    assert fig53a(1e5, -1.) == fig53a(1e5, 10.)


def test_fig60b_inverts_source_table_with_both_inputs():
    np.testing.assert_allclose(fig60b(2., .637), 1.)
    np.testing.assert_allclose(fig60b(2., (.637 + .732)/2), 1.1)
    np.testing.assert_allclose(fig60b(2.25, (.637 + .873)/2), 1.)
    assert fig60b(2., 10.) == 2.4
    assert fig60b(2., -1.) == 0.
    assert fig60b(30., 1.667) == 1.
    assert fig60b(2.0005, .637) == 1.  # GLOOK relative snapping.
    # Lower BETA extrapolation follows TLINEX mode 1, without snapping.
    np.testing.assert_allclose(fig60b(1.2, .047 + (.104-.047)*(-.05/.25)), .2)
    assert get_table_manager().lookup('FIG60B', 2., .637) == 1.
