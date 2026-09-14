"""
Regression tests for the IDEAL translation.

Source of truth: datcom-legacy/datcom_2000/ideal.f.

IDEAL builds the Weber coefficients SLOPE consumes. Beyond the source
structure, the derived section parameters are checked against thin-airfoil
theory and against published NACA 2412 data, which is an independent
constraint rather than a restatement of the source.
"""

import numpy as np
import pytest

from pydatcom.geometry.ideal import (
    calculate_ideal, _weber_grid, _STATIONS,
)


def _naca(thickness=0.12, camber=0.02, position=0.4, points=61):
    """A NACA four-digit style section."""
    x = np.linspace(0.0, 1.0, points)
    t = thickness / 0.2 * (0.2969 * np.sqrt(np.clip(x, 0.0, None)) -
                           0.1260 * x - 0.3516 * x**2 +
                           0.2843 * x**3 - 0.1015 * x**4)
    if camber == 0.0:
        c = np.zeros_like(x)
    else:
        c = np.where(x < position,
                     camber / position**2 * (2 * position * x - x**2),
                     camber / (1 - position)**2 *
                     ((1 - 2 * position) + 2 * position * x - x**2))
    return x, t, c


def _run(**kwargs):
    x, t, c = _naca(**{k: v for k, v in kwargs.items()
                       if k in ('thickness', 'camber', 'position', 'points')})
    rest = {k: v for k, v in kwargs.items()
            if k not in ('thickness', 'camber', 'position', 'points')}
    return calculate_ideal(x, t, c, leading_edge_radius=0.0158, **rest)


# --------------------------------------------------------------------------
# The Weber grid
# --------------------------------------------------------------------------

def test_grid_is_thirty_two_cosine_spaced_stations():
    theta_nu, theta_mu, x = _weber_grid()
    assert len(theta_nu) == _STATIONS == 32
    assert len(theta_mu) == _STATIONS - 1
    assert x[-1] == pytest.approx(0.0)          # leading edge
    assert x[0] < 1.0 and x[0] > 0.99           # near the trailing edge
    assert np.all(np.diff(x) < 0.0)             # runs TE to LE


def test_leading_edge_parameter():
    """A0 = sqrt(2*RHO)."""
    assert _run()['a0'] == pytest.approx(np.sqrt(2.0 * 0.0158))


# --------------------------------------------------------------------------
# Weber coefficients
# --------------------------------------------------------------------------

def test_all_five_coefficient_sets_are_returned():
    result = _run()
    for key in ('st1', 'st2', 'st3', 'st4', 'st5'):
        assert result[key].shape == (_STATIONS,)


def test_coefficients_dividing_by_sin_are_nan_at_the_leading_edge():
    """ST2, ST4 and ST5 divide by sin(THNU), which vanishes at station N."""
    result = _run()
    for key in ('st2', 'st4', 'st5'):
        assert np.isnan(result[key][-1]), key
        assert np.all(np.isfinite(result[key][:-1])), key


def test_coefficients_not_dividing_by_sin_are_finite_throughout():
    result = _run()
    for key in ('st1', 'st3'):
        assert np.all(np.isfinite(result[key])), key


def test_symmetric_section_has_no_camber_coefficients():
    """ST4 and ST5 are built from the mean line, so they vanish without it."""
    result = _run(camber=0.0)
    assert np.allclose(result['st4'][:-1], 0.0, atol=1e-12)
    assert np.allclose(result['st5'][:-1], 0.0, atol=1e-12)


def test_thickness_coefficients_survive_a_symmetric_section():
    """ST1, ST2 and ST3 come from thickness and are unaffected by camber."""
    plain = _run(camber=0.0)
    cambered = _run(camber=0.04)
    for key in ('st1', 'st3'):
        assert np.allclose(plain[key], cambered[key])


# --------------------------------------------------------------------------
# Derived section parameters
# --------------------------------------------------------------------------

def test_symmetric_section_has_zero_derived_parameters():
    """The strongest correctness check available without the original."""
    result = _run(camber=0.0)
    assert result['alpha_zero_lift'] == pytest.approx(0.0, abs=1e-10)
    assert result['cli'] == pytest.approx(0.0, abs=1e-10)
    assert result['cm_c4'] == pytest.approx(0.0, abs=1e-10)
    assert result['ideal_alpha'] == pytest.approx(0.0, abs=1e-10)


def test_derived_parameters_scale_almost_linearly_with_camber():
    """Thin-airfoil theory is linear in the mean line, to within ASMINT.

    The analytical relationship is exactly linear, but the mean line is
    resampled onto the Weber stations by ASMINT, whose equal-axis-scaling
    step makes the spline shape depend on the data range. Doubling the
    camber therefore changes the fit slightly: the result is linear to about
    3e-7 relative rather than to machine precision. That nonlinearity is a
    property of the source's ASMINT, not of this routine.
    """
    single = _run(camber=0.02)
    double = _run(camber=0.04)
    for key in ('alpha_zero_lift', 'cli', 'cm_c4', 'ideal_alpha'):
        assert double[key] == pytest.approx(2.0 * single[key], rel=1e-5), key


def test_asmint_axis_scaling_is_the_source_of_that_nonlinearity():
    """Pin the cause, so the loose tolerance above is not read as slop."""
    from pydatcom.utils.interpolation import asmint
    x = np.linspace(0.0, 1.0, 61)
    position = 0.4

    def mean_line(amount):
        return np.where(x < position,
                        amount / position**2 * (2 * position * x - x**2),
                        amount / (1 - position)**2 *
                        ((1 - 2 * position) + 2 * position * x - x**2))

    stations = 0.5 * (np.cos(np.arange(1, 33) * np.pi / 32) + 1.0)
    single = np.asarray(asmint(x, mean_line(0.02), stations))
    double = np.asarray(asmint(x, mean_line(0.04), stations))
    departure = np.abs(double - 2.0 * single).max() / np.abs(double).max()
    assert 1e-9 < departure < 1e-5


def test_naca_2412_matches_published_section_data():
    """Independent check: published NACA 2412 has alpha0 near -2.0 degrees
    and a quarter-chord moment near -0.05."""
    result = _run(thickness=0.12, camber=0.02, position=0.4)
    assert -2.5 < result['alpha_zero_lift'] < -1.7
    assert -0.07 < result['cm_c4'] < -0.03
    assert 0.15 < result['cli'] < 0.40


def test_camber_makes_the_zero_lift_angle_negative():
    for camber in (0.01, 0.02, 0.04):
        assert _run(camber=camber)['alpha_zero_lift'] < 0.0


def test_ideal_lift_follows_from_the_two_angles():
    """CLI = 2*pi/RAD * (AI - ALO)."""
    from pydatcom.utils.constants import PI, RAD
    result = _run()
    assert result['cli'] == pytest.approx(
        2.0 * PI / RAD * (result['ideal_alpha'] - result['alpha_zero_lift']))


# --------------------------------------------------------------------------
# Branches and validation
# --------------------------------------------------------------------------

def test_supplied_cli_skips_the_derivation():
    """The source jumps past the AI/ALO/CLI block when CLI is already set."""
    result = _run(cli=0.30, alpha_zero_lift=-2.0)
    assert result['cli'] == pytest.approx(0.30)
    assert result['alpha_zero_lift'] == pytest.approx(-2.0)
    assert result['ideal_alpha'] is None
    # The moment is still computed, from the supplied zero-lift angle.
    assert np.isfinite(result['cm_c4'])


def test_supplied_cli_without_alpha_zero_is_rejected():
    with pytest.raises(ValueError, match="zero-lift angle"):
        _run(cli=0.30)


def test_supersonic_section_uses_the_other_resampler():
    """The source switches from the ASMINT spline to TBFUNX."""
    spline = _run(supersonic_section=False)
    linear = _run(supersonic_section=True)
    assert not np.allclose(spline['thickness'], linear['thickness'])
    assert np.all(np.isfinite(linear['thickness']))


def test_rejects_mismatched_section_arrays():
    x, t, c = _naca()
    with pytest.raises(ValueError):
        calculate_ideal(x, t[:10], c, 0.0158)
    with pytest.raises(ValueError):
        calculate_ideal([0.5], [0.01], [0.0], 0.0158)
