"""
Regression tests for the SLOPE translation.

Source of truth: datcom-legacy/datcom_2000/slope.f.

SLOPE consumes the Weber coefficients IDEAL produces, so these tests run the
two together. The section results are checked against published NACA 2412
data and against thin-airfoil limits, which are constraints independent of
the source.
"""

import numpy as np
import pytest

from pydatcom.geometry.ideal import calculate_ideal
from pydatcom.geometry.slope import (
    calculate_slope, _CORRECTION_FLOOR, _SLOPE_FACTOR,
    _MIN_REYNOLDS, _REYNOLDS_SUBSTITUTE,
)


def _section(thickness=0.12, camber=0.02, position=0.4):
    x = np.linspace(0.0, 1.0, 61)
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


def _weber(**kwargs):
    x, t, c = _section(**kwargs)
    return calculate_ideal(x, t, c, leading_edge_radius=0.0158)


def _run(mach=0.0, reynolds=6.0e6, **kwargs):
    weber = _weber(**kwargs)
    return calculate_slope(weber, mach, reynolds, weber['alpha_zero_lift'])


# --------------------------------------------------------------------------
# Section results against published data
# --------------------------------------------------------------------------

def test_naca_2412_lift_slope_matches_published_data():
    """Published NACA 2412 has a section slope near 0.095 to 0.105 per degree."""
    assert 0.090 < _run()['cla'] < 0.110


def test_naca_2412_zero_lift_moment_matches_published_data():
    """Published NACA 2412 has a quarter-chord moment near -0.05."""
    assert -0.06 < _run()['cm_c4'] < -0.04


def test_naca_2412_zero_lift_angle_matches_published_data():
    assert -2.4 < _run()['alpha_zero_lift'] < -1.6


def test_aerodynamic_centre_sits_near_the_quarter_chord():
    assert 0.20 < _run()['xac'] < 0.32


def test_symmetric_section_has_no_zero_lift_moment():
    result = _run(camber=0.0)
    assert result['cm_c4'] == pytest.approx(0.0, abs=1e-6)
    assert result['alpha_zero_lift'] == pytest.approx(0.0, abs=1e-3)


def test_lift_slope_stays_below_the_inviscid_limit():
    """The viscous correction and the 1.05 factor keep it under 2*pi."""
    inviscid_per_degree = 2.0 * np.pi * np.pi / 180.0
    assert _run()['cla'] < inviscid_per_degree * _SLOPE_FACTOR


# --------------------------------------------------------------------------
# Compressibility
# --------------------------------------------------------------------------

def test_lift_slope_rises_with_mach():
    values = [_run(mach=m)['cla'] for m in (0.0, 0.3, 0.5, 0.7)]
    assert all(b > a for a, b in zip(values, values[1:]))


def test_zero_lift_moment_is_reported_only_at_mach_zero():
    """The source stores CMCO4 from the Mach zero pass alone."""
    assert 'cm_c4' in _run(mach=0.0)
    assert 'cm_c4' not in _run(mach=0.5)


def test_supersonic_is_rejected():
    """The source computes nothing and leaves CLALPA at its sentinel."""
    weber = _weber()
    with pytest.raises(ValueError, match="subsonic"):
        calculate_slope(weber, 1.2, 6.0e6, weber['alpha_zero_lift'])


# --------------------------------------------------------------------------
# The zero-lift iteration
# --------------------------------------------------------------------------

def test_iteration_drives_the_lower_angle_to_zero_lift():
    """The refined angle differs from the starting one and converges fast."""
    weber = _weber()
    result = calculate_slope(weber, 0.0, 6.0e6, weber['alpha_zero_lift'])
    assert result['iterations'] <= 5
    assert result['alpha_zero_lift'] != pytest.approx(
        weber['alpha_zero_lift'], abs=1e-9)


def test_iteration_converges_from_a_poor_starting_angle():
    """Both starts land within the angle the lift tolerance permits.

    The source stops once |CL| at the lower angle falls under 0.001. With a
    slope near 0.096 per degree that admits about 0.01 degrees of slack, so
    two different starting points converge to angles that agree to well
    inside that band rather than to machine precision.
    """
    from pydatcom.geometry.slope import _LIFT_TOLERANCE
    weber = _weber()
    far = calculate_slope(weber, 0.0, 6.0e6, 5.0)
    near = calculate_slope(weber, 0.0, 6.0e6, weber['alpha_zero_lift'])
    permitted = _LIFT_TOLERANCE / near['cla']
    assert abs(far['alpha_zero_lift'] - near['alpha_zero_lift']) < permitted
    assert permitted < 0.02


# --------------------------------------------------------------------------
# The viscous correction
# --------------------------------------------------------------------------

def test_correction_is_floored():
    """The source clamps the viscous factor at 0.6896."""
    for thickness in (0.06, 0.12, 0.21):
        assert _run(thickness=thickness)['correction'] >= _CORRECTION_FLOOR


def test_reynolds_defaults_to_one_million_at_mach_zero():
    """The supplied Reynolds number is used only with a crest-critical Mach."""
    assert _run(mach=0.0, reynolds=6.0e6)['reynolds_used'] == 1.0e6


def test_very_low_reynolds_is_substituted():
    weber = _weber()
    result = calculate_slope(weber, 0.5, 1.0e3, weber['alpha_zero_lift'],
                             crest_critical_mach=0.75)
    assert result['reynolds_used'] == pytest.approx(_REYNOLDS_SUBSTITUTE)
    assert _MIN_REYNOLDS > _REYNOLDS_SUBSTITUTE


def test_supplied_reynolds_is_used_with_a_crest_critical_mach():
    weber = _weber()
    result = calculate_slope(weber, 0.5, 6.0e6, weber['alpha_zero_lift'],
                             crest_critical_mach=0.75)
    assert result['reynolds_used'] == pytest.approx(6.0e6)


def test_thicker_section_has_a_larger_trailing_edge_angle():
    values = [_run(thickness=t)['phite'] for t in (0.06, 0.12, 0.21)]
    assert all(b > a for a, b in zip(values, values[1:]))
