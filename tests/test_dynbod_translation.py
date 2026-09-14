"""
Regression tests for the DYNBOD translation.

Source of truth: datcom-legacy/datcom_2000/dynbod.f.

Both branches are closed-form in the source -- the hypersonic figures
7.2.1.1-9A and 7.2.1.2-12 are algebraic, not tabulated -- so these tests
check the expressions themselves, the branch structure, and the damping
signs that make the derivatives physically meaningful.
"""

import numpy as np
import pytest

from pydatcom.aerodynamics.dynbod import (
    calculate_dynbod_subsonic, calculate_dynbod_hypersonic,
)


def _body():
    x = np.array([0., 1., 2., 3., 4., 6., 8., 9., 9.5, 10.])
    r = np.array([0., .35, .6, .75, .8, .8, .8, .72, .62, .55])
    return x, np.pi * r**2


def _subsonic(**kwargs):
    x, s = _body()
    params = dict(cla_body=0.0274, cma_body=0.0109, body_length=10.0,
                  body_area=float(np.pi * 0.8**2),
                  base_area=float(np.pi * 0.55**2),
                  xcg=5.0, sref=float(np.pi * 0.8**2), cbar=10.0)
    params.update(kwargs)
    return calculate_dynbod_subsonic(x, s, **params)


def _hypersonic(**kwargs):
    params = dict(theta_nose=np.deg2rad(15.0),
                  theta_afterbody=np.deg2rad(5.0),
                  theta_flare=np.deg2rad(8.0),
                  d_nose=1.0, d1=1.6, d2=2.0,
                  nose_length=3.0, afterbody_length=8.0,
                  cna=[0.02, 0.01, 0.005], cma=[0.01, 0.004, 0.002],
                  xcg=5.0, sref=2.0, cbar=10.0)
    params.update(kwargs)
    return calculate_dynbod_hypersonic(**params)


# --------------------------------------------------------------------------
# Subsonic / transonic / supersonic branch
# --------------------------------------------------------------------------

def test_subsonic_matches_source_expressions():
    """CLQ, CMQ, CLAD and CMAD against an independent form of the source."""
    x, s = _body()
    result = _subsonic()

    volume = float(np.trapezoid(s, x)) if hasattr(np, 'trapezoid') else \
        float(np.trapz(s, x))
    centroid = (float(np.trapezoid(s * x, x)) if hasattr(np, 'trapezoid')
                else float(np.trapz(s * x, x))) / volume
    body_area = float(np.pi * 0.8**2)
    h = 5.0 / 10.0
    v = volume / (body_area * 10.0)
    arm = centroid / 10.0 - h
    denom = 1.0 - h - v
    area_ratio = float(np.pi * 0.55**2) / body_area

    assert result['volume'] == pytest.approx(volume, rel=1e-9)
    assert result['centroid'] == pytest.approx(centroid, rel=1e-9)
    assert result['clq'] == pytest.approx(
        2.0 * 0.0274 * (1.0 - h) * area_ratio * 1.0**2, rel=1e-9)
    assert result['clad'] == pytest.approx(
        2.0 * 0.0274 * v * area_ratio * 1.0, rel=1e-9)
    assert result['cmad'] == pytest.approx(
        2.0 * 0.0109 * (v * arm) / denom * area_ratio, rel=1e-9)


def test_acceleration_derivative_uses_first_power_of_length_ratio():
    """CLAD scales with LB/CBARR, the others with its square.

    The asymmetry is in the source and is dimensionally correct, so it must
    survive translation rather than being tidied away.
    """
    short = _subsonic(cbar=10.0)
    long = _subsonic(cbar=5.0)
    assert long['clq'] / short['clq'] == pytest.approx(4.0)
    assert long['cmq'] / short['cmq'] == pytest.approx(4.0)
    assert long['cmad'] / short['cmad'] == pytest.approx(4.0)
    assert long['clad'] / short['clad'] == pytest.approx(2.0)


def test_pitch_damping_is_negative():
    """CMQ and CMAD must damp for a statically stable body slope."""
    result = _subsonic(cma_body=0.0109)
    assert result['cmq'] < 0.0
    assert result['cmad'] < 0.0


def test_clq_falls_as_cg_moves_aft():
    """CLQ carries (1 - XCG/LB), so an aft CG reduces it."""
    forward = _subsonic(xcg=2.0)['clq']
    aft = _subsonic(xcg=8.0)['clq']
    assert forward > aft > 0.0


def test_subsonic_rejects_degenerate_configurations():
    x, s = _body()
    with pytest.raises(ValueError):
        calculate_dynbod_subsonic(x[:1], s[:1], 0.03, 0.01, 10.0, 2.0, 1.0,
                                  5.0, 2.0, 10.0)
    with pytest.raises(ValueError):
        _subsonic(sref=0.0)
    with pytest.raises(ValueError, match="denominator"):
        # Contrive 1 - XCG/LB - VB/(SB*LB) == 0.
        result = _subsonic()
        v = result['volume_ratio']
        _subsonic(xcg=10.0 * (1.0 - v))


# --------------------------------------------------------------------------
# Hypersonic branch: Figures 7.2.1.1-9A and 7.2.1.2-12
# --------------------------------------------------------------------------

def test_hypersonic_segment_sums():
    """CNQ and CMQ are the sums of their three segment contributions."""
    result = _hypersonic()
    assert result['cnq'] == pytest.approx(sum(result['cnq_segments']))
    assert result['cmq'] == pytest.approx(sum(result['cmq_segments']))


def test_nose_taper_is_fixed_at_zero():
    """The source pins LAMN to zero: the nose closes to a point."""
    assert _hypersonic()['taper_ratios'][0] == 0.0


def test_boattail_disables_the_flare_term():
    """A nonpositive flare angle sets LAMF to 1, zeroing that segment."""
    result = _hypersonic(theta_flare=-1.0)
    assert result['boat_tailed']
    assert result['taper_ratios'][2] == 1.0


def test_unit_taper_segment_contributes_nothing():
    """A segment of unit taper has no projected area change."""
    # d1 == d2 makes the flare taper unity.
    result = _hypersonic(d1=2.0, d2=2.0)
    assert result['taper_ratios'][2] == pytest.approx(1.0)
    # Its CNQ then reduces to the moment-arm term alone.
    arm = (5.0 - 8.0) / 10.0
    assert result['cnq_segments'][2] == pytest.approx(-2.0 * arm * 0.005)


def test_figure_7211_9a_expression():
    """One segment of Figure 7.2.1.1-9A against an independent form."""
    theta, taper, diameter = np.deg2rad(5.0), 0.625, 1.6
    scale = np.pi / (4.0 * 2.0 * 10.0)
    expected = (0.66667 / np.tan(theta) *
                (2.0 * (1.0 - taper**3) -
                 3.0 * taper * np.cos(theta)**2 * (1.0 - taper**2)) *
                scale * diameter**3)
    result = _hypersonic()
    arm = (5.0 - 3.0) / 10.0
    assert result['cnq_segments'][1] == pytest.approx(
        expected - 2.0 * arm * 0.01, rel=1e-9)


def test_figure_72112_12_expression():
    """One segment of Figure 7.2.1.2-12 against an independent form."""
    theta, taper, diameter = np.deg2rad(5.0), 0.625, 1.6
    scale = np.pi / (4.0 * 2.0 * 10.0) / 10.0
    a = 6.0 * taper**2 * (1.0 - taper**2)
    b = -8.0 * taper * (1.0 - taper**3)
    c = 3.0 * (1.0 - taper**4)
    cmqp = (-(a * np.cos(theta)**4 + b * np.cos(theta)**2 + c) /
            (6.0 * np.sin(theta)**2) - taper**4 / 2.0) * scale * diameter**4
    cnqp = (0.66667 / np.tan(theta) *
            (2.0 * (1.0 - taper**3) -
             3.0 * taper * np.cos(theta)**2 * (1.0 - taper**2)) *
            np.pi / (4.0 * 2.0 * 10.0) * diameter**3)
    arm = (5.0 - 3.0) / 10.0
    expected = (cmqp - 2.0 * arm * 0.004 + arm * cnqp -
                2.0 * arm**2 * 0.01)
    assert _hypersonic()['cmq_segments'][1] == pytest.approx(expected,
                                                             rel=1e-9)


def test_hypersonic_rejects_bad_input():
    with pytest.raises(ValueError):
        _hypersonic(sref=0.0)
    with pytest.raises(ValueError):
        _hypersonic(cna=[0.01, 0.01])
    with pytest.raises(ValueError, match="tan"):
        # A zero half angle with a non-unit taper is inconsistent.
        _hypersonic(theta_afterbody=0.0)
