"""
Regression tests for the VRTCDO / VFCDO translation.

Source of truth: datcom-legacy/datcom_2000/vrtcdo.f, vfcdo.f.

The two source routines are identical apart from the COMMON offsets that
select which surface they read, so one translation and one test file cover
both the vertical tail and the ventral fin.
"""

import numpy as np
import pytest

from pydatcom.aerodynamics.vertical_panel import calculate_vertical_panel_cdo

_SWEEP = np.deg2rad(35.0)


def _panel(**kwargs):
    params = dict(mach=2.0, reynolds_per_length=2.0e6, sref=135.0,
                  mac_inboard=4.0, area_inboard=20.0,
                  root_chord=5.0, tip_chord=3.0, semispan=5.0,
                  tan_le=np.tan(_SWEEP), cos_le=np.cos(_SWEEP),
                  thickness_ratio=0.05, leading_edge_radius=0.01)
    params.update(kwargs)
    return calculate_vertical_panel_cdo(**params)


def _cranked(**kwargs):
    outboard = np.deg2rad(45.0)
    params = dict(straight=False, mac_outboard=3.0, area_outboard=8.0,
                  break_chord=4.0, span_outboard=2.0,
                  tan_le_outboard=np.tan(outboard),
                  cos_le_outboard=np.cos(outboard),
                  tan_te_outboard=np.tan(np.deg2rad(10.0)),
                  leading_edge_radius_outboard=0.012)
    params.update(kwargs)
    return _panel(**params)


# --------------------------------------------------------------------------
# Skin friction
# --------------------------------------------------------------------------

def test_straight_friction_uses_two_faces_of_the_panel_area():
    """CDF = CF * SRSTAR/SR * 2, the factor of two being both faces."""
    result = _panel()
    assert result['cd_friction'] == pytest.approx(
        result['cf_inboard'] * 20.0 / 135.0 * 2.0, rel=1e-12)
    assert result['cf_outboard'] is None


def test_cranked_friction_area_weights_the_two_panels():
    """CDF = (CFI*SISTAR + CFO*SOSTAR)/SR * 2."""
    result = _cranked()
    assert result['cd_friction'] == pytest.approx(
        (result['cf_inboard'] * 20.0 + result['cf_outboard'] * 8.0) /
        135.0 * 2.0, rel=1e-12)


def test_friction_falls_with_mach():
    values = [_panel(mach=m)['cd_friction'] for m in (1.2, 2.0, 3.0)]
    assert all(b < a for a, b in zip(values, values[1:]))


def test_roughness_cutoff_limits_the_reynolds_number():
    smooth = _panel(roughness=0.4e-4)
    rough = _panel(roughness=6.4e-4)
    assert (rough['roughness_cutoff_reynolds'] <
            smooth['roughness_cutoff_reynolds'])
    assert rough['cf_inboard'] >= smooth['cf_inboard']


def test_zero_roughness_skips_the_cutoff():
    """The source's RUFF=0 path takes the raw Reynolds number."""
    result = _panel(roughness=0.0)
    assert result['roughness_cutoff_reynolds'] is None
    assert result['reynolds_used'] == pytest.approx(4.0 * 2.0e6)


def test_mach_lookup_is_capped_at_three():
    """Above Mach 3 the friction lookup reuses the Mach 3 values."""
    assert _panel(mach=3.0)['cf_inboard'] == pytest.approx(
        _panel(mach=4.0)['cf_inboard'])


# --------------------------------------------------------------------------
# Wave drag
# --------------------------------------------------------------------------

def test_sharp_leading_edge_uses_ksharp():
    """With KSHARP supplied there is no leading-edge bluntness term."""
    result = _panel(ksharp=1.2)
    assert result['leading_edge'] == 'sharp'
    assert result['cd_leading_edge'] == 0.0
    beta = np.sqrt(3.0)
    area = (5.0 + 3.0) / 2.0 * 5.0
    assert result['cd_wave'] == pytest.approx(
        1.2 * 0.05**2 * area / 135.0 / beta, rel=1e-12)


def test_round_leading_edge_adds_a_bluntness_term():
    result = _panel()
    assert result['leading_edge'] == 'round'
    assert result['cd_leading_edge'] > 0.0
    assert result['cd_wave'] > result['cd_leading_edge']


def test_sharp_leading_edge_has_less_wave_drag_than_round():
    assert _panel(ksharp=1.2)['cd_wave'] < _panel()['cd_wave']


def test_sonic_leading_edge_switches_the_denominator():
    """beta/tan(sweep) >= 1 divides by beta, otherwise by tan(sweep)."""
    # tan(35 deg) = 0.7002, so the switch is at M = sqrt(1 + 0.7002^2).
    switch = np.sqrt(1.0 + np.tan(_SWEEP)**2)
    assert not _panel(mach=switch - 0.02)['sonic_leading_edge']
    assert _panel(mach=switch + 0.02)['sonic_leading_edge']


def test_wave_drag_falls_with_mach_beyond_the_switch():
    values = [_panel(mach=m)['cd_wave'] for m in (1.5, 2.0, 2.5, 3.0)]
    assert all(b < a for a, b in zip(values, values[1:]))


def test_cranked_wave_drag_rebuilds_the_outboard_trapezoid():
    """S = (CB + SPANIN*(tanLEo - tanTEo) + CT) * SPANS / 2."""
    result = _cranked()
    span_inboard = 5.0 - 2.0
    root_equivalent = 4.0 + span_inboard * (np.tan(np.deg2rad(45.0)) -
                                            np.tan(np.deg2rad(10.0)))
    assert result['wave_drag_area'] == pytest.approx(
        (root_equivalent + 3.0) * 5.0 * 0.5, rel=1e-12)


def test_straight_wave_drag_uses_the_panel_trapezoid():
    assert _panel()['wave_drag_area'] == pytest.approx(
        (5.0 + 3.0) / 2.0 * 5.0)


# --------------------------------------------------------------------------
# Validation
# --------------------------------------------------------------------------

def test_subsonic_is_rejected():
    with pytest.raises(ValueError, match="supersonic"):
        _panel(mach=0.8)


def test_unswept_panel_is_rejected():
    """The source floors tan(sweep) at UNUSED rather than allowing zero."""
    with pytest.raises(ValueError, match="tan"):
        _panel(tan_le=0.0)


def test_cranked_without_outboard_geometry_is_rejected():
    with pytest.raises(ValueError, match="outboard"):
        _panel(straight=False)


def test_bad_references_are_rejected():
    with pytest.raises(ValueError):
        _panel(sref=0.0)
    with pytest.raises(ValueError):
        _panel(roughness=-1.0)
