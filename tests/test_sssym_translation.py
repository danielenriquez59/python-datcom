"""
Regression tests for the SSSYM translation.

Source of truth: datcom-legacy/datcom_2000/sssym.f.

SSSYM is closed-form with no figure tables, so these tests cover the three
moment contributions, the flap-type gate, and the physical signs a
trailing-edge flap must produce.
"""

import numpy as np
import pytest

from pydatcom.aerodynamics.sssym import calculate_sssym, TRAILING_EDGE


def _run(**kwargs):
    params = dict(deflections=[-10., 0., 10., 20.],
                  root_chord=6.0, sref=135.0,
                  span_inboard=3.0, span_outboard=9.0,
                  flap_chord=1.2, flap_area=14.0,
                  k3=0.9, taper_ratio_flap=0.5,
                  tan_hinge_line=np.tan(np.deg2rad(10.0)),
                  tan_le=np.tan(np.deg2rad(30.0)),
                  x_surface=10.0, chord_inboard=5.0, xcg=12.0,
                  bcmd1=-0.45, bcld1=0.55, bcld2=0.30)
    params.update(kwargs)
    return calculate_sssym(**params)


# --------------------------------------------------------------------------
# The flap-type gate
# --------------------------------------------------------------------------

def test_trailing_edge_flaps_are_computed():
    assert _run(flap_type=TRAILING_EDGE)['applicable']


@pytest.mark.parametrize("flap_type", [0.0, 2.0, 3.0])
def test_other_flap_types_produce_nothing(flap_type):
    """The source falls straight to its RETURN for anything but type 1."""
    result = _run(flap_type=flap_type)
    assert not result['applicable']
    assert 'delta_cm' not in result
    assert 'trailing edge' in result['reason']


# --------------------------------------------------------------------------
# The three moment contributions
# --------------------------------------------------------------------------

def test_total_moment_is_the_sum_of_its_three_parts():
    result = _run()
    assert result['cmd_total'] == pytest.approx(
        result['cmd1'] + result['cmd2'] + result['cmd3'])


def test_first_moment_term_carries_the_k1_taper_factor():
    """K1 = K3*(1 + TRTOFL + TRTOFL^2)."""
    result = _run(k3=0.9, taper_ratio_flap=0.5)
    assert result['k1'] == pytest.approx(0.9 * (1.0 + 0.5 + 0.25))
    assert result['cmd1'] == pytest.approx(
        result['k1'] / 3.0 * result['flap_span'] * 1.2 * -0.45 /
        (6.0 * 135.0))


def test_second_moment_term_vanishes_with_an_unswept_hinge_line():
    """K2 = K3*TANHL, so a zero hinge sweep removes that contribution."""
    result = _run(tan_hinge_line=0.0)
    assert result['k2'] == 0.0
    assert result['cmd2'] == 0.0
    assert result['cmd_total'] == pytest.approx(
        result['cmd1'] + result['cmd3'])


def test_third_moment_term_vanishes_when_the_flap_sits_at_the_cg():
    """CMD3 carries the flap's arm ahead of the moment reference."""
    base = _run()
    at_cg = _run(xcg=base['moment_arm'] + 12.0)
    assert at_cg['moment_arm'] == pytest.approx(0.0)
    assert at_cg['cmd3'] == pytest.approx(0.0)


def test_moment_arm_matches_the_source_expression():
    """DELXF = XW + ALOCI*TANLE + CI - CFI - XCG."""
    result = _run()
    assert result['moment_arm'] == pytest.approx(
        10.0 + 3.0 * np.tan(np.deg2rad(30.0)) + 5.0 - 1.2 - 12.0)


def test_flap_span_is_twice_the_station_difference():
    """BEF = 2*(ALOCO - ALOCI), both panels."""
    assert _run()['flap_span'] == pytest.approx(2.0 * (9.0 - 3.0))


# --------------------------------------------------------------------------
# Increments and physical signs
# --------------------------------------------------------------------------

def test_increments_are_linear_in_deflection():
    result = _run()
    for key in ('delta_cm', 'delta_cl'):
        steps = np.diff(result[key])
        assert np.allclose(steps, steps[0])


def test_zero_deflection_gives_zero_increment():
    result = _run(deflections=[0.0])
    assert result['delta_cm'][0] == 0.0
    assert result['delta_cl'][0] == 0.0


def test_trailing_edge_down_raises_lift_and_lowers_the_nose():
    """A positive deflection must add lift and a nose-down moment."""
    result = _run(deflections=[10.0])
    assert result['cld'] > 0.0
    assert result['delta_cl'][0] > 0.0
    assert result['cmd_total'] < 0.0
    assert result['delta_cm'][0] < 0.0


def test_lift_effectiveness_scales_with_flap_area():
    """CLD = K3*BCLD1*SF/SR."""
    small = _run(flap_area=7.0)['cld']
    large = _run(flap_area=14.0)['cld']
    assert large == pytest.approx(2.0 * small)


def test_a_larger_flap_span_increases_the_moment_terms():
    narrow = _run(span_outboard=6.0)
    wide = _run(span_outboard=9.0)
    assert abs(wide['cmd1']) > abs(narrow['cmd1'])


# --------------------------------------------------------------------------
# Validation
# --------------------------------------------------------------------------

def test_rejects_nonpositive_references():
    with pytest.raises(ValueError):
        _run(root_chord=0.0)
    with pytest.raises(ValueError):
        _run(sref=0.0)
