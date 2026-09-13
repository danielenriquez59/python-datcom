"""Physics invariants for the high-priority DATCOM translation fixes."""

import numpy as np
import pytest

from pydatcom.aerodynamics.body_alone import calculate_body_alone_subsonic
from pydatcom.aerodynamics.drag import calculate_total_drag
from pydatcom.aerodynamics.hypersonic import calculate_hypersonic_coefficients
from pydatcom.aerodynamics.lift import (
    calculate_lift_curve_slope_compressible,
    calculate_wing_lift_subsonic,
)
from pydatcom.aerodynamics.moment import (
    calculate_axial_force_coefficient,
    calculate_normal_force_coefficient,
    calculate_wing_moment_coefficient,
    calculate_total_pitching_moment,
)
from pydatcom.aerodynamics.stability import calculate_static_stability_margin
from pydatcom.aerodynamics.supersonic import (
    calculate_supersonic_coefficients,
    calculate_supersonic_wave_drag,
)
from pydatcom.aerodynamics.subsonic import calculate_subsonic_coefficients
from pydatcom.aerodynamics.transonic import calculate_transonic_coefficients


def _wing_state(sref=100.0):
    return {
        'wing_aspect_ratio': 6.0,
        'wing_taper_ratio': 0.5,
        'wing_savsi': 0.0,
        'wing_tovc': 0.08,
        'wing_area': 50.0,
        'wing_mac': 4.0,
        'options_sref': sref,
        'options_cbarr': 4.0,
        'synths_xw': 10.0,
        'synths_xcg': 12.0,
        'body_swet_sref': 0.0,
        'drag_misc': 0.0,
    }


def _body_state():
    x = np.array([0.0, 2.0, 5.0, 8.0, 10.0])
    r = np.array([0.0, 1.0, 1.0, 1.0, 0.5])
    return {
        'body_nx': len(x),
        'body_x': x,
        'body_r': r,
        'body_s': np.pi * r**2,
        'body_p': 2.0 * np.pi * r,
        'options_sref': np.pi,
        'options_cbarr': 10.0,
        'options_rougfc': 1.0e-4,
        'synths_xcg': 5.0,
    }


def test_trsoni_finite_wing_slope_and_degree_units():
    ar, mach, sweep = 6.0, 0.6, 15.0
    beta2 = 1.0 - mach**2
    expected = 2.0 * np.pi * ar / (
        2.0 + np.sqrt(4.0 + ar**2 *
                      (beta2 + np.tan(np.deg2rad(sweep))**2))
    )
    assert calculate_lift_curve_slope_compressible(ar, 0.5, mach, sweep) == pytest.approx(expected)

    result = calculate_wing_lift_subsonic(_wing_state(50.0), 5.0, mach)
    assert result['cla_per_deg'] == pytest.approx(result['cla'] * np.pi / 180.0)
    assert result['cl'] == pytest.approx(result['cla_per_deg'] * 5.0)


def test_cli_alphai_and_wing_incidence_define_body_zero_lift_angle():
    state = _wing_state(50.0)
    state.update(wing_cli=0.12, wing_alphai=1.0, wing_clalpa=[0.1], synths_aliw=2.0)
    result = calculate_wing_lift_subsonic(state, -2.2, 0.6)
    assert result['section_alpha_zero'] == pytest.approx(-0.2)
    assert result['alpha_zero'] == pytest.approx(-2.2)
    assert result['cl'] == pytest.approx(0.0, abs=1e-14)


def test_reference_area_rescales_wing_lift_and_induced_drag_only():
    first = _wing_state(100.0)
    second = _wing_state(200.0)
    lift1 = calculate_wing_lift_subsonic(first, 5.0, 0.6)
    lift2 = calculate_wing_lift_subsonic(second, 5.0, 0.6)
    assert lift1['cl'] * 100.0 == pytest.approx(lift2['cl'] * 200.0)

    drag1 = calculate_total_drag(first, lift1['cl'], 0.6, 1.0e7)
    drag2 = calculate_total_drag(second, lift2['cl'], 0.6, 1.0e7)
    assert drag1['cd_induced'] * 100.0 == pytest.approx(drag2['cd_induced'] * 200.0)


def test_wind_body_force_rotation_is_invertible():
    cl, cd, alpha = 0.4, 0.07, 13.0
    cn = calculate_normal_force_coefficient(cl, cd, alpha)
    ca = calculate_axial_force_coefficient(cl, cd, alpha)
    angle = np.deg2rad(alpha)
    assert cn * np.cos(angle) - ca * np.sin(angle) == pytest.approx(cl)
    assert ca * np.cos(angle) + cn * np.sin(angle) == pytest.approx(cd)


def test_bodyrt_subsonic_is_symmetric_in_angle_of_attack():
    state = _body_state()
    positive = calculate_body_alone_subsonic(state, 10.0, 0.5, 1.0e7)
    negative = calculate_body_alone_subsonic(state, -10.0, 0.5, 1.0e7)
    assert negative['cn'] == pytest.approx(-positive['cn'])
    assert negative['cl'] == pytest.approx(-positive['cl'])
    assert negative['cm'] == pytest.approx(-positive['cm'])
    assert negative['cd'] == pytest.approx(positive['cd'])
    assert negative['ca'] == pytest.approx(positive['ca'])


def test_hypbod_is_symmetric_and_already_uses_sref():
    state = _body_state()
    positive = calculate_hypersonic_coefficients(state, 10.0, 6.0)
    negative = calculate_hypersonic_coefficients(state, -10.0, 6.0)
    assert positive['method'] == 'legacy_hypbod_nasa_tn_d176'
    assert negative['cl'] == pytest.approx(-positive['cl'])
    assert negative['cn'] == pytest.approx(-positive['cn'])
    assert negative['cm'] == pytest.approx(-positive['cm'])
    assert negative['cd'] == pytest.approx(positive['cd'])
    assert negative['ca'] == pytest.approx(positive['ca'])


def test_supersonic_moment_uses_dimensional_half_mac_location():
    state = _wing_state(50.0)
    result = calculate_supersonic_coefficients(state, 5.0, 2.0, 1.0e7)
    assert result['xac'] == pytest.approx(12.0)
    assert result['cm'] == pytest.approx(0.0, abs=1e-14)
    assert calculate_wing_moment_coefficient(0.4, 12.0, 12.0, 4.0, 4.0, -0.02) == pytest.approx(-0.02)


def test_subsonic_default_xac_uses_actual_mac_not_reference_chord():
    state = _wing_state(50.0)
    state.update(options_cbarr=10.0, synths_xcg=11.0)
    result = calculate_total_pitching_moment(state, 0.4, 5.0, 0.6)
    assert result['xac'] == pytest.approx(11.0)
    assert result['cm_wing'] == pytest.approx(0.0)


def test_supersonic_lift_drag_recovers_two_dimensional_limit():
    beta = np.sqrt(2.0**2 - 1.0)
    result = calculate_supersonic_wave_drag(2.0, 0.0, 1.0e12, 0.2)
    assert result['cd_wave_lift'] == pytest.approx(beta * 0.2**2 / 4.0, rel=1e-10)


def test_component_neutral_point_does_not_move_with_cg():
    state = _wing_state(100.0)
    state.update(
        wing_xac_abs=11.0,
        wing_cla=4.0,
        htail_area=20.0,
        htail_cla=3.0,
        synths_xh=24.0,
    )
    first = calculate_static_stability_margin(state, 0.6)
    state['synths_xcg'] = 14.0
    second = calculate_static_stability_margin(state, 0.6)
    assert first['xnp'] == pytest.approx(second['xnp'])
    assert first['stable']
    assert not second['stable']


def test_subsonic_to_transonic_boundary_is_continuous():
    state = _wing_state(50.0)
    below = calculate_subsonic_coefficients(state, 5.0, 0.9 - 1.0e-8, 1.0e7)
    boundary = calculate_transonic_coefficients(state, 5.0, 0.9, 1.0e7)
    for coefficient in ('cl', 'cd', 'cm'):
        assert below[coefficient] == pytest.approx(boundary[coefficient], rel=1e-6, abs=1e-9)
