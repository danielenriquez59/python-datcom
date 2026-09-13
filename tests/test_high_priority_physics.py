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
    calculate_cmalph_zero_lift_moment,
)
from pydatcom.aerodynamics.stability import calculate_static_stability_margin
from pydatcom.aerodynamics.supersonic import (
    calculate_supersonic_coefficients,
    calculate_supdrg_straight_wing,
    calculate_supdrg_skin_friction,
    calculate_supersonic_wave_drag,
)
from pydatcom.aerodynamics.subsonic import calculate_subsonic_coefficients
from pydatcom.aerodynamics.transonic import calculate_transonic_coefficients
from pydatcom.geometry.wing import calculate_straight_exposed_geometry


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


def test_cmalph_constant_section_zero_lift_moment():
    state = {
        'wing_area': 4.0,
        'wing_mac': 1.0,
        'wing_chrdr': 1.0,
        'wing_chrdtp': 1.0,
        'wing_sspn': 2.0,
        'wing_aspect_ratio': 4.0,
        'wing_savsi': 0.0,
        'wing_chstat': 0.25,
        'wing_cmo': -0.03,
        'wing_twista': 0.0,
        'options_sref': 4.0,
        'options_cbarr': 1.0,
    }
    # At M=0 CALM=1; unswept multiplier is AR/(/(AR+2)=2/3.
    assert calculate_cmalph_zero_lift_moment(state, 0.0) == pytest.approx(-0.02)
    # CMALPH clamps the Mach correction at its M=0.9 endpoint.
    assert calculate_cmalph_zero_lift_moment(state, 2.0) == pytest.approx(-0.02 * 1.445)


def test_wtgeom_straight_wing_uses_exposed_planform():
    state = {
        'wing_type': 1.0,
        'wing_chrdr': 4.0,
        'wing_chrdtp': 2.0,
        'wing_sspn': 5.0,
        'wing_sspne': 4.0,
        'wing_savsi': 0.0,
        'wing_chstat': 0.0,
        'synths_xw': 10.0,
    }
    result = calculate_straight_exposed_geometry(state)
    exposed_root = 3.6
    taper = 2.0 / exposed_root
    expected_mac = (2.0 * exposed_root * (1.0 + taper + taper**2) /
                    (3.0 * (1.0 + taper)))
    assert result['root_chord'] == pytest.approx(exposed_root)
    assert result['area'] == pytest.approx(4.0 * (exposed_root + 2.0))
    assert result['aspect_ratio'] == pytest.approx(4.0 * 4.0**2 / result['area'])
    assert result['mac'] == pytest.approx(expected_mac)
    assert result['exposed_root_x'] == pytest.approx(10.0)


def test_wtgeom_converts_sweep_reference_station():
    from pydatcom.geometry.wing import WingGeometry

    state = {
        'wing_type': 1.0,
        'wing_chrdr': 10.0,
        'wing_chrdtp': 5.0,
        'wing_sspn': 25.0,
        'wing_savsi': 30.0,
        'wing_chstat': 0.25,
    }
    wing = WingGeometry(state)
    expected_le = np.rad2deg(np.arctan(
        np.tan(np.deg2rad(30.0)) + 0.25 * (10.0 - 5.0) / 25.0))
    expected_half = np.rad2deg(np.arctan(
        np.tan(np.deg2rad(30.0)) + (0.5 - 0.25) * (5.0 - 10.0) / 25.0))
    assert wing.calculate_sweep_at_station(0.0) == pytest.approx(expected_le)
    assert wing.calculate_sweep_at_station(0.5) == pytest.approx(expected_half)


def test_supersonic_lift_drag_recovers_two_dimensional_limit():
    beta = np.sqrt(2.0**2 - 1.0)
    result = calculate_supersonic_wave_drag(2.0, 0.0, 1.0e12, 0.2)
    assert result['cd_wave_lift'] == pytest.approx(beta * 0.2**2 / 4.0, rel=1e-10)


def test_supdrg_straight_wing_figure_58_and_reference_basis():
    # beta*b/2/RLW = 0.4 lands exactly on the third source table entry.
    state = {
        'wing_type': 1.0,
        'wing_chrdr': 1.0,
        'wing_chrdtp': 1.0,
        'wing_sspn': 2.0,
        'wing_area': 4.0,
        'wing_aspect_ratio': 4.0,
        'wing_savsi': 0.0,
        'wing_chstat': 0.0,
        'wing_tceff': 0.05,
        'wing_leri': 0.0,
        'options_sref': 4.0,
    }
    result = calculate_supdrg_straight_wing(state, np.sqrt(1.04), 0.2)
    assert result['dragc'] == pytest.approx(0.593)
    assert result['p'] == pytest.approx(1.0)
    assert result['cd_wave_lift'] == pytest.approx(0.593 / (2.0 * np.pi) * 0.2**2)
    assert result['cd_wave_volume'] == pytest.approx(16.0 * 0.05**2 / (3.0 * 0.2))

    state['options_sref'] = 8.0
    scaled = calculate_supdrg_straight_wing(state, np.sqrt(1.04), 0.2)
    assert scaled['cd_wave_lift'] * 8.0 == pytest.approx(result['cd_wave_lift'] * 4.0)
    assert scaled['cd_wave_volume'] * 8.0 == pytest.approx(result['cd_wave_volume'] * 4.0)


def test_supdrg_straight_wing_skin_friction_and_roughness_cutoff():
    state = {
        'wing_type': 1.0,
        'wing_chrdr': 1.0,
        'wing_chrdtp': 1.0,
        'wing_sspn': 2.0,
        'options_sref': 8.0,
        'options_rougfc': 1.6e-4,
    }
    smooth = calculate_supdrg_skin_friction(state, 2.0, 1.0e7)
    assert smooth['cd_friction'] == pytest.approx(smooth['cf'])
    assert smooth['reynolds_used'] == pytest.approx(1.0e7)
    assert smooth['cf'] == pytest.approx(0.002304082970, rel=2e-9)

    state['options_rougfc'] = 1.0e-3
    rough = calculate_supdrg_skin_friction(state, 2.0, 1.0e9)
    expected_cutoff = (12.0 / 1.0e-3)**1.0482 * 10.0**1.98509
    assert rough['roughness_cutoff_reynolds'] == pytest.approx(expected_cutoff)
    assert rough['reynolds_used'] == pytest.approx(expected_cutoff)
    assert rough['cf'] == pytest.approx(0.00305908123105, rel=2e-9)

    capped = calculate_supdrg_skin_friction(state, 4.0, 1.0e9)
    at_three = calculate_supdrg_skin_friction(state, 3.0, 1.0e9)
    assert capped['mach_lookup'] == 3.0
    assert capped['roughness_cutoff_reynolds'] == pytest.approx(
        at_three['roughness_cutoff_reynolds'])


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
