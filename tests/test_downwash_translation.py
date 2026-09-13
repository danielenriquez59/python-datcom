"""
Regression tests for the DWASH / INFTGM / CLWBT / CDWBT translations.

Source of truth: datcom-legacy/datcom_2000/dwash.f, inftgm.f, clwbt.f,
cdwbt.f.  The figure tests reproduce the source DATA statements exactly; the
gradient tests check the translated result against an independently written
form of the published DATCOM Section 4.4.1 equation.
"""

import numpy as np
import pytest

from pydatcom.aerodynamics.downwash import (
    fig4417_68a, fig4417_68b, calculate_downwash_geometry,
    calculate_downwash_gradient_441, calculate_downwash,
)
from pydatcom.aerodynamics.wing_body_tail import (
    calculate_clwbt, calculate_cdwbt, calculate_tail_load,
)
from pydatcom.aerodynamics.moment import calculate_total_pitching_moment
from pydatcom.geometry.wing import calculate_straight_exposed_geometry


def _aircraft_state():
    """A straight-tapered wing-tail aircraft with complete geometry."""
    return {
        'wing_type': 1.0, 'wing_chrdr': 6.0, 'wing_chrdtp': 3.0,
        'wing_sspn': 15.0, 'wing_sspne': 13.5,
        'wing_savsi': 0.0, 'wing_chstat': 0.25,
        'htail_type': 1.0, 'htail_chrdr': 3.0, 'htail_chrdtp': 1.5,
        'htail_sspn': 6.0, 'htail_sspne': 5.4,
        'htail_savsi': 0.0, 'htail_chstat': 0.25,
        'synths_xw': 10.0, 'synths_zw': 0.0, 'synths_aliw': 0.0,
        'synths_xh': 35.0, 'synths_zh': 2.0, 'synths_alih': 0.0,
        'synths_xcg': 12.0,
        'options_sref': 135.0, 'options_cbarr': 4.6667,
        'flight_mach': 0.3,
    }


# --------------------------------------------------------------------------
# Figure 4.4.1-68A: dwash.f DATA X4157A / Y4157A
# --------------------------------------------------------------------------

_FIG68A_SOURCE = [
    (0., 1.5), (19., 1.5), (20., 1.5), (25., 1.57), (30., 1.68),
    (35., 1.81), (40., 2.0), (45., 2.25), (50., 2.75), (54., 3.0),
]


@pytest.mark.parametrize("sweep,expected", _FIG68A_SOURCE)
def test_fig4417_68a_source_coordinates(sweep, expected):
    """Every source table coordinate is reproduced exactly."""
    assert fig4417_68a(sweep) == pytest.approx(expected, abs=1e-12)


def test_fig4417_68a_is_monotonic_in_sweep():
    """Separation onset requires a larger DELTAY as sweep increases."""
    sweeps = np.linspace(0.0, 54.0, 100)
    values = [fig4417_68a(s) for s in sweeps]
    assert np.all(np.diff(values) >= -1e-12)


# --------------------------------------------------------------------------
# Figure 4.4.1-68B: dwash.f DATA X157B1 / X157B2 / Y4157B
# --------------------------------------------------------------------------

_FIG68B_AOBV = [0.0, 0.2, 0.6, 1.0]
_FIG68B_BHBV = [.2, .3, .4, .5, .6, .7, .8, .833, 1.0, 1.1, 1.2, 1.3,
                1.4, 1.5]
# Y4157B(14,4), column by column exactly as the DATA statement lists it.
_FIG68B_COLUMNS = [
    [1., 1.03, 1.068, 1.11, 1.176, 1.26, 1.36, 1.34, 1.20, 1.11, 1.0,
     .88, .74, .6],
    [.96, .96, .965, .976, .98, 1.016, 1.036, 1.04, .946, .88, .8,
     .69, .58, .46],
    [.74, .73, .72, .71, .69, .67, .644, .638, .582, .56, .52, .48,
     .44, .4],
    [.5, .5, .5, .49, .476, .46, .448, .44, .4, .372, .352, .32, .3, .27],
]


def test_fig4417_68b_reproduces_every_source_entry():
    """All 56 original table entries are returned at their own grid point."""
    for col, aobv in enumerate(_FIG68B_AOBV):
        for row, bhbv in enumerate(_FIG68B_BHBV):
            assert fig4417_68b(aobv, bhbv) == pytest.approx(
                _FIG68B_COLUMNS[col][row], abs=1e-12), (
                    f"mismatch at a/(bv/2)={aobv}, bh/bv={bhbv}")


def test_fig4417_68b_decreases_with_tail_height():
    """Raising the tail out of the wake reduces the average downwash."""
    for bhbv in (0.4, 0.8, 1.2):
        values = [fig4417_68b(a, bhbv) for a in _FIG68B_AOBV]
        assert np.all(np.diff(values) < 0.0)


def test_fig4417_68b_is_clamped_at_zero():
    """DWASH clamps DEBODE at zero; extrapolation must not go negative."""
    assert fig4417_68b(1.0, 3.0) >= 0.0


# --------------------------------------------------------------------------
# INFTGM downwash synthesizing dimensions
# --------------------------------------------------------------------------

def test_inftgm_geometry_matches_hand_calculation():
    """A(193), A(194), A(24) and A(12) for zero incidence and no hinge axis."""
    state = _aircraft_state()
    geometry = calculate_downwash_geometry(state)
    tail = calculate_straight_exposed_geometry(state, component='htail')

    # A(193) = XH - XW - CHRDR (both incidences zero)
    assert geometry['a193'] == pytest.approx(35.0 - 10.0 - 6.0)
    # DXBH = XBRSTH + DXSTAR: the exposed MAC quarter chord aft of the
    # exposed root leading edge, carried inboard to the centerline.  SAVSI=0
    # with CHSTAT=0.25 is an unswept quarter chord, so the leading edge is
    # still swept by tan_le = 0.25*(CHRDR-CHRDTP)/SSPN.
    buried = 6.0 - 5.4
    assert tail['tan_le'] == pytest.approx(0.25 * (3.0 - 1.5) / 6.0)
    assert geometry['dxbh'] == pytest.approx(
        tail['mac'] / 4.0 + tail['mac_span_location'] * tail['tan_le'] +
        buried * tail['tan_le'])
    assert geometry['a194'] == pytest.approx(geometry['a193'] +
                                             geometry['dxbh'])
    # With zero incidence the arm is A(194) and the height is ZH - ZW.
    assert geometry['tail_arm'] == pytest.approx(geometry['a194'])
    assert geometry['tail_height'] == pytest.approx(2.0)
    assert geometry['tail_angle'] == pytest.approx(
        np.arctan2(2.0, geometry['a194']))


def test_inftgm_hinge_axis_rotation():
    """HINAX rotates the tail reference point before the arm is formed."""
    state = _aircraft_state()
    base = calculate_downwash_geometry(state)
    state['synths_alih'] = 5.0
    state['synths_hinax'] = 34.0
    rotated = calculate_downwash_geometry(state)
    alih = np.deg2rad(5.0)
    assert rotated['zh_rotated'] == pytest.approx(
        2.0 + np.sin(alih) * (34.0 - 35.0))
    assert rotated['xh_rotated'] == pytest.approx(
        34.0 * (1.0 - np.cos(alih)) + 35.0 * np.cos(alih))
    assert rotated['tail_arm'] != base['tail_arm']


def test_inftgm_wing_incidence_changes_tail_height():
    """Wing incidence tilts the chord plane the tail height is measured from."""
    state = _aircraft_state()
    level = calculate_downwash_geometry(state)
    state['synths_aliw'] = 3.0
    tilted = calculate_downwash_geometry(state)
    assert tilted['tail_height'] != pytest.approx(level['tail_height'])


# --------------------------------------------------------------------------
# DATCOM Section 4.4.1 downwash gradient
# --------------------------------------------------------------------------

def test_gradient_matches_published_section_441_equation():
    """Independent form of de/da = 4.44*(KA*KL*KH*sqrt(cos Lc/4))**1.19."""
    state = _aircraft_state()
    result = calculate_downwash_gradient_441(state)
    wing = calculate_straight_exposed_geometry(state, component='wing')

    ar = wing['aspect_ratio']
    taper = wing['taper_ratio']
    sspne = wing['semispan']
    ka = 1.0 / ar - 1.0 / (1.0 + ar**1.7)
    kl = (10.0 - 3.0 * taper) / 7.0
    # Source XLH: reference separation less each MAC quarter-chord offset.
    xlh = (35.0 - wing_tail_offset(state, 'htail')) - (
        10.0 - wing_tail_offset(state, 'wing'))
    height = calculate_downwash_geometry(state)['tail_height']
    kh = (1.0 - abs(0.5 * height / sspne)) / (xlh / sspne)**(1.0 / 3.0)
    cos_c4 = 1.0 / np.sqrt(1.0 + wing['tan_c4']**2)
    expected = 4.44 * (ka * kl * kh * np.sqrt(cos_c4))**1.19

    assert result['k_a'] == pytest.approx(ka, rel=1e-12)
    assert result['k_lambda'] == pytest.approx(kl, rel=1e-12)
    assert result['k_h'] == pytest.approx(kh, rel=1e-12)
    assert result['deda'] == pytest.approx(expected, rel=1e-12)


def wing_tail_offset(state, component):
    """Theoretical MAC quarter-chord offset used by the source XLH."""
    return calculate_straight_exposed_geometry(
        state, component=component)['mac_c4_theoretical']


def test_gradient_is_physically_reasonable():
    """A conventional tail-aft layout has de/da between 0.2 and 0.6."""
    result = calculate_downwash_gradient_441(_aircraft_state())
    assert 0.2 < result['deda'] < 0.6


def test_gradient_falls_with_aspect_ratio():
    """Higher aspect ratio produces weaker downwash at the tail."""
    gradients = []
    for sspn, sspne in ((9.0, 8.1), (15.0, 13.5), (21.0, 18.9)):
        state = _aircraft_state()
        state['wing_sspn'] = sspn
        state['wing_sspne'] = sspne
        gradients.append(calculate_downwash_gradient_441(state)['deda'])
    assert gradients[0] > gradients[1] > gradients[2]


def test_gradient_falls_as_tail_is_raised():
    """Moving the tail away from the wake reduces the gradient."""
    low = _aircraft_state()
    low['synths_zh'] = 0.0
    high = _aircraft_state()
    high['synths_zh'] = 6.0
    assert (calculate_downwash_gradient_441(high)['deda'] <
            calculate_downwash_gradient_441(low)['deda'])


def test_gradient_falls_with_longer_tail_arm():
    """A longer arm puts the tail in weaker induced flow."""
    near = _aircraft_state()
    near['synths_xh'] = 28.0
    far = _aircraft_state()
    far['synths_xh'] = 45.0
    assert (calculate_downwash_gradient_441(far)['deda'] <
            calculate_downwash_gradient_441(near)['deda'])


def test_gradient_rejects_nonpositive_tail_arm():
    """A tail ahead of the wing is not a valid Section 4.4.1 configuration."""
    state = _aircraft_state()
    state['synths_xh'] = 2.0
    with pytest.raises(ValueError, match="positive tail arm"):
        calculate_downwash_gradient_441(state)


def test_downwash_angle_is_linear_in_alpha():
    """The 4.4.1 gradient is constant, so eps is linear through alpha_0L."""
    state = _aircraft_state()
    deda = calculate_downwash_gradient_441(state)['deda']
    for alpha in (-6.0, 0.0, 4.0, 10.0):
        result = calculate_downwash(state, alpha)
        assert result['eps_deg'] == pytest.approx(deda * alpha)
    assert calculate_downwash(state, 0.0)['eps_deg'] == pytest.approx(0.0)


def test_downwash_angle_is_less_than_alpha():
    """Downwash cannot exceed the angle of attack that produced it."""
    state = _aircraft_state()
    for alpha in (2.0, 5.0, 10.0):
        assert 0.0 < calculate_downwash(state, alpha)['eps_deg'] < alpha


# --------------------------------------------------------------------------
# CLWBT / CDWBT
# --------------------------------------------------------------------------

def test_clwbt_matches_source_expression():
    """CLBWH = ((CLH-CLI)*(KHB+KBH) + CLI*(KKHB+KKBH))*QOQI + CLBW."""
    khb, kbh, kkhb, kkbh = 1.15, 0.22, 1.05, 0.18
    clah, alih, qoqi = 0.062, 2.0, 0.92
    clh, clbw = 0.31, 0.48
    result = calculate_clwbt(clbw, clh, clah, alih, 1.7, qoqi,
                             khb, kbh, kkhb, kkbh)
    cli = clah * alih
    expected = ((clh - cli) * (khb + kbh) + cli * (kkhb + kkbh)) * qoqi + clbw
    assert result['cl_total'] == pytest.approx(expected, rel=1e-12)


def test_clwbt_slope_matches_source_expression():
    """CLABWH = CLABW + (KHB+KBH)*CLAH*(1-DEODA)*QOQI."""
    result = calculate_clwbt(0.4, 0.3, 0.062, 0.0, 1.5, 0.9,
                             khb=1.1, kbh=0.2,
                             cla_wing_body=0.081, deda=0.35)
    expected = 0.081 + (1.1 + 0.2) * 0.062 * (1.0 - 0.35) * 0.9
    assert result['cla_total'] == pytest.approx(expected, rel=1e-12)


def test_clwbt_no_carryover_reduces_to_isolated_tail():
    """KHB=1, KBH=0 leaves the tail load scaled only by the q ratio."""
    result = calculate_clwbt(0.5, 0.3, 0.06, 0.0, 1.0, 1.0)
    assert result['cl_total'] == pytest.approx(0.8)


def test_cdwbt_matches_source_expression():
    """WBTCD = CDOV + CDWB + QOQI*(CDH*cos(EPS) + CLH*sin(EPS))."""
    eps = 2.4
    result = calculate_cdwbt(0.031, 0.0042, 0.26, eps, 0.93,
                             cdo_vertical=0.0019)
    rad = np.deg2rad(eps)
    expected = 0.0019 + 0.031 + 0.93 * (0.0042 * np.cos(rad) +
                                        0.26 * np.sin(rad))
    assert result['cd_total'] == pytest.approx(expected, rel=1e-12)


def test_cdwbt_tail_lift_adds_drag_through_downwash():
    """A lifting tail in downwash contributes a positive drag increment."""
    lifting = calculate_cdwbt(0.03, 0.004, 0.25, 3.0, 1.0)
    unloaded = calculate_cdwbt(0.03, 0.004, 0.0, 3.0, 1.0)
    assert lifting['cd_total'] > unloaded['cd_total']


def test_cdwbt_zero_downwash_is_plain_summation():
    """With EPS=0 the tail contributes only its own drag."""
    result = calculate_cdwbt(0.03, 0.004, 0.25, 0.0, 1.0)
    assert result['cd_total'] == pytest.approx(0.034)


def test_cdwbt_rejects_negative_dynamic_pressure_ratio():
    with pytest.raises(ValueError):
        calculate_cdwbt(0.03, 0.004, 0.25, 1.0, -0.1)


# --------------------------------------------------------------------------
# Tail load and the aircraft moment
# --------------------------------------------------------------------------

def test_tail_lift_slope_is_per_degree():
    """The tail slope must be per degree, matching the degree angles used.

    calculate_lift_curve_slope_compressible returns a per-radian slope.
    Using it directly against a degree angle inflates the tail load by
    180/pi, which is how a 3965% static margin reached an example run.
    """
    state = _aircraft_state()
    downwash = calculate_downwash(state, 10.0)
    load = calculate_tail_load(state, 10.0, downwash)
    # A subsonic tail of AR ~5 has a slope near 2*pi*AR/(2+AR) per radian,
    # i.e. well under 0.12 per degree.
    assert 0.02 < load['cla_tail'] < 0.12
    assert 0.0 < load['cl_tail'] < 0.5


def test_tail_moment_magnitude_is_physical():
    """Cm_tail must stay in a plausible range for a conventional layout."""
    state = _aircraft_state()
    for alpha in (-4.0, 0.0, 5.0, 10.0):
        result = calculate_total_pitching_moment(state, 0.5, alpha, 0.3)
        assert abs(result['cm_tail']) < 1.5, (
            f"Cm_tail={result['cm_tail']:.4f} at alpha={alpha} is outside "
            "any physical range for this configuration")


def test_static_margin_is_physical():
    """dCm/dCL must give a static margin of a few tens of percent MAC."""
    state = _aircraft_state()
    from pydatcom.aerodynamics.lift import calculate_wing_lift_subsonic
    alphas = np.array([-4.0, -2.0, 0.0, 2.0, 4.0, 6.0])
    cls, cms = [], []
    for alpha in alphas:
        cl = calculate_wing_lift_subsonic(state, alpha, 0.3)['cl']
        cls.append(cl)
        cms.append(calculate_total_pitching_moment(
            state, cl, alpha, 0.3)['cm_total'])
    dcm_dcl = np.polyfit(cls, cms, 1)[0]
    assert -2.0 < dcm_dcl < 0.0, (
        f"dCm/dCL={dcm_dcl:.4f} implies a {-dcm_dcl*100:.0f}% static margin")


def test_tail_load_uses_local_flow_angle():
    """ALPT = alpha - eps, and the load is normalized to SREF."""
    state = _aircraft_state()
    downwash = calculate_downwash(state, 6.0)
    load = calculate_tail_load(state, 6.0, downwash)
    assert load['alpt_deg'] == pytest.approx(6.0 - downwash['eps_deg'])
    assert load['cla_tail_sref'] == pytest.approx(
        load['cla_tail'] * load['tail_area'] / 135.0)
    assert load['cl_tail'] == pytest.approx(
        load['cla_tail_sref'] * load['alpt_deg'])


def test_tail_moment_is_no_longer_zero():
    """PHYSICS_REVIEW item 11: the tail must contribute to aircraft Cm."""
    state = _aircraft_state()
    result = calculate_total_pitching_moment(state, 0.5, 5.0, 0.3)
    assert result['tail_supported']
    assert result['tail_method'] == 'legacy_clwbt_with_carryover'
    assert abs(result['cm_tail']) > 1e-6


def test_tail_is_stabilizing():
    """dCm/dalpha must be more negative with the tail than without it."""
    state = _aircraft_state()
    no_tail = {k: v for k, v in state.items() if not k.startswith('htail_')}

    def slope(config):
        low = calculate_total_pitching_moment(config, 0.30, 2.0, 0.3)
        high = calculate_total_pitching_moment(config, 0.45, 8.0, 0.3)
        return (high['cm_total'] - low['cm_total']) / 6.0

    assert not calculate_total_pitching_moment(
        no_tail, 0.5, 5.0, 0.3)['tail_supported']
    assert slope(state) < slope(no_tail)


def test_tail_moment_grows_with_tail_arm():
    """A longer arm produces a larger restoring moment."""
    near = _aircraft_state()
    near['synths_xh'] = 28.0
    far = _aircraft_state()
    far['synths_xh'] = 45.0
    cm_near = calculate_total_pitching_moment(near, 0.5, 5.0, 0.3)['cm_tail']
    cm_far = calculate_total_pitching_moment(far, 0.5, 5.0, 0.3)['cm_tail']
    assert cm_far < cm_near < 0.0


def test_supplied_tail_lift_still_wins():
    """An upstream translated tail load overrides the internal buildup."""
    state = _aircraft_state()
    state['aero_cl_tail'] = 0.0
    result = calculate_total_pitching_moment(state, 0.5, 5.0, 0.3)
    assert result['tail_method'] == 'supplied_cl_tail'
    assert result['cm_tail'] == pytest.approx(0.0)


def test_incomplete_tail_geometry_is_reported_not_guessed():
    """Without tail geometry the result must say so rather than invent a load."""
    state = _aircraft_state()
    del state['htail_chrdr']
    result = calculate_total_pitching_moment(state, 0.5, 5.0, 0.3)
    assert not result['tail_supported']
    assert result['tail_method'] == 'unsupported_configuration'
    assert result['cm_tail'] == 0.0
