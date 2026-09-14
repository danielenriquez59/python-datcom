"""GRDEFF translation tests.

The figure tables are checked against ``grdeff.f`` itself: the test parses
the FORTRAN DATA statements with its own reader and compares, so a mistyped
constant cannot pass by agreeing with a mistyped expectation.  A wrong
ordinate in a ground-effect figure produces a plausible number with no
error anywhere, which is exactly the failure a self-consistent test misses.
"""

import re
import pathlib

import numpy as np
import pytest

from pydatcom.aerodynamics.ground_effect import (
    calculate_grdeff, ground_effect_geometry, ground_effect_incidence,
    figure_4711_14, figure_4711_15, figure_4711_17, figure_4711_18a,
    figure_4711_21,
    _F14_HWOB2, _F14_DXOB2, _F14_X, _F15_HWCOCR, _F15_CLOCOS, _F15_LOLOM1,
    _F18A_INVERSE_TAPER, _F18A_ASPECT_RATIO, _F18A_BWOB,
    _F21_HWOCBR, _F21_CL, _F21_BW,
    _F17_HMAC_OVER_CR, _F17_CURVE_A, _F17_CURVE_B,
)

_SOURCE = (pathlib.Path(__file__).resolve().parents[1] /
           'datcom-legacy' / 'datcom_2000' / 'grdeff.f')


def _logical_lines(text):
    """Join fixed-form FORTRAN continuation lines into logical statements."""
    joined = []
    for raw in text.splitlines():
        if not raw.strip() or raw[0] in 'Cc*':
            continue
        body = raw[6:] if len(raw) > 6 else ''
        continued = len(raw) > 5 and raw[5] not in (' ', '0')
        if continued and joined:
            joined[-1] += body
        else:
            joined.append(body)
    return joined


def _expand(field):
    """Expand the FORTRAN repeat form ``n*value`` into individual values."""
    values = []
    for item in field.split(','):
        item = item.strip()
        if not item:
            continue
        if '*' in item:
            count, _, value = item.partition('*')
            values.extend([float(value)] * int(count))
        else:
            values.append(float(item))
    return np.array(values)


def _parse_data(name):
    """Read one named DATA array out of grdeff.f.

    Spaces are removed first, so the pattern anchors on either ``DATA<name>/``
    or ``,<name>/`` -- a DATA statement may declare several arrays in turn.
    """
    text = _SOURCE.read_text(errors='replace')
    pattern = re.compile(r'(?:DATA|,)' + re.escape(name) + r'/([^/]*)/')
    for line in _logical_lines(text):
        squeezed = line.strip().replace(' ', '')
        if not squeezed.startswith('DATA'):
            continue
        match = pattern.search(squeezed)
        if match:
            return _expand(match.group(1))
    raise AssertionError(f"DATA {name} not found in grdeff.f")


# --------------------------------------------------------------------------
# The figure tables, against the source
# --------------------------------------------------------------------------

@pytest.mark.parametrize('name,expected', [
    ('X218', _F14_HWOB2), ('X118', _F14_DXOB2),
    ('X219', _F15_HWCOCR), ('X119', _F15_CLOCOS),
    ('X222A', _F18A_INVERSE_TAPER), ('X122A', _F18A_ASPECT_RATIO),
    ('X225', _F21_HWOCBR), ('X125', _F21_CL),
    ('X4717A', _F17_HMAC_OVER_CR),
    ('Y4717A', _F17_CURVE_A), ('Y4717B', _F17_CURVE_B),
])
def test_grid_matches_source(name, expected):
    np.testing.assert_allclose(_parse_data(name), expected, rtol=0, atol=0)


@pytest.mark.parametrize('name,table,shape', [
    ('Y18', _F14_X, (11, 7)),
    ('Y19', _F15_LOLOM1, (12, 9)),
    ('Y22A', _F18A_BWOB, (6, 4)),
    ('Y25', _F21_BW, (11, 9)),
])
def test_dependent_table_matches_source(name, table, shape):
    """Values and FORTRAN Y(NX2,NX1) column-major layout both checked."""
    flat = _parse_data(name)
    assert flat.size == shape[0] * shape[1], (
        f"{name} holds {flat.size} values, not {shape[0]}x{shape[1]}")
    np.testing.assert_allclose(table, flat.reshape(shape, order='F'),
                               rtol=0, atol=0)


def test_every_figure_is_rectangular():
    """Each dependent table exactly fills its two grids -- no short DATA."""
    for dependent, first, second in (
            (_F14_X, _F14_HWOB2, _F14_DXOB2),
            (_F15_LOLOM1, _F15_HWCOCR, _F15_CLOCOS),
            (_F18A_BWOB, _F18A_INVERSE_TAPER, _F18A_ASPECT_RATIO),
            (_F21_BW, _F21_HWOCBR, _F21_CL)):
        assert dependent.shape == (first.size, second.size)
    for curve in (_F17_CURVE_A, _F17_CURVE_B):
        assert curve.shape == _F17_HMAC_OVER_CR.shape


# --------------------------------------------------------------------------
# Figure lookups land on their tabulated nodes
# --------------------------------------------------------------------------

def test_figures_return_exact_nodes():
    # Figure 4.7.1-14 at DX/(b/2)=0.0 (X118[3]) and h/(b/2)=0.5 (X218[5]).
    assert figure_4711_14(0.0, 0.5) == pytest.approx(_F14_X[5, 3])
    assert figure_4711_15(10.0, 0.6) == pytest.approx(_F15_LOLOM1[2, 2])
    assert figure_4711_18a(6.0, 2.0) == pytest.approx(_F18A_BWOB[2, 1])
    assert figure_4711_21(0.6, 0.5) == pytest.approx(_F21_BW[3, 3])
    assert figure_4711_17(3, 0.5) == pytest.approx(_F17_CURVE_A[5])
    assert figure_4711_17(5, 0.5) == pytest.approx(_F17_CURVE_B[5])


def test_figure_4711_17_only_applies_to_flap_types_three_to_five():
    """The source's IF(IFTYPE.LT.3 .OR. IFTYPE.GT.5) leaves DDCLF zero."""
    for kind in (0, 1, 2, 6, 7, 8):
        assert figure_4711_17(kind, 0.5) == 0.0
    assert figure_4711_17(4, 0.5) == figure_4711_17(3, 0.5)
    assert figure_4711_17(5, 0.5) != figure_4711_17(3, 0.5)


def test_figure_4711_14_decays_with_height():
    """X falls monotonically as the wing climbs, at every DX station."""
    for dxob2 in _F14_DXOB2:
        values = [figure_4711_14(float(dxob2), float(h))
                  for h in np.linspace(0.0, 1.0, 11)]
        assert all(a >= b for a, b in zip(values, values[1:]))


# --------------------------------------------------------------------------
# Geometry
# --------------------------------------------------------------------------

def _wing(**overrides):
    base = dict(sspn=15.0, sspnop=0.0, chrdr=6.0, sspndd=0.0,
                dhdadi=3.0, dhdado=None, planform_type=1.0, tovc=0.12)
    base.update(overrides)
    return base


_SYNTH = dict(xcg=20.0, zw=2.0, aliw=1.0, xh=45.0, zh=4.0, alih=0.0)
_SWEEP = dict(tan_c4=0.20, cos_c4=0.98, tan_c4_inboard=0.20,
              tan_c4_outboard=0.25, cos_c4_inboard=0.98,
              cos_c4_outboard=0.97, zero_sweep_station=0.5)
_WT = dict(aspect_ratio=7.0, taper_ratio=0.5, mac=6.2, y_mac=6.4, mac_c4=3.1)


def test_straight_planform_dx_uses_quarter_chord_tangent():
    geometry = ground_effect_geometry(5.0, _wing(), _SYNTH, _SWEEP, _WT)
    assert geometry['dx'] == pytest.approx(0.5 * 6.0 - 0.20 * 0.75 * 15.0)
    assert geometry['dxob2'] == pytest.approx(geometry['dx'] / 15.0)
    assert geometry['dx_uses_zero_sweep_station'] is False


def test_cranked_planform_splits_dx_across_both_panels():
    """With SSPNOP above a quarter semispan both panel tangents contribute."""
    wing = _wing(planform_type=3.0, sspnop=6.0)
    geometry = ground_effect_geometry(5.0, wing, _SYNTH, _SWEEP, _WT)
    expected = (0.5 * 6.0 - 0.25 * (6.0 - 0.25 * 15.0)
                - 0.20 * (15.0 - 6.0))
    assert geometry['dx'] == pytest.approx(expected)
    assert geometry['dx_uses_zero_sweep_station'] is False


def test_small_outboard_panel_dx_multiplies_span_by_a_chord_fraction():
    """Source defect: A(8) is a chord station, used where a tangent belongs.

    INFTGM sets A(8) under "DETERMINE LOCATION OF ZERO SWEEP ANGLE", giving
    it values like 1.0, .25*CP and A(110)/A(13)+CHSTAT -- all chord
    fractions.  The two sibling DX branches multiply the same span by a
    sweep tangent.  The source form is preserved; the flag records it and
    the result reports what the parallel branches would have given.
    """
    wing = _wing(planform_type=3.0, sspnop=2.0)
    geometry = ground_effect_geometry(5.0, wing, _SYNTH, _SWEEP, _WT)
    assert geometry['dx_uses_zero_sweep_station'] is True
    assert geometry['dx'] == pytest.approx(0.5 * 6.0 - 0.75 * 15.0 * 0.5)
    # What A(68), the inboard tangent, would have produced instead.
    assert geometry['dx_with_inboard_tangent'] == pytest.approx(
        0.5 * 6.0 - 0.75 * 15.0 * 0.20)
    assert geometry['dx'] != pytest.approx(
        geometry['dx_with_inboard_tangent'])


@pytest.mark.parametrize('dhdadi,dhdado,sspndd', [
    (None, None, 0.0),        # neither dihedral given
    (3.0, None, 0.0),         # inboard only
    (None, 5.0, 12.0),        # outboard only, past the quarter semispan
    (None, 5.0, 2.0),         # outboard only, inside the quarter semispan
    (3.0, 5.0, 12.0),         # both, past the quarter semispan
    (3.0, 5.0, 2.0),          # both, inside the quarter semispan
])
def test_every_dihedral_branch_is_reachable_and_finite(dhdadi, dhdado, sspndd):
    wing = _wing(dhdadi=dhdadi, dhdado=dhdado, sspndd=sspndd)
    geometry = ground_effect_geometry(5.0, wing, _SYNTH, _SWEEP, _WT)
    assert np.isfinite(geometry['hw'])
    assert np.isfinite(geometry['hwmac4'])
    # HW sits at or above the 75% root chord elevation for positive dihedral.
    assert geometry['hw'] >= geometry['h75cr'] - 1e-9


def test_dihedral_raises_the_wing_above_the_flat_case():
    flat = ground_effect_geometry(5.0, _wing(dhdadi=None), _SYNTH, _SWEEP,
                                  _WT)
    dihedral = ground_effect_geometry(5.0, _wing(dhdadi=6.0), _SYNTH, _SWEEP,
                                      _WT)
    assert dihedral['hw'] > flat['hw']


def test_height_parameters_follow_the_source_definitions():
    geometry = ground_effect_geometry(5.0, _wing(), _SYNTH, _SWEEP, _WT)
    hwob2 = geometry['hwob2']
    assert geometry['r'] == pytest.approx((1.0 + hwob2**2)**0.5 - hwob2)
    assert geometry['sigma'] == pytest.approx(np.exp(-2.48 * hwob2**0.768))
    assert geometry['hwocbr'] == pytest.approx(geometry['hw'] / 6.2)


def test_ground_proximity_raises_sigma_and_r():
    """Both fall off with height; that is what makes the effect vanish."""
    low = ground_effect_geometry(1.0, _wing(), _SYNTH, _SWEEP, _WT)
    high = ground_effect_geometry(40.0, _wing(), _SYNTH, _SWEEP, _WT)
    assert low['sigma'] > high['sigma']
    assert low['r'] > high['r']


def test_hwmac4_is_computed_even_without_a_tail():
    """Divergence: the source leaves HWMAC4 unset when HTPL is false.

    Figure 4.7.1-17 reads HWMAC4 with no tail test, so without a tail the
    source reads GR(9) -- zero on the first case, the previous case's value
    afterwards.  The formula uses only wing quantities, so it is evaluated
    unconditionally here and agrees wherever the source defines it.
    """
    without = ground_effect_geometry(5.0, _wing(), _SYNTH, _SWEEP, _WT)
    assert without['hwmac4_computed_without_tail'] is True
    assert np.isfinite(without['hwmac4'])

    tail = dict(sspn=6.0, sspndd=0.0, dhdadi=0.0, dhdado=None)
    with_tail = ground_effect_geometry(
        5.0, _wing(), _SYNTH, _SWEEP, _WT, tail=tail,
        tail_theoretical=dict(y_mac=2.6, mac_c4=1.3))
    assert with_tail['hwmac4_computed_without_tail'] is False
    # The tail branch changes nothing about the wing MAC elevation.
    assert with_tail['hwmac4'] == pytest.approx(without['hwmac4'])
    assert np.isfinite(with_tail['htmac4'])


def test_nonpositive_geometry_is_rejected():
    with pytest.raises(ValueError, match="root chord"):
        ground_effect_geometry(5.0, _wing(chrdr=0.0), _SYNTH, _SWEEP, _WT)
    with pytest.raises(ValueError, match="MAC"):
        ground_effect_geometry(5.0, _wing(), _SYNTH, _SWEEP,
                               dict(_WT, mac=0.0))


# --------------------------------------------------------------------------
# The full buildup
# --------------------------------------------------------------------------

_ALPHA = np.array([-2., 0., 2., 4., 6., 8., 10.])


def _case(aspect_ratio=7.0, ground_height=5.0, with_tail=False, **kwargs):
    cl = 0.08 * _ALPHA + 0.15
    cd = 0.02 + cl**2 / (np.pi * aspect_ratio * 0.85)
    cm = -0.05 - 0.01 * _ALPHA
    cl_tail = cl + 0.03 * _ALPHA
    blocks = {}
    for name in ('wing_body', 'wing_body_vertical', 'wing_body_tail',
                 'wing_body_tail_vertical'):
        tailish = name in ('wing_body_tail', 'wing_body_tail_vertical')
        blocks[name] = dict(cd=cd.copy(),
                            cl=(cl_tail if tailish else cl).copy(),
                            cm=cm.copy(),
                            cla=np.full_like(_ALPHA, 0.08),
                            cma=np.full_like(_ALPHA, -0.01))
    extra = {}
    if with_tail:
        extra = dict(
            tail=dict(sspn=6.0, sspndd=0.0, dhdadi=0.0, dhdado=None),
            tail_theoretical=dict(y_mac=2.6, mac_c4=1.3),
            downwash=dict(slope=np.full_like(_ALPHA, 0.35),
                          angle=0.35 * _ALPHA + 0.5))
    extra.update(kwargs)
    return calculate_grdeff(
        ground_height, _ALPHA, _wing(), _SYNTH, _SWEEP,
        dict(_WT, aspect_ratio=aspect_ratio),
        dict(cd=cd, cl=cl, flap_dcl=0.0), blocks, **extra)


def test_ground_effect_raises_lift_and_cuts_induced_drag():
    """The defining physics: closer to the ground, more lift, less drag."""
    result = _case(has_vertical_panel=True)
    positive = _ALPHA > 0
    assert np.all(result['dclwbg'][positive] > 0.0)
    assert np.all(result['dcdlwg'][positive] < 0.0)
    # The incidence increment is what drives it, and it is a reduction.
    assert np.all(result['incidence']['dalpha'][positive] < 0.0)


def test_effect_weakens_as_the_aircraft_climbs():
    near = _case(ground_height=2.0, has_vertical_panel=True)
    far = _case(ground_height=30.0, has_vertical_panel=True)
    positive = _ALPHA > 0
    assert np.all(near['dclwbg'][positive] > far['dclwbg'][positive])
    assert np.all(np.abs(near['dcdlwg'][positive]) >
                  np.abs(far['dcdlwg'][positive]))


def test_low_and_high_aspect_ratio_paths_are_selected_at_three():
    assert _case(aspect_ratio=3.0)['incidence']['path'] == 'high_aspect_ratio'
    assert _case(aspect_ratio=2.99)['incidence']['path'] == 'low_aspect_ratio'
    low = _case(aspect_ratio=2.0)['incidence']
    assert 'k' in low and 'bw' in low
    assert low['bw'].shape == _ALPHA.shape


def test_moment_increment_follows_the_free_air_centre_of_pressure():
    """DXCP is BWI(121)/BWI(101), taken once at the first angle."""
    result = _case(has_vertical_panel=True)
    assert result['dxcp'] == pytest.approx(-0.01 / 0.08)
    np.testing.assert_allclose(result['dcmwbg'],
                               result['dxcp'] * result['dclwbg'])


def test_normal_and_axial_force_are_consistent_with_lift_and_drag():
    result = _case(has_vertical_panel=True)
    block = result['buildup']['wing_body']
    sin_a, cos_a = np.sin(np.deg2rad(_ALPHA)), np.cos(np.deg2rad(_ALPHA))
    # Rotating CN/CA back through alpha must return CL/CD.
    np.testing.assert_allclose(block['cn'] * cos_a - block['ca'] * sin_a,
                               block['cl'], rtol=1e-6, atol=1e-9)
    np.testing.assert_allclose(block['cn'] * sin_a + block['ca'] * cos_a,
                               block['cd'], rtol=1e-6, atol=1e-9)


def test_tail_path_reduces_downwash_near_the_ground():
    """The image-vortex ratio is below one, so the tail sees less downwash."""
    result = _case(with_tail=True, has_vertical_panel=True)
    assert result['has_horizontal_tail'] is True
    assert 0.0 < result['tail']['ratio'] < 1.0
    supplied = 0.35 * _ALPHA + 0.5
    np.testing.assert_allclose(result['tail']['ddwash'],
                               supplied * result['tail']['ratio'])
    # Less downwash means a higher effective tail angle.
    positive = supplied > 0
    assert np.all(result['tail']['alphat'][positive] >
                  (_ALPHA - supplied)[positive])


def test_tail_lift_without_a_vertical_panel_is_flagged():
    """CLG is built on BWV(J+20), which the restore block never refreshes."""
    with_panel = _case(with_tail=True, has_vertical_panel=True)
    without = _case(with_tail=True, has_vertical_panel=False)
    assert with_panel['tail_lift_read_stale_vertical'] is False
    assert without['tail_lift_read_stale_vertical'] is True


def test_free_air_inputs_are_never_mutated():
    """The source's CIOM save/restore exists to guarantee this."""
    cl = 0.08 * _ALPHA + 0.15
    cd = 0.02 + cl**2 / (np.pi * 7.0 * 0.85)
    cm = -0.05 - 0.01 * _ALPHA
    blocks = {name: dict(cd=cd.copy(), cl=cl.copy(), cm=cm.copy(),
                         cla=np.full_like(_ALPHA, 0.08),
                         cma=np.full_like(_ALPHA, -0.01))
              for name in ('wing_body', 'wing_body_vertical',
                           'wing_body_tail', 'wing_body_tail_vertical')}
    snapshot = {name: {key: value.copy() for key, value in block.items()}
                for name, block in blocks.items()}
    calculate_grdeff(5.0, _ALPHA, _wing(), _SYNTH, _SWEEP, _WT,
                     dict(cd=cd, cl=cl, flap_dcl=0.0), blocks,
                     has_vertical_panel=True)
    for name, block in blocks.items():
        for key, value in block.items():
            np.testing.assert_array_equal(value, snapshot[name][key])


def test_repeated_heights_are_independent():
    """Two heights run in either order give the same answer at each."""
    first = _case(ground_height=3.0, has_vertical_panel=True)
    second = _case(ground_height=12.0, has_vertical_panel=True)
    again = _case(ground_height=3.0, has_vertical_panel=True)
    np.testing.assert_allclose(first['dclwbg'], again['dclwbg'])
    assert not np.allclose(first['dclwbg'], second['dclwbg'])


def test_flap_increment_only_reaches_the_wing_without_a_tail():
    """The source adds WING(L+200) to CLWF only when HTPL is false."""
    geometry = ground_effect_geometry(5.0, _wing(), _SYNTH, _SWEEP, _WT)
    cl = 0.08 * _ALPHA + 0.15
    common = dict(geometry=geometry, alpha_deg=_ALPHA, wing=_wing(),
                  sweep=_SWEEP, wing_theoretical=_WT,
                  wing_body=dict(cl=cl, cla=np.full_like(_ALPHA, 0.08)))
    with_flap = ground_effect_incidence(
        wing_alone=dict(cl=cl, flap_dcl=0.4), has_horizontal_tail=False,
        **common)
    without_flap = ground_effect_incidence(
        wing_alone=dict(cl=cl, flap_dcl=0.0), has_horizontal_tail=False,
        **common)
    assert not np.allclose(with_flap['dalpha'], without_flap['dalpha'])

    tailed_with = ground_effect_incidence(
        wing_alone=dict(cl=cl, flap_dcl=0.4), has_horizontal_tail=True,
        **common)
    tailed_without = ground_effect_incidence(
        wing_alone=dict(cl=cl, flap_dcl=0.0), has_horizontal_tail=True,
        **common)
    np.testing.assert_allclose(tailed_with['dalpha'],
                               tailed_without['dalpha'])


def test_zero_taper_is_rejected_by_the_tail_path():
    """Figure 4.7.1-18A is entered on 1/A(118), so taper cannot be zero."""
    cl = 0.08 * _ALPHA + 0.15
    blocks = {name: dict(cd=np.zeros_like(_ALPHA), cl=cl.copy(),
                         cm=np.zeros_like(_ALPHA),
                         cla=np.full_like(_ALPHA, 0.08),
                         cma=np.full_like(_ALPHA, -0.01))
              for name in ('wing_body', 'wing_body_vertical',
                           'wing_body_tail', 'wing_body_tail_vertical')}
    with pytest.raises(ValueError, match="taper"):
        calculate_grdeff(
            5.0, _ALPHA, _wing(), _SYNTH, _SWEEP,
            dict(_WT, taper_ratio=0.0),
            dict(cd=np.zeros_like(_ALPHA), cl=cl), blocks,
            tail=dict(sspn=6.0, sspndd=0.0, dhdadi=0.0, dhdado=None),
            tail_theoretical=dict(y_mac=2.6, mac_c4=1.3),
            downwash=dict(slope=np.full_like(_ALPHA, 0.35),
                          angle=0.35 * _ALPHA + 0.5))


def test_mismatched_schedule_is_rejected():
    cl = 0.08 * _ALPHA + 0.15
    blocks = {name: dict(cd=cl[:3], cl=cl[:3], cm=cl[:3],
                         cla=cl[:3], cma=cl[:3])
              for name in ('wing_body', 'wing_body_vertical',
                           'wing_body_tail', 'wing_body_tail_vertical')}
    with pytest.raises(ValueError, match="match the angle schedule"):
        calculate_grdeff(5.0, _ALPHA, _wing(), _SYNTH, _SWEEP, _WT,
                         dict(cd=cl, cl=cl), blocks)
