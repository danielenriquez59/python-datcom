"""CDRAG translation tests.

As with the other table-heavy routines, the figure constants are checked
against ``cdrag.f`` by a parser this test carries itself, so a mistyped
ordinate cannot pass by agreeing with a mistyped expectation.  Two of the
tables are assembled in the source through EQUIVALENCE -- ``Y42`` from four
quarters and ``Y48`` from three thirds -- and that assembly is checked too.
"""

import re
import pathlib

import numpy as np
import pytest

from pydatcom.aerodynamics.cdrag import (
    calculate_cdrag, STRAIGHT_TAPERED, CRANKED, DOUBLE_DELTA, CURVED,
    _F28B_COS_SWEEP, _F28B_MACH, _F28B_RLS,
    _F28BD_COS_SWEEP, _F28BD_MACH, _F28BD_RLS,
    _F42_ANGLE, _F42_DRAG_FACTOR, _F42_TAPER, _F42,
    _F48_ANGLE, _F48_DRAG_FACTOR, _F48_TAPER, _F48,
    _F53B_X, _F53B_Y, _F54_ASPECT_RATIO, _F54_CL_OVER_AR, _F54,
    _F27_MACH, _F27_INTERCEPT,
    _F4312_10A_DIAMETER_RATIO, _F4312_10A_KWB,
    _suction_parameter,
)

_SOURCE = (pathlib.Path(__file__).resolve().parents[1] /
           'datcom-legacy' / 'datcom_2000' / 'cdrag.f')


def _parse_all():
    """Every DATA array in cdrag.f, by name."""
    text = _SOURCE.read_text(errors='replace')
    lines = []
    for raw in text.splitlines():
        if not raw.strip() or raw[0] in 'Cc*':
            continue
        body = raw[6:] if len(raw) > 6 else ''
        if len(raw) > 5 and raw[5] not in (' ', '0') and lines:
            lines[-1] += body
        else:
            lines.append(body)
    tables = {}
    for line in lines:
        squeezed = line.strip().replace(' ', '')
        if not squeezed.startswith('DATA'):
            continue
        for match in re.finditer(r'(?:DATA|,)([A-Z0-9]+)/([^/]*)/', squeezed):
            values = []
            try:
                for item in match.group(2).split(','):
                    item = item.strip()
                    if not item:
                        continue
                    if '*' in item:
                        count, _, value = item.partition('*')
                        values.extend([float(value)] * int(count))
                    else:
                        values.append(float(item))
            except ValueError:
                continue
            tables[match.group(1)] = np.array(values)
    return tables


_TABLES = _parse_all()


@pytest.mark.parametrize('name,expected', [
    ('X228B', _F28B_COS_SWEEP), ('X128B', _F28B_MACH),
    ('X228BD', _F28BD_COS_SWEEP), ('X128BD', _F28BD_MACH),
    ('X142', _F42_ANGLE), ('X242', _F42_DRAG_FACTOR), ('X342', _F42_TAPER),
    ('X148', _F48_ANGLE), ('X248', _F48_DRAG_FACTOR), ('X348', _F48_TAPER),
    ('X53B', _F53B_X), ('Y53B', _F53B_Y),
    ('X254', _F54_ASPECT_RATIO), ('X154', _F54_CL_OVER_AR),
    ('X27M', _F27_MACH), ('X27I', _F27_INTERCEPT),
    ('X1OA', _F4312_10A_DIAMETER_RATIO), ('Y1OA', _F4312_10A_KWB),
])
def test_grid_matches_source(name, expected):
    np.testing.assert_allclose(_TABLES[name], expected, rtol=0, atol=0)


@pytest.mark.parametrize('names,table,shape', [
    (('Y28B',), _F28B_RLS, (11, 4)),
    (('Y28BD',), _F28BD_RLS, (9, 4)),
    (('Y42A', 'Y42B', 'Y42C', 'Y42D'), _F42, (12, 4, 9)),
    (('Y48A', 'Y48B', 'Y48C'), _F48, (10, 4, 9)),
    (('Y54',), _F54, (14, 12)),
])
def test_dependent_table_matches_source(names, table, shape):
    """Values, the EQUIVALENCE assembly, and the column-major layout."""
    flat = np.concatenate([_TABLES[name] for name in names])
    assert flat.size == int(np.prod(shape)), (
        f"{'+'.join(names)} holds {flat.size} values, not {shape}")
    np.testing.assert_allclose(table, flat.reshape(shape, order='F'),
                               rtol=0, atol=0)


def test_every_figure_is_rectangular():
    assert _F28B_RLS.shape == (_F28B_COS_SWEEP.size, _F28B_MACH.size)
    assert _F28BD_RLS.shape == (_F28BD_COS_SWEEP.size, _F28BD_MACH.size)
    assert _F42.shape == (_F42_DRAG_FACTOR.size, _F42_ANGLE.size,
                          _F42_TAPER.size)
    assert _F48.shape == (_F48_DRAG_FACTOR.size, _F48_ANGLE.size,
                          _F48_TAPER.size)
    assert _F54.shape == (_F54_ASPECT_RATIO.size, _F54_CL_OVER_AR.size)
    assert _F53B_X.shape == _F53B_Y.shape
    assert _F27_MACH.shape == _F27_INTERCEPT.shape
    assert _F4312_10A_DIAMETER_RATIO.shape == _F4312_10A_KWB.shape


# --------------------------------------------------------------------------
# A representative configuration
# --------------------------------------------------------------------------

_ALPHA = np.array([0., 2., 4., 6., 8., 10.])
_CL = 0.08 * _ALPHA + 0.05


def _geometry(**overrides):
    base = dict(area=180.0, area_inboard=110.0, area_outboard=70.0,
                mac=6.2, mac_inboard=7.4, mac_outboard=4.6,
                aspect_ratio=6.0, aspect_ratio_outboard=3.2,
                taper_ratio=0.45, taper_ratio_outboard=0.5,
                span_inboard=9.0, sspne=16.0, sspn=18.0)
    base.update(overrides)
    return base


def _section(**overrides):
    base = dict(tovc=0.12, xovc=0.35, tovco=0.10, xovco=0.28,
                leri=0.008, lero=0.006, twista=-2.0)
    base.update(overrides)
    return base


def _sweep(**overrides):
    base = dict(cos_le=0.90, tan_le=0.48, cos_le_inboard=0.88,
                tan_le_inboard=0.54, cos_le_outboard=0.93,
                tan_le_outboard=0.40, tan_c4=0.36,
                cos_max_thickness=0.94, cos_max_thickness_inboard=0.92,
                cos_max_thickness_outboard=0.96)
    base.update(overrides)
    return base


def _run(mach=0.6, beta=0.8, cl=None, **kwargs):
    return calculate_cdrag(
        mach, _ALPHA, _CL if cl is None else cl, 0.085,
        _geometry(**kwargs.pop('geometry', {})),
        _section(**kwargs.pop('section', {})),
        _sweep(**kwargs.pop('sweep', {})),
        300.0, 0.00016, 2.0e6, beta, **kwargs)


# --------------------------------------------------------------------------
# Structure and physics
# --------------------------------------------------------------------------

def test_total_drag_is_zero_lift_plus_lift_dependent():
    result = _run()
    np.testing.assert_allclose(result['cd'], result['cdo'] + result['cdl'])


def test_lift_dependent_drag_grows_with_lift():
    result = _run()
    assert np.all(np.diff(result['cdl']) > 0.0)
    assert result['cdo'] > 0.0


def test_lift_dependent_drag_is_quadratic_in_lift_on_the_straight_path():
    """With zero twist the only CL term left is the CL^2 induced one."""
    result = _run(section=dict(twista=0.0))
    ratio = result['cdl'][1:] / _CL[1:]**2
    np.testing.assert_allclose(ratio, ratio[0], rtol=1e-9)


def test_zero_lift_drag_rises_with_mach():
    values = [_run(mach=m, beta=b)['cdo']
              for m, b in ((0.3, 0.954), (0.6, 0.8), (0.85, 0.527))]
    assert values[0] < values[1] < values[2]


def test_thicker_section_raises_zero_lift_drag():
    thin = _run(section=dict(tovc=0.06))['cdo']
    thick = _run(section=dict(tovc=0.18))['cdo']
    assert thick > thin


def test_zero_lift_drag_scales_with_area_over_reference():
    """CDO carries a 2*S/SREF wetted-area factor."""
    base = _run()['cdo']
    doubled = _run(geometry=dict(area=360.0))['cdo']
    assert doubled == pytest.approx(2.0 * base)


def test_rougher_surface_caps_the_reynolds_number():
    """A rough surface lowers the Figure 4.1.5.1-27 cutoff."""
    smooth = calculate_cdrag(0.6, _ALPHA, _CL, 0.085, _geometry(),
                             _section(), _sweep(), 300.0, 1.0e-6, 2.0e6, 0.8)
    rough = calculate_cdrag(0.6, _ALPHA, _CL, 0.085, _geometry(),
                            _section(), _sweep(), 300.0, 1.0e-3, 2.0e6, 0.8)
    assert rough['cutoff_reynolds'] < smooth['cutoff_reynolds']
    assert rough['panels']['whole']['reynolds_was_capped'] is True
    assert rough['cdo'] > smooth['cdo']


# --------------------------------------------------------------------------
# Path selection
# --------------------------------------------------------------------------

def test_planform_type_selects_the_path():
    assert _run()['path'] == 'straight_tapered'
    assert _run(planform_type=CRANKED)['path'] == 'cranked'
    for kind in (DOUBLE_DELTA, CURVED):
        assert _run(planform_type=kind)['path'] == 'double_delta_or_curved'


def test_non_straight_planforms_split_zero_lift_drag_across_panels():
    straight = _run()
    cranked = _run(planform_type=CRANKED)
    assert set(straight['panels']) == {'whole'}
    assert set(cranked['panels']) == {'inboard', 'outboard'}
    assert cranked['cdo'] == pytest.approx(
        cranked['panels']['inboard']['cdo'] +
        cranked['panels']['outboard']['cdo'])


def test_inboard_panel_uses_the_dashed_correlation_curve():
    """The two Figure 4.1.5.1-28B variants differ below cos sweep 0.65."""
    low = _run(planform_type=CRANKED,
               sweep=dict(cos_max_thickness_inboard=0.50,
                          cos_max_thickness_outboard=0.50))
    inboard = low['panels']['inboard']['rls']
    outboard = low['panels']['outboard']['rls']
    # The dashed table starts at 0.45 and is flat to 0.65; the solid one
    # starts at 0.50 and is still climbing there.
    assert inboard != pytest.approx(outboard)


def test_double_delta_drag_is_the_vortex_term_in_local_angle():
    result = _run(planform_type=DOUBLE_DELTA)
    expected = 0.95 * _CL * np.tan(np.deg2rad(_ALPHA))
    np.testing.assert_allclose(result['cdl'], expected, rtol=1e-9)


def test_form_factor_coefficient_switches_at_thirty_percent_chord():
    forward = _run(section=dict(xovc=0.25))['panels']['whole']
    aft = _run(section=dict(xovc=0.35))['panels']['whole']
    assert forward['capl'] == 2.00
    assert aft['capl'] == 1.20
    # Exactly 0.30 takes the aft branch, matching the source's .GE. test.
    assert _run(section=dict(xovc=0.30))['panels']['whole']['capl'] == 1.20


# --------------------------------------------------------------------------
# Body carryover
# --------------------------------------------------------------------------

def test_body_carryover_scales_lift_before_the_drag_buildup():
    """Figure 4.3.1.2-10A is entered on (SSPN-SSPNE)/SSPN."""
    without = _run()
    with_body = _run(body_present=True)
    assert without['carryover'] == 1.0
    # d/b = (18-16)/18 = 0.1111, between the 0.1 and 0.2 table entries.
    assert with_body['carryover'] == pytest.approx(
        1.08 + (0.1111111 - 0.1) / 0.1 * (1.16 - 1.08), rel=1e-4)
    assert np.all(with_body['cdl'] >= without['cdl'])


def test_carryover_is_one_when_the_wing_reaches_the_centreline():
    result = _run(body_present=True,
                  geometry=dict(sspne=18.0, sspn=18.0))
    assert result['carryover'] == pytest.approx(1.0)


# --------------------------------------------------------------------------
# Source defects, pinned
# --------------------------------------------------------------------------

def test_inboard_suction_argument_is_hardcoded_zero_on_the_cranked_path():
    """Source: TEMPI=0.0 where the outboard forms A(6)*A(28)/A(85).

    The inboard analogue is never built, so the inboard leading-edge
    suction is always read at a zero abscissa.  Changing the inboard
    aspect ratio, taper or sweep therefore cannot move it, while the
    matching outboard quantities do move the outboard value.
    """
    base = _run(planform_type=CRANKED)
    assert base['inboard_suction_argument_is_zero'] is True

    # The inboard value is exactly what the suction curve gives at zero.
    reynolds_in = base['detail']['suction_inboard']
    assert reynolds_in == pytest.approx(
        _suction_parameter(
            2.0e6 * _section()['leri'] * _geometry()['mac_inboard'] /
            abs(_sweep()['tan_le_inboard']) *
            np.sqrt(1.0 - (0.6 * _sweep()['cos_le_inboard'])**2),
            0.0))

    # Had the inboard analogue A(5)*A(26)/A(61) been formed, it would land
    # somewhere else on the same curve.
    intended = _suction_parameter(
        2.0e6 * _section()['leri'] * _geometry()['mac_inboard'] /
        abs(_sweep()['tan_le_inboard']) *
        np.sqrt(1.0 - (0.6 * _sweep()['cos_le_inboard'])**2),
        4.0 * 0.6 / _sweep()['cos_le_inboard'])
    assert intended != pytest.approx(base['detail']['suction_inboard'])

    # The outboard argument genuinely responds to its own geometry, which
    # is what makes the inboard's insensitivity a defect rather than a
    # property of the curve.
    other = _run(planform_type=CRANKED,
                 geometry=dict(aspect_ratio_outboard=6.0))
    assert (other['detail']['suction_outboard'] !=
            pytest.approx(base['detail']['suction_outboard']))


def test_straight_path_oswald_carries_a_factor_the_cranked_one_does_not():
    """Source: D(30)=1.1*TEMP/(...) on the straight path only.

    The consequence is visible -- the straight path can report an "Oswald
    efficiency" above 1.0, which the identical cranked expression cannot.
    """
    straight = _run()
    cranked = _run(planform_type=CRANKED)
    assert straight['oswald_carries_source_1_1_factor'] is True
    assert cranked['oswald_carries_source_1_1_factor'] is False
    assert straight['oswald_efficiency'] > 1.0
    assert cranked['oswald_efficiency'] < 1.0
    # Dividing the 1.1 back out brings it into range.
    assert straight['oswald_efficiency'] / 1.1 < 1.0


def test_cutoff_reynolds_is_built_from_the_exposed_mac_for_every_panel():
    """Source: D(2) uses A(16) once and both panels reuse the result."""
    cranked = _run(planform_type=CRANKED)
    straight = _run()
    assert cranked['cutoff_reynolds'] == pytest.approx(
        straight['cutoff_reynolds'])
    # Changing a panel MAC cannot move the shared cutoff.
    moved = _run(planform_type=CRANKED, geometry=dict(mac_inboard=20.0))
    assert moved['cutoff_reynolds'] == pytest.approx(
        cranked['cutoff_reynolds'])


# --------------------------------------------------------------------------
# Rejections
# --------------------------------------------------------------------------

@pytest.mark.parametrize('kwargs', [
    dict(sref=0.0), dict(roughness=0.0), dict(area=0.0),
])
def test_nonpositive_scales_are_rejected(kwargs):
    sref = kwargs.get('sref', 300.0)
    roughness = kwargs.get('roughness', 0.00016)
    geometry = _geometry(area=kwargs['area']) if 'area' in kwargs \
        else _geometry()
    with pytest.raises(ValueError, match="positive"):
        calculate_cdrag(0.6, _ALPHA, _CL, 0.085, geometry, _section(),
                        _sweep(), sref, roughness, 2.0e6, 0.8)


def test_mismatched_lift_schedule_is_rejected():
    with pytest.raises(ValueError, match="matching the schedule"):
        calculate_cdrag(0.6, _ALPHA, _CL[:3], 0.085, _geometry(),
                        _section(), _sweep(), 300.0, 0.00016, 2.0e6, 0.8)
