"""HINGE translation tests.

The figure constants are checked against ``hinge.f`` by a parser this test
carries itself.  The *shapes* are checked just as carefully: TLINEX declares
its dependent table ``Y(NX2,NX1)``, so a table reshaped in the other order
still holds every correct value and still interpolates without error, while
returning the wrong number.  Eleven of these eighteen figures are
non-square, so the orientation is checked against each call site's own
declared ``NX1``/``NX2``.
"""

import re
import pathlib

import numpy as np
import pytest

from pydatcom.aerodynamics import hinge as H
from pydatcom.aerodynamics.hinge import (
    calculate_hinge, ROUND_NOSE, ELLIPTIC_NOSE, SHARP_NOSE)

_SOURCE = (pathlib.Path(__file__).resolve().parents[1] /
           'datcom-legacy' / 'datcom_2000' / 'hinge.f')


def _parse_all():
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

# (python grid, FORTRAN name)
_GRIDS = [
    (H._F6111_39A_TC, 'X1125A'), (H._F6111_39A_CFOCA, 'X2125A'),
    (H._F6111_39B_CLOCLT, 'X1125B'), (H._F6111_39B_CFOCA, 'X2125B'),
    (H._F4112_8A_LOG_RF, 'X1128A'), (H._F4112_8A_TANPHP, 'X2128A'),
    (H._F6131_11A_TC, 'X1317A'), (H._F6131_11A_CFOCA, 'X2317A'),
    (H._F6131_11B_CLACLT, 'X1317B'), (H._F6131_11B_CFOCA, 'X2317B'),
    (H._F6131_12A_BALANCE, 'X1318A'),
    (H._F6131_12A_SHARP, 'Y318A1'), (H._F6131_12A_ELLIPTIC, 'Y318A2'),
    (H._F6131_12A_ROUND, 'Y318A3'),
    (H._F6132_12A_TC, 'X1327A'), (H._F6132_12A_CFOCA, 'X2327A'),
    (H._F6132_12B_CLACLT, 'X1327B'), (H._F6132_12B_CFOCA, 'X2327B'),
    (H._F6132_13A_BALANCE, 'X1328A'), (H._F6132_13A, 'Y1328A'),
    (H._F6132_13B_TC, 'X1328B'), (H._F6132_13B_BALANCE, 'X2328B'),
    (H._F6132_13C_TC, 'X1328C'), (H._F6132_13C_BALANCE, 'X2328C'),
    (H._F6161_19A_AR, 'X6115A'), (H._F6161_19A, 'Y6115A'),
    (H._F6161_19B_ETA, 'X6115B'), (H._F6161_19B, 'Y6115B'),
    (H._F6161_19C_CBOCF, 'X16116'), (H._F6161_19C_CFOCAP, 'X26116'),
    (H._F6162_15A_CFOCAP, 'X1629A'), (H._F6162_15A_AR, 'X2629A'),
    (H._F6162_15B_ETA, 'X1629B'), (H._F6162_15B, 'Y1629B'),
]

# (python table, FORTRAN name, first grid, second grid) -- the source's
# NX1 and NX2, in that order, at the routine's own TLINEX call site.
_SURFACES = [
    (H._F6111_39A, 'Y1125A', H._F6111_39A_TC, H._F6111_39A_CFOCA),
    (H._F6111_39B, 'Y1125B', H._F6111_39B_CLOCLT, H._F6111_39B_CFOCA),
    (H._F4112_8A, 'Y1128A', H._F4112_8A_LOG_RF, H._F4112_8A_TANPHP),
    (H._F6131_11A, 'Y1317A', H._F6131_11A_TC, H._F6131_11A_CFOCA),
    (H._F6131_11B, 'Y1317B', H._F6131_11B_CLACLT, H._F6131_11B_CFOCA),
    (H._F6132_12A, 'Y1327A', H._F6132_12A_TC, H._F6132_12A_CFOCA),
    (H._F6132_12B, 'Y1327B', H._F6132_12B_CLACLT, H._F6132_12B_CFOCA),
    (H._F6132_13B, 'Y1328B', H._F6132_13B_TC, H._F6132_13B_BALANCE),
    (H._F6132_13C, 'Y1328C', H._F6132_13C_TC, H._F6132_13C_BALANCE),
    (H._F6161_19C, 'Y16116', H._F6161_19C_CBOCF, H._F6161_19C_CFOCAP),
    (H._F6162_15A, 'Y1629A', H._F6162_15A_CFOCAP, H._F6162_15A_AR),
]


@pytest.mark.parametrize('grid,name', _GRIDS, ids=[n for _, n in _GRIDS])
def test_grid_matches_source(grid, name):
    np.testing.assert_allclose(_TABLES[name], grid, rtol=0, atol=0)


@pytest.mark.parametrize('table,name,first,second', _SURFACES,
                         ids=[s[1] for s in _SURFACES])
def test_surface_matches_source_in_the_tlinex_orientation(table, name, first,
                                                          second):
    """Values and the Y(NX2,NX1) layout TLINEX declares.

    A table reshaped the other way round still contains every correct
    value and still interpolates without error, so only comparing against
    the source in the right orientation catches the mistake.
    """
    flat = _TABLES[name]
    assert flat.size == first.size * second.size
    assert table.shape == (second.size, first.size), (
        f"{name} must be (NX2, NX1) = ({second.size}, {first.size})")
    np.testing.assert_allclose(
        table, flat.reshape((second.size, first.size), order='F'),
        rtol=0, atol=0)


def test_non_square_figures_would_catch_a_transposed_reshape():
    """Guards the test above: most of these tables are not square."""
    non_square = [name for _, name, first, second in _SURFACES
                  if first.size != second.size]
    assert len(non_square) >= 10


# --------------------------------------------------------------------------
# A representative flap
# --------------------------------------------------------------------------

def _surface(**overrides):
    base = dict(tovc=0.12, tovco=0.10, semispan=16.0, semispan_exposed=14.0,
                sspnop=0.0, aspect_ratio=6.0, sin_c4=0.34, cos_c4=0.94,
                sweep_c4_deg=20.0, sweep_le_deg=24.0, cos_le=0.913,
                tan_te=0.05, mac_exposed=6.2)
    base.update(overrides)
    return base


def _flap(**overrides):
    base = dict(span_inboard=4.0, span_outboard=11.0, chord_inboard=1.5,
                chord_outboard=1.1, chord_ratio=0.25, chord_balance=None,
                thickness_at_hinge=0.10, tan_te_angle=0.10,
                tan_te_angle_effective=0.12)
    base.update(overrides)
    return base


_DELTA = np.array([-10., -5., 5., 10., 15.])


def _run(mach=0.4, **kwargs):
    return calculate_hinge(
        mach, kwargs.pop('deflections', _DELTA),
        _surface(**kwargs.pop('surface', {})),
        _flap(**kwargs.pop('flap', {})),
        kwargs.pop('slope_ratio', 0.85), 2.0e6, 0.10,
        section_dcl=kwargs.pop('section_dcl', 0.045 * _DELTA), **kwargs)


# --------------------------------------------------------------------------
# Behaviour
# --------------------------------------------------------------------------

def test_both_derivatives_are_restoring():
    """A plain trailing-edge flap has negative CHA and CHD."""
    result = _run()
    assert result['cha'] < 0.0
    assert np.all(result['chd'] < 0.0)
    assert result['cha_section'] < 0.0
    assert result['chd_section'] < 0.0


def test_section_values_carry_the_prandtl_glauert_factor():
    """Both section derivatives are divided by sqrt(1-M^2)."""
    low = _run(mach=0.2)
    high = _run(mach=0.7)
    assert low['beta'] == pytest.approx(np.sqrt(1.0 - 0.2**2))
    assert abs(high['cha_section']) > abs(low['cha_section'])
    # The incompressible value behind the factor is the same.
    assert (low['cha_section'] * low['beta'] ==
            pytest.approx(high['cha_section'] * high['beta']))


def test_nose_balance_reduces_the_hinge_moment():
    """That is what an aerodynamic balance is for."""
    plain = _run()
    balanced = _run(flap=dict(chord_balance=0.45))
    assert plain['has_nose_balance'] is False
    assert balanced['has_nose_balance'] is True
    assert abs(balanced['cha']) < abs(plain['cha'])
    assert np.all(np.abs(balanced['chd']) < np.abs(plain['chd']))


def test_nose_shape_orders_the_balance_effectiveness():
    """Round balances most, sharp least, elliptic between."""
    factors = {}
    for kind in (SHARP_NOSE, ELLIPTIC_NOSE, ROUND_NOSE):
        result = _run(flap=dict(chord_balance=0.45), nose_type=kind)
        factors[kind] = result['alpha_nose_factor']
    assert factors[ROUND_NOSE] < factors[ELLIPTIC_NOSE] < factors[SHARP_NOSE]


def test_out_of_range_nose_type_falls_through_to_sharp():
    """A FORTRAN computed GO TO with an out-of-range index falls through.

    Both of this routine's computed GO TOs are followed immediately by
    their sharp-nose branch.
    """
    sharp = _run(flap=dict(chord_balance=0.45), nose_type=SHARP_NOSE)
    for kind in (0, 4, 99):
        other = _run(flap=dict(chord_balance=0.45), nose_type=kind)
        assert other['alpha_nose_factor'] == pytest.approx(
            sharp['alpha_nose_factor'])
        assert other['delta_nose_factor'] == pytest.approx(
            sharp['delta_nose_factor'])


def test_balance_ratio_is_dropped_when_the_radicand_goes_negative():
    """A balance chord smaller than half the hinge thickness has no ratio."""
    result = _run(flap=dict(chord_balance=0.01, thickness_at_hinge=0.90))
    assert result['has_nose_balance'] is False
    assert result['balance_ratio'] is None
    assert result['alpha_nose_factor'] is None


def test_thickness_switches_to_outboard_past_the_planform_break():
    """Source: IF(BIF.GE.(BO2-SSPNOP)) TC=TOVCO."""
    inboard = _run(surface=dict(sspnop=2.0), flap=dict(span_inboard=4.0))
    outboard = _run(surface=dict(sspnop=13.0), flap=dict(span_inboard=4.0))
    assert inboard['thickness_used'] == pytest.approx(0.12)
    assert inboard['thickness_is_outboard'] is False
    assert outboard['thickness_used'] == pytest.approx(0.10)
    assert outboard['thickness_is_outboard'] is True


def test_supplied_and_averaged_section_increments_agree():
    """SDCL used directly must match the four-per-deflection averages."""
    supplied = 0.045 * _DELTA
    alpha_delta = -supplied / (_DELTA * 0.10)
    # Four identical entries per deflection average back to the same value.
    dcl = np.repeat(supplied, 4)
    aldag = np.repeat(alpha_delta, 4)
    direct = _run(section_dcl=supplied)
    averaged = calculate_hinge(
        0.4, _DELTA, _surface(), _flap(), 0.85, 2.0e6, 0.10,
        section_dcl_table=dcl, alpha_delta_table=aldag)
    np.testing.assert_allclose(direct['chd'], averaged['chd'], rtol=1e-12)
    np.testing.assert_allclose(direct['section_increments'],
                               averaged['section_increments'], rtol=1e-12)


def test_span_weighting_uses_both_flap_ends():
    """K is read at the inboard and outboard stations and differenced."""
    narrow = _run(flap=dict(span_inboard=7.0, span_outboard=8.0))
    wide = _run(flap=dict(span_inboard=2.0, span_outboard=14.0))
    assert narrow['k_alpha'] != pytest.approx(wide['k_alpha'])
    assert narrow['k_delta'] != pytest.approx(wide['k_delta'])


def test_normal_chord_ratio_reduces_to_the_flap_chord_ratio_unswept():
    """With no sweep anywhere the normal conversion is the identity."""
    result = _run(surface=dict(sweep_c4_deg=0.0, sweep_le_deg=0.0,
                               cos_le=1.0, cos_c4=1.0, sin_c4=0.0,
                               tan_te=0.0),
                  flap=dict(chord_inboard=1.3, chord_outboard=1.3))
    assert result['sweep_hinge_deg'] == pytest.approx(0.0)
    assert result['normal_chord_ratio'] == pytest.approx(0.25)


# --------------------------------------------------------------------------
# Rejections
# --------------------------------------------------------------------------

@pytest.mark.parametrize('mach', [1.0, 1.5, np.nan])
def test_non_subsonic_mach_is_rejected(mach):
    with pytest.raises(ValueError, match="subsonic"):
        _run(mach=mach)


def test_zero_deflection_is_rejected():
    """The source divides DCHD by each deflection."""
    with pytest.raises(ValueError, match="deflection"):
        _run(deflections=np.array([-5.0, 0.0, 5.0]),
             section_dcl=np.array([-0.2, 0.0, 0.2]))


def test_degenerate_flap_span_is_rejected():
    with pytest.raises(ValueError, match="span"):
        _run(flap=dict(span_inboard=6.0, span_outboard=6.0))


def test_missing_section_increment_tables_are_rejected():
    with pytest.raises(ValueError, match="SDCL"):
        calculate_hinge(0.4, _DELTA, _surface(), _flap(), 0.85, 2.0e6, 0.10)


def test_short_increment_tables_are_rejected():
    with pytest.raises(ValueError, match="four"):
        calculate_hinge(0.4, _DELTA, _surface(), _flap(), 0.85, 2.0e6, 0.10,
                        section_dcl_table=np.zeros(4),
                        alpha_delta_table=np.zeros(4))
