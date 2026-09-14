"""
Regression tests for the MAXCL translation.

Source of truth: datcom-legacy/datcom_2000/maxcl.f.

This is the first translated caller of INTERX's three-variable branch, so
the DEL1 lookup exercises the TLIN3X memory layout that was settled against
a compiled probe.
"""

import numpy as np
import pytest

from pydatcom.aerodynamics.maxcl import (
    calculate_maxcl, _PARM58, _EVAL58, _PARM59, _EVAL59,
    _PARM60, _EVAL60, _PARM61, _EVAL61, _PARM62, _EVAL62,
)


def _run(**kwargs):
    params = dict(delta_y=2.5, thickness_station=0.35,
                  reynolds_per_length=[2.0e6], chord=4.0)
    params.update(kwargs)
    return calculate_maxcl(**params)


# --------------------------------------------------------------------------
# Table shapes
# --------------------------------------------------------------------------

def test_table_shapes_match_the_source_calls():
    assert len(_PARM58) == 26 and len(_EVAL58) == 52   # LGH (13,4), LIND 13
    assert len(_PARM59) == 30 and len(_EVAL59) == 160  # LGH (10,4,4), LIND 10
    assert len(_PARM60) == 20 and len(_EVAL60) == 30   # LGH (10,3), LIND 10
    assert len(_PARM61) == 16 and len(_EVAL61) == 32   # LGH (8,4), LIND 8
    assert len(_PARM62) == 8 and len(_EVAL62) == 8


def test_three_variable_grid_offsets():
    """PARM59 packs three grids at stride 10, so the third starts at 20."""
    assert _PARM59[:10] == [0., 1., 1.5, 2., 2.5, 3., 3.5, 4., 4.5, 5.]
    assert _PARM59[10:14] == [0., .02, .04, .06]
    assert _PARM59[20:24] == [.15, .30, .40, .50]


@pytest.mark.parametrize("delta_y,station,expected", [
    (0.0, 0.30, 0.8), (3.0, 0.30, 1.58), (5.0, 0.30, 1.42),
    (0.0, 0.45, 0.8), (5.0, 0.45, 1.35), (3.0, 0.45, 1.35),
])
def test_base_clmax_corners(delta_y, station, expected):
    """Values readable straight from the source EVAL58 DATA."""
    assert _run(delta_y=delta_y, thickness_station=station)['clbase'] == \
        pytest.approx(expected, abs=1e-12)


# --------------------------------------------------------------------------
# Increments and their gates
# --------------------------------------------------------------------------

def test_uncambered_section_has_no_camber_increments():
    result = _run(camber=False)
    assert result['del1'] == 0.0
    assert result['del2'] == 0.0


def test_camber_enables_both_camber_increments():
    result = _run(camber=True, camber_ratio=0.04, camber_station=0.30)
    assert result['del1'] != 0.0
    assert result['del2'] != 0.0


def test_thickness_station_of_030_skips_del2():
    """The source tests XOVC against 0.30 exactly before computing DEL2."""
    result = _run(thickness_station=0.30, camber=True,
                  camber_ratio=0.04, camber_station=0.30)
    assert result['del1'] != 0.0
    assert result['del2'] == 0.0


def test_near_030_station_is_snapped_before_the_test():
    """A station within 1e-5 of 0.30 is snapped, so DEL2 is also skipped."""
    result = _run(thickness_station=0.300001, camber=True,
                  camber_ratio=0.04, camber_station=0.30)
    assert result['thickness_station_used'] == 0.30
    assert result['del2'] == 0.0


def test_low_reynolds_number_substitutes_a_nominal_value():
    """The source uses 9.0e6 when the section Reynolds number is below one."""
    tiny = _run(reynolds_per_length=[1.0e-9])['del3'][0]
    nominal = _run(reynolds_per_length=[9.0e6 / 4.0])['del3'][0]
    assert tiny == pytest.approx(nominal)


def test_clmax_is_the_sum_of_base_and_increments():
    result = _run(camber=True, camber_ratio=0.04, camber_station=0.30)
    assert result['clmax'][0] == pytest.approx(
        result['clbase'] + result['del1'] + result['del2'] +
        result['del3'][0])


def test_del4_is_computed_but_excluded():
    """The source evaluates DEL4 and never adds it; the dead path is kept."""
    result = _run(delta_y=2.5)
    assert result['del4_unused'] == pytest.approx(-0.46)
    assert result['clmax'][0] != pytest.approx(
        result['clbase'] + result['del3'][0] + result['del4_unused'])


# --------------------------------------------------------------------------
# Schedule and validation
# --------------------------------------------------------------------------

def test_one_result_per_flight_condition():
    result = _run(reynolds_per_length=[1.0e6, 2.0e6, 4.0e6])
    assert result['clmax'].shape == (3,)
    assert result['del3'].shape == (3,)


def test_clmax_is_physically_plausible():
    for delta_y in (1.0, 2.5, 4.0):
        value = _run(delta_y=delta_y)['clmax'][0]
        assert 0.5 < value < 2.0


def test_rejects_bad_input():
    with pytest.raises(ValueError):
        _run(chord=0.0)
    with pytest.raises(ValueError):
        _run(reynolds_per_length=[])
