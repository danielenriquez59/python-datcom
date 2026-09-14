"""
Regression tests for the SDDVC translation.

Source of truth: datcom-legacy/datcom_2000/sddvc.f.

The figure lookups go through the already-tested INTERX, so these cover the
table shapes SDDVC owns, the ICASE branch selection, the case-3 blend, and
the taper interpolation that combines the four figure parts.
"""

import numpy as np
import pytest

from pydatcom.aerodynamics.vortex_core import (
    calculate_sddvc, _TAPER_GRID, _XA, _YA, _XB, _YB1, _YB2,
    _XC, _YC1, _YC2, _XD, _YD,
)


def _run(**kwargs):
    params = dict(x=[1.5, 3.0], beta_aspect=4.0, taper_ratio=0.4,
                  icase=2, sweep_te=0.3, sweep_reference=1.0)
    params.update(kwargs)
    return calculate_sddvc(**params)


# --------------------------------------------------------------------------
# Table shapes: the LIND and LDEP of the source INTERX calls
# --------------------------------------------------------------------------

def test_figure_part_shapes_match_the_source_calls():
    assert len(_XA) == 7 and len(_YA) == 10        # LGH (5,2), LIND 5
    assert len(_XB) == 12 and len(_YB1) == 32      # LGH (8,4), LIND 8
    assert len(_YB2) == 32
    assert len(_XC) == 8 and len(_YC1) == 15       # LGH (5,3), LIND 5
    assert len(_YC2) == 15
    assert len(_XD) == 9 and len(_YD) == 20        # LGH (5,4), LIND 5


def test_taper_grid_matches_the_four_figure_parts():
    """74A is taper 0, 74B/74C cover 0.25 and 0.50, 74D is taper 1."""
    assert list(_TAPER_GRID) == [0.0, 0.25, 0.50, 1.0]


def test_figure_74b_curves_cross_at_their_flat_tail():
    """YB2 sits below YB1 almost everywhere but crosses it at the tail.

    YB1 flattens at 2.51 while YB2 reaches 2.52, so a blanket "YB2 is
    always lower" assertion would be asserting an assumption rather than
    the source data.
    """
    crossings = [i for i, (a, b) in enumerate(zip(_YB1, _YB2)) if b > a]
    assert crossings == [14, 15]
    assert _YB1[14] == 2.51 and _YB2[14] == 2.52


# --------------------------------------------------------------------------
# Taper interpolation
# --------------------------------------------------------------------------

@pytest.mark.parametrize("index,taper", list(enumerate(_TAPER_GRID)))
def test_taper_grid_points_return_their_figure_value(index, taper):
    """At a grid taper the result is exactly that part's value."""
    result = _run(taper_ratio=float(taper))
    for station in range(2):
        assert result['dhb'][station] == pytest.approx(
            result['per_taper'][station][index])


def test_displacement_falls_with_taper_ratio():
    """A more rectangular wing sheds its vortex closer to the chord plane."""
    values = [_run(taper_ratio=t)['dhb'][0] for t in (0.25, 0.5, 1.0)]
    assert all(b < a for a, b in zip(values, values[1:]))


# --------------------------------------------------------------------------
# ICASE branches
# --------------------------------------------------------------------------

def test_case_two_uses_figure_74b():
    """Case 2 keeps the half-chord-unswept values untouched."""
    result = _run(icase=2)
    assert result['icase'] == 2
    assert np.all(np.isfinite(result['per_taper']))


def test_case_three_blends_between_the_two_figures():
    """Case 3 lies between the case 1 and case 2 results."""
    case1 = _run(icase=1)['dhb']
    case2 = _run(icase=2)['dhb']
    case3 = _run(icase=3)['dhb']
    for a, b, c in zip(case1, case2, case3):
        assert min(a, b) <= c <= max(a, b)


def test_case_three_with_zero_trailing_edge_sweep_matches_case_two():
    """The blend weight is SWEPTE/SWEPR, so zero sweep keeps the 74B value."""
    blended = _run(icase=3, sweep_te=0.0)['dhb']
    assert blended == pytest.approx(_run(icase=2)['dhb'])


def test_case_three_with_full_sweep_ratio_matches_case_one():
    """A unit weight moves the mid-taper values all the way to 74C."""
    blended = _run(icase=3, sweep_te=1.0, sweep_reference=1.0)['dhb']
    assert blended == pytest.approx(_run(icase=1)['dhb'])


# --------------------------------------------------------------------------
# Behaviour and validation
# --------------------------------------------------------------------------

def test_displacement_grows_downstream():
    result = _run(x=[1.0, 6.0], taper_ratio=0.25)
    assert result['dhb'][1] > result['dhb'][0]


def test_both_stations_are_evaluated_independently():
    same = _run(x=[2.0, 2.0])
    assert same['dhb'][0] == pytest.approx(same['dhb'][1])


def test_rejects_bad_input():
    with pytest.raises(ValueError, match="two streamwise"):
        _run(x=[1.0, 2.0, 3.0])
    with pytest.raises(ValueError, match="ICASE"):
        _run(icase=4)
    with pytest.raises(ValueError, match="reference sweep"):
        _run(icase=3, sweep_reference=0.0)
