"""
Regression tests for the CALCA0 translation.

Source of truth: datcom-legacy/datcom_2000/calca0.f.

Covers the shared bracket search, the DA0OT column compression, the twist
and camber corrections, and the gates that skip each of them.
"""

import numpy as np
import pytest

from pydatcom.aerodynamics.calca0 import (
    calculate_calca0, _bracket, _twist_correction, _camber_factor,
    _DA0OT, _TAPER_GRID, _ASPECT_GRID, _SWEEP_GRID,
    _CAMBER_TOC, _CAMBER_NPT, _CAMBER_LOCX, _CAMBER_LOCY,
    _CAMBER_CX, _CAMBER_CY, _EXACT_TOLERANCE,
)


def _run(**kwargs):
    params = dict(section_alpha_zero=-2.0, taper_ratio=0.5,
                  aspect_ratio=6.0, sweep_c4_deg=20.0,
                  cos_sweep_c4=np.cos(np.deg2rad(20.0)), mach=0.4,
                  twist_deg=-3.0, thickness_ratio=0.12, camber=False)
    params.update(kwargs)
    return calculate_calca0(**params)


# --------------------------------------------------------------------------
# Tables
# --------------------------------------------------------------------------

def test_da0ot_shape_and_corners():
    """DA0OT(22,10): 22 sweeps by 10 columns."""
    assert _DA0OT.shape == (22, 10)
    assert _DA0OT[0, 0] == pytest.approx(-.399)
    assert _DA0OT[21, 0] == pytest.approx(-.372)
    assert _DA0OT[0, 6] == pytest.approx(-.419)
    assert _DA0OT[21, 9] == pytest.approx(-.425)


def test_da0ot_columns_overlap_between_taper_groups():
    """The index is 3*(taper-1)+aspect, so 3 groups of 4 share two columns.

    Twelve taper/aspect combinations address only ten columns: column 4 is
    both (taper 1, aspect 4) and (taper 2, aspect 1), and column 7 likewise.
    That compression is in the source, not introduced by the translation.
    """
    combinations = [3 * (taper - 1) + aspect
                    for taper in (1, 2, 3) for aspect in (1, 2, 3, 4)]
    assert len(combinations) == 12
    assert max(combinations) == 10
    assert sorted(set(combinations)) == list(range(1, 11))
    assert combinations.count(4) == 2 and combinations.count(7) == 2


def test_camber_curve_offsets_cover_every_value():
    """NPT, LOCX and LOCY slice CX and CY into seven ragged curves."""
    assert sum(_CAMBER_NPT) == len(_CAMBER_CY) == 33
    assert len(_CAMBER_CX) == 12
    for position in range(7):
        count = _CAMBER_NPT[position]
        x_start = _CAMBER_LOCX[position] - 1
        y_start = _CAMBER_LOCY[position] - 1
        assert x_start + count <= len(_CAMBER_CX)
        assert y_start + count <= len(_CAMBER_CY)


def test_camber_thickness_grid_descends():
    """TOC runs 16 down to 7, so the bracket search reverses its test."""
    assert list(_CAMBER_TOC) == sorted(_CAMBER_TOC, reverse=True)


# --------------------------------------------------------------------------
# The shared bracket search
# --------------------------------------------------------------------------

def test_bracket_exact_hit_skips_interpolation():
    index, skip, fraction = _bracket(_TAPER_GRID, 0.5)
    assert (index, skip, fraction) == (2, True, 0.0)


def test_bracket_within_tolerance_counts_as_exact():
    """A query within 2e-2 of a grid point is treated as landing on it."""
    _, skip, _ = _bracket(_ASPECT_GRID, 6.0 + _EXACT_TOLERANCE / 2.0)
    assert skip


def test_bracket_interior_gives_the_standard_fraction():
    index, skip, fraction = _bracket(_TAPER_GRID, 0.3)
    assert index == 2 and not skip
    assert fraction == pytest.approx((0.3 - 0.0) / (0.5 - 0.0))


def test_bracket_first_point_skips_interpolation():
    _, skip, _ = _bracket(_TAPER_GRID, 0.0)
    assert skip


def test_bracket_above_the_grid_clamps():
    """The loop runs to completion and falls into the exact-hit label."""
    index, skip, _ = _bracket(_SWEEP_GRID, 70.0)
    assert index == len(_SWEEP_GRID) and skip


def test_bracket_handles_a_descending_grid():
    index, skip, fraction = _bracket(_CAMBER_TOC, 13.0, descending=True)
    assert index == 3 and not skip
    assert fraction == pytest.approx((13.0 - 14.0) / (12.0 - 14.0))


# --------------------------------------------------------------------------
# Twist correction
# --------------------------------------------------------------------------

def test_twist_below_half_a_degree_is_skipped():
    """The source tests ABS(TWISTA) against 0.5 before doing any work."""
    assert _run(twist_deg=0.4)['twist_term'] is None
    assert _run(twist_deg=-0.4)['twist_term'] is None
    assert _run(twist_deg=0.6)['twist_term'] is not None


def test_twist_enters_linearly():
    """alpha0 = TWISTA * DA0OT + section alpha0."""
    result = _run(twist_deg=-3.0)
    assert result['alpha_zero_lift'] == pytest.approx(
        -3.0 * result['twist_term'] + -2.0)


def test_twist_term_is_negative_across_the_table():
    """Every DA0OT entry is negative, so the correction has one sign."""
    assert np.all(_DA0OT < 0.0)
    for taper in (0.0, 0.25, 0.5, 1.0):
        for aspect in (1.5, 4.0, 10.0):
            assert _twist_correction(taper, aspect, 10.0)['value'] < 0.0


def test_low_aspect_ratio_column_clamp_is_unreachable():
    """The source's IF(IAR.EQ.4 .AND. A(7).LE.0.5) can never fire.

    The bracket search returns index 4 only for an aspect ratio above 6.0,
    which cannot simultaneously be at or below 0.5. The guard is carried in
    the translation for fidelity; this pins its dead-ness so its absence
    from the behaviour is not mistaken for an omission.
    """
    for aspect in (0.1, 0.4, 0.5, 1.0, 3.0):
        index, _, _ = _bracket(_ASPECT_GRID, aspect)
        assert index < 4
    for aspect in (6.5, 9.0, 10.0, 50.0):
        index, _, _ = _bracket(_ASPECT_GRID, aspect)
        assert index == 4 and aspect > 0.5


# --------------------------------------------------------------------------
# Camber correction
# --------------------------------------------------------------------------

def test_camber_is_skipped_when_not_requested():
    assert _run(camber=False)['camber_factor'] is None


def test_camber_scales_the_zero_lift_angle():
    result = _run(camber=True, twist_deg=0.0)
    assert result['alpha_zero_lift'] == pytest.approx(
        result['camber_factor'] * -2.0)


def test_camber_factor_falls_as_mach_rises():
    """The curves run from 1.0 down through zero as cos(sweep)*Mach grows."""
    values = [_camber_factor(12.0, 1.0, m)['factor']
              for m in (0.4, 0.6, 0.7)]
    assert all(b <= a for a, b in zip(values, values[1:]))


def test_camber_thickness_selects_different_curves():
    thin = _camber_factor(7.0, 1.0, 0.8)
    thick = _camber_factor(16.0, 1.0, 0.8)
    assert thin['thickness_index'] != thick['thickness_index']


def test_both_corrections_compose():
    result = _run(twist_deg=-3.0, camber=True)
    expected = result['camber_factor'] * (
        -3.0 * result['twist_term'] + -2.0)
    assert result['alpha_zero_lift'] == pytest.approx(expected)
