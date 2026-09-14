"""
Regression tests for the SUBPAW / SUBPAH translation.

Source of truth: datcom-legacy/datcom_2000/subpaw.f, subpah.f.

The two source routines are identical apart from whitespace and COMMON
offsets, so one translation and one test file cover both the wing and the
horizontal tail.
"""

import numpy as np
import pytest

from pydatcom.aerodynamics.subpaw import (
    calculate_subpaw, _FIG_71410_6_X, _FIG_71410_6_Y,
    _FIG_71420_8_X, _FIG_71420_8_Y,
    _COMPRESSIBILITY_MACH, _BETA_ASPECT_LIMIT,
)

_SWEEP = np.deg2rad(25.0)


def _run(mach=0.3, **kwargs):
    params = dict(cla=0.075, section_cla_compressible=0.10,
                  cos_sweep_c4=np.cos(_SWEEP), tan_sweep_c4=np.tan(_SWEEP),
                  aspect_ratio=4.0, area=120.0, sref=135.0,
                  mac=4.5, cbarr=4.6667, dcm_dcl=-0.10, dcm_dclq=-0.30,
                  xac_root_fraction=0.25, dxcg=-0.5, taper_ratio=0.0)
    params.update(kwargs)
    return calculate_subpaw(mach, **params)


# --------------------------------------------------------------------------
# Figures
# --------------------------------------------------------------------------

def test_figure_table_shapes():
    assert len(_FIG_71410_6_X) == len(_FIG_71410_6_Y) == 8
    assert len(_FIG_71420_8_X) == len(_FIG_71420_8_Y) == 9


def test_figures_start_at_zero_and_rise():
    for x, y in ((_FIG_71410_6_X, _FIG_71410_6_Y),
                 (_FIG_71420_8_X, _FIG_71420_8_Y)):
        assert x[0] == 0.0 and y[0] == 0.0
        assert all(b > a for a, b in zip(y, y[1:]))
        assert x[-1] == _BETA_ASPECT_LIMIT


# --------------------------------------------------------------------------
# Pitching derivatives
# --------------------------------------------------------------------------

def test_clq_matches_the_source_expression():
    result = _run()
    assert result['clq'] == pytest.approx(
        (0.5 + 2.0 * -0.30) * 0.075 * 4.5 / 4.6667)


def test_clq_is_independent_of_mach():
    values = [_run(m)['clq'] for m in (0.1, 0.5, 0.9)]
    assert all(v == pytest.approx(values[0]) for v in values)


def test_compressibility_correction_applies_only_above_mach_02():
    at_limit = _run(_COMPRESSIBILITY_MACH)
    above = _run(_COMPRESSIBILITY_MACH + 0.01)
    assert at_limit['cmq'] == pytest.approx(at_limit['cmq_incompressible'])
    assert above['cmq'] != pytest.approx(above['cmq_incompressible'])


def test_pitch_damping_grows_in_magnitude_with_mach():
    values = [abs(_run(m)['cmq']) for m in (0.3, 0.6, 0.9)]
    assert all(b > a for a, b in zip(values, values[1:]))


def test_cmq_is_negative_for_a_stable_surface():
    assert _run()['cmq'] < 0.0


# --------------------------------------------------------------------------
# Acceleration derivative gates
# --------------------------------------------------------------------------

def test_supersonic_blocks_the_acceleration_derivatives():
    # A highly swept surface keeps M*cos below one so the CMQ guard passes.
    sweep = np.deg2rad(70.0)
    result = _run(1.2, cos_sweep_c4=np.cos(sweep), tan_sweep_c4=np.tan(sweep))
    assert not result['acceleration_available']
    assert result['acceleration_gate'] == 'supersonic'
    assert 'clad' not in result


def test_taper_blocks_the_acceleration_derivatives():
    """The source computes them only for an untapered planform."""
    result = _run(taper_ratio=0.4)
    assert not result['acceleration_available']
    assert result['acceleration_gate'] == 'tapered_surface'


def test_high_beta_aspect_ratio_blocks_the_acceleration_derivatives():
    result = _run(0.3, aspect_ratio=8.0)
    assert not result['acceleration_available']
    assert result['acceleration_gate'] == 'beta_aspect_ratio_above_four'


def test_pitching_derivatives_survive_every_gate():
    """CLQ and CMQ are produced even when the acceleration path is blocked."""
    for kwargs in (dict(taper_ratio=0.4), dict(aspect_ratio=8.0)):
        result = _run(0.3, **kwargs)
        assert np.isfinite(result['clq']) and np.isfinite(result['cmq'])


# --------------------------------------------------------------------------
# Acceleration derivatives
# --------------------------------------------------------------------------

def test_acceleration_derivatives_are_produced_when_allowed():
    result = _run(0.3)
    assert result['acceleration_available']
    assert np.isfinite(result['clad']) and np.isfinite(result['cmad'])


def test_cmad_carries_the_cg_shift():
    """CMAD = CMADPP + (DXCG/CBARR)*CLAD."""
    result = _run(0.3, dxcg=-0.5)
    assert result['cmad'] == pytest.approx(
        result['cmad_before_cg_shift'] + (-0.5 / 4.6667) * result['clad'])


def test_zero_cg_offset_leaves_cmad_unshifted():
    result = _run(0.3, dxcg=0.0)
    assert result['cmad'] == pytest.approx(result['cmad_before_cg_shift'])


def test_acceleration_derivatives_grow_with_mach():
    values = [abs(_run(m)['clad']) for m in (0.1, 0.5, 0.9)]
    assert all(b > a for a, b in zip(values, values[1:]))


# --------------------------------------------------------------------------
# Validation
# --------------------------------------------------------------------------

def test_fast_swept_case_is_rejected_rather_than_returning_nan():
    """The source's sqrt(1 - M^2 cos^2) runs before its supersonic test."""
    with pytest.raises(ValueError, match=r"M\*cos"):
        _run(1.2)


def test_rejects_bad_references():
    with pytest.raises(ValueError):
        _run(sref=0.0)
    with pytest.raises(ValueError):
        _run(cbarr=0.0)
