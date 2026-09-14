"""
Regression tests for the WINGCL translation.

Source of truth: datcom-legacy/datcom_2000/wingcl.f.

Covers the transonic wing CL curve and the CLB interpolation. The CDL
section is not translated because its figure tables are four elements short
in the source; the guard that reports this is tested here.
"""

import numpy as np
import pytest

from pydatcom.aerodynamics.wingcl import (
    calculate_wingcl, calculate_wingcl_clb, calculate_wingcl_cdl,
    _CDL_REQUIRED, _CDL_SOURCE_LENGTH, _PARM,
)

_ALPHA = np.array([-4., 0., 4., 8., 10., 12., 14., 16.])


def _lift(**kwargs):
    params = dict(alpha_schedule=_ALPHA, alpha_zero_lift=-2.0,
                  alpha_stall_onset=10.0, alpha_clmax=14.0,
                  cla=0.08, clmax=1.25)
    params.update(kwargs)
    return calculate_wingcl(**params)


# --------------------------------------------------------------------------
# Transonic wing CL
# --------------------------------------------------------------------------

def test_lift_is_linear_below_stall_onset():
    """CL = CLA*(alpha - ALPHAO) up to and including ALPHAS."""
    result = _lift()
    for index in range(5):          # -4 through 10 degrees
        assert result['cl'][index] == pytest.approx(
            0.08 * (_ALPHA[index] + 2.0))


def test_lift_reaches_clmax_exactly_at_aclmax():
    result = _lift()
    assert result['cl'][6] == pytest.approx(1.25)


def test_lift_between_onset_and_clmax_is_above_neither_bound():
    result = _lift()
    assert 0.96 < result['cl'][5] < 1.25


def test_angles_beyond_clmax_are_left_unset():
    """The source breaks out of its loop past ACLMAX."""
    result = _lift()
    assert np.isnan(result['cl'][7])


def test_weak_clmax_is_corrected_up_to_the_linear_curve():
    """A stated CLMAX below the linear value at ALPHAS is raised."""
    result = _lift(clmax=0.5)
    assert result['stall_corrected']
    assert result['clmax'] == pytest.approx(0.08 * (10.0 + 2.0))
    assert result['alpha_clmax'] == pytest.approx(10.0)


def test_aclmax_below_onset_is_corrected():
    result = _lift(alpha_clmax=8.0)
    assert result['stall_corrected']
    assert result['alpha_clmax'] == pytest.approx(10.0)


def test_corrected_case_is_purely_linear():
    """With ACLMAX collapsed onto ALPHAS there is no nonlinear region."""
    result = _lift(clmax=0.5)
    assert not result['nonlinear']
    assert result['exponent'] is None


def test_supplied_lift_values_are_kept():
    supplied = [None, 0.99, None, None, None, None, None, None]
    result = _lift(cl_supplied=supplied)
    assert result['cl'][1] == pytest.approx(0.99)
    assert result['cl'][0] == pytest.approx(0.08 * (-4.0 + 2.0))


def test_lift_rejects_inconsistent_stall_angles():
    with pytest.raises(ValueError, match="zero-lift"):
        _lift(alpha_stall_onset=-3.0)
    with pytest.raises(ValueError):
        calculate_wingcl([], -2.0, 10.0, 14.0, 0.08, 1.25)


# --------------------------------------------------------------------------
# Transonic wing CLB
# --------------------------------------------------------------------------

def _clb(**kwargs):
    params = dict(cl=[0.2, 0.4, 0.6], mach=0.9, cla=0.08,
                  clb_subsonic=-0.0012, clb_supersonic=-0.0008,
                  cla_mach06=0.075, cla_mach14=0.06)
    params.update(kwargs)
    return calculate_wingcl_clb(**params)


def test_clb_is_proportional_to_lift():
    result = _clb()
    assert result['clb'] == pytest.approx(
        result['clb_per_cl'] * np.array([0.2, 0.4, 0.6]))


def test_clb_matches_the_source_equation():
    """Linear on Mach between anchors normalised by their own slopes."""
    subsonic = -0.0012 / 0.075**2
    supersonic = -0.0008 / 0.06**2
    expected = ((supersonic - subsonic) * (0.9 - 0.6) / 0.8 + subsonic) * 0.08**2
    assert _clb()['clb_per_cl'] == pytest.approx(expected, rel=1e-12)


def test_clb_reduces_to_each_anchor_at_its_own_mach():
    at_subsonic = _clb(mach=0.6)['clb_per_cl']
    assert at_subsonic == pytest.approx(-0.0012 / 0.075**2 * 0.08**2)
    at_supersonic = _clb(mach=1.4)['clb_per_cl']
    assert at_supersonic == pytest.approx(-0.0008 / 0.06**2 * 0.08**2)


def test_clb_skips_unset_lift_entries():
    result = _clb(cl=[0.2, np.nan, 0.6])
    assert np.isnan(result['clb'][1])
    assert np.all(np.isfinite(result['clb'][[0, 2]]))


def test_clb_rejects_zero_anchor_slopes():
    with pytest.raises(ValueError, match="anchor"):
        _clb(cla_mach06=0.0)


# --------------------------------------------------------------------------
# CDL: the defective source tables
# --------------------------------------------------------------------------

def test_cdl_grid_offsets_use_a_stride_of_seven():
    """PARM packs three grids at LIND=7, so they start at 0, 7 and 14."""
    assert _PARM[:7] == [-4., -3., -2., -1., 0., 1., 2.]
    assert _PARM[7:13] == [.5, .75, 1., 1.5, 1.75, 2.]
    assert _PARM[14:18] == [0., .2, .5, 1.]


def test_cdl_required_length_exceeds_the_source_arrays():
    """7*6*4 = 168, but DEP55A and DEP55B hold 164 each."""
    assert _CDL_REQUIRED == 168
    assert _CDL_SOURCE_LENGTH == 164
    assert _CDL_REQUIRED - _CDL_SOURCE_LENGTH == 4


def test_cdl_rejects_the_source_length_with_an_explanation():
    with pytest.raises(ValueError, match="168 elements"):
        calculate_wingcl_cdl(0.95, 0.06, 4.0, 0.5, 0.5,
                             [0.0] * 164, [0.0] * 164)


def test_cdl_accepts_a_completed_table():
    """With both tables completed the chain runs end to end."""
    table_a = list(np.linspace(0.2, 1.2, _CDL_REQUIRED))
    table_b = list(np.linspace(0.3, 1.3, _CDL_REQUIRED))
    result = calculate_wingcl_cdl(0.95, 0.06, 4.0, 0.5, 0.5,
                                  table_a, table_b)
    assert np.isfinite(result['cdl_per_cl2'])
    assert result['anchor_sweep0'] != result['anchor_sweep3']


def test_cdl_rejects_nonpositive_thickness():
    table = [0.5] * _CDL_REQUIRED
    with pytest.raises(ValueError, match="thickness"):
        calculate_wingcl_cdl(0.95, 0.0, 4.0, 0.5, 0.5, table, table)
