"""
Regression tests for the CNCA translation.

Source of truth: datcom-legacy/datcom_2000/cnca.f.

The axis rotation is shared with pydatcom.aerodynamics.moment, so these
cover the derivative pass, the UNUSED masking, and the source's asymmetric
first-angle guard.
"""

import numpy as np
import pytest

from pydatcom.aerodynamics.cnca import calculate_cnca
from pydatcom.aerodynamics.moment import (
    calculate_normal_force_coefficient, calculate_axial_force_coefficient,
)
from pydatcom.utils.constants import UNUSED

_ALPHA = np.array([-4., -2., 0., 2., 4., 8.])
_SLOPE = 0.08
_CM_SLOPE = -0.01


def _case():
    cl = _SLOPE * _ALPHA + 0.05
    cd = 0.02 + 0.05 * cl**2
    cm = _CM_SLOPE * _ALPHA
    return cl, cd, cm


def test_rotation_agrees_with_the_existing_moment_helpers():
    """CNCA's rotation must match the one already in moment.py."""
    cl, cd, _ = _case()
    result = calculate_cnca(_ALPHA, cl, cd)
    for index, alpha in enumerate(_ALPHA):
        assert result['cn'][index] == pytest.approx(
            calculate_normal_force_coefficient(cl[index], cd[index], alpha))
        assert result['ca'][index] == pytest.approx(
            calculate_axial_force_coefficient(cl[index], cd[index], alpha))


def test_lift_slope_recovers_a_linear_curve():
    """A linear CL schedule must differentiate back to its own slope."""
    cl, cd, _ = _case()
    result = calculate_cnca(_ALPHA, cl, cd)
    assert np.allclose(result['cla'][1:], _SLOPE, atol=1e-9)


def test_moment_slope_recovers_a_linear_curve():
    cl, cd, cm = _case()
    result = calculate_cnca(_ALPHA, cl, cd, cm)
    assert np.allclose(result['cma'], _CM_SLOPE, atol=1e-9)


def test_first_angle_has_no_lift_slope():
    """The source guards its CLA call with J >= 2."""
    cl, cd, _ = _case()
    result = calculate_cnca(_ALPHA, cl, cd)
    assert np.isnan(result['cla'][0])
    assert np.all(np.isfinite(result['cla'][1:]))


def test_moment_slope_carries_no_first_angle_guard():
    """CMA is evaluated at every angle; the asymmetry is in the source."""
    cl, cd, cm = _case()
    result = calculate_cnca(_ALPHA, cl, cd, cm)
    assert np.all(np.isfinite(result['cma']))


def test_unused_entries_are_masked_not_computed():
    """A sentinel CL leaves CN and CA absent at that angle only."""
    cl, cd, _ = _case()
    cl = cl.copy()
    cl[2] = UNUSED
    result = calculate_cnca(_ALPHA, cl, cd)
    assert np.isnan(result['cn'][2]) and np.isnan(result['ca'][2])
    others = [0, 1, 3, 4, 5]
    assert np.all(np.isfinite(result['cn'][others]))


def test_unused_drag_also_masks_the_pair():
    cl, cd, _ = _case()
    cd = cd.copy()
    cd[4] = UNUSED
    result = calculate_cnca(_ALPHA, cl, cd)
    assert np.isnan(result['cn'][4]) and np.isnan(result['ca'][4])


def test_zero_alpha_leaves_the_coefficients_unrotated():
    result = calculate_cnca([0.0, 2.0], [0.5, 0.6], [0.02, 0.03])
    assert result['cn'][0] == pytest.approx(0.5)
    assert result['ca'][0] == pytest.approx(0.02)


def test_moment_is_optional():
    cl, cd, _ = _case()
    assert 'cma' not in calculate_cnca(_ALPHA, cl, cd)


def test_rejects_mismatched_input():
    cl, cd, cm = _case()
    with pytest.raises(ValueError):
        calculate_cnca(_ALPHA, cl[:3], cd)
    with pytest.raises(ValueError):
        calculate_cnca([], [], [])
    with pytest.raises(ValueError):
        calculate_cnca(_ALPHA, cl, cd, cm[:3])
