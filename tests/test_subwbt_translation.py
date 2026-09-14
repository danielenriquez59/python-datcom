"""
Regression tests for the SUBWBT translation.

Source of truth: datcom-legacy/datcom_2000/subwbt.f.

Two source defects in the ventral fin block are reproduced deliberately;
these tests pin both so neither is mistaken for a translation slip, and so
that correcting them later is a visible decision.
"""

import numpy as np
import pytest

from pydatcom.aerodynamics.subwbt import calculate_subwbt

_ALPHA = np.array([0., 5., 10., 15.])
_VT = dict(arm_x=20.0, arm_z=6.0, cyb=-0.008)
_VF = dict(arm_x=18.0, arm_z=-3.0, cyb=-0.003)


def _run(**kwargs):
    zero = np.zeros_like(_ALPHA)
    params = dict(alpha_deg=_ALPHA, cyp_wing_body=zero,
                  cnp_wing_body=zero, cnr_wing_body=zero,
                  blref=30.0, vertical_tail=_VT, ventral_fin=_VF)
    params.update(kwargs)
    return calculate_subwbt(**params)


# --------------------------------------------------------------------------
# Panel geometry
# --------------------------------------------------------------------------

def test_vertical_tail_cyp_matches_the_source_expression():
    result = _run(ventral_fin=None)
    cos_a = np.cos(_ALPHA / 180.0 * np.pi)
    sin_a = np.sin(_ALPHA / 180.0 * np.pi)
    effective = 6.0 * cos_a - 20.0 * sin_a
    assert result['vertical_tail']['cyp'] == pytest.approx(
        2.0 * (effective - 6.0) * -0.008 / 30.0)


def test_vertical_tail_cyp_vanishes_at_zero_alpha():
    """At alpha=0 the effective arm equals ZP, so 2*(ZEE-ZP) is zero."""
    result = _run(ventral_fin=None)
    assert result['vertical_tail']['cyp'][0] == pytest.approx(0.0, abs=1e-15)


def test_yaw_damping_is_negative():
    """CNR must damp for panels aft of the CG."""
    assert np.all(_run()['cnr'] < 0.0)


def test_increments_add_to_the_wing_body_values():
    base = np.array([0.01, 0.02, 0.03, 0.04])
    result = _run(cyp_wing_body=base)
    assert result['cyp'] == pytest.approx(
        base + result['vertical_tail']['cyp'] + result['ventral_fin']['cyp'])


# --------------------------------------------------------------------------
# Source defect 1: the ventral CYP grouping
# --------------------------------------------------------------------------

def test_ventral_cyp_uses_the_source_ungrouped_form():
    """The source writes (2.*ZEE-ZPF) where the tail writes 2.*(ZEE-ZP).

    A missing pair of parentheses is the obvious reading, but correcting it
    would change results, so the source form stands.
    """
    result = _run(vertical_tail=None)
    cos_a = np.cos(_ALPHA / 180.0 * np.pi)
    sin_a = np.sin(_ALPHA / 180.0 * np.pi)
    effective = -3.0 * cos_a - 18.0 * sin_a
    ungrouped = (2.0 * effective - -3.0) * -0.003 / 30.0
    grouped = 2.0 * (effective - -3.0) * -0.003 / 30.0
    assert result['ventral_fin']['cyp'] == pytest.approx(ungrouped)
    assert not np.allclose(ungrouped, grouped)


def test_ventral_cyp_does_not_vanish_at_zero_alpha():
    """The consequence of the grouping defect, stated directly.

    The vertical tail's CYP is zero at alpha=0; the ventral fin's is not,
    because 2*ZEE - ZPF leaves ZPF behind where 2*(ZEE - ZPF) would cancel.
    """
    result = _run()
    assert result['vertical_tail']['cyp'][0] == pytest.approx(0.0, abs=1e-15)
    assert result['ventral_fin']['cyp'][0] == pytest.approx(
        -3.0 * -0.003 / 30.0)
    assert result['ventral_fin']['cyp'][0] != 0.0


# --------------------------------------------------------------------------
# Source defect 2: the ventral CNR sideslip derivative
# --------------------------------------------------------------------------

def test_ventral_cnr_uses_the_vertical_tail_sideslip_derivative():
    """The source multiplies the ventral CNR by DYBV, not DYBF.

    Every other ventral term uses DYBF. The same substitution appears in
    CLRDER, so it is a recurring transcription slip in the original.
    """
    result = _run()
    cos_a = np.cos(_ALPHA / 180.0 * np.pi)
    sin_a = np.sin(_ALPHA / 180.0 * np.pi)
    moment_arm = 18.0 * cos_a + -3.0 * sin_a
    with_vt_cyb = 2.0 * moment_arm**2 * _VT['cyb'] / 30.0**2
    with_vf_cyb = 2.0 * moment_arm**2 * _VF['cyb'] / 30.0**2
    assert result['ventral_fin']['cnr'] == pytest.approx(with_vt_cyb)
    assert not np.allclose(with_vt_cyb, with_vf_cyb)


def test_ventral_alone_falls_back_to_its_own_derivative():
    """With no vertical tail there is no DYBV to substitute."""
    result = _run(vertical_tail=None)
    assert result['ventral_cnr_fell_back_to_own_cyb']
    cos_a = np.cos(_ALPHA / 180.0 * np.pi)
    sin_a = np.sin(_ALPHA / 180.0 * np.pi)
    moment_arm = 18.0 * cos_a + -3.0 * sin_a
    assert result['ventral_fin']['cnr'] == pytest.approx(
        2.0 * moment_arm**2 * _VF['cyb'] / 30.0**2)


def test_substitution_flag_is_clear_when_both_panels_are_present():
    assert not _run()['ventral_cnr_fell_back_to_own_cyb']


# --------------------------------------------------------------------------
# Configuration bookkeeping
# --------------------------------------------------------------------------

def test_horizontal_tail_contributes_nothing():
    result = _run()
    for key in ('cyp', 'cnp', 'cnr'):
        assert np.allclose(result['horizontal_tail'][key], 0.0)


def test_wing_body_tail_equals_wing_body():
    base = np.array([0.01, 0.02, 0.03, 0.04])
    result = _run(cyp_wing_body=base)
    assert result['wing_body_tail']['cyp'] == pytest.approx(base)


def test_wing_body_vertical_equals_the_total():
    result = _run()
    for key in ('cyp', 'cnp', 'cnr'):
        assert result['wing_body_vertical'][key] == pytest.approx(result[key])


def test_no_panels_leaves_the_wing_body_values_untouched():
    base = np.array([0.01, 0.02, 0.03, 0.04])
    result = _run(vertical_tail=None, ventral_fin=None, cyp_wing_body=base)
    assert result['cyp'] == pytest.approx(base)


def test_rejects_bad_input():
    with pytest.raises(ValueError):
        _run(blref=0.0)
    with pytest.raises(ValueError):
        _run(cyp_wing_body=np.zeros(2))
