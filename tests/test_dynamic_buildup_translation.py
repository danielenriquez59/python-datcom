"""
Regression tests for the DNPAWB and DNPWBT translations.

Source of truth: datcom-legacy/datcom_2000/dnpawb.f, dnpwbt.f.

Both routines are pure combination, so these tests cover the combination
expressions, the regime and geometry gates the source applies, and the
branch split on wing-to-tail span ratio.
"""

import numpy as np
import pytest

from pydatcom.aerodynamics.dynamic_buildup import (
    calculate_dnpawb, calculate_dnpwbt,
)

_WING = dict(clq_wing=-8.0, cmq_wing=-12.0, clad_wing=-2.0, cmad_wing=-3.0)
_BODY = dict(clq_body=-0.5, cmq_body=-0.8, clad_body=-0.1, cmad_body=-0.2)
_K = dict(khb=1.08, kbh=0.13)


def _wing_body(**kwargs):
    params = dict(_WING, **_BODY)
    params.update(_K)
    params.setdefault('beta_aspect_ratio', 2.0)
    params.setdefault('taper_ratio', 0.0)
    params.update(kwargs)
    return calculate_dnpawb(**params)


def _wing_body_tail(**kwargs):
    base = _wing_body()
    params = dict(clq_wing_body=base['clq'], cmq_wing_body=base['cmq'],
                  clad_wing_body=base['clad'], cmad_wing_body=base['cmad'],
                  qoqi=[0.95, 0.93, 0.90], deda=[0.31, 0.31, 0.31],
                  cla_tail=0.06, khb=1.08, kbh=0.13,
                  dxac=-23.0, cbar=4.6667,
                  wing_span=30.0, tail_span=12.0)
    params.update(kwargs)
    return calculate_dnpwbt(**params)


# --------------------------------------------------------------------------
# DNPAWB
# --------------------------------------------------------------------------

def test_dnpawb_applies_carryover_to_the_wing_only():
    """CLQWB = (KWB+KBW)*CLQW + CLQB, body added unscaled."""
    result = _wing_body()
    carryover = 1.08 + 0.13
    assert result['clq'] == pytest.approx(carryover * -8.0 + -0.5)
    assert result['cmq'] == pytest.approx(carryover * -12.0 + -0.8)
    assert result['clad'] == pytest.approx(carryover * -2.0 + -0.1)
    assert result['cmad'] == pytest.approx(carryover * -3.0 + -0.2)


def test_dnpawb_hypersonic_zeroes_body_acceleration_terms():
    """The source sets CLADB and CMADB to zero in hypersonic flow."""
    result = _wing_body(hypersonic=True, supersonic=True)
    carryover = 1.08 + 0.13
    assert result['clad'] == pytest.approx(carryover * -2.0)
    assert result['cmad'] == pytest.approx(carryover * -3.0)


@pytest.mark.parametrize("beta_ar", [-0.5, 4.5])
def test_dnpawb_gates_acceleration_outside_beta_aspect_range(beta_ar):
    """Outside supersonic flow, BETA*AR must lie in [0, 4]."""
    result = _wing_body(beta_aspect_ratio=beta_ar)
    assert not result['acceleration_available']
    assert result['acceleration_gate'] == 'beta_aspect_ratio_out_of_range'
    assert 'clad' not in result
    # The pitching derivatives are still produced.
    assert 'clq' in result and 'cmq' in result


def test_dnpawb_gate_does_not_apply_supersonically():
    """The source skips both gates when SUPERS is set."""
    result = _wing_body(beta_aspect_ratio=9.0, taper_ratio=0.4,
                        supersonic=True)
    assert result['acceleration_available']


def test_dnpawb_gates_acceleration_for_a_tapered_wing():
    result = _wing_body(taper_ratio=0.4)
    assert not result['acceleration_available']
    assert result['acceleration_gate'] == 'tapered_wing'


def test_dnpawb_transonic_cmad_sentinel():
    """The source marks an unavailable transonic CMAD with a literal 1000."""
    result = _wing_body(cmad_wing=1000.0, transonic=True)
    assert not result['acceleration_available']
    assert result['acceleration_gate'] == 'transonic_cmad_unavailable'
    # CLAD is still produced before the source returns.
    assert 'clad' in result
    assert 'cmad' not in result


# --------------------------------------------------------------------------
# DNPWBT
# --------------------------------------------------------------------------

def test_dnpwbt_wide_wing_branch_expressions():
    """BW/BH >= 1.5 folds the arm into the increment before applying it."""
    result = _wing_body_tail()
    assert result['branch'] == 'wide_wing'
    base = _wing_body()
    arm = 23.0 / 4.6667
    qoqi = np.array([0.95, 0.93, 0.90])
    increment = 2.0 * (1.08 + 0.13) * qoqi * arm * 0.06
    assert result['clq'] == pytest.approx(base['clq'] + increment)
    assert result['cmq'] == pytest.approx(base['cmq'] - increment * arm)
    assert result['clad'] == pytest.approx(base['clad'] + increment * 0.31)


def test_dnpwbt_narrow_wing_branch_expressions():
    """BW/BH < 1.5 applies the arm outside and adds the jet term."""
    result = _wing_body_tail(wing_span=15.0, tail_span=12.0,
                             qoqi=[0.95], deda=[0.31], jet_term=[0.02])
    assert result['branch'] == 'narrow_wing'
    base = _wing_body()
    arm = 23.0 / 4.6667
    increment = 2.0 * (1.08 + 0.13) * 0.95 * 0.06
    assert result['clq'] == pytest.approx(base['clq'] + (increment + 0.02) * arm)
    assert result['cmq'] == pytest.approx(base['cmq'] - (increment + 0.02) * arm**2)


def test_dnpwbt_narrow_acceleration_uses_the_jet_term_with_doubled_arm():
    """CLDWBT = CLADWB - 2*DXOCB*DTJ; the source doubles the arm first."""
    result = _wing_body_tail(wing_span=15.0, tail_span=12.0,
                             qoqi=[0.95], deda=[0.31], jet_term=[0.02])
    base = _wing_body()
    doubled = 2.0 * 23.0 / 4.6667
    assert result['clad'] == pytest.approx(base['clad'] - doubled * 0.02)
    assert result['cmad'] == pytest.approx(
        base['cmad'] + (doubled**2 / 2.0) * 0.02)


def test_dnpwbt_branch_splits_at_a_span_ratio_of_1_5():
    wide = _wing_body_tail(wing_span=18.0, tail_span=12.0)
    narrow = _wing_body_tail(wing_span=17.9, tail_span=12.0)
    assert wide['span_ratio'] == pytest.approx(1.5)
    assert wide['branch'] == 'wide_wing'
    assert narrow['branch'] == 'narrow_wing'


def test_tail_increases_pitch_damping():
    """Adding the tail must make CMQ more negative."""
    base = _wing_body()
    result = _wing_body_tail()
    assert np.all(result['cmq'] < base['cmq'])


def test_dnpwbt_narrow_branch_gates_acceleration():
    """Transonic flow and a tapered wing both suppress the branch."""
    for kwargs, gate in ((dict(transonic=True), 'transonic'),
                         (dict(taper_ratio=0.4), 'tapered_wing')):
        result = _wing_body_tail(wing_span=15.0, tail_span=12.0,
                                 qoqi=[0.95], deda=[0.31], **kwargs)
        assert not result['acceleration_available']
        assert result['acceleration_gate'] == gate
        assert 'clad' not in result


def test_dnpwbt_rejects_bad_input():
    with pytest.raises(ValueError):
        _wing_body_tail(qoqi=[0.9, 0.9], deda=[0.3])
    with pytest.raises(ValueError):
        _wing_body_tail(tail_span=0.0)
    with pytest.raises(ValueError):
        _wing_body_tail(jet_term=[0.01])
