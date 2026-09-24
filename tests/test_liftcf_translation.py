"""
Regression tests for the ANGLES and LIFTCF translations.

Source of truth: datcom-legacy/datcom_2000/angles.f, liftcf.f.

Both are checked against compiled probes of the legacy routines
(tools/probes/angles.py and tools/probes/liftcf.py), and every LIFTCF table
against a re-parse of the FORTRAN DATA statements.
"""

import json
import math
import pathlib
import sys

import numpy as np
import pytest

from pydatcom.aerodynamics import liftcf as module
from pydatcom.aerodynamics.liftcf import calculate_liftcf
from pydatcom.aerodynamics.wtlift import calculate_wtlift
from pydatcom.utils.legacy_numeric import angles, zerang

_ROOT = pathlib.Path(__file__).resolve().parent.parent
sys.path.insert(0, str(_ROOT / 'tools'))
from fortran_data import parse  # noqa: E402

_FIXTURES = _ROOT / 'tests' / 'fixtures' / 'probes'
_ANGLES = json.loads((_FIXTURES / 'angles.json').read_text())
_LIFTCF = json.loads((_FIXTURES / 'liftcf.json').read_text())


# --------------------------------------------------------------------------
# ANGLES
# --------------------------------------------------------------------------

@pytest.mark.parametrize("case", range(len(_ANGLES)))
def test_angles_matches_compiled_routine(case):
    probe = _ANGLES[case]
    np.testing.assert_allclose(angles(probe['entry'], probe['arg']),
                               probe['result'], rtol=1e-12, atol=1e-14)


def test_angles_probe_reaches_every_entry_and_the_early_return():
    assert {p['entry'] for p in _ANGLES} == {1, 2, 3, 4, 5, 6}
    # At least one case returns with a stale 9. still in the record.
    assert any(9.0 in p['result'] for p in _ANGLES)


def test_angles_early_return_keeps_the_previous_values():
    """An angle within 2**-16 rad of the stored one changes nothing else."""
    previous = [10.0, math.radians(10.0), 9.0, 9.0, 9.0, math.radians(10.0)]
    nudged = angles(1, [10.0 + 1e-4] + previous[1:])
    assert nudged[1:] == previous[1:]
    moved = angles(1, [10.1] + previous[1:])
    assert moved[2] == pytest.approx(math.sin(math.radians(10.1)))


def test_angles_wraps_and_zeroes():
    assert angles(1, [270.0] + zerang()[1:])[0] == pytest.approx(-90.0)
    assert angles(1, [1e-4] + [5.0] * 5) == zerang()
    quarter = angles(1, [90.0] + zerang()[1:])
    assert quarter[3] == 0.0 and quarter[4] == 2.0**20


def test_angles_does_not_modify_its_argument():
    arg = [30.0, 0.0, 0.0, 1.0, 0.0, 0.0]
    angles(1, arg)
    assert arg == [30.0, 0.0, 0.0, 1.0, 0.0, 0.0]


# --------------------------------------------------------------------------
# LIFTCF against execution
# --------------------------------------------------------------------------

@pytest.mark.parametrize("case", range(len(_LIFTCF)))
def test_liftcf_matches_compiled_routine(case):
    """CL, CN, and for the straight method B(45) and the angle state."""
    probe = _LIFTCF[case]
    result = calculate_liftcf(**probe['inputs'])
    out = probe['outputs']
    np.testing.assert_allclose(result['cl'], out['CL'], rtol=1e-9, atol=1e-12)
    np.testing.assert_allclose(result['cn'], out['CN'], rtol=1e-9, atol=1e-12)
    if 'b45' in result:
        assert result['b45'] == pytest.approx(out['B45'][0], rel=1e-9)
        np.testing.assert_allclose(result['angle_state'], out['STATE'],
                                   rtol=1e-12, atol=1e-14)
    # The double delta flips negative angles in place and restores them.
    assert out['ALPHA'] == probe['inputs']['alpha_deg']


def _run(probe):
    return calculate_liftcf(**probe['inputs'])


def test_liftcf_probe_reaches_every_branch():
    """Guard the fixture: each planform and branch stays covered."""
    by_type = {}
    for p in _LIFTCF:
        by_type.setdefault(p['inputs']['planform_type'], []).append(p)
    assert set(by_type) == {1.0, 2.0, 3.0, 4.0}
    straight = by_type[1.0]
    results = [_run(p) for p in straight]
    assert {r['low_aspect_ratio'] for r in results} == {True, False}
    assert any(p['inputs']['geometry']['aspect_ratio'] <= 1.0
               for p in straight)
    # Every straight case crosses the stall angle.
    assert all(min(p['inputs']['alpha_deg']) < p['inputs']['lift']['alpha_clmax']
               < max(p['inputs']['alpha_deg']) for p in straight)
    # The stale stall record is actually reused.
    assert any(r['angle_state'][2] == 0.3 for r in results)
    for kind in (2.0, 3.0):
        flags = {_run(p)['dashed_region'] for p in by_type[kind]}
        assert flags == {True, False}
    breaks = [_run(p)['alpha_break'] for p in by_type[3.0]]
    assert any(1.0 / p['inputs']['geometry']['tan_le'] > 1.0
               for p in by_type[3.0])
    assert max(breaks) <= 7.0


@pytest.mark.parametrize("names", [
    ['AJ', 'TRAT', 'DC', 'A58', 'CLJ58', 'TIR', 'D', 'C90I', 'C90',
     'X13356', 'X23356', 'Y13356', 'X33356', 'Y33356', 'X13357', 'Y13357'],
])
def test_every_table_value_matches_the_source(names):
    source = parse('liftcf')
    for name in names:
        np.testing.assert_array_equal(getattr(module, f'_{name}'),
                                      source[name], err_msg=name)


# --------------------------------------------------------------------------
# Physical behaviour of the straight method, through WTLIFT
# --------------------------------------------------------------------------

def _wing(aspect_ratio=6.0, sref=150.0, alpha=None):
    geometry = dict(area=150.0, aspect_ratio=aspect_ratio, taper_ratio=0.5,
                    sweep_le_deg=10.0, cos_le=math.cos(math.radians(10.0)),
                    tan_le=math.tan(math.radians(10.0)), tan_c2=0.1,
                    arclss_classified=1.3, arclss_ratio=3.0)
    section = dict(deltay=2.0, xovc=0.3, cla=0.1, clmax=1.4)
    flight = dict(mach=0.3, beta=math.sqrt(1 - 0.09), alpha_zero_lift=-2.0)
    lift = calculate_wtlift(1.0, geometry, section, flight, sref)
    geometry.update(arclss_factor=0.3, a159=lift['a159'], a160=lift['a160'])
    if alpha is None:
        alpha = list(np.arange(-6.0, 30.0, 1.0))
    return lift, calculate_liftcf(1.0, alpha, geometry, section, lift,
                                  flight, sref)


def test_lift_is_zero_at_the_zero_lift_angle():
    _, curve = _wing(alpha=[-2.0])
    assert curve['cl'][0] == 0.0


def test_initial_slope_matches_wtlift():
    """Near zero lift the curve's slope is WTLIFT's CLa, plus the small
    nonlinear CNaa*|sin|*sin term that vanishes there."""
    lift, curve = _wing(alpha=[-2.1, -1.9])
    slope = (curve['cl'][1] - curve['cl'][0]) / 0.2
    assert slope == pytest.approx(lift['cla'], rel=1e-3)


def test_lift_peaks_near_the_predicted_stall():
    lift, curve = _wing()
    alpha = np.arange(-6.0, 30.0, 1.0)
    peak = alpha[np.argmax(curve['cl'])]
    assert abs(peak - lift['alpha_clmax']) <= 1.5
    assert max(curve['cl']) == pytest.approx(lift['clmax'], rel=0.05)
    # And it falls beyond it.
    assert curve['cl'][-1] < max(curve['cl'])


def test_lift_force_is_independent_of_reference_area():
    """CL*SREF is a force, so it cannot depend on the reference chosen."""
    _, a = _wing(sref=150.0)
    _, b = _wing(sref=300.0)
    np.testing.assert_allclose(a['cl'] * 150.0, b['cl'] * 300.0, rtol=1e-12)


def test_unscaled_cnaa90_below_unit_aspect_ratio_is_preserved():
    """At A <= 1 the source leaves Figure 4.1.3.3-55B's CNAA90 off the SREF
    basis.  Below the stall the force is still reference-independent; past
    it, it is not, which is how the omission shows.
    """
    alpha = [0.0, 4.0, 50.0, 70.0]
    lift, a = _wing(aspect_ratio=0.8, sref=150.0, alpha=alpha)
    _, b = _wing(aspect_ratio=0.8, sref=300.0, alpha=alpha)
    assert lift["alpha_clmax"] < 50.0
    np.testing.assert_allclose(a['cl'][:2] * 150.0, b['cl'][:2] * 300.0,
                               rtol=1e-12)
    assert abs(a['cl'][3] * 150.0 - b['cl'][3] * 300.0) > 1e-3 * abs(
        a['cl'][3] * 150.0)
