"""
Regression tests for CLMCH0 and the M15O17 lift pass (CALCA0, WTLIFT,
LIFTCF in sequence).

Source of truth: datcom-legacy/datcom_2000/clmch0.f, m15o17.f, calca0.f.

Checked against a compiled probe (tools/probes/clmch0.py) that runs CLMCH0
and then M15O17 at a flight Mach number, as the main program does; the
second pass inherits LIFTCF's angle state from the first.  This is also the
first execution-backed check of CALCA0's twist and camber branches.
"""

import json
import pathlib

import numpy as np
import pytest

from pydatcom.aerodynamics.surface_lift import (
    calculate_clmch0, calculate_surface_lift, section_alpha_zero,
)
from pydatcom.utils.constants import UNUSED

_ROOT = pathlib.Path(__file__).resolve().parent.parent
_PROBE = json.loads(
    (_ROOT / 'tests' / 'fixtures' / 'probes' / 'clmch0.json').read_text())


def _passes(inputs):
    mach0 = calculate_clmch0(inputs['planform_type'], inputs['alpha_deg'],
                             inputs['geometry'], inputs['section'],
                             inputs['sref'])
    state = mach0['surface']['liftcf'].get('angle_state')
    flight = calculate_surface_lift(inputs['planform_type'],
                                    inputs['alpha_deg'], inputs['geometry'],
                                    inputs['section'], inputs['flight'],
                                    inputs['sref'], state)
    return mach0, flight


def _close(actual, expected, what):
    np.testing.assert_allclose(actual, expected, rtol=1e-9, atol=1e-12,
                               err_msg=what)


@pytest.mark.parametrize("case", range(len(_PROBE)))
def test_matches_compiled_routines(case):
    probe = _PROBE[case]
    out = probe['outputs']
    mach0, flight = _passes(probe['inputs'])
    if 'cl' in mach0:
        _close(mach0['cl'], out['B3'], 'B(3)')
        _close([mach0['cla'], mach0['alpha_zero_lift'], mach0['alpha_clmax']],
               out['M0'][:3], 'B(48), A(126), A(127)')
    _close(mach0['surface']['calca0']['alpha_zero_lift'], out['M0'][1],
           'A(126)')
    _close(flight['liftcf']['cl'], out['CL'], 'CL')
    _close(flight['liftcf']['cn'], out['CN'], 'CN')
    wtlift = flight['wtlift']
    if wtlift['computed']:
        _close([wtlift['cla'], wtlift['alpha_clmax'], wtlift['clmax'],
                flight['alpha_zero_lift']], out['M1'][:4],
               'WING(101), B(43), B(44), B(49)')
    if 'b45' in flight['liftcf']:
        _close(mach0['surface']['liftcf']['angle_state'], out['S0'],
               'state after the Mach-zero pass')
        _close(flight['liftcf']['b45'], out['M1'][4], 'B(45)')
        _close(flight['liftcf']['angle_state'], out['S1'],
               'state after the flight pass')


def test_probe_reaches_every_branch():
    passes = [_passes(p['inputs'])[1] for p in _PROBE]
    calca0 = [f['calca0'] for f in passes]
    assert any(c['twist_term'] is not None for c in calca0)
    assert any(c['camber_factor'] is not None for c in calca0)
    assert any(c['twist_term'] is not None and c['camber_factor'] is not None
               for c in calca0)
    assert any(p['inputs']['section']['swafp'] != UNUSED for p in _PROBE)
    assert {p['inputs']['planform_type'] for p in _PROBE} == {1.0, 2.0,
                                                               3.0, 4.0}
    assert {f['wtlift'].get('low_aspect_ratio') for f in passes
            if f['wtlift']['computed']} == {True, False}


def test_section_alpha_zero_prefers_the_supplied_value():
    assert section_alpha_zero(-2.2, 1.1, 0.15, 0.1) == -2.2
    assert section_alpha_zero(UNUSED, 1.1, 0.15, 0.1) == pytest.approx(-0.4)


def test_mach_zero_pass_is_independent_of_the_flight_section_data():
    """CLMCH0 uses CLALPA(1) and CLMAXL, never the flight Mach's values."""
    inputs = json.loads(json.dumps(_PROBE[0]['inputs']))
    a = calculate_clmch0(inputs['planform_type'], inputs['alpha_deg'],
                         inputs['geometry'], inputs['section'], inputs['sref'])
    inputs['section'].update(cla=0.2, clmax=3.0)
    b = calculate_clmch0(inputs['planform_type'], inputs['alpha_deg'],
                         inputs['geometry'], inputs['section'], inputs['sref'])
    np.testing.assert_array_equal(a['cl'], b['cl'])
    assert a['alpha_clmax'] == b['alpha_clmax']
