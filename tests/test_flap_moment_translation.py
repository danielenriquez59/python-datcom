"""
Regression tests for FLAPCM and M37O45, the flap pitching moment.

Source of truth: datcom-legacy/datcom_2000/flapcm.f, m37o45.f.

Checked against a compiled probe (tools/probes/flapcm.py) on the moment
increments and every ``/SUPWH/`` and ``TCD`` word.  The probe's cases run
in one program, so the replay carries FLAPCM's saved span stations.
"""

import json
import pathlib

import pytest

from pydatcom.aerodynamics.flap_moment import (ET_DATA, calculate_flapcm,
                                               m37o45)

_ROOT = pathlib.Path(__file__).resolve().parent.parent
_PROBE = json.loads(
    (_ROOT / 'tests' / 'fixtures' / 'probes' / 'flapcm.json').read_text())


def _replay():
    state = {'et': ET_DATA, 'kinbd': 0, 'koutbd': 0}
    out = []
    for p in _PROBE:
        r = calculate_flapcm(dict(p['inputs'], state=state))
        state = r['state']
        out.append(r)
    return out


_RESULTS = _replay()


@pytest.mark.parametrize("case", range(len(_PROBE)))
def test_matches_compiled_routine(case):
    c, o, r = _PROBE[case]['inputs'], _PROBE[case]['outputs'], \
        _RESULTS[case]
    n = int(c['f'][15] + .5)
    assert r['delcm'][:n] == pytest.approx(o['D'][:n], rel=1e-11)
    assert r['fcm'] == pytest.approx(o['FCM'], rel=1e-10, abs=1e-14)
    assert r['tcd'] == pytest.approx(o['TCD'], rel=1e-10, abs=1e-14)
    assert r['flp'][0] == o['ETA'][0] and r['flp'][4] == o['ETA'][1]


def test_probe_reaches_every_branch():
    types = {p['inputs']['f'][16] for p in _PROBE}
    assert {1.0, 2.0, 5.0, 6.0} <= types
    assert any(p['inputs']['f'][18] != 1e-30 for p in _PROBE)   # SDCL
    assert any(p['inputs']['f'][28] != 1e-30 for p in _PROBE)   # SCMD
    assert any(p['inputs']['htpl'] for p in _PROBE)
    assert any(r['flp'][0] != p['inputs']['flp'][0]
               for r, p in zip(_RESULTS, _PROBE))              # nudged


def test_station_edit_carries_to_the_next_call():
    """Case 9 repeats case 0's inputs, but case 0 overwrote a standard
    span station with its inboard flap edge."""
    assert _PROBE[9]['inputs'] == _PROBE[0]['inputs']
    assert _RESULTS[9]['delcm'][0] != pytest.approx(_RESULTS[0]['delcm'][0],
                                                    rel=1e-3)
    assert _RESULTS[0]['state']['et'] != ET_DATA


def test_m37o45_routes():
    calls = []
    r = m37o45(False, 5.0, 0.5, 0.3, lambda: calls.append(1) or 'flapcm',
               None)
    assert r == 'flapcm' and calls == [1]
    tail = {'tante': -0.1, 'tanle': 0.7, 'bsto2': 5.0, 'crh': 4.0}
    r = m37o45(True, 5.0, 0.5, 0.3, lambda: calls.append(2), tail)
    assert calls == [1] and len(r['gdh']) == 4
