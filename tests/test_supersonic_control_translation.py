"""
Regression tests for DFLCON, the supersonic control derivatives.

Source of truth: datcom-legacy/datcom_2000/dflcon.f.

Checked against a compiled probe (tools/probes/dflcon.py) on its four
``/POWR/`` outputs over tapered and untapered, inboard and tip controls.
"""

import json
import pathlib

import pytest

from pydatcom.aerodynamics.supersonic_control import calculate_dflcon

_ROOT = pathlib.Path(__file__).resolve().parent.parent
_PROBE = json.loads(
    (_ROOT / 'tests' / 'fixtures' / 'probes' / 'dflcon.json').read_text())


@pytest.mark.parametrize("case", range(len(_PROBE)))
def test_matches_compiled_routine(case):
    c, o = _PROBE[case]['inputs'], _PROBE[case]['outputs']['R']
    r = calculate_dflcon(c['a'], c['d'], c['fl'], c['beta'], c['inbord'],
                         c['taperd'])
    if r is None:
        assert o == [-9.0] * 4                  # outputs left alone
        return
    assert [r['cld'], r['clld'], r['cmd'], r['chd']] == pytest.approx(
        o, rel=1e-9, abs=1e-10)        # cancellation in the sums


def test_probe_covers_each_form():
    kinds = {(p['inputs']['taperd'], p['inputs']['inbord']) for p in _PROBE}
    assert kinds == {(True, True), (True, False), (False, True),
                     (False, False)}
    assert any(abs(p['inputs']['d']) > 1 for p in _PROBE)
    assert any(p['inputs']['fl'] == 0.0 for p in _PROBE)
