"""
Regression tests for M41O51, M53O65 and M56O70.

Source of truth: datcom-legacy/datcom_2000/m41o51.f, m53o65.f, m56o70.f.

M56O70 is checked against a compiled probe (tools/probes/m56o70.py) that
runs VTAREA for both panels and BDAREA through the overlay.
"""

import json
import pathlib

import pytest

from pydatcom.aerodynamics.supersonic_overlays import m41o51, m53o65, m56o70
from pydatcom.utils.constants import UNUSED

_ROOT = pathlib.Path(__file__).resolve().parent.parent
_PROBE = json.loads(
    (_ROOT / 'tests' / 'fixtures' / 'probes' / 'm56o70.json').read_text())
_I = 2


@pytest.mark.parametrize("case", range(len(_PROBE)))
def test_m56o70_matches_compiled_overlay(case):
    c, o = _PROBE[case]['inputs'], _PROBE[case]['outputs']
    r = m56o70(c)
    ours = [r[b][k + _I] for b in ('vtin', 'vfin', 'htin')
            for k in (94, 114, 134)] + [r['syna'][6], r['syna'][7]]
    assert ours == pytest.approx(o['R'], rel=1e-11, abs=1e-14)
    assert r['ig'] == int(o['IG'][0])


def test_absent_panel_still_runs_when_a_later_word_is_unset():
    """PL .AND. A .OR. B .OR. C: with no vertical tail but VTIN(114+I)
    unset, VTAREA still fills the words."""
    c = _PROBE[1]['inputs']
    assert not c['vtpl'] and c['vtin'][str(114 + _I)] == UNUSED
    assert m56o70(c)['vtin'][114 + _I] != UNUSED
    c2 = _PROBE[2]['inputs']
    assert m56o70(c2)['vtin'][114 + _I] == 3.3


def test_tail_position_is_restored():
    for p in _PROBE:
        c = p['inputs']
        r = m56o70(c)
        assert (r['syna'][6], r['syna'][7]) == (c['syna']['6'],
                                                c['syna']['7'])


def test_m41o51_and_m53o65():
    calls = []
    r, words = m41o51(5.0, lambda: calls.append(1))
    assert r is None and not calls
    assert words == {250 + j: -UNUSED for j in range(2, 11)}
    m41o51(1.0, lambda: calls.append(1))
    assert calls == [1]
    dl, dr = list(range(10)), [0.5 * v for v in range(10)]
    r, words = m53o65(3.0, dl, dr, lambda: calls.append(2))
    assert r is None and calls == [1]
    assert words[201] == 0.0 and words[210] == 4.5
