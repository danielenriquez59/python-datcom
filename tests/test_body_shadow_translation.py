"""
Regression tests for PTINT2 and BDAREA, the body area in the horizontal
tail's Mach zone.

Source of truth: datcom-legacy/datcom_2000/ptint2.f, bdarea.f.

Checked against a compiled probe (tools/probes/bdarea.py) on the three
``HTIN`` words set and the abort flag.
"""

import json
import pathlib

import pytest

import pydatcom.aerodynamics.body_shadow as module
from pydatcom.aerodynamics.body_shadow import calculate_bdarea

_ROOT = pathlib.Path(__file__).resolve().parent.parent
_PROBE = json.loads(
    (_ROOT / 'tests' / 'fixtures' / 'probes' / 'bdarea.json').read_text())


def _run(c):
    return calculate_bdarea(c['xb'], c['rb'], c['mach'], c['syna'],
                            c['htin'], c['aht'])


@pytest.mark.parametrize("case", range(len(_PROBE)))
def test_matches_compiled_routine(case):
    c, o = _PROBE[case]['inputs'], _PROBE[case]['outputs']['R']
    r = _run(c)
    assert r['abort'] == bool(o[0])
    if r['abort']:
        assert o[1:] == [-2.0, -3.0, -1.0]          # words left alone
    else:
        assert [r['sb'], r['s'], r['xbar']] == pytest.approx(
            o[1:], rel=1e-12)


def test_probe_reaches_every_area_shape(monkeypatch):
    seen = set()
    calls = []
    original = module.ptint2

    def spy(*args):
        r = original(*args)
        calls.append(r)
        return r
    monkeypatch.setattr(module, 'ptint2', spy)
    for p in _PROBE:
        calls.clear()
        if _run(p['inputs'])['abort']:
            seen.add('abort')
            continue
        le, te = calls
        for side, a, b in (('u', 'indxui', 'inxuie'),
                           ('l', 'indxli', 'inxlie')):
            kind = 1 if te[b] == 0 else (2 if le[b] == 0 else 3)
            seen.add((side, kind, le[a] == te[a] if kind == 1 else None))
    for side in 'ul':
        assert {(side, 1, True), (side, 1, False), (side, 2, None),
                (side, 3, None)} <= seen
    assert 'abort' in seen


def test_shadow_is_at_most_the_extended_area():
    for p in _PROBE:
        r = _run(p['inputs'])
        if not r['abort']:
            assert 0.0 < r['sb'] <= r['s'] + 1e-12
