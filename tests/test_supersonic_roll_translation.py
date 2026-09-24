"""
Regression tests for SPRYAW, the supersonic control roll and yaw.

Source of truth: datcom-legacy/datcom_2000/spryaw.f.

Checked against a compiled probe (tools/probes/spryaw.py) on all 59
``/POWR/`` words and the ``HT``, ``BODY`` and ``WING`` curves it sets.
"""

import json
import pathlib

import pytest

from pydatcom.aerodynamics.supersonic_roll import calculate_spryaw

_ROOT = pathlib.Path(__file__).resolve().parent.parent
_PROBE = json.loads(
    (_ROOT / 'tests' / 'fixtures' / 'probes' / 'spryaw.json').read_text())


def _blocks(r):
    ht, body, wing = [0.0] * 20, [0.0] * 30, [0.0] * 30
    for i, v in enumerate(r.get('clrlal', []) + r.get('clrlsp', [])):
        ht[i] = v
    for i, v in enumerate(r.get('cnywsp', [])):
        ht[10 + i] = v
    for i, v in r.get('cnywal', {}).items():
        body[i] = v
    for i, v in r.get('clrlht', {}).items():
        wing[i] = v
    return ht, body, wing


@pytest.mark.parametrize("case", range(len(_PROBE)))
def test_matches_compiled_routine(case):
    c, o = _PROBE[case]['inputs'], _PROBE[case]['outputs']
    r = calculate_spryaw(c)
    ht, body, wing = _blocks(r)
    assert r['spr'] == pytest.approx(o['SPR'], rel=1e-9, abs=1e-14)
    assert ht == pytest.approx(o['HT'][:20], rel=1e-9, abs=1e-14)
    assert body == pytest.approx(o['BODY'][:30], rel=1e-9, abs=1e-14)
    assert wing == pytest.approx(o['WING'][:30], rel=1e-9, abs=1e-14)


def test_probe_reaches_every_branch():
    rs = [calculate_spryaw(p['inputs']) for p in _PROBE]
    assert any('clrlal' in r for r in rs)
    assert any('clrlsp' in r for r in rs)
    assert any('clrlht' in r for r in rs)
    assert sum(r['returned'] for r in rs) == 2


def test_inboard_flag_is_reversed_from_sshing():
    """Case 0's flap ends at 78% of the exposed span, inboard of the tip,
    yet SPRYAW hands DFLCON INBORD=.FALSE.: the compiled derivatives are
    the tip-control ones."""
    from pydatcom.aerodynamics.supersonic_control import calculate_dflcon
    c = _PROBE[0]['inputs']
    assert c['f']['15'] / c['win']['3'] < 0.99
    s = _PROBE[0]['outputs']['SPR']
    trtofl = c['f']['13'] / c['f']['12']
    tip = calculate_dflcon(s[32], s[33], trtofl, s[0], False, True)
    inboard = calculate_dflcon(s[32], s[33], trtofl, s[0], True, True)
    assert s[11] == pytest.approx(tip['cld'], rel=1e-9)
    assert s[11] != pytest.approx(inboard['cld'], rel=1e-6)
