"""
Regression tests for SSHING, the supersonic hinge moments.

Source of truth: datcom-legacy/datcom_2000/sshing.f.

Checked against a compiled probe (tools/probes/sshing.py, SSSYM stubbed)
on all 59 ``/POWR/`` words and the ``WING`` words set.  The probe's cases
run in one program, so the replay carries SSHING's saved locals.
"""

import json
import math
import pathlib

import pytest

from pydatcom.aerodynamics.ptcp import calculate_ptcp
from pydatcom.aerodynamics.supersonic_hinge import calculate_sshing

_ROOT = pathlib.Path(__file__).resolve().parent.parent
_PROBE = json.loads(
    (_ROOT / 'tests' / 'fixtures' / 'probes' / 'sshing.json').read_text())


def _replay():
    loc = {k: 0.0 for k in ('lamhl', 'k', 'mu', 'xt', 'xr', 'ci', 'ct',
                            'tovca', 'kaseno')}
    out = []
    for p in _PROBE:
        r = calculate_sshing(dict(p['inputs'], stale=loc))
        loc = r['locals']
        out.append(r)
    return out


_RESULTS = _replay()


@pytest.mark.parametrize("case", range(len(_PROBE)))
def test_matches_compiled_routine(case):
    o, r = _PROBE[case]['outputs'], _RESULTS[case]
    assert r['spr'] == pytest.approx(o['SPR'], rel=1e-9, abs=1e-10)
    words = [-7.0] * 5
    if r.get('clalds'):
        words[:len(r['clalds'])] = r['clalds']
    if r.get('chabc') is not None:
        words[3:] = [r['chabc'], r['chdelr']]
    assert words == pytest.approx(o['W'], rel=1e-9, abs=1e-10)


def test_probe_reaches_every_branch():
    assert any(r['returned'] for r in _RESULTS)
    assert any(r.get('clalds') for r in _RESULTS)
    assert {r['kaseno'] for r in _RESULTS} >= {1, 2, 3, 5, 6}
    assert any(p['inputs']['htpl'] for p in _PROBE)
    assert any(p['inputs']['a']['62'] == 0.0 for p in _PROBE)


def test_hinge_moment_drops_the_tip():
    """CHAT sums SPAMT+SPAMR after the loop, when SPAMT is the root's."""
    i = next(k for k, r in enumerate(_RESULTS)
             if r.get('chabc') is not None and not r['returned'])
    c, spr = _PROBE[i]['inputs'], _RESULTS[i]['spr']
    f = c['f']
    aloci, aloco, cfi, cfo = f['14'], f['15'], f['12'], f['13']
    root = sum(spr[44:52])
    tip = sum(spr[36:44])
    base = -2. / (57.2957795 * spr[0] * math.sqrt(1. - spr[31]**2))
    denom = (aloco - aloci) * (cfi**2 + cfi * cfo + cfo**2)
    assert spr[52] == pytest.approx(base * (1. - 2. * root / denom))
    assert tip != root


def test_unswept_root_repeats_the_tip_moments():
    i = next(k for k, p in enumerate(_PROBE) if p['inputs']['a']['62'] == 0.)
    spr = _RESULTS[i]['spr']
    assert spr[44:52] == spr[36:44]


def test_ptcp_last_station_follows_the_source_loop():
    """The ten-point loops build EN by adding 0.1, reaching 0.99999... where
    an exact 1.0 would put ARCCOS at its -1 defect (0 instead of pi)."""
    p = calculate_ptcp(0.92661, 3, 1.0, 0.3, math.sqrt(1.05**2 - 1), -0.1)
    assert p['pp'][9] == pytest.approx(1.0, abs=1e-6)
    assert p['pressure_ratio'] == pytest.approx(0.3400525976526264,
                                                rel=1e-9)
