"""
Regression tests for FG6115, JETPWE and M30O36, the jet power effects.

Source of truth: datcom-legacy/datcom_2000/fg6115.f, jetpwe.f, m30o36.f.

Checked against a compiled probe of the overlay (tools/probes/jetpwe.py)
on every ``/POWR/`` and ``/IPOWER/`` word it sets.  The probe's cases run
in one program, so the replay carries JETPWE's saved ``COSAIH``.
"""

import json
import pathlib
import sys

import numpy as np
import pytest

from pydatcom.aerodynamics import jet_power as module
from pydatcom.aerodynamics.jet_power import (calculate_jetpwe, fg6115,
                                             m30o36_block)

_ROOT = pathlib.Path(__file__).resolve().parent.parent
sys.path.insert(0, str(_ROOT / 'tools'))
from fortran_data import parse  # noqa: E402

_PROBE = json.loads(
    (_ROOT / 'tests' / 'fixtures' / 'probes' / 'jetpwe.json').read_text())

_WORDS = {'atp': 1, 'deuda': 24, 'epslon': 25, 'atj': 26, 'xep': 47,
          'zjp': 48, 'xjp': 49, 'xhp': 50, 'ain': 51, 'vin': 52,
          'tinotj': 53, 'vjpovi': 54, 'zjporj': 55, 'de': 56, 'zjpobh': 57,
          'ytob2h': 58, 'debode': 59, 'zjpxhp': 60, 'srtpco': 61,
          'zjdexh': 62, 'comp1': 63, 'pteopi': 64, 'rjporj': 65, 'rjp': 66,
          'dxporj': 67, 'dxp': 68, 'xepc': 69, 'xhpc': 70, 'ztp': 71,
          'zbart': 93, 'dcmt1': 94, 'xl': 114, 'dlh': 135}
_CURVES = {'dclt': 2, 'dclnj': 27, 'dclhe': 73, 'dcmnj': 115, 'dcme': 136}
_CLOSE = dict(rtol=1e-11, atol=1e-14)


def _replay():
    cosaih, out = 0.0, []
    for p in _PROBE:
        c = p['inputs']
        r = calculate_jetpwe(c['alpha'], c['mach'], c['jet'], c['win'],
                             c['a'], c['position'], c['htpl'], c['tail'],
                             c['sref'], c['cbarr'],
                             {'cosaih': cosaih, 'comp1': c['comp1']})
        cosaih = r['cosaih']
        out.append(r)
    return out


_RESULTS = _replay()


@pytest.mark.parametrize("case", range(len(_PROBE)))
def test_matches_compiled_overlay(case):
    c, o = _PROBE[case]['inputs'], _PROBE[case]['outputs']
    r = _RESULTS[case]
    pw, power = o['PW'], o['POWER']
    na = len(c['alpha'])
    for key, index in _WORDS.items():
        if key in r:
            assert r[key] == pytest.approx(pw[index - 1], rel=1e-11,
                                           abs=1e-14), key
    for key, start in _CURVES.items():
        np.testing.assert_allclose(r[key], pw[start - 1:start - 1 + na],
                                   err_msg=key, **_CLOSE)
    for key, start in [('dcdpow', 1), ('dclpow', 21), ('dcmpow', 41)]:
        np.testing.assert_allclose(r[key], power[start - 1:start - 1 + na],
                                   err_msg=key, **_CLOSE)
    s = m30o36_block(c['alpha'], r['dcdpow'], r['dclpow'], r['dcmpow'])
    for key, start in [('cn', 61), ('ca', 81), ('cla', 101), ('cma', 121)]:
        np.testing.assert_allclose(s[key], power[start - 1:start - 1 + na],
                                   rtol=1e-10, atol=1e-14, err_msg=key)
    assert r['win1'] == o['WIN1'][0]


def test_probe_reaches_every_branch():
    kases = {r['kase'] for r in _RESULTS}
    assert kases == {None, 1, 2, 3}
    v = [r['vjpovi'] for r in _RESULTS if r.get('kase') in (1, 3)]
    assert any(x <= 2 for x in v) and any(x >= 8 for x in v)
    assert any(2 < x < 4 for x in v) and any(4 < x < 8 for x in v)
    assert any(abs(x - 4) < 1e-2 for x in v)
    wins = [(p['inputs']['win'], r['win1']) for p, r in zip(_PROBE,
                                                          _RESULTS)]
    assert any(w['2'] != 1e-30 and r == w['5'] for w, r in wins)


def test_inboard_jet_overwrites_the_tip_chord():
    """A jet on the inboard panel of a cranked wing leaves WINGIN(1), the
    tip chord, set to the break chord WINGIN(5)."""
    hits = 0
    for p, r in zip(_PROBE, _RESULTS):
        c = p['inputs']
        inboard = (c['jet']['nengsj'] != 1.0 and c['win']['2'] != 1e-30 and
                   c['jet']['jelloc'] <= c['win']['4'] - c['win']['2'])
        if inboard:
            hits += 1
            assert p['outputs']['WIN1'][0] == c['win']['5'] != c['win']['1']
    assert hits


def test_no_tail_reads_the_previous_incidence():
    """With no tail, DLH uses COSAIH from the previous call."""
    k = next(i for i, p in enumerate(_PROBE) if not p['inputs']['htpl'])
    assert _RESULTS[k]['cosaih'] == _RESULTS[k - 1]['cosaih']
    c = _PROBE[k]['inputs']
    fresh = calculate_jetpwe(c['alpha'], c['mach'], c['jet'], c['win'],
                             c['a'], c['position'], False, c['tail'],
                             c['sref'], c['cbarr'],
                             {'cosaih': 0.0, 'comp1': c['comp1']})
    assert fresh['dlh'] != _RESULTS[k]['dlh']


def test_fg6115_blends_parts():
    a = fg6115(10.0, 4.0, 2.0)
    b = fg6115(10.0, 4.0, 4.0)
    c = fg6115(10.0, 4.0, 8.0)
    assert fg6115(10.0, 4.0, 3.0) == pytest.approx((a + b) / 2)
    assert fg6115(10.0, 4.0, 6.0) == pytest.approx((b + c) / 2)
    assert fg6115(10.0, 4.0, 4.005) == b


def test_tables_match_the_source():
    fg, jp = parse('fg6115'), parse('jetpwe')
    for source in (fg, jp):
        for name, values in source.items():
            if name.startswith('I') and len(values) == 1:
                continue
            np.testing.assert_array_equal(getattr(module, '_' + name),
                                          values, err_msg=name)
