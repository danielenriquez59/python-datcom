"""
Regression tests for the NACA section coordinate routines: COORD1, COORD4,
COORD5, COORD6, CORD4M, CORD5M and XYCORD.

Source of truth: datcom-legacy/datcom_2000/coord*.f, cord4m.f, cord5m.f,
xycord.f.

Checked against a compiled probe (tools/probes/naca_sections.py) on every
surface, mean-line and thickness ordinate and the scalar words set.
"""

import json
import pathlib

import pytest

import pydatcom.geometry.naca_sections as module

_ROOT = pathlib.Path(__file__).resolve().parent.parent
_PROBE = json.loads(
    (_ROOT / 'tests' / 'fixtures' / 'probes' / 'naca_sections.json')
    .read_text())

_ROUTINES = {'COORD1': module.coord1, 'COORD4': module.coord4,
             'COORD5': module.coord5, 'COORD6': module.coord6,
             'CORD4M': module.cord4m, 'CORD5M': module.cord5m}


def _run(c):
    if c['routine'] == 'XYCORD':
        return module.xycord(c['x'], c['yu'], c['yl'])
    return _ROUTINES[c['routine']](c['digits'], c['x'])


@pytest.mark.parametrize("case", range(len(_PROBE)))
def test_matches_compiled_routine(case):
    c, o = _PROBE[case]['inputs'], _PROBE[case]['outputs']
    r = _run(c)
    for key in ('xu', 'xl', 'yun', 'yln', 'thn', 'cam'):
        assert r[key] == pytest.approx(o[key.upper()], abs=1e-14), key
    s = o['S']
    for key, index in (('rho', 0), ('t', 1), ('zm', 2), ('zp', 3),
                       ('alphai', 4), ('alphao', 5)):
        if key in r:
            assert r[key] == pytest.approx(s[index], rel=1e-12,
                                           abs=1e-15), key


def test_every_routine_is_probed():
    assert {p['inputs']['routine'] for p in _PROBE} == set(_ROUTINES) | {
        'XYCORD'}


def test_coord1_zeroes_the_lower_leading_edge_twice():
    """COORD1's closing block writes XL(1) twice where the other routines
    write XU(1): the upper leading-edge station keeps what was there."""
    c = next(p['inputs'] for p in _PROBE if p['inputs']['routine'] == 'COORD1')
    prev = {k: [7.0] * len(c['x']) for k in ('xu', 'xl', 'yun', 'yln',
                                             'thn', 'cam')}
    r = module.coord1(c['digits'], c['x'], prev)
    assert r['xu'][0] == 7.0 and r['yun'][0] == 7.0
    assert r['xl'][0] == 0.0 and r['yln'][0] == 0.0


def test_six_series_trailing_edge_is_straight():
    """For a subscripted 6-series section, stations past 0.8 chord lie on
    the straight line from the first such station to the trailing edge."""
    c = next(p['inputs'] for p in _PROBE if p['inputs']['routine'] == 'COORD6'
             and p['inputs']['digits']['ii'] > 0 and
             p['inputs']['digits']['jj'] > 0)
    r = module.coord6(c['digits'], c['x'])
    pts = [(x, y) for x, y in zip(r['xu'], r['yun']) if x >= 0.8]
    (x0, y0), (x1, y1) = pts[0], pts[1]
    slope = -y0 / (1. - x0)
    assert y1 == pytest.approx(y0 + slope * (x1 - x0), abs=1e-14)


_CORDSP = json.loads(
    (_ROOT / 'tests' / 'fixtures' / 'probes' / 'cordsp.json').read_text())


@pytest.mark.parametrize("case", range(len(_CORDSP)))
def test_cordsp_matches_compiled_routine(case):
    c, o = _CORDSP[case]['inputs'], _CORDSP[case]['outputs']
    r = module.cordsp(c['digits'], c['x'], c['surface_in'])
    for key, tags in (('xu', ('XU',)), ('xl', ('XL',)),
                      ('yu', ('YU', 'YUU')), ('yl', ('YL', 'YLL'))):
        for tag in tags:
            assert r[key] == pytest.approx(o[tag], abs=1e-15), tag
    s = r['surface_in']
    assert [r['rho'], r['toc'], s[16], s[18], s[62], s[63], s[70],
            s[71]] == pytest.approx(o['S'], abs=1e-15)


def test_cordsp_zeroes_the_outboard_radius_of_a_straight_wing():
    """The planform test is a mixed-mode comparison that never matches, so
    WGIN(63) is set even for a straight tapered wing (planform 1)."""
    c, o = _CORDSP[0]['inputs'], _CORDSP[0]['outputs']
    assert c['surface_in']['15'] == 1 and c['surface_in']['63'] == 1e-30
    assert o['S'][5] == 0.0


_AIRFOL = json.loads(
    (_ROOT / 'tests' / 'fixtures' / 'probes' / 'airfol.json').read_text())
_U = 1e-30


@pytest.mark.parametrize("case", range(len(_AIRFOL)))
def test_decode_and_airfol_match_compiled_routines(case):
    card, o = _AIRFOL[case]['inputs']['card'], _AIRFOL[case]['outputs']
    d = module.decode(card)
    assert [d['na'], d['l']] + list(d['digits'].values()) == [
        int(v) for v in o['D']]
    r = module.airfol(card, {15: 1.0, 16: _U, 18: _U, 62: _U, 63: _U,
                             70: _U, 71: _U})
    for key in ('x', 'xu', 'xl', 'yun', 'yln', 'thn', 'cam'):
        assert r[key] == pytest.approx(o[key.upper()], abs=1e-14), key


def test_decode_reaches_every_family():
    families = {module.decode(p['inputs']['card'])['na'] for p in _AIRFOL}
    assert families == {1, 2, 3, 4, 5, 6, 7}
