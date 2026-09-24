"""
Regression tests for SDWASH, the supersonic downwash and dynamic pressure
at the horizontal tail.

Source of truth: datcom-legacy/datcom_2000/sdwash.f.

Checked against a compiled probe (tools/probes/sdwash.py) on every
``/SUPDW/`` and ``/IDWASH/`` word the routine sets, and on ``NALPHA``,
``NF`` and ``JDETCH``.  INFTGM is stubbed in the probe; its ``A`` words are
inputs here.
"""

import json
import pathlib

import numpy as np
import pytest

from pydatcom.aerodynamics.supersonic_flow_field import calculate_sdwash

_ROOT = pathlib.Path(__file__).resolve().parent.parent
_PROBE = json.loads(
    (_ROOT / 'tests' / 'fixtures' / 'probes' / 'sdwash.json').read_text())

_CLOSE = dict(rtol=1e-10, atol=1e-14)


def _run(c):
    return calculate_sdwash(c, user_epsilon=c['user'])


@pytest.mark.parametrize("case", range(len(_PROBE)))
def test_matches_compiled_routine(case):
    c, o = _PROBE[case]['inputs'], _PROBE[case]['outputs']
    r = _run(c)
    dwa, idw = o['DWA'], o['IDW']
    na = len(c['alpha'])
    assert (r['nalpha'], r['nf'], r['jdetch']) == tuple(
        int(v) for v in o['INT'])
    assert r['beta'] == pytest.approx(dwa[1], rel=1e-15)
    if 'x' in r:
        for key, start in [('x', 3), ('y', 5), ('z', 7), ('dhb', 69),
                           ('zeff', 29), ('depx', 129), ('depavg', 109)]:
            np.testing.assert_allclose(
                r[key], dwa[start - 1:start - 1 + len(r[key])],
                err_msg=key, **_CLOSE)
    for key, start in [('alpha', 9), ('clanl', 169), ('m', 189),
                       ('zc', 210)]:
        np.testing.assert_allclose(r[key], dwa[start - 1:start - 1 + na],
                                   err_msg=key, **_CLOSE)
    assert r['zwakec'] == pytest.approx(dwa[208], rel=1e-12)
    assert r['delqo'] == pytest.approx(dwa[229], rel=1e-12)
    d = r['dpresr']
    if d is not None:
        np.testing.assert_allclose(
            [d['dle'], d['deltaz'], d['xsur'], d.get('theta1', 0.0),
             d.get('delte', 0.0), d.get('thete', 0.0)], dwa[230:236],
            **_CLOSE)
    for key, start in [('qqinfy', 0), ('dwangl', 20), ('depda', 40)]:
        np.testing.assert_allclose(r[key], idw[start:start + na],
                                   err_msg=key, **_CLOSE)


def test_probe_reaches_every_branch():
    results = [(p['inputs'], _run(p['inputs'])) for p in _PROBE]
    assert {r.get('icase') for _, r in results} >= {1, 2, 3, 4}
    assert any('x' not in r and not c['user'] for c, r in results)  # VISDW
    assert any(c['user'] and 'x' not in r for c, r in results)
    assert any(r['jdetch'] > 0 and 'x' in r for _, r in results)   # RETURN
    assert any(r['jdetch'] > 0 and 'x' not in r for _, r in results)
    assert any(r['nalpha'] == 0 for _, r in results)                # JDETCH 0
    assert any(r['nf'] == c['nf'] != 0 for c, r in results)
    assert any(any(abs(z / r['zwakec']) <= 1 for z in r['zc'])
               for _, r in results)                                  # wake
    assert any(r['dpresr'] is not None for _, r in results)


def test_simple_downwash_for_supersonic_edges():
    """With a supersonic leading edge the downwash angle is the simple
    1.62*CL/(pi*A) estimate on the wing's own area, in degrees."""
    c, r = next((p['inputs'], _run(p['inputs'])) for p in _PROBE
                if not p['inputs']['user'] and 'x' not in _run(p['inputs']))
    a = c['a']
    for j in range(r['nalpha']):
        expect = (1.62 * c['wing']['cl'][j] / (3.141592654 * a['120']) *
                  57.2957795 * c['sref'] / a['3'])
        assert r['dwangl'][j] == pytest.approx(expect, rel=1e-14)
