"""
Regression tests for SUPWB, SUPHB and their overlay M20O24.

Source of truth: datcom-legacy/datcom_2000/supwb.f, suphb.f, m20o24.f.

Checked against a compiled probe of the whole overlay
(tools/probes/supwb.py), including SUPCM0's write over the first moment.
"""

import json
import pathlib
import sys

import numpy as np
import pytest

from pydatcom.aerodynamics import supersonic_wing_body as module
from pydatcom.aerodynamics import transonic_buildup as tb
from pydatcom.aerodynamics.supersonic_wing_body import (
    calculate_supwb, m20o24_slopes,
)
from pydatcom.aerodynamics.transonic_buildup import calculate_supcm0
from pydatcom.utils.constants import UNUSED

_ROOT = pathlib.Path(__file__).resolve().parent.parent
sys.path.insert(0, str(_ROOT / 'tools'))
from fortran_data import parse  # noqa: E402

_PROBE = json.loads(
    (_ROOT / 'tests' / 'fixtures' / 'probes' / 'supwb.json').read_text())


def _blocks(c, key):
    s = c[key]
    return ({int(k): v for k, v in s['in'].items()},
            {int(k): v for k, v in s['a'].items()}, s['aero'])


def _surface_run(c, key, moment_block=None):
    win, a, aero = _blocks(c, key)
    p = c['position']
    tail = key == 'tail'
    st = c['stale']
    return calculate_supwb(
        c['alpha'], c['mach'],
        {'spans': win[3], 'span': win[4], 'cr': win[6], 'type': win[15]}, a,
        {'x': p['xh'] if tail else p['xw'],
         'incidence': p['alih'] if tail else p['aliw'],
         'zw': p['zw'], 'zcg': p['zcg']},
        c['body'], aero, c['sref'], c['cbarr'], tail=tail,
        moment_block=moment_block,
        stale={'kkwb': st['hkkwb' if tail else 'kkwb'],
               'kkbw': st['hkkbw' if tail else 'kkbw']})


def _cm0(c, key, x, stale):
    win, a, _ = _blocks(c, key)
    surface = {'sspn': win[4], 'sspne': win[3], 'chrdr': win[6],
               'tovc': win[16], 'ler': win[62], 'twista': win[11],
               'ycm': win[93], 'cld': win[94], 'a38': a[38], 'a80': a[80],
               'a118': a[118], 'a120': a[120], 'a122': a[122]}
    return calculate_supcm0(surface, x, c['position']['zw'],
                            c['body']['x'][-1], c['body']['sbd120'],
                            c['mach'], c['tr'], stale)['cm0']


def _replay(c):
    wing = _surface_run(c, 'wing')
    out = {'wing': wing}
    if c['htpl']:
        block = {k: wing[k] for k in ('kwb', 'kkwb', 'kbw', 'kkbw', 'xacbw',
                                      'dd', 'ivbw', 'gamma', 'xacw')}
        block['incidence'] = c['position']['aliw']
        out['tail'] = _surface_run(c, 'tail', block)
    out['wing_cm0'] = _cm0(c, 'wing', c['position']['xw'], wing['cm'][0])
    if c['htpl']:
        out['tail_cm0'] = _cm0(c, 'tail', c['position']['xh'],
                               out['tail']['cm'][0])
    return out


def _check_block(r, cm0, arr, swb, alpha):
    na = len(alpha)
    close = dict(rtol=1e-11, atol=1e-16)
    np.testing.assert_allclose(arr[0:na], r['cd'], **close)
    np.testing.assert_allclose(arr[20:20 + na], r['cl'], **close)
    np.testing.assert_allclose(arr[40:40 + na], [cm0] + r['cm'][1:], **close)
    s = m20o24_slopes(alpha, r)
    np.testing.assert_allclose(arr[60:60 + na], s['cn'], **close)
    np.testing.assert_allclose(arr[80:80 + na], s['ca'], **close)
    np.testing.assert_allclose(arr[100:100 + na], s['cla'], **close)
    np.testing.assert_allclose(arr[120:120 + na], s['cma'], **close)
    for index, key in [(2, 'kkwb'), (3, 'xacn'), (4, 'cd0'), (5, 'dd'),
                       (6, 'beta'), (7, 'clabw'), (8, 'xacbw'), (9, 'fa'),
                       (10, 'cli'), (11, 'kbw'), (32, 'rkbw'),
                       (33, 'clawb'), (34, 'fn'), (35, 'kwb'), (36, 'xac'),
                       (37, 'kkbw'), (38, 'rlap'), (39, 'xaca'),
                       (60, 'trino'), (61, 'xcpln')]:
        if r.get(key) is not None:
            # INTKBW's ACOS near one amplifies rounding: 1e-9 for its words.
            assert swb[index - 1] == pytest.approx(r[key], rel=1e-9,
                                                   abs=1e-16), key
    np.testing.assert_allclose(swb[11:11 + na], r['ivbw'], **close)
    np.testing.assert_allclose(swb[39:39 + na], r['gamma'], **close)


@pytest.mark.parametrize("case", range(len(_PROBE)))
def test_matches_compiled_overlay(case):
    c, o = _PROBE[case]['inputs'], _PROBE[case]['outputs']
    r = _replay(c)
    _check_block(r['wing'], r['wing_cm0'], o['BW'], o['SWB'], c['alpha'])
    na = len(c['alpha'])
    np.testing.assert_allclose(o['BWV'][:60], o['BW'][:60])
    slg134, stg134, delxw, delxh, _ = o['X']
    assert slg134 == pytest.approx(r['wing']['xacw'], rel=1e-12)
    assert delxw == pytest.approx(r['wing']['delxw'], rel=1e-12)
    np.testing.assert_allclose(o['AB'], r['wing']['alpha_body'], rtol=1e-14)
    assert o['CM0'][0] == pytest.approx(r['wing_cm0'], rel=1e-11)
    if c['htpl']:
        _check_block(r['tail'], r['tail_cm0'], o['BH'], o['SHB'], c['alpha'])
        assert stg134 == pytest.approx(r['tail']['xacw'], rel=1e-12)
        assert delxh == pytest.approx(r['tail']['delxw'], rel=1e-12)
    else:
        assert o['BH'][:na] == [0.0] * na


def test_probe_reaches_every_branch():
    results = [_replay(p['inputs'])['wing'] for p in _PROBE]
    assert any('rkbw' in r and r['rkbw'] > 0 for r in results)
    assert any('rkbw' not in r for r in results)
    assert any(r.get('rkbw') == 0.0 for r in results)       # DX <= -CR
    assert any('trino' not in r for r in results)           # untapered
    assert any(r['cli'] is None for r in results)           # not straight
    assert any(r['rlap'] == 0.0 for r in results)
    assert any(p['inputs']['body']['bnose'] == 1.0 for p in _PROBE)
    assert any(p['inputs']['position']['aliw'] == 0.0 for p in _PROBE)


def test_tail_moment_reads_the_wing_blocks():
    """SUPHB's moment loop takes the wing's factors: replacing them with the
    tail's own changes the compiled result the translation reproduces."""
    c = _PROBE[0]['inputs']
    r = _replay(c)
    own = {k: r['tail'][k] for k in ('kwb', 'kkwb', 'kbw', 'kkbw', 'xacbw',
                                     'dd', 'ivbw', 'gamma', 'xacw')}
    own['incidence'] = c['position']['alih']
    alt = _surface_run(c, 'tail', own)
    assert not np.allclose(alt['cm'][1:], r['tail']['cm'][1:])
    np.testing.assert_allclose(_PROBE[0]['outputs']['BH'][41:46],
                               r['tail']['cm'][1:], rtol=1e-11)


def test_supcm0_overwrites_the_first_moment():
    """SUPCM0 stores its regression CM0 in BW(41), SUPWB's moment at the
    first angle, whenever the regression is in range."""
    hits = 0
    for p in _PROBE:
        r = _replay(p['inputs'])
        if r['wing_cm0'] != r['wing']['cm'][0]:
            hits += 1
            assert p['outputs']['BW'][40] == pytest.approx(r['wing_cm0'])
    assert hits


def test_tables_match_the_source():
    source = parse('supwb')
    assert source == parse('suphb')
    for name, values in source.items():
        if name.startswith('I') and len(values) == 1:
            continue
        mod = module if hasattr(module, '_' + name) else tb
        np.testing.assert_array_equal(getattr(mod, '_' + name), values,
                                      err_msg=name)
