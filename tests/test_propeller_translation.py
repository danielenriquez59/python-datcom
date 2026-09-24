"""
Regression tests for PRPWEF and M13O15, the propeller power effects.

Source of truth: datcom-legacy/datcom_2000/prpwef.f, m13o15.f.

Checked against a compiled probe (tools/probes/prpwef.py) on every
``/POWR/`` word the routine sets; the tables are re-parsed from source.
"""

import json
import math
import pathlib
import sys

import numpy as np
import pytest

from pydatcom.aerodynamics import propeller as module
from pydatcom.aerodynamics.propeller import (
    POWR_LAYOUT, calculate_m13o15, calculate_prpwef,
)

_ROOT = pathlib.Path(__file__).resolve().parent.parent
sys.path.insert(0, str(_ROOT / 'tools'))
from fortran_data import parse  # noqa: E402

_PROBE = json.loads(
    (_ROOT / 'tests' / 'fixtures' / 'probes' / 'prpwef.json').read_text())


def _inputs(c):
    w, t, v = dict(c['wing']), dict(c['tail']), dict(c['vertical'])
    w['cf'] = w['d10'] + w['d11']
    t['cf'] = t['dht10'] + t['dht11']
    v['cf'] = v['dvt10'] + v['dvt11']
    return {'alpha': c['alpha'], 'sref': c['sref'], 'cbarr': c['cbarr'],
            'htpl': c['htpl'], 'power': c['power'], 'wing': w,
            'a': {int(k): x for k, x in c['a'].items()},
            'position': c['position'], 'tail': t, 'vertical': v,
            'body': c['body'], 'dwash': c['dwash']}


def _run(c, routine=calculate_m13o15):
    return routine(_inputs(c), c['stale'])


def _close(ours, theirs, name):
    if math.isnan(theirs):
        assert math.isnan(ours), name
    else:
        assert ours == pytest.approx(theirs, rel=1e-10, abs=1e-14), name


@pytest.mark.parametrize("case", range(len(_PROBE)))
def test_matches_compiled_routine(case):
    c, o = _PROBE[case]['inputs'], _PROBE[case]['outputs']
    r = _run(c)
    pw, na = o['PW'], len(c['alpha'])
    for name, start, length in POWR_LAYOUT:
        if length == 20:
            for j in range(na):
                _close(r[name][j], pw[start - 1 + j], f'{name}[{j}]')
        elif name == 'argcs':
            for k in range(6):
                _close(r['argcs'][k], pw[start - 1 + k], f'argcs[{k}]')
        elif name in r['powr']:
            _close(r['powr'][name], pw[start - 1], name)
    for tag, key, start in [('POWER', 'cdpow', 0), ('POWER', 'dclpon', 20),
                            ('POWER', 'dcm', 40), ('POWER', 'cn', 60),
                            ('POWER', 'ca', 80), ('POWER', 'cla', 100),
                            ('POWER', 'cma', 120)]:
        for j in range(na):
            _close(r[key][j], o[tag][start + j], f'{key}[{j}]')
    kn, cosaiw, _, nalpha = o['OUT']
    assert r['kn'] == pytest.approx(kn, rel=1e-12)
    assert r['cosaiw'] == pytest.approx(cosaiw, rel=1e-14)
    assert r['nalpha'] == nalpha


def test_probe_reaches_every_branch():
    results = [_run(p['inputs'], calculate_prpwef) for p in _PROBE]
    powr = [r['powr'] for r in results]
    assert any(s['deuda'] == -1.0 for s in powr)
    assert any(s['deuda'] != -1.0 for s in powr)
    assert any('scapi' in s for s in powr)           # past the break
    assert any('tri' in s for s in powr)             # inner panel
    assert any(p['inputs']['power']['nengsp'] == 2.0 for p in _PROBE)
    assert any(p['inputs']['power']['nengsp'] > 2.0 for p in _PROBE)
    assert any(s['cm0i'] != s['cm0te0'] for s in powr)
    assert any(r['kn'] != p['inputs']['power']['kn']
               for r, p in zip(results, _PROBE))


def test_centreline_sweep_is_radians_used_as_degrees():
    """A single propeller takes the sweep from A(69), the ANGLES test word
    of the quarter-chord record (radians), and applies COS(DEG*SWEEPA)."""
    c = _PROBE[0]['inputs']
    r = _run(c, calculate_prpwef)['powr']
    assert r['sweepa'] == c['a']['69']
    assert r['cossw'] == pytest.approx(math.cos(0.01745329 * c['a']['69']))
    assert r['cossw'] > 0.9999


def test_upwash_lookup_is_dimensional():
    """Figure 4.4.1-61 is read at XBARP in feet, not XBARP/CRP."""
    c = _PROBE[0]['inputs']
    r = _run(c, calculate_prpwef)['powr']
    assert r['xbarp'] > 2.0          # beyond the chart's 0.25..2.0 grid
    assert r['deuda'] != -1.0


def test_tables_match_the_source():
    source = parse('prpwef')
    for name, values in source.items():
        if name in ('Y415A', 'Y415B', 'Y415C'):
            continue
        if name == 'PISQRD':
            assert module._PISQRD == values[0]
            continue
        if name.startswith('I') and len(values) == 1:
            continue
        np.testing.assert_array_equal(getattr(module, '_' + name), values,
                                      err_msg=name)
    np.testing.assert_array_equal(
        module._Y41412, source['Y415A'] + source['Y415B'] + source['Y415C'])


def test_nan_lookup_returns_the_last_node():
    """GLOOK's comparisons all fail for NaN, so the lookup lands on the last
    node; the probe's two-engine cases (ASTARI = 0/0) depend on it."""
    from pydatcom.utils.legacy_tables import tlin1x
    assert tlin1x([1.0, 2.0, 4.0], [3.0, 5.0, 9.0], math.nan, 2, 2) == 9.0
    assert tlin1x([4.0, 2.0, 1.0], [3.0, 5.0, 9.0], math.nan) == 9.0
    with pytest.raises(ValueError):
        tlin1x([1.0, 2.0], [3.0, 5.0], math.inf)
    two_engine = [p for p in _PROBE if p['inputs']['power']['nengsp'] == 2.0]
    assert two_engine and all(math.isnan(p['outputs']['PW'][210])
                              for p in two_engine)
