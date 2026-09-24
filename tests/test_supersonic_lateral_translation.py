"""
Regression tests for SUPLAT and SUPLAH, the supersonic sideslip
derivatives of the wing, tail and body.

Source of truth: datcom-legacy/datcom_2000/suplat.f, suplah.f.

Checked against a compiled probe (tools/probes/suplat.py) on every
``SBETA`` word and every surface and body word the routines set.
"""

import json
import math
import pathlib
import sys

import numpy as np
import pytest

from pydatcom.aerodynamics import supersonic_lateral as module
from pydatcom.aerodynamics.supersonic_lateral import calculate_suplat

_ROOT = pathlib.Path(__file__).resolve().parent.parent
sys.path.insert(0, str(_ROOT / 'tools'))
from fortran_data import parse  # noqa: E402

_PROBE = json.loads(
    (_ROOT / 'tests' / 'fixtures' / 'probes' / 'suplat.json').read_text())


def _same(a, b):
    if math.isnan(a) or math.isnan(b):
        return math.isnan(a) and math.isnan(b)
    if math.isinf(a) or math.isinf(b):
        return a == b
    return a == pytest.approx(b, rel=1e-10, abs=1e-14)


@pytest.mark.parametrize("case", range(len(_PROBE)))
def test_matches_compiled_routine(case):
    c, o = _PROBE[case]['inputs'], _PROBE[case]['outputs']
    r = calculate_suplat(c, tail=c['tail'])
    for k, v in r['sla'].items():
        assert _same(v, o['SLA'][k - 1]), ('sla', k)
    for key, off in (('cybw', 0), ('cnbw', 20), ('clbw', 40)):
        for j, v in enumerate(r[key]):
            assert _same(v, o['SURF'][off + j]), (key, j)
    for block, tag in (('bwi', 'BODY'), ('bwh', 'BWH'), ('bwv', 'BWV'),
                       ('bwhv', 'BWHV')):
        for k, v in r.get(block, {}).items():
            assert _same(v, o[tag][k - 141]), (block, k)
    w = r['win']
    assert [w[12], w[13], w[14]] == o['W'][:3]
    if not c['tail']:
        assert r['lf'] == int(o['LF'][0])


def test_probe_reaches_every_branch():
    rs = [(p['inputs'], calculate_suplat(p['inputs'],
                                         tail=p['inputs']['tail']))
          for p in _PROBE]
    assert any(5 in r['sla'] for c, r in rs if not c['tail'])     # delta
    assert any(3 in r['sla'] and 5 not in r['sla'] for c, r in rs)
    assert any(not c['bo'] for c, _ in rs)
    assert any(r.get('lf') == 1 for _, r in rs)
    assert any(21 in r['sla'] for _, r in rs)                     # K_H(B)
    assert any(c['tail'] for c, _ in rs)


def test_delta_branch_boundaries_differ():
    """At beta/tan(LE) = 1 SUPLAH takes the delta branch (<= 1), where
    1/Q = 0 makes its side force infinite; SUPLAT (< 0.998) does not."""
    tail = next(p for p in _PROBE if p['inputs']['tail'] and
                p['inputs']['a']['118'] == 0.0)
    wing = _PROBE[-1]
    assert wing['inputs']['a']['62'] == tail['inputs']['a']['62']
    assert math.isinf(tail['outputs']['SURF'][0])
    assert wing['outputs']['SURF'][0] == 0.0


def test_suplah_swept_forward_power_is_nan():
    p = next(p for p in _PROBE if p['inputs']['tail'] and
             p['inputs']['a']['62'] < 0)
    assert all(math.isnan(v) for v in p['outputs']['SURF'][40:44])


def test_tables_match_the_source():
    source = parse('suplat')
    assert parse('suplah') == {k: v for k, v in source.items()
                               if k not in ('XA25OO', 'XB25OO', 'YA25OO',
                                            'YB25OO')}
    for name, values in source.items():
        if name.startswith('Y2225'):
            continue
        np.testing.assert_array_equal(getattr(module, '_' + name), values,
                                      err_msg=name)
    np.testing.assert_array_equal(
        module._Y12225, sum((source['Y2225' + p] for p in 'ABCDE'), []))


def test_m23o27_words():
    from pydatcom.aerodynamics.supersonic_lateral import m23o27_words
    from pydatcom.utils.constants import UNUSED
    body = {k: 0.1 * k for k in range(141, 190)}
    vt = {k: 0.01 * k for k in range(141, 190)}
    vf = {k: 0.001 * k for k in range(141, 190)}
    r = m23o27_words(3, body, vt, vf)
    assert r['surface'] == {142: -UNUSED, 143: -UNUSED, 162: -UNUSED,
                            163: -UNUSED}
    assert r['bv'][141] == pytest.approx(0.111 * 141)
    assert r['bv'][183] == pytest.approx(0.111 * 183)
    assert r['bv'][142] == -UNUSED and 184 not in r['bv']
