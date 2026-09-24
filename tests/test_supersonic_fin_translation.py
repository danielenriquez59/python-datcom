"""
Regression tests for MASRAT, SUPLAV and SUPLAF, the supersonic
vertical-panel sideslip increments.

Source of truth: datcom-legacy/datcom_2000/masrat.f, suplav.f, suplaf.f.

Checked against a compiled probe (tools/probes/suplav.py) on every
``SBETA``, panel-block, ``BWV`` and ``BWHV`` word the routines set.
"""

import json
import pathlib
import sys

import numpy as np
import pytest

from pydatcom.aerodynamics import supersonic_fin as module
from pydatcom.aerodynamics.supersonic_fin import calculate_suplav, masrat

_ROOT = pathlib.Path(__file__).resolve().parent.parent
sys.path.insert(0, str(_ROOT / 'tools'))
from fortran_data import parse  # noqa: E402

_PROBE = json.loads(
    (_ROOT / 'tests' / 'fixtures' / 'probes' / 'suplav.json').read_text())


def _run(c):
    return calculate_suplav(c['alpha'], c['ventral'], c)


@pytest.mark.parametrize("case", range(len(_PROBE)))
def test_matches_compiled_routine(case):
    c, o = _PROBE[case]['inputs'], _PROBE[case]['outputs']
    r = _run(c)
    off = 31 if c['ventral'] else 0
    sla = o['SLA']
    for key, index in [('rkvwb', 24), ('rkvb', 25), ('rkpvwb', 26),
                       ('dcybv', 27), ('rkvhb', 28), ('zp', 29),
                       ('rlp', 30)]:
        if key in r:
            assert r[key] == pytest.approx(sla[off + index - 1], rel=1e-12,
                                           abs=1e-15), key
    panel, na = o['P'], len(c['alpha'])
    assert r['vt141'] == pytest.approx(panel[0], rel=1e-12)
    assert r['vt161'] == pytest.approx(panel[20], rel=1e-12)
    vt181 = r['vt181'] if r['vt181'] is not None else [0.0] * na
    np.testing.assert_allclose(vt181, panel[40:40 + na], rtol=1e-12,
                               atol=1e-16)
    for block, tag in (('bwv', 'BWV'), ('bwhv', 'BWHV')):
        for index, value in r[block].items():
            assert value == pytest.approx(o[tag][index - 141], rel=1e-12,
                                          abs=1e-16), (block, index)


def test_probe_reaches_every_branch():
    heights = set()
    for p in _PROBE:
        c = p['inputs']
        r1 = c['vtin']['4'] - c['vtin']['3']
        z = -c['position']['zw'] / r1
        flip = (c['position']['vertup'] if c['ventral']
                else not c['position']['vertup'])
        heights.add(-z if flip else z)
        if c['htpl']:
            zt = c['position']['zh'] / r1          # ALIH = 0 cases
            if c['position']['alih'] == 0.0:
                heights.add(-zt if flip else zt)
    assert {0.0, 1.0, -1.0} <= heights
    assert any(abs(h) < 1 and h not in (0.0,) for h in heights)
    assert any(abs(h) > 1 for h in heights)
    for v in (False, True):
        mine = [p['inputs'] for p in _PROBE if p['inputs']['ventral'] == v]
        assert any(c['htpl'] and c['straight'] for c in mine)
        assert any(not c['straight'] for c in mine)
        assert any(not c['htpl'] for c in mine)


def test_ventral_fin_tail_on_repeats_tail_off():
    """SUPLAF's DCVWHB shares VF(141), which is reset to the tail-off DCYBV
    before the wing-body-tail sums read it."""
    for p in _PROBE:
        c = p['inputs']
        if c['ventral'] and c['htpl'] and c['straight']:
            r = _run(c)
            assert r['dcvwhb'] == r['dcybv']
            added = p['outputs']['BWHV'][0] - c['bwhv']['141']
            assert added == pytest.approx(r['dcybv'], rel=1e-12)
            own = dict(c, ventral=False)
            assert _run(own)['dcvwhb'] != _run(own)['dcybv']
            return
    pytest.fail('no ventral tail-on case')


def test_masrat_interpolates_between_the_three_positions():
    lo, mid, hi = (masrat(0.3, 0.5, z) for z in (-1.0, 0.0, 1.0))
    assert masrat(0.3, 0.5, 0.5) == pytest.approx((mid + hi) / 2)
    assert masrat(0.3, 0.5, -0.5) == pytest.approx((lo + mid) / 2)
    assert masrat(0.3, 0.5, 2.0) == hi


def test_tables_match_the_source():
    source = parse('masrat')
    for name, values in source.items():
        if name.startswith('I') and len(values) == 1:
            continue
        np.testing.assert_array_equal(getattr(module, '_' + name), values,
                                      err_msg=name)
    for stem in ('suplav', 'suplaf'):
        s = parse(stem)
        assert s['XAMF'] == module._XAMF and s['YAMF'] == module._YAMF
