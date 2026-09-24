"""
Regression tests for SUPWBT, the supersonic wing-body-tail.

Source of truth: datcom-legacy/datcom_2000/supwbt.f.

Checked against a compiled probe (tools/probes/supwbt.py) on every
``STP``, ``BWH``, ``BWHV``, ``FACT`` and ``BD`` word the routine sets.
"""

import json
import pathlib
import sys

import numpy as np
import pytest

from pydatcom.aerodynamics import supersonic_wbt as module
from pydatcom.aerodynamics import supersonic_wing_body as swb
from pydatcom.aerodynamics import transonic_buildup as tb
from pydatcom.aerodynamics.supersonic_wbt import calculate_supwbt

_ROOT = pathlib.Path(__file__).resolve().parent.parent
sys.path.insert(0, str(_ROOT / 'tools'))
from fortran_data import parse  # noqa: E402

_PROBE = json.loads(
    (_ROOT / 'tests' / 'fixtures' / 'probes' / 'supwbt.json').read_text())

_CLOSE = dict(rtol=1e-10, atol=1e-14)


def _run(c):
    return calculate_supwbt(c, user_downwash=c['user'])


@pytest.mark.parametrize("case", range(len(_PROBE)))
def test_matches_compiled_routine(case):
    c, o = _PROBE[case]['inputs'], _PROBE[case]['outputs']
    r = _run(c)
    if r is None:
        assert not any(o['BWH']) and not any(o['BWHV'])
        return
    na, stp = r['nalpha'], o['STP']
    for key, start in [('cmah', 2), ('cltb', 22), ('cdawb', 42),
                       ('gamma', 111), ('ivbh', 133)]:
        np.testing.assert_allclose(r[key], stp[start - 1:start - 1 + na],
                                   err_msg=key, **_CLOSE)
    if r['ivwh']:
        np.testing.assert_allclose(r['ivwh'], stp[70:70 + na], **_CLOSE)
        np.testing.assert_allclose(r['deltat'], stp[90:90 + na], **_CLOSE)
        np.testing.assert_allclose(r['fact102'], o['FACT'][0:na], **_CLOSE)
        np.testing.assert_allclose(r['fact122'], o['FACT'][20:20 + na],
                                   **_CLOSE)
    for key, index in [('dd', 62), ('trino', 63), ('rkbw', 64), ('kbw', 65),
                       ('kwb', 66), ('clahb', 67), ('clabh', 68),
                       ('yt', 69), ('rcreo2', 70), ('kkbw', 131),
                       ('kkwb', 132), ('dxacwb', 153), ('cd0wbt', 154),
                       ('cd0wbv', 155)]:
        if key in r:
            # INTKBW's ACOS near one amplifies rounding: 1e-9 for its words.
            assert r[key] == pytest.approx(stp[index - 1], rel=1e-9,
                                           abs=1e-14), key
    for block, keys in [('BWH', ('cd', 'cl', 'cm', 'cn', 'ca', 'cla',
                                 'cma')),
                        ('BWHV', ('cdv', 'cl', 'cm', 'cnv', 'cav', 'cla',
                                  'cma'))]:
        for k, key in enumerate(keys):
            np.testing.assert_allclose(r[key], o[block][20 * k:20 * k + na],
                                       err_msg=f'{block} {key}', **_CLOSE)
    np.testing.assert_allclose(
        [r['bd'][k] for k in (2, 3, 58, 63, 64, 84, 761, 762)], o['BD'],
        **_CLOSE)
    # The zero leading-edge sweep snapped to 1e-5 is restored on exit.
    assert o['AHT'][0] == c['tail']['a']['62']


def test_probe_reaches_every_branch():
    results = [_run(p['inputs']) for p in _PROBE]
    live = [r for r in results if r is not None]
    assert any(r is None for r in results)                     # JDETCH = 0
    assert any(r['nalpha'] < len(p['inputs']['alpha'])
               for r, p in zip(results, _PROBE) if r)          # shortened
    assert any('trino' not in r for r in live)                 # triangular
    assert any('rkbw' not in r for r in live)                  # Fig. 10
    assert any(r.get('rkbw') == 0.0 for r in live)             # DX <= -CR
    assert any(r['fact102'] for r in live)                     # canard
    assert any(0.0 in r['fact102'][1:] for r in live)          # ALPAHT = 0
    assert any(r['kkbw'] == 0.33 for r in live)                # stale KK
    assert any(p['inputs']['user'] for p in _PROBE)


def test_canard_reads_dedalp_as_epsilon():
    """FACT(J+101) multiplies by EPSLON(J+20), which the DWASH layout makes
    DEDALP(J): changing the downwash angles leaves that factor's first
    term alone, changing DEDALP moves it."""
    c = json.loads(json.dumps(_PROBE[4]['inputs']))
    base = _run(c)['fact102']
    c['dwash']['dedalp'] = [d + 0.1 for d in c['dwash']['dedalp']]
    moved = _run(c)['fact102']
    assert any(a != b for a, b in zip(base, moved) if a)


def test_tables_match_the_source():
    source = parse('supwbt')
    np.testing.assert_array_equal(tb._TFIG10, source['TFIG10'])
    np.testing.assert_array_equal(tb._DKWB10, source['DKWB10'])
    np.testing.assert_array_equal(tb._DKBW10, source['DKBW10'])
    np.testing.assert_array_equal(swb._T4312A, source['T4312A'])
    np.testing.assert_array_equal(swb._D4312A, source['D4312A'])
    for name in ('T4312B', 'D4312B', 'T4467', 'D4467'):
        np.testing.assert_array_equal(getattr(module, '_' + name),
                                      source[name], err_msg=name)
    # SUPWB moves Figure 4.3.1.2-12B's end abscissae inward; SUPWBT does not.
    assert swb._T4312B[0] == 0.015 and module._T4312B[0] == 0.0
