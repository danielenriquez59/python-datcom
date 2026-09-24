"""
Regression tests for the transonic buildup: TRANCM/TRHTCM with
WBTRAN/HBTRAN and WBCM1, TRACM0, TRANCD, and WBTRA with TRAWBT.

Source of truth: datcom-legacy/datcom_2000/trancm.f, trhtcm.f, wbtran.f,
hbtran.f, wbcm1.f, tracm0.f, trancd.f, trawbt.f, wbtra.f.

Checked against a compiled probe (tools/probes/transonic_buildup.py) that
runs the routines in overlay order on shared COMMON, so each pass sees
what the one before it left.
"""

import json
import pathlib
import sys

import numpy as np
import pytest

from pydatcom.aerodynamics import transonic_buildup as module
from pydatcom.aerodynamics.transonic_buildup import (
    calculate_trancd, calculate_tracm0, calculate_trancm, calculate_trawbt,
    calculate_wbcm1, m24o30_body_words,
)
from pydatcom.aerodynamics.wing_body import _X21C, _X38B, _Y21C, _Y38B
from pydatcom.utils.constants import PI, UNUSED

_ROOT = pathlib.Path(__file__).resolve().parent.parent
sys.path.insert(0, str(_ROOT / 'tools'))
from fortran_data import parse  # noqa: E402

_PROBE = json.loads((_ROOT / 'tests' / 'fixtures' / 'probes' /
                     'transonic_buildup.json').read_text())


def _blocks(case, key):
    s = case[key]
    return ({int(k): v for k, v in s['in'].items()},
            {int(k): v for k, v in s['a'].items()})


def _surface(win, a):
    return {'type': win[15], 'sspn': win[4], 'sspne': win[3],
            'chrdr': win[6], 'tovc': win[16], 'ler': win[62],
            'twista': win[11], 'ycm': win[93], 'cld': win[94],
            'a38': a[38], 'a56': a[56], 'a80': a[80], 'a118': a[118],
            'a120': a[120], 'a122': a[122]}


def _replay(c):
    """The probe's call sequence: M25O31, TRANCD, then WBTRA."""
    body = c['body']
    out = {}
    for key, x, cla, stale in [('wing', c['xw'], c['cla_wing'], 'wb'),
                               ('tail', c['xh'], c['cla_tail'], 'hb')]:
        if not c['wgpl' if key == 'wing' else 'htpl']:
            continue
        win, a = _blocks(c, key)
        body_args = None
        if c['bo']:
            body_args = {
                'alpha_deg': c['alpha'],
                'surface': {'sspn': win[4], 'sspne': win[3],
                            'chrdr': win[6]},
                'body_length': body['length'],
                'alpha0_body': body['alpha0'], 'cla_body': body['cla'],
                'cma_body': body['cma'], 'sref': c['sref'],
                'bd87': body['bd87'],
                'stale_wb13': c['stale'][stale + '13'],
                'stale_wb14': c['stale'][stale + '14']}
            if key == 'tail':
                # TRHTCM reads the wing's SWB(8); zero if WBTRAN never ran.
                body_args['xacbw_14'] = (out['wing']['wbtran']['xacbw']
                                         if 'wing' in out else 0.0)
        out[key] = calculate_trancm(c['mach'], c['mfb'], win[16], a, cla,
                                    c['cbarr'], c['xcg'], x, c['wgpl'],
                                    body_args)
        if c['bo']:
            out[key + '_cm0'] = calculate_tracm0(
                _surface(win, a), x, c['zw'], body['length'],
                body['sbd120'], c['mach'], c['tr'],
                c['stale'][('bw' if key == 'wing' else 'bh') + '41'])
    flight = {'mach': c['mach'], 'reynolds_per_length': c['reynolds'],
              'tr': c['tr']}
    rbody = {'length': body['length'], 'x_max_area': body['x_max_area'],
             'max_diameter': body['max_diameter']}
    for key, x, cd0, gate in [
            ('wing', c['xw'], c['cd0_wb'], c['bo'] and c['wgpl']),
            ('tail', c['xh'], c['cd0_hb'], c['htpl'])]:
        win, a = _blocks(c, key)
        if gate:
            out[key + '_cd'] = calculate_trancd(c['alpha'], _surface(win, a),
                                                cd0, rbody, x, flight)
    win, a = _blocks(c, 'wing')
    hin, aht = _blocks(c, 'tail')
    wing_body_cla = (out['wing']['wbtran']['cla'] if 'wing' in out and
                     c['bo'] else c['bw101'])
    tail_wbtran = out['tail']['wbtran'] if c['bo'] and 'tail' in out else {}
    kq, _, kd = c['kdwash']
    out['trawbt'] = calculate_trawbt(
        {'sspn': win[4], 'span_break': win[12], 'dihedral_in': win[13],
         'dihedral_out': win[14]}, a, c['b48'],
        {'type': hin[15], 'sspn': hin[4]}, aht,
        {'xcg': c['xcg'], 'xw': c['xw'], 'xh': c['xh'], 'aliw': c['aliw']},
        c['cla_wing'], c['cd0_wing'], wing_body_cla,
        tail_wbtran.get('clawb', 0.0), tail_wbtran.get('clabw', 0.0),
        c['cd0_tail'], out['wing']['xac'] if 'wing' in out else 0.0,
        c['sref'], c['cbarr'],
        c['dwash']['deda'] if kd else None,
        c['dwash']['q_ratio'] if kq else None)
    out['wing_body_cla'] = wing_body_cla
    return out


def _check_trancm(r, tra, w, swb):
    close = dict(rel=1e-11, abs=1e-14)
    assert tra[83 - 71:89 - 71] == pytest.approx(r['xmv'], **close)
    assert tra[89 - 71:95 - 71] == pytest.approx(r['xacv'], **close)
    assert tra[95 - 71] == pytest.approx(r['xacw'], **close)
    assert tra[105 - 71] == pytest.approx(r['xac'], **close)
    if 'delxac' in r:
        assert tra[96 - 71] == pytest.approx(r['delxac'], **close)
        assert tra[97 - 71:105 - 71] == pytest.approx(r['xacp'], **close)
    else:
        assert tra[96 - 71:105 - 71] == [0.0] * 9
    assert w[0] == pytest.approx(r['cma'], **close)
    assert w[6] == pytest.approx(r['dxcg'], **close)
    if 'wbtran' not in r:
        return
    t = r['wbtran']
    assert tra[0] == pytest.approx(t['clawb'], **close)
    assert tra[1] == pytest.approx(t['clabw'], **close)
    assert tra[106 - 71] == pytest.approx(r['xacbw'], **close)
    assert tra[107 - 71] == pytest.approx(r['xacwb'], **close)
    assert w[1] == pytest.approx(t['cla'], **close)
    assert w[2] == pytest.approx(r['cma_wing_body'], **close)
    assert w[3:6] == pytest.approx([r['wbcm1'][k] for k in
                                    ('wb13', 'wb14', 'wb15')], **close)
    for index, key in [(5, 'dd'), (8, 'xacbw'), (11, 'kbw'), (32, 'rkbw'),
                       (35, 'kwb'), (39, 'xaca'), (60, 'trino')]:
        if key in t:
            assert swb[index - 1] == pytest.approx(t[key], **close), key


@pytest.mark.parametrize("case", range(len(_PROBE)))
def test_matches_compiled_routines(case):
    c, o = _PROBE[case]['inputs'], _PROBE[case]['outputs']
    r = _replay(c)
    if 'wing' in r:
        _check_trancm(r['wing'], o['WTRA'], o['W'], o['WSWB'])
        if c['bo']:
            np.testing.assert_allclose(o['AB'], r['wing']['wbtran'][
                'alpha_body'], rtol=1e-14)
    if 'tail' in r:
        _check_trancm(r['tail'], o['HTRA'], o['H'], o['HSWB'])
    for key, tra, w in [('wing', o['WTRA'], o['W']),
                        ('tail', o['HTRA'], o['H'])]:
        if key + '_cm0' in r:
            cm0 = r[key + '_cm0']['cm0']
            assert w[7] == pytest.approx(cm0, rel=1e-11)
            assert tra[74 - 71] == pytest.approx(cm0, rel=1e-11)
    for key, cd0, tag in [('wing', c['cd0_wb'], 'BWCD'),
                          ('tail', c['cd0_hb'], 'BHCD')]:
        cd = r.get(key + '_cd')
        if cd is None:
            assert o[tag][0] == pytest.approx(cd0, rel=1e-14)
        else:
            np.testing.assert_allclose(o[tag], cd, rtol=1e-11)
    q, deda, clawbt, cmawbt, bwh101, bwh121, zwc, zc, dj, dqoq = o['WBT']
    t = r['trawbt']
    if t is None:
        assert (q, deda) == (c['dwash']['q_ratio'], c['dwash']['deda'])
        assert (clawbt, cmawbt, bwh121) == (0.0, 0.0, 0.0)
        return
    close = dict(rel=1e-11, abs=1e-15)
    assert q == pytest.approx(t['q_ratio'], **close)
    assert deda == pytest.approx(t['deda'], **close)
    assert clawbt == pytest.approx(t['cla'], **close)
    assert bwh101 == pytest.approx(t['cla'], **close)
    assert cmawbt == pytest.approx(t['cma'], **close)
    assert bwh121 == pytest.approx(t['cma'], **close)
    assert (zwc, zc, dj, dqoq) == pytest.approx(
        (t['zwc'], t['zc'], t['dj'], t['dqoq']), **close)
    assert o['HTRA'][108 - 71] == pytest.approx(t['cd0_tail'], **close)


def test_probe_reaches_every_branch():
    results = [_replay(p['inputs']) for p in _PROBE]
    wings = [r['wing'] for r in results if 'wing' in r]
    trans = [w['wbtran'] for w in wings if 'wbtran' in w]
    assert any('delxac' in w for w in wings)
    assert any('delxac' not in w for w in wings)
    assert any('rkbw' in t for t in trans) and any('rkbw' not in t
                                                   for t in trans)
    assert any('trino' not in t for t in trans)   # untapered or Mach 1
    assert any(r['wing_cd'] is None for r in results if 'wing_cd' in r)
    assert any(r['wing_cd'] is not None for r in results if 'wing_cd' in r)
    assert any(not r['wing_cm0']['computed'] for r in results
               if 'wing_cm0' in r)
    assert any(r['wing_cm0']['computed'] for r in results if 'wing_cm0' in r)
    assert any(r['trawbt'] is None for r in results)
    assert any(r['trawbt'] is not None for r in results)
    assert any(p['inputs']['wgpl'] and not p['inputs']['bo'] for p in _PROBE)


def test_tables_match_the_source():
    cm, wt, tw = parse('trancm'), parse('wbtran'), parse('trawbt')
    assert cm == parse('trhtcm') and wt == parse('hbtran')
    np.testing.assert_array_equal(module._XM, cm['XM'])
    np.testing.assert_array_equal(module._T422AF, cm['T422AF'])
    for name, parts in [('_SUBAF1', 'ABC'), ('_SUBAF2', 'DEF'),
                        ('_SUPAF1', 'GHI'), ('_SUPAF2', 'JKL')]:
        np.testing.assert_array_equal(
            getattr(module, name), sum((cm['DUMY' + p] for p in parts), []))
    np.testing.assert_array_equal(module._T425AD, cm['T425AD'])
    np.testing.assert_array_equal(
        module._D425AD, sum((cm[f'DUMY{i}'] for i in range(1, 13)), []))
    np.testing.assert_array_equal(module._T428, cm['T428'])
    np.testing.assert_array_equal(module._D428, cm['D428'])
    np.testing.assert_array_equal(module._D4311A, wt['DUMYA'] + wt['DUMYB'])
    for name in ('TFIG10', 'DKWB10', 'DKBW10', 'T4311A', 'T4311B', 'D4311B',
                 'T4337A', 'D4337A', 'T4337B', 'D4337B'):
        np.testing.assert_array_equal(getattr(module, '_' + name), wt[name])
    for name, values in tw.items():
        np.testing.assert_array_equal(getattr(module, '_' + name), values)
    w1 = parse('wbcm1')
    for name, ours in [('X38B', _X38B), ('Y38B', _Y38B), ('X21C', _X21C),
                       ('Y21C', _Y21C)]:
        np.testing.assert_array_equal(ours, w1[name])


def test_tracm0_reynolds_number_is_the_mach_number_times_the_mac():
    """TRACM0 passes ``FLC(I+2)*A(122)``, Mach number times the MAC, so the
    regression's Reynolds term always sits on its 8e5 clamp."""
    c = _PROBE[0]['inputs']
    win, a = _blocks(c, 'wing')
    r = calculate_tracm0(_surface(win, a), c['xw'], c['zw'],
                         c['body']['length'], c['body']['sbd120'],
                         c['mach'], c['tr'])
    assert r['computed']
    assert r['reynolds'] == pytest.approx(c['mach'] * a[122])
    assert r['reynolds'] < 8.0e5


def test_wbcm1_keeps_stale_words_when_the_ellipse_fails():
    """A zero ellipse half-height (the two anchors equal) has no real root;
    the source prints and returns without touching WB(13)/WB(14)."""
    a = {7: 4.0, 10: 7.0, 27: 0.0, 38: 0.5, 44: 0.0}
    # TEMP0 = 0.5 gives WB(15) = 0.25, and A(44) = 0 gives TEMP4 = 0.25.
    r = calculate_wbcm1(a, 15.0, 4.0, 5.5, stale_wb13=0.7, stale_wb14=0.6)
    assert r['wb15'] == pytest.approx(0.25) and r['temp4'] == 0.25
    assert r['ellipse_failed']
    assert (r['wb13'], r['wb14']) == (0.7, 0.6)


def test_carryover_centre_always_moves_aft_with_mach():
    """``ABS(XACBW4-XACBW6)`` makes the 0.6-to-1.4 interpolation slope
    nonnegative whichever anchor is further aft."""
    for p in _PROBE:
        r = _replay(p['inputs']).get('wing')
        if r and 'wbtran' in r:
            slope = (r['xacbw'] - r['wbcm1']['wb13']) / (p['inputs']['mach']
                                                         - 0.6)
            assert slope >= 0.0


def test_tail_reads_the_wings_mach_14_carryover_centre():
    """HBTRAN writes the tail's SWB(8) into the second /SUPWBB/ block;
    TRHTCM reads the first, the wing's."""
    for p in _PROBE:
        c, o = p['inputs'], p['outputs']
        if not (c['bo'] and c['wgpl'] and c['htpl']):
            continue
        r = _replay(c)['tail']
        own = r['wbtran']['xacbw']
        wing = o['WSWB'][7]
        assert o['HSWB'][7] == pytest.approx(own, rel=1e-12)
        x6 = r['wbcm1']['wb13']
        assert o['HTRA'][106 - 71] == pytest.approx(
            x6 + abs(wing - x6) / 0.8 * (c['mach'] - 0.6), rel=1e-12)
        assert abs(own - wing) > 1e-3


def test_trancd_marks_angles_without_data():
    c = _PROBE[0]['inputs']
    win, a = _blocks(c, 'wing')
    cd = calculate_trancd([0.0, 4.0, 13.0], _surface(win, a), 0.02,
                          {'length': 44.0, 'x_max_area': 20.0,
                           'max_diameter': 4.2}, 16.0,
                          {'mach': 0.95, 'reynolds_per_length': 2e6,
                           'tr': 0.4})
    # 13 degrees is past the transonic tables' 11-degree limit.
    assert cd[2] == -UNUSED
    assert np.all(np.abs(cd[:2] - 0.02) < 0.05) and cd[1] > cd[0]


def test_m24o30_body_words():
    w = m24o30_body_words(0.004, 0.02, 5.5, 30.0)
    assert w == {'body_141': -0.004, 'body_161': pytest.approx(-0.02 * 5.5
                                                               / 30.0)}


# --- Second level: WBCLB, CLBCLC and SETUP2 (tools/probes/second_level.py)

from pydatcom.aerodynamics.transonic_buildup import (  # noqa: E402
    calculate_clbclc, calculate_wbclb, setup2_step,
)

_SECOND = json.loads((_ROOT / 'tests' / 'fixtures' / 'probes' /
                      'second_level.json').read_text())


def _wbclb(c):
    return calculate_wbclb(
        c['alpha'], c['alpha_body'], {'sspn': c['sspn']}, c['body_x'],
        c['body_s'], c['xc'], c['taper'], c['aliw'], c['cl_w'], c['cla_w'],
        c['cl_b'], c['cla_b'], c['cd_w'], c['cd_b'], c['kwb'], c['kbw'],
        c['cla_wb'], {'clb_14': c['clbn14'], 'cna_14': c['cnam14'],
                      'clb_mfb': c['clblfb'], 'cla_mfb': c['clamfb']},
        c['mach'], c['mfb'], c['clwb'], c['cdwb'], c['clbb'])


@pytest.mark.parametrize("case", range(len(_SECOND['wbclb'])))
def test_wbclb_matches_compiled_routine(case):
    c, o = _SECOND['wbclb'][case]['inputs'], _SECOND['wbclb'][case]['outputs']
    r = _wbclb(c)
    np.testing.assert_allclose(r['cl'], o['CLWB'], rtol=1e-12, atol=1e-40)
    np.testing.assert_allclose(r['cd'], o['CDWB'], rtol=1e-12, atol=1e-40)
    np.testing.assert_allclose(r['clb'], o['CLBB'], rtol=1e-12, atol=1e-40)
    assert (r['kkwb'], r['kkbw'], r['clbcl'], r['ratio']) == pytest.approx(
        o['K'], rel=1e-12)
    np.testing.assert_allclose(r['ivbw'], o['IV'], rtol=1e-12, atol=1e-15)
    np.testing.assert_allclose(r['go2pav'], o['GO'], rtol=1e-12, atol=1e-15)


def test_wbclb_ratio_is_the_body_radius_bodowg_writes_back():
    """BODOWG assigns its radius argument, which is WBCLB's YB, so the
    ratio is sqrt(Smax/pi)/SSPN whatever the exposed-root offset."""
    c = _SECOND['wbclb'][0]['inputs']
    r = _wbclb(c)
    assert r['ratio'] == pytest.approx(
        np.sqrt(max(c["body_s"]) / PI) / c["sspn"], rel=1e-14)
    assert _SECOND['wbclb'][0]['outputs']['K'][3] == pytest.approx(
        r['ratio'], rel=1e-12)


def test_wbclb_probe_reaches_every_branch():
    results = [_wbclb(w['inputs']) for w in _SECOND['wbclb']]
    assert any(r['ratio'] >= 1 / 3 for r in results)
    assert any(r['ratio'] < 1 / 3 for r in results)
    assert any(any(v != 0.0 for v in r['ivbw']) for r in results)
    drags = [w['inputs']['cdwb'][1] for w in _SECOND['wbclb']]
    assert UNUSED in drags and any(d != UNUSED for d in drags)


def _block(step, name):
    cl, clb = _SECOND['blocks'][f'{step}{name}']
    block = {21 + k: v for k, v in enumerate(cl)}
    block.update({181 + k: v for k, v in enumerate(clb)})
    return block


@pytest.mark.parametrize("case", range(len(_SECOND['setup2'])))
def test_setup2_matches_compiled_routine(case):
    c = _SECOND['setup2'][case]['inputs']
    calls = _SECOND['setup2'][case]['calls']
    state = {
        'mach': c['mach'], 'flc2': c['flc2'], 'nalpha': c['nalpha'],
        'bo': c['bo'], 'wgpl': c['wgpl'], 'htpl': c['htpl'],
        'subson': False, 'transn': True, 'supers': False, 'done': False,
        'sec': {k: -5.0 for k in range(1, 24)},
        'tra6': c['mfbw'], 'trah6': c['mfbh'], 'b': [0.0, 0.0],
        'bht': [0.0, 0.0],
        'wingin': {68: c['win68'], 69: c['win69'], 21: 0.0, 41: 0.0},
        'htin': {68: c['hin68'], 69: c['hin69'], 21: 0.0, 41: 0.0},
    }
    nf = -1
    for step, call in enumerate(calls, start=1):
        state.update({
            'wbt67': 0.01 * step, 'stp155': -0.02 * step,
            'bw101': 0.08 + 0.001 * step, 'bh101': 0.05 + 0.001 * step,
            **{name.lower(): _block(step, name)
               for name in ('WING', 'HT', 'BW', 'BH')}})
        nf = setup2_step(nf, state)
        assert (nf, state['nalpha']) == (call['nf'], call['nalpha'])
        ours = [state['mach'], *state['b'], *state['bht'],
                state['wingin'][21], state['wingin'][41],
                state['htin'][21], state['htin'][41]]
        assert ours == pytest.approx(call['S'], rel=1e-12)
        sec = state['sec']
        assert [sec[k] for k in range(1, 15)] == pytest.approx(
            call['SEC'], rel=1e-12)
        assert [sec[k] for k in range(17, 24)] == pytest.approx(
            call['SEX'], rel=1e-12)
        assert [state['subson'], state['transn'], state['supers']] == \
            call['LOG']
    assert nf == -8 and state['done']


def test_clbclc_skips_zero_lift_and_unused_ratios():
    data = {21: 0.0, 22: 0.5, 23: 0.4, 181: 9.0, 182: -UNUSED, 183: 0.02}
    assert calculate_clbclc(data, 3) == pytest.approx(0.05)
    assert calculate_clbclc({21: 0.0, 181: 1.0}, 1) == UNUSED
