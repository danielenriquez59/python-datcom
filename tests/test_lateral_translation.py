"""
Regression tests for the subsonic lateral pass: M29O35, SUBLAT, M17O21.

Source of truth: datcom-legacy/datcom_2000/m29o35.f, sublat.f, m17o21.f.

Checked against a compiled probe (tools/probes/lateral.py) that runs M29O35
and then M17O21, which calls SUBLAT for the wing and again for the
horizontal tail, as the main program does.  Every derivative block, both
``/SBETA/`` arrays and the dihedral inputs M29O35 writes back are compared.
"""

import json
import pathlib
import sys

import numpy as np
import pytest

from pydatcom.aerodynamics import lateral as module
from pydatcom.aerodynamics.lateral import (
    calculate_m17o21, calculate_m29o35_block, calculate_sublat,
)
from pydatcom.utils.constants import UNUSED

_ROOT = pathlib.Path(__file__).resolve().parent.parent
sys.path.insert(0, str(_ROOT / 'tools'))
from fortran_data import parse  # noqa: E402

_PROBE = json.loads(
    (_ROOT / 'tests' / 'fixtures' / 'probes' / 'lateral.json').read_text())


def _ints(mapping):
    return {int(k): float(v) for k, v in mapping.items()}


def _full(stb):
    block = {k: 0.0 for k in range(1, 136)}
    block.update(stb)
    return block


def _run(c):
    """M29O35's two passes, SUBLAT twice, then M17O21."""
    f = c['flags']
    wingin, htin = _ints(c['wingin']), _ints(c['htin'])
    a, aht = _ints(c['a']), _ints(c['aht'])
    avt, avf = _ints(c['avt']), _ints(c['avf'])
    vtin, vfin = _ints(c['vtin']), _ints(c['vfin'])
    syna = _ints(c['syna'])
    vtin['type'], vfin['type'] = vtin[15], vfin[15]
    wing_geo = calculate_m29o35_block(
        wingin, a, avt, vtin, aht, syna, c['body'], c['bd1'], c['bd66'],
        a[10], syna[9], f['vertup'], f['wgpl'], f['bo'], f['vtpl'],
        f['htpl'], tail_pass=False)
    stb = _full(wing_geo['stb'])
    tail_geo = calculate_m29o35_block(
        htin, aht, avf, vfin, aht, syna, c['body'], c['bd1'], c['bd66'],
        a[10], syna[12], f['vertup'], f['htpl'], f['bo'], f['vfpl'],
        f['htpl'], tail_pass=True, stb123_wing=stb[123])
    stbh = _full(tail_geo['stb'])
    wingin.update(wing_geo['dihedral'])
    htin.update(tail_geo['dihedral'])
    wingin['type'], htin['type'] = wingin[15], htin[15]
    flight = {'mach': c['mach'], 'reynolds_per_length': c['rn'],
              'alpha': c['alpha']}
    common = dict(aht=aht, htin=htin, tvtin=c['tvtin'], syna=syna,
                  flight=flight, bd1=c['bd1'], sref=c['sref'],
                  cbarr=c['cbarr'], blref=c['blref'])
    b = {2: c['b']['2']}
    b.update({3 + j: v for j, v in enumerate(c['wing_cl0'])})
    bht = {2: c['bht']['2']}
    bht.update({3 + j: v for j, v in enumerate(c['ht_cl0'])})
    wing = None
    if f['wgpl'] or f['vtpl'] or f['tvtpan']:
        wing = calculate_sublat(
            stb, c['wing'], wingin, a, b, c['bw_cl'],
            c['body_lat']['cla'], avt=avt, vtin=vtin, vt_cla=vtin[21],
            flags={k: f[k] for k in ('wgpl', 'bo', 'vtpl', 'tvtpan', 'htpl',
                                     'transn')}, ity=0, **common)
    tail = None
    if f['htpl'] or f['vfpl']:
        tail = calculate_sublat(
            stbh, c['ht'], htin, aht, bht, c['bh_cl'], c['body_lat']['cla'],
            avt=avf, vtin=vfin, vt_cla=vfin[21],
            flags={'wgpl': f['htpl'], 'bo': f['bo'], 'vtpl': f['vfpl'],
                   'tvtpan': False, 'htpl': False, 'transn': f['transn']},
            ity=1, **common)
    combos = calculate_m17o21(len(c['alpha']), c['body_lat'], wing, tail,
                              f['transn'])
    return wing, tail, combos, wingin


def _row(block, n):
    if block is None:
        return np.zeros(2 + n)
    cyb = np.asarray(block['cyb'], dtype=float).ravel()[0]
    cnb = np.asarray(block['cnb'], dtype=float).ravel()[0]
    clb = np.asarray(block.get('clb', np.zeros(n)), dtype=float)
    return np.concatenate([[cyb, cnb], clb])


def _get(result, key):
    return None if result is None else result.get(key)


@pytest.mark.parametrize("case", range(len(_PROBE)))
def test_matches_compiled_routines(case):
    probe = _PROBE[case]
    c, out = probe['inputs'], probe['outputs']
    wing, tail, combos, wingin = _run(c)
    n = len(c['alpha'])
    transn = c['flags']['transn']
    expected = {
        'WING': _row(_get(wing, 'surface'), n),
        'HT': _row(_get(tail, 'surface'), n),
        'VT': _row(_get(wing, 'vertical'), n),
        'VF': _row(_get(tail, 'vertical'), n),
        'BW': _row(_get(wing, 'combination'), n),
        'BH': _row(_get(tail, 'combination'), n),
        'BV': _row(combos['bv'], n), 'BWH': _row(combos['bwh'], n),
        'BWV': _row(combos['bwv'], n), 'BWHV': _row(combos['bwhv'], n),
    }
    if transn:
        # M17O21's closing overwrite of the first CY_beta and Cn_beta.
        for key in ('WING', 'HT', 'VT', 'VF'):
            expected[key][:2] = UNUSED
    for key, mine in expected.items():
        # Blocks a transonic SUBLAT leaves unset keep the source's zeros.
        np.testing.assert_allclose(mine, out[key], rtol=1e-9, atol=1e-14,
                                   err_msg=key)
    if wing is not None:
        _compare_stb(wing['stb'], out['STB'], 'STB')
    if tail is not None:
        _compare_stb(tail['stb'], out['STBH'], 'STBH')
    np.testing.assert_allclose([wingin[12], wingin[13], wingin[14]],
                               out['DIH'], rtol=1e-12, atol=1e-30)


def _compare_stb(mine, reference, name):
    for index in range(1, 136):
        value = mine.get(index, 0.0)
        assert value == pytest.approx(reference[index - 1], rel=1e-9,
                                      abs=1e-14), f'{name}({index})'


def test_probe_reaches_every_branch():
    runs = [(p['inputs'], _run(p['inputs'])) for p in _PROBE]
    flags = [c['flags'] for c, _ in runs]
    assert any(not f['bo'] for f in flags)
    assert any(f['transn'] for f in flags)
    assert any(f['tvtpan'] for f in flags)
    assert any(not f['htpl'] for f in flags)
    inputs = [c for c, _ in runs]
    assert any(c['wingin']['12'] != UNUSED for c in inputs)       # non-uniform
    assert any('zu' in c['body'] for c in inputs)
    assert any(c['a']['120'] < 1.0 for c in inputs)
    assert {c['wingin']['15'] for c in inputs} >= {1.0, 2.0, 3.0}
    assert any(c['vtin']['15'] != 1.0 for c in inputs)
    assert any(r[0] is not None and 'twin_correction' in r[0]
               for _, r in runs)
    assert any(c['syna']['7'] < 0.0 for c in inputs)


def test_every_table_value_matches_the_source():
    source = parse('sublat')
    names = ['X327', 'X127', 'X227', 'Y27', 'X128A', 'X228A', 'Y28A',
             'X128B', 'X228B', 'Y28B', 'X329', 'X129', 'X229', 'Y29',
             'X130A', 'X230A', 'Y30A', 'X130B', 'X230B', 'Y30B', 'X431',
             'X331', 'X131', 'X231', 'X1526', 'X2526', 'Y526', 'X158A',
             'X258A', 'Y58A', 'X158B', 'X258B', 'Y58B', 'X158C', 'X258C',
             'Y58C', 'X122A', 'X222A', 'Y22A', 'X122B', 'X222B', 'Y22B',
             'X5322C', 'Y5322C', 'X5322D', 'Y5322D', 'X5324A', 'Y5324A',
             'X124B', 'X224B', 'Y24B', 'X124C', 'X224C', 'Y24C']
    for name in names:
        np.testing.assert_array_equal(getattr(module, f'_{name}'),
                                      source[name], err_msg=name)
    np.testing.assert_array_equal(
        module._Y31, source['T31A'] + source['T31B'] + source['T31C'])


def test_no_body_skips_the_vertical_tail():
    probe = next(p for p in _PROBE if not p['inputs']['flags']['bo'])
    wing, _, _, _ = _run(probe['inputs'])
    assert 'vertical' not in wing and 'surface' in wing


def test_lone_phif_doubles_the_vertical_tail():
    """PHIF set and PHIV not: the wing pass scales by 2*cos(UNUSED)**2."""
    probe = next(p for p in _PROBE
                 if p['inputs']['syna']['19'] != UNUSED
                 and p['inputs']['syna']['18'] == UNUSED)
    wing, _, _, _ = _run(probe['inputs'])
    assert wing['twin_correction'] == pytest.approx(2.0)


def test_directional_stability_signs():
    """A vertical tail gives positive Cn_beta and negative CY_beta; the
    body alone is directionally unstable."""
    probe = _PROBE[0]
    wing, _, combos, _ = _run(probe['inputs'])
    assert wing['vertical']['cyb'] < 0.0
    assert wing['vertical']['cnb'] > 0.0
    assert wing['combination']['cnb'] < 0.0
    assert combos['bwv']['cnb'] > wing['combination']['cnb']
