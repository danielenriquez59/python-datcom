"""
Regression tests for the block set-up routines (INIZ, INITZ1, INITZ2,
CLEARA, SYNDIM, SECI, SECO) and the executive sequencing (MAIN05, MAIN06,
MAIN07).

Source of truth: the matching files in datcom-legacy/datcom_2000.

SYNDIM, SECI and SECO are checked against compiled probes
(tools/probes/section_setup.py) on every word of their blocks, including
the words they leave alone.
"""

import json
import pathlib

import pytest

from pydatcom.core.executive import main05, main06, main07
from pydatcom.io.block_setup import (INITZ1_BLOCKS, cleara, initz1, initz2,
                                     iniz, seci, seco, syndim)
from pydatcom.utils.constants import UNUSED

_PROBE = json.loads((pathlib.Path(__file__).resolve().parent / 'fixtures' /
                     'probes' / 'section_setup.json').read_text())
_APPROX = dict(rel=1e-12, abs=1e-300)


@pytest.mark.parametrize("case", range(len(_PROBE['syndim']['cases'])))
def test_syndim_matches_compiled_routine(case):
    c, r = _PROBE['syndim']['cases'][case], _PROBE['syndim']['records'][case]
    bd = [0.0] + [0.03 * k for k in range(1, 101)]
    a = [0.0] + [0.01 * k for k in range(1, 196)]
    aht = [0.0] + [0.02 * k for k in range(1, 196)]
    bd[77], bd[33], bd[65], bd[74] = (c['bd77'], c['bd33'], c['bd65'],
                                      c['bd74'])
    a[62], aht[62] = c['a62'], c['aht62']
    syna = [0.0] * 20
    syna[1], syna[6], syna[8] = c['xcg'], c['xh'], c['alih']
    htin = [0.0, 0.0, 0.0, c['htin3'], c['htin4']]
    syndim(bd, a, aht, c['sspn'], c['sspne'], syna, htin)
    assert bd[1:101] == pytest.approx(r['BD'], **_APPROX)
    assert [a[173], aht[173]] == pytest.approx(r['ARM'], **_APPROX)


@pytest.mark.parametrize("case", range(len(_PROBE['seci']['cases'])))
def test_seci_matches_compiled_routine(case):
    c, r = _PROBE['seci']['cases'][case], _PROBE['seci']['records'][case]
    typein = [0.0] + [0.5 * k for k in range(1, 163)]
    for k, v in c['typein'].items():
        typein[int(k)] = v
    o = seci([0.0] + c['a'], typein, c['straight'])
    assert [o['atype'], o['l'], o['l']] == pytest.approx(r['OUT'])
    assert o['naca'] == pytest.approx(r['NACA'], **_APPROX)
    assert o['cla'] == pytest.approx(r['CLA'], **_APPROX)
    assert o['cbar'] == pytest.approx(r['CBAR'][0], **_APPROX)
    for key, name, mark in (('x', 'X', 1.), ('yu', 'YU', 4.),
                            ('yl', 'YL', 5.), ('thn', 'THN', 6.),
                            ('cam', 'CAM', 7.)):
        expected = [-k * mark for k in range(1, 61)]
        if key in o:
            expected[:len(o[key])] = o[key]
        assert expected == pytest.approx(r[name], **_APPROX), name


@pytest.mark.parametrize("case", range(len(_PROBE['seco']['cases'])))
def test_seco_matches_compiled_routine(case):
    c, r = _PROBE['seco']['cases'][case], _PROBE['seco']['records'][case]
    a = [0.0] + [0.25 * k for k in range(1, 163)]
    for k in c['unset']:
        a[k] = UNUSED
    seco(a, c['camber'], c['atype'], c['nmach'], _PROBE['seco']['section'])
    assert a[1:] == pytest.approx(r['A'], **_APPROX)


def test_initialisers():
    bwh, bwhv, bh = [0.0] * 381, [0.0] * 381, [0.0] * 381
    iniz(bwh, bwhv, bh)
    assert bh[144] == bwh[47] == bwhv[7] == UNUSED
    assert bh[145] == bwh[48] == bwhv[8] == 0.0
    blocks = initz1()
    assert {k: len(v) - 1 for k, v in blocks.items()} == INITZ1_BLOCKS
    assert all(set(v[1:]) == {UNUSED} for v in blocks.values())


def test_initz2_leaves_the_vertical_tail():
    blocks = [[1.0] * 401 for _ in range(8)]
    initz2(*blocks)
    powr, fcm, supdw, body, wing, ht, vt, vf = blocks
    assert powr[315] == fcm[282] == supdw[93] == UNUSED
    assert body[201] == wing[400] == ht[380] == vf[201] == UNUSED
    assert body[200] == 1.0 and set(vt[1:]) == {1.0}


def test_cleara_clears_only_the_last_ldum_word():
    out = cleara()
    assert out['ldum121'] == 0 and out['iovly'] == 999
    assert len(out['ifig']) == 20 and out['ifig'][0][0] == '    '


def test_main05_adds_a_zero_deflection_and_leaves_it():
    f = [0.0] * 117
    f[1], f[2], f[16] = 10.0, 20.0, 2.0
    calls = []
    ran = main05(True, False, False, False, f, calls.append)
    assert ran == calls == ['M36O44', 'M37O45', 'M38O46']
    assert f[16] == 2.0 and f[3] == 0.0
    f[74] = 1.0
    assert main05(True, True, True, True, f, lambda n: None) == [
        'M55O67', 'M52O64']
    assert main05(False, False, True, True, f, lambda n: None) == ['M38O46']


def test_main06_and_main07_sequences():
    assert main06(True, True, 5.0, 0.8, lambda n: None) == [
        'M36O44', 'M52O64', 'M37O45', 'M40O50']
    assert main06(False, True, 5.0, 1.2, lambda n: None) == [
        'M52O64', 'M40O50']
    assert main07(True, True, lambda n: None) == ['M41O51', 'M53O65']
