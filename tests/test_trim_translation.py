"""
Regression tests for DRAGFP, TRIMRT and TRIMR2 (overlay M38O46).

Source of truth: datcom-legacy/datcom_2000/dragfp.f, trimrt.f, trimr2.f.

Checked against a compiled probe (tools/probes/trim.py); the tables are
re-parsed from source.
"""

import json
import pathlib
import sys

import numpy as np
import pytest

from pydatcom.aerodynamics import trim as module
from pydatcom.aerodynamics.trim import (
    calculate_dragfp, calculate_trimr2, calculate_trimrt,
)

_ROOT = pathlib.Path(__file__).resolve().parent.parent
sys.path.insert(0, str(_ROOT / 'tools'))
from fortran_data import parse  # noqa: E402

_PROBE = json.loads(
    (_ROOT / 'tests' / 'fixtures' / 'probes' / 'trim.json').read_text())


def _dragfp(c):
    t = c['tail']
    return calculate_dragfp(c['wing']['b23'], c['flap'], c['wing'],
                            {'a3': t['a3'], 'a7': t['a7'],
                             'alpha': t['alpha'], 'epsilon': c['epsilon']},
                            c['sref'])


def _trimrt(c, dragfp):
    """TRIMRT after DRAGFP, whose DELCDM is TRIMRT's DCDMIN (both
    ``WING(231..)``)."""
    fl, b = c['flags'], c['blocks']
    flap = dict(c['flap'])
    if 'delcdm' in dragfp:
        flap['dcdmin'] = dragfp['delcdm']
    if fl['wgpl'] and fl['htpl']:
        untrimmed = dict(b['BWH'])
        if fl['vtpl']:
            untrimmed['cd'] = b['BWHV']['cd']
        aclmax = c['tail']['bht43']
    else:
        untrimmed = b['BW'] if fl['bo'] else b['WING']
        aclmax = c['wing']['b43']
    return calculate_trimrt(c['alpha'], c['epsilon'], untrimmed, aclmax,
                            flap, dragfp['dcdi'])


def _trimr2(c):
    t2 = dict(c['t2'], alpha_clmax=c['tail']['bht43'])
    return calculate_trimr2(c['alpha'], c['epsilon'], c['q_ratio'],
                            c['blocks']['BW'], t2, c['wbt'], c['position'],
                            c['sref'], c['cbarr'])


@pytest.mark.parametrize("case", range(len(_PROBE)))
def test_dragfp_matches_compiled_routine(case):
    c, o = _PROBE[case]['inputs'], _PROBE[case]['outputs']
    r = _dragfp(c)
    np.testing.assert_allclose(r['dcdi'], o['DCDI'], rtol=1e-11, atol=1e-16)
    if c['flap']['ftype'] <= 5.0:
        np.testing.assert_allclose(r['delcdm'], o['DELCDM'], rtol=1e-11)
        np.testing.assert_allclose(r['delcdf'], o['DELCDF'], rtol=1e-11,
                                   atol=1e-16)
        assert r['kprm'] == pytest.approx(o['KPRM'][0], rel=1e-12)
    else:
        # Not written: WING(231..) keeps the flap routines' DCDMIN.
        np.testing.assert_allclose(o['DELCDM'], c['flap']['dcdmin'])


@pytest.mark.parametrize("case", range(len(_PROBE)))
def test_trimrt_matches_compiled_routine(case):
    c, o = _PROBE[case]['inputs'], _PROBE[case]['outputs']
    r = _trimrt(c, _dragfp(c))
    na, n = len(c['alpha']), r['ntrim']
    assert o['TRM'][20] == n
    assert o['TRM'][21] == (-9.0 if r['tstop'] is None else r['tstop'])
    np.testing.assert_allclose(r['alpha'], o['TRM'][:na], rtol=1e-14)
    for key, tag in [('utcl', 'UTCL'), ('utcm', 'UTCM'), ('utcd', 'UTCD')]:
        np.testing.assert_allclose(r[key], o[tag], rtol=1e-14)
    for key, tag in [('deltat', 'DELTAT'), ('dclt', 'DCLT'),
                     ('clmaxt', 'CLMAXT'), ('cdmint', 'CDMINT'),
                     ('chdt', 'CHDT'), ('cdit', 'CDIT')]:
        values = r[key]
        assert len(values) in (0, n)
        np.testing.assert_allclose(values, o[tag][:len(values)], rtol=1e-11,
                                   atol=1e-15, err_msg=key)
        assert o[tag][len(values):] == [0.0] * (na - len(values))


@pytest.mark.parametrize("case", range(len(_PROBE)))
def test_trimr2_matches_compiled_routine(case):
    c, o = _PROBE[case]['inputs'], _PROBE[case]['outputs']
    r = _trimr2(c)
    assert o['TRM2'][20] == r['ntrim'] and o['TRM2'][21] == r['tstop']
    for key, tag in [('aliht', 'ALIHT'), ('cdhtrm', 'CDHTRM'),
                     ('cmhtrm', 'CMHTRM'), ('hmtrm', 'HMTRM'),
                     ('hmunt', 'HMUNT'), ('clwbt', 'CLWBT'),
                     ('cdwbt', 'CDWBT')]:
        for j, value in enumerate(r[key]):
            if value is not None:
                assert value == pytest.approx(o[tag][j], rel=1e-11,
                                              abs=1e-15), f"{key}[{j}]"
    for j, value in enumerate(r['clhtrm']):
        assert value == pytest.approx(o['TRM2'][j], rel=1e-11, abs=1e-15)
        if value != -1000.0:
            assert o['HT261'][j] == pytest.approx(value, rel=1e-11)


def test_probe_reaches_every_branch():
    rt = [_trimrt(p['inputs'], _dragfp(p['inputs'])) for p in _PROBE]
    assert {r['tstop'] for r in rt} >= {None, 1.0, 2.0}
    r2 = [_trimr2(p['inputs']) for p in _PROBE]
    assert {r['tstop'] for r in r2} >= {1.0, 2.0, 3.0}
    assert any(-1000.0 in r['aliht'] for r in r2)
    types = {p['inputs']['flap']['ftype'] for p in _PROBE}
    assert types >= {1.0, 2.0, 3.0, 5.0, 6.0}
    dcm = [p['inputs']['flap']['dcm'] for p in _PROBE]
    assert any(d[0] > d[-2] for d in dcm) and any(d[0] < d[-2] for d in dcm)


def test_trimrt_ignores_the_last_deflection():
    """NDELTA is reduced by one: changing the last deflection's moment
    increment changes nothing."""
    c = _PROBE[0]['inputs']
    dragfp = _dragfp(c)
    base = _trimrt(c, dragfp)
    flap = dict(c['flap'], dcm=c['flap']['dcm'][:-1] + [5.0])
    assert _trimrt(dict(c, flap=flap), dragfp)['deltat'] == base['deltat']


def test_dragfp_uses_the_tails_geometry_when_there_is_one():
    c = next(p['inputs'] for p in _PROBE if p['inputs']['tail']['a3'] != 1e-30)
    r = _dragfp(c)
    assert r['scale'] == c['tail']['a3'] / c['sref']
    assert r['aspect'] == c['tail']['a7']


def test_tables_match_the_source():
    source = parse('dragfp')
    for name, values in source.items():
        np.testing.assert_array_equal(getattr(module, '_' + name), values,
                                      err_msg=name)
