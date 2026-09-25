"""
Regression tests for the supersonic rate derivatives: SUPCLD/SUPHLD
(CL-alpha-dot).

Source of truth: datcom-legacy/datcom_2000/supcld.f, suphld.f.

Checked against compiled probes (tools/probes/supcld.py), both routines
on their own COMMON blocks, on every word of ``DYN`` and ``A`` and the
result word, with CALCA linked.
"""

import json
import pathlib
import sys

import numpy as np
import pytest

from pydatcom.aerodynamics import supersonic_rate as module
from pydatcom.aerodynamics.supersonic_rate import calculate_supcld

_ROOT = pathlib.Path(__file__).resolve().parent.parent
sys.path.insert(0, str(_ROOT / 'tools'))
from fortran_data import parse  # noqa: E402

_CLD = json.loads((_ROOT / 'tests' / 'fixtures' / 'probes' /
                   'supcld.json').read_text())


def _inputs(c):
    a = {k: 0.03 * k for k in range(1, 196)}
    a.update({int(k): v for k, v in c['a'].items()})
    wingin = {k: 0.04 * k for k in range(1, 102)}
    wingin.update({int(k): v for k, v in c['wingin'].items()})
    slg = {k: 0.02 * k for k in range(1, 142)}
    slg.update({int(k): v for k, v in c['slg'].items()})
    dyn = {k: 0.01 * k for k in range(1, 214)}
    dyn[9], dyn[10] = 0.8, 1.3
    return dict(c, a=a, wingin=wingin, slg=slg, dyn=dyn,
                result={241: -7.0})


@pytest.mark.parametrize("name", ['supcld', 'suphld'])
@pytest.mark.parametrize("case", range(len(_CLD['cases'])))
def test_supcld_matches_compiled_routine(name, case):
    r = calculate_supcld(_inputs(_CLD['cases'][case]))
    o = _CLD['records'][name][case]
    np.testing.assert_allclose([r['dyn'][k] for k in range(1, 214)],
                               o['DYN'], rtol=1e-9, atol=1e-14)
    np.testing.assert_allclose([r['a'][k] for k in range(1, 196)], o['A'],
                               rtol=1e-9, atol=1e-14)
    assert r['result'][241] == pytest.approx(o['CLAD'][0], rel=1e-9)


def test_low_taper_is_left_at_a_quarter():
    k = next(i for i, c in enumerate(_CLD['cases'])
             if 0 < c['a']['27'] < 0.25 and c['slg']['1'] / c['a']['62'] < 1)
    assert calculate_supcld(_inputs(_CLD['cases'][k]))['a'][27] == 0.25


def test_tables_match_the_source():
    for name, values in parse('supcld').items():
        if name.startswith('I'):
            continue
        np.testing.assert_array_equal(getattr(module, '_' + name), values,
                                      err_msg=name)


# --- SUPCMD / SUPHMD (tools/probes/supcmd.py) ------------------------------

from pydatcom.aerodynamics.supersonic_rate import calculate_supcmd  # noqa

_CMD = json.loads((_ROOT / 'tests' / 'fixtures' / 'probes' /
                   'supcmd.json').read_text())


@pytest.mark.parametrize("name", ['supcmd', 'suphmd'])
@pytest.mark.parametrize("case", range(len(_CMD['cases'])))
def test_supcmd_matches_compiled_routine(name, case):
    c = _CMD['cases'][case]
    a = {k: 0.03 * k for k in range(1, 196)}
    a.update({int(k): v for k, v in c['a'].items()})
    a[173] = 1.25
    slg = {k: 0.02 * k for k in range(1, 142)}
    slg.update({int(k): v for k, v in c['slg'].items()})
    dyn = {k: 0.01 * k for k in range(1, 214)}
    r = calculate_supcmd(dict(c, a=a, slg=slg, dyn=dyn,
                              result={241: -7.0, 261: 3.5}))
    o = _CMD['records'][name][case]
    np.testing.assert_allclose([r['dyn'][k] for k in range(1, 214)],
                               o['DYN'], rtol=1e-9, atol=1e-14)
    np.testing.assert_allclose([r['a'][k] for k in range(1, 196)], o['A'],
                               rtol=1e-9, atol=1e-14)
    assert [r['result'][241], r['result'][261]] == pytest.approx(
        o['RES'], rel=1e-9)


def test_supcmd_tables_match_the_source():
    for name, values in parse('supcmd').items():
        if name.startswith('I'):
            continue
        np.testing.assert_array_equal(getattr(module, '_MD_' + name), values,
                                      err_msg=name)


# --- SUPCMQ / SUPHMQ (tools/probes/supcmq.py) ------------------------------

from pydatcom.aerodynamics.supersonic_rate import calculate_supcmq  # noqa

_CMQ = json.loads((_ROOT / 'tests' / 'fixtures' / 'probes' /
                   'supcmq.json').read_text())


def _cmq_inputs(c):
    a = {k: 0.03 * k for k in range(1, 196)}
    a.update({3: 280.0, 7: c['ar'], 10: 11.0, 16: 7.8, 27: c['taper'],
              29: 12.0, 50: 0.3, 62: c['tanle'], 67: 0.9, 68: 0.5,
              173: 1.25})
    slg = {k: 0.02 * k for k in range(1, 142)}
    slg.update({7: 3.1, 134: 0.42, 135: -0.12})
    tra = {k: 0.005 * k for k in range(1, 109)}
    tra.update({int(k): v for k, v in _CMQ['tra'].items()})
    return {'straight': c['straight'], 'transn': c['transn'],
            'mach': c['mach'], 'sr': 320.0, 'cbarr': 8.5, 'a': a,
            'wingin': {69: 0.1}, 'slg': slg,
            'dyn': {k: 0.01 * k for k in range(1, 214)}, 'tra': tra,
            'result': {101: 0.065, 201: 1.7, 221: -9.0}}


@pytest.mark.parametrize("name", ['supcmq', 'suphmq'])
@pytest.mark.parametrize("case", range(len(_CMQ['cases'])))
def test_supcmq_matches_compiled_routine(name, case):
    r = calculate_supcmq(_cmq_inputs(_CMQ['cases'][case]),
                         tail=(name == 'suphmq'))
    o = _CMQ['records'][name][case]
    np.testing.assert_allclose([r['dyn'][k] for k in range(1, 214)],
                               o['DYN'], rtol=1e-9, atol=1e-14)
    np.testing.assert_allclose([r['a'][k] for k in range(1, 196)], o['A'],
                               rtol=1e-9, atol=1e-14)
    assert r['result'][221] == pytest.approx(o['CMQ'][0], rel=1e-9)


def test_supcmq_leaves_a_zero_sweep_nudged():
    """On the supersonic path the zero leading-edge tangent stays 1E-5."""
    k = next(i for i, c in enumerate(_CMQ['cases'])
             if c['tanle'] == 0.0 and not c['transn'])
    assert calculate_supcmq(_cmq_inputs(_CMQ['cases'][k]))['a'][62] == \
        .00001


def test_supcmq_tables_match_the_source():
    for name, values in parse('supcmq').items():
        if name.startswith('I'):
            continue
        np.testing.assert_array_equal(getattr(module, '_MQ_' + name), values,
                                      err_msg=name)


# --- SUPPAW / SUPPAH (tools/probes/suppaw.py) ------------------------------

from pydatcom.aerodynamics.supersonic_rate import calculate_suppaw  # noqa

_PAW = json.loads((_ROOT / 'tests' / 'fixtures' / 'probes' /
                   'suppaw.json').read_text())


def _paw_inputs(c):
    a = {k: 0.03 * k for k in range(1, 196)}
    a.update({3: 280.0, 7: c['ar'], 10: 11.0, 16: 7.8, 27: c['taper'],
              29: 12.0, 62: c['tanle']})
    slg = {k: 0.02 * k for k in range(1, 142)}
    slg.update({1: c['beta'] or 0.7, 7: 3.1, 134: 0.42, 135: -0.12})
    return {'straight': c['straight'], 'transn': c['transn'], 'sr': 320.0,
            'cbarr': 8.5, 'a': a, 'slg': slg,
            'dyn': {k: 0.01 * k for k in range(1, 214)},
            'result': {201: 1.7}}


@pytest.mark.parametrize("name", ['suppaw', 'suppah'])
@pytest.mark.parametrize("case", range(len(_PAW['cases'])))
def test_suppaw_matches_compiled_routine(name, case):
    r = calculate_suppaw(_paw_inputs(_PAW['cases'][case]))
    o = _PAW['records'][name][case]
    for block, size, tag in (('dyn', 213, 'DYN'), ('a', 195, 'A'),
                             ('slg', 141, 'SLG')):
        np.testing.assert_allclose([r[block][k] for k in range(1, size + 1)],
                                   o[tag], rtol=1e-9, atol=1e-14,
                                   err_msg=tag)
    assert r['result'][201] == pytest.approx(o['CLQ'][0], rel=1e-9)


def test_suppaw_transonic_overwrites_beta():
    k = next(i for i, c in enumerate(_PAW['cases']) if c['transn'])
    r = calculate_suppaw(_paw_inputs(_PAW['cases'][k]))
    assert r['slg'][1] == pytest.approx(0.663324958)


def test_suppaw_tables_match_the_source():
    for name, values in parse('suppaw').items():
        if name.startswith('I'):
            continue
        np.testing.assert_array_equal(getattr(module, '_QW_' + name), values,
                                      err_msg=name)
