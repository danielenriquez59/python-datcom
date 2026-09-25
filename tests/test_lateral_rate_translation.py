"""
Regression tests for the lateral rate derivatives: SUPHYW (supersonic
horizontal tail).

Source of truth: datcom-legacy/datcom_2000/suphyw.f.

Checked against a compiled probe (tools/probes/suphyw.py) on ``DYN``
204-213 and ``HT`` 280-340, over leading-edge regimes and taper ratios in
one program (INTEP3 keeps its chart pair).
"""

import json
import pathlib
import sys

import numpy as np
import pytest

from pydatcom.aerodynamics import lateral_rate as module
from pydatcom.aerodynamics.lateral_rate import calculate_suphyw

_ROOT = pathlib.Path(__file__).resolve().parent.parent
sys.path.insert(0, str(_ROOT / 'tools'))
from fortran_data import parse  # noqa: E402

_HYW = json.loads((_ROOT / 'tests' / 'fixtures' / 'probes' /
                   'suphyw.json').read_text())


def _replay():
    state, out = {}, []
    for c in _HYW['cases']:
        r = calculate_suphyw({
            'mach': c['mach'], 'alpha': _HYW['alpha'], 'nalpha': 4,
            'sr': 320.0, 'blref': 36.0, 'htpl': c['htpl'], 'bo': c['bo'],
            'a': {4: 60.0, 58: c['lamle'], 62: c['tanle'], 74: 0.5,
                  118: c['taper'], 120: c['aw']},
            'htin': {3: c['sspne'], 4: 6.5, 6: 5.0},
            'slg': {1: c['beta'], 3: 0.93}, 'syna': {1: 22.0},
            'dyn': {k: 0.002 * k for k in range(1, 214)},
            'ht': {k: 0.001 * k for k in range(1, 381)}, 'state': state})
        state = r['state']
        out.append(r)
    return out


_RESULTS = _replay()


@pytest.mark.parametrize("case", range(len(_HYW['cases'])))
def test_suphyw_matches_compiled_routine(case):
    r, o = _RESULTS[case], _HYW['records'][case]
    np.testing.assert_allclose([r['dyn'][k] for k in range(204, 214)],
                               o['DYN'], rtol=1e-9, atol=1e-14)
    np.testing.assert_allclose([r['ht'][k] for k in range(280, 341)],
                               o['HT'], rtol=1e-9, atol=1e-14)


def test_subsonic_edge_cnp_uses_the_stale_word():
    """CN-p is built from DYN(212), which SUPHYW never sets."""
    k = next(i for i, c in enumerate(_HYW['cases'])
             if c['beta'] / c['tanle'] < 1 and not
             _RESULTS[i]['returned_early'])
    r, alpha = _RESULTS[k], _HYW['alpha']
    assert r['dyn'][212] == pytest.approx(0.002 * 212)
    assert r['ht'][323] == pytest.approx(
        (r['dyn'][212] - r['ht'][281]) * alpha[2] / 57.2957795)


def test_tables_match_the_source():
    for name, values in parse('suphyw').items():
        if name.startswith('I') or name in ('X13170', 'X23170', 'Y23170'):
            continue
        np.testing.assert_array_equal(getattr(module, '_' + name), values,
                                      err_msg=name)


# --- SUBHYW (tools/probes/subhyw.py) ---------------------------------------

import math  # noqa: E402

from pydatcom.aerodynamics.lateral_rate import calculate_subhyw  # noqa

_HB = json.loads((_ROOT / 'tests' / 'fixtures' / 'probes' /
                  'subhyw.json').read_text())


def _subhyw_replay():
    state, out = {}, []
    for c in _HB['cases']:
        ht = {k: 0.001 * k for k in range(1, 381)}
        for k, v in enumerate(_HB['cl'][c['cl']]):
            ht[21 + k] = v
        for k, v in enumerate(_HB['cd']):
            ht[1 + k] = v
        r = calculate_subhyw({
            'mach': c['mach'], 'alpha': _HB['alpha'], 'nalpha': 5,
            'sr': 320.0, 'blref': 36.0, 'htpl': c['htpl'], 'bo': c['bo'],
            'im': 1,
            'a': {4: 60.0, 27: c['taper'], 64: c['sweep'],
                  68: math.tan(math.radians(c['sweep'])), 120: c['ar'],
                  197: c['beta']},
            'htin': {3: c['sspne'], 4: 6.5, 21: 0.105},
            'syna': {5: 1.2, 7: 2.0}, 'stb': {122: 3.0},
            'dyn': {k: 0.002 * k for k in range(1, 214)}, 'ht': ht,
            'state': state})
        state = r['state']
        out.append(r)
    return out


_HB_RESULTS = _subhyw_replay()


@pytest.mark.parametrize("case", range(len(_HB['cases'])))
def test_subhyw_matches_compiled_routine(case):
    r, o = _HB_RESULTS[case], _HB['records'][case]
    np.testing.assert_allclose([r['dyn'][k] for k in range(40, 100)],
                               o['DYN'], rtol=1e-9, atol=1e-14)
    np.testing.assert_allclose([r['ht'][k] for k in range(281, 301)],
                               o['HT'], rtol=1e-9, atol=1e-14)


# --- SUBRYW (tools/probes/subryw.py) ---------------------------------------

from pydatcom.aerodynamics.lateral_rate import calculate_subryw  # noqa

_RB = json.loads((_ROOT / 'tests' / 'fixtures' / 'probes' /
                  'subryw.json').read_text())


def _subryw_replay():
    state, out = {}, []
    for c in _RB['cases']:
        wing = {k: 0.001 * k for k in range(1, 401)}
        for k, v in enumerate(_RB['cl'][c['cl']]):
            wing[21 + k] = v
        for k, v in enumerate(_RB['cd']):
            wing[1 + k] = v
        wing[26] = c['next_cl']
        dyn = {k: 0.002 * k for k in range(1, 214)}
        dyn[21] = -0.08
        sweep = math.radians(c['sweep'])
        r = calculate_subryw({
            'mach': c['mach'], 'alpha': _RB['alpha'], 'nalpha': 5,
            'sr': 320.0, 'blref': 36.0, 'cbarr': 8.5, 'wgpl': c['htpl'],
            'bo': c['bo'], 'im': 1,
            'a': {4: 60.0, 16: 7.5, 27: c['taper'], 64: c['sweep'],
                  67: math.cos(sweep), 68: math.tan(sweep), 120: c['ar'],
                  197: c['beta']},
            'htin': {3: c['sspne'], 4: 6.5, 11: -2.0, 21: 0.105},
            'syna': {3: 2.0, 5: 1.2}, 'stb': {122: 3.0}, 'dyn': dyn,
            'wing': wing, 'state': state})
        state = r['state']
        out.append(r)
    return out


_RB_RESULTS = _subryw_replay()


@pytest.mark.parametrize("case", range(len(_RB['cases'])))
def test_subryw_matches_compiled_routine(case):
    r, o = _RB_RESULTS[case], _RB['records'][case]
    np.testing.assert_allclose([r['dyn'][k] for k in range(40, 210)],
                               o['DYN'], rtol=1e-9, atol=1e-14)
    np.testing.assert_allclose([r['wing'][k] for k in range(281, 361)],
                               o['HT'], rtol=1e-9, atol=1e-14)


def test_subryw_zero_lift_term_follows_the_word_after_the_curve():
    k = next(i for i, c in enumerate(_RB['cases']) if c['next_cl'] == 0.0)
    assert _RB_RESULTS[k]['dyn'][46] == 0.0


# --- SUPRYW, SUPHYW's wing twin (tools/probes/suphyw.py) -------------------

from pydatcom.aerodynamics.lateral_rate import calculate_supryw  # noqa


def test_supryw_matches_compiled_routine():
    state = {}
    for c, o in zip(_HYW['cases'], _HYW['wing_records']):
        r = calculate_supryw({
            'mach': c['mach'], 'alpha': _HYW['alpha'], 'nalpha': 4,
            'sr': 320.0, 'blref': 36.0, 'wgpl': c['htpl'], 'bo': c['bo'],
            'a': {4: 60.0, 58: c['lamle'], 62: c['tanle'], 74: 0.5,
                  118: c['taper'], 120: c['aw']},
            'htin': {3: c['sspne'], 4: 6.5, 6: 5.0},
            'slg': {1: c['beta'], 3: 0.93}, 'syna': {1: 22.0},
            'dyn': {k: 0.002 * k for k in range(1, 214)},
            'ht': {k: 0.001 * k for k in range(1, 401)}, 'state': state})
        state = r['state']
        np.testing.assert_allclose(
            [r['dyn'][k] for k in range(204, 214)] +
            [r['ht'][k] for k in range(280, 341)], o['DYN'] + o['HT'],
            rtol=1e-9, atol=1e-14)
