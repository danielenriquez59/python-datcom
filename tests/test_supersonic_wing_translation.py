"""
Regression tests for SUPLNG, the supersonic wing lift, pitching moment and
drag, its horizontal-tail twin SUPLTG, and their overlays M27O33 and
M22O26.

Source of truth: datcom-legacy/datcom_2000/suplng.f, supltg.f, m27o33.f
and m22o26.f.

Checked against compiled probes (tools/probes/suplng.py and supltg.py) on
every word of ``/SUPWH/`` (bar the LOGICAL ``DETACH``, checked on its
own), the surface's geometry block and its result words 1-200.  The
SUPLTG cases run in one program, so the replay carries the saved ``RACH``.
"""

import json
import pathlib
import sys

import numpy as np
import pytest

from pydatcom.aerodynamics import supersonic_wing as module
from pydatcom.aerodynamics.supersonic_wing import (
    calculate_suplng, calculate_supltg, m22o26, m27o33)

_ROOT = pathlib.Path(__file__).resolve().parent.parent
sys.path.insert(0, str(_ROOT / 'tools'))
from fortran_data import parse  # noqa: E402


def _load(name):
    return json.loads((_ROOT / 'tests' / 'fixtures' / 'probes' /
                       f'{name}.json').read_text())


_SUPLNG, _M27O33 = _load('suplng'), _load('m27o33')
_SUPLTG, _M22O26 = _load('supltg'), _load('m22o26')


def _replay(fn, probe):
    state, out = {'rach': 0.0}, []
    for p in probe:
        r = fn(dict(p['inputs'], state=state))
        state = r['state']
        out.append(r)
    return out


_TAIL = _replay(calculate_supltg, _SUPLTG)
_TAIL_OVERLAY = _replay(m22o26, _M22O26)
_UNUSED_TABLES = {'D455', 'T455', 'I27'}


def _check(r, o):
    slg = np.delete(np.array(r['slg'][1:]), 93)
    np.testing.assert_allclose(slg, np.delete(np.array(o['SLG']), 93),
                               rtol=1e-9, atol=1e-14, err_msg='SLG')
    np.testing.assert_allclose(r['a'][1:], o['A'], rtol=1e-9, atol=1e-14,
                               err_msg='A')
    np.testing.assert_allclose(r['wing'][1:201], o['HT'], rtol=1e-9,
                               atol=1e-14, err_msg='WING')


@pytest.mark.parametrize("case", range(len(_SUPLNG)))
def test_suplng_matches_compiled_routine(case):
    p = _SUPLNG[case]
    r = calculate_suplng(p['inputs'])
    _check(r, p['outputs'])
    if p['outputs']['DET'][0]:
        assert r['detach']


@pytest.mark.parametrize("case", range(len(_M27O33)))
def test_m27o33_matches_compiled_overlay(case):
    p = _M27O33[case]
    _check(m27o33(p['inputs']), p['outputs'])


@pytest.mark.parametrize("case", range(len(_SUPLTG)))
def test_supltg_matches_compiled_routine(case):
    _check(_TAIL[case], _SUPLTG[case]['outputs'])
    assert _TAIL[case]['detach'] == bool(_SUPLTG[case]['outputs']['DET'][0])


@pytest.mark.parametrize("case", range(len(_M22O26)))
def test_m22o26_matches_compiled_overlay(case):
    _check(_TAIL_OVERLAY[case], _M22O26[case]['outputs'])


def test_cranked_tail_reuses_the_wave_drag_arg():
    """A cranked tail's drag due to lift is formed with the wave-drag
    ``ARG``, not ``(1+P)/(PI*A*P)``, and ``P``/``DRAGC`` stand as given."""
    k = next(i for i, p in enumerate(_SUPLTG)
             if not p['inputs']['straight'])
    c, r = _SUPLTG[k]['inputs'], _TAIL[k]
    assert r['slg'][81] == c['slg'][80] and r['slg'][82] == c['slg'][81]
    textbook = (1. + r['slg'][82]) / (3.141592654 * r['a'][7] *
                                      r['slg'][82])
    cl = r['wing'][21]
    implied = r['slg'][53] / (r['slg'][81] * c['sw'] / r['a'][3] * cl ** 2)
    assert implied != pytest.approx(textbook)


def test_probe_reaches_every_branch():
    runs = [(p['inputs'], calculate_suplng(p['inputs'])) for p in _SUPLNG]
    assert any(not c['straight'] for c, _ in runs)
    assert any(r['detach'] for _, r in runs)
    # The detachment search succeeds at one degree in one case and gives
    # up (90 degrees) in another; the transition band is populated.
    detang = {round(r['slg'][115], 6) for c, r in runs
              if c['straight'] and not r['detach'] and
              r['slg'][2] > 1. and r['slg'][93] >= 1.}
    assert round(1. / 57.2957795, 6) in detang
    assert round(90. / 57.2957795, 6) in detang


def test_detachment_search_tries_one_degree():
    """The search counts degrees against a limit in radians, so with
    angles of attack under a radian it stops after one degree."""
    c = dict(_SUPLNG[3]['inputs'])
    r = calculate_suplng(c)
    assert r['slg'][115] == pytest.approx(90. / 57.2957795)


def test_tables_match_the_source():
    for name, values in parse('suplng').items():
        if name in _UNUSED_TABLES:
            assert not hasattr(module, '_' + name)
            continue
        np.testing.assert_array_equal(getattr(module, '_' + name), values,
                                      err_msg=name)
