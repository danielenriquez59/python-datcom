"""
Regression tests for overlay M10O12: WGEOTL, WBTAIL and its closing pass.

Source of truth: datcom-legacy/datcom_2000/m10o12.f, wgeotl.f, wbtail.f.

Checked against a compiled probe of the overlay itself
(tools/probes/m10o12.py), so the chaining of the three routines is tested
as the program runs it, and every table against a re-parse of the FORTRAN.
"""

import json
import pathlib
import sys

import numpy as np
import pytest

from pydatcom.aerodynamics import wbtail as module
from pydatcom.aerodynamics.wbtail import (
    NOT_AVAILABLE, calculate_m10o12, calculate_wbtail, calculate_wgeotl,
)

_ROOT = pathlib.Path(__file__).resolve().parent.parent
sys.path.insert(0, str(_ROOT / 'tools'))
from fortran_data import parse  # noqa: E402

_PROBE = json.loads(
    (_ROOT / 'tests' / 'fixtures' / 'probes' / 'm10o12.json').read_text())

_COMPONENTS = ('cd', 'cl', 'cm', 'cn', 'ca', 'cla', 'cma')


def _run(c):
    """The overlay's sequence: WGEOTL, the canard span override, WBTAIL,
    then the closing pass."""
    v = c['vertical']
    wing_body = {k: c['wing_body'][k] for k in ('cd', 'cl', 'cm', 'cla', 'cma')}
    if not c['htpl']:
        return None, None, calculate_m10o12(c['alpha_deg'], wing_body,
                                            v['vt1'], v['vf1'])
    wing = dict(c['wing'], twash=c['twash'])
    geo = calculate_wgeotl(c['alpha_deg'], wing, c['tail'],
                           c['synthesis']['aliw'])
    downwash = dict(c['downwash'], ali=geo['ali'])
    if geo['vortex_span'] is not None:
        downwash['vortex_span'] = [geo['vortex_span']] * len(c['alpha_deg'])
    bw = dict(c['wing_body'], cla=c['wing_body']['cla'][0],
              cma=c['wing_body']['cma'][0])
    tail = calculate_wbtail(c['alpha_deg'], c['tail'], c['tail_alone'], bw,
                            downwash, wing, c['synthesis'], c['body'],
                            c['sref'], c['cbarr'], v['dvt20'] + v['dvf20'])
    return geo, tail, calculate_m10o12(c['alpha_deg'], wing_body, v['vt1'],
                                       v['vf1'], tail)


def _close(actual, expected, what):
    np.testing.assert_allclose(actual, expected, rtol=1e-9, atol=1e-12,
                               err_msg=what)


@pytest.mark.parametrize("case", range(len(_PROBE)))
def test_matches_compiled_overlay(case):
    """BWH, BWHV, BWV, the WGEOTL factors, WBT, and the tail geometry."""
    probe = _PROBE[case]
    out = probe['outputs']
    geo, tail, result = _run(probe['inputs'])
    for block in ('bwh', 'bwhv', 'bwv'):
        if block in result:
            for comp in _COMPONENTS:
                _close(result[block][comp], out[block.upper() + comp.upper()],
                       f'{block} {comp}')
    if geo is None:
        return
    n = len(probe['inputs']['alpha_deg'])
    _close(geo['ali'], out['FACT42'], 'FACT(42)')
    if tail['canard']:
        _close(tail['fact101'], out['FACT102'], 'FACT(102)')
        _close(tail['fact121'], out['FACT122'], 'FACT(122)')
        _close([geo['vortex_span']] * n, out['FACT82'], 'FACT(82)')
    g = tail['geometry']
    _close([g[k] for k in ('bd8', 'bd30', 'bd31', 'bd58', 'bd63', 'bd64',
                           'bd84', 'bd761', 'bd762')], out['GEOM'], 'BD')
    _close(tail['cd_slope'], out['BD95'], 'BD(95)')
    wbt = out['WBT']
    for key, slot in [('khb', 1), ('kbh', 2), ('wbt3', 3), ('wbt4', 4),
                      ('wbt66', 66), ('wbt67', 67), ('body_radius', 108),
                      ('khb_incidence', 150), ('kbh_incidence', 151)]:
        assert tail[key] == pytest.approx(wbt[slot - 1], rel=1e-9), key
    for key, start in [('wbt87', 87), ('tail_lift', 109),
                       ('vortex_lift', 129), ('ivbw', 67), ('go2pav', 45)]:
        _close(tail[key], wbt[start:start + n], key)


def test_probe_reaches_every_branch():
    results = [_run(p['inputs']) for p in _PROBE]
    tails = [t for _, t, _ in results if t is not None]
    assert {t['canard'] for t in tails} == {True, False}
    assert any(not p['inputs']['htpl'] for p in _PROBE)
    assert any(NOT_AVAILABLE in r['bwh']['cma'] for _, _, r in results
               if 'bwh' in r)
    assert any(-module.UNUSED in r['bwv']['cd'] for _, _, r in results)
    assert any(np.any(t['vortex_lift'] != 0.0) for t in tails)
    assert any(p['inputs']['synthesis']['alih'] != 0.0 for p in _PROBE)


def test_every_table_value_matches_the_source():
    source = parse('wgeotl')
    for name in ('X41602', 'X41601', 'X41603', 'Y44160'):
        np.testing.assert_array_equal(getattr(module, f'_{name}'),
                                      source[name], err_msg=name)
    # The carryover tables are shared with WBLIFT; pin them to WBTAIL's own.
    source = parse('wbtail')
    for name in ('X10A', 'Y10A', 'X10B', 'Y10B', 'X12A1', 'Y12A1', 'X12A2',
                 'Y12A2'):
        np.testing.assert_array_equal(getattr(module, f'_{name}'),
                                      source[name], err_msg=name)


def test_tail_vortex_term_uses_the_maximum_body_radius():
    """BODOWG overwrites WBT(108) with the maximum body radius before the
    vortex term reads it; the term scales with that radius, not SSPN-SSPNE."""
    probe = next(p for p in _PROBE if p['inputs']['htpl'])
    _, tail, _ = _run(probe['inputs'])
    inputs = probe['inputs']
    radius = np.sqrt(max(inputs['body']['s']) / np.pi)
    assert tail['body_radius'] == pytest.approx(radius)
    assert tail['body_radius'] != pytest.approx(
        inputs['tail']['sspn'] - inputs['tail']['sspne'])


def test_tail_is_stabilising():
    """The tail makes the wing-body-tail moment slope more negative."""
    probe = _PROBE[0]
    _, tail, result = _run(probe['inputs'])
    bw_cm = np.array(probe['inputs']['wing_body']['cm'])
    alpha = np.array(probe['inputs']['alpha_deg'])
    assert (np.polyfit(alpha[:4], result['bwh']['cm'][:4], 1)[0] <
            np.polyfit(alpha[:4], bw_cm[:4], 1)[0])


def test_first_slope_is_wbtail_s_analytic_value():
    """M10O12 recomputes CLa and CMa from the second angle only."""
    probe = _PROBE[0]
    _, tail, result = _run(probe['inputs'])
    assert result['bwh']['cla'][0] == tail['cla'][0]
    assert result['bwh']['cma'][0] == tail['cma'][0]
    assert result['bwh']['cla'][1] != tail['cla'][1]
