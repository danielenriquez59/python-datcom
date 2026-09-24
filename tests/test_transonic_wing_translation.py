"""
Regression tests for TRANWG and CLMXB1, TRSONI's two anchors.

Source of truth: datcom-legacy/datcom_2000/tranwg.f, clmxb1.f.

Checked against a compiled probe (tools/probes/tranwg.py); the shared
tables are pinned to each routine's own DATA.
"""

import json
import math
import pathlib
import sys

import numpy as np
import pytest

from pydatcom.aerodynamics.transonic_wing import (
    calculate_clmxb1, calculate_tranwg,
)

_ROOT = pathlib.Path(__file__).resolve().parent.parent
sys.path.insert(0, str(_ROOT / 'tools'))
from fortran_data import parse  # noqa: E402

_PROBE = json.loads(
    (_ROOT / 'tests' / 'fixtures' / 'probes' / 'tranwg.json').read_text())


def _run(probe):
    i = probe['inputs']
    a = {int(k): v for k, v in i['a'].items()}
    m = i['clmxb1']
    return (calculate_tranwg(a, i['deltay'], i['sref']),
            calculate_clmxb1(m['bu4'], m['a160'], i['deltay'], m['xovc'],
                             a[3], i['sref']))


@pytest.mark.parametrize("case", range(len(_PROBE)))
def test_matches_compiled_routines(case):
    probe = _PROBE[case]
    tranwg, clmxb1 = _run(probe)
    cna, dcna, a62, a86 = probe['outputs']['OUT']
    assert tranwg['cna'] == pytest.approx(cna, rel=1e-12)
    assert tranwg['dcna'] == pytest.approx(dcna, rel=1e-10, abs=1e-15)
    assert (tranwg['a62'], tranwg['a86']) == (a62, a86)
    assert clmxb1['clmax'] == pytest.approx(probe['outputs']['CLS'][0],
                                            rel=1e-12)


def test_probe_reaches_every_branch():
    branches = {b for p in _PROBE for b in _run(p)[0]['branches']}
    assert branches == {'56a', '56g', 'rectangular'}
    regimes = set()
    for p in _PROBE:
        a = p['inputs']['a']
        for mach in (1.3, 1.4, 1.5):
            regimes.add(math.sqrt(mach**2 - 1) / (a['62'] or 1e-5) > 1.0)
    assert regimes == {True, False}


def test_shared_tables_are_the_routines_own():
    """TRANWG's figures are VTLIFT's and CLMXB1's are CLMXBS's."""
    tranwg, vtlift = parse('tranwg'), parse('vtlift')
    for name in ('A1350', 'DA50', 'B1350', 'DB50', 'G13246', 'DG3246',
                 'T13246', 'DUMY1', 'DUMY2', 'DUMY3', 'DUMY4', 'DUMY5',
                 'DUMY6'):
        assert tranwg[name] == vtlift[name], name
    clmxb1, clmxbs = parse('clmxb1'), parse('clmxbs')
    for name in ('C1ABC', 'DYAG', 'CBASE', 'C2A', 'AMN', 'DE'):
        assert clmxb1[name] == clmxbs[name], name


def test_slope_derivative_matches_a_finite_difference_of_the_fit():
    """DCNA is the Mach derivative at 1.4 of the quadratic in beta through
    the three anchors, CN*beta**2 = AA*beta**2 + B*beta + C."""
    tranwg, _ = _run(_PROBE[0])
    betas = [math.sqrt(m * m - 1.0) for m in (1.3, 1.4, 1.5)]
    coeffs = np.polyfit(betas, [c * b * b for c, b in
                                zip(tranwg['anchor_slopes'], betas)], 2)

    def cn(mach):
        beta = math.sqrt(mach * mach - 1.0)
        return np.polyval(coeffs, beta) / beta**2

    h = 1e-6
    assert tranwg['dcna'] == pytest.approx((cn(1.4 + h) - cn(1.4 - h)) / (2 * h),
                                           rel=1e-5)


# --------------------------------------------------------------------------
# TRSONI, run over a Mach sequence as the main loop runs it
# --------------------------------------------------------------------------

from pydatcom.aerodynamics import transonic_wing as module  # noqa: E402
from pydatcom.aerodynamics.transonic_wing import calculate_trsoni  # noqa: E402

_TRSONI = json.loads(
    (_ROOT / 'tests' / 'fixtures' / 'probes' / 'trsoni.json').read_text())


def _sequence(config):
    """Every Mach of a configuration, the stale TRA words carried over."""
    i = config['inputs']
    a = {int(k): v for k, v in i['a'].items()}
    stale, out = None, []
    for run in config['runs']:
        r = calculate_trsoni(run['mach'], i['alpha'], i['wing'], a,
                             i['sref'], i['roughness'], i['body'], stale,
                             i['nf'])
        out.append((run, r))
        stale = [run['outputs']['TRA'][42 + k] for k in range(15)]
    return out


@pytest.mark.parametrize("config", range(len(_TRSONI)))
def test_trsoni_matches_compiled_routine(config):
    i = _TRSONI[config]['inputs']
    for run, r in _sequence(_TRSONI[config]):
        out = run['outputs']
        for index, value in r['tra'].items():
            assert value == pytest.approx(out['TRA'][index - 1], rel=1e-9,
                                          abs=1e-14), f"TRA({index})"
        cla, cd0, a160, a62, a86 = out['W']
        if 'cla' in r:
            assert r['cla'] == pytest.approx(cla, rel=1e-10)
            assert r['cd0'] == pytest.approx(cd0, rel=1e-10, abs=1e-15)
            assert (r['a62'], r['a86']) == pytest.approx((a62, a86))
        if 'a160' in r:
            assert r['a160'] == pytest.approx(a160, rel=1e-12)
        if i['body'] and not r.get('early_return'):
            body_cla, body_cma, cd0_wb, db = out['B']
            assert r['body_cla'] == pytest.approx(body_cla, rel=1e-12)
            assert r['body_cma'] == pytest.approx(body_cma, rel=1e-12)
            assert r['cd0_wing_body'] == pytest.approx(cd0_wb, rel=1e-10)
            assert r['base_diameter'] == pytest.approx(db, rel=1e-12)
            np.testing.assert_allclose(r['body_cd'], out['CDJ'], rtol=1e-10)


def test_trsoni_probe_reaches_every_branch():
    results = [r for c in _TRSONI for _, r in _sequence(c)]
    assert any(r.get('early_return') for r in results)
    assert any('cla' not in r and not r.get('early_return') for r in results)
    lowar = [r for r in results if r.get('low_aspect_ratio')]
    assert any(r['a160'] <= 4.5 for r in lowar)
    assert any(r['a160'] > 4.5 for r in lowar)
    assert any(r.get('low_aspect_ratio') is False for r in results)
    assert any('body_cd' in r for r in results)
    machs = {run['mach'] for c in _TRSONI for run in c['runs']}
    assert min(machs) < 1.0 < 1.2 < max(machs)


def test_trsoni_tables_match_the_source():
    source = parse('trsoni')
    for name in ('X', 'Y', 'TR', 'DR', 'X27M', 'X27I', 'T43A', 'D43A',
                 'T43B', 'D43B', 'T44A', 'D44A', 'T44B', 'D44B', 'T44C',
                 'D44C', 'T418A', 'D418A', 'T18B1', 'D18B1', 'C18B1',
                 'T18B2', 'D18B2', 'C18B2', 'T419A', 'D419A', 'T419B',
                 'D419B', 'T419C', 'D419C', 'T429L', 'D429L', 'T429R',
                 'D429R', 'T424', 'D424', 'T426', 'D426'):
        np.testing.assert_array_equal(getattr(module, f'_{name}'),
                                      source[name], err_msg=name)


def test_supersonic_wave_drag_points_are_never_stored():
    """The four points above normal Mach one keep their stale (zero) values,
    so the faired wave drag falls back to zero past about Mach 1.1."""
    i = _TRSONI[0]['inputs']
    a = {int(k): v for k, v in i['a'].items()}
    r = calculate_trsoni(1.3, i['alpha'], i['wing'], a, i['sref'],
                         i['roughness'], None)
    assert r['stale_wave_points'] == [11, 12, 13, 14]
    assert [r['tra'][43 + k] for k in (11, 12, 13, 14)] == [0.0] * 4
    assert r['wave_drag'] == 0.0
    peak = calculate_trsoni(1.05, i['alpha'], i['wing'], a, i['sref'],
                            i['roughness'], None)['wave_drag']
    assert peak > 0.01


def test_lift_slope_meets_its_anchors():
    """TRANF passes through the subsonic slope below the force break and
    TRANWG's slope at Mach 1.4."""
    i = _TRSONI[0]['inputs']
    a = {int(k): v for k, v in i['a'].items()}
    at14 = calculate_trsoni(1.4, i['alpha'], i['wing'], a, i['sref'],
                            i['roughness'], None)
    assert at14['cla'] == pytest.approx(at14['tra'][1], rel=1e-12)
    low = calculate_trsoni(0.5, i['alpha'], i['wing'], a, i['sref'],
                           i['roughness'], None)
    assert low['tra'][16] == 0.5                  # MT(1) = XM = MACH
    assert low['cla'] == pytest.approx(low['tra'][21], rel=1e-12)
