"""
Regression tests for overlay M06O06: BODYRT, BODYJM, and the body slope and
lateral pass.

Source of truth: datcom-legacy/datcom_2000/m06o06.f, bodyrt.f, bodyjm.f.

Checked against a compiled probe of the overlay (tools/probes/m06o06.py),
which was also BODYRT's first check against execution and found its
crossflow planform integral starting at the wrong station.
"""

import json
import pathlib

import numpy as np
import pytest

from pydatcom.aerodynamics.bodyrt import calculate_bodyrt, calculate_m06o06
from pydatcom.utils.constants import UNUSED

_ROOT = pathlib.Path(__file__).resolve().parent.parent
_PROBE = json.loads(
    (_ROOT / 'tests' / 'fixtures' / 'probes' / 'm06o06.json').read_text())

_CURVES = [('CD', 'cd'), ('CL', 'cl'), ('CM', 'cm'), ('CN', 'cn'),
           ('CA', 'ca'), ('CLA', 'cla'), ('CMA', 'cma'), ('CYB', 'cyb'),
           ('CNB', 'cnb'), ('CLB', 'clb')]


def _run(c):
    return calculate_m06o06(c['x'], c['s'], c['p'], c['r'], c['alpha'],
                            c['alpha_zero'], c['mach'], c['rn'], c['sref'],
                            c['cbar'], c['blref'], c['bd33'], c['xcg'],
                            c['roughness'], c['method'], c['ellip'],
                            c['kbody'])


@pytest.mark.parametrize("case", range(len(_PROBE)))
def test_matches_compiled_overlay(case):
    probe = _PROBE[case]
    c, out = probe['inputs'], probe['outputs']
    result = _run(c)
    for tag, key in _CURVES:
        np.testing.assert_allclose(result[key], out[tag], rtol=1e-9,
                                   atol=1e-13, err_msg=tag)
    rt = result['bodyrt']
    base, cd0, cm0, eta, volume, centroid = out['BD']
    assert rt['base_area'] == pytest.approx(base, rel=1e-12)
    assert rt['cd_zero_lift'] == pytest.approx(cd0, rel=1e-12)
    assert rt['bd76_eta'] == pytest.approx(eta, rel=1e-12)
    if c['kbody']:
        assert result['cm0'] == pytest.approx(cm0, abs=1e-12)
    jm = result['bodyjm']
    if jm is not None:
        assert jm['volume'] == pytest.approx(volume, rel=1e-12)
        assert jm['centroid'] == pytest.approx(centroid, rel=1e-12)
        assert jm['planform_area'] == pytest.approx(out['MORE'][0], rel=1e-12)
        assert jm['ellipticity'] == pytest.approx(out['MORE'][1])


def test_probe_reaches_every_branch():
    inputs = [p['inputs'] for p in _PROBE]
    assert {c['method'] > 1.5 for c in inputs} == {True, False}
    ellip = {c['ellip'] for c in inputs if c['method'] > 1.5}
    assert UNUSED in ellip and 1.0 in ellip
    assert any(e < 1.0 and e != UNUSED for e in ellip)
    assert any(e > 1.0 for e in ellip)
    assert any(c['kbody'] for c in inputs)
    assert any(c['alpha_zero'] != 0.0 for c in inputs)


def test_crossflow_planform_starts_at_the_substituted_station():
    """BODYRT integrates the planform from TMP5, the substituted station;
    an earlier translation started at the original one, which changed the
    viscous force by several percent on a body with no station at TMP5."""
    probe = next(p for p in _PROBE
                 if p['inputs']['r'][0] == 0.0 and p['inputs']['mach'] == 0.7)
    c = probe['inputs']
    rt = calculate_bodyrt(c['x'], c['s'], c['p'], c['r'],
                          np.array(c['alpha']) + c['alpha_zero'], c['mach'],
                          c['rn'], c['sref'], c['cbar'], c['bd33'],
                          c['roughness'])
    np.testing.assert_allclose(rt['cn'], probe['outputs']['CL'], rtol=1e-9,
                               atol=1e-13)


def test_body_is_directionally_unstable_and_rolls_neutrally():
    """CY_b = -CLa, Cn_b = -(c/b)*CMa, Cl_b = 0 for the axisymmetric body."""
    result = _run(_PROBE[0]['inputs'])
    c = _PROBE[0]['inputs']
    np.testing.assert_array_equal(result['cyb'], -result['cla'])
    np.testing.assert_allclose(result['cnb'],
                               -(c['cbar'] / c['blref']) * result['cma'])
    assert np.all(result['clb'] == 0.0)
    assert np.all(result['cnb'] < 0.0)


def test_normal_force_is_odd_and_drag_even():
    """For an uncambered body at zero incidence."""
    probe = next(p for p in _PROBE if p['inputs']['alpha_zero'] == 0.0)
    result = _run(probe['inputs'])
    alpha = np.array(probe['inputs']['alpha'])
    for a in (2.0,):
        i, k = list(alpha).index(-a), list(alpha).index(a)
        assert result['cl'][i] == pytest.approx(-result['cl'][k])
        assert result['cd'][i] == pytest.approx(result['cd'][k])
