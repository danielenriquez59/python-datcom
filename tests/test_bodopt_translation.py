"""
Regression tests for BODOPT and overlay M04O04: the asymmetric body.

Source of truth: datcom-legacy/datcom_2000/bodopt.f, m04o04.f.

Checked against a compiled probe of the overlay (tools/probes/m04o04.py),
and every table against a re-parse of the FORTRAN.
"""

import json
import pathlib
import sys

import numpy as np
import pytest

from pydatcom.aerodynamics import bodopt as module
from pydatcom.aerodynamics.bodopt import calculate_m04o04

_ROOT = pathlib.Path(__file__).resolve().parent.parent
sys.path.insert(0, str(_ROOT / 'tools'))
from fortran_data import parse  # noqa: E402

_PROBE = json.loads(
    (_ROOT / 'tests' / 'fixtures' / 'probes' / 'm04o04.json').read_text())


def _run(c):
    return calculate_m04o04(c['x'], c['s'], c['p'], c['r'], c['zu'],
                            c['zl'], c['alpha'], c['mach'], c['rn'],
                            c['sref'], c['cbar'], c['blref'], c['xcg'],
                            c['roughness'], c['kbody'])


@pytest.mark.parametrize("case", range(len(_PROBE)))
def test_matches_compiled_overlay(case):
    probe = _PROBE[case]
    c, out = probe['inputs'], probe['outputs']
    result = _run(c)
    for tag, key in [('CD', 'cd'), ('CL', 'cl'), ('CM', 'cm'), ('CN', 'cn'),
                     ('CA', 'ca'), ('CLA', 'cla'), ('CMA', 'cma'),
                     ('CYB', 'cyb'), ('CNB', 'cnb'), ('CLB', 'clb')]:
        np.testing.assert_allclose(result[key], out[tag], rtol=1e-9,
                                   atol=1e-13, err_msg=tag)
    opt = result['bodopt']
    alpha_zero, cm0, cd0, cd10 = out['BD']
    assert opt['alpha_zero_lift'] == pytest.approx(alpha_zero, rel=1e-12)
    assert opt['cd_zero_lift'] == pytest.approx(cd0, rel=1e-12)
    assert cd10 == pytest.approx(cd0)
    # With experimental data M04O04 replaces BD(62) with the curves' CM0.
    expected = result['cm0'] if c['kbody'] else opt['cm0']
    assert expected == pytest.approx(cm0, rel=1e-10, abs=1e-15)
    np.testing.assert_allclose(opt['cdl'], out['CDL'], rtol=1e-10,
                               atol=1e-15)
    # BODOPT also writes the wing's local angles B(23) onward.
    np.testing.assert_allclose(np.array(c['alpha']) + c['aliw'], out['BLOC'])


def test_probe_reaches_every_branch():
    inputs = [p['inputs'] for p in _PROBE]
    results = [_run(c)['bodopt'] for c in inputs]
    assert any(c['s'][-1] == 0.0 for c in inputs)     # default RATIO, FR
    assert any(c['s'][-1] > 0.0 for c in inputs)
    assert any(c['kbody'] for c in inputs)
    # The viscous force switches on inside the schedule.
    assert any(r['onset_angle'] < max(c['alpha'])
               and r['onset_angle'] > min(c['alpha'])
               for r, c in zip(results, inputs))


def test_every_table_value_matches_the_source():
    source = parse('bodopt')
    for name in ('XBA1', 'YBA1', 'XBA2', 'YBA2', 'X1BA3', 'X2BA3', 'YBA3',
                 'X27M', 'X27I'):
        np.testing.assert_array_equal(getattr(module, f'_{name}'),
                                      source[name], err_msg=name)


def test_suspect_kv_entries_are_preserved():
    """Figure BA-2 breaks its rise at 3.626 and 3.9; kept as chart data."""
    assert list(module._YBA2[4:8]) == [3.21, 3.626, 3.31, 3.9]


def test_camber_gives_a_nonzero_zero_lift_angle_and_symmetric_body_none():
    c = dict(_PROBE[0]['inputs'])
    cambered = _run(c)['bodopt']['alpha_zero_lift']
    half = (np.array(c['zu']) - np.array(c['zl'])) / 2.0
    symmetric = dict(c, zu=list(half), zl=list(-half))
    assert cambered != 0.0
    assert _run(symmetric)['bodopt']['alpha_zero_lift'] == pytest.approx(
        0.0, abs=1e-12)
