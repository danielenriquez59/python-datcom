"""
Regression tests for CTABS, the control-tab hinge moments and control
force, and the overlay M36O44.

Source of truth: datcom-legacy/datcom_2000/ctabs.f and m36o44.f.

CTABS is checked against a compiled probe (tools/probes/ctabs.py) on
``BW``, ``BH``, ``BV``, ``BWH`` and ``BWHV`` words 201-380.
"""

import json
import pathlib

import numpy as np
import pytest

from pydatcom.aerodynamics.control_tab import calculate_ctabs, m36o44
from pydatcom.utils.constants import UNUSED

_PROBE = json.loads((pathlib.Path(__file__).resolve().parent / 'fixtures' /
                     'probes' / 'ctabs.json').read_text())


@pytest.mark.parametrize("case", range(len(_PROBE)))
def test_matches_compiled_routine(case):
    r = calculate_ctabs(_PROBE[case]['inputs'])
    for name, values in _PROBE[case]['outputs'].items():
        ours = [r[name.lower()][201 + i] for i in range(len(values))]
        np.testing.assert_allclose(ours, values, rtol=1e-9, atol=1e-14,
                                   err_msg=name)


def test_probe_reaches_every_branch():
    c = [p['inputs'] for p in _PROBE]
    assert {p['f']['117'] for p in c} == {1.0, 2.0, 3.0}
    rl = {p['f']['135'] for p in c}
    assert min(rl) < 0.0 and 0.0 in rl and max(rl) > 0.0
    assert any(p['pinf'] == UNUSED for p in c)


def test_force_coefficient_overwrites_the_force():
    """CFC and FC share BW 201: with a known pressure only the
    coefficient survives."""
    c = _PROBE[0]['inputs']
    r = calculate_ctabs(c)
    no_q = calculate_ctabs(dict(c, pinf=UNUSED))
    q = 0.7 * c['pinf'] * c['mach'] ** 2
    assert r['bw'][201] == pytest.approx(
        no_q['bw'][201] / (q * c['sref'] * c['cbarr']))


def test_m36o44_sequence():
    calls = []

    def run(name):
        return lambda: calls.append(name) or name

    ran, wing = m36o44(False, 1.0, True, run('LIFTFP'), run('HINGE'),
                       run('CTABS'))
    assert calls == ['LIFTFP', 'HINGE', 'CTABS']
    assert wing == {250 + j: -UNUSED for j in range(2, 11)}
    calls.clear()
    ran, wing = m36o44(False, 2.0, False, run('LIFTFP'), run('HINGE'),
                       run('CTABS'))
    assert calls == ['LIFTFP'] and list(ran) == ['LIFTFP']
    calls.clear()
    assert m36o44(True, 1.0, True, run('LIFTFP'), run('HINGE'),
                  run('CTABS'))[1] == {}
    assert calls == ['LIFTFP']
