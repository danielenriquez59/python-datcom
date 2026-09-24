"""
Regression tests for WBCD, WBCDL and TABLES: the wing-body and tail-body
regression drag, and M08O10's vertical-panel output setup.

Source of truth: datcom-legacy/datcom_2000/wbcd.f, wbcdl.f, tables.f,
m08o10.f.

TABLES is checked directly against a compiled probe at 232 Mach and angle
points, and WBCD for both halves (tools/probes/wbcd.py).
"""

import json
import pathlib

import numpy as np
import pytest

from pydatcom.aerodynamics.vertical_panel import m08o10_panel_block
from pydatcom.aerodynamics.wing_body import (
    UNUSED, calculate_tables, calculate_wbcd, calculate_wbcdl,
)

_ROOT = pathlib.Path(__file__).resolve().parent.parent
_PROBE = json.loads(
    (_ROOT / 'tests' / 'fixtures' / 'probes' / 'wbcd.json').read_text())


@pytest.mark.parametrize("point", range(len(_PROBE['tables'])))
def test_tables_matches_compiled_routine(point):
    t = _PROBE['tables'][point]
    b = calculate_tables(t['mach'], t['alpha'])
    assert (b is None) == t['ndm']
    if b is not None:
        np.testing.assert_allclose(b, t['b'], rtol=1e-10, atol=1e-15)


def _half(inputs, key):
    comb, arr, x = (('wing_body', 'BW', inputs['xw']) if key == 'wing' else
                    ('tail_body', 'BH', inputs['xh']))
    flight = {'mach': inputs['mach'], 'reynolds_per_length': inputs['rn'],
              'tr': inputs['tr']}
    return arr, calculate_wbcd(inputs['alpha'], inputs[key], inputs[comb],
                               inputs['body'], x, flight)


@pytest.mark.parametrize("case", range(len(_PROBE['wbcd'])))
def test_wbcd_matches_compiled_routine(case):
    probe = _PROBE['wbcd'][case]
    for key in ('wing', 'tail'):
        arr, result = _half(probe['inputs'], key)
        if result is None:
            # The source leaves the buildup drag untouched: zero here.
            for tag in ('CD', 'CN', 'CA'):
                assert np.all(np.array(probe['outputs'][arr + tag]) == 0.0)
            continue
        for tag, name in [('CD', 'cd'), ('CN', 'cn'), ('CA', 'ca')]:
            np.testing.assert_allclose(result[name],
                                       probe['outputs'][arr + tag],
                                       rtol=1e-10, atol=1e-15)


def test_probe_reaches_every_branch():
    ndm = [t['ndm'] for t in _PROBE['tables']]
    assert any(ndm) and not all(ndm)
    machs = {t['mach'] for t in _PROBE['tables']}
    assert {0.9, 1.0, 1.1, 2.5} <= machs and any(m > 2.5 for m in machs)
    results = [(_half(p['inputs'], 'wing')[1], _half(p['inputs'], 'tail')[1])
               for p in _PROBE['wbcd']]
    assert any(w is None for w, _ in results)
    assert any(t is not None for _, t in results)
    assert any(w is not None and UNUSED in w['cdl'] for w, _ in results)


def test_angle_limits_by_speed_range():
    """18, 11, 12 and 15 degrees below 0.9, to 1.0, to 1.1, and above."""
    for mach, limit in ((0.5, 18.0), (0.95, 11.0), (1.05, 12.0), (1.5, 15.0)):
        assert calculate_tables(mach, limit) is not None
        assert calculate_tables(mach, limit + 0.01) is None
    assert calculate_tables(2.6, 2.0) is None


def test_tables_is_continuous_across_the_integer_angles():
    for alpha in (3.0, 7.0, 10.0):
        below = calculate_tables(0.5, alpha - 1e-9)
        above = calculate_tables(0.5, alpha + 1e-9)
        np.testing.assert_allclose(below, above, atol=1e-7)


def test_wbcdl_twist_is_checked_in_degrees_and_applied_in_radians():
    args = dict(aspect_ratio=4.0, tan_le=0.5, tovc=0.06, nose_length=4.0,
                afterbody_length=4.0, taper_ratio=0.4,
                leading_edge_radius=0.008, ycm=0.01, cld=0.2,
                reynolds=5e6, tr=0.4, mach=0.5, alpha_deg=[4.0])
    assert calculate_wbcdl(twist_deg=-9.0, **args) is not None
    assert calculate_wbcdl(twist_deg=-9.5, **args) is None
    b = calculate_tables(0.5, 4.0)
    diff = (calculate_wbcdl(twist_deg=-2.0, **args)[0] -
            calculate_wbcdl(twist_deg=0.0, **args)[0])
    assert diff == pytest.approx(b[12] * 2.0 / 57.2957795)


def test_m08o10_panel_block():
    block = m08o10_panel_block(0.004, 4)
    assert block['cd'][0] == 0.004 and block['cl'][0] == 0.0
    assert np.all(block['cd'][1:] == -1.0e-30)
