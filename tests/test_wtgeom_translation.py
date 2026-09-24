"""
Regression tests for the complete WTGEOM translation and SETUP1's sweep
records.

Source of truth: datcom-legacy/datcom_2000/wtgeom.f, setup1.f, m02o02.f.

Every one of the 195 A-block words is checked against a compiled probe
(tools/probes/wtgeom.py), for single- and two-panel surfaces and for a
repeated call on a block that already holds values.
"""

import json
import math
import pathlib

import numpy as np
import pytest

from pydatcom.geometry.wing import calculate_straight_exposed_geometry
from pydatcom.geometry.wtgeom import (
    calculate_wtgeom, sweep_records, vertical_panel_adjustments,
)
from pydatcom.utils.constants import UNUSED

_ROOT = pathlib.Path(__file__).resolve().parent.parent
_PROBE = json.loads(
    (_ROOT / 'tests' / 'fixtures' / 'probes' / 'wtgeom.json').read_text())


def _run(probe):
    inputs, out = probe['inputs'], probe['outputs']
    before = {k + 1: v for k, v in enumerate(out['BEFORE'])}
    records = sweep_records(inputs['savsi'], inputs['savso'], before)
    a_in = dict(before)
    a_in.update(records)
    a_in[174] = inputs['xovc']
    ain = {int(k): v for k, v in inputs['ain'].items()}
    return records, calculate_wtgeom(ain, a_in)


@pytest.mark.parametrize("case", range(len(_PROBE)))
def test_matches_compiled_routine(case):
    probe = _PROBE[case]
    records, result = _run(probe)
    np.testing.assert_allclose([records[106 + k] for k in range(12)],
                               probe['outputs']['REC'], rtol=1e-12,
                               atol=1e-14)
    np.testing.assert_allclose([result['a'][k] for k in range(1, 196)],
                               probe['outputs']['A'], rtol=1e-9, atol=1e-13)
    assert result['ain'][5] == pytest.approx(probe['outputs']['AIN5'][0])


def test_probe_reaches_every_branch():
    inputs = [p['inputs'] for p in _PROBE]
    results = [_run(p)[1] for p in _PROBE]
    assert {r['single_panel'] for r in results} == {True, False}
    assert {i['ain']['9'] for i in inputs} >= {0.0, 0.25, 0.5, 1.0}
    assert any(i['ain']['9'] == i['xovc'] for i in inputs)
    assert any(i['repeat'] for i in inputs)
    assert any(i['savsi'] < 0 for i in inputs)


def test_straight_subset_agrees_with_the_named_view():
    """calculate_straight_exposed_geometry is the named-key view of the
    single-panel subset; the two translations must agree."""
    state = {'wing_type': 1.0, 'wing_chrdr': 7.0, 'wing_chrdtp': 3.0,
             'wing_sspn': 15.0, 'wing_sspne': 13.0, 'wing_savsi': 25.0,
             'wing_chstat': 0.25}
    named = calculate_straight_exposed_geometry(state)
    a_in = sweep_records(25.0, 25.0)
    a_in[174] = 0.3
    a = calculate_wtgeom({1: 3.0, 2: 0.0, 3: 13.0, 4: 15.0, 5: 3.0, 6: 7.0,
                          9: 0.25, 66: 0.3}, a_in)['a']
    pairs = [('area', 3), ('aspect_ratio', 7), ('root_chord', 10),
             ('taper_ratio', 27), ('mac', 16), ('mac_span_location', 31),
             ('tan_le', 38), ('tan_c4', 44), ('area_theoretical', 4),
             ('aspect_ratio_theoretical', 120), ('mac_theoretical', 122),
             ('mac_c4_theoretical', 161), ('mac_le_theoretical', 195)]
    for key, index in pairs:
        assert named[key] == pytest.approx(a[index], rel=1e-8), key


def test_two_panel_exposed_sweep_loses_its_sign():
    """The exposed record of a forward-swept two-panel surface comes back
    positive through ANGLES(4); a single-panel one keeps its sign."""
    probe = next(p for p in _PROBE
                 if p['inputs']['savsi'] < 0 and p['inputs']['ain']['2'] > 0)
    _, result = _run(probe)
    a = result['a']
    assert a[62] < 0.0 and a[86] < 0.0           # both panels forward
    assert a[38] > 0.0                            # exposed LE: sign lost
    single = next(p for p in _PROBE
                  if p['inputs']['savsi'] < 0 and p['inputs']['ain']['2'] == 0)
    assert _run(single)[1]['a'][44] < 0.0


def test_sweep_stations_follow_the_taper_relation():
    """tan(sweep_x) = tan(sweep_ref) + 4(1-l)/(A(1+l)) * (x_ref - x) on
    each panel, for the LE, c/4, c/2 and TE records."""
    a_in = sweep_records(30.0, 30.0)
    a_in[174] = 0.3
    a = calculate_wtgeom({1: 2.0, 2: 0.0, 3: 12.0, 4: 14.0, 5: 2.0, 6: 8.0,
                          9: 0.25, 66: 0.3}, a_in)['a']
    factor = 4.0 * (1.0 - a[26]) / (a[5] * (1.0 + a[26]))
    for station, start in ((0.0, 58), (0.5, 70), (1.0, 76)):
        assert a[start + 4] == pytest.approx(
            math.tan(math.radians(30.0)) + factor * (0.25 - station))


def test_vertical_panel_halving_and_offset():
    a = {1: 10.0, 3: 20.0, 7: 4.0, 120: 5.0, 5: UNUSED, 130: 1.0,
         133: 2.0, 136: 3.0, 16: 9.0}
    out = vertical_panel_adjustments(a, 0.5)
    assert (out[1], out[3], out[7], out[120]) == (5.0, 10.0, 2.0, 2.5)
    assert out[5] == UNUSED and out[16] == 9.0
    assert (out[130], out[133], out[136]) == (1.5, 2.5, 3.5)
