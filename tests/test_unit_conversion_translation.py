"""
Regression tests for CONV, the input unit conversion and scaling.

Source of truth: datcom-legacy/datcom_2000/conv.f.

Checked against a compiled probe (tools/probes/conv.py) on every word of
the arrays and scalars CONV converts, the body-origin shift, and the line
it prints.
"""

import json
import pathlib

import numpy as np
import pytest

from pydatcom.io.unit_conversion import calculate_conv

_PROBE = json.loads((pathlib.Path(__file__).resolve().parent / 'fixtures' /
                     'probes' / 'conv.json').read_text())


@pytest.mark.parametrize("case", range(len(_PROBE)))
def test_matches_compiled_routine(case):
    c, o = _PROBE[case]['inputs'], _PROBE[case]['outputs']
    r = calculate_conv(c)
    for name, values in o.items():
        if name == '_text':
            assert values == ['0 ' + r['message']]
            continue
        if name == 'SC':
            ours = [r[k] for k in ('wt', 'rougfc', 'sref', 'cbarr', 'blref')]
        elif name == 'BD':
            ours = [r['bd'][k] for k in (11, 33, 65, 74, 82)]
        else:
            ours = r[name.lower()][1:]
        np.testing.assert_allclose(ours, values, rtol=1e-9, atol=1e-14,
                                   err_msg=name)


def test_probe_reaches_every_unit_and_flag():
    c = [p['inputs'] for p in _PROBE]
    assert {p['idim'] for p in c} == {1, 2, 3, 4}
    for flag in ('symfp', 'asyfp', 'trajet', 'hypef'):
        assert any(p[flag] for p in c)
    assert any(p['x0'] == 0.0 for p in c)
