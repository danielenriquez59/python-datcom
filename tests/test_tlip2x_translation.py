"""
Regression tests for TLIP2X, TLINEX over a packed table.

Source of truth: datcom-legacy/datcom_2000/tlip2x.f.

Checked against a compiled probe (tools/probes/tlip2x.py) on interior,
grid-point and extrapolated points, for a whole table and for a slice of a
three-variable table selected through the shape array.
"""

import json
import pathlib

import pytest

from pydatcom.utils.packed_tables import tlip2x, unpack_table

_PROBE = json.loads((pathlib.Path(__file__).resolve().parent / 'fixtures' /
                     'probes' / 'tlip2x.json').read_text())


@pytest.mark.parametrize("case", range(len(_PROBE['cases'])))
def test_matches_compiled_routine(case):
    p = _PROBE
    q1, q2, l1, l2, u1, u2 = p['cases'][case]
    rec = p['records'][case]
    plain = tlip2x(p['x1'], p['x2'], p['packed'], q1, q2, l1, l2, u1, u2)
    sliced = tlip2x(p['x1'], p['x2'], p['packed3'], q1, q2, l1, l2, u1, u2,
                    shape=[-4, 3, 2, 0, 0, 2, 0])
    assert plain == pytest.approx(rec['PLAIN'][0], rel=1e-12, abs=1e-14)
    assert sliced == pytest.approx(rec['SLICE'][0], rel=1e-12, abs=1e-14)


def test_probe_tables_pack_as_intended():
    assert list(unpack_table(_PROBE['packed'], 12)) == \
        pytest.approx(_PROBE['y'])


# --- TLIP3X and INTEP3 (tools/probes/tlip3x.py) ----------------------------

from pydatcom.utils.packed_tables import intep3, tlip3x  # noqa: E402

_P3 = json.loads((pathlib.Path(__file__).resolve().parent / 'fixtures' /
                  'probes' / 'tlip3x.json').read_text())
_N_T3 = len(_P3['points']) * len(_P3['modes'])


@pytest.mark.parametrize("case", range(_N_T3))
def test_tlip3x_matches_compiled_routine(case):
    q = _P3['points'][case // len(_P3['modes'])]
    m = _P3['modes'][case % len(_P3['modes'])]
    value = tlip3x(_P3['x1'], _P3['x2'], _P3['x3'], _P3['packed3'], *q, *m)
    assert value == pytest.approx(_P3['records'][case]['T3'][0], rel=1e-12,
                                  abs=1e-14)


def _intep3_replay():
    out, n = [], _N_T3
    for n2d in (4, 1):
        state = {}
        for lam in _P3['lamdas']:
            charts = {k: {'x1': _P3['x1'], 'x2': _P3['x2'], 'packed': w,
                          'nxx': [3, 0, 0, 0, 0, 0, 0], 'nx2': 4}
                      for k, w in _P3['charts'].items()}
            charts['d']['nx2'] = n2d
            out.append((intep3(0.8, 0.3, lam, charts, state),
                        _P3['records'][n]['I3'][0]))
            n += 1
    return out


@pytest.mark.parametrize("case", range(2 * len(_P3['lamdas'])))
def test_intep3_matches_compiled_routine(case):
    ours, expected = _intep3_replay()[case]
    assert ours == pytest.approx(expected, rel=1e-12, abs=1e-14)
