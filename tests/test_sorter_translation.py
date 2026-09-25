"""
Regression tests for SORTER, the INTEGER table sort.

Source of truth: datcom-legacy/datcom_2000/sorter.f.

Checked against a compiled probe (tools/probes/sorter.py) on ``A``, ``B``,
the row and column orders and ``IFLAG``.
"""

import json
import pathlib

import pytest

from pydatcom.io.sorter import sorter

_PROBE = json.loads((pathlib.Path(__file__).resolve().parent / 'fixtures' /
                     'probes' / 'sorter.json').read_text())
_MARK = -7


def _flat(table):
    return [table[i][j] for j in range(len(table[0]))
            for i in range(len(table))]


@pytest.mark.parametrize("case", range(len(_PROBE)))
def test_matches_compiled_routine(case):
    c, o = _PROBE[case]['inputs'], _PROBE[case]['outputs']
    a, b, ir, ic, iflag = sorter(c['a'], c['irow'], c['icol'], _MARK)
    ours = {'A': _flat(a), 'B': _flat(b) if b else [_MARK] * len(_flat(a)),
            'IR': ir, 'IC': ic, 'FLAG': [iflag]}
    for key, values in ours.items():
        assert [int(v) for v in o[key]] == values, key


def test_sorts_rows_then_columns():
    a, b, ir, ic, iflag = sorter([[5, 3], [2, 8]], 1, 1)
    assert ir == [2, 1] and a == [[2, 8], [5, 3]]
    assert ic == [1, 2] and b == [[2, 8], [5, 3]] and iflag == 0


def test_overflowing_sentinel_repeats_the_first_choice():
    """With BMAX+BMAX wrapping negative, the chosen entry stays the
    minimum and is chosen every time."""
    ir = sorter([[2000000000], [1500000000], [-4]], 0, 1)[2]
    assert ir == [3, 3, 3]
