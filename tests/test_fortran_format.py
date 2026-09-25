"""
Tests for the FORTRAN FORMAT interpreter, against records written by
gfortran (tools/probes/fortran_format.py).
"""

import json
import pathlib

import pytest

from pydatcom.io.fortran_format import fortran_write

_PROBE = json.loads((pathlib.Path(__file__).resolve().parent / 'fixtures' /
                     'probes' / 'fortran_format.json').read_text())


@pytest.mark.parametrize("case", range(len(_PROBE)))
def test_matches_gfortran(case):
    c = _PROBE[case]
    assert fortran_write(c['format'], c['items']) == c['records']


def test_literals_after_the_last_item_are_written():
    assert fortran_write('(1X,I2,3H ok,I2)', [5]) == ['  5 ok']


def test_scale_factor_shifts_decimal_digits_exactly():
    assert fortran_write('(1P,F10.3)', [9.99995]) == ['   100.000']
