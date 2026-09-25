"""
Regression tests for the namelist reader internals: FINDCH, SKIPBL,
EXTRST, TODEC, FINDVN, FORINT, FORLOG, FORREA, SUBINT/SUBLOG/SUBREA and
READCD.

Source of truth: the matching files in datcom-legacy/datcom_2000.

Checked against the records the compiled routines print
(tools/probes/namelist_io.py).  READCD is checked on the cards before end
of file: under gfortran its retry after END stops the run with "Read past
ENDFILE record", where the source expects END again.
"""

import io
import json
import pathlib

import pytest

from pydatcom.io import namelist_io as m
from pydatcom.io.fortran_format import fortran_write

_PROBE = json.loads((pathlib.Path(__file__).resolve().parent / 'fixtures' /
                     'probes' / 'namelist_io.json').read_text())


def _card(text):
    return list(text.ljust(80)[:80])


def _ours():
    p, out = _PROBE, []
    for text, ch, k in p['finds']:
        out.append(fortran_write('(I4)', [m.findch(_card(text), ch, k)]))
    for text, k in p['blanks']:
        out.append(fortran_write('(I4)', [m.skipbl(_card(text), k)]))
    out.append([''.join(m.extrst(_card('ABCDEFGHIJ'), 3, 7))])
    for text in p['numbers']:
        out.append(fortran_write('(1P,E25.16,I3)', m.todec(_card(text))))
    vname = list(''.join(p['names']))
    lens = [len(n) for n in p['names']]
    for text in p['lookups']:
        out.append(fortran_write('(I4,L2)',
                                 m.findvn(lens, _card(text), vname)))
    name = list('VARNAME ')
    out += [m.forint(name, v) for v in p['ints']]
    out += [m.forlog(name, v) for v in p['logs']]
    out += [m.forrea(name, v) for v in p['reals']]
    unit = io.StringIO('\n'.join(p['cards']) + '\n')
    kol = _card('PREVIOUS')
    for _ in p['cards']:
        kol, eof = m.readcd(unit, kol)
        out.append(fortran_write('(80A1,L2)', kol + [eof]))
    return out


_OURS = _ours()


@pytest.mark.parametrize("case", range(len(_PROBE['records'])))
def test_matches_compiled_routines(case):
    assert _OURS[case] == _PROBE['records'][case]


def test_readcd_reports_end_of_file():
    kol, eof = m.readcd(io.StringIO(''), _card('KEEP'))
    assert eof and ''.join(kol).strip() == 'KEEP'


def test_todec_digit_overflow_wraps():
    """The mantissa digits accumulate in a 32-bit INTEGER."""
    assert m.todec(_card('12345678901'))[0] == -539222987.0


def test_assignment_helpers():
    assert m.subint(3) == 3 and m.sublog(True) and m.subrea(2.5) == 2.5
