"""
Regression tests for NAMEW (namelist printing), TOINT (integer decoding)
and the dump overlays M11O13 and M42O52.

Source of truth: datcom-legacy/datcom_2000/namew.f, toint.f, m11o13.f,
m42o52.f.

NAMEW and TOINT are checked against records the compiled routines print
(tools/probes/namelist_write.py).  NAMEW runs in the source's own word
size, where a LOGICAL word printed through FORREA shows its bits as a
REAL.
"""

import json
import pathlib

import pytest

from pydatcom.io.dump_overlays import m11o13, m42o52
from pydatcom.io.namelist_io import namew, toint
from pydatcom.io.printers import dmpary
from pydatcom.utils.constants import UNUSED

_PROBE = json.loads((pathlib.Path(__file__).resolve().parent / 'fixtures' /
                     'probes' / 'namelist_write.json').read_text())


def _block():
    cb = [v if v is not None else 0.0 for v in _PROBE['reals']]
    for k, v in _PROBE['logicals'].items():
        cb[int(k) - 1] = v
    return cb


@pytest.mark.parametrize("case", range(len(_PROBE['lists'])))
def test_namew_matches_compiled_routine(case):
    c = _PROBE['lists'][case]
    lines, _ = namew('$', list('FLTCON'), list(''.join(_PROBE['names'])),
                     _PROBE['lenvn'], c['vdime'], _block(), c['loc'])
    assert lines == _PROBE['namew'][case]


@pytest.mark.parametrize("case", range(len(_PROBE['ints'])))
def test_toint_matches_compiled_routine(case):
    ians, ierr = toint(list(_PROBE['ints'][case].ljust(80)[:80]))
    assert [f'{ians:12d}{ierr:12d}'] == _PROBE['toint'][case]


def test_namew_zero_length_prints_the_name_alone():
    """A zero length prints no values and keeps the saved VTYPE."""
    lines, vtype = namew('$', list('X'), list('AB'), [1, 1], [-1, 0],
                         [True, False], [1, 2])
    assert vtype == 0 and lines[2] == '0B      ='


def test_m11o13_starts_from_unused_and_dumps():
    seen = {}

    def grdeff(gr):
        seen['gr'] = list(gr)
        return [1.0] * 303

    gr, lines = m11o13(grdeff, dump=True)
    assert seen['gr'] == [UNUSED] * 303
    assert lines == dmpary([1.0] * 303, 'GR', 2)
    assert m11o13(grdeff, dump=False)[1] == []


def test_m42o52_dumps_six_blocks_in_order():
    blocks = {k: [float(n)] for n, k in
              enumerate(('f', 'body', 'ht', 'hyp', 'vt', 'wing'))}
    _, lines = m42o52(lambda: blocks, dump=True)
    names = [ln[2:6] for ln in lines if '(' in ln]
    assert names == ['   F', 'BODY', '  HT', ' HYP', '  VT', 'WING']
