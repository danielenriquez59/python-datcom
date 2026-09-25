"""
Regression tests for DMPARY, PRCSID and SWRITE, the small printers.

Source of truth: datcom-legacy/datcom_2000/dmpary.f, prcsid.f, swrite.f.

Checked against the records the compiled routines print
(tools/probes/printers.py); SWRITE runs on AUXOUT's own formats.
"""

import json
import pathlib

import pytest

from pydatcom.io.printers import dmpary, prcsid, swrite

_PROBE = json.loads((pathlib.Path(__file__).resolve().parent / 'fixtures' /
                     'probes' / 'printers.json').read_text())
_RECORDS = _PROBE['records']
_N_DMP, _N_TITLE = len(_PROBE['dmpary']), len(_PROBE['titles'])


@pytest.mark.parametrize("case", range(_N_DMP))
def test_dmpary_matches_compiled_routine(case):
    c = _PROBE['dmpary'][case]
    assert dmpary(c['array'], c['name'], c['nlet']) == _RECORDS[case]


@pytest.mark.parametrize("case", range(_N_TITLE))
def test_prcsid_matches_compiled_routine(case):
    title = list(_PROBE['titles'][case].ljust(74))
    assert prcsid(True, title) == _RECORDS[_N_DMP + case]


def test_prcsid_is_silent_without_head():
    assert _RECORDS[_N_DMP + _N_TITLE] == []
    assert prcsid(False, list('TITLE'.ljust(74))) == []


@pytest.mark.parametrize("case", range(len(_PROBE['swrite'])))
def test_swrite_matches_compiled_routine(case):
    c = _PROBE['swrite'][case]
    lines, ndmf, naf = swrite(c['last'], c['form'], c['icont'],
                              c['columns'], len(c['columns'][0]))
    flags = '@@FLAGS' + (' T' if ndmf else ' F') + (' T' if naf else ' F')
    assert lines + [flags] == _RECORDS[_N_DMP + _N_TITLE + 1 + case]


# --- MESSGE, in the source's word size (tools/probes/messge.py) -----------

from pydatcom.io.printers import messge  # noqa: E402

_MESSGE = json.loads((pathlib.Path(__file__).resolve().parent / 'fixtures' /
                      'probes' / 'messge.json').read_text())


@pytest.mark.parametrize("case", range(len(_MESSGE['cases'])))
def test_messge_matches_compiled_routine(case):
    c = _MESSGE['cases'][case]
    m = list(_MESSGE['msscl'][case])
    m[1], m[2] = 'MSG1', 'MSG2'
    lines = messge(['TLIN', 'EX  '], c['mess'], _MESSGE['x'], m, c['nf'],
                   c['iovly'])
    assert [ln.rstrip() for ln in lines] == _MESSGE['records'][case]
