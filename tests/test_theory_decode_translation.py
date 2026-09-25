"""
Regression tests for TOLOG, REPTCT (namelist value decoding), THEORY (the
airfoil-section report) and M51O63.

Source of truth: datcom-legacy/datcom_2000/tolog.f, reptct.f, theory.f,
m51o63.f.

TOLOG, REPTCT and THEORY are checked against compiled probes
(tools/probes/theory_decode.py); THEORY runs with IDEAL and SLOPE stubbed,
and with its CALL EXIT redirected, since gfortran binds that call to its
own intrinsic EXIT, which ends the run.
"""

import json
import pathlib

import pytest

from pydatcom.geometry.airfoil_report import theory
from pydatcom.io.block_setup import m51o63
from pydatcom.io.fortran_format import fortran_write
from pydatcom.io.namelist_io import reptct, tolog
from pydatcom.utils.constants import UNUSED

_PROBE = json.loads((pathlib.Path(__file__).resolve().parent / 'fixtures' /
                     'probes' / 'theory_decode.json').read_text())
_N_LOG = len(_PROBE['logs'])


@pytest.mark.parametrize("case", range(_N_LOG))
def test_tolog_matches_compiled_routine(case):
    card = list(_PROBE['logs'][case].ljust(80))
    assert fortran_write('(L2,I3)', tolog(card)) == _PROBE['decode'][case]


@pytest.mark.parametrize("case", range(len(_PROBE['repts'])))
def test_reptct_matches_compiled_routine(case):
    kol = list(_PROBE['repts'][case].ljust(80))
    irept, ierr = reptct(kol)
    assert fortran_write('(2I6,1X,80A1)', [irept, ierr] + kol) == \
        _PROBE['decode'][_N_LOG + case]


def _slope(mach, renn):
    return (UNUSED if mach > 0.7 else 0.1 + 0.02 * mach), 0.25 + 0.01 * mach


@pytest.mark.parametrize("case", range(len(_PROBE['theory'])))
def test_theory_matches_compiled_routine(case):
    c, rec = _PROBE['theory'][case], _PROBE['records'][case]
    s = dict(_PROBE['section'], surface=c['surface'], mcc=c['mcc'],
             cla=c['cla'], xac=[-1.0] * len(c['machs']))
    r = theory(s, c['machs'], [1e6 * (i + 1) for i in range(len(c['machs']))],
               4.0, lambda: None, _slope)
    k = rec.index('@@VALUES')
    assert r['lines'] == [ln for ln in rec[:k] if ln != '@@EXIT']
    assert r['exit'] == ('@@EXIT' in rec[:k])
    values = [float(v) for v in rec[k + 1].split()]
    assert [r['cla0']] + r['cla'] + r['xac'] == pytest.approx(values,
                                                              rel=1e-12)


def test_m51o63_dispatch():
    powr = [0.0] * 316
    assert m51o63(1, lambda: None, lambda: None, powr) == ['INITZ1']
    assert m51o63(2, lambda: None, lambda: None, powr) == ['INITZ2']
    assert m51o63(3, lambda: None, lambda: None, powr) == []
    assert set(powr[1:]) == {UNUSED}
