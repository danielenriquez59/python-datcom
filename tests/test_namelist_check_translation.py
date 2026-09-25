"""
Regression tests for VNAME, LVALUE and RVALUE, the namelist card syntax
checks.

Source of truth: datcom-legacy/datcom_2000/vname.f, lvalue.f, rvalue.f.

Checked against a compiled probe (tools/probes/namelist_check.py) on the
column reached, the value count and the fault count (and for VNAME the
name matched, its subscript and whether it was subscripted), on 30 cards.
"""

import json
import pathlib

import pytest

from pydatcom.io.namelist_check import lvalue, rvalue, vname

_PROBE = json.loads((pathlib.Path(__file__).resolve().parent / 'fixtures' /
                     'probes' / 'namelist_check.json').read_text())


@pytest.mark.parametrize("case", range(len(_PROBE)))
def test_matches_compiled_routine(case):
    c = _PROBE[case]['inputs']
    expected = [int(v) for v in _PROBE[case]['outputs']['OUT']]
    if c['routine'] == 'VNAME':
        l, i, found, array, ndms, nf = vname(c['card'], c['l'], c['names'],
                                             c['nf'])
        assert [l, i, int(found), int(array), ndms, nf] == expected
    else:
        check = lvalue if c['routine'] == 'LVALUE' else rvalue
        assert list(check(c['card'], c['l'], 0, c['nf'])) == expected


def test_clean_cards_have_no_faults():
    assert rvalue('1.0,2.5,-3.E2,4.5E-1$', 1, 0, 0) == (21, 4, 0)
    assert lvalue('3*.TRUE.,.FALSE.$', 1, 0, 0) == (17, 4, 0)
    assert vname('  ALSCHD(3)=1.,', 3, ['ALSCHD'], 0)[:5] == \
        (12, 1, True, True, 2)


def test_integer_values_count_as_faults():
    """RVALUE wants a decimal point or an exponent on every value."""
    assert rvalue('1,2$', 1, 0, 0)[2] == 2


# --- TESTOR, over a card sequence (its END and counts are saved) ----------

from pydatcom.io.namelist_check import testor as run_testor  # noqa: E402

_TESTOR = json.loads((pathlib.Path(__file__).resolve().parent / 'fixtures' /
                      'probes' / 'testor.json').read_text())


def _testor_replay():
    names = [n for n, _ in _TESTOR['names']]
    ldm = [d for _, d in _TESTOR['names']]
    state, out = {}, []
    for text, col, nam, k, ier in _TESTOR['cards']:
        out.append(run_testor(text.ljust(80)[:80] + '####', col, nam, k, ier,
                          names, ldm, state))
    return out


_TESTOR_RESULTS = _testor_replay()


@pytest.mark.parametrize("case", range(len(_TESTOR['cards'])))
def test_testor_matches_compiled_routine(case):
    r, rec = _TESTOR_RESULTS[case], _TESTOR['records'][case]
    assert r['unit6'] == rec['unit6']
    assert [ln.rstrip() for ln in r['unit11'] if ln.strip()] == rec['unit11']
    assert r['l'] == rec['l']


def test_testor_adds_a_missing_termination():
    k = next(i for i, c in enumerate(_TESTOR['cards'])
             if c[0].startswith(' $FLTCON NMACH=1.0'))
    assert 'MISSING NAMELIST TERMINATION' in _TESTOR_RESULTS[k]['unit6'][0]
    assert _TESTOR_RESULTS[k]['unit11'][0].startswith(' $END')
