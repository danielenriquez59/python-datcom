"""
Regression tests for XPERNM, M34O42 and TEST, the experimental-namelist
card counts.

Source of truth: datcom-legacy/datcom_2000/xpernm.f, m34o42.f, tbtrn.f.

Checked against a compiled probe (tools/probes/xpernm.py) that writes each
deck to unit 8 and runs XPERNM on it.
"""

import json
import pathlib

import pytest

from pydatcom.io.experimental_cards import m34o42, xpernm
from pydatcom.io.experimental_cards import test_name as name_matches

_ROOT = pathlib.Path(__file__).resolve().parent.parent
_PROBE = json.loads(
    (_ROOT / 'tests' / 'fixtures' / 'probes' / 'xpernm.json').read_text())


@pytest.mark.parametrize("case", range(len(_PROBE)))
def test_matches_compiled_routine(case):
    c, o = _PROBE[case]['inputs'], _PROBE[case]['outputs']
    r = xpernm(c['deck'], c['nlist'])
    assert [r['klist']] + r['nlist'] == [int(v) for v in o['R']]
    assert r['warning'] == ('_text' in o)


def test_name_needs_a_blank_after_it():
    assert name_matches('EXPR A=1', 'EXPR')
    assert not name_matches('EXPR01 A=1', 'EXPR')
    assert name_matches('EXPR', 'EXPR')


def test_count_runs_to_the_next_header():
    deck = [' $EXPR A=1$', ' $SYNTHS X=1$', ' $EXPR B=2$', 'CASEID']
    r = xpernm(deck, [0] * 100)
    assert r['nlist'][:2] == [2000, 2000]


def test_overlay_skips_without_a_list():
    assert m34o42(0, [' $EXPR A=1$'], [0] * 100) is None
    assert m34o42(1, [' $EXPR A=1$'], [0] * 100)['klist'] == 1
