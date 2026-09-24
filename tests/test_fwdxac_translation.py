"""
Regression tests for the FWDXAC translation.

Source of truth: datcom-legacy/datcom_2000/fwdxac.f.

The four tables were extracted from the source DATA statements by parsing;
these tests re-parse the FORTRAN independently and compare every value.  The
lookups are checked against a compiled probe of the legacy routine, see
tools/probes/fwdxac.py.
"""

import json
import pathlib
import re

import numpy as np
import pytest

from pydatcom.aerodynamics.fwdxac import calculate_fwdxac, _TABLES

_ROOT = pathlib.Path(__file__).resolve().parent.parent
_SOURCE = _ROOT / 'datcom-legacy' / 'datcom_2000' / 'fwdxac.f'
_PROBE = json.loads(
    (_ROOT / 'tests' / 'fixtures' / 'probes' / 'fwdxac.json').read_text())


def _parse_source_block(name):
    """Pull one DATA array out of the FORTRAN and expand it."""
    collected, active = [], False
    for line in _SOURCE.read_text(errors='replace').splitlines():
        if re.match(r'\s*DATA\s+' + name + r'\s*/', line):
            active = True
            line = line.split('/', 1)[1]
        elif active:
            line = line[6:] if len(line) > 6 else ''
        if active:
            if '/' in line:
                collected.append(line.split('/')[0])
                break
            collected.append(line)
    return [float(t) for t in ''.join(collected).split(',') if t.strip()]


def _source_table(prefix, first):
    """A table in the source's flat order, from its six EQUIVALENCEd parts."""
    values = []
    for k in range(6):
        values.extend(_parse_source_block(f'{prefix}{first + k}'))
    return values


def _flat(table):
    """Undo the TLIN3X reshape, back to the source's storage order."""
    return _TABLES[table].transpose(2, 1, 0).ravel()


@pytest.mark.parametrize("table,prefix,first", [
    ('TYSUBL', 'SUBT', 1), ('TYSUBR', 'SUBT', 7),
    ('TYSUPL', 'SUPT', 1), ('TYSUPR', 'SUPT', 7),
])
def test_every_table_value_matches_the_source(table, prefix, first):
    """All 216 entries of each table, in source storage order."""
    expected = _source_table(prefix, first)
    assert len(expected) == 216
    np.testing.assert_array_equal(_flat(table), expected)


@pytest.mark.parametrize("case", range(len(_PROBE)))
def test_matches_compiled_routine(case):
    """Every probe point, including extrapolation and clamping."""
    probe = _PROBE[case]
    assert calculate_fwdxac(**probe['inputs'])['xac'] == pytest.approx(
        probe['xac'], rel=1e-9, abs=1e-12)


def test_probe_reaches_every_branch():
    """Guard the fixture: all four regime/abscissa branches are probed."""
    reached = {(p['inputs']['mach'] < 1.0, abs(p['inputs']['factor']) > 1.0)
               for p in _PROBE}
    assert reached == {(True, True), (True, False), (False, True),
                       (False, False)}


def test_abscissa_turns_over_at_one():
    """tan/beta up to one, then beta/tan, so the abscissa never exceeds one."""
    assert calculate_fwdxac(3., .3, .4, .6)['abscissa'] == pytest.approx(.4)
    assert calculate_fwdxac(3., .3, 2.5, .6)['abscissa'] == pytest.approx(.4)
    assert calculate_fwdxac(3., .3, 1.0, .6)['abscissa'] == pytest.approx(1.)


def test_tysubr_is_the_table_that_joins_its_neighbours():
    """The data shows which table the subsonic beta/tan branch should read.

    TYSUBR meets TYSUBL exactly at an abscissa of one (tan/beta = 1), on
    every curve and taper, and TYSUPL at zero (Mach one) to chart-reading
    precision.  TYSUPR, which the source actually reads there, misses both
    by far more.
    """
    subr, subl = _TABLES['TYSUBR'], _TABLES['TYSUBL']
    supl, supr = _TABLES['TYSUPL'], _TABLES['TYSUPR']
    np.testing.assert_array_equal(subr[-1], subl[-1])
    np.testing.assert_allclose(subr[0], supl[0], atol=0.02 + 1e-12)
    assert np.max(np.abs(supr[-1] - subl[-1])) > 0.15
    assert np.max(np.abs(supr[0] - supl[0])) > 0.25


def test_source_defect_is_preserved_and_reported():
    """Subsonic tan/beta > 1 reads TYSUPR, as executed; TYSUBR is reported."""
    result = calculate_fwdxac(3.0, .3, 2.0, .6)
    assert result['table'] == 'TYSUPR'
    assert result['source_defect'] == 'subsonic_beta_over_tan_reads_TYSUPR'
    assert result['xac_tysubr'] != pytest.approx(result['xac'], abs=0.1)
    # The defect leaves a step at tan/beta = 1 that TYSUBR would not.
    below = calculate_fwdxac(3.0, .3, 1.0, .6)
    above = calculate_fwdxac(3.0, .3, 1.0 + 1e-9, .6)
    assert abs(above['xac'] - below['xac']) > 0.05
    assert above['xac_tysubr'] == pytest.approx(below['xac'], abs=1e-8)
    # No other branch carries the flag.
    for factor, mach in ((.5, .6), (.5, 1.6), (2.0, 1.6)):
        assert 'source_defect' not in calculate_fwdxac(3.0, .3, factor, mach)


def test_suspect_supt6_entry_is_preserved():
    """SUPT6 row five keeps its -.56, which breaks an otherwise smooth row."""
    row = _source_table('SUPT', 1)[5 * 36 + 4 * 6:5 * 36 + 5 * 6]
    assert row == [-.71, -.71, -.68, -.56, -.64, -.63]
    np.testing.assert_array_equal(_TABLES['TYSUPL'][:, 4, 5], row)
