"""
Regression tests for the TABLEC translation.

Source of truth: datcom-legacy/datcom_2000/tablec.f.

The 266 coefficient values were extracted from the source DATA statements by
parsing rather than transcribed by hand. These tests re-parse the FORTRAN
directly and compare, so the embedded table is checked against the file it
came from rather than against itself.
"""

import pathlib
import re

import numpy as np
import pytest

from pydatcom.aerodynamics.tablec import (
    calculate_tablec, _CCM, _MACH, _COEFFICIENTS, _MACH_POINTS,
)

_SOURCE = (pathlib.Path(__file__).resolve().parent.parent /
           'datcom-legacy' / 'datcom_2000' / 'tablec.f')


def _parse_source_block(name):
    """Pull one DATA array out of the FORTRAN and expand it."""
    text = _SOURCE.read_text(errors='replace')
    collected, active = [], False
    for line in text.splitlines():
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
    values = []
    for token in ''.join(collected).split(','):
        token = token.strip()
        if token:
            values.append(float(token))
    return values


pytestmark = pytest.mark.skipif(
    not _SOURCE.exists(),
    reason="legacy tablec.f not present")


# --------------------------------------------------------------------------
# The embedded table against the source file
# --------------------------------------------------------------------------

def test_mach_grid_matches_the_source():
    assert _MACH == pytest.approx(_parse_source_block('ZM'))


def test_every_coefficient_value_matches_the_source():
    """All 266 entries, re-parsed from the FORTRAN and compared in order."""
    expected = []
    for block in ('C1', 'C2', 'C3', 'C4'):
        expected.extend(_parse_source_block(block))
    assert len(expected) == _COEFFICIENTS * _MACH_POINTS == 266
    assert _CCM.ravel() == pytest.approx(expected, rel=1e-12)


def test_block_lengths_match_the_source_declarations():
    for block, declared in (('C1', 70), ('C2', 70), ('C3', 70), ('C4', 56)):
        assert len(_parse_source_block(block)) == declared


def test_table_shape():
    assert _CCM.shape == (_COEFFICIENTS, _MACH_POINTS) == (19, 14)


# --------------------------------------------------------------------------
# Interpolation
# --------------------------------------------------------------------------

@pytest.mark.parametrize("index", range(14))
def test_each_grid_mach_returns_its_own_column(index):
    result = calculate_tablec(_MACH[index])
    assert result['c'] == pytest.approx(_CCM[:, index], abs=1e-12)


def test_all_nineteen_coefficients_are_returned():
    assert calculate_tablec(1.0)['c'].shape == (19,)


def test_interior_mach_lies_between_its_neighbours():
    """Linear interpolation keeps each coefficient inside its bracket."""
    low, high = _MACH[4], _MACH[5]          # 0.90 and 0.95
    middle = calculate_tablec((low + high) / 2.0)['c']
    for index in range(_COEFFICIENTS):
        bounds = sorted((_CCM[index, 4], _CCM[index, 5]))
        assert bounds[0] - 1e-12 <= middle[index] <= bounds[1] + 1e-12


def test_mach_outside_the_grid_is_clamped():
    """End mode 0 at both ends, so no extrapolation."""
    assert calculate_tablec(0.05)['c'] == pytest.approx(_CCM[:, 0])
    assert calculate_tablec(9.0)['c'] == pytest.approx(_CCM[:, -1])


def test_mach_grid_is_increasing():
    assert np.all(np.diff(_MACH) > 0)
    assert _MACH[0] == 0.40 and _MACH[-1] == 2.50
