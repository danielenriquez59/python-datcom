"""
Regression tests for the TBSUP translation.

Source of truth: datcom-legacy/datcom_2000/tbsup.f.

The 992 coefficients were extracted by parsing rather than transcribed, and
this test re-parses the FORTRAN independently to check them.
"""

import pathlib
import re

import numpy as np
import pytest

from pydatcom.aerodynamics.tbsup import (
    calculate_tbsup, _BSUP, _BLOCK, _ANGLE_GROUPS, _MACH_STRIDE,
    _MACH_BASE, _BUMP_ABOVE, _BUMP, _UNREACHABLE,
)

_SOURCE = (pathlib.Path(__file__).resolve().parent.parent /
           'datcom-legacy' / 'datcom_2000' / 'tbsup.f')

pytestmark = pytest.mark.skipif(
    not _SOURCE.exists(), reason="legacy tbsup.f not present")


def _parse_source():
    text = _SOURCE.read_text(errors='replace')
    values = []
    for index in range(16):
        name = f'B{index}'
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
        raw = re.sub(r'E\s+(\d)', r'E+\1', ''.join(collected))
        values.extend(float(t.strip()) for t in raw.split(',') if t.strip())
    return values


def test_every_coefficient_matches_the_source():
    expected = _parse_source()
    assert len(expected) == 992
    assert _BSUP == pytest.approx(expected, rel=1e-12)


def test_offset_follows_the_source_expression():
    assert calculate_tbsup(11, 2)['offset'] == 0
    assert calculate_tbsup(12, 2)['offset'] == _MACH_STRIDE
    assert calculate_tbsup(13, 5)['offset'] == (
        (5 - 2) * _BLOCK + (13 - _MACH_BASE) * _MACH_STRIDE + _BUMP)


def test_the_bump_applies_only_above_mach_index_twelve():
    below = calculate_tbsup(_BUMP_ABOVE, 2)['offset']
    above = calculate_tbsup(_BUMP_ABOVE + 1, 2)['offset']
    assert above - below == _MACH_STRIDE + _BUMP


def test_thirty_two_coefficients_are_unreachable():
    """The +32 bump skips a whole block of the stored table.

    Mach indices 11 and 12 cover positions 0 to 479, and 13 and 14 resume at
    512, so 480 to 511 are never read. This is the mirror image of WINGCL's
    defect, where the table is shorter than its indexing needs; here it is
    longer than the indexing uses.
    """
    reachable = set()
    for mach in range(_MACH_BASE, _MACH_BASE + 4):
        for angle in range(2, 2 + _ANGLE_GROUPS):
            offset = calculate_tbsup(mach, angle)['offset']
            reachable.update(range(offset, offset + _BLOCK))
    gap = sorted(set(range(_BSUP.size)) - reachable)
    assert len(reachable) == 960
    assert len(gap) == _BUMP == 32
    assert gap == list(_UNREACHABLE)


def test_reachable_blocks_do_not_overlap():
    seen = set()
    for mach in range(_MACH_BASE, _MACH_BASE + 4):
        for angle in range(2, 2 + _ANGLE_GROUPS):
            offset = calculate_tbsup(mach, angle)['offset']
            assert offset not in seen
            seen.add(offset)
    assert len(seen) == 4 * _ANGLE_GROUPS


def test_first_angle_index_returns_a_zero_block():
    for mach in range(_MACH_BASE, _MACH_BASE + 4):
        result = calculate_tbsup(mach, 1)
        assert result['zeroed']
        assert np.all(result['coefficients'] == 0.0)


def test_each_selection_returns_sixteen_coefficients():
    assert calculate_tbsup(13, 10)['coefficients'].shape == (_BLOCK,)


def test_returned_block_is_a_copy():
    block = calculate_tbsup(13, 5)['coefficients']
    original = block[0]
    block[0] = 999.0
    assert calculate_tbsup(13, 5)['coefficients'][0] == pytest.approx(original)


def test_out_of_range_selection_is_rejected():
    with pytest.raises(ValueError, match="outside"):
        calculate_tbsup(20, 16)
    with pytest.raises(ValueError, match="outside"):
        calculate_tbsup(_MACH_BASE, 0)
