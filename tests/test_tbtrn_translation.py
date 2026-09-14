"""
Regression tests for the TBTRN translation.

Source of truth: datcom-legacy/datcom_2000/tbtrn.f.

The 1376 coefficients were extracted by parsing rather than transcribed, and
this test re-parses the FORTRAN independently to check them.
"""

import pathlib
import re

import numpy as np
import pytest

from pydatcom.aerodynamics.tbtrn import (
    calculate_tbtrn, _BTRN, _BLOCK, _MACH_OFFSET, _MACH_GROUPS, _offset,
)

_SOURCE = (pathlib.Path(__file__).resolve().parent.parent /
           'datcom-legacy' / 'datcom_2000' / 'tbtrn.f')

pytestmark = pytest.mark.skipif(
    not _SOURCE.exists(), reason="legacy tbtrn.f not present")


def _parse_source():
    text = _SOURCE.read_text(errors='replace')
    values = []
    for index in range(23):
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
    assert len(expected) == 1376
    assert _BTRN == pytest.approx(expected, rel=1e-12)


def test_offset_chain_matches_the_source():
    """L = (IA-2)*16, then +288 at IM=6, +464 at IM=7, +(IM-8)*240+656 above."""
    assert _offset(5, 2) == 0
    assert _offset(6, 2) == 288
    assert _offset(7, 2) == 464
    assert _offset(8, 2) == 656
    assert _offset(9, 2) == 896
    assert _offset(10, 2) == 1136
    assert _offset(8, 5) == (5 - 2) * _BLOCK + 656


def test_mach_indices_one_through_five_address_the_same_block():
    """None of them triggers an offset test, so all leave L at zero.

    The same kind of aliasing TBSUB shows between its first two indices.
    """
    offsets = {calculate_tbtrn(mach, 4)['offset'] for mach in range(1, 6)}
    assert len(offsets) == 1


def test_mach_blocks_are_unequal():
    """The offset chain gives blocks of 18, 11, 12 and 15 angle groups."""
    assert _MACH_GROUPS[1] == 18
    assert _MACH_GROUPS[6] == 11
    assert _MACH_GROUPS[7] == 12
    assert all(_MACH_GROUPS[mach] == 15 for mach in (8, 9, 10))


def test_every_stored_coefficient_is_reachable():
    """Unlike TBSUP, this table has no dead region."""
    reachable = set()
    for mach, groups in _MACH_GROUPS.items():
        for angle in range(2, 2 + groups):
            offset = calculate_tbtrn(mach, angle)['offset']
            reachable.update(range(offset, offset + _BLOCK))
    assert len(reachable) == _BTRN.size == 1376


def test_block_sizes_match_their_declared_starts():
    """Each Mach block runs exactly up to the next block's start."""
    starts = sorted(set(_MACH_OFFSET.values()))
    for index, start in enumerate(starts[:-1]):
        span = starts[index + 1] - start
        mach = next(m for m, o in _MACH_OFFSET.items() if o == start)
        assert span == _MACH_GROUPS[mach] * _BLOCK


def test_first_angle_index_returns_a_zero_block():
    for mach in (1, 6, 10):
        result = calculate_tbtrn(mach, 1)
        assert result['zeroed']
        assert np.all(result['coefficients'] == 0.0)


def test_each_selection_returns_sixteen_coefficients():
    assert calculate_tbtrn(8, 10)['coefficients'].shape == (_BLOCK,)


def test_returned_block_is_a_copy():
    block = calculate_tbtrn(8, 5)['coefficients']
    original = block[0]
    block[0] = 999.0
    assert calculate_tbtrn(8, 5)['coefficients'][0] == pytest.approx(original)


def test_out_of_range_selection_is_rejected():
    with pytest.raises(ValueError, match="outside"):
        calculate_tbtrn(10, 40)
    with pytest.raises(ValueError, match="outside"):
        calculate_tbtrn(1, 0)
