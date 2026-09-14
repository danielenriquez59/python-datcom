"""
Regression tests for the TBSUB translation.

Source of truth: datcom-legacy/datcom_2000/tbsub.f.

As with TABLEC, the 864 coefficients were extracted by parsing rather than
transcribed, and this test re-parses the FORTRAN independently to check
them. The source writes some exponents as "E 00" with a space where the
sign would be, which the parser has to accept.
"""

import pathlib
import re

import numpy as np
import pytest

from pydatcom.aerodynamics.tbsub import (
    calculate_tbsub, _BSUB, _BLOCK, _ANGLE_GROUPS, _MACH_STRIDE,
)

_SOURCE = (pathlib.Path(__file__).resolve().parent.parent /
           'datcom-legacy' / 'datcom_2000' / 'tbsub.f')

pytestmark = pytest.mark.skipif(
    not _SOURCE.exists(), reason="legacy tbsub.f not present")


def _parse_source():
    text = _SOURCE.read_text(errors='replace')
    values = []
    for index in range(14):
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


# --------------------------------------------------------------------------
# The embedded table against the source file
# --------------------------------------------------------------------------

def test_every_coefficient_matches_the_source():
    expected = _parse_source()
    assert len(expected) == 864
    assert _BSUB == pytest.approx(expected, rel=1e-12)


def test_table_dimensions_factor_as_the_source_indexing_implies():
    """864 = 3 Mach blocks of 18 angle groups of 16 coefficients."""
    assert _BSUB.size == 864
    assert _BLOCK == 16 and _ANGLE_GROUPS == 18
    assert _MACH_STRIDE == _BLOCK * _ANGLE_GROUPS == 288
    assert _BSUB.size == 3 * _MACH_STRIDE


def test_space_exponent_format_is_parsed():
    """The source writes at least one value as '-4.5359688E 00'."""
    text = _SOURCE.read_text(errors='replace')
    assert re.search(r'E\s+\d\d', text)
    assert np.all(np.isfinite(_BSUB))


# --------------------------------------------------------------------------
# Selection
# --------------------------------------------------------------------------

def test_first_angle_index_returns_a_zero_block():
    """The source fills BT with zeros when IA is one."""
    for mach in (1, 2, 3, 4):
        result = calculate_tbsub(mach, 1)
        assert result['zeroed']
        assert np.all(result['coefficients'] == 0.0)
        assert result['offset'] is None


def test_offset_follows_the_source_expression():
    assert calculate_tbsub(3, 5)['offset'] == (5 - 2) * 16 + (3 - 2) * 288


def test_mach_indices_one_and_two_address_the_same_block():
    """The source adds its Mach offset only from IM=2, and (2-2)*288 is zero.

    So IM=1 and IM=2 alias. That is a property of the source's indexing,
    not of the translation.
    """
    first = calculate_tbsub(1, 7)
    second = calculate_tbsub(2, 7)
    assert first['offset'] == second['offset']
    assert first['coefficients'] == pytest.approx(second['coefficients'])


def test_each_selection_returns_sixteen_coefficients():
    assert calculate_tbsub(3, 10)['coefficients'].shape == (16,)


def test_the_full_table_is_reachable():
    """IM=2..4 with IA=2..19 covers every block exactly once."""
    seen = set()
    for mach in (2, 3, 4):
        for angle in range(2, 2 + _ANGLE_GROUPS):
            offset = calculate_tbsub(mach, angle)['offset']
            assert offset not in seen
            seen.add(offset)
    assert len(seen) == 3 * _ANGLE_GROUPS
    assert min(seen) == 0 and max(seen) == _BSUB.size - _BLOCK


def test_returned_block_is_a_copy():
    """Mutating the result must not corrupt the table."""
    block = calculate_tbsub(3, 5)['coefficients']
    original = block[0]
    block[0] = 999.0
    assert calculate_tbsub(3, 5)['coefficients'][0] == pytest.approx(original)


def test_out_of_range_selection_is_rejected():
    with pytest.raises(ValueError, match="outside"):
        calculate_tbsub(5, 19)
    with pytest.raises(ValueError, match="outside"):
        calculate_tbsub(1, 0)
