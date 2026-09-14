"""
Regression tests for the YUP and TLIP1X translations.

Source of truth: datcom-legacy/datcom_2000/yup.f, tlip1x.f.

YUP unpacks a 1970s storage scheme that fits two three-significant-digit
floats into one integer word. These tests build packed words from a
reference packer written directly from the documented digit layout, so the
unpacker is checked against an independent construction rather than against
itself.
"""

import numpy as np
import pytest

from pydatcom.utils.packed_tables import (
    yup, unpack_table, tlip1x,
    _BOTH_EXPONENTS_NEGATIVE, _SECOND_EXPONENT_NEGATIVE,
)


def _pack(digits1, exponent1, digits2, exponent2,
          negative1=False, negative2=False, exponent_signs=0):
    """Build a packed word from the documented digit positions."""
    word = (exponent_signs * 10**8 + exponent1 * 10**7 +
            exponent2 * 10**6 + digits1 * 1000 + digits2)
    if negative2:
        word += 10**9
    return -word if negative1 else word


# --------------------------------------------------------------------------
# The packing layout
# --------------------------------------------------------------------------

def test_two_values_share_one_word():
    word = _pack(123, 0, 456, 0)
    assert yup(1, [word]) == pytest.approx(123.0)
    assert yup(2, [word]) == pytest.approx(456.0)


def test_odd_indices_take_the_first_value():
    words = [_pack(100, 0, 200, 0), _pack(300, 0, 400, 0)]
    assert yup(1, words) == pytest.approx(100.0)
    assert yup(3, words) == pytest.approx(300.0)


def test_even_indices_take_the_second_value():
    words = [_pack(100, 0, 200, 0), _pack(300, 0, 400, 0)]
    assert yup(2, words) == pytest.approx(200.0)
    assert yup(4, words) == pytest.approx(400.0)


# --------------------------------------------------------------------------
# Signs
# --------------------------------------------------------------------------

def test_the_first_value_takes_the_sign_of_the_word():
    word = _pack(123, 0, 456, 0, negative1=True)
    assert yup(1, [word]) == pytest.approx(-123.0)
    assert yup(2, [word]) == pytest.approx(456.0)


def test_the_second_value_has_its_own_sign_digit():
    word = _pack(123, 0, 456, 0, negative2=True)
    assert yup(1, [word]) == pytest.approx(123.0)
    assert yup(2, [word]) == pytest.approx(-456.0)


def test_both_values_can_be_negative_at_once():
    word = _pack(123, 0, 456, 0, negative1=True, negative2=True)
    assert yup(1, [word]) == pytest.approx(-123.0)
    assert yup(2, [word]) == pytest.approx(-456.0)


# --------------------------------------------------------------------------
# Exponents and the NSPS selector
# --------------------------------------------------------------------------

def test_positive_exponents_scale_both_values():
    word = _pack(123, 2, 456, 3)
    assert yup(1, [word]) == pytest.approx(12300.0)
    assert yup(2, [word]) == pytest.approx(456000.0)


def test_nsps_one_makes_both_exponents_negative():
    word = _pack(123, 2, 456, 3, exponent_signs=_BOTH_EXPONENTS_NEGATIVE)
    assert yup(1, [word]) == pytest.approx(1.23)
    assert yup(2, [word]) == pytest.approx(0.456)


def test_nsps_two_makes_only_the_second_exponent_negative():
    word = _pack(123, 2, 456, 3, exponent_signs=_SECOND_EXPONENT_NEGATIVE)
    assert yup(1, [word]) == pytest.approx(12300.0)
    assert yup(2, [word]) == pytest.approx(0.456)


def test_nsps_above_two_makes_only_the_first_exponent_negative():
    """The source's fall-through branch covers every remaining value."""
    for selector in (3, 4, 9):
        word = _pack(123, 2, 456, 3, exponent_signs=selector)
        assert yup(1, [word]) == pytest.approx(1.23)
        assert yup(2, [word]) == pytest.approx(456000.0)


def test_nsps_zero_leaves_both_exponents_positive():
    word = _pack(123, 2, 456, 3, exponent_signs=0)
    assert yup(1, [word]) == pytest.approx(12300.0)
    assert yup(2, [word]) == pytest.approx(456000.0)


# --------------------------------------------------------------------------
# Multidimensional offsets
# --------------------------------------------------------------------------

def test_negative_first_shape_entry_adds_the_block_offset():
    """M(1) < 0 selects the source's I1/I3/I4 offset calculation."""
    words = [_pack(100, 0, 200, 0), _pack(300, 0, 400, 0),
             _pack(500, 0, 600, 0), _pack(700, 0, 800, 0)]
    # NX2=-2, NX1=2, NX3=0, NX4=0, I1=2 -> offset of 2 logical positions.
    shape = [-2, 2, 0, 0, 2, 1, 1]
    assert yup(1, words, shape) == pytest.approx(yup(3, words))
    assert yup(2, words, shape) == pytest.approx(yup(4, words))


def test_positive_first_shape_entry_is_a_flat_index():
    words = [_pack(100, 0, 200, 0), _pack(300, 0, 400, 0)]
    shape = [2, 2, 0, 0, 2, 1, 1]
    assert yup(3, words, shape) == pytest.approx(300.0)


# --------------------------------------------------------------------------
# Bulk unpacking and TLIP1X
# --------------------------------------------------------------------------

def test_unpack_table_recovers_the_logical_order():
    words = [_pack(100, 0, 200, 0), _pack(300, 0, 400, 0),
             _pack(500, 0, 600, 0)]
    assert unpack_table(words, 6) == pytest.approx(
        [100., 200., 300., 400., 500., 600.])


def test_tlip1x_interpolates_the_unpacked_values():
    words = [_pack(100, 0, 200, 0), _pack(300, 0, 400, 0),
             _pack(500, 0, 600, 0)]
    x = [0., 1., 2., 3., 4., 5.]
    assert tlip1x(x, words, 2.5) == pytest.approx(350.0)
    assert tlip1x(x, words, 2.0) == pytest.approx(300.0)


def test_tlip1x_matches_tlin1x_on_the_unpacked_data():
    """TLIP1X is TLIN1X with an indirection, so the two must agree."""
    from pydatcom.utils.legacy_tables import tlin1x
    words = [_pack(100, 0, 250, 0), _pack(310, 0, 400, 0),
             _pack(505, 0, 600, 0)]
    x = [0., 1., 2., 3., 4., 5.]
    y = unpack_table(words, 6)
    for query in (0.5, 1.7, 3.3, 4.9):
        assert tlip1x(x, words, query) == pytest.approx(tlin1x(x, y, query))


def test_out_of_range_index_is_rejected():
    with pytest.raises(ValueError, match="outside"):
        yup(99, [_pack(100, 0, 200, 0)])


def test_tlip1x_rejects_a_short_abscissa_array():
    with pytest.raises(ValueError):
        tlip1x([1.0], [_pack(100, 0, 200, 0)], 1.0)
