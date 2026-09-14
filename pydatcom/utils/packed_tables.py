"""
Packed table storage: YUP and TLIP1X.

``YUP`` unpacks the compression scheme the legacy tables use to halve their
storage: two floating-point numbers, each carrying three significant digits,
share one integer word along with their signs and decimal exponents.

Reading the packed word as a positive integer, its decimal digits hold, from
the right:

===========  ==========================================
positions     contents
===========  ==========================================
1 to 3        the second number's three digits
4 to 6        the first number's three digits
7             the second number's decimal exponent
8             the first number's decimal exponent
9             which exponents are negative
10            the second number's sign
===========  ==========================================

The first number's sign is the sign of the word itself.  Odd one-based
indices select the first number of a word, even indices the second.

``TLIP1X`` is then ``TLIN1X`` reading through that indirection: identical
bracket, snapping and extrapolation behaviour, with every ordinate fetched
by ``YUP`` instead of indexed directly.

Reference: datcom-legacy/datcom_2000/yup.f, tlip1x.f
"""

import numpy as np
from typing import Dict, Optional, Sequence
import logging

from pydatcom.utils.legacy_tables import tlin1x

logger = logging.getLogger(__name__)

# Which exponents the NSPS digit marks negative.
_BOTH_EXPONENTS_NEGATIVE = 1
_SECOND_EXPONENT_NEGATIVE = 2


def yup(index: int, packed: Sequence[int],
        shape: Optional[Sequence[int]] = None) -> float:
    """Translate YUP: unpack one ordinate from a packed table.

    Args:
        index: One-based position in the logical (unpacked) table, the
            source's ``I``.
        packed: The packed integer words, the source's ``NY``.
        shape: The source's ``M`` array, holding ``NX2, NX1, NX3, NX4, I1,
            I3, I4``.  A negative first entry selects the multidimensional
            offset calculation; omitting it treats ``index`` as a flat
            position.

    Returns:
        The unpacked value.

    Raises:
        ValueError: If the index falls outside the packed array.
    """
    raw = int(index)
    if shape is not None and len(shape) >= 7 and shape[0] < 0:
        first = abs(int(shape[0]))
        second = first * int(shape[1])
        third = second * int(shape[2])
        if shape[1] != 0:
            raw += first * (int(shape[4]) - 1)
        if shape[2] != 0:
            raw += second * (int(shape[5]) - 1)
        if shape[3] != 0:
            raw += third * (int(shape[6]) - 1)

    # Odd positions take the first number of a word, even the second.
    word_index = raw // 2
    if raw != 2 * word_index:
        word_index += 1
        take_first = True
    else:
        take_first = False
    if not 1 <= word_index <= len(packed):
        raise ValueError(
            f"YUP index {index} resolves to packed word {word_index}, "
            f"outside the {len(packed)}-word array")

    word = int(packed[word_index - 1])
    sign_first = -1 if word < 0 else 1
    word = abs(word)

    digits_second = word % 1000
    rest = word // 1000
    digits_first = rest % 1000
    rest //= 1000
    exponent_second = rest % 10
    rest //= 10
    exponent_first = rest % 10
    rest //= 10
    exponent_signs = rest % 10
    rest //= 10
    sign_second = rest % 10

    value_first = float(digits_first)
    if sign_first < 0:
        value_first = -value_first
    value_second = float(digits_second)
    if sign_second == 1:
        value_second = -value_second

    if exponent_signs == _BOTH_EXPONENTS_NEGATIVE:
        exponent_first = -exponent_first
        exponent_second = -exponent_second
    elif exponent_signs == _SECOND_EXPONENT_NEGATIVE:
        exponent_second = -exponent_second
    elif exponent_signs != 0:
        exponent_first = -exponent_first

    value_first *= 10.0**exponent_first
    value_second *= 10.0**exponent_second
    return float(value_first if take_first else value_second)


def unpack_table(packed: Sequence[int], count: int,
                 shape: Optional[Sequence[int]] = None) -> np.ndarray:
    """Unpack ``count`` consecutive ordinates into a plain array.

    Args:
        packed: The packed integer words.
        count: How many logical values to recover.
        shape: The source's ``M`` array, passed through to :func:`yup`.

    Returns:
        The unpacked values, in logical order.
    """
    return np.array([yup(position, packed, shape)
                     for position in range(1, int(count) + 1)])


def tlip1x(x, packed: Sequence[int], query: float,
           lower: int = 0, upper: int = 0,
           shape: Optional[Sequence[int]] = None) -> float:
    """Translate TLIP1X: TLIN1X over a packed table.

    The interpolation itself is identical to ``TLIN1X`` -- same bracket,
    snapping, clamping and linear or quadratic extrapolation -- so this
    unpacks the ordinates and defers to the existing translation rather than
    duplicating that logic.

    Args:
        x: Abscissas.
        packed: Packed ordinate words.
        query: Where to interpolate.
        lower: Lower end mode, the source's ``LX1L``.
        upper: Upper end mode, ``LX1U``.
        shape: The source's ``M`` array, passed through to :func:`yup`.

    Returns:
        The interpolated value.

    Raises:
        ValueError: If the abscissa array is unusable.
    """
    x = np.asarray(x, dtype=float)
    if x.ndim != 1 or x.size < 2:
        raise ValueError("TLIP1X needs at least two abscissas")
    y = unpack_table(packed, x.size, shape)
    return float(tlin1x(x, y, query, lower, upper))
