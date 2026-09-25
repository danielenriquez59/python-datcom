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


def tlip2x(x1, x2, packed: Sequence[int], query1: float, query2: float,
           lower1: int = 0, lower2: int = 0, upper1: int = 0,
           upper2: int = 0, shape: Optional[Sequence[int]] = None,
           nx2: Optional[int] = None) -> float:
    """Translate TLIP2X: TLINEX over a packed table ``Y(NX2, NX1)``.

    Each X1 column is interpolated in X2 by :func:`tlip1x`, its ordinates
    offset by YUP's shape array with the column index in word 5; the
    columns are then interpolated in X1 with TLIN1X's rules.

    Args:
        x1, x2: Abscissas.  packed: The packed ordinate words.
        query1, query2: The point.  lower1..upper2: The end modes
            ``LX1L``, ``LX2L``, ``LX1U``, ``LX2U``.
        shape: The caller's ``NXX``: with a first word of zero or more,
            that word is ``NX1``; with a negative first word (a slice of a
            larger table) ``NX1`` is its second word and its offsets are
            kept.  Omitted, ``NX1`` is the length of ``x1``.
        nx2: ``NX2``, the ordinates per column; omitted, the length of
            ``x2``.

    Returns:
        The interpolated value.  The source's MESSGE diagnostics on
        extrapolation are not reproduced, as for TLINEX.
    """
    x1 = np.asarray(x1, dtype=float)
    nx2 = len(x2) if nx2 is None else int(nx2)
    if shape is not None and shape[0] < 0:
        nx1 = int(shape[1])
        base = [int(v) for v in shape[:7]]
    else:
        nx1 = x1.size if shape is None else int(shape[0])
        base = [-nx2, nx1, 0, 0, 0, 0, 0]
    columns = []
    for column in range(1, nx1 + 1):
        np1 = list(base)
        np1[4] = column
        columns.append(tlip1x(x2[:nx2], packed, query2, lower2, upper2,
                              np1))
    return float(tlin1x(x1[:nx1], columns, query1, lower1, upper1))


def tlip3x(x1, x2, x3, packed: Sequence[int], query1: float,
           query2: float, query3: float, lower1: int = 0, lower2: int = 0,
           lower3: int = 0, upper1: int = 0, upper2: int = 0,
           upper3: int = 0) -> float:
    """Translate TLIP3X: TLIN3X over a packed table ``Y(NX2, NX1, NX3)``.

    Each X3 slice is interpolated by :func:`tlip2x`, reached through YUP's
    shape array ``(-NX2, NX1, NX3, 0, 0, k, 0)``, and the slices are then
    interpolated in X3 with TLIN1X's rules.  The MESSGE diagnostics are not
    reproduced.
    """
    nx1, nx2, nx3 = len(x1), len(x2), len(x3)
    slices = [tlip2x(x1, x2, packed, query1, query2, lower1, lower2,
                     upper1, upper2,
                     shape=[-nx2, nx1, nx3, 0, 0, x3_slot, 0])
              for x3_slot in range(1, nx3 + 1)]
    return float(tlin1x(np.asarray(x3, dtype=float), slices, query3,
                        lower3, upper3))


def intep3(arg1: float, arg2: float, lamda: float, charts: dict,
           state: Optional[dict] = None) -> float:
    """Translate INTEP3: interpolate between packed charts drawn at
    ``LAMDA`` = 0, 0.25, 0.5, 0.75 and 1.

    Args:
        arg1, arg2: The point on each chart.  lamda: The chart parameter.
        charts: ``'a'`` to ``'e'`` (for 0 to 1), each a dict with ``x1``,
            ``x2``, ``packed``, ``nxx`` (the source's ``N1``, ``NXX``
            words) and ``nx2``, or ``None`` for a caller's dummy chart
            (``0,0,0,1,1,0``), which evaluates to 0.
        state: The saved local ``it``, the chart pair last used; a
            ``LAMDA`` outside [0, 1] reuses it.

    Returns:
        ``ANS``, linear in ``LAMDA`` between the two charts bracketing it.
        When chart D has a single X2 point, a ``LAMDA`` above 0.5 pairs C
        with E instead.
    """
    st = state if state is not None else {}
    it = st.get('it', 1)
    if lamda == 0.:
        it = 1
    if 0.0 < lamda <= 0.25:
        it = 2
    if 0.25 < lamda <= 0.50:
        it = 3
    if 0.50 < lamda <= 0.75:
        it = 4
    if 0.75 < lamda <= 1.00:
        it = 5
    if it in (4, 5) and charts['d']['nx2'] == 1:
        it = 6
    st['it'] = it

    def look(key):
        c = charts[key]
        if c is None:
            # A dummy chart (``0,0,0,1,1,0`` in the call): one point, packed
            # word 0, which TLIP2X evaluates to 0.
            return 0.0
        return tlip2x(c['x1'], c['x2'], c['packed'], arg1, arg2,
                      shape=c['nxx'], nx2=c['nx2'])

    # (first chart, its LAMDA), (second chart, its LAMDA) per IT.
    pairs = {1: (('a', 0.0), ('a', 1.0)), 2: (('a', 0.0), ('b', 0.25)),
             3: (('b', 0.25), ('c', 0.50)), 4: (('c', 0.50), ('d', 0.75)),
             5: (('d', 0.75), ('e', 1.00)), 6: (('c', 0.50), ('e', 1.00))}
    (k1, x1), (k2, x2) = pairs[it]
    aa1, aa2 = look(k1), look(k2)
    arg = (lamda - x1) / (x2 - x1)
    return aa1 + (aa2 - aa1) * arg
