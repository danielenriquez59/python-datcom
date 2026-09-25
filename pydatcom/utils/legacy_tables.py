"""Numerical contracts of DATCOM GLOOK, SWITCH, TLIN1X, and TLINEX.

Unlike generic interpolation, legacy table lookup snaps to grid points with
relative error below 0.001. Extrapolation controls refer to the first (lower)
and last (upper) *table ends*, including for descending grids. Modes <= 0
clamp, 1 extrapolates linearly, and > 1 extrapolates quadratically when three
points exist. Mode 0 also requests a legacy diagnostic, as do positive modes.

Python adaptations: indices are zero based; unused interpolation fractions
are None; invalid grids/shapes raise ValueError instead of out-of-bounds
Fortran access. Message printing and Hollerith ROUT/MESS bookkeeping are not
translated. SWITCH exposes the flags so callers can supply diagnostics.

The public functions check their grids and tables once, then hand plain
Python floats to :func:`_interp1`.  Tables here hold a handful of points, and
per-call NumPy overhead on arrays that small cost far more than the
arithmetic.  :func:`_interp1` reads ordinates through a getter, so a
multi-variable lookup interpolates only the columns its outer bracket needs,
as the FORTRAN label sequence does.
"""

import math
from typing import Callable, List, NamedTuple, Optional, Sequence

import numpy as np

from .legacy_numeric import _parabola


class Lookup(NamedTuple):
    index: int
    fraction: Optional[float]


class Switches(NamedTuple):
    no_interpolation: bool
    after_last: bool
    before_first: bool
    message_requested: bool
    extrapolate: bool
    ascending: bool
    use_extrapolation: bool


def _grid(x) -> List[float]:
    x = np.asarray(x, dtype=float)
    if x.ndim != 1 or not len(x):
        raise ValueError("A grid must be a nonempty finite one-dimensional array")
    values = x.tolist()
    if not all(map(math.isfinite, values)):
        raise ValueError("A grid must be a nonempty finite one-dimensional array")
    # A repeated coordinate is accepted, as GLOOK accepts it: a query can
    # only reach the zero-width interval by snapping to its first node.
    # TRANJT's Figure 6.3.2-30 grid repeats 4.04.
    if len(values) > 1:
        pairs = zip(values, values[1:])
        if values[0] < values[-1]:
            monotonic = all(a <= b for a, b in pairs)
        else:
            monotonic = (values[0] != values[-1] and
                         all(a >= b for a, b in pairs))
        if not monotonic:
            raise ValueError("Grid coordinates must be monotonic")
    return values


def _table(y, shape, message: str) -> list:
    """Check a table's shape and finiteness; return it as nested lists."""
    y = np.asarray(y, dtype=float)
    if y.shape != shape or not np.isfinite(y).all():
        raise ValueError(message)
    return y.tolist()


def _query(query) -> float:
    query = float(query)
    if math.isinf(query):
        raise ValueError("The query must be finite")
    return query


def _finite(value: float) -> float:
    """An inner lookup's result, checked as the next level checks its Y."""
    if not math.isfinite(value):
        raise ValueError("Y must be a finite array with the same shape as X")
    return value


def _divide(numerator: float, denominator: float) -> float:
    """Divide as NumPy does: a zero-width end interval gives inf or NaN."""
    if denominator != 0.0:
        return numerator / denominator
    return float(np.float64(numerator) / denominator)


def _glook(x: Sequence[float], query: float, ascending: bool) -> Lookup:
    previous = 0.0
    for index, coordinate in enumerate(x):
        difference = query - coordinate
        denominator = coordinate if query == 0.0 else query
        if abs(denominator) <= 0.0001:
            denominator = 1.0
        if abs(difference / denominator) < 1.e-3:
            return Lookup(index, None)
        if (ascending and difference < 0.0) or (not ascending and difference > 0.0):
            if index == 0:
                return Lookup(index, None)
            return Lookup(index, previous / (previous - difference))
        previous = difference
    return Lookup(len(x)-1, None)


def _interp1(x: Sequence[float], y: Callable[[int], float], query: float,
             lower: int, upper: int) -> float:
    """TLIN1X on a checked grid, reading ordinate ``i`` as ``y(i)``."""
    first, last = x[0], x[-1]
    ascending = first <= last
    if ascending:
        after, before = query > last, query < first
    else:
        after, before = query < last, query > first
    mode = upper if after else lower
    if not ((after or before) and mode > 0):
        index, fraction = _glook(x, query, ascending)
        if fraction is None:
            return y(index)
        low = y(index-1)
        return low + fraction*(y(index)-low)

    n = len(x)
    if n < 2:
        y(0)    # an inner lookup's own error takes precedence, as it did
        raise ValueError("Extrapolation requires at least two grid points")
    if mode > 1 and n > 2:
        points = (0, 1, 2) if before else (n-3, n-2, n-1)
        return _parabola([x[i] for i in points], [y(i) for i in points],
                         query)
    if before:
        y0, y1 = y(0), y(1)
        fraction = _divide(query-x[0], x[1]-x[0])
        return y0 + fraction*(y1-y0)
    y_before, y_last = y(n-2), y(n-1)
    fraction = _divide(query-x[-1], x[-1]-x[-2])
    return y_last + fraction*(y_last-y_before)


def glook(x, query: float, ascending: Optional[bool] = None) -> Lookup:
    """Translate GLOOK with its caller-initialized NOING=False contract.

    A None fraction means use Y[index] directly (snap or endpoint clamp).
    Otherwise interpolate between index-1 and index with the given fraction.

    A NaN query fails every comparison, so the source's search runs off the
    end and returns the last node, exactly; that is kept.  Infinite queries
    are rejected.
    """
    x = _grid(x)
    query = _query(query)
    if ascending is None:
        ascending = x[0] <= x[-1]
    return _glook(x, query, ascending)


def switch(x, query: float, lower: int = 0, upper: int = 0) -> Switches:
    """Translate SWITCH's seven LG flags, preserving table-end semantics."""
    x = _grid(x)
    query = _query(query)
    ascending = x[0] <= x[-1]
    after = query > x[-1] if ascending else query < x[-1]
    before = query < x[0] if ascending else query > x[0]
    mode = upper if after else lower
    outside = after or before
    message = outside and mode >= 0
    extrapolate = outside and mode > 0
    return Switches(False, after, before, message, extrapolate, ascending, extrapolate)


def tlin1x(x, y, query: float, lower: int = 0, upper: int = 0) -> float:
    """TLIN1X: snapped linear lookup with table-end extrapolation controls.

    SWITCH chooses extrapolation before GLOOK, so a positive mode bypasses
    snapping outside the grid, even arbitrarily close to an endpoint.
    Singleton tables support direct lookup/clamping; extrapolation needs
    at least two points (the legacy routine otherwise accesses outside X/Y).
    """
    x = _grid(x)
    y = _table(y, (len(x),),
               "Y must be a finite array with the same shape as X")
    return _interp1(x, y.__getitem__, _query(query), lower, upper)


def tlinex(x1, x2, y, query1: float, query2: float,
           lower1: int = 0, lower2: int = 0,
           upper1: int = 0, upper2: int = 0) -> float:
    """TLINEX: interpolate X2 columns, then interpolate their values in X1.

    Y must have shape ``(len(x2), len(x1))``, matching Fortran Y(NX2,NX1).
    Each axis uses TLIN1X's snapping and independent extrapolation modes.
    Only the X1 columns the outer lookup reads are interpolated, as in the
    legacy label sequence.
    """
    x1, x2 = _grid(x1), _grid(x2)
    rows = _table(y, (len(x2), len(x1)),
                  "Y must be finite with shape (len(x2), len(x1))")
    query1, query2 = _query(query1), _query(query2)

    def column(i):
        return _finite(_interp1(x2, lambda j: rows[j][i], query2,
                                lower2, upper2))

    return _interp1(x1, column, query1, lower1, upper1)


def tlinex_flat(x1, x2, flat, query1: float, query2: float,
                lower1: int = 0, lower2: int = 0,
                upper1: int = 0, upper2: int = 0) -> float:
    """TLINEX on a table given flat, in the source's ``Y(NX2,NX1)`` order.

    The X2 index varies fastest, so ``Y(j,i)`` is ``flat[i*len(x2) + j]``:
    each run of ``len(x2)`` values is one X1 column, the order DATA
    statements list a figure's curves in.
    """
    n1, n2 = len(x1), len(x2)
    grid = np.asarray(flat, dtype=float)
    if grid.size != n1*n2:
        raise ValueError(
            f"A flat table for {n1}x{n2} grids needs {n1*n2} values, "
            f"not {grid.size}")
    return tlinex(x1, x2, grid.reshape(n1, n2).T, query1, query2,
                  lower1, lower2, upper1, upper2)


def tlin3x(x1, x2, x3, y, query1: float, query2: float, query3: float,
           lower1: int = 0, lower2: int = 0, lower3: int = 0,
           upper1: int = 0, upper2: int = 0, upper3: int = 0) -> float:
    """TLIN3X: linear interpolation of Y = F(X1, X2, X3).

    ``Y`` must have shape ``(len(x2), len(x1), len(x3))``.  The source
    declares ``Y(NX1,NX2,NX3)`` but only ever touches it by handing slices to
    ``TLINEX``, which declares ``Y(NX2,NX1)``, so the declaration is dead.
    The effective layout was established by compiling and running the legacy
    routine over an asymmetric table, and is confirmed independently by
    ``tlin4x.f``, whose header states outright that the structure is
    ``Y(NX2,NX1,NX3,NX4)``.

    Each X3 slice the outer lookup reads is interpolated as :func:`tlinex`
    does, then those results are interpolated along X3.
    """
    x1, x2, x3 = _grid(x1), _grid(x2), _grid(x3)
    cube = _table(y, (len(x2), len(x1), len(x3)),
                  "Y must be finite with shape (len(x2), len(x1), len(x3))")
    query1, query2, query3 = _query(query1), _query(query2), _query(query3)

    def plane(k):
        def column(i):
            return _finite(_interp1(x2, lambda j: cube[j][i][k], query2,
                                    lower2, upper2))
        return _finite(_interp1(x1, column, query1, lower1, upper1))

    return _interp1(x3, plane, query3, lower3, upper3)


def tlin4x(x1, x2, x3, x4, y, query1: float, query2: float,
           query3: float, query4: float,
           lower1: int = 0, lower2: int = 0, lower3: int = 0, lower4: int = 0,
           upper1: int = 0, upper2: int = 0, upper3: int = 0,
           upper4: int = 0) -> float:
    """TLIN4X: linear interpolation of Y = F(X1, X2, X3, X4).

    ``Y`` must have shape ``(len(x2), len(x1), len(x3), len(x4))``, which the
    source header states directly as ``Y(NX2,NX1,NX3,NX4)``.

    ``INTERX`` cannot reach this routine -- its four-variable branch is
    commented out with "TLIN4X CALL DELETED TO SAVE CORE" -- but the routine
    itself is live, called directly by ``latflp.f``, ``sublat.f``,
    ``trancm.f`` and ``trhtcm.f``.

    Each X4 slice the outer lookup reads is interpolated as :func:`tlin3x`
    does, then those results are interpolated along X4.
    """
    x1, x2, x3, x4 = _grid(x1), _grid(x2), _grid(x3), _grid(x4)
    hyper = _table(
        y, (len(x2), len(x1), len(x3), len(x4)),
        "Y must be finite with shape (len(x2), len(x1), len(x3), len(x4))")
    query1, query2 = _query(query1), _query(query2)
    query3, query4 = _query(query3), _query(query4)

    def volume(m):
        def plane(k):
            def column(i):
                return _finite(_interp1(x2, lambda j: hyper[j][i][k][m],
                                        query2, lower2, upper2))
            return _finite(_interp1(x1, column, query1, lower1, upper1))
        return _finite(_interp1(x3, plane, query3, lower3, upper3))

    return _interp1(x4, volume, query4, lower4, upper4)
