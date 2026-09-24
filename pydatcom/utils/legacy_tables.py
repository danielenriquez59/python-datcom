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
"""

from typing import NamedTuple, Optional

import numpy as np

from .legacy_numeric import quad


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


def _grid(x):
    x = np.asarray(x, dtype=float)
    if x.ndim != 1 or not len(x) or not np.all(np.isfinite(x)):
        raise ValueError("A grid must be a nonempty finite one-dimensional array")
    if len(x) > 1 and not (np.all(np.diff(x) > 0) or np.all(np.diff(x) < 0)):
        raise ValueError("Grid coordinates must be strictly monotonic")
    return x


def glook(x, query: float, ascending: Optional[bool] = None) -> Lookup:
    """Translate GLOOK with its caller-initialized NOING=False contract.

    A None fraction means use Y[index] directly (snap or endpoint clamp).
    Otherwise interpolate between index-1 and index with the given fraction.

    A NaN query fails every comparison, so the source's search runs off the
    end and returns the last node, exactly; that is kept.  Infinite queries
    are rejected.
    """
    x = _grid(x)
    if np.isinf(query):
        raise ValueError("The query must be finite")
    if ascending is None:
        ascending = x[0] <= x[-1]
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
            return Lookup(index, float(previous / (previous - difference)))
        previous = difference
    return Lookup(len(x)-1, None)


def switch(x, query: float, lower: int = 0, upper: int = 0) -> Switches:
    """Translate SWITCH's seven LG flags, preserving table-end semantics."""
    x = _grid(x)
    if np.isinf(query):
        raise ValueError("The query must be finite")
    ascending = bool(x[0] <= x[-1])
    after = bool(query > x[-1] if ascending else query < x[-1])
    before = bool(query < x[0] if ascending else query > x[0])
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
    y = np.asarray(y, dtype=float)
    if y.shape != x.shape or not np.all(np.isfinite(y)):
        raise ValueError("Y must be a finite array with the same shape as X")
    flags = switch(x, query, lower, upper)
    if not flags.use_extrapolation:
        index, fraction = glook(x, query, flags.ascending)
        if fraction is None:
            return float(y[index])
        return float(y[index-1] + fraction*(y[index]-y[index-1]))

    if len(x) < 2:
        raise ValueError("Extrapolation requires at least two grid points")
    mode = lower if flags.before_first else upper
    if mode > 1 and len(x) > 2:
        selected = slice(0, 3) if flags.before_first else slice(-3, None)
        return quad(x[selected], y[selected], query)
    if flags.before_first:
        fraction = (query-x[0])/(x[1]-x[0])
        return float(y[0] + fraction*(y[1]-y[0]))
    fraction = (query-x[-1])/(x[-1]-x[-2])
    return float(y[-1] + fraction*(y[-1]-y[-2]))


def tlinex(x1, x2, y, query1: float, query2: float,
           lower1: int = 0, lower2: int = 0,
           upper1: int = 0, upper2: int = 0) -> float:
    """TLINEX: interpolate X2 columns, then interpolate their values in X1.

    Y must have shape ``(len(x2), len(x1))``, matching Fortran Y(NX2,NX1).
    Each axis uses TLIN1X's snapping and independent extrapolation modes.
    The legacy label sequence evaluates only needed X1 columns; this version
    evaluates all columns before the same outer linear/quadratic operation.
    The numerical result is equivalent for finite valid tables.
    """
    x1, x2 = _grid(x1), _grid(x2)
    y = np.asarray(y, dtype=float)
    if y.shape != (len(x2), len(x1)) or not np.all(np.isfinite(y)):
        raise ValueError("Y must be finite with shape (len(x2), len(x1))")
    columns = [tlin1x(x2, y[:, index], query2, lower2, upper2)
               for index in range(len(x1))]
    return tlin1x(x1, columns, query1, lower1, upper1)


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

    Each X3 slice is interpolated by :func:`tlinex`, then those results are
    interpolated along X3 by :func:`tlin1x`.  The legacy label sequence
    evaluates only the needed slices; this evaluates all of them before the
    same outer operation, which is numerically equivalent for finite tables.
    """
    x1, x2, x3 = _grid(x1), _grid(x2), _grid(x3)
    y = np.asarray(y, dtype=float)
    if y.shape != (len(x2), len(x1), len(x3)) or not np.all(np.isfinite(y)):
        raise ValueError(
            "Y must be finite with shape (len(x2), len(x1), len(x3))")
    slices = [tlinex(x1, x2, y[:, :, k], query1, query2,
                     lower1, lower2, upper1, upper2)
              for k in range(len(x3))]
    return tlin1x(x3, slices, query3, lower3, upper3)


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

    Each X4 slice is interpolated by :func:`tlin3x`, then those results are
    interpolated along X4.
    """
    x1, x2, x3, x4 = _grid(x1), _grid(x2), _grid(x3), _grid(x4)
    y = np.asarray(y, dtype=float)
    expected = (len(x2), len(x1), len(x3), len(x4))
    if y.shape != expected or not np.all(np.isfinite(y)):
        raise ValueError(
            "Y must be finite with shape (len(x2), len(x1), len(x3), len(x4))")
    slices = [tlin3x(x1, x2, x3, y[:, :, :, l], query1, query2, query3,
                     lower1, lower2, lower3, upper1, upper2, upper3)
              for l in range(len(x4))]
    return tlin1x(x4, slices, query4, lower4, upper4)
