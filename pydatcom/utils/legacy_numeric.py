"""Numerical kernels translated from datcom-legacy/datcom_2000.

These expose FORTRAN routine contracts explicitly. Existing generic NumPy
wrappers are not substitutes for the legacy integration and lookup modes.
"""

import numpy as np

from .constants import PI, UNUSED


def tranf(x, y, left_slope: float, right_slope: float,
          query: float) -> float:
    """TRANF: nonnegative piecewise-cubic transonic interpolation.

    Endpoint derivatives are supplied by the caller. Interior derivatives
    are the tangent of the mean angle of the adjacent secants, with the
    source's flattening rules at nonpositive ordinates and slope reversals.
    The first or last cubic is also used for extrapolation.
    """
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    if (x.ndim != 1 or x.shape != y.shape or len(x) < 2 or
            not np.all(np.isfinite(x)) or not np.all(np.isfinite(y)) or
            not np.all(np.diff(x) > 0.0)):
        raise ValueError("TRANF requires finite matching arrays and increasing X")
    if not all(np.isfinite(value) for value in (left_slope, right_slope, query)):
        raise ValueError("TRANF slopes and query must be finite")

    right = int(np.searchsorted(x[1:-1], query, side='right')) + 1
    left = right - 1

    def point_slope(index: int) -> float:
        if index == 0:
            return float(left_slope)
        if index == len(x) - 1:
            return float(right_slope)
        slope_left = (y[index] - y[index - 1]) / (x[index] - x[index - 1])
        slope_right = (y[index + 1] - y[index]) / (x[index + 1] - x[index])
        slope = np.tan((np.arctan(slope_left) + np.arctan(slope_right)) / 2.0)
        if y[index] <= 0.0:
            slope = 0.0
        if (len(x) > 10 and
                (abs(slope_left) <= UNUSED or abs(slope_right) <= UNUSED)):
            slope = 0.0
        if slope_right != 0.0 and slope_left / slope_right < 0.0:
            slope = 0.0
        return float(slope)

    width = x[right] - x[left]
    parameter = (query - x[left]) / width
    slope_left = point_slope(left)
    slope_right = point_slope(right)
    value = (
        (2.0 * parameter**3 - 3.0 * parameter**2 + 1.0) * y[left]
        + (parameter**3 - 2.0 * parameter**2 + parameter) * width * slope_left
        + (-2.0 * parameter**3 + 3.0 * parameter**2) * y[right]
        + (parameter**3 - parameter**2) * width * slope_right
    )
    return max(0.0, float(value))


def quad(x, y, query: float, derivative: bool = False) -> float:
    """QUAD: evaluate the parabola through three points, or its derivative.

    ``derivative=True`` replaces FORTRAN's YA='DERI' input flag. Divided
    differences are algebraically equivalent to the source determinants,
    but avoid powers of large absolute coordinates. Distinct X is required.
    """
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    if x.shape != (3,) or y.shape != (3,) or len(np.unique(x)) != 3:
        raise ValueError("QUAD requires three points with distinct X coordinates")
    first = (y[1] - y[0]) / (x[1] - x[0])
    second = ((y[2] - y[1]) / (x[2] - x[1]) - first) / (x[2] - x[0])
    if derivative:
        return float(first + second * ((query - x[0]) + (query - x[1])))
    return float(y[0] + (query - x[0]) * (first + (query - x[1]) * second))


def trapz(f, y, ntest: int = 1) -> np.ndarray:
    """TRAPZ: total/cumulative trapezoids or volume of circular frustums.

    ``ntest=1`` returns a one-element total integral; other nonnegative
    values return cumulative integrals starting at zero. Negative NTEST
    interprets F as radii and returns a one-element volume using the exact
    conical-frustum rule, as in the source. Coordinate direction is retained.
    Unlike the FORTRAN work array, only initialized results are returned.
    """
    f = np.asarray(f, dtype=float)
    y = np.asarray(y, dtype=float)
    if f.ndim != 1 or f.shape != y.shape or len(f) == 0:
        raise ValueError("TRAPZ requires nonempty matching one-dimensional arrays")
    width = np.diff(y)
    if ntest < 0:
        increments = PI / 3 * width * (f[1:]**2 + f[:-1]**2 + f[1:]*f[:-1])
    else:
        increments = .5 * width * (f[1:] + f[:-1])
    cumulative = np.concatenate(([0.], np.cumsum(increments)))
    return cumulative[-1:] if ntest == 1 or ntest < 0 else cumulative


def tbfunx(x, y, query: float, lower: int = 0, upper: int = 0):
    """TBFUNX: return (value, derivative) for an increasing legacy table.

    Interior values are LINEAR, while derivatives come from a local QUAD
    fit. End modes <=0 clamp the value but retain the quadratic derivative;
    mode 1 extrapolates linearly; modes >1 extrapolate quadratically.
    One point returns a constant; two points always extrapolate linearly.
    Source message/work-array arguments are omitted. Duplicate or unordered
    X is rejected instead of continuing through the source warning/sentinel.
    """
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    if x.ndim != 1 or x.shape != y.shape or len(x) == 0:
        raise ValueError("TBFUNX requires nonempty matching one-dimensional arrays")
    if not np.all(np.isfinite(x)) or not np.all(np.diff(x) > 0):
        raise ValueError("TBFUNX requires finite, strictly increasing X coordinates")
    n = len(x)
    if n == 1:
        return float(y[0]), 0.0
    if n == 2:
        slope = (y[1] - y[0]) / (x[1] - x[0])
        return float(y[0] + slope * (query - x[0])), float(slope)

    if x[0] < query < x[-1]:
        left = max(1, int(np.searchsorted(x, query, side='right')) - 1)
        window = slice(left - 1, left + 2)
        # Labels 1000-1010: leftmost interval uses the first pair.
        segment = 0 if query < x[1] else left
        slope = (y[segment+1] - y[segment]) / (x[segment+1] - x[segment])
        value = y[segment] + slope * (query - x[segment])
        return float(value), quad(x[window], y[window], query, derivative=True)

    at_upper = query >= x[-1]
    start = n - 3 if at_upper else 0
    endpoint = n - 1 if at_upper else 0
    mode = upper if at_upper else lower
    if mode == 1:
        # Preserve labels 1020-1040, including NP=3: LE=NP3=0
        # makes LEXU=1 override the lower-end slope in that special case.
        segment = n - 2 if upper == 1 and start == n - 3 else 0
        slope = (y[segment+1] - y[segment]) / (x[segment+1] - x[segment])
        anchor = n - 1 if segment == n - 2 else 0
        return float(y[anchor] + slope * (query - x[anchor])), float(slope)
    window = slice(start, start + 3)
    value = y[endpoint] if mode <= 0 else quad(x[window], y[window], query)
    return float(value), quad(x[window], y[window], query, derivative=True)


def arccos(value: float) -> float:
    """Translate ARCCOS, the source's "standard FORTRAN only" inverse cosine.

    Inside ``[-1, 1]`` this is the ordinary inverse cosine, built from
    ``atan(sqrt(1-a^2)/a)`` with ``pi`` added when that lands negative, and
    ``pi/2`` at exactly zero.

    Outside ``[-1, 1]`` it is *not* an inverse cosine at all: the source
    returns ``log|a + sqrt(a^2-1)|``, the inverse hyperbolic cosine.  For
    ``a > 1`` the true ``arccos`` is imaginary and this is its magnitude, so
    the continuation is deliberate.  For ``a < -1`` the same expression is
    evaluated on a quantity that shrinks toward zero, so the result turns
    negative; that is the source's behaviour and is preserved.

    At exactly ``a = -1`` the routine returns zero where the true inverse
    cosine is ``pi``.  ``sqrt(1-a^2)/a`` is negative zero there, and the
    source gates its ``pi`` correction on ``X .LT. 0.0``, which negative
    zero does not satisfy.  Every other argument in ``(-1, 0)`` picks the
    correction up normally, so the defect is confined to that one point.
    It is reproduced rather than corrected.

    Reference: datcom-legacy/datcom_2000/arccos.f
    """
    value = float(value)
    if not np.isfinite(value):
        raise ValueError("ARCCOS requires a finite argument")
    if value == 0.0:
        return float(np.pi / 2.0)
    if abs(value) <= 1.0:
        angle = float(np.arctan(np.sqrt(1.0 - value**2) / value))
        return angle if angle >= 0.0 else float(np.pi + angle)
    return float(np.log(abs(value + np.sqrt(value**2 - 1.0))))
