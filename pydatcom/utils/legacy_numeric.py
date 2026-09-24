"""Numerical kernels translated from datcom-legacy/datcom_2000.

These expose FORTRAN routine contracts explicitly. Existing generic NumPy
wrappers are not substitutes for the legacy integration and lookup modes.
"""

import math

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


def tbfunx(x, y, query: float, lower: int = 0, upper: int = 0,
           ordered: bool = True):
    """TBFUNX: return (value, derivative) for an increasing legacy table.

    Interior values are LINEAR, while derivatives come from a local QUAD
    fit. End modes <=0 clamp the value but retain the quadratic derivative;
    mode 1 extrapolates linearly; modes >1 extrapolate quadratically.
    One point returns a constant; two points always extrapolate linearly.
    Source message/work-array arguments are omitted. Duplicate or unordered
    X is rejected instead of continuing through the source warning/sentinel.

    ``ordered=False`` accepts an unordered X, as the source silently does,
    for callers that search a curve which is not monotonic, such as WBCM
    inverting a lift curve past the stall.  Only the first and last points
    then decide between interpolation and extrapolation, and the interior
    search takes the last point at or below the query, exactly as the
    source's ``DO 1000`` loop does.  Adjacent duplicates still raise.
    """
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    if x.ndim != 1 or x.shape != y.shape or len(x) == 0:
        raise ValueError("TBFUNX requires nonempty matching one-dimensional arrays")
    if not np.all(np.isfinite(x)):
        raise ValueError("TBFUNX requires finite X coordinates")
    if ordered and not np.all(np.diff(x) > 0):
        raise ValueError("TBFUNX requires finite, strictly increasing X coordinates")
    if not ordered and np.any(np.diff(x) == 0):
        raise ValueError("TBFUNX divides by zero on adjacent duplicate X")
    n = len(x)
    if n == 1:
        return float(y[0]), 0.0
    if n == 2:
        slope = (y[1] - y[0]) / (x[1] - x[0])
        return float(y[0] + slope * (query - x[0])), float(slope)

    if x[0] < query < x[-1]:
        if ordered:
            left = max(1, int(np.searchsorted(x, query, side='right')) - 1)
        else:
            left = max(1, max(i for i in range(n - 1) if query >= x[i]))
        window = slice(left - 1, left + 2)
        # Labels 1000-1010: below XA(2) the source interpolates on the
        # window's first pair, XA(L-1) and XA(L).  For an ordered table that
        # is the table's first pair.
        segment = left - 1 if query < x[1] else left
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


# ANGLES' own DATA constants.  CON differs from the COMMON block's DEG in
# the eighth digit, so it is kept separately.
_ANGLES_CON = 1.74532925199433E-2
_ANGLES_HPI = 1.57079632679489E+0
_ANGLES_OPI = 3.14159265358979E+0
_ANGLES_TPI = 6.28318530717959E+0
_ANGLES_EPS = 1.52587890625000E-5
_ANGLES_SPE = 1.04857600000000E+6


def zerang() -> list:
    """Translate ZERANG: the zero-angle record ``[0, 0, 0, 1, 0, 0]``."""
    return [0.0, 0.0, 0.0, 1.0, 0.0, 0.0]


def _sign(a: float, b: float) -> float:
    """FORTRAN SIGN(A,B), with gfortran's treatment of a negative zero."""
    return math.copysign(abs(a), b)


def angles(entry: int, arg) -> list:
    """Translate ANGLES: complete a six-element angle record from one part.

    The record is ``[degrees, radians, sin, cos, tan, test]``, and ``entry``
    names which element was just set: 1 degrees, 2 radians, 3 sine,
    4 cosine, 5 tangent, 6 a sine-cosine pair to be normalised.  The sign
    of ``entry`` is ignored, as in the source.

    The routine is stateful through ``test``, the last angle it resolved,
    in radians.  If the new angle is within ``EPS`` (2**-16 radians) of it,
    the source returns at once and every element except the one the caller
    set keeps its previous value.  An angle within ``EPS`` of zero resets
    the whole record to zero.  Callers that keep a record across calls,
    as LIFTCF does, must therefore pass the previous record back in.

    Angles are wrapped into ``(-pi, pi]``; a cosine of zero gives a tangent
    of ``+-2**20``.  Entries 3 and 5 reconstruct a non-negative cosine, so
    they resolve into ``[-pi/2, pi/2]``.

    Returns:
        The new record as a list; ``arg`` is not modified.

    Reference: datcom-legacy/datcom_2000/angles.f
    """
    out = [float(v) for v in arg]
    eps = _ANGLES_EPS
    entry = abs(int(entry))
    if not 1 <= entry <= 6:
        raise ValueError("ANGLES entry must be 1 to 6")

    def resolve(a2: float, derive: bool) -> list:
        # Labels 1040 to 1070.
        a2 = math.fmod(a2 + _ANGLES_TPI, _ANGLES_TPI)
        if a2 > _ANGLES_OPI:
            a2 -= _ANGLES_TPI
        if abs(a2 - out[5]) <= eps:
            return out
        if abs(a2) <= eps:
            return zerang()
        out[0] = a2 / _ANGLES_CON
        out[1] = a2
        if derive:
            sine = math.sin(a2)
            out[2] = 0.0 if abs(sine) <= eps else sine
            cosine = math.cos(a2)
            if abs(cosine) > eps:
                out[3] = cosine
                out[4] = out[2] / cosine
            else:
                out[3] = 0.0
                out[4] = _sign(_ANGLES_SPE, a2)
        out[5] = a2
        return out

    def quarter(a3: float) -> list:
        # Label 1160: a cosine of zero.
        out[2] = _sign(1.0, a3)
        out[3] = 0.0
        out[4] = _sign(_ANGLES_SPE, a3)
        return resolve(_sign(_ANGLES_HPI, a3), False)

    def from_pair(a3: float, a4: float) -> list:
        # Label 1150.
        if abs(a4) <= eps:
            return quarter(a3)
        if abs(a3) <= eps:
            a3 = 0.0
        out[2] = a3
        out[3] = a4
        out[4] = a3 / a4
        return resolve(math.atan2(a3, a4), False)

    def from_sine(a3: float) -> list:
        # Label 1090.
        if abs(a3) <= eps:
            return zerang()
        if abs(a3) >= 1.0 - eps:
            return quarter(a3)
        a4 = math.sqrt(1.0 - a3**2)
        if abs(a4 - out[3]) <= eps and a3 * out[2] >= 0.0:
            return out
        return from_pair(a3, a4)

    if entry == 1:
        return resolve(out[0] * _ANGLES_CON, True)
    if entry == 2:
        return resolve(out[1], True)
    if entry == 3:
        return from_sine(out[2])
    if entry == 4:
        a4 = out[3]
        if abs(a4) >= 1.0 - eps:
            a4 = _sign(1.0, a4)
        a3 = math.sqrt(1.0 - a4**2)
        if abs(a4) <= eps:
            return quarter(a3)
        if abs(a3 - out[2]) <= eps and a4 * out[3] >= 0.0:
            return out
        return from_pair(a3, a4)
    if entry == 5:
        a1 = out[4]
        if abs(a1) > _ANGLES_SPE:
            a1 = _sign(_ANGLES_SPE, a1)
        return from_sine(_sign(a1 / math.sqrt(1.0 + a1**2), a1))
    a3, a4 = out[2], out[3]
    norm = math.sqrt(a3**2 + a4**2)
    if norm == 0.0:
        return zerang()
    return from_pair(a3 / norm, a4 / norm)
