"""Numerical kernels translated from datcom-legacy/datcom_2000.

These expose FORTRAN routine contracts explicitly. Existing generic NumPy
wrappers are not substitutes for the legacy integration and lookup modes.
"""

import math

import numpy as np

from .constants import PI, RAD, UNUSED


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
    if derivative:
        first = (y[1] - y[0]) / (x[1] - x[0])
        second = ((y[2] - y[1]) / (x[2] - x[1]) - first) / (x[2] - x[0])
        return float(first + second * ((query - x[0]) + (query - x[1])))
    return float(_parabola(x.tolist(), y.tolist(), query))


def _parabola(x, y, query: float) -> float:
    """QUAD's value on three plain-float points, for the table kernels.

    Raises:
        ValueError: For repeated X, as :func:`quad` does.
    """
    x0, x1, x2 = x
    y0, y1, y2 = y
    if x0 == x1 or x1 == x2 or x0 == x2:
        raise ValueError("QUAD requires three points with distinct X coordinates")
    first = (y1 - y0) / (x1 - x0)
    second = ((y2 - y1) / (x2 - x1) - first) / (x2 - x0)
    return y0 + (query - x0) * (first + (query - x1) * second)


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
            left = max(1, max(bracket_index for bracket_index in range(n - 1)
                              if query >= x[bracket_index]))
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


def sleq(a, b):
    """Translate SLEQ: Gauss-Jordan solution of ``A x = b`` without pivoting.

    On a zero pivot the source rotates every row down by one (the last row
    to the top, on the partly reduced matrix), counts the rotation, and
    restarts the elimination from the first row; after ``N`` rotations it
    prints a warning and returns without setting ``X``.  Both are kept.

    Args:
        a: ``N`` by ``N`` coefficients (the source's ``A(N,M)`` with
            ``M >= N+1``; the right-hand side goes in column ``N+1``).
        b: Right-hand side.

    Returns:
        ``(x, ok)``: the solution, or ``None`` with ``ok`` False after the
        warning.

    Reference: datcom-legacy/datcom_2000/sleq.f
    """
    b = np.asarray(b, dtype=float)
    n = len(b)
    work = np.zeros((n, n + 1))
    work[:, :n] = np.asarray(a, dtype=float)
    work[:, n] = b
    rotations = 0
    while True:
        for pivot_row in range(n):
            if work[pivot_row, pivot_row] == 0.0:
                break
            work[pivot_row, pivot_row + 1:] = (
                work[pivot_row, pivot_row + 1:] / work[pivot_row, pivot_row])
            work[pivot_row, pivot_row] = 1.0
            for row_index in range(n):
                if row_index == pivot_row:
                    continue
                work[row_index, pivot_row + 1:] = (
                    work[row_index, pivot_row + 1:]
                    - work[row_index, pivot_row]
                    * work[pivot_row, pivot_row + 1:])
                work[row_index, pivot_row] = 0.0
        else:
            return work[:, n].copy(), True
        rotations += 1
        if rotations > n:
            return None, False
        work = np.roll(work, 1, axis=0)


def quadin(y, h: float) -> float:
    """Translate QUADIN: integrate equally spaced ordinates.

    Five-point Newton-Cotes (Boole) panels from the start, then the one to
    four points left over by the trapezoid, Simpson or three-eighths rule
    over the *last* points.  Nonpositive ``h`` or no points give zero.

    Reference: datcom-legacy/datcom_2000/quadin.f
    """
    y = np.asarray(y, dtype=float)
    n = len(y)
    if h <= 0.0 or n <= 0:
        return 0.0
    ans = 0.0
    k = n
    if n >= 4:
        for panel_start in range(1, n + 1, 4):
            k = n - panel_start + 1
            if k >= 5:
                w = y[panel_start - 1:panel_start + 4]
                ans += 7.*w[0] + 32.*w[1] + 12.*w[2] + 32.*w[3] + 7.*w[4]
        ans = ans * h / 22.5
    if k == 2:
        ans += h * (y[n - 1] + y[n - 2]) / 2.
    elif k == 3:
        ans += h * (y[n - 3] + 4. * y[n - 2] + y[n - 1]) / 3.
    elif k == 4:
        ans += 3. * h * (y[n - 4] + 3. * y[n - 3] + 3. * y[n - 2] +
                         y[n - 1]) / 8.
    return float(ans)


def mach2(nu_deg: float):
    """Translate MACH2: the Mach number of a Prandtl-Meyer angle.

    Newton iteration on ``sqrt(6) atan(b/sqrt(6)) - atan(b) = nu`` from
    ``b = 1``, to 2e-6, at most 30 steps.  Returns ``(mach, ier)``: ``ier``
    1 for a negative angle (Mach 1), 2 above 130 degrees (Mach 1000), 3
    when the iteration does not converge (the last iterate is used).

    Reference: datcom-legacy/datcom_2000/mach2.f
    """
    if nu_deg < 0.0:
        return 1.0, 1
    if nu_deg == 0.0:
        return 1.0, 0
    if nu_deg > 130.0:
        return 1000.0, 2
    sqr6 = 2.44948974
    rnu = nu_deg / RAD
    b = 1.0
    ier = 3
    for _ in range(30):
        f = sqr6 * math.atan(b / sqr6) - math.atan(b) - rnu
        if abs(f) <= 2.0e-6:
            ier = 0
            break
        b = b - f * (1. + b * b) * (6. + b * b) / (5. * b * b)
    return math.sqrt(b * b + 1.0), ier


def simul2(x, c1, c2):
    """Translate SIMUL2: where two tabulated curves cross.

    Scans for the first exact touch or sign change of ``c2 - c1``, then
    halves the step from the left bracket, moving the bracket whenever the
    difference keeps the left point's sign, until the difference is under
    0.1% of ``c1`` (of ``c2`` where ``c1`` is zero), exactly zero, or 101
    halvings have passed.  Both curves are read with TBFUNX.  Returns
    ``(x, c1)`` at the crossing, or ``(-1000, -1000)`` when they do not
    cross.

    Reference: datcom-legacy/datcom_2000/simul2.f
    """
    x = [float(v) for v in x]
    c1 = [float(v) for v in c1]
    c2 = [float(v) for v in c2]
    signp = None
    for point_index in range(len(x)):
        cross = c2[point_index] - c1[point_index]
        if cross == 0.0:
            return x[point_index], c1[point_index]
        sign = math.copysign(1.0, cross)
        if point_index > 0 and sign != signp:
            break
        signp = sign
    else:
        return -1000.0, -1000.0
    xd2 = x[point_index] - x[point_index - 1]
    xsrt = x[point_index - 1]
    kount = 0
    while True:
        xd2 = xd2 / 2.
        xtst = xsrt + xd2
        cd1 = tbfunx(x, c1, xtst, 0, 0)[0]
        cd2 = tbfunx(x, c2, xtst, 0, 0)[0]
        cr = cd2 - cd1
        if cr == 0.0:
            return xtst, cd1
        cdv = cr / cd2 if cd1 == 0.0 else cr / cd1
        if abs(cdv) < 0.001 or kount > 100:
            return xtst, cd1
        kount += 1
        if math.copysign(1.0, cr) == signp:
            xsrt = xtst


def tlinvs(x1, x2, y, xa2: float, za: float) -> float:
    """Translate TLINVS: the ``X1`` at which TLINEX's surface reaches ``za``.

    ``y`` has TLINEX's ``(len(x2), len(x1))`` layout, and is taken to fall
    with ``X1``: at the ``X2`` bracket, a value at or above the first
    column returns ``X1(1)`` and one at or below the last returns
    ``X1(NX1)``.  Otherwise a relative-step search from ``X1(NX1/2)``
    runs TLINEX until within 1e-3 of ``za`` or for 21 evaluations, clamped
    to the grid.

    Raises:
        ValueError: Where the source would read past ``Y``: ``xa2`` beyond
            the last ``X2`` with ``za`` between that row's end values.

    Reference: datcom-legacy/datcom_2000/tlinvs.f
    """
    from .legacy_tables import tlinex
    x1 = np.asarray(x1, dtype=float)
    x2 = np.asarray(x2, dtype=float)
    y = np.asarray(y, dtype=float)
    nx1, nx2 = len(x1), len(x2)

    def search():
        dgss = x1[nx1 // 2 - 1]
        kount = 0
        while True:
            zgss = float(tlinex(x1, x2, y, dgss, xa2, 0, 0, 0, 0))
            kount += 1
            if abs(zgss - za) < 1.e-3 or kount > 20:
                return float(dgss)
            din = dgss * abs(zgss - za) / zgss
            if zgss > za:
                dgss = min(dgss + din, x1[-1])
            else:
                dgss = max(dgss - din, x1[0])

    def row(row_index):              # labels 1050-1070 at an exact row
        if not za < y[row_index, 0]:
            return float(x1[0])
        if not za > y[row_index, -1]:
            return float(x1[-1])
        return search()

    if not xa2 > x2[0]:
        if not za < y[0, 0]:
            return float(x1[0])
        if not za > y[0, -1]:
            return float(x1[-1])
    if not xa2 < x2[-1]:
        if not za < y[-1, 0]:
            return float(x1[0])
        if not za > y[-1, -1]:
            return float(x1[-1])
    for row_index in range(1, nx2):
        if x2[row_index] == xa2:
            return row(row_index)
        if x2[row_index] > xa2:
            rat = ((xa2 - x2[row_index - 1])
                   / (x2[row_index] - x2[row_index - 1]))
            if za >= (y[row_index - 1, 0]
                      + (y[row_index, 0] - y[row_index - 1, 0]) * rat):
                return float(x1[0])
            if za <= (y[row_index - 1, -1]
                      + (y[row_index, -1] - y[row_index - 1, -1]) * rat):
                return float(x1[-1])
            return search()
    raise ValueError("TLINVS would read past its table beyond the last X2")


def inter3(arg1: float, arg2: float, rl: float, tables) -> float:
    """Translate INTER3: TLINEX in five tables bracketed by Reynolds number.

    ``tables`` holds five ``(x1, x2, y)`` triples, ``y`` in TLINEX's
    ``(len(x2), len(x1))`` layout, for Reynolds numbers 1e5, 1e6, 1e7, 1e8
    and 1e9.  The two tables bracketing ``rl`` are read and interpolated
    *linearly in Reynolds number*; below 1e5 the first table is used (its
    two reads are identical), above 1e9 the last.

    Reference: datcom-legacy/datcom_2000/inter3.f
    """
    from .legacy_tables import tlinex

    def read(table_index):
        x1, x2, y = tables[table_index]
        return float(tlinex(x1, x2, np.asarray(y, dtype=float), arg1, arg2,
                            0, 0, 0, 0))

    decades = [1e5, 1e6, 1e7, 1e8, 1e9]
    if rl > 1e9:
        return read(4)
    if rl <= 1e5:
        return read(0)
    it = next(decade_index for decade_index in range(1, 5)
              if rl <= decades[decade_index])
    low, high = read(it - 1), read(it)
    x1, x2 = decades[it - 1], decades[it]
    return low + (high - low) * (rl - x1) / (x2 - x1)


def simul4(coff, eq):
    """Translate SIMUL4: four simultaneous equations by Cramer's rule.

    ``coff`` is the source's 16-word ``COFF``, equation ``i``'s
    coefficient on unknown ``m`` at ``COFF(4(i-1)+m)``; ``eq`` is the
    right-hand side.  Each unknown is DET4 with words ``m, m+4, m+8,
    m+12`` replaced by ``eq``, over DET4 of ``coff``.  A singular system divides by zero as
    the source does.

    Reference: datcom-legacy/datcom_2000/simul4.f
    """
    from .math_utils import det4
    coff = np.asarray(coff, dtype=float).reshape(-1)
    d = det4(coff)
    unk = []
    for unknown_index in range(4):
        de = coff.copy()
        de[unknown_index::4] = np.asarray(eq, dtype=float)
        g = det4(de)
        unk.append(g / d if d != 0.0 else
                   (math.nan if g == 0.0 else math.copysign(math.inf, g)))
    return unk
