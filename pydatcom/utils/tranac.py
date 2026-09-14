"""
TRANAC: transonic planform CL by non-linear interpolation.

A sibling of ``TRANF``.  Both define interior slopes as the tangent of the
mean angle of the adjacent secants and fit a cubic per interval, so the
curve and its first derivative are continuous, with the end slopes supplied
by the caller.

The two differ in what surrounds that core.  ``TRANF`` clamps its output
non-negative and flattens the slope at non-positive ordinates.  ``TRANAC``
instead carries a ``DELY`` offset that is applied to the fifth point of an
eight-point table, and solves the cubic coefficients explicitly rather than
through a Hermite form.

Reference: datcom-legacy/datcom_2000/tranac.f
"""

import numpy as np
from typing import Dict, Sequence
import logging

from pydatcom.utils.constants import UNUSED

logger = logging.getLogger(__name__)

# The DELY mechanism applies only to this table shape and point.
_DELY_POINTS = 8
_DELY_INDEX = 5          # one-based, as the source writes it


def tranac(x: Sequence[float], y: Sequence[float],
           slope_left: float, slope_right: float,
           query: float, dely: float = 0.0) -> Dict[str, float]:
    """Translate TRANAC: cubic interpolation with angular-average slopes.

    Args:
        x: Abscissas, strictly increasing.
        y: Ordinates.
        slope_left: ``DYL``, the slope imposed at the first point.
        slope_right: ``DYR``, the slope imposed at the last point.
        query: ``XVAL``, where to interpolate.
        dely: ``DELY``, an offset applied to the fifth ordinate of an
            eight-point table only.

    Returns:
        Dictionary with ``value`` and the cubic coefficients and interval
        used.

    Raises:
        ValueError: If the arrays are mismatched, shorter than two points,
            not increasing, or the cubic solve is degenerate.

    Notes:
        The source mutates the caller's ``Y`` array, adding ``DELY`` before
        the solve and subtracting it afterwards.  This translation copies
        instead, so the caller's data is never touched; the numerical result
        is identical because the source restores the value before returning.

        A negative ``DELY`` also forces the slope at the fifth point of an
        eight-point table to zero, independently of the offset itself.  The
        source tests ``DELY .LT. -UNUSED``, which for its 1e-30 sentinel is
        simply "negative".
    """
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float).copy()
    count = x.size
    if x.ndim != 1 or x.shape != y.shape or count < 2:
        raise ValueError("TRANAC needs matching arrays of at least two points")
    if not np.all(np.diff(x) > 0.0):
        raise ValueError("TRANAC requires strictly increasing abscissas")

    # Bracket: the source scans interior points, falling through to the last.
    bracket = count
    for position in range(2, count):        # one-based J = 2 .. NPT-1
        if query < x[position - 1]:
            bracket = position
            break
    left = bracket - 1                      # one-based J after J = J-1
    right = left + 1                        # one-based K

    dely_active = count == _DELY_POINTS
    slopes = []
    for offset in range(2):
        point = left + offset               # one-based MP
        if point == 1:
            slopes.append(float(slope_left))
            continue
        if point == count:
            slopes.append(float(slope_right))
            continue
        below = y[point - 2]
        here = y[point - 1]
        above = y[point]
        secant_left = (below - here) / (x[point - 2] - x[point - 1])
        secant_right = (here - above) / (x[point - 1] - x[point])
        mean_angle = (np.arctan(secant_left) + np.arctan(secant_right)) / 2.0
        slope = np.sin(mean_angle) / np.cos(mean_angle)
        # A negative DELY flattens the slope at the fifth point.
        if dely_active and point == _DELY_INDEX and dely < -UNUSED:
            slope = 0.0
        # A sign reversal between the two secants also flattens it.
        if secant_right != 0.0 and secant_left / secant_right < 0.0:
            slope = 0.0
        slopes.append(float(slope))

    # The offset applies only while the cubic is being solved.
    work = y.copy()
    if dely_active:
        for point in (left, right):
            if point == _DELY_INDEX:
                work[point - 1] += dely

    xj, xk = x[left - 1], x[right - 1]
    yj, yk = work[left - 1], work[right - 1]
    xj2, xk2 = xj * xj, xk * xk
    dx = xj - xk
    dy = yj - yk
    dx2 = xj2 - xk2
    dx3 = xj2 * xj - xk2 * xk
    dslope = slopes[0] - slopes[1]

    red = 2.0 * xk * dx - dx2
    green = 3.0 * xk2 * dx - dx3
    yellow = slopes[1] * dx - dy
    denominator = 3.0 * dx2 * red - 2.0 * dx * green
    if denominator == 0.0:
        raise ValueError("TRANAC's cubic solve is degenerate on this interval")

    a = (dslope * red - 2.0 * dx * yellow) / denominator
    b = (3.0 * dx2 * yellow - dslope * green) / denominator
    c = (dy - a * dx3 - b * dx2) / dx
    d = yk - a * xk2 * xk - b * xk2 - c * xk

    value = ((a * query + b) * query + c) * query + d
    return {
        'value': float(value),
        'a': float(a), 'b': float(b), 'c': float(c), 'd': float(d),
        'interval': (left, right),
        'slopes': tuple(slopes),
        'dely_applied': bool(dely_active and dely != 0.0),
        'method': 'legacy_tranac',
    }
