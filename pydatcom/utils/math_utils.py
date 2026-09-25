"""
Mathematical utility functions for PyDATCOM.

Ported from DATCOM FORTRAN subroutines with bounds checking and error handling.
References: datcom.f lines 390-451 (ARCSIN), 452-515 (AREA1, AREA2).
ARCCOS and SIGN are in ``legacy_numeric``.
"""

import math
import numpy as np
from typing import Tuple
import logging

from pydatcom.utils.constants import PI

logger = logging.getLogger(__name__)


def arcsin(a: float) -> float:
    """
    Translate ARCSIN: arc sine in radians.

    Exactly +-1 returns +-pi/2 with the source's ``PI``.  Outside [-1, 1]
    the source prints an error and returns 1000; that sentinel is kept,
    with a logged error, because callers test for it.

    Reference: datcom-legacy/datcom_2000/arcsin.f
    """
    if abs(a) == 1.0:
        return PI / 2.0 * a / abs(a)
    if abs(a) > 1.0:
        logger.error(f"ARCSIN of {a:.5e} is out of bounds")
        return 1000.0
    return math.atan(a / math.sqrt(1.0 - a**2))


def area1(x: np.ndarray, y: np.ndarray, nsum: int) -> float:
    """
    Calculate incremental area of vertical tail shadowed by Mach line.
    
    Uses Heron's formula for triangle area.
    Reference: FORTRAN AREA1, datcom.f line 452
    
    Args:
        x: Array of x-coordinates (length 4)
        y: Array of y-coordinates (length 4)
        nsum: Number of points to use (3, 4, or 6)
        
    Returns:
        Computed area
    """
    # Calculate triangle area using first 3 points
    a = math.sqrt((x[1] - x[0])**2 + (y[1] - y[0])**2)
    b = math.sqrt((x[2] - x[1])**2 + (y[2] - y[1])**2)
    c = math.sqrt((x[0] - x[2])**2 + (y[0] - y[2])**2)
    s = (a + b + c) / 2.0
    area = math.sqrt(s * (s - a) * (s - b) * (s - c))
    
    # Add second triangle if nsum is 4 or 6
    if nsum == 4 or nsum == 6:
        a = math.sqrt((x[3] - x[0])**2 + (y[3] - y[0])**2)
        b = math.sqrt((x[2] - x[3])**2 + (y[2] - y[3])**2)
        c = math.sqrt((x[0] - x[2])**2 + (y[0] - y[2])**2)
        s = (a + b + c) / 2.0
        area2 = math.sqrt(s * (s - a) * (s - b) * (s - c))
        area += area2
    
    return area


def area2(x: np.ndarray, y: np.ndarray, inum: int) -> Tuple[float, float, float]:
    """
    Calculate body-shadow area and first moments as in FORTRAN AREA2.

    Reference: datcom-legacy/datcom_2000/area2.f and its BDAREA callers.
    Triangle (1, 2, 3) is always included; INUM=2 adds (1, 3, 4).
    Areas are unsigned, including for the lower body profile. This is a
    sum of triangle areas, not a signed polygon-area calculation.

    BDAREA's quadrilaterals have Y(1)=Y(4)=0. AREA2's second-triangle
    centroid formula assumes this geometry, which is required here too.
    Determinants replace Heron's formula and vertex averages replace
    intersecting medians, avoiding singular slopes for vertical medians
    and giving zero moments for collapsed triangles.

    Args:
        x: One-dimensional x-coordinates, at least 3 or 4 entries.
        y: Matching y-coordinates.
        inum: Legacy selector: 3 for a triangle, 2 for a quadrilateral.

    Returns:
        (area, AX, AY), where AX=integral(x dA) and AY=integral(y dA).
        These are first moments, not centroid coordinates. Divide each
        moment by nonzero area to obtain the centroid.

    Raises:
        ValueError: For an unsupported selector or body-shadow geometry.
    """
    if inum not in (2, 3):
        raise ValueError("AREA2 inum must be 3 (triangle) or 2 (quadrilateral)")
    count = 3 if inum == 3 else 4
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    if x.ndim != 1 or y.ndim != 1 or len(x) < count or len(y) < count:
        raise ValueError(f"AREA2 requires at least {count} x and y coordinates")
    if not np.all(np.isfinite(x[:count])) or not np.all(np.isfinite(y[:count])):
        raise ValueError("AREA2 coordinates must be finite")
    if inum == 2 and (y[0] != 0.0 or y[3] != 0.0):
        raise ValueError("AREA2 quadrilateral requires y[0] = y[3] = 0")

    area = ax = ay = 0.0
    triangles = ((0, 1, 2),) if inum == 3 else ((0, 1, 2), (0, 2, 3))
    for v0, v1, v2 in triangles:
        da = abs((x[v1] - x[v0]) * (y[v2] - y[v0])
                 - (x[v2] - x[v0]) * (y[v1] - y[v0])) / 2.0
        area += da
        ax += da * (x[v0] + x[v1] + x[v2]) / 3.0
        ay += da * (y[v0] + y[v1] + y[v2]) / 3.0
    return float(area), float(ax), float(ay)


def det4(a: np.ndarray) -> float:
    """
    Translate DET4: a 4x4 determinant by cofactor expansion.

    ``a`` is the source's 16-word array (column-major ``A(4,4)``); a 4x4
    array is taken in the same element order.  The expansion runs along
    ``A(1..4)`` with 3x3 minors formed by skipping every fourth word, as
    the source does, rather than through a factorisation.  An exactly
    singular matrix gives exactly zero, which SIMUL4 tests for; an LU
    factorisation (``np.linalg.det``) leaves round-off there instead.

    Reference: datcom-legacy/datcom_2000/det4.f
    """
    flat = np.asarray(a, dtype=float)
    if flat.ndim == 2:
        flat = flat.reshape(-1, order='F')
    p = 0.0
    for cofactor_row in range(1, 5):
        a3 = [flat[word - 1] for word in range(5, 17)
              if (word - cofactor_row) % 4 != 0]
        pp = (a3[0] * (a3[4] * a3[8] - a3[5] * a3[7]) -
              a3[1] * (a3[3] * a3[8] - a3[5] * a3[6]) +
              a3[2] * (a3[3] * a3[7] - a3[4] * a3[6]))
        if cofactor_row in (2, 4):
            pp = -pp
        p += flat[cofactor_row - 1] * pp
    return float(p)
