"""
Mathematical utility functions for PyDATCOM.

Ported from DATCOM FORTRAN subroutines with bounds checking and error handling.
References: datcom.f lines 390-451 (ARCSIN, ARCCOS), 452-515 (AREA1, AREA2)
"""

import numpy as np
from typing import Tuple, List
import logging

from pydatcom.utils.constants import PI

logger = logging.getLogger(__name__)


def arcsin(a: float) -> float:
    """
    Arc sine with bounds checking.
    
    Handles values outside [-1, 1] gracefully.
    Reference: FORTRAN ARCSIN function, datcom.f line 434
    
    Args:
        a: Input value
        
    Returns:
        Arc sine in radians
        
    Raises:
        ValueError: If abs(a) > 1.0
    """
    # Handle exact ±1 case
    if abs(a) == 1.0:
        return PI / 2.0 * a / abs(a)
    
    # Check bounds
    if abs(a) > 1.0:
        logger.error(f"ARCSIN of {a:.5e} is out of bounds")
        raise ValueError(f"arcsin argument {a} is out of range [-1, 1]")
    
    # Standard calculation
    return np.arctan(a / np.sqrt(1.0 - a**2))


def arccos(a: float) -> float:
    """
    Arc cosine with extended domain handling.
    
    Uses inverse cosh for |a| > 1.0 to provide smooth extension.
    Reference: FORTRAN ARCCOS function, datcom.f line 390
    
    Args:
        a: Input value
        
    Returns:
        Arc cosine in radians (or extended value for |a| > 1)
    """
    # Handle a = 0 case
    if a == 0.0:
        return PI / 2.0
    
    # Handle |a| > 1 using inverse cosh
    if abs(a) > 1.0:
        x = np.log(abs(a + np.sqrt(a**2 - 1.0)))
        return x
    
    # Standard calculation
    x = np.arctan(np.sqrt(1.0 - a**2) / a)
    
    # Adjust for negative values
    if x < 0.0:
        x = PI + x
    
    return x


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
    a = np.sqrt((x[1] - x[0])**2 + (y[1] - y[0])**2)
    b = np.sqrt((x[2] - x[1])**2 + (y[2] - y[1])**2)
    c = np.sqrt((x[0] - x[2])**2 + (y[0] - y[2])**2)
    s = (a + b + c) / 2.0
    area = np.sqrt(s * (s - a) * (s - b) * (s - c))
    
    # Add second triangle if nsum is 4 or 6
    if nsum == 4 or nsum == 6:
        a = np.sqrt((x[3] - x[0])**2 + (y[3] - y[0])**2)
        b = np.sqrt((x[2] - x[3])**2 + (y[2] - y[3])**2)
        c = np.sqrt((x[0] - x[2])**2 + (y[0] - y[2])**2)
        s = (a + b + c) / 2.0
        area2 = np.sqrt(s * (s - a) * (s - b) * (s - c))
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
    for i, j, k in triangles:
        da = abs((x[j] - x[i]) * (y[k] - y[i])
                 - (x[k] - x[i]) * (y[j] - y[i])) / 2.0
        area += da
        ax += da * (x[i] + x[j] + x[k]) / 3.0
        ay += da * (y[i] + y[j] + y[k]) / 3.0
    return float(area), float(ax), float(ay)


def det4(a: np.ndarray) -> float:
    """
    Calculate determinant of 4x4 matrix.
    
    Reference: FORTRAN DET4 subroutine
    
    Args:
        a: 4x4 matrix as 1D array (16 elements, row-major)
        
    Returns:
        Determinant value
    """
    # Reshape to 2D if needed
    if a.shape == (16,):
        matrix = a.reshape((4, 4))
    else:
        matrix = a
    
    return np.linalg.det(matrix)


def solve_linear(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    """
    Solve linear system Ax = b.
    
    Args:
        a: Coefficient matrix
        b: Right-hand side vector
        
    Returns:
        Solution vector x
    """
    return np.linalg.solve(a, b)


# DEPRECATED: Use NumPy directly
# These wrappers kept for backward compatibility only

def trapz_integrate(x: np.ndarray, y: np.ndarray) -> float:
    """
    DEPRECATED: Use np.trapz() or np.trapezoid() directly.
    
    Trapezoidal integration.
    
    Args:
        x: Independent variable array
        y: Dependent variable array
        
    Returns:
        Integral value
    """
    # Try new name first (NumPy 2.0+), fall back to old name
    try:
        return np.trapezoid(y, x)
    except AttributeError:
        return np.trapz(y, x)


def linear_interp(x: float, x_data: np.ndarray, y_data: np.ndarray) -> float:
    """
    DEPRECATED: Use np.interp() directly.
    
    Linear interpolation.
    
    Args:
        x: Point to interpolate at
        x_data: Known x values
        y_data: Known y values
        
    Returns:
        Interpolated y value
    """
    return np.interp(x, x_data, y_data)


def sign(a: float, b: float) -> float:
    """
    FORTRAN SIGN function: returns abs(a) with sign of b.
    
    Args:
        a: Magnitude value
        b: Sign value
        
    Returns:
        abs(a) * sign(b)
    """
    return np.abs(a) * np.sign(b) if b != 0 else np.abs(a)

