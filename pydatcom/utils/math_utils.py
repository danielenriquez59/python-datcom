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
    Calculate area using trapezoidal rule and centroid.
    
    Reference: FORTRAN AREA2, datcom.f line 471
    
    Args:
        x: Array of x-coordinates
        y: Array of y-coordinates  
        inum: Number of points
        
    Returns:
        Tuple of (area, x_centroid, y_centroid)
    """
    area = 0.0
    ax = 0.0
    ay = 0.0
    
    for i in range(inum - 1):
        # Trapezoidal area increment
        da = (y[i] + y[i + 1]) * (x[i + 1] - x[i]) / 2.0
        area += da
        
        # Centroid contributions
        dx = x[i + 1] - x[i]
        dy = y[i + 1] - y[i]
        
        # X centroid contribution
        dax = (x[i] + dx / 3.0) * (y[i] * dx + dy * dx / 2.0)
        ax += dax
        
        # Y centroid contribution  
        day = (y[i] / 2.0 + dy / 3.0) * (y[i] * dx + dy * dx / 2.0)
        ay += day
    
    # Finalize centroids
    if abs(area) > 1e-10:
        ax = ax / area
        ay = ay / area
    
    return area, ax, ay


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

