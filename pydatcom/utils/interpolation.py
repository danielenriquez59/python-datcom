"""
Interpolation functions for PyDATCOM.

Holds ASMINT and a plain bilinear lookup.  The legacy table routines
(GLOOK, TLIN1X, TLINEX, INTER3, INTEP3) live in ``legacy_tables``,
``legacy_numeric`` and ``packed_tables``.

Reference: datcom.f line 516 (ASMINT)
"""

import numpy as np
import logging

logger = logging.getLogger(__name__)


def asmint(x_data: np.ndarray, y_data: np.ndarray, x_vals: np.ndarray) -> np.ndarray:
    """
    Asymmetric interpolation with continuous derivatives.
    
    Uses cubic polynomials between points with parabolic end sections.
    Method produces smooth curves with continuous first derivatives and
    allows extrapolation beyond data range.
    
    Reference: FORTRAN ASMINT subroutine, datcom.f line 516
    
    Algorithm:
    - End intervals: 2nd order polynomial (parabola)
    - Interior intervals: 3rd order polynomial (cubic)
    - Interior slopes: tangent of the average of the two secant angles
    - Continuous derivatives everywhere
    
    Args:
        x_data: Input X data points (must be monotonic)
        y_data: Input Y data points
        x_vals: X values where interpolation is needed
        
    Returns:
        Interpolated Y values at x_vals
    """
    x_data = np.asarray(x_data, dtype=float)
    y_data = np.asarray(y_data, dtype=float)
    x_vals = np.asarray(x_vals, dtype=float)
    if x_data.ndim != 1 or y_data.shape != x_data.shape:
        raise ValueError("ASMINT requires matching one-dimensional X and Y arrays")
    if len(x_data) < 2 or not np.all(np.isfinite(x_data)) or not np.all(np.diff(x_data) > 0):
        raise ValueError("ASMINT requires at least two finite, strictly increasing X values")
    npt = len(x_data)
    if npt < 3:
        # Fall back to linear interpolation for < 3 points
        return np.interp(x_vals, x_data, y_data)
    
    y_vals = np.zeros_like(x_vals, dtype=float)
    
    for query_index, xval in enumerate(x_vals):
        # FORTRAN branches directly to 1080 at a data point. Do not
        # overwrite a knot value by falling through to another interval.
        knot = np.searchsorted(x_data, xval)
        if knot < npt and xval == x_data[knot]:
            y_vals[query_index] = y_data[knot]
            continue
        # Determine location
        if xval < x_data[1]:
            # Left end parabola (extrapolation or first interval)
            # Labels 1020/1040 use the slope at X(2), not X(1).
            left_knot = 1
            right_knot = 2
            locate = 1
        elif xval > x_data[-2]:
            # Right end parabola (extrapolation or last interval)
            left_knot = npt - 2
            right_knot = npt - 1
            locate = 3
        else:
            # Interior cubic
            # Find interval
            for idx in range(1, npt - 2):
                if xval < x_data[idx + 1]:
                    left_knot = idx
                    right_knot = idx + 1
                    locate = 2
                    break
            else:
                # If we didn't break, use last interval
                left_knot = npt - 2
                right_knot = npt - 1
                locate = 2
        
        # Calculate slopes at end points
        yp = np.zeros(2)
        
        for slope_end in range(2):
            lp = (left_knot + slope_end - 1 if left_knot + slope_end - 1 >= 0
                  else 0)
            mp = left_knot + slope_end
            rp = (right_knot + slope_end if right_knot + slope_end < npt
                  else npt - 1)
            
            # Left and right slopes
            if x_data[mp] != x_data[lp]:
                sl = (y_data[lp] - y_data[mp]) / (x_data[lp] - x_data[mp])
            else:
                sl = 0.0
            
            if x_data[rp] != x_data[mp]:
                sr = (y_data[mp] - y_data[rp]) / (x_data[mp] - x_data[rp])
            else:
                sr = 0.0
            
            # Average angle method
            angl = np.arctan(sl)
            angr = np.arctan(sr)
            angav = (angl + angr) / 2.0
            yp[slope_end] = np.tan(angav)
            
            # For end points, only calculate one slope
            if locate != 2:
                break
        
        # Calculate polynomial coefficients
        if locate == 2:
            # Interior: cubic polynomial
            # Hermite form is algebraically identical to labels 1050-1060,
            # but avoids cancellation in powers of absolute X coordinates.
            # In particular, do not zero the cubic for small airfoil intervals.
            width = x_data[right_knot] - x_data[left_knot]
            t = (xval - x_data[left_knot]) / width
            y_vals[query_index] = (
                (2*t**3 - 3*t**2 + 1) * y_data[left_knot]
                + (t**3 - 2*t**2 + t) * width * yp[0]
                + (-2*t**3 + 3*t**2) * y_data[right_knot]
                + (t**3 - t**2) * width * yp[1]
            )
            continue
        else:
            # End parabola
            j_idx = 0 if locate == 1 else npt - 2
            k_idx = j_idx + 1
            width = x_data[k_idx] - x_data[j_idx]
            secant = (y_data[k_idx] - y_data[j_idx]) / width
            # Labels 1070 constrain the derivative at the inner endpoint.
            curvature = (yp[0] - secant) / width
            if locate == 3:
                curvature = -curvature
            offset = xval - x_data[j_idx]
            y_vals[query_index] = (y_data[j_idx] + secant * offset
                                   + curvature * offset * (offset - width))

    return y_vals


def bilinear_interp(x: float, y: float,
                    x_data: np.ndarray, y_data: np.ndarray,
                    z_table: np.ndarray) -> float:
    """
    Bilinear interpolation on 2D grid.
    
    Args:
        x: X value to interpolate at
        y: Y value to interpolate at
        x_data: X grid points
        y_data: Y grid points
        z_table: 2D table of Z values [len(x_data), len(y_data)]
        
    Returns:
        Interpolated Z value
    """
    # Find bracketing indices
    if x <= x_data[0]:
        i1, i2 = 0, 0
        fx = 0.0
    elif x >= x_data[-1]:
        i1, i2 = len(x_data) - 1, len(x_data) - 1
        fx = 0.0
    else:
        i2 = np.searchsorted(x_data, x)
        i1 = i2 - 1
        fx = (x - x_data[i1]) / (x_data[i2] - x_data[i1])
    
    if y <= y_data[0]:
        j1, j2 = 0, 0
        fy = 0.0
    elif y >= y_data[-1]:
        j1, j2 = len(y_data) - 1, len(y_data) - 1
        fy = 0.0
    else:
        j2 = np.searchsorted(y_data, y)
        j1 = j2 - 1
        fy = (y - y_data[j1]) / (y_data[j2] - y_data[j1])
    
    # Bilinear interpolation
    if i1 == i2 and j1 == j2:
        return z_table[i1, j1]
    elif i1 == i2:
        return z_table[i1, j1] * (1 - fy) + z_table[i1, j2] * fy
    elif j1 == j2:
        return z_table[i1, j1] * (1 - fx) + z_table[i2, j1] * fx
    else:
        # Full bilinear
        z11 = z_table[i1, j1]
        z12 = z_table[i1, j2]
        z21 = z_table[i2, j1]
        z22 = z_table[i2, j2]
        
        return (z11 * (1 - fx) * (1 - fy) +
                z21 * fx * (1 - fy) +
                z12 * (1 - fx) * fy +
                z22 * fx * fy)
