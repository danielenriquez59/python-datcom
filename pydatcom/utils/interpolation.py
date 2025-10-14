"""
Interpolation functions for PyDATCOM.

Provides various interpolation methods used throughout DATCOM:
- Linear interpolation
- Cubic spline interpolation
- Asymmetric interpolation (ASMINT)
- Multi-dimensional table lookup (INTER3, INTEP3)
- General lookup (GLOOK)

Reference: datcom.f lines 516 (ASMINT), 15275 (INTEP3), 15344 (INTER3), 10761 (GLOOK)
"""

import numpy as np
from scipy import interpolate
from typing import Tuple, Optional, Callable
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
    - Slopes at each point: average of left and right linear slopes
    - Continuous derivatives everywhere
    
    Args:
        x_data: Input X data points (must be monotonic)
        y_data: Input Y data points
        x_vals: X values where interpolation is needed
        
    Returns:
        Interpolated Y values at x_vals
    """
    npt = len(x_data)
    if npt < 3:
        # Fall back to linear interpolation for < 3 points
        return np.interp(x_vals, x_data, y_data)
    
    y_vals = np.zeros_like(x_vals, dtype=float)
    
    for i, xval in enumerate(x_vals):
        # Determine location
        if xval < x_data[1]:
            # Left end parabola (extrapolation or first interval)
            j = 0
            k = 1
            locate = 1
        elif xval > x_data[-2]:
            # Right end parabola (extrapolation or last interval)
            j = npt - 2
            k = npt - 1
            locate = 3
        else:
            # Interior cubic
            # Find interval
            for idx in range(1, npt - 2):
                if xval == x_data[idx]:
                    y_vals[i] = y_data[idx]
                    continue
                if xval < x_data[idx + 1]:
                    j = idx
                    k = idx + 1
                    locate = 2
                    break
            else:
                # If we didn't break, use last interval
                j = npt - 2
                k = npt - 1
                locate = 2
        
        # Calculate slopes at end points
        yp = np.zeros(2)
        
        for n in range(2):
            lp = j + n - 1 if j + n - 1 >= 0 else 0
            mp = j + n
            rp = k + n if k + n < npt else npt - 1
            
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
            yp[n] = np.tan(angav)
            
            # For end points, only calculate one slope
            if locate != 2:
                break
        
        # Calculate polynomial coefficients
        if locate == 2:
            # Interior: cubic polynomial
            x1s = x_data[j]**2
            x2s = x_data[k]**2
            x12f = x_data[j] - x_data[k]
            y12f = y_data[j] - y_data[k]
            x12s = x1s - x2s
            x12c = x1s * x_data[j] - x2s * x_data[k]
            y12p = yp[0] - yp[1]
            
            tw = 2.0
            th = 3.0
            red = tw * x_data[k] * x12f - x12s
            grn = th * x2s * x12f - x12c
            yel = yp[1] * x12f - y12f
            e = th * x12s * red - tw * x12f * grn
            
            if abs(e) > 1e-10:
                a = (y12p * red - tw * x12f * yel) / e
                b = (th * x12s * yel - y12p * grn) / e
            else:
                a = 0.0
                b = 0.0
            
            c = (y12f - a * x12c - b * x12s) / x12f if abs(x12f) > 1e-10 else 0.0
            d = y_data[k] - a * x2s * x_data[k] - b * x2s - c * x_data[k]
        else:
            # End parabola
            j_idx = 0 if locate == 1 else npt - 2
            k_idx = 1 if locate == 1 else npt - 1
            l_idx = 1 if locate == 1 else npt - 2
            
            if x_data[k_idx] != x_data[j_idx]:
                z = (y_data[j_idx] - y_data[k_idx]) / (x_data[j_idx] - x_data[k_idx])
            else:
                z = 0.0
            
            a = 0.0
            denom = 2.0 * x_data[l_idx] - x_data[j_idx] - x_data[k_idx]
            if abs(denom) > 1e-10:
                b = (yp[0] - z) / denom
            else:
                b = 0.0
            c = yp[0] - 2.0 * b * x_data[l_idx]
            d = y_data[j_idx] - ((b * x_data[j_idx] + c) * x_data[j_idx])
        
        # Evaluate polynomial
        y_vals[i] = (((a * xval + b) * xval) + c) * xval + d
    
    return y_vals


def linear_interpolation_2d(arg1: float, arg2: float,
                            x1_data: np.ndarray, x2_data: np.ndarray,
                            y_table: np.ndarray) -> float:
    """
    2D linear interpolation (bilinear).
    
    Args:
        arg1: First argument value
        arg2: Second argument value
        x1_data: First dimension data points
        x2_data: Second dimension data points
        y_table: 2D table of values [len(x1_data), len(x2_data)]
        
    Returns:
        Interpolated value
    """
    from scipy.interpolate import RectBivariateSpline
    
    # Use scipy for robust 2D interpolation
    interp_func = RectBivariateSpline(x1_data, x2_data, y_table, kx=1, ky=1)
    return float(interp_func(arg1, arg2))


def glook(x_grid: np.ndarray, x_alpha: np.ndarray, 
          nas_grid: np.ndarray, table_grid: np.ndarray) -> np.ndarray:
    """
    General table lookup with interpolation.
    
    Reference: FORTRAN GLOOK subroutine, datcom.f line 10761
    
    Args:
        x_grid: Grid X values
        x_alpha: Alpha schedule values
        nas_grid: Number of alpha points per grid station
        table_grid: Table of values to interpolate
        
    Returns:
        Interpolated values at x_alpha points
    """
    # Simplified implementation using linear interpolation
    # Full FORTRAN version has complex multi-dimensional logic
    
    nxg = len(x_grid)
    nval = len(x_alpha)
    result = np.zeros(nval)
    
    for i in range(nval):
        # Simple 1D interpolation
        # Full implementation would handle 2D grid
        result[i] = np.interp(x_alpha[i], x_grid, table_grid[:nxg])
    
    return result


class TableInterpolator:
    """
    Multi-dimensional table interpolation for DATCOM figures.
    
    Handles complex table lookups with multiple dimensions and
    taper ratio (lambda) variations.
    
    Reference: INTER3, INTEP3 subroutines
    """
    
    def __init__(self, table_data: dict):
        """
        Initialize with table data.
        
        Args:
            table_data: Dictionary containing table arrays and metadata
        """
        self.table_data = table_data
        self._cache = {}
    
    def lookup_2d(self, arg1: float, arg2: float, table_name: str) -> float:
        """
        Perform 2D table lookup.
        
        Args:
            arg1: First argument (e.g., Mach number)
            arg2: Second argument (e.g., angle of attack)
            table_name: Name of table to lookup
            
        Returns:
            Interpolated value
        """
        if table_name not in self.table_data:
            logger.error(f"Table {table_name} not found")
            return 0.0
        
        table = self.table_data[table_name]
        
        # Extract table dimensions
        x1_data = np.array(table.get('x1', []))
        x2_data = np.array(table.get('x2', []))
        y_data = np.array(table.get('y', []))
        
        if len(x1_data) == 0 or len(x2_data) == 0:
            return 0.0
        
        # Perform 2D interpolation
        return linear_interpolation_2d(arg1, arg2, x1_data, x2_data, y_data)
    
    def lookup_3d(self, arg1: float, arg2: float, lambda_val: float,
                  table_set: str) -> float:
        """
        Perform 3D table lookup with taper ratio (lambda).
        
        This corresponds to INTEP3/INTER3 subroutines which handle
        tables at different taper ratios.
        
        Args:
            arg1: First argument
            arg2: Second argument
            lambda_val: Taper ratio (0.0 to 1.0)
            table_set: Set of tables at different lambda values
            
        Returns:
            Interpolated value
        """
        # Determine which lambda tables to use
        if lambda_val <= 0.0:
            it = 1  # Use lambda=0 table only
        elif lambda_val <= 0.25:
            it = 2  # Interpolate between 0 and 0.25
        elif lambda_val <= 0.50:
            it = 3  # Interpolate between 0.25 and 0.50
        elif lambda_val <= 0.75:
            it = 4  # Interpolate between 0.50 and 0.75
        else:
            it = 5  # Interpolate between 0.75 and 1.0
        
        # Get table data for appropriate lambda values
        # Simplified: use linear interpolation on lambda
        lambda_tables = {
            0.0: f"{table_set}_lam00",
            0.25: f"{table_set}_lam25",
            0.50: f"{table_set}_lam50",
            0.75: f"{table_set}_lam75",
            1.0: f"{table_set}_lam100",
        }
        
        # Find bracketing lambda values
        lambda_points = [0.0, 0.25, 0.50, 0.75, 1.0]
        idx = np.searchsorted(lambda_points, lambda_val)
        
        if idx == 0:
            # Below first point, use first table
            return self.lookup_2d(arg1, arg2, lambda_tables[0.0])
        elif idx >= len(lambda_points):
            # Above last point, use last table
            return self.lookup_2d(arg1, arg2, lambda_tables[1.0])
        else:
            # Interpolate between two tables
            lam1 = lambda_points[idx - 1]
            lam2 = lambda_points[idx]
            
            val1 = self.lookup_2d(arg1, arg2, lambda_tables.get(lam1, table_set))
            val2 = self.lookup_2d(arg1, arg2, lambda_tables.get(lam2, table_set))
            
            # Linear interpolation on lambda
            frac = (lambda_val - lam1) / (lam2 - lam1) if lam2 != lam1 else 0.0
            return val1 + (val2 - val1) * frac


def cubic_spline_interp(x_data: np.ndarray, y_data: np.ndarray, 
                        x_vals: np.ndarray, extrapolate: bool = True) -> np.ndarray:
    """
    Cubic spline interpolation with optional extrapolation.
    
    Args:
        x_data: Known X values
        y_data: Known Y values
        x_vals: X values to interpolate at
        extrapolate: Allow extrapolation beyond data range
        
    Returns:
        Interpolated Y values
    """
    if extrapolate:
        spline = interpolate.interp1d(x_data, y_data, kind='cubic', 
                                     fill_value='extrapolate', bounds_error=False)
    else:
        spline = interpolate.interp1d(x_data, y_data, kind='cubic', 
                                     fill_value=np.nan, bounds_error=False)
    
    return spline(x_vals)


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


def find_nearest_index(value: float, array: np.ndarray) -> int:
    """
    Find index of nearest value in array.
    
    Args:
        value: Value to find
        array: Array to search
        
    Returns:
        Index of nearest element
    """
    idx = np.abs(array - value).argmin()
    return int(idx)


def interpolate_along_curve(s: np.ndarray, x: np.ndarray, y: np.ndarray,
                            s_vals: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """
    Interpolate X and Y coordinates along a curve parameterized by arc length s.
    
    Args:
        s: Arc length parameter at known points
        x: X coordinates at known points
        y: Y coordinates at known points
        s_vals: Arc length values where interpolation is needed
        
    Returns:
        Tuple of (x_vals, y_vals) interpolated coordinates
    """
    x_interp = np.interp(s_vals, s, x)
    y_interp = np.interp(s_vals, s, y)
    
    return x_interp, y_interp

