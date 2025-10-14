"""
DATCOM figure and table lookup functions.

Implements empirical correlation figures from DATCOM charts.
These are critical for drag, lift, and moment calculations.

Reference: datcom.f lines 9425 (FIG26), 9467 (FIG53A), 9503 (FIG60B), 9571 (FIG68)
"""

import numpy as np
from typing import Tuple
import logging

logger = logging.getLogger(__name__)


def fig26(reynolds: float, mach: float) -> float:
    """
    Compute skin friction coefficient from Figure 4.1.5.1-26.
    
    Turbulent flat plate skin friction coefficient as function of
    Reynolds number and Mach number.
    
    Reference: FORTRAN FIG26 subroutine, datcom.f line 9425
    
    Args:
        reynolds: Reynolds number
        mach: Mach number
        
    Returns:
        Skin friction coefficient (CF)
    """
    # Polynomial coefficients for different Mach numbers
    # These come from curve fits to DATCOM Figure 4.1.5.1-26
    a_coef = np.array([
        4.12963E-6, 3.92725E-6, 4.55853E-6, 4.49735E-6, 4.11442E-6,
        4.51587E-6, 4.61971E-6, 4.53836E-6, 3.86772E-6
    ])
    b_coef = np.array([
        -1.36204E-4, -1.3037E-4, -1.48715E-4, -1.47407E-4, -1.36505E-4,
        -1.47222E-4, -1.48773E-4, -1.44676E-4, -1.23287E-4
    ])
    c_coef = np.array([
        1.7162E-3, 1.65388E-3, 1.85005E-3, 1.83955E-3, 1.72383E-3,
        1.8252E-3, 1.81973E-3, 1.74944E-3, 1.49335E-3
    ])
    d_coef = np.array([
        -9.88935E-3, -9.59519E-3, -1.0503E-2, -1.04587E-2, -9.91294E-3,
        -1.02911E-2, -1.01084E-2, -9.59421E-3, -8.22227E-3
    ])
    e_coef = np.array([
        2.23641E-2, 2.18366E-2, 2.33437E-2, 2.32398E-2, 2.22626E-2,
        2.2622E-2, 2.18584E-2, 2.04571E-2, 1.76472E-2
    ])
    mach_points = np.array([0.0, 0.3, 0.7, 0.9, 1.0, 1.5, 2.0, 2.5, 3.0])
    
    # Find Mach number range
    m_idx = 8  # Default to last range
    for m in range(8):
        if mach >= mach_points[m] and mach < mach_points[m + 1]:
            m_idx = m
            break
    
    # Log10 of Reynolds number
    x = np.log10(reynolds)
    
    # Polynomial evaluation
    def eval_poly(idx):
        return x * (e_coef[idx] + x * (d_coef[idx] + x * (
            c_coef[idx] + x * (b_coef[idx] + x * a_coef[idx]))))
    
    # Interpolate between Mach numbers if needed
    if m_idx < 8 and abs(mach - mach_points[m_idx]) > 0.02:
        cf_m = eval_poly(m_idx)
        cf_n = eval_poly(m_idx + 1)
        
        # Linear interpolation on Mach
        frac = (mach - mach_points[m_idx]) / (mach_points[m_idx + 1] - mach_points[m_idx])
        cf = cf_m + frac * (cf_n - cf_m)
    else:
        cf = eval_poly(m_idx)
    
    return cf


def fig53a(rv: float, z: float) -> float:
    """
    Compute from Figure 4.1.5.2-53A.
    
    Reference: FORTRAN FIG53A subroutine, datcom.f line 9467
    
    Args:
        rv: Reynolds number based on volume
        z: Z parameter
        
    Returns:
        R value from figure
    """
    # Polynomial coefficients for different Z values
    a_coef = np.array([
        8.98425E-03, 4.50351E-03, 6.0128E-03, 1.07637E-02, 8.48758E-03
    ])
    b_coef = np.array([
        -0.138262, -0.064239, -0.08858, -0.167056, -0.130342
    ])
    c_coef = np.array([
        0.718213, 0.263365, 0.410583, 0.893775, 0.67405
    ])
    d_coef = np.array([
        -1.49813, -0.478074, -0.771195, -1.79298, -1.35093
    ])
    e_coef = np.array([
        1.00929, 0.162566, 0.423931, 1.27184, 0.96068
    ])
    z_points = np.array([0.0, 0.25, 0.50, 0.75, 1.0])
    
    # Find Z range
    z_idx = 4  # Default to last
    for i in range(4):
        if z >= z_points[i] and z < z_points[i + 1]:
            z_idx = i
            break
    
    # Log10 of RV
    x = np.log10(rv)
    
    # Polynomial evaluation
    def eval_poly(idx):
        return x * (e_coef[idx] + x * (d_coef[idx] + x * (
            c_coef[idx] + x * (b_coef[idx] + x * a_coef[idx]))))
    
    # Interpolate between Z values
    if z_idx < 4:
        r_z1 = eval_poly(z_idx)
        r_z2 = eval_poly(z_idx + 1)
        
        frac = (z - z_points[z_idx]) / (z_points[z_idx + 1] - z_points[z_idx])
        r = r_z1 + frac * (r_z2 - r_z1)
    else:
        r = eval_poly(z_idx)
    
    return r


def fig60b(beta: float) -> Tuple[float, float]:
    """
    Compute BTANA and CNAA from Figure 4.3.1.1-60B.
    
    Reference: FORTRAN FIG60B subroutine, datcom.f line 9503
    
    Args:
        beta: Prandtl-Glauert parameter sqrt(|M²-1|)
        
    Returns:
        Tuple of (btana, cnaa)
    """
    # Data from DATCOM Figure 60B
    beta_data = np.array([
        0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8,
        2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4, 3.6, 3.8, 4.0
    ])
    
    btana_data = np.array([
        0.0, 0.0, 0.0, 0.011, 0.045, 0.095, 0.16, 0.229, 0.30, 0.369,
        0.435, 0.494, 0.548, 0.596, 0.640, 0.679, 0.715, 0.748, 0.778,
        0.806, 0.832
    ])
    
    cnaa_data = np.array([
        6.28, 6.28, 6.29, 6.29, 6.30, 6.31, 6.32, 6.33, 6.34, 6.35,
        6.36, 6.37, 6.38, 6.39, 6.40, 6.40, 6.41, 6.42, 6.42, 6.43, 6.43
    ])
    
    # Interpolate
    btana = np.interp(beta, beta_data, btana_data)
    cnaa = np.interp(beta, beta_data, cnaa_data)
    
    return btana, cnaa


def fig68(mach: float, delta: float) -> Tuple[float, int]:
    """
    Compute shock wave angle from Figure 4.4.1.1-68.
    
    Oblique shock relations for given Mach number and deflection angle.
    
    Reference: FORTRAN FIG68 subroutine, datcom.f line 9571
    
    Args:
        mach: Mach number
        delta: Flow deflection angle (degrees)
        
    Returns:
        Tuple of (theta_shock_angle_deg, error_code)
        error_code: 0=success, 1=no solution exists
    """
    # Simplified implementation using shock relations
    # Full implementation would use tabulated data
    
    if mach <= 1.0:
        logger.warning(f"FIG68 called with subsonic Mach {mach}")
        return 0.0, 1
    
    # Approximate shock angle using theta-beta-M relation
    # For small deflections: theta ≈ delta * sqrt((M²-1)/M²)
    beta_prandtl = np.sqrt(mach**2 - 1.0)
    delta_rad = np.deg2rad(delta)
    
    # Iterative solution for shock angle (simplified)
    # Full version would use Newton-Raphson on oblique shock equation
    theta_approx = np.rad2deg(np.arcsin(1.0 / mach))  # Mach angle
    theta = theta_approx + delta  # Approximate shock angle
    
    # Validate solution exists
    max_deflection = np.rad2deg(np.arcsin(1.0 / mach)) * 2.0
    if delta > max_deflection:
        logger.warning(f"Deflection {delta}° exceeds max for M={mach}")
        return 0.0, 1
    
    return theta, 0


class DatcomTableManager:
    """
    Manager for DATCOM table data.
    
    Loads and caches DATCOM figure data for efficient lookup.
    Tables are stored in YAML/JSON format in data/tables/ directory.
    """
    
    def __init__(self):
        """Initialize table manager with empty cache."""
        self._cache = {}
        self._tables_loaded = False
    
    def load_table(self, table_name: str) -> dict:
        """
        Load table data from file.
        
        Args:
            table_name: Name of table (e.g., 'FIG26', 'FIG53A')
            
        Returns:
            Dictionary with table data
        """
        if table_name in self._cache:
            return self._cache[table_name]
        
        # For now, return empty dict
        # Full implementation would load from pydatcom/data/tables/
        logger.warning(f"Table {table_name} not yet loaded from file")
        self._cache[table_name] = {}
        return {}
    
    def lookup(self, table_name: str, *args) -> float:
        """
        Perform table lookup with interpolation.
        
        Args:
            table_name: Name of table
            *args: Lookup arguments (depends on table)
            
        Returns:
            Interpolated value
        """
        # Delegate to specific figure functions
        if table_name == 'FIG26' and len(args) == 2:
            return fig26(args[0], args[1])
        elif table_name == 'FIG53A' and len(args) == 2:
            return fig53a(args[0], args[1])
        elif table_name == 'FIG60B' and len(args) == 1:
            btana, cnaa = fig60b(args[0])
            return btana  # Return first value
        elif table_name == 'FIG68' and len(args) == 2:
            theta, ierr = fig68(args[0], args[1])
            return theta
        else:
            logger.error(f"Unknown table {table_name} or wrong arguments")
            return 0.0


# Create global table manager instance
_table_manager = DatcomTableManager()


def get_table_manager() -> DatcomTableManager:
    """Get the global table manager instance."""
    return _table_manager

