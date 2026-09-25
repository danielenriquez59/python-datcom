"""
DATCOM figure and table lookup functions.

Implements empirical correlation figures from DATCOM charts.
These are critical for drag, lift, and moment calculations.

Reference: datcom.f lines 9425 (FIG26), 9467 (FIG53A), 9503 (FIG60B), 9571 (FIG68)
"""

import math
import numpy as np
from typing import Tuple
import logging

from pydatcom.utils.constants import RAD, UNUSED

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
    for mach_bracket in range(8):
        if (mach >= mach_points[mach_bracket]
                and mach < mach_points[mach_bracket + 1]):
            m_idx = mach_bracket
            break
    
    # Log10 of Reynolds number
    x = math.log10(reynolds)
    
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
        -1.32167, -9.04994E-02, -0.489542, -1.80535, -1.22853
    ])
    e_coef = np.array([
        0.493821, -0.735445, -0.317567, 1.02609, 0.471355
    ])
    z_points = np.array([0.0, 1.0, 2.0, 4.0, 10.0])

    # Label 1070 avoids LOG10(0); every other return clips negative R.
    if rv == 0.0:
        return 0.0
    if rv < 0.0:
        raise ValueError("FIG53A requires a nonnegative Reynolds number")
    
    # Find Z range
    z_idx = 4  # Default to last
    for z_bracket in range(4):
        if z >= z_points[z_bracket] and z < z_points[z_bracket + 1]:
            z_idx = z_bracket
            break
    
    # Log10 of RV
    x = math.log10(rv)
    
    # Polynomial evaluation
    def eval_poly(idx):
        return x * (e_coef[idx] + x * (d_coef[idx] + x * (
            c_coef[idx] + x * (b_coef[idx] + x * a_coef[idx]))))
    
    # Interpolate between Z values
    if z_idx < 4 and abs(z - z_points[z_idx]) > 0.001:
        r_z1 = eval_poly(z_idx)
        r_z2 = eval_poly(z_idx + 1)
        
        frac = (z - z_points[z_idx]) / (z_points[z_idx + 1] - z_points[z_idx])
        r = r_z1 + frac * (r_z2 - r_z1)
    else:
        r = eval_poly(z_idx)
    
    return max(float(r), 0.0)


def fig60b(beta: float, btana: float) -> float:
    """FIG60B: fin CNAA for input BETA and BETA*TAN(ALPHA).

    Translates Figure 4.1.3.3-60B, for straight tapered fins with a
    supersonic leading edge and attached shock. BTANA is an INPUT, not
    an output. Returns CN per SIN(ALPHA)**2 as specified by the source.

    First interpolate the source's 9-by-13 table in BETA, then invert
    the resulting BTANA curve. TLINEX modes (1,1,0,0) extrapolate below
    the BETA table and clamp above it. TBFUNX's value is piecewise linear;
    its quadratic derivative is unused here.
    """
    from pydatcom.utils.legacy_tables import tlinex

    cnaa_grid = np.array([0., .2, .4, .6, .8, 1., 1.2, 1.4, 1.6, 1.8, 2., 2.2, 2.4])
    beta_grid = np.array([1.25, 1.5, 1.75, 2., 2.5, 3., 4., 5., 20.])
    # Consecutive groups of 13 in FORTRAN DATA Z60B are columns of
    # Y(NX2,NX1), i.e. a whole CNAA sweep at a fixed BETA.
    table = np.array([
        [0., .047, .089, .127, .160, .188, .211, .231, .247, .258, .266, .269, .270],
        [0., .104, .192, .259, .311, .357, .395, .426, .452, .474, .492, .505, .515],
        [0., .112, .215, .308, .405, .499, .572, .624, .667, .704, .735, .760, .780],
        [0., .132, .270, .400, .521, .637, .732, .820, .878, .925, .962, .986, 1.000],
        [0., .183, .353, .534, .700, .873, 1.040, 1.190, 1.300, 1.388, 1.461, 1.520, 1.563],
        [0., .223, .447, .652, .865, 1.096, 1.325, 1.525, 1.697, 1.842, 1.978, 2.105, 2.222],
        [0., .270, .576, .845, 1.106, 1.447, 1.803, 2.222, 2.511, 2.707, 2.875, 3.014, 3.125],
        [0., .270, .575, .890, 1.170, 1.563, 1.960, 2.473, 2.921, 3.333, 3.744, 4.148, 4.545],
        [0., .306, .629, .925, 1.232, 1.667, 2.262, 3.004, 4.065, 5.814, 9.259, 20.83, 1000.],
    ]).T
    curve = np.zeros(13)
    for cnaa_slot in range(1, 13):
        curve[cnaa_slot] = tlinex(
            beta_grid, cnaa_grid, table, beta, cnaa_grid[cnaa_slot],
            lower1=1, lower2=1, upper1=0, upper2=0)

    query = min(btana, curve[-1])
    # TBFUNX labels 1020-1040 give upper endpoint precedence.
    if query >= curve[-1]:
        return float(cnaa_grid[-1])
    if query <= curve[0]:
        return float(cnaa_grid[0])
    # Value-only translation of TBFUNX labels 1000-1010. Do not use
    # searchsorted: lower-BETA extrapolation can produce a nonmonotonic
    # curve, and the legacy routine scans the entire interior table.
    left = max(bracket_index for bracket_index in range(12)
               if query >= curve[bracket_index])
    left = max(left, 1)
    if query < curve[1]:
        left = 0
    width = curve[left+1] - curve[left]
    if width == 0.:
        width = UNUSED
    return float(cnaa_grid[left] + (cnaa_grid[left+1] - cnaa_grid[left])
                 * (query - curve[left]) / width)


def fig68(mach: float, delta: float) -> Tuple[float, int]:
    """
    Compute shock wave angle from Figure 4.4.1.1-68.
    
    Translate the FIG68 cubic solution for a weak oblique shock (gamma=1.4).
    
    Reference: FORTRAN FIG68 subroutine, datcom.f line 9571
    
    Args:
        mach: Mach number
        delta: Flow deflection angle (degrees)
        
    Returns:
        Tuple of (theta_shock_angle_deg, error_code)
        error_code: 0=attached shock, 1=negative deflection (Mach angle),
            2=detached shock (maximum attached *deflection* angle),
            3=subsonic Mach number (zero angle).
    """
    # DR is initialized to this value in the legacy /CONSNT/ block.
    dr = RAD
    if mach < 1.0:
        return 0.0, 3
    tmin = dr * math.asin(1.0 / mach)
    if mach == 1.0 and delta > 0.0:
        return 0.0, 2
    if delta < 0.0:
        return float(tmin), 1
    if delta == 0.0:
        return float(tmin), 0

    xm2 = mach * mach
    sin2_tmax = (3.0*xm2 - 5.0 + math.sqrt(9.0*xm2*xm2 + 12.0*xm2 + 60.0)) / (7.0*xm2)
    tmax = math.asin(math.sqrt(sin2_tmax))
    dmax = dr * math.atan(1.0 / (math.tan(tmax) * (1.2*xm2 / (xm2*sin2_tmax - 1.0) - 1.0)))
    if delta > dmax:
        # Legacy label 1030 returns wedge angle DMAX, not shock angle TMAX.
        return float(dmax), 2
    if delta == dmax:
        return float(tmax * dr), 0

    # X = sin(shock angle)**2; depress the cubic with X = Y - P/3.
    s2d = math.sin(delta / dr)**2
    p = -(xm2 + 2.0)/xm2 - 1.4*s2d
    q = (2.0*xm2 + 1.0)/(xm2*xm2) + (1.44 + 0.4/xm2)*s2d
    r = -(1.0 - s2d)/(xm2*xm2)
    a = q - p*p/3.0
    b = (2.0*p*p*p - 9.0*p*q + 27.0*r)/27.0
    cosp = -b / (2.0*math.sqrt(-a*a*a/27.0))
    # Clip roundoff at a repeated root (zero deflection or detachment).
    phi3 = math.acos(np.clip(cosp, -1.0, 1.0))/3.0
    s2t = 2.0*math.sqrt(-a/3.0)*math.cos(phi3 + 240.0/dr) - p/3.0
    return float(math.asin(math.sqrt(np.clip(s2t, 0.0, 1.0))) * dr), 0


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
        elif table_name == 'FIG60B' and len(args) == 2:
            return fig60b(args[0], args[1])
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


def angdet(mach: float) -> float:
    """Translate ANGDET: the wedge turn angle at shock detachment.

    Solves NACA TR 1135 equation 168 for the shock angle at maximum
    deflection, then substitutes it into equation 138 to get that
    deflection.  Complements :func:`fig68`, which returns the same limit
    from its own tabulated cubic.

    Args:
        mach: Free-stream Mach number, above one.

    Returns:
        The detachment turn angle in radians.

    Raises:
        ValueError: For a subsonic or sonic Mach number.
    """
    import numpy as _np
    gamma = 1.4
    if mach <= 1.0:
        raise ValueError("ANGDET requires supersonic Mach")

    # Equation 168: sin^2 of the shock angle at maximum deflection.
    sin2 = (((1.0 + gamma) * mach**2 - 4.0 +
             _np.sqrt((gamma + 1.0) * ((gamma + 1.0) * mach**4 +
                                       8.0 * (gamma - 1.0) * mach**2 + 16.0)))
            / (4.0 * gamma * mach**2))
    shock = _np.arcsin(_np.sqrt(sin2))

    # Equation 138: the deflection that shock angle corresponds to.
    denominator = mach**2 * _np.sin(shock)**2 - 1.0
    if denominator == 0.0:
        raise ValueError("ANGDET equation 138 divides by zero at this Mach")
    cot_delta = _np.tan(shock) * (((gamma + 1.0) * mach**2) /
                                  (2.0 * denominator) - 1.0)
    return float(_np.arctan(1.0 / cot_delta))
