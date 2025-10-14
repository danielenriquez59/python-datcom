"""
Physical and mathematical constants used throughout DATCOM.

Corresponds to FORTRAN COMMON /CONSNT/ block.
Reference: datcom.f line 6

Note: Uses numpy constants directly. Import as:
    from pydatcom.utils.constants import UNUSED
    import numpy as np
    # Use np.pi, np.deg2rad(), np.rad2deg() directly
"""

import numpy as np
from typing import Dict

# DEPRECATED: Use np.pi directly
# Kept for backward compatibility with state dict
PI = np.pi
DEG = np.rad2deg(1.0)  # Radians to degrees conversion (use np.rad2deg() instead)
RAD = np.deg2rad(1.0)  # Degrees to radians conversion (use np.deg2rad() instead)

# DATCOM-specific constants
UNUSED = -999.0   # Sentinel value for unused parameters
KAND = 0          # Additional constant from COMMON block


def get_constants_dict() -> Dict[str, float]:
    """
    Return constants as a dictionary with state manager naming convention.
    
    Returns:
        Dictionary with constants_* prefixed keys
    """
    return {
        'constants_pi': np.pi,
        'constants_deg': np.rad2deg(1.0),
        'constants_rad': np.deg2rad(1.0),
        'constants_unused': UNUSED,
        'constants_kand': KAND,
    }

