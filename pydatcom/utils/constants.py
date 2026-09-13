"""
Physical and mathematical constants used throughout DATCOM.

Corresponds to FORTRAN COMMON /CONSNT/ block.
Reference: datcom.f line 6

Note: Uses numpy constants directly. Import as:
    from pydatcom.utils.constants import UNUSED
    import numpy as np
    # Use np.pi, np.deg2rad(), np.rad2deg() directly
"""

from typing import Dict, Union

# DEPRECATED: Use np.pi directly
# Kept for backward compatibility with state dict
PI = 3.141592654
DEG = 0.01745329  # Multiply degrees by DEG to obtain radians.
RAD = 57.2957795  # Multiply radians by RAD to obtain degrees.

# DATCOM-specific constants
# Source: BLOCKD, DATA CONST. Keep its precision and sentinel convention
# for translated routines; use np.deg2rad/rad2deg for generic math.
UNUSED = 1.e-30
KAND = '$'  # Portable representation of the source Hollerith 4H$ delimiter.


def get_constants_dict() -> Dict[str, Union[float, str]]:
    """
    Return constants as a dictionary with state manager naming convention.
    
    Returns:
        Dictionary with constants_* prefixed keys
    """
    return {
        'constants_pi': PI,
        'constants_deg': DEG,
        'constants_rad': RAD,
        'constants_unused': UNUSED,
        'constants_kand': KAND,
    }

