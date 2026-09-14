"""
Supersonic downwash, Figure 4.4.1-76.

The figure is split across five source routines by wing taper ratio, each a
chain of table lookups in which every stage feeds the next:

- ``SDWA``  Figure 4.4.1-76A, taper 0.00
- ``SDWB``  Figure 4.4.1-76B, taper 1.00
- ``SDWC``  Figure 4.4.1-76C, taper 0.25 and 0.50
- ``SDWD``  Figure 4.4.1-76D, taper 0.00, 0.25 and 0.50
- ``SDWE``  Figure 4.4.1-76E, general taper

Each stage is an ``INTERX`` two-variable lookup whose dependent tables use
the flat ``TABLE(LIND,4)`` packing, so the grid for variable ``k`` starts at
offset ``k*LIND`` regardless of how many entries the previous grid used.
The source end modes are carried through per call.

Reference: datcom-legacy/datcom_2000/sdwa.f, sdwb.f, sdwc.f, sdwd.f, sdwe.f
"""

import numpy as np
from typing import Dict
import logging

from pydatcom.utils.legacy_interp import interx

logger = logging.getLogger(__name__)

# --------------------------------------------------------------------------
# Figure 4.4.1-76A, taper 0.00.  sdwa.f
# --------------------------------------------------------------------------
_XA1 = [.90, 1.0, 1.2, 1.4, 1.6, 2.4, 0.0, 0.1, 0.3, 0.5]
_YA1 = [.90, 1.18, 1.65, 1.95, 2.10, 2.34,
        .55, 0.84, 1.30, 1.60, 1.74, 2.00,
        0.0, 0.00, 0.50, 0.87, 1.10, 1.50,
        0.0, 0.00, 0.00, 0.30, 0.61, 1.00]
_XA2 = [0.0, 3.0, 0.0, 4.0, 8.0, 12.]
_YA2 = [0.0, 4.45, 0.0, 3.00, 0.0, 2.50]
_XA3 = [0.0, 3.0, 0.0, 0.0, 0.0, 0.0, 0.3, 0.4, 0.6, 0.8]
_YA3 = [0.00, 0.60, -.04, 0.54, -.09, 0.44, -.17, 0.35, -.28, 0.26]

# --------------------------------------------------------------------------
# Figure 4.4.1-76B, taper 1.00.  sdwb.f
# --------------------------------------------------------------------------
_XB1 = [1.0, 1.2, 1.4, 1.5, 1.6, 2.0, 2.2, 2.4, 0.0, 0.1, 0.3, 0.5]
_YB1 = [0.32, 0.74, 1.21, 1.39, 1.44, 1.56, 1.62, 1.77,
        0.11, 0.60, 1.10, 1.25, 1.30, 1.42, 1.50, 1.62,
        0.00, 0.20, 0.76, 1.00, 1.10, 1.34, 1.42, 1.47,
        0.00, 0.00, 0.00, 0.50, 0.75, 1.06, 1.12, 1.20]
_XB2 = [0.00, .8, 1.70, 3.0, 2.00, 4., 8.00, 12.]
_YB2 = [0.00, .60, 4.43, 12.7,
        0.00, .80, 1.70, 3.00,
        0.00, .40, 0.80, 1.40,
        0.00, .27, 0.50, 0.94]
_XB3 = [0.0, 5.0, 0.0, 0.0, .15, .30]
_YB3 = [0.0, 1.00, 0.0, 1.10, 0.0, 1.20]


def calculate_sdwa(x: float, y: float, z: float,
                   beta_aspect: float) -> Dict[str, float]:
    """Translate SDWA: Figure 4.4.1-76A downwash, taper 0.00.

    Three chained lookups: 76A1 on ``(x, z)``, 76A2 on that result and the
    compressible aspect-ratio parameter, then 76A3 on that result and ``y``.

    Args:
        x: Streamwise station.
        y: Spanwise station.
        z: Vertical station.
        beta_aspect: ``AB``, the compressible aspect-ratio parameter.

    Returns:
        Dictionary with ``sdw`` and the two intermediate stage values.
    """
    stage1 = interx(2, _XA1, [x, z], [6, 4], _YA1, lind=6,
                    lx1l=2, lx2l=0, lx1u=1, lx2u=2)
    stage2 = interx(2, _XA2, [stage1, beta_aspect], [2, 3], _YA2, lind=3,
                    lx1l=0, lx2l=1, lx1u=1, lx2u=1)
    sdw = interx(2, _XA3, [stage2, y], [2, 5], _YA3, lind=5,
                 lx1l=0, lx2l=0, lx1u=1, lx2u=2)
    return {
        'sdw': float(sdw),
        'stage1': float(stage1),
        'stage2': float(stage2),
        'taper': 0.00,
        'method': 'legacy_sdwa',
    }


def calculate_sdwb(x: float, y: float, z: float,
                   beta_aspect: float) -> Dict[str, float]:
    """Translate SDWB: Figure 4.4.1-76B downwash, taper 1.00.

    Args:
        x: Streamwise station.
        y: Spanwise station.
        z: Vertical station.
        beta_aspect: ``AB``, the compressible aspect-ratio parameter.

    Returns:
        Dictionary with ``sdw`` and the two intermediate stage values.
    """
    stage1 = interx(2, _XB1, [x, z], [8, 4], _YB1, lind=8,
                    lx1l=2, lx2l=0, lx1u=1, lx2u=2)
    stage2 = interx(2, _XB2, [stage1, beta_aspect], [4, 4], _YB2, lind=4,
                    lx1l=0, lx2l=2, lx1u=1, lx2u=2)
    sdw = interx(2, _XB3, [stage2, y], [2, 3], _YB3, lind=3,
                 lx1l=0, lx2l=0, lx1u=1, lx2u=2)
    return {
        'sdw': float(sdw),
        'stage1': float(stage1),
        'stage2': float(stage2),
        'taper': 1.00,
        'method': 'legacy_sdwb',
    }
