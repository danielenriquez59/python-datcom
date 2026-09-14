"""
SDDVC: downward displacement of the vortex core, Figure 4.4.1-74.

Feeds the supersonic downwash routines.  The figure comes in four parts,
each covering a different taper ratio, and the result is interpolated across
taper at the end:

- 74A, taper 0.0
- 74B, taper 0.25 and 0.50, with the half chord unswept
- 74C, taper 0.25 and 0.50, with the trailing edge unswept
- 74D, taper 1.0

Which of 74B and 74C applies is selected by ``icase``; case 3 uses both and
blends them on the trailing-edge sweep.

Reference: datcom-legacy/datcom_2000/sddvc.f
"""

import numpy as np
from typing import Dict, Sequence
import logging

from pydatcom.utils.legacy_interp import interx
from pydatcom.utils.legacy_numeric import tbfunx

logger = logging.getLogger(__name__)

# The four taper ratios the figure parts correspond to.
_TAPER_GRID = np.array([0.0, 0.25, 0.50, 1.0])

# Figure 4.4.1-74A, taper 0.0.  sddvc.f DATA XA / YA, LIND=5.
_XA = [.5, 1.5, 3.0, 5.0, 6.0, 8., 12.]
_YA = [.5, .81, 1.50, 2.47, 2.93,
       .3, .61, 1.21, 2.00, 2.38]

# Figure 4.4.1-74B, half chord unswept.  LIND=8.  YB1 is taper 0.25 and
# YB2 taper 0.50.
_XB = [.5, .75, 1.50, 2.00, 2.50, 3.00, 3.50, 6.00, 2., 4., 8., 12.]
_YB1 = [1.25, 1.60, 1.60, 2.07, 2.96, 3.81, 4.67, 9.03,
        .700, .750, 1.12, 1.49, 1.94, 2.51, 2.51, 2.51,
        .410, .460, .700, .900, 1.11, 1.32, 1.52, 2.58,
        .310, .330, .500, .630, .800, .960, 1.11, 1.89]
_YB2 = [1.10, 1.50, 1.50, 1.96, 2.75, 3.56, 4.38, 8.33,
        .600, .670, 1.00, 1.30, 1.63, 2.08, 2.52, 2.52,
        .310, .340, .540, .700, .860, 1.02, 1.18, 1.94,
        .220, .221, .320, .450, .580, .700, .810, 1.40]

# Figure 4.4.1-74C, trailing edge unswept.  LIND=5.
_XC = [.5, 1.0, 2.0, 3.5, 6.0, 3.2, 6.4, 12.8]
_YC1 = [.910, 1.03, 1.46, 2.41, 4.16,
        .500, .600, .940, 1.58, 2.66,
        .250, .320, .620, 1.08, 1.81]
_YC2 = [.910, 1.03, 1.40, 2.31, 3.99,
        .500, .560, .880, 1.46, 2.41,
        .250, .320, .570, .960, 1.59]

# Figure 4.4.1-74D, taper 1.0.  LIND=5.
_XD = [1.0, 1.5, 2.0, 3.75, 6.0, 2., 4., 8., 12.]
_YD = [1.00, 1.41, 1.88, 4.00, 6.77,
       .510, .600, .710, 1.37, 2.20,
       .290, .300, .370, .660, 1.02,
       .130, .120, .200, .400, .690]


def calculate_sddvc(x: Sequence[float], beta_aspect: float,
                    taper_ratio: float, icase: int,
                    sweep_te: float = 0.0,
                    sweep_reference: float = 1.0) -> Dict[str, object]:
    """Translate SDDVC: vortex core downward displacement.

    Args:
        x: Two streamwise stations, the source's ``X(2)``.
        beta_aspect: ``ABETA``, the compressible aspect-ratio parameter.
        taper_ratio: ``TAPR``, the wing taper ratio to interpolate to.
        icase: Selects which mid-taper figure applies.  1 uses Figure
            4.4.1-74C, 2 uses 74B, and 3 uses both and blends them on the
            trailing-edge sweep.
        sweep_te: ``SWEPTE``, used only by case 3.
        sweep_reference: ``SWEPR``, the blend denominator for case 3.

    Returns:
        Dictionary with ``dhb``, the displacement at each station, and
        ``per_taper``, a 2x4 array of the four figure-part values behind
        each station's result.

    Raises:
        ValueError: If ``x`` is not two stations, ``icase`` is not 1, 2 or 3,
            or case 3 is requested with a zero reference sweep.
    """
    x = np.asarray(x, dtype=float)
    if x.shape != (2,):
        raise ValueError("SDDVC expects exactly two streamwise stations")
    if icase not in (1, 2, 3):
        raise ValueError("SDDVC ICASE must be 1, 2 or 3")
    if icase == 3 and sweep_reference == 0.0:
        raise ValueError("SDDVC case 3 divides by the reference sweep")

    displacement = np.empty(2)
    per_taper = np.empty((2, 4))
    for station in range(2):
        query = [float(x[station]), float(beta_aspect)]
        values = np.empty(4)

        # Figure 4.4.1-74A: taper 0.
        values[0] = interx(2, _XA, query, [5, 2], _YA, lind=5,
                           lx1l=1, lx2l=1, lx1u=1, lx2u=1)

        # Figure 4.4.1-74B: applied for cases 2 and 3.
        blend_b = None
        if icase in (2, 3):
            b1 = interx(2, _XB, query, [8, 4], _YB1, lind=8,
                        lx1l=2, lx2l=2, lx1u=2, lx2u=2)
            b2 = interx(2, _XB, query, [8, 4], _YB2, lind=8,
                        lx1l=2, lx2l=2, lx1u=2, lx2u=2)
            values[1], values[2] = b1, b2
            if icase == 3:
                blend_b = (b1, b2)

        # Figure 4.4.1-74C: applied for every case except 2.
        if icase != 2:
            c1 = interx(2, _XC, query, [5, 3], _YC1, lind=5,
                        lx1l=2, lx2l=2, lx1u=1, lx2u=2)
            c2 = interx(2, _XC, query, [5, 3], _YC2, lind=5,
                        lx1l=2, lx2l=2, lx1u=1, lx2u=2)
            if blend_b is not None:
                # Case 3 interpolates from the 74B value towards the 74C
                # value in proportion to the trailing-edge sweep.
                values[1] = blend_b[0] + sweep_te * (c1 - blend_b[0]) / sweep_reference
                values[2] = blend_b[1] + sweep_te * (c2 - blend_b[1]) / sweep_reference
            else:
                values[1], values[2] = c1, c2

        # Figure 4.4.1-74D: taper 1.
        values[3] = interx(2, _XD, query, [5, 4], _YD, lind=5,
                           lx1l=2, lx2l=2, lx1u=1, lx2u=2)

        displacement[station] = tbfunx(_TAPER_GRID, values, float(taper_ratio),
                                       lower=0, upper=0)[0]
        per_taper[station] = values

    return {
        'dhb': displacement,
        'per_taper': per_taper,
        'taper_grid': _TAPER_GRID,
        'icase': int(icase),
        'method': 'legacy_sddvc',
    }
