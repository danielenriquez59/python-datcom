"""
CSLOPE: wing surface slope at six chord stations.

Produces the ``WGIN(95)`` to ``WGIN(100)`` slope angles that the downwash
path uses, by differentiating the wing surface coordinates at the 0, 20, 40,
60, 80 and 100 percent chord stations.

Which surface is read depends on where the horizontal tail sits: a tail
above the wing reads the upper surface, a tail at or below it reads the
lower surface with the sign reversed.  Stations the user has already
supplied are left alone, and the whole routine is skipped when there is no
horizontal tail.

Reference: datcom-legacy/datcom_2000/cslope.f
"""

import numpy as np
from typing import Dict, Optional, Sequence
import logging

from pydatcom.utils.constants import RAD
from pydatcom.utils.legacy_numeric import tbfunx

logger = logging.getLogger(__name__)

# The source evaluates exactly these six stations.
_STATIONS = (0.0, 0.2, 0.4, 0.6, 0.8, 1.0)


def calculate_cslope(x_upper: Sequence[float], y_upper: Sequence[float],
                     x_lower: Sequence[float], y_lower: Sequence[float],
                     z_wing: float, z_tail: float,
                     has_horizontal_tail: bool = True,
                     supplied: Optional[Sequence[Optional[float]]] = None
                     ) -> Dict[str, object]:
    """Translate CSLOPE: wing surface slope angles at six chord stations.

    Args:
        x_upper: Upper surface chord stations.
        y_upper: Upper surface ordinates.
        x_lower: Lower surface chord stations.
        y_lower: Lower surface ordinates.
        z_wing: Wing vertical location, ``ZW``.
        z_tail: Horizontal tail vertical location, ``ZH``.
        has_horizontal_tail: When false the source returns immediately and
            every station stays unset.
        supplied: Per-station user values.  A station whose entry is not
            ``None`` is left untouched, matching the source's ``UNUSED``
            test against ``WGIN(95)`` to ``WGIN(100)``.

    Returns:
        Dictionary with ``slope_angles_deg`` and ``slopes``, both six
        elements with NaN where the source leaves a station unset, plus
        which surface was read.

    Raises:
        ValueError: If a needed surface has mismatched or too-short arrays.

    Notes:
        The tail-above test is strict: ``ZH > ZW`` selects the upper
        surface, so a tail exactly level with the wing reads the lower
        surface like one below it.
    """
    stations = np.full(len(_STATIONS), np.nan)
    angles = np.full(len(_STATIONS), np.nan)
    if supplied is None:
        supplied = [None] * len(_STATIONS)
    if len(supplied) != len(_STATIONS):
        raise ValueError("CSLOPE takes six per-station supplied values")

    if not has_horizontal_tail:
        return {
            'slopes': stations,
            'slope_angles_deg': angles,
            'surface': None,
            'method': 'legacy_cslope',
        }

    # ZH > ZW reads the upper surface; level or below reads the lower one.
    tail_above = z_tail > z_wing
    if tail_above:
        x, y, sign, surface = x_upper, y_upper, 1.0, 'upper'
    else:
        x, y, sign, surface = x_lower, y_lower, -1.0, 'lower'

    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    if x.ndim != 1 or x.shape != y.shape or x.size < 2:
        raise ValueError(
            f"CSLOPE needs matching {surface}-surface coordinates with at "
            "least two points")

    for station_slot, station in enumerate(_STATIONS):
        if supplied[station_slot] is not None:
            continue
        slope = sign * tbfunx(x, y, station, lower=0, upper=0)[1]
        stations[station_slot] = slope
        angles[station_slot] = np.arctan(slope) * RAD

    return {
        'slopes': stations,
        'slope_angles_deg': angles,
        'surface': surface,
        'chord_stations': _STATIONS,
        'method': 'legacy_cslope',
    }
