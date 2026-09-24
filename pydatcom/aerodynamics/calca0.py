"""
CALCA0: lifting-surface zero-lift angle of attack.

Starts from the section zero-lift angle and applies two corrections:

- **Twist**, through Figure 4.1.3.1-5's ``DA0OT`` table, interpolated over
  taper ratio, aspect ratio and quarter-chord sweep.
- **Camber**, through a set of seven ragged curves selected by thickness
  ratio and read at ``cos(sweep) * Mach``.

Both corrections use the same bracket-search idiom, which is factored out
here: walk the grid until the query is matched within 2e-2 or falls below a
point, and skip interpolation entirely at an exact hit, at or below the
first point, or above the last.

Reference: datcom-legacy/datcom_2000/calca0.f
"""

import numpy as np
from typing import Dict, Optional, Sequence, Tuple
import logging

from pydatcom.utils.legacy_numeric import tbfunx

logger = logging.getLogger(__name__)

# The source treats a query within this distance of a grid point as exact.
_EXACT_TOLERANCE = 2.0e-2

# DA0OT independent grids.
_TAPER_GRID = (0.0, 0.5, 1.0)
_ASPECT_GRID = (1.5, 3.5, 6.0, 10.0)
_SWEEP_GRID = (-45., -40., -35., -30., -25., -20., -15., -10., -5., 0.,
               5., 10., 15., 20., 25., 30., 35., 40., 45., 50., 55., 60.)

# DA0OT(22,10): 22 sweeps by 10 columns.  Columns are addressed as
# 3*(taper-1) + aspect, so the three taper groups overlap at columns 4 and
# 7 rather than occupying 12 distinct slots.  That compression is in the
# source; see the module tests.
_DA0OT_COLUMNS = 10
_DA0OT = np.array([
    # Column 1
    [-.399, -.3995, -.400, -.400, -.4005, -.4005, -.4005, -.4005, -.4005,
     -.400, -.400, -.400, -.3995, -.3990, -.3985, -.398, -.397, -.3955,
     -.3940, -.385, -.380, -.372],
    # Column 2
    [-.384, -.385, -.3855, -.386, -.386, -.3865, -.386, -.3855, -.3845,
     -.384, -.3825, -.3815, -.380, -.378, -.375, -.372, -.369, -.364,
     -.358, -.350, -.3435, -.335],
    # Column 3
    [-.375, -.375, -.375, -.3745, -.374, -.373, -.372, -.371, -.370,
     -.3685, -.367, -.365, -.362, -.359, -.355, -.3515, -.347, -.342,
     -.336, -.331, -.325, -.318],
    # Column 4
    [-.417, -.4155, -.414, -.413, -.412, -.411, -.4105, -.410, -.4095,
     -.409, -.408, -.4075, -.407, -.4065, -.406, -.405, -.404, -.4025,
     -.401, -.396, -.393, -.387],
    # Column 5
    [-.430, -.427, -.424, -.422, -.420, -.4175, -.4155, -.414, -.412,
     -.410, -.4085, -.407, -.405, -.403, -.401, -.399, -.396, -.393,
     -.390, -.385, -.381, -.375],
    # Column 6
    [-.437, -.434, -.431, -.428, -.4245, -.422, -.419, -.417, -.414,
     -.412, -.409, -.407, -.404, -.402, -.399, -.3965, -.394, -.391,
     -.3885, -.3875, -.386, -.380],
    # Column 7
    [-.419, -.417, -.416, -.414, -.413, -.413, -.412, -.411, -.4105,
     -.410, -.4095, -.409, -.408, -.4075, -.407, -.406, -.405, -.4035,
     -.402, -.400, -.395, -.390],
    # Column 8
    [-.4405, -.4365, -.433, -.430, -.427, -.425, -.422, -.4205, -.419,
     -.4175, -.416, -.415, -.413, -.412, -.410, -.409, -.407, -.406,
     -.406, -.4055, -.405, -.405],
    # Column 9
    [-.456, -.451, -.447, -.442, -.439, -.436, -.433, -.431, -.4285,
     -.426, -.424, -.422, -.420, -.419, -.417, -.416, -.415, -.414,
     -.415, -.416, -.417, -.418],
    # Column 10
    [-.469, -.465, -.460, -.456, -.452, -.449, -.445, -.442, -.439,
     -.437, -.434, -.432, -.429, -.428, -.426, -.425, -.424, -.423,
     -.423, -.423, -.424, -.425],
]).T   # -> shape (22, 10), indexed [sweep, column]

# Camber curves, selected by thickness ratio in percent.  TOC descends.
_CAMBER_TOC = (16., 14., 12., 10., 9., 8., 7.)
_CAMBER_NPT = (7, 6, 4, 4, 3, 4, 5)
_CAMBER_LOCX = (1, 3, 5, 6, 8, 8, 8)
_CAMBER_LOCY = (1, 8, 14, 18, 22, 25, 29)
_CAMBER_CX = (.4, .45, .5, .55, .6, .65, .7, .75, .8, .85, .9, .95)
_CAMBER_CY = (1., .95, .85, .65, .35, -.15, -1.,
              1., .95, .75, .45, -.06, -1.25,
              1., .9, .6, -.5,
              1., 1., .85, .1,
              1., .8, 0.,
              1., .95, .65, 0.,
              1., 1., .95, .85, .55)


def _bracket(grid: Sequence[float], value: float,
             descending: bool = False) -> Tuple[int, bool, float]:
    """The source's shared bracket search.

    Walks the grid until the query matches a point within ``2e-2`` or falls
    past it.  Returns the one-based index the source's ``II`` holds, whether
    interpolation is skipped, and the bracket fraction.

    Skipping happens on an exact hit, at or below the first grid point, and
    when the query lies beyond the last: the loop runs to completion and
    falls into the exact-hit label.
    """
    previous = 0.0
    index = len(grid)
    skip = True

    for position, point in enumerate(grid, start=1):
        index = position
        residual = value - point

        if abs(residual) < _EXACT_TOLERANCE:
            return index, True, 0.0

        past = residual > 0.0 if descending else residual < 0.0
        if past:
            skip = False
            break

        previous = residual
    else:
        # Loop completed without finding a bracket: clamp to the last point.
        return index, True, 0.0

    if index == 1:
        return index, True, 0.0

    fraction = previous / (previous - residual)
    return index, skip, float(fraction)


def _twist_correction(taper_ratio: float, aspect_ratio: float,
                      sweep_c4_deg: float) -> Dict[str, float]:
    """Interpolate DA0OT over taper, aspect ratio and sweep."""
    taper_index, taper_skip, taper_fraction = _bracket(
        _TAPER_GRID, taper_ratio,
    )
    aspect_index, aspect_skip, aspect_fraction = _bracket(
        _ASPECT_GRID, aspect_ratio,
    )
    sweep_index, sweep_skip, sweep_fraction = _bracket(
        _SWEEP_GRID, sweep_c4_deg,
    )

    # The source carries IF(IAR.EQ.4 .AND. A(7).LE.0.5) IAR=3 here. The
    # condition is unreachable: the bracket search only returns index 4
    # for an aspect ratio above 6.0, which cannot also be at or below
    # 0.5. It is kept for fidelity and pinned by a test asserting it
    # never fires, so its absence from the behaviour is deliberate.
    if aspect_index == 4 and aspect_ratio <= 0.5:  # pragma: no cover
        aspect_index = 3

    def column_value(taper_offset: int, aspect: int) -> float:
        column = taper_offset + aspect
        row = sweep_index - 1
        value = _DA0OT[row, column - 1]
        if not sweep_skip:
            below = _DA0OT[row - 1, column - 1]
            value = below + sweep_fraction * (value - below)
        return value

    def at_taper(taper_offset: int, aspect: int) -> float:
        value = column_value(taper_offset, aspect)
        if aspect_skip:
            return value
        lower = column_value(taper_offset, aspect - 1)
        return lower + aspect_fraction * (value - lower)

    taper_offset = 3 * (taper_index - 1)
    result = at_taper(taper_offset, aspect_index)
    if not taper_skip:
        # Second pass at the taper group below, then blend.
        lower_aspect = 3 if aspect_index == 4 else aspect_index
        lower = at_taper(taper_offset - 3, lower_aspect)
        result = lower + taper_fraction * (result - lower)

    return {
        'value': float(result),
        'taper_index': taper_index,
        'aspect_index': aspect_index,
        'sweep_index': sweep_index,
    }


def _camber_factor(thickness_percent: float, cos_sweep_c4: float,
                   mach: float) -> Dict[str, float]:
    """Read the camber curves at cos(sweep)*Mach, blending on thickness."""
    index, skip, fraction = _bracket(
        _CAMBER_TOC, thickness_percent, descending=True,
    )
    compressibility_arg = cos_sweep_c4 * mach

    def curve(position: int) -> float:
        count = _CAMBER_NPT[position - 1]
        x_start = _CAMBER_LOCX[position - 1] - 1
        y_start = _CAMBER_LOCY[position - 1] - 1
        x = np.asarray(_CAMBER_CX[x_start:x_start + count], dtype=float)
        y = np.asarray(_CAMBER_CY[y_start:y_start + count], dtype=float)
        return float(
            tbfunx(x, y, compressibility_arg, lower=-1, upper=0)[0],
        )

    value = curve(index)
    if not skip and index > 1:
        thicker_curve = curve(index - 1)
        value = value + (thicker_curve - value) * fraction

    return {
        'factor': float(value),
        'thickness_index': index,
        'query': float(compressibility_arg),
    }


def calculate_calca0(section_alpha_zero: float,
                     taper_ratio: float, aspect_ratio: float,
                     sweep_c4_deg: float, cos_sweep_c4: float,
                     mach: float,
                     twist_deg: float = 0.0,
                     thickness_ratio: float = 0.0,
                     camber: bool = False) -> Dict[str, object]:
    """Translate CALCA0: lifting-surface zero-lift angle.

    Args:
        section_alpha_zero: Section zero-lift angle, the source's initial
            ``A(134)``.
        taper_ratio: Exposed taper ratio.
        aspect_ratio: Exposed aspect ratio.
        sweep_c4_deg: Quarter-chord sweep, degrees.
        cos_sweep_c4: Cosine of that sweep.
        mach: Free-stream Mach number.
        twist_deg: Wing twist ``TWISTA``.  Twist below half a degree in
            magnitude skips the correction entirely, as the source does.
        thickness_ratio: Section thickness ratio as a fraction of chord.
        camber: Whether the section is cambered.

    Returns:
        Dictionary with ``alpha_zero_lift`` and the two correction terms.

    Notes:
        The ``DA0OT`` table is addressed as ``3*(taper-1) + aspect``, so its
        three taper groups overlap at columns 4 and 7 instead of occupying
        twelve distinct slots.  That compression is in the source and is
        preserved; a test pins the resulting column map.

        The source's low-aspect-ratio column clamp is dead code; see the
        comment at its site.
    """
    alpha_zero = float(section_alpha_zero)
    twist_term = None
    camber_factor = None

    # The source skips the whole twist path below half a degree of twist.
    if abs(twist_deg) >= 0.5:
        twist = _twist_correction(taper_ratio, aspect_ratio, sweep_c4_deg)
        twist_term = twist['value']
        alpha_zero = twist_deg * twist_term + alpha_zero

    if camber:
        factor = _camber_factor(thickness_ratio * 100.0, cos_sweep_c4, mach)
        camber_factor = factor['factor']
        alpha_zero = camber_factor * alpha_zero

    return {
        'alpha_zero_lift': float(alpha_zero),
        'twist_term': twist_term,
        'camber_factor': camber_factor,
        'method': 'legacy_calca0',
    }
