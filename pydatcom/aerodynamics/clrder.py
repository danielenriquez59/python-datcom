"""
CLRDER: rolling moment due to yaw rate, the lateral dynamic derivative CLR.

The first translated piece of the lateral-directional axis.  Three
contributions are summed:

- A lift-dependent term, ``CL * CLR/CL``, built from Figure 7.1.3.2-10 and
  the compressibility correction of DATCOM Section 7.1.3.2.
- A dihedral term, linear in the geometric dihedral angle.
- A twist term from Figure 7.1.3.2-11.

Vertical tail and ventral fin increments are then applied through the
sideslip derivative and the panel's moment arms.

Reference: datcom-legacy/datcom_2000/clrder.f
"""

import numpy as np
from typing import Dict, Optional, Sequence
import logging

from pydatcom.utils.constants import PI, RAD
from pydatcom.utils.legacy_interp import interx

logger = logging.getLogger(__name__)

# Figure 7.1.3.2-10: clrder.f DATA X73210 / Y73210.  LIND=10, so the first
# ten entries are the aspect-ratio grid and the next four the taper grid.
_FIG_71320_10_AR = [1., 2., 3., 4., 5., 6., 7., 8., 9., 10.]
_FIG_71320_10_TAPER = [0., .25, .5, 1.]
_FIG_71320_10_DEP = [
    2., 3.3, 3.9, 4.2, 4.5, 4.8, 4.9, 4.95, 5., 5.,
    3., 4.15, 4.95, 5.4, 5.8, 6., 6.2, 6.4, 6.45, 6.5,
    3.95, 5.5, 6.25, 6.8, 7.15, 7.45, 7.75, 7.8, 7.9, 7.95,
    5., 6.8, 7.5, 8.05, 8.5, 8.85, 9., 9.15, 9.2, 9.25,
]

# Figure 7.1.3.2-11: clrder.f DATA X73211 / Y73211.  LIND=9.  The final
# taper entry is the source sentinel 99., whose column repeats the 0.4
# column, flattening the table above that taper ratio.
_FIG_71320_11_AR = [2., 3., 4., 5., 6., 7., 8., 9., 10.]
_FIG_71320_11_TAPER = [0., .2, .4, 99.]
_FIG_71320_11_DEP = [
    .0014, .00175, .00205, .0022, .00215, .0021, .00215, .00225, .00255,
    .0017, .0022, .0025, .00285, .0032, .00323, .00325, .00335, .0036,
    .002, .00245, .0028, .00315, .0035, .0036, .00375, .00397, .0043,
    .002, .00245, .0028, .00315, .0035, .0036, .00375, .00397, .0043,
]

# The sweep interpolation grid of Section 7.1.3.2 and its two coefficient
# sets, at a fixed 15-degree spacing.
_SWEEP_GRID = [0., 15., 30., 45., 60.]
_UNITS = [.025, .03, .03571428571, .04444444444, .056]
_UNITI = [.05, .06, .07857142857, .1, .12]
_SWEEP_STEP = 15.0


def calculate_clr_wing(cl: Sequence[float], aspect_ratio: float,
                       taper_ratio: float, sweep_c4_deg: float,
                       mach: float, dihedral_deg: float = 0.0,
                       twist_deg: float = 0.0) -> Dict[str, object]:
    """Translate CLRDER's subsonic wing rolling moment due to yaw rate.

    ``CLR = (CL*(CLR/CL) + dCLR/dGamma * Gamma + dCLR/dTwist * twist)/RAD``

    Args:
        cl: Lift coefficient at each angle of attack.
        aspect_ratio: Exposed aspect ratio.
        taper_ratio: Exposed taper ratio.
        sweep_c4_deg: Quarter-chord sweep, degrees.
        mach: Free-stream Mach number, below one.
        dihedral_deg: Geometric dihedral, degrees.
        twist_deg: Wing twist, degrees.

    Returns:
        Dictionary with the per-angle ``clr`` and the three factors behind
        it.

    Raises:
        ValueError: For nonpositive aspect ratio, negative sweep, or a Mach
            number at or above one.

    Notes:
        The source skips the whole lift-dependent path when the quarter-
        chord sweep is negative, leaving CLR at whatever it held.  A forward-
        swept wing is rejected here instead of silently returning stale data.

        Sweep is taken in degrees, which is how DATCOM's input specifies it
        as SAVSI.  The source holds it in radians and converts back with its
        truncated RAD = 57.2957795 for the sweep bracket, so a nominal 15
        degrees arrives as 14.999999997 and falls into the interval below.
        Working in degrees throughout avoids that knife edge; it changes the
        selected interval only for a sweep landing exactly on a grid point.
    """
    cl = np.atleast_1d(np.asarray(cl, dtype=float))
    if aspect_ratio <= 0.0:
        raise ValueError("CLRDER requires a positive aspect ratio")
    if mach >= 1.0:
        raise ValueError("this CLRDER branch is subsonic; Mach must be < 1")
    if sweep_c4_deg < 0.0:
        raise ValueError(
            "CLRDER's lift-dependent path is not defined for a forward-swept "
            "wing; the source leaves CLR untouched in that case")

    sweep_deg = float(sweep_c4_deg)
    sweep_rad = np.deg2rad(sweep_deg)
    cos_sweep = np.cos(sweep_rad)
    tan_sweep = np.tan(sweep_rad)
    beta = np.sqrt(1.0 - (mach * cos_sweep)**2)
    if beta <= 0.0:
        raise ValueError("CLRDER compressibility factor vanished")
    beta_ar = aspect_ratio * beta

    # Section 7.1.3.2 compressibility correction.
    tan_sweep_sq = tan_sweep ** 2
    numerator = (
        1.0
        + aspect_ratio * (1.0 - beta ** 2)
        / (2.0 * beta * (beta_ar + 2.0 * cos_sweep))
        + (beta_ar + 2.0 * cos_sweep) / (beta_ar + 4.0 * cos_sweep)
        * tan_sweep_sq / 8.0
    )
    denominator = (
        1.0
        + (aspect_ratio + 2.0 * cos_sweep)
        / (aspect_ratio + 4.0 * cos_sweep) * tan_sweep_sq / 8.0
    )
    correction = numerator / denominator

    dihedral_factor = (
        PI * aspect_ratio * np.sin(sweep_rad)
        / (12.0 * (aspect_ratio + 4.0 * cos_sweep))
    )

    unit = interx(
        2, [_FIG_71320_10_AR, _FIG_71320_10_TAPER],
        [aspect_ratio, taper_ratio], [10, 4], _FIG_71320_10_DEP,
        lx1l=1, lx2l=1, lx1u=1, lx2u=1,
    )

    # Sweep bracket: last grid point at or below the wing sweep.
    sweep_index = 0
    for edge in _SWEEP_GRID:
        if sweep_deg >= edge:
            sweep_index += 1
    sweep_index = min(sweep_index - 1, len(_SWEEP_GRID) - 2)
    sweep_index = max(sweep_index, 0)

    upper_sweep = _SWEEP_GRID[sweep_index + 1]
    offset_from_upper = upper_sweep - sweep_deg
    clr_cl_zero = (
        (_UNITI[sweep_index + 1] - _UNITI[sweep_index]) / _SWEEP_STEP
        * offset_from_upper
        - _UNITI[sweep_index + 1]
        + (
            (_UNITS[sweep_index + 1] - _UNITS[sweep_index]) / _SWEEP_STEP
            * offset_from_upper
            - _UNITS[sweep_index + 1]
        ) * unit
    )

    twist_factor = interx(
        2, [_FIG_71320_11_AR, _FIG_71320_11_TAPER],
        [aspect_ratio, taper_ratio], [9, 4], _FIG_71320_11_DEP,
    )

    clr_cl = -clr_cl_zero * correction
    clr = (
        cl * clr_cl
        + dihedral_factor * (dihedral_deg / RAD)
        + twist_factor * twist_deg
    ) / RAD

    return {
        'clr': clr,
        'clr_per_cl': float(clr_cl),
        'compressibility_correction': float(correction),
        'dihedral_factor': float(dihedral_factor),
        'twist_factor': float(twist_factor),
        'figure_10_unit': float(unit),
        'sweep_index': sweep_index,
        'method': 'legacy_clrder_wing',
    }


def calculate_clr_panel_increment(alpha_deg: Sequence[float],
                                  cyb_panel: float,
                                  arm_x: float, arm_z: float,
                                  blref: float) -> Dict[str, object]:
    """Translate CLRDER's vertical tail or ventral fin increment.

    ``dCLR = -2*dCYb*(lp*cos a + zp*sin a)*(zp*cos a - lp*sin a)/b^2``
    and the companion roll-damping term
    ``dCLP = 2*dCYb*(zp*cos a - lp*sin a)*(zp*cos a - lp*sin a - zp)/b^2``.

    Args:
        alpha_deg: Angle-of-attack schedule, degrees.
        cyb_panel: The panel's sideslip side-force derivative.
        arm_x: Longitudinal arm from the CG, the source's ``LP``.
        arm_z: Vertical arm from the CG, ``ZP``.
        blref: Lateral reference length, ``BLREF``.

    Returns:
        Dictionary with per-angle ``dclr`` and ``dclp``.

    Raises:
        ValueError: If the reference length is nonpositive.
    """
    if blref <= 0.0:
        raise ValueError("CLRDER requires a positive lateral reference length")
    alpha = np.atleast_1d(np.asarray(alpha_deg, dtype=float))
    sin_alpha = np.sin(alpha / RAD)
    cos_alpha = np.cos(alpha / RAD)
    blref_sq = blref ** 2

    arm_along_wind = arm_x * cos_alpha + arm_z * sin_alpha
    arm_across_wind = arm_z * cos_alpha - arm_x * sin_alpha

    return {
        'dclr': -2.0 * cyb_panel * arm_along_wind * arm_across_wind / blref_sq,
        'dclp': (
            2.0 * cyb_panel * arm_across_wind * (arm_across_wind - arm_z)
            / blref_sq
        ),
        'method': 'legacy_clrder_panel',
    }
