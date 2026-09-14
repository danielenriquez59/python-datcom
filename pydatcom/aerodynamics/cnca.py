"""
CNCA: normal and axial force coefficients, and the slope pass.

The source applies the same three operations to each of nine component
configurations -- body, wing, horizontal tail, and the wing-body,
body-tail, body-vertical, wing-body-tail and wing-body-vertical-tail
buildups:

- rotate ``CL`` and ``CD`` into ``CN`` and ``CA``
- differentiate ``CL`` over the angle schedule to give ``CLA``
- for the wing only, differentiate ``CM`` to give ``CMA``

The rotation itself already exists in ``pydatcom.aerodynamics.moment`` as
``calculate_normal_force_coefficient`` and
``calculate_axial_force_coefficient``, and agrees with the source.  What
this module adds is the derivative pass, the ``UNUSED`` masking that skips
components with no data, and the source's rule that no ``CLA`` is produced
at the first angle.

Reference: datcom-legacy/datcom_2000/cnca.f
"""

import numpy as np
from typing import Dict, Optional, Sequence
import logging

from pydatcom.utils.constants import RAD, UNUSED
from pydatcom.utils.legacy_numeric import tbfunx

logger = logging.getLogger(__name__)


def _is_set(values: np.ndarray) -> np.ndarray:
    """Which entries hold real data rather than the legacy sentinel."""
    return np.abs(np.abs(values) - UNUSED) > 1.0e-40


def calculate_cnca(alpha_deg: Sequence[float], cl: Sequence[float],
                   cd: Sequence[float],
                   cm: Optional[Sequence[float]] = None
                   ) -> Dict[str, np.ndarray]:
    """Translate CNCA for one component configuration.

    ``CN = CL*cos(a) + CD*sin(a)`` and ``CA = CD*cos(a) - CL*sin(a)``, with
    ``CLA`` and optionally ``CMA`` taken as TBFUNX slopes over the angle
    schedule.

    Args:
        alpha_deg: Angle-of-attack schedule, degrees.
        cl: Lift coefficient at each angle.  Entries equal to the legacy
            ``UNUSED`` sentinel are treated as absent.
        cd: Drag coefficient at each angle, same convention.
        cm: Pitching moment at each angle.  The source computes ``CMA`` for
            the wing only, so this is optional.

    Returns:
        Dictionary with ``cn``, ``ca`` and ``cla``, plus ``cma`` when ``cm``
        is supplied.  Entries the source would leave untouched come back as
        NaN rather than a stale value.

    Raises:
        ValueError: If the arrays do not match or are empty.

    Notes:
        The source guards its ``CLA`` call with ``J .GE. 2``, so no lift
        slope is produced at the first angle of the schedule.  That is
        preserved: ``cla[0]`` is NaN.  The wing ``CMA`` call carries no such
        guard and is evaluated at every angle, an asymmetry that is in the
        source rather than introduced here.
    """
    alpha = np.atleast_1d(np.asarray(alpha_deg, dtype=float))
    lift = np.atleast_1d(np.asarray(cl, dtype=float))
    drag = np.atleast_1d(np.asarray(cd, dtype=float))
    if alpha.size == 0 or lift.shape != alpha.shape or drag.shape != alpha.shape:
        raise ValueError("CNCA requires matching nonempty angle, CL and CD arrays")

    sin_a = np.sin(alpha / RAD)
    cos_a = np.cos(alpha / RAD)

    lift_set = _is_set(lift)
    drag_set = _is_set(drag)
    both = lift_set & drag_set

    normal = np.full(alpha.shape, np.nan)
    axial = np.full(alpha.shape, np.nan)
    normal[both] = lift[both] * cos_a[both] + drag[both] * sin_a[both]
    axial[both] = drag[both] * cos_a[both] - lift[both] * sin_a[both]

    # CLA: a TBFUNX slope over the schedule, skipped at the first angle.
    slope = np.full(alpha.shape, np.nan)
    if alpha.size >= 2 and np.all(np.diff(alpha) > 0) and np.all(lift_set):
        for index in range(1, alpha.size):
            slope[index] = tbfunx(alpha, lift, alpha[index],
                                  lower=0, upper=0)[1]

    result = {
        'cn': normal,
        'ca': axial,
        'cla': slope,
        'method': 'legacy_cnca',
    }

    if cm is not None:
        moment = np.atleast_1d(np.asarray(cm, dtype=float))
        if moment.shape != alpha.shape:
            raise ValueError("CNCA CM array must match the angle schedule")
        moment_slope = np.full(alpha.shape, np.nan)
        # The source also tests element 2 of the array before proceeding.
        if (alpha.size >= 2 and np.all(np.diff(alpha) > 0) and
                np.all(_is_set(moment))):
            for index in range(alpha.size):
                moment_slope[index] = tbfunx(alpha, moment, alpha[index],
                                             lower=0, upper=0)[1]
        result['cma'] = moment_slope

    return result
