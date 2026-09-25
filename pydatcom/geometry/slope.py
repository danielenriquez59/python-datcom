"""
SLOPE: section lift-curve slope, zero-lift moment and aerodynamic centre.

Uses the Weber pressure distribution that :func:`pydatcom.geometry.ideal`
builds, including compressibility.  Pressure is integrated at two angles of
attack one degree apart; the zero-lift angle is iterated until the lower
angle produces no lift, then the slope follows from the pair.

A viscous correction derived from the trailing-edge angle and the Reynolds
number is applied to the inviscid slope, and at Mach zero the routine also
returns the crest-critical Mach number and its lift coefficient.

Reference: datcom-legacy/datcom_2000/slope.f
"""

import math
import numpy as np
from typing import Dict, Optional
import logging

from pydatcom.utils.constants import PI, DEG

logger = logging.getLogger(__name__)

# The source iterates until the lower angle's lift falls under this.
_LIFT_TOLERANCE = 0.001
# Its guard against a runaway iteration is implicit; this bounds it.
_MAX_ITERATIONS = 50

# The viscous correction is floored here and scaled by this factor.
_CORRECTION_FLOOR = 0.6896
_SLOPE_FACTOR = 1.05

# Reynolds numbers below this are replaced, with a warning, as the source
# does.
_MIN_REYNOLDS = 2.71828e5
_REYNOLDS_SUBSTITUTE = 2.7182e5

# The trailing-edge angle is formed over this chord interval.
_PHITE_INTERVAL = 0.089272624


def _pressure_distribution(alpha_deg: float, weber: Dict, mach: float,
                           beta: float, tmach: float):
    """Upper and lower surface pressure coefficients at one angle."""
    st1, st2, st3, st4, st5 = (weber['st1'], weber['st2'], weber['st3'],
                               weber['st4'], weber['st5'])
    x = weber['x']
    count = len(x)
    last = count - 1

    cos_a = math.cos(alpha_deg * DEG)
    sin_a = math.sin(alpha_deg * DEG)

    cp_upper = np.empty(count)
    cp_lower = np.empty(count)

    for station in range(last):
        cpi = 1.0 - (1.0 + st1[station])**2 / (1.0 + st2[station]**2)
        bo2 = 1.0 - mach**2 * (1.0 - mach * cpi)
        if bo2 < 0.0:
            raise ValueError(
                "SLOPE's local Mach term went negative; the analysis is not "
                "possible for this section at this Mach number")
        bo = math.sqrt(bo2)
        chord = x[station]
        root = math.sqrt((1.0 - chord) / chord) if chord > 0.0 else 0.0

        upper = ((cos_a * (1.0 + st1[station] / bo + st4[station] / beta) +
                  sin_a / beta * (1.0 + st3[station] / bo) * root) /
                 math.sqrt(1.0 + ((st2[station] + st5[station]) / bo)**2))
        lower = ((cos_a * (1.0 + st1[station] / bo - st4[station] / beta) -
                  sin_a / beta * (1.0 + st3[station] / bo) * root) /
                 math.sqrt(1.0 + ((st2[station] - st5[station]) / bo)**2))

        if mach == 0.0:
            cp_upper[station] = 1.0 - upper**2
            cp_lower[station] = 1.0 - lower**2
        else:
            for surface_velocity, surface_label in (
                    (upper, 'upper'), (lower, 'lower')):
                if 1.0 + 0.2 * mach**2 * (1.0 - surface_velocity**2) < 0.0:
                    raise ValueError(
                        f"SLOPE's {surface_label}-surface isentropic term "
                        "went negative; the analysis is not possible")
            cp_upper[station] = tmach * ((1.0 + 0.2 * mach**2 *
                                          (1.0 - upper**2))**3.5 - 1.0)
            cp_lower[station] = tmach * ((1.0 + 0.2 * mach**2 *
                                          (1.0 - lower**2))**3.5 - 1.0)

    # The leading-edge station is handled separately.
    a0 = weber['a0']
    le_velocity = sin_a * ((1.0 + st3[last]) / (0.5 * a0)) if a0 != 0.0 else 0.0
    if tmach == 0.0:
        cp_le = 1.0 - le_velocity**2
    else:
        term = 1.0 + 0.2 * mach**2 * (1.0 - (le_velocity / beta)**2)
        if term < 0.0:
            raise ValueError(
                "SLOPE's leading-edge isentropic term went negative")
        cp_le = tmach * (term**3.5 - 1.0)
    cp_upper[last] = cp_le
    cp_lower[last] = cp_le
    return cp_upper, cp_lower, cos_a, le_velocity


def _integrate(cp_upper, cp_lower, weber: Dict, cos_a: float):
    """Section lift and quarter-chord moment from the pressure difference."""
    theta = weber['theta_nu']
    x = weber['x']
    count = len(x)
    last = count - 1

    load = cp_upper[:last] - cp_lower[:last]
    sin_theta = np.sin(theta[:last]) / 2.0
    cl = float(np.sum(-load * sin_theta)) * PI / (cos_a * count)
    cm = float(np.sum(load * (x[:last] - 0.25) * sin_theta)) * PI / count
    return cl, cm


def calculate_slope(weber: Dict, mach: float, reynolds: float,
                    alpha_zero_lift: float,
                    crest_critical_mach: float = 0.0) -> Dict[str, object]:
    """Translate SLOPE: section CLA, CM0 and aerodynamic centre.

    Args:
        weber: The dictionary :func:`pydatcom.geometry.ideal.calculate_ideal`
            returns, supplying the five coefficient sets and the grid.
        mach: Free-stream Mach number; must be below one.
        reynolds: Section Reynolds number.
        alpha_zero_lift: Starting zero-lift angle, refined by the iteration.
        crest_critical_mach: ``MCC`` from a previous Mach-zero pass; the
            source uses it only to decide whether the supplied Reynolds
            number is usable.

    Returns:
        Dictionary with ``cla`` per degree, ``cm_c4``, ``xac``, the refined
        zero-lift angle, and at Mach zero the crest-critical results.

    Raises:
        ValueError: For Mach at or above one, or when the pressure analysis
            fails, which the source reports and abandons.

    Notes:
        The source computes nothing at all for Mach at or above one, leaving
        ``CLALPA`` at its UNUSED sentinel.  That is an error here rather than
        a silent sentinel.
    """
    if mach >= 1.0:
        raise ValueError(
            "SLOPE is subsonic only; the source leaves CLALPA unset at "
            "Mach 1 and above")

    beta = math.sqrt(1.0 - mach**2)
    tmach = 0.0 if mach == 0.0 else 1.0 / (0.7 * mach**2)

    # Iterate the zero-lift angle until the lower angle carries no lift.
    alpha_zero = float(alpha_zero_lift)
    lift_pair = moment_pair = None
    iterations = 0
    for iterations in range(1, _MAX_ITERATIONS + 1):
        trial_angles = (alpha_zero, alpha_zero + 1.0)
        results = []
        for angle_deg in trial_angles:
            cp_upper, cp_lower, cos_a, _ = _pressure_distribution(
                angle_deg, weber, mach, beta, tmach)
            results.append(_integrate(cp_upper, cp_lower, weber, cos_a))
        lift_pair = [item[0] for item in results]
        moment_pair = [item[1] for item in results]
        if abs(lift_pair[0]) <= _LIFT_TOLERANCE:
            break
        lift_slope = lift_pair[1] - lift_pair[0]
        if lift_slope == 0.0:
            raise ValueError("SLOPE's iteration stalled: the two angles "
                             "produced identical lift")
        alpha_zero = alpha_zero - lift_pair[0] / lift_slope
    else:
        raise ValueError(
            f"SLOPE's zero-lift iteration did not converge in "
            f"{_MAX_ITERATIONS} passes")

    # Viscous correction from the trailing-edge angle and Reynolds number.
    thickness = weber['thickness']
    used_reynolds = reynolds
    if not (reynolds >= 10.0 and crest_critical_mach != 0.0 and mach != 0.0):
        used_reynolds = 1.0e6
    if used_reynolds < _MIN_REYNOLDS:
        logger.warning(
            "Reynolds number too low for the airfoil section module; "
            "section characteristics based on %.4g", _REYNOLDS_SUBSTITUTE)
        used_reynolds = _REYNOLDS_SUBSTITUTE

    phite = ((thickness[5] + thickness[6]) / 2.0 - thickness[1]) / \
        _PHITE_INTERVAL
    exponent = -1.0 + 5.0 * phite / 2.0
    correction = 1.0 - (math.log(used_reynolds / 1.0e5))**exponent * \
        (0.232 + 1.785 * phite - 2.950 * phite**2)
    correction = max(correction, _CORRECTION_FLOOR)

    cla = (lift_pair[1] - lift_pair[0]) * correction * _SLOPE_FACTOR
    cma = moment_pair[1] - moment_pair[0]
    if cla == 0.0:
        raise ValueError("SLOPE produced a zero lift-curve slope")
    xac = 0.25 - cma / cla

    result = {
        'cla': float(cla),
        'cma': float(cma),
        'xac': float(xac),
        'alpha_zero_lift': float(alpha_zero),
        'correction': float(correction),
        'phite': float(phite),
        'reynolds_used': float(used_reynolds),
        'iterations': iterations,
        'method': 'legacy_slope',
    }
    # The source stores the Mach-zero moment as the section CM0.
    if mach == 0.0:
        result['cm_c4'] = float(moment_pair[0])
    return result
