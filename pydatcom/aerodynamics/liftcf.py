"""
LIFTCF: lifting-surface lift and normal force at every angle of attack.

Runs after WTLIFT, whose slope, maximum lift and stall angle it consumes,
and fills the lift curve ``AOUT(21)`` onward and the normal-force curve
``AOUT(61)`` onward.  Each planform type has its own method:

- **Straight tapered** (Section 4.1.3.3): below the stall angle, a
  normal-force slope built from ``B(45)`` and the Figure 4.1.3.3-55A
  increment; above it, a fairing toward the Figure 4.1.3.3-55B value at 90
  degrees.  CN is ``(CNa*cos + CNaa*|sin|)*sin`` and CL is ``CN*cos``.
- **Curved**: Figure 4.1.3.3-58, scaled by ``sqrt(b/2 / length)``.
- **Double delta**: Figure 4.1.3.3-56, the solid curve, or the dashed
  low-aspect-ratio curve near Mach one below 12 degrees.
- **Cranked**: linear up to the Figure 4.1.3.3-57 break angle, then the
  Figure 4.1.3.3-56 increment above it.

The tables were extracted from the source DATA statements by parsing, with
``tools/fortran_data.py``, rather than by hand.

Reference: datcom-legacy/datcom_2000/liftcf.f
"""

import math
import numpy as np
from typing import Dict, Optional, Sequence
import logging

from pydatcom.aerodynamics.cdrag import (
    CRANKED, CURVED, DOUBLE_DELTA, STRAIGHT_TAPERED,
)
from pydatcom.aerodynamics.wtlift import calculate_clmxbs
from pydatcom.utils.constants import DEG, RAD
from pydatcom.utils.legacy_numeric import angles, tbfunx
from pydatcom.utils.legacy_tables import tlinex, tlin3x

logger = logging.getLogger(__name__)

# Figure 4.1.3.3-55A: the CN_alpha increment below the stall angle.  DC is
# DCNAA(7,12): taper ratio fastest, one run of seven per J value.
_AJ = np.array([
    -4.0, -2.0, -1.0, 0.0, 0.15, 0.25, 0.5,
    1.0, 1.5, 2.0, 2.5, 3.0,
])
_TRAT = np.array([
    0.2, 0.4, 0.46, 0.5, 0.55, 0.6, 1.0,
])
_DC = np.array([
    -0.9, -0.676, -0.607, -0.564, -0.506, -0.452, 0.0,
    -0.4, -0.3, -0.27, -0.25, -0.225, -0.2, 0.0,
    -0.1, -0.075, -0.065, -0.0625, -0.0562, -0.05, 0.0,
    0.25, 0.25, 0.225, 0.2083, 0.1874, 0.1666, 0.0,
    0.6, 0.6, 0.54, 0.5, 0.45, 0.4, 0.0,
    0.75, 0.75, 0.675, 0.625, 0.562, 0.5, 0.0,
    1.1, 1.1, 0.99, 0.917, 0.825, 0.734, 0.0,
    1.6, 1.6, 1.44, 1.334, 1.2, 1.068, 0.0,
    2.0, 2.0, 2.0, 1.852, 1.666, 1.48, 0.0,
    2.35, 2.35, 2.35, 2.35, 2.115, 1.88, 0.0,
    2.65, 2.65, 2.65, 2.65, 2.65, 2.36, 0.0,
    2.85, 2.85, 2.85, 2.85, 2.85, 2.85, 0.0,
])

# Figure 4.1.3.3-55B: CN_alpha at 90 degrees.  AR, AI and TIR are
# EQUIVALENCEd onto one array, so both of its abscissas are the TIR grid.
_TIR = np.array([
    0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0,
])
_C90 = np.array([
    2.0, 1.32, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2,
])
_C90I = np.array([
    2.0, 1.4, 1.26, 1.21, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2,
])

# The stall-parameter factor D, over the same TIR grid.
_D = np.array([
    0.0, -0.35, -0.7, -1.0, -1.27, -1.46, -1.55, -1.4, -1.0, -0.21, 0.0,
])

# Figure 4.1.3.3-58: curved planforms, CL over sqrt(b/2) against angle.
_A58 = np.array([
    0.0, 2.0, 4.0, 6.0, 8.0, 10.0, 12.0, 14.0, 16.0, 18.0, 20.0,
])
_CLJ58 = np.array([
    0.0, 0.11, 0.238, 0.375, 0.521, 0.68, 0.85, 1.025, 1.205, 1.39, 1.575,
])

# Figure 4.1.3.3-56 solid: angle of attack and beta*tan(LE sweep).  Y13356
# is (12,6), the beta*tan axis fastest.
_X13356 = np.array([
    0.0, 4.0, 8.0, 12.0, 16.0, 20.0,
])
_X23356 = np.array([
    0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0,
])
_Y13356 = np.array([
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.165, 0.157, 0.12, 0.095, 0.082, 0.074, 0.066, 0.061, 0.058, 0.055, 0.055, 0.055,
    0.297, 0.276, 0.227, 0.19, 0.167, 0.151, 0.142, 0.135, 0.13, 0.127, 0.125, 0.125,
    0.425, 0.409, 0.35, 0.3, 0.268, 0.25, 0.235, 0.225, 0.217, 0.212, 0.209, 0.208,
    0.568, 0.545, 0.48, 0.42, 0.377, 0.348, 0.327, 0.311, 0.301, 0.296, 0.293, 0.293,
    0.701, 0.675, 0.602, 0.533, 0.482, 0.443, 0.422, 0.404, 0.392, 0.385, 0.381, 0.38,
])

# Figure 4.1.3.3-56 dashed: aspect ratio.  Y33356 is (8,4,5) by
# declaration and read as (8,5,4): beta*tan fastest, then aspect ratio,
# then the first four angles.
_X33356 = np.array([
    1.0, 1.5, 2.0, 2.5, 3.0,
])
_Y33356 = np.array([
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.07, 0.07, 0.07, 0.07, 0.07, 0.07, 0.066, 0.061,
    0.087, 0.087, 0.087, 0.085, 0.081, 0.074, 0.066, 0.061,
    0.117, 0.113, 0.105, 0.095, 0.082, 0.074, 0.066, 0.061,
    0.145, 0.142, 0.12, 0.095, 0.082, 0.074, 0.066, 0.061,
    0.165, 0.157, 0.12, 0.095, 0.082, 0.074, 0.066, 0.061,
    0.145, 0.145, 0.145, 0.145, 0.145, 0.142, 0.138, 0.135,
    0.18, 0.18, 0.178, 0.172, 0.163, 0.151, 0.142, 0.135,
    0.215, 0.214, 0.206, 0.19, 0.167, 0.155, 0.142, 0.135,
    0.285, 0.273, 0.227, 0.19, 0.167, 0.151, 0.142, 0.135,
    0.297, 0.276, 0.227, 0.19, 0.167, 0.151, 0.142, 0.135,
    0.425, 0.409, 0.35, 0.3, 0.268, 0.25, 0.235, 0.225,
    0.425, 0.409, 0.35, 0.3, 0.268, 0.25, 0.235, 0.225,
    0.425, 0.409, 0.35, 0.3, 0.268, 0.25, 0.235, 0.225,
    0.425, 0.409, 0.35, 0.3, 0.268, 0.25, 0.235, 0.225,
    0.425, 0.409, 0.35, 0.3, 0.268, 0.25, 0.235, 0.225,
])

# Figure 4.1.3.3-57: the cranked-wing break angle.
_X13357 = np.array([
    0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0,
])
_Y13357 = np.array([
    7.0, 7.0, 7.0, 7.0, 7.0, 6.9, 6.75, 6.5, 6.2, 5.8, 5.3,
])


_FIG_55A = _DC.reshape(12, 7).T
_FIG_56_SOLID = _Y13356.reshape(6, 12).T
_FIG_56_DASHED = _Y33356.reshape(4, 5, 8).transpose(2, 1, 0)


def _fig56(alpha: float, btanle: float, aspect_ratio: float,
           dashed: bool, solid_points: int = 6) -> float:
    """Figure 4.1.3.3-56, solid or dashed, as the double-delta and cranked
    paths call it.  ``solid_points`` is the solid call's ``NX1``."""
    if dashed:
        return tlin3x(
            _X33356, _X23356[:8], _X13356[:4], _FIG_56_DASHED,
            aspect_ratio, btanle, alpha, 2, 2, 2, 2, 1, 2,
        )
    return tlinex(
        _X13356[:solid_points], _X23356, _FIG_56_SOLID[:, :solid_points],
        alpha, btanle, 2, 2, 2, 1,
    )


def _straight(alpha, geometry, section, lift, flight, sref, state):
    """Labels 1000 to 1070, the straight tapered method."""
    area = float(geometry['area'])
    aspect_ratio = float(geometry['aspect_ratio'])
    factor = float(geometry['arclss_factor']) + 1.0
    cla = float(lift['cla'])
    clmax = float(lift['clmax'])
    alpha_clmax = float(lift['alpha_clmax'])
    alpha_zero = float(flight['alpha_zero_lift'])
    beta = float(flight['beta'])
    a160 = float(geometry['a160'])

    ajay = (.3 * factor * aspect_ratio * float(geometry['cos_le']) *
            (factor * (float(geometry['a159']) + 1.0) - (a160 / 7.0)**3) /
            beta)
    stall = angles(1, [alpha_clmax - alpha_zero] + state[1:6])
    local = list(state[6:12])
    cna = cla * RAD
    b45 = (clmax / (stall[2] * stall[3]) - cna * stall[3]) / abs(stall[2])

    low_aspect_ratio = aspect_ratio < float(geometry['arclss_ratio'])
    clsmax = clmax if low_aspect_ratio else None
    cn, cl = [], []

    for a in alpha:
        if a >= alpha_clmax:
            if clsmax is None:
                # The first angle past the stall on a high-aspect-ratio
                # surface evaluates CLMXBS once at BU4 = 4/beta.
                clsmax = calculate_clmxbs(
                    4.0 / beta, float(flight['mach']), a160,
                    float(section['deltay']), float(section['xovc']),
                    area, sref)['clmax']
            temp = alpha_zero / (90.0 - alpha_clmax)
            local = angles(1, [a * (1.0 + temp) - 90.0 * temp] + local[1:])
            angle_ratio = abs(stall[4] / local[4])
            dj, _ = tbfunx(_TIR, _D, angle_ratio, 0, 0)
            if aspect_ratio > 1.0:
                cnaa90, _ = tbfunx(_TIR, _C90I, 1.0 / aspect_ratio, 0, 0)
                cnaa90 = cnaa90 * area / sref
            else:
                cnaa90, _ = tbfunx(_TIR, _C90, aspect_ratio, 0, 0)
            cnaaj = (
                b45 + (cnaa90 - b45) * (1.0 - angle_ratio)
                + RAD * dj * cla / 2.3 * (beta * clmax / clsmax) ** 2
            )
        else:
            local = angles(1, [a - alpha_zero] + local[1:])
            angle_ratio = abs(local[4] / stall[4])
            increment = tlinex(
                _AJ, _TRAT, _FIG_55A, ajay, angle_ratio, 2, 1, -1, 0,
            )
            cnaaj = b45 + increment * area / sref
        normal = (RAD * cla * local[3] + cnaaj * abs(local[2])) * local[2]
        cn.append(normal)
        cl.append(normal * math.cos(DEG * a))
    return {
        'cl': np.array(cl), 'cn': np.array(cn), 'b45': float(b45),
        'ajay': float(ajay), 'low_aspect_ratio': bool(low_aspect_ratio),
        'angle_state': stall + local,
    }


def calculate_liftcf(planform_type: float,
                     alpha_deg: Sequence[float],
                     geometry: Dict[str, float],
                     section: Dict[str, float],
                     lift: Dict[str, float],
                     flight: Dict[str, float],
                     sref: float,
                     angle_state: Optional[Sequence[float]] = None
                     ) -> Dict[str, object]:
    """Translate LIFTCF: CL and CN of one lifting surface at every angle.

    Args:
        planform_type: ``AIN(15)``; one of the ``cdrag`` WTYPE constants.
        alpha_deg: ``B(23)`` onward, the surface's local angles.
        geometry: ``A`` block entries.  Straight tapered: ``area`` A(3),
            ``aspect_ratio`` A(7), ``cos_le`` A(37), ARCLSS's
            ``arclss_factor`` A(123) and ``arclss_ratio`` A(125), and
            WTLIFT's ``a159`` and ``a160``.  Double delta and cranked:
            ``area``, ``aspect_ratio``, ``inboard_span`` A(23),
            ``aspect_ratio_inboard`` A(5) and ``tan_le`` A(62).  Curved:
            ``area`` and ``planform_length`` A(29).
        section: ``deltay`` AIN(17) and ``xovc`` AIN(18) for the straight
            method; ``sspne`` AIN(3) for the others.
        lift: WTLIFT's ``cla`` AOUT(101), ``clmax`` B(44) and
            ``alpha_clmax`` B(43).  The curved method needs none of them.
        flight: ``mach`` B(1), ``beta`` B(2) and ``alpha_zero_lift`` B(49).
        sref: ``SREF``.
        angle_state: ``A(147)`` to ``A(158)``, the two ANGLES records the
            straight method keeps, as the previous call left them.  Zero
            by default, as in a freshly cleared COMMON block.  See
            :func:`pydatcom.utils.legacy_numeric.angles`.

    Returns:
        Dictionary with ``cl`` (AOUT(21) onward) and ``cn`` (AOUT(61)
        onward).  The straight method adds ``b45`` and the final
        ``angle_state``; the cranked method adds ``alpha_break``.

    Raises:
        ValueError: For a planform type outside the four the source knows.

    Notes:
        Above the stall, ``CNAA90`` from Figure 4.1.3.3-55B is put on the
        SREF basis with ``A(3)/SREF`` when the aspect ratio exceeds one but
        not when it is one or less, although every other term it is combined
        with is on the SREF basis on both branches.  The source form is kept;
        the two agree only when SREF equals the exposed area.

        The cranked method's break-angle lookup passes the solid curve with
        ``NX1=5`` where every other call passes 6, dropping the 20-degree
        column.  Figure 4.1.3.3-57 keeps the break at or below 7 degrees
        across its grid, so the column is never reached; the source form
        is kept.
    """
    kind = float(planform_type)
    alpha = [float(a) for a in alpha_deg]
    if angle_state is None:
        state = [0.0] * 12
    else:
        state = [float(v) for v in angle_state]
    if len(state) != 12:
        raise ValueError("LIFTCF's angle state is the twelve words "
                         "A(147) to A(158)")

    if kind == STRAIGHT_TAPERED:
        result = _straight(alpha, geometry, section, lift, flight,
                           float(sref), state)
        result['method'] = 'legacy_liftcf_straight'
        return result

    area = float(geometry['area'])
    if kind == CURVED:
        root = math.sqrt(float(section['sspne']) /
                         float(geometry['planform_length']))
        cl = np.array([root * tbfunx(_A58, _CLJ58, a, 0, 1)[0] * area / sref
                       for a in alpha])
        cn = cl / np.cos(DEG * np.array(alpha))
        return {'cl': cl, 'cn': cn, 'method': 'legacy_liftcf_curved'}

    if kind not in (DOUBLE_DELTA, CRANKED):
        raise ValueError(f"LIFTCF: method not applicable to type {kind:g}")

    aspect_ratio = float(geometry['aspect_ratio'])
    cla = float(lift['cla'])
    con = (float(geometry['inboard_span']) /
           (float(section['sspne']) * float(geometry['aspect_ratio_inboard']))
           * RAD)
    tan_le = float(geometry['tan_le'])
    btanle = float(flight['beta']) * tan_le
    dashed_region = not (aspect_ratio >= 3.0 or float(flight['mach']) < 0.7
                         or btanle > 7.0)

    if kind == DOUBLE_DELTA:
        cl = []
        for a in alpha:
            magnitude = abs(a)
            dashed = dashed_region and magnitude < 12.0
            value = _fig56(magnitude, btanle, aspect_ratio, dashed) * cla * con
            cl.append(-value if a < 0.0 else value)
        cl = np.array(cl)
        cn = cl / np.cos(DEG * np.array(alpha))
        return {'cl': cl, 'cn': cn, 'dashed_region': dashed_region,
                'method': 'legacy_liftcf_double_delta'}

    alpha_break, _ = tbfunx(_X13357, _Y13357, 1.0 / tan_le, 0, 2)
    cl_break = cla * alpha_break
    dashed = dashed_region and alpha_break < 12.0
    cl_break_nonlinear = _fig56(alpha_break, btanle, aspect_ratio, dashed,
                                solid_points=5) * cla * con
    cl = []
    for a in alpha:
        if a <= alpha_break:
            cl.append(cla * a)
            continue
        dashed = dashed_region and a < 12.0
        value = _fig56(a, btanle, aspect_ratio, dashed) * cla * con
        cl.append(cl_break + value - cl_break_nonlinear)
    cl = np.array(cl)
    cn = cl / np.cos(DEG * np.array(alpha))
    return {'cl': cl, 'cn': cn, 'alpha_break': float(alpha_break),
            'dashed_region': dashed_region, 'method': 'legacy_liftcf_cranked'}
