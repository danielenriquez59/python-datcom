"""
WTLIFT and CLMXBS: lifting-surface lift-curve slope and maximum lift.

``WTLIFT`` runs once per lifting surface before the lift curve is built,
and fills the quantities the rest of the program reads:

- ``AOUT(101)``: the lift-curve slope on the SREF basis, from the
  Section 4.1.3.2 Helmbold-Diederich form with the Figure 4.1.3.2-52 ratio
  for cranked wings.
- ``B(44)`` and ``B(43)``: the maximum lift coefficient and the angle it
  occurs at.  High-aspect-ratio surfaces take Figures 4.1.3.4-21A/B and the
  Figure 4.1.3.4-22 Mach correction; low-aspect-ratio surfaces take
  ``CLMXBS`` (Figures 4.1.3.4-23 and 24A) and Figures 4.1.3.4-25A/B.
- ``A(144)``, ``A(145)``, ``A(146)``, ``A(159)``, ``A(160)``, and for
  non-straight planforms the panel slopes ``A(171)`` and ``A(172)``.

``CLMXBS`` is called only from here and from LIFTCF, so it lives with its
caller.

The tables were extracted from the source DATA statements by parsing, with
``tools/fortran_data.py``, rather than by hand.

Reference: datcom-legacy/datcom_2000/wtlift.f, clmxbs.f
"""

import numpy as np
from typing import Dict
import logging

from pydatcom.aerodynamics.cdrag import CRANKED, CURVED, STRAIGHT_TAPERED
from pydatcom.utils.constants import DEG, PI
from pydatcom.utils.legacy_numeric import tbfunx
from pydatcom.utils.legacy_tables import tlinex, tlin3x

logger = logging.getLogger(__name__)

# CLMXBS: Figure 4.1.3.4-23, base CLMAX.  C1ABC, DYAG, and CBASE declared
# BASE(19,12): the first six columns are part (A), the last six part (B).

_C1ABC = np.array([
    0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8,
    2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4, 3.6,
])

_DYAG = np.array([
    0.0, 0.25, 0.5, 0.75, 1.0, 1.35,
])

_CBASE = np.array([
    0.9, 1.38, 1.58, 1.65, 1.65, 1.58, 1.45, 1.32, 1.21, 1.13,
    1.06, 1.02, 0.99, 0.96, 0.94, 0.92, 0.91, 0.9, 0.9, 0.75,
    1.28, 1.47, 1.55, 1.56, 1.5, 1.38, 1.27, 1.17, 1.08, 1.02,
    0.98, 0.94, 0.92, 0.9, 0.89, 0.88, 0.87, 0.87, 0.65, 1.21,
    1.39, 1.47, 1.47, 1.4, 1.27, 1.17, 1.08, 1.02, 0.97, 0.94,
    0.92, 0.9, 0.88, 0.87, 0.86, 0.85, 0.84, 0.6, 1.11, 1.3,
    1.37, 1.38, 1.33, 1.22, 1.13, 1.06, 1.01, 0.96, 0.92, 0.9,
    0.88, 0.86, 0.85, 0.84, 0.83, 0.83, 0.4, 1.02, 1.2, 1.29,
    1.31, 1.26, 1.18, 1.11, 1.04, 0.98, 0.94, 0.9, 0.88, 0.86,
    0.84, 0.83, 0.82, 0.82, 0.81, 0.4, 1.02, 1.2, 1.29, 1.31,
    1.26, 1.18, 1.09, 1.02, 0.96, 0.91, 0.87, 0.85, 0.84, 0.82,
    0.81, 0.8, 0.8, 0.8, 0.68, 1.18, 1.37, 1.43, 1.44, 1.37,
    1.24, 1.17, 1.1, 1.06, 1.02, 0.99, 0.96, 0.94, 0.92, 0.91,
    0.9, 0.89, 0.88, 0.65, 1.11, 1.29, 1.36, 1.37, 1.31, 1.2,
    1.13, 1.07, 1.03, 0.99, 0.96, 0.94, 0.92, 0.9, 0.89, 0.88,
    0.87, 0.86, 0.6, 1.07, 1.23, 1.29, 1.3, 1.25, 1.16, 1.1,
    1.05, 1.01, 0.97, 0.94, 0.92, 0.9, 0.88, 0.87, 0.86, 0.85,
    0.84, 0.55, 0.99, 1.15, 1.23, 1.25, 1.21, 1.14, 1.09, 1.04,
    0.99, 0.95, 0.92, 0.9, 0.88, 0.86, 0.85, 0.84, 0.83, 0.83,
    0.5, 0.91, 1.08, 1.17, 1.2, 1.17, 1.13, 1.07, 1.03, 0.98,
    0.94, 0.9, 0.88, 0.86, 0.84, 0.83, 0.82, 0.82, 0.82, 0.5,
    0.91, 1.08, 1.17, 1.2, 1.17, 1.11, 1.06, 1.01, 0.96, 0.91,
    0.88, 0.85, 0.83, 0.82, 0.81, 0.8, 0.8, 0.8,
])


# CLMXBS: Figure 4.1.3.4-24A, the CLMAX increment.  DE declared DV(15,3).

_C2A = np.array([
    0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 13.0, 14.0,
])

_AMN = np.array([
    0.2, 0.4, 0.6,
])

_DE = np.array([
    -0.11, -0.1, -0.085, -0.045, 0.0, 0.055, 0.115, 0.175, 0.23, 0.285, 0.335, 0.365, 0.36, 0.315, 0.225,
    -0.11, -0.1, -0.085, -0.045, 0.0, 0.05, 0.1, 0.15, 0.2, 0.255, 0.3, 0.33, 0.32, 0.275, 0.185,
    -0.11, -0.1, -0.085, -0.045, 0.0, 0.04, 0.08, 0.115, 0.15, 0.18, 0.21, 0.225, 0.21, 0.16, 0.08,
])


# WTLIFT: Figure 4.1.3.4-24B, the C2 factor.

_TR = np.array([
    0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0,
])

_C2 = np.array([
    0.0, 0.21, 0.5, 0.9, 1.08, 1.05, 1.0, 0.94, 0.9, 0.86, 0.85,
])


# WTLIFT: Figure 4.1.3.4-21A, CLMAX/cl_max.  CLL declared CLOCL(13,7).

_SALE = np.array([
    0.0, 5.0, 10.0, 15.0, 20.0, 25.0, 30.0, 35.0, 40.0, 45.0, 50.0, 55.0, 60.0,
])

_DELTAY = np.array([
    1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.5,
])

_CLL = np.array([
    0.9, 0.915, 0.93, 0.95, 0.975, 1.0, 1.03, 1.07, 1.1, 1.15, 1.19, 1.25, 1.3,
    0.9, 0.91, 0.92, 0.94, 0.96, 0.98, 1.0, 1.02, 1.05, 1.08, 1.11, 1.15, 1.19,
    0.9, 0.91, 0.92, 0.93, 0.94, 0.95, 0.96, 0.965, 0.975, 0.99, 1.0, 1.01, 1.03,
    0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.89, 0.89, 0.885, 0.88, 0.88, 0.875, 0.87,
    0.9, 0.895, 0.89, 0.88, 0.87, 0.86, 0.85, 0.835, 0.82, 0.8, 0.78, 0.755, 0.73,
    0.9, 0.89, 0.88, 0.87, 0.85, 0.83, 0.81, 0.79, 0.76, 0.725, 0.695, 0.65, 0.59,
    0.9, 0.89, 0.87, 0.86, 0.84, 0.82, 0.795, 0.77, 0.73, 0.7, 0.65, 0.59, 0.52,
])


# WTLIFT: Figure 4.1.3.4-21B, the stall-angle increment.  DACLL declared
# DACLMX(13,4).

_DY = np.array([
    1.2, 2.0, 3.0, 4.0,
])

_DACLL = np.array([
    1.75, 1.9, 2.2, 2.7, 3.4, 4.15, 5.1, 6.1, 7.3, 8.7, 10.15, 11.75, 13.3,
    0.1, 0.5, 1.05, 1.65, 2.3, 3.1, 3.9, 4.7, 5.7, 6.7, 7.7, 8.75, 9.8,
    1.2, 1.4, 1.7, 2.0, 2.4, 2.85, 3.35, 3.7, 4.25, 4.7, 5.3, 5.9, 6.65,
    2.2, 2.1, 2.0, 2.0, 2.1, 2.15, 2.3, 2.4, 2.55, 2.7, 2.9, 3.05, 3.3,
])


# WTLIFT: Figure 4.1.3.4-22, the Mach correction to CLMAX.  DCAR declared
# DCLTB(5,6,4).

_DYA = np.array([
    2.0, 2.25, 2.5, 3.0, 4.0, 4.5,
])

_AMACH = np.array([
    0.2, 0.3, 0.4, 0.5, 0.6,
])

_SALE4 = np.array([
    0.0, 20.0, 40.0, 60.0,
])

_DCAR = np.array([
    0.0, -0.02, -0.02, 0.0, 0.0, 0.0, -0.13, -0.19, -0.2, -0.19,
    0.0, -0.185, -0.32, -0.4, -0.445, 0.0, -0.21, -0.36, -0.45, -0.5,
    0.0, -0.24, -0.42, -0.545, -0.64, -0.0, -0.24, -0.455, -0.605, -0.72,
    0.0, -0.04, -0.045, -0.02, 0.0, 0.0, -0.095, -0.165, -0.22, -0.25,
    0.0, -0.105, -0.19, -0.26, -0.29, 0.0, -0.125, -0.225, -0.3, -0.345,
    0.0, -0.15, -0.265, -0.37, -0.45, 0.0, -0.15, -0.29, -0.41, -0.51,
    0.0, -0.04, -0.07, -0.095, -0.1, 0.0, -0.043, -0.08, -0.108, -0.123,
    0.0, -0.045, -0.085, -0.12, -0.145, 0.0, -0.05, -0.095, -0.14, -0.19,
    0.0, -0.07, -0.125, -0.19, -0.26, 0.0, -0.07, -0.14, -0.215, -0.3,
    0.0, 0.0, 0.0, 0.0, -0.02, 0.0, 0.0, 0.0, 0.0, -0.022,
    0.0, 0.0, 0.0, 0.0, -0.03, 0.0, 0.0, 0.0, -0.02, -0.07,
    0.0, 0.0, -0.02, -0.055, -0.085, 0.0, -0.04, -0.085, -0.14, -0.2,
])


# WTLIFT: Figure 4.1.3.4-25A, the low-aspect-ratio stall angle.

_C1ABCS = np.array([
    0.0, 0.4, 0.8, 1.2, 1.6, 2.0, 2.4, 2.8, 3.2, 3.6,
])

_ACLMX = np.array([
    35.0, 35.0, 35.0, 32.0, 28.0, 25.0, 23.2, 22.0, 21.5, 21.0,
])


# WTLIFT: Figure 4.1.3.4-25B, its increment above and below C1 = 4.5.
# DACL is 20 by 3 and DACLO 10 by 10, the first grid fastest in each.

_DMN = np.array([
    0.2, 0.4, 0.6,
])

_C1TAB = np.array([
    4.5, 5.0, 5.5, 6.0, 6.5, 7.0, 7.5, 8.0, 8.5, 9.0,
    9.5, 10.0, 10.5, 11.0, 11.5, 12.0, 12.5, 13.0, 13.5, 14.0,
])

_DACL = np.array([
    0.0, 0.5, 0.9, 1.4, 1.9, 2.5, 3.3, 4.0, 4.6, 5.6,
    6.4, 7.3, 8.2, 9.2, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0,
    0.0, 0.2, 0.4, 0.7, 1.2, 1.7, 2.4, 3.0, 3.7, 4.6,
    5.2, 6.0, 6.9, 7.8, 8.6, 9.5, 10.4, 11.4, 12.3, 13.5,
    0.0, 0.0, 0.1, 0.2, 0.5, 0.7, 1.0, 1.3, 1.6, 2.0,
    2.5, 3.0, 3.6, 4.3, 4.9, 5.5, 6.2, 7.0, 7.6, 8.5,
])

_ACLE = np.array([
    0.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 30.0,
])

_C1TABO = np.array([
    0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5,
])

_DACLO = np.array([
    10.0, 8.5, 6.9, 5.5, 4.0, 2.6, 1.5, 0.7, 0.1, 0.0,
    8.7, 7.3, 5.3, 4.2, 2.6, 1.4, 0.5, -0.2, -0.5, 0.0,
    7.5, 5.9, 4.2, 2.5, 1.2, 0.0, -0.7, -1.1, -0.8, 0.0,
    5.5, 3.4, 1.6, 0.0, -1.3, -2.1, -2.5, -2.0, -0.8, 0.0,
    3.0, 0.7, -1.4, -3.3, -4.3, -4.3, -3.1, -2.0, -0.8, 0.0,
    0.3, -2.5, -4.7, -5.8, -5.3, -4.3, -3.1, -2.0, -0.8, 0.0,
    -2.2, -5.0, -6.7, -6.3, -5.3, -4.3, -3.1, -2.0, -0.8, 0.0,
    -3.3, -6.6, -7.2, -6.3, -5.3, -4.3, -3.1, -2.0, -0.8, 0.0,
    -4.2, -7.0, -7.2, -6.3, -5.3, -4.3, -3.1, -2.0, -0.8, 0.0,
    -8.5, -7.9, -7.2, -6.3, -5.3, -4.3, -3.1, -2.0, -0.8, 0.0,
])


# WTLIFT: Figure 4.1.3.2-52, the cranked-wing lift-slope ratio.

_BA = np.array([
    1.15, 1.4, 2.0, 2.2, 3.0, 3.6, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0,
])

_CLOVCL = np.array([
    1.25, 1.2, 1.12, 1.1, 1.04, 1.0, 0.985, 0.96, 0.95, 0.94, 0.94, 0.94,
])


def _columns(flat, fastest: int, count: int) -> np.ndarray:
    """Source storage, ``count`` runs of ``fastest`` values, as TLINEX's
    ``(len(x2), len(x1))`` layout."""
    return np.asarray(flat, dtype=float).reshape(count, fastest).T


def _slope(aspect_ratio: float, section_cla: float, beta: float,
           tan_c2: float) -> float:
    """The Helmbold-Diederich lift-curve slope per degree, as WTLIFT forms it.

    ``2*pi*A / (2 + sqrt((2*pi*A/cla)**2 * (1 + tan^2/beta^2) + 4))``, with
    the section slope ``cla`` per degree and the ``2*pi`` carried in degrees
    through the source's ``DEG``.
    """
    numerator = 2.0 * PI * aspect_ratio * DEG
    aspect_over_cla_squared = (aspect_ratio * DEG * 2.0 * PI / section_cla) ** 2
    sweep_compressibility = 1.0 + tan_c2 ** 2 / beta ** 2
    denominator = 2.0 + np.sqrt(aspect_over_cla_squared * sweep_compressibility + 4.0)
    return numerator / denominator


def calculate_clmxbs(c1p1ac: float, mach: float, a160: float,
                     deltay: float, xovc: float, area: float,
                     sref: float) -> Dict[str, float]:
    """Translate CLMXBS: low-aspect-ratio maximum lift coefficient.

    Args:
        c1p1ac: ``BU4``, ``(C1+1)*A*cos(LE sweep)/beta``, the Figure
            4.1.3.4-23 abscissa.
        mach: ``B(1)``.
        a160: ``A(160)``, ``(C2+1)*A*tan(LE sweep)``.
        deltay: ``AIN(17)``, the leading-edge sharpness parameter.
        xovc: ``AIN(18)``, the maximum-thickness chord station, which picks
            Figure 4.1.3.4-23 part (A) at or below 0.35 and part (B) above.
        area: ``A(3)``, the exposed area.
        sref: ``SREF``.

    Returns:
        Dictionary with ``clmax`` (``CLSMAX``, on the SREF basis), the base
        and increment it is built from, and the ``part`` used.
    """
    increment = tlinex(_AMN, _C2A, _columns(_DE, 15, 3), mach, a160,
                       -1, 2, 0, 2)
    table_part = 'A' if xovc <= 0.35 else 'B'
    cbase_block = _CBASE[:114] if table_part == 'A' else _CBASE[114:]
    base = tlinex(_DYAG, _C1ABC, _columns(cbase_block, 19, 6), deltay, c1p1ac,
                  0, 0, -1, 0)
    return {
        'clmax': float((base + increment) * area / sref),
        'clmax_base': float(base),
        'clmax_increment': float(increment),
        'part': table_part,
        'method': 'legacy_clmxbs',
    }


def calculate_wtlift(planform_type: float,
                     geometry: Dict[str, float],
                     section: Dict[str, float],
                     flight: Dict[str, float],
                     sref: float) -> Dict[str, object]:
    """Translate WTLIFT: lift-curve slope, maximum lift and its angle.

    Args:
        planform_type: ``AIN(15)``; one of the ``cdrag`` WTYPE constants.
        geometry: ``A`` block entries ``area`` A(3), ``aspect_ratio`` A(7),
            ``taper_ratio`` A(27), ``sweep_le_deg`` A(34), ``cos_le`` A(37),
            ``tan_le`` A(38), ``tan_c2`` A(50), and ARCLSS's
            ``arclss_classified`` A(124) and ``arclss_ratio`` A(125).  A
            non-straight planform also needs the inboard and outboard panel
            values ``aspect_ratio_inboard`` A(5), ``tan_c2_inboard`` A(74),
            ``aspect_ratio_outboard`` A(168) and ``tan_c2_outboard`` A(98).
        section: ``deltay`` AIN(17), ``xovc`` AIN(18), and the section
            ``cla`` A(131) per degree and ``clmax`` A(132).
        flight: ``mach`` B(1), ``beta`` B(2), and ``alpha_zero_lift`` B(49).
        sref: ``SREF``.

    Returns:
        Dictionary with ``cla`` (AOUT(101), per degree, SREF basis),
        ``clmax`` (B(44)), ``alpha_clmax`` (B(43)), ``low_aspect_ratio``,
        and the A-block quantities WTLIFT writes.  ``computed`` is False
        for a curved planform: the source prints "NO CLALPHA COMPUTATION"
        and returns having set only the panel slopes ``a171``/``a172``.

    Notes:
        On the low-aspect-ratio path ``B(43)`` comes straight from Figure
        4.1.3.4-25 without the zero-lift angle that the high-aspect-ratio
        path adds, since that figure gives the stall angle itself.
    """
    planform_kind = float(planform_type)
    area = float(geometry['area'])
    aspect_ratio = float(geometry['aspect_ratio'])
    taper = float(geometry['taper_ratio'])
    sweep_le = float(geometry['sweep_le_deg'])
    tan_c2 = float(geometry['tan_c2'])
    tan_le = float(geometry['tan_le'])
    cla_section = float(section['cla'])
    beta = float(flight['beta'])
    mach = float(flight['mach'])
    deltay = float(section['deltay'])

    result = {'method': 'legacy_wtlift', 'computed': planform_kind != CURVED}
    if planform_kind != STRAIGHT_TAPERED:
        result['a172'] = float(_slope(
            float(geometry['aspect_ratio_outboard']),
            cla_section, beta,
            float(geometry['tan_c2_outboard']),
        ))
        result['a171'] = float(_slope(
            float(geometry['aspect_ratio_inboard']),
            cla_section, beta,
            float(geometry['tan_c2_inboard']),
        ))
    if planform_kind == CURVED:
        logger.warning("WTLIFT: no CLALPHA computation for a curved planform")
        return result

    low_aspect_ratio = aspect_ratio < float(geometry['arclss_ratio'])
    a159, _ = tbfunx(_TR, _C2, taper, 0, 0)
    a160 = (a159 + 1.0) * tan_le * aspect_ratio
    a145 = tlinex(_DELTAY, _SALE, _columns(_CLL, 13, 7), deltay, sweep_le,
                  -1, 0, -1, 2)
    a144 = tlinex(_DY, _SALE, _columns(_DACLL, 13, 4), deltay, sweep_le,
                  -1, 0, -1, 2)
    cla = _slope(aspect_ratio, cla_section, beta, tan_c2) * area / sref
    if planform_kind == CRANKED:
        ratio, _ = tbfunx(_BA, _CLOVCL, aspect_ratio * beta, 2, 0)
        cla = ratio * cla
        result['cranked_ratio'] = float(ratio)

    result.update({
        'cla': float(cla),
        'low_aspect_ratio': bool(low_aspect_ratio),
        'a144': float(a144),
        'a145': float(a145),
        'a159': float(a159),
        'a160': float(a160),
    })
    if not low_aspect_ratio:
        # DCAR is DCLTB(5,6,4): Mach fastest, then DELTAY, then sweep.
        mach_correction = tlin3x(
            _DYA, _AMACH, _SALE4,
            _DCAR.reshape(4, 6, 5).transpose(2, 1, 0),
            deltay, mach, sweep_le, 0, -1, 0, 2, 2, 2)
        clmax_section = float(section['clmax'])
        clmax = (a145 * clmax_section + mach_correction) * area / sref
        result.update({
            'a146': float(a145 * clmax_section),
            'mach_correction': float(mach_correction),
            'clmax': float(clmax),
            'alpha_clmax': float(clmax / cla + float(flight['alpha_zero_lift'])
                                 + a144),
        })
        return result

    c1p1ac = aspect_ratio * float(geometry['arclss_classified']) / beta
    base = calculate_clmxbs(c1p1ac, mach, a160, deltay, float(section['xovc']),
                            area, sref)
    alpha_base, _ = tbfunx(_C1ABCS, _ACLMX, c1p1ac, 0, 0)
    if a160 <= 4.5:
        low_ar_abscissa = (aspect_ratio * float(geometry['cos_le']) *
                           (1.0 + 4.0 * taper ** 2))
        increment = tlinex(_ACLE, _C1TABO, _columns(_DACLO, 10, 10),
                           low_ar_abscissa, a160, 0, 0, 2, 0)
    else:
        increment = tlinex(_DMN, _C1TAB, _columns(_DACL, 20, 3),
                           mach, a160, -1, 0, 2, 0)
    result.update({
        'c1p1ac': float(c1p1ac),
        'clmax': base['clmax'],
        'clmxbs': base,
        'alpha_clmax_base': float(alpha_base),
        'alpha_clmax_increment': float(increment),
        'alpha_clmax': float(alpha_base + increment),
    })
    return result
