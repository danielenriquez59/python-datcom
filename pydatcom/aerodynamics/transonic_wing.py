"""
Transonic wing support routines for TRSONI: CLMXB1 and TRANWG.

- ``CLMXB1`` is CLMXBS evaluated at Mach 0.6, the lower anchor of the
  transonic maximum-lift fairing.  Its tables are the same digitisations as
  CLMXBS's (the source labels them Figure 4.1.3.4-16/17 rather than 23/24),
  which a test pins.
- ``TRSONI`` itself: the transonic lift-curve slope, faired by TRANF from
  the subsonic slope below the force-break Mach (Figures 4.1.3.2-53A/B,
  -54A-C) up to TRANWG's value at Mach 1.4; for a low aspect ratio, CLMAX
  and its angle carried from Mach 0.6 by Figures 4.1.3.4-25 and -26; the
  zero-lift drag, friction at Mach 0.6 plus the Figure 4.1.5.1-29 wave drag
  faired over fifteen Mach numbers; and the body's transonic slopes and
  drag (Figures 4.2.3.1-24 and -26).
- ``TRANWG`` is the upper anchor of the lift-curve-slope fairing: the
  supersonic normal-force slope of Section 4.1.3.2 at Mach 1.3, 1.4 and
  1.5, from Figures 4.1.3.2-56A-G and -60A/B, and the slope and its Mach
  derivative at 1.4 from a quadratic in beta through the three.  Its
  figures are the same as VTLIFT's, with the same end modes, and are
  shared with ``vertical_lift_figures``.

Reference: datcom-legacy/datcom_2000/clmxb1.f, tranwg.f, trsoni.f
"""

import math
import numpy as np
from typing import Dict, Mapping, Optional, Sequence
import logging

from pydatcom.aerodynamics.cdrag import STRAIGHT_TAPERED
from pydatcom.aerodynamics.vertical_lift_figures import (
    fig4132_56a, fig4132_56g, fig4132_60a, fig4132_60b,
)
from pydatcom.aerodynamics.wtlift import calculate_clmxbs
from pydatcom.utils.constants import PI, RAD
from pydatcom.utils.legacy_interp import interx
from pydatcom.utils.legacy_numeric import tbfunx, tranf
from pydatcom.utils.legacy_tables import tlinex
from pydatcom.utils.table_lookup import fig26

logger = logging.getLogger(__name__)

# TRANWG's three anchor Mach numbers.
_ANCHOR_MACH = (1.3, 1.4, 1.5)

# The source replaces a zero leading-edge tangent with this, in COMMON.
_ZERO_TANGENT = 0.00001


def calculate_clmxb1(bu4: float, a160: float, deltay: float, xovc: float,
                     area: float, sref: float) -> Dict[str, float]:
    """Translate CLMXB1: CLMXBS at Mach 0.6."""
    return dict(calculate_clmxbs(bu4, 0.60, a160, deltay, xovc, area, sref),
                method='legacy_clmxb1')


def calculate_tranwg(a: Mapping[int, float], deltay: float,
                     sref: float) -> Dict[str, object]:
    """Translate TRANWG: the Mach 1.4 normal-force slope and its derivative.

    Args:
        a: The wing's ``A`` block, one-based: ``3`` exposed area, ``7``
            exposed aspect ratio, ``25`` break-to-root chord ratio, ``27``
            exposed taper, ``58`` inboard LE sweep (degrees), ``61`` its
            cosine, ``62`` and ``86`` the inboard and outboard LE tangents.
        deltay: ``WINGIN(17)``, the leading-edge sharpness parameter.
        sref: ``SREF``.

    Returns:
        Dictionary with ``cna`` (the slope at Mach 1.4, per degree, SREF
        basis), ``dcna`` (its Mach derivative there), the three anchor
        slopes, and the ``a62`` and ``a86`` the routine leaves in COMMON
        (a zero tangent is replaced by ``1e-5``).

    Notes:
        The rectangular-wing branch is taken only for a taper of exactly
        one *and* a leading-edge sweep of exactly zero; any other wing,
        however nearly rectangular, goes through Figure 4.1.3.2-56A.
    """
    tan_out = float(a[86]) or _ZERO_TANGENT
    tan_le = float(a[62]) or _ZERO_TANGENT
    aspect_ratio = float(a[7])
    cos_le = float(a[61])
    rectangular = float(a[25]) == 1.0 and float(a[58]) == 0.0
    slopes, betas, branches = [], [], []
    for mach in _ANCHOR_MACH:
        beta = math.sqrt(mach**2 - 1.0)
        betas.append(beta)
        bovert = beta / tan_le
        supersonic_le = bovert > 1.0
        if not supersonic_le:
            ratio = fig4132_60a(bovert, deltay / cos_le)
        else:
            deltdt = math.atan(deltay / (5.85 * cos_le)) * RAD
            ratio = fig4132_60b(1.0 / bovert, deltdt)
        if not rectangular:
            bcna = fig4132_56a(bovert, aspect_ratio * tan_le, float(a[27]))
            theory = bcna / beta if supersonic_le else bcna / tan_le
            branch = '56a'
        elif aspect_ratio * beta <= 1.0:
            theory = fig4132_56g(aspect_ratio * beta) * aspect_ratio
            branch = '56g'
        else:
            theory = (4.0 - 2.0 * (1.0 / (aspect_ratio * beta))) / beta
            branch = 'rectangular'
        slopes.append(theory * ratio * float(a[3]) / (sref * RAD))
        branches.append(branch)

    cn, bb = slopes, betas
    f1 = (cn[1] * bb[1]**2 - cn[0] * bb[0]**2) / (bb[1] - bb[0])
    f2 = (cn[2] * bb[2]**2 - cn[0] * bb[0]**2) / (bb[2] - bb[0])
    aa = (f2 - f1) / (bb[2] - bb[1])
    b = f1 - aa * (bb[1] + bb[0])
    c = cn[1] * bb[1]**2 - aa * bb[1]**2 - b * bb[1]
    dcna = (-b * _ANCHOR_MACH[1] / bb[1]**3 -
            2.0 * c * _ANCHOR_MACH[1] / bb[1]**4)
    return {'cna': float(cn[1]), 'dcna': float(dcna),
            'anchor_slopes': [float(v) for v in cn], 'branches': branches,
            'a62': tan_le, 'a86': tan_out, 'method': 'legacy_tranwg'}


# TRSONI's tables, extracted from the source DATA by tools/fortran_data.py.
_X = np.array([
    0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9,
    1.0,
])
_Y = np.array([
    0.0, 0.225, 0.47, 0.496, 0.43, 0.32, 0.21, 0.125, 0.075, 0.0475,
    0.0,
])
_TR = np.array([
    0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9,
    1.0,
])
_DR = np.array([
    0.0, 0.21, 0.5, 0.9, 1.08, 1.05, 1.0, 0.94, 0.9, 0.86,
    0.85,
])
_X27M = np.array([
    0.0, 1.0, 2.0, 3.0,
])
_X27I = np.array([
    1.5778, 1.67221, 1.98509, 2.28874,
])
_T43A = np.array([
    0.0, 2.0, 4.0, 6.0, 8.0, 10.0, 12.0, 14.0, 0.0, 1.0,
    2.0, 3.0, 4.0, 6.0, 8.0,
])
_D43A = np.array([
    1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
    1.0, 1.0, 1.0, 1.0, 0.98, 0.9, 1.0, 1.0, 1.0, 1.0,
    1.0, 0.983, 0.89, 0.855, 1.0, 1.0, 1.0, 1.0, 0.939, 0.888,
    0.857, 0.835, 1.0, 1.0, 0.958, 0.904, 0.869, 0.842, 0.82, 0.802,
    1.0, 0.952, 0.908, 0.872, 0.842, 0.816, 0.8, 0.787, 1.0, 0.952,
    0.908, 0.872, 0.842, 0.816, 0.8, 0.787,
])
_T43B = np.array([
    0.0, 0.17453, 0.34907, 0.5236, 0.69813, 1.0472, 1.22173, 1.5708, 0.799, 0.85,
    0.9, 0.95, 1.0,
])
_D43B = np.array([
    0.8, 0.807, 0.82, 0.841, 0.867, 0.918, 0.945, 1.0, 0.85, 0.855,
    0.865, 0.88, 0.9, 0.938, 0.958, 1.0, 0.9, 0.903, 0.91, 0.92,
    0.932, 0.958, 0.974, 1.0, 0.95, 0.952, 0.955, 0.96, 0.967, 0.98,
    0.987, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
])
_T44A = np.array([
    0.0, 2.0, 4.0, 6.0, 8.0, 12.0, 14.0, 1.0, 2.0, 3.0,
    4.0, 6.0, 8.0,
])
_D44A = np.array([
    1.09, 1.088, 1.056, 1.008, 0.962, 0.872, 0.828, 0.98, 1.088, 1.1,
    1.05, 0.992, 0.882, 0.83, 0.85, 1.0, 1.075, 1.06, 1.005, 0.89,
    0.835, 0.74, 0.9, 1.015, 1.07, 1.025, 0.895, 0.837, 1.15, 1.15,
    1.15, 1.15, 1.07, 0.915, 0.915, 1.12, 1.12, 1.12, 1.12, 1.12,
    0.94, 0.94,
])
_T44B = np.array([
    0.0, 4.0, 6.0, 8.0, 10.0, 12.0, 14.0, 1.0, 2.0, 3.0,
    4.0, 6.0, 8.0,
])
_D44B = np.array([
    -0.14, -0.14, -0.14, -0.14, -0.14, 0.0, 0.44, -0.07, -0.07, -0.07,
    -0.07, 0.0, 0.31, 0.64, 0.05, 0.05, 0.05, 0.06, 0.23, 0.49,
    0.75, 0.06, 0.06, 0.06, 0.15, 0.35, 0.58, 0.79, 0.08, 0.08,
    0.14, 0.29, 0.48, 0.66, 0.84, 0.09, 0.15, 0.24, 0.39, 0.555,
    0.72, 0.89,
])
_T44C = np.array([
    0.0, 2.0, 4.0, 6.0, 8.0, 16.0,
])
_D44C = np.array([
    -0.04, 0.01, 0.075, 0.13, 0.15, 0.15,
])
_T418A = np.array([
    0.0, 0.4, 0.8, 1.2, 1.6, 2.0, 2.4, 2.8, 3.2, 3.6,
])
_D418A = np.array([
    35.0, 35.0, 35.0, 32.0, 28.0, 25.0, 23.2, 22.0, 21.5, 21.0,
])
_T18B1 = np.array([
    0.2, 0.4, 0.6,
])
_D18B1 = np.array([
    4.5, 5.0, 5.5, 6.0, 6.5, 7.0, 7.5, 8.0, 8.5, 9.0,
    9.5, 10.0, 10.5, 11.0, 11.5, 12.0, 12.5, 13.0, 13.5, 14.0,
])
_C18B1 = np.array([
    0.0, 0.5, 0.9, 1.4, 1.9, 2.5, 3.3, 4.0, 4.6, 5.6,
    6.4, 7.3, 8.2, 9.2, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0,
    0.0, 0.2, 0.4, 0.7, 1.2, 1.7, 2.4, 3.0, 3.7, 4.6,
    5.2, 6.0, 6.9, 7.8, 8.6, 9.5, 10.4, 11.4, 12.3, 13.5,
    0.0, 0.0, 0.1, 0.2, 0.5, 0.7, 1.0, 1.3, 1.6, 2.0,
    2.5, 3.0, 3.6, 4.3, 4.9, 5.5, 6.2, 7.0, 7.6, 8.5,
])
_T18B2 = np.array([
    0.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 30.0,
])
_D18B2 = np.array([
    0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5,
])
_C18B2 = np.array([
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
_T419A = np.array([
    0.8, 1.0, 1.2, 1.6, 2.0, 2.4, 3.0, 0.0, 0.0, 0.0,
    0.0, 0.6, 0.7, 0.8, 0.9, 0.95, 1.0, 1.05, 1.1, 1.15,
    1.2, 1.4,
])
_D419A = np.array([
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.05, -0.05, -0.05,
    -0.03, -0.015, -0.01, 0.0, -0.115, -0.1, -0.083, -0.047, -0.018, 0.0,
    0.025, -0.145, -0.11, -0.08, -0.02, 0.022, 0.05, 0.097, -0.12, -0.08,
    -0.03, 0.05, 0.1, 0.14, 0.204, -0.06, 0.0, 0.06, 0.15, 0.21,
    0.25, 0.325, 0.09, 0.13, 0.18, 0.257, 0.31, 0.33, 0.36, 0.113,
    0.165, 0.21, 0.28, 0.32, 0.345, 0.36, 0.108, 0.16, 0.205, 0.284,
    0.32, 0.344, 0.36, 0.096, 0.145, 0.19, 0.273, 0.32, 0.343, 0.355,
    0.02, 0.05, 0.1, 0.2, 0.3, 0.33, 0.34,
])
_T419B = np.array([
    0.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0,
    13.0, 14.0, 16.0,
])
_D419B = np.array([
    1.0, 1.0, 0.93, 0.8, 0.55, 0.4, 0.3, 0.24, 0.2, 0.16,
    0.14, 0.12, 0.12,
])
_T419C = np.array([
    0.0, 0.4, 0.8, 1.0, 1.2, 1.6, 2.0, 2.4, 3.0, 5.0,
    0.0, 0.0, 0.0, 0.6, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95,
    1.0, 1.05, 1.1, 1.15, 1.2, 1.4,
])
_D419C = np.array([
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0,
    -1.25, -1.22, -1.2, -1.2, -1.2, -1.1, -1.1, -1.0, -1.0, -0.85,
    -2.05, -1.75, -1.5, -1.4, -1.3, -1.2, -1.1, -1.0, -1.0, -0.95,
    -2.2, -1.9, -1.6, -1.5, -1.4, -1.2, -1.1, -1.0, -0.9, -0.63,
    -1.85, -1.6, -1.4, -1.3, -1.2, -1.0, -0.8, -0.7, -0.6, -0.4,
    -1.8, -1.38, -1.0, -0.9, -0.7, -0.5, -0.3, -0.1, 0.0, 0.1,
    -0.6, -0.2, 0.2, 0.3, 0.5, 0.8, 1.0, 1.2, 1.4, 2.0,
    1.9, 2.4, 3.0, 3.2, 3.5, 3.9, 4.3, 4.6, 5.0, 6.2,
    2.4, 4.05, 5.5, 6.2, 7.0, 8.0, 8.9, 10.0, 10.9, 13.9,
    1.6, 4.04, 6.9, 8.0, 9.1, 10.6, 11.9, 13.0, 14.3, 18.2,
    -1.2, 3.9, 7.6, 8.9, 10.2, 12.1, 13.5, 14.7, 15.9, 18.8,
    -11.0, -0.6, 8.1, 10.5, 12.0, 14.4, 16.2, 17.5, 18.8, 21.3,
])
_T429L = np.array([
    0.0, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.475, 0.0, 0.5,
    1.0, 1.5, 2.0, 3.0, 4.0,
])
_D429L = np.array([
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.9, 0.71,
    0.55, 0.39, 0.21, 0.1, 0.0, 0.0, 2.34, 2.03, 1.7, 1.05,
    0.51, 0.21, 0.01, 0.0, 3.01, 2.92, 2.68, 1.99, 1.0, 0.32,
    0.01, 0.0, 3.3, 3.27, 3.12, 2.5, 1.26, 0.42, 0.03, 0.0,
    3.49, 3.43, 3.38, 2.8, 1.52, 0.51, 0.06, 0.0, 3.61, 3.58,
    3.52, 3.14, 1.83, 0.52, 0.09, 0.0,
])
_T429R = np.array([
    0.0, 0.4, 0.8, 1.2, 1.4, 0.0, 0.0, 0.0, 0.5, 1.0,
    1.5, 2.0, 3.0, 4.0,
])
_D429R = np.array([
    0.0, 0.0, 0.0, 0.0, 0.0, 0.9, 1.08, 1.18, 1.28, 1.31,
    2.34, 2.53, 2.68, 2.75, 2.78, 3.01, 3.08, 3.1, 3.12, 3.12,
    3.3, 3.32, 3.32, 3.3, 3.28, 3.49, 3.49, 3.48, 3.42, 3.39,
    3.61, 3.61, 3.6, 3.56, 3.51,
])
_T424 = np.array([
    0.75, 0.85, 0.9, 0.95, 1.0, 1.015, 1.03, 1.05, 1.075, 1.125,
    1.15, 1.2, 1.3, 0.011, 0.067, 0.122,
])
_D424 = np.array([
    0.011, 0.014, 0.018, 0.031, 0.069, 0.04, 0.02, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.067, 0.071, 0.078, 0.092, 0.14, 0.15, 0.12,
    0.086, 0.067, 0.05, 0.047, 0.045, 0.05, 0.122, 0.13, 0.139, 0.155,
    0.204, 0.22, 0.236, 0.21, 0.168, 0.11, 0.101, 0.1, 0.1,
])
_T426 = np.array([
    6.0, 7.0, 8.0, 10.0, 12.0, 14.0, 16.0, 18.0, 20.0, 24.0,
    26.0, 1.0, 1.025, 1.05, 1.1, 1.2,
])
_D426 = np.array([
    0.118, 0.084, 0.065, 0.041, 0.029, 0.02, 0.016, 0.012, 0.01, 0.007,
    0.006, 0.147, 0.11, 0.088, 0.061, 0.043, 0.033, 0.026, 0.022, 0.019,
    0.015, 0.013, 0.166, 0.12, 0.097, 0.068, 0.05, 0.039, 0.031, 0.026,
    0.022, 0.015, 0.013, 0.186, 0.138, 0.11, 0.076, 0.057, 0.044, 0.034,
    0.028, 0.022, 0.015, 0.013, 0.2, 0.152, 0.12, 0.08, 0.059, 0.046,
    0.036, 0.028, 0.022, 0.015, 0.013,
])


def calculate_trsoni(mach: float, alpha_deg: Sequence[float],
                     wing: Mapping, a: Mapping[int, float], sref: float,
                     roughness: float, body: Optional[Mapping] = None,
                     stale_wave_drag: Optional[Sequence[float]] = None,
                     nf: int = 0) -> Dict[str, object]:
    """Translate TRSONI: transonic wing and body lift, maximum lift and drag.

    Args:
        mach: ``FLC(NZ+2)``.
        alpha_deg: ``FLC(23)`` onward.
        wing: ``type`` (``WINGIN(15)``), ``tovc`` (16), ``deltay`` (17),
            ``xovc`` (18), ``clamo`` (69), and ``cd0`` (``WING(1)``), which
            stands in the wing-body sum when the planform is not straight.
        a: The wing's ``A`` block, one-based: 3, 7, 16, 25, 27, 58, 61, 62,
            67, 71 (half-chord sweep, radians), 74, 86, 128, 129.
        sref: ``SREF``.
        roughness: ``RUFF``, the third word of ``/OPTION/``.
        body: With a body (``BO``): ``cla_06`` and ``cma_06`` (``BODY(101)``,
            ``BODY(121)``), ``cf_06`` ``BD(92)``, ``cd0_06`` ``BD(61)``,
            ``cd_base_06`` ``BD(60)``, ``wetted_area`` ``BD(93)``,
            ``base_area`` ``BD(57)``, and the Mach 1.4 body's ``length``
            ``SBD(2)``, ``base_diameter`` ``SBD(6)``, ``cla_14``
            ``SBD(18)``, ``cma_14`` ``SBD(110)``, ``max_diameter``
            ``SBD(120)`` and ``cd0_14`` ``SBD(124)``.
        stale_wave_drag: ``TRA(43)`` to ``TRA(57)`` as the previous pass left
            them; zero by default.  See Notes.
        nf: ``NF``; a negative value returns after the force-break Mach.

    Returns:
        Dictionary with ``cla`` (``WING(101)``), ``cd0`` (``WING(1)``),
        ``tra`` (the ``TRA`` words defined, one-based), and with a body
        ``body_cla``, ``body_cma`` (the ``BODY(101)``/``(121)`` it
        overwrites), ``body_cd`` (``BODY(1)`` onward), ``cd0_wing_body``
        (``BW(1)``) and ``base_diameter`` (``SBD(6)`` after the floor).

    Notes:
        The wave-drag loop's supersonic branch (normal Mach above one) jumps
        past the statement that scales the Figure 4.1.5.1-29 value and
        stores it in ``CDW2``, so those points of the TRANF fairing hold
        whatever ``TRA(43..57)`` held: zero on the first pass, the previous
        Mach's values after.  Kept; the stale values are an input.

        The Mach 0.6 stall-angle anchor reads Figure 4.1.3.4-25B at
        ``MT(1)``, the lift fairing's lower Mach (up to 0.75, extrapolated
        quadratically past the 0.6 grid), not at 0.6.  Kept.
    """
    return _transonic_surface(mach, alpha_deg, wing, a, sref, roughness, body,
                              stale_wave_drag, nf, False)


def _transonic_surface(mach, alpha_deg, wing, a, sref, roughness, body,
                       stale_wave_drag, nf, stores_supersonic_points):
    """TRSONI and TRSONJ; they differ only in whether the wave-drag loop
    stores its supersonic points."""
    tra: Dict[int, float] = {}
    alpha = np.asarray(alpha_deg, dtype=float)
    result: Dict[str, object] = {'tra': tra, 'method': 'legacy_trsoni'}
    cd0w = float(wing.get('cd0', 0.0))
    tra[4] = float(mach)
    if float(wing['type']) == STRAIGHT_TAPERED:
        g = {int(k): float(v) for k, v in a.items()}
        toc, arstar, srstar = float(wing['tovc']), g[7], g[3]
        tranwg = calculate_tranwg(g, float(wing['deltay']), sref)
        g[62], g[86] = tranwg['a62'], tranwg['a86']
        tanle, tanc2, cosle = g[62], g[74], g[61]
        tra[1], tra[82] = tranwg['cna'], tranwg['dcna']
        k = float(wing['clamo']) * RAD / (2.0 * PI)
        tra[3] = k
        beta6 = 0.80
        mfb0 = interx(2, _T43A, [toc * 100.0, arstar], [8, 7], _D43A,
                      lind=8, lx1l=2, lx2l=2, lx1u=2, lx2u=-1)
        mfb = interx(2, _T43B, [g[71], mfb0], [8, 5], _D43B, lind=8,
                     lx1l=0, lx2l=2, lx1u=2, lx2u=0)
        tra[5], tra[6] = mfb0, mfb
        if nf < 0:
            result['early_return'] = True
            return result
        aoc = interx(2, _T44B, [toc * 100.0, arstar], [7, 6], _D44B,
                     lind=7, lx1l=0, lx2l=2, lx1u=1, lx2u=2)
        cfbct = interx(2, _T44A, [toc * 100.0, arstar], [7, 6], _D44A,
                       lind=7, lx1l=2, lx2l=1, lx1u=2, lx2u=2)
        betafb = 0.0 if mfb > 0.98 else math.sqrt(1.0 - mfb**2)
        arg1 = 2.0 * PI * arstar / RAD
        arg2 = (arstar / k)**2
        clafbt = arg1 / (2.0 + math.sqrt(arg2 * (betafb**2 + tanc2**2) + 4.0))
        arg2_6 = (0.8 * arstar / k)**2
        arg3_6 = 1.0 + (tanc2 / 0.8)**2 if arg2_6 > 0.0 else 0.0
        claw6 = arg1 * srstar / ((2.0 + math.sqrt(arg2_6 * arg3_6 + 4.0)) * sref)
        xm = 0.75
        if xm > mfb - 0.1:
            xm = mfb - 0.1
        if xm > mach:
            xm = mach
        bb = math.sqrt(1.0 - xm**2)
        arg3 = bb**2 + tanc2**2
        arg4 = 2.0 + math.sqrt(arg2 * arg3 + 4.0)
        claw7 = arg1 * srstar / (arg4 * sref)
        dcla7 = (xm * arg1 * arg2 * srstar / sref /
                 (arg4**2 * math.sqrt(arg2 * arg3 + 4.0)))
        clafb = clafbt * cfbct * srstar / sref
        claa = (1.0 - aoc) * clafb
        boc = interx(1, _T44C, [toc * 100.0], [6], _D44C, lind=6,
                     lx1l=2, lx1u=1)
        clab = (1.0 - boc) * clafb
        mt = [xm, mfb, mfb + .07, mfb + .14, 1.40]
        clamt = [claw7, clafb, claa, clab, tranwg['cna']]
        cla = tranf(mt, clamt, dcla7, tranwg['dcna'], mach)
        tra.update({7: aoc, 8: cfbct, 9: betafb, 10: clafbt, 12: clafb,
                    13: claa, 14: boc, 15: clab, 70: claw6})
        tra.update({16 + i: v for i, v in enumerate(mt)})
        tra.update({21 + i: v for i, v in enumerate(clamt)})
        result['cla'] = float(cla)

        # Aspect-ratio classification and, for a low one, CLMAX.
        c1, _ = tbfunx(_X, _Y, g[27], 0, 0)
        aratio = g[128] / ((c1 + 1.0) * cosle)
        c2 = interx(1, _TR, [g[27]], [11], _DR, lind=11)
        tra[27], tra[28] = c1, aratio
        if arstar < aratio:
            bu4 = (c1 + 1.0) * arstar * cosle / beta6
            arg = (c2 + 1.0) * arstar * tanle
            clmax6 = calculate_clmxb1(bu4, arg, float(wing['deltay']),
                                      float(wing['xovc']), srstar,
                                      sref)['clmax']
            aclbas, _ = tbfunx(_T418A, _D418A, bu4, 0, 0)
            if arg <= 4.5:
                dacma6 = tlinex(_T18B2, _D18B2, _C18B2.reshape(10, 10).T,
                                arstar * cosle * (1.0 + 4.0 * g[27]**2),
                                arg, 0, 2, 2, 0)
            else:
                dacma6 = tlinex(_T18B1, _D18B1, _C18B1.reshape(3, 20).T,
                                mt[0], arg, 0, 0, 2, 1)
            alclm6 = aclbas + dacma6
            c3 = interx(1, _T419B, [arg], [13], _D419B, lind=13)
            var1 = bu4 * beta6
            dalcm = interx(2, _T419C, [var1, mach], [10, 13], _D419C,
                           lind=13, lx1l=2, lx1u=2, lx2u=2)
            dclmax = interx(2, _T419A, [c3 * var1, mach], [7, 11], _D419A,
                            lind=11, lx1l=2, lx1u=2, lx2u=2)
            tra.update({29: bu4, 30: clmax6, 31: aclbas, 32: dacma6, 33: c3,
                        34: dalcm, 35: dclmax, 36: alclm6,
                        37: alclm6 + dalcm,
                        38: clmax6 + dclmax * (srstar / sref)})
            result.update({'low_aspect_ratio': True, 'a160': float(arg),
                           'clmax': tra[38], 'alpha_clmax': tra[37]})
        else:
            result['low_aspect_ratio'] = False

        # Zero-lift drag: friction at Mach 0.6 and faired wave drag.
        cept, _ = tbfunx(_X27M, _X27I, 0.6, 0, 0)
        rlcoff = (12.0 * g[16] / roughness)**1.0482 * 10.0**cept
        rnn = min(g[16] * g[129], rlcoff)
        cf = fig26(rnn, 0.60)
        rl = 2.0 if float(wing['xovc']) < 0.30 else 1.2
        cdf = cf * (1.0 + rl * toc) * 2.0 * srstar / sref
        root = math.sqrt(g[67])
        xmtd = np.array([0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.90, 0.925,
                         0.950, 0.975, 1.000, 1.050, 1.100, 1.200,
                         1.400]) / root
        arg2 = toc**.3333
        var2 = arstar * arg2
        arg3 = g[67]**2.5
        cdw2 = (np.zeros(15) if stale_wave_drag is None
                else np.array(stale_wave_drag, dtype=float))
        stale_points = []
        for i, xmi in enumerate(xmtd):
            var1 = math.sqrt(abs((xmi * root)**2 - 1.0)) / arg2
            if xmi * root <= 1.0:
                cdw1 = interx(2, _T429L, [var1, var2], [8, 7], _D429L,
                              lind=8, lx1u=-1, lx2u=2)
                cdw2[i] = cdw1 * toc**1.666 * arg3 * srstar / sref
            else:
                cdw1 = interx(2, _T429R, [var1, var2], [5, 7], _D429R,
                              lind=7, lx2u=2)
                if stores_supersonic_points:
                    cdw2[i] = cdw1 * toc**1.666 * arg3 * srstar / sref
                else:
                    # TRSONI computes CDW1 here and never stores it.
                    stale_points.append(i)
        cdw = tranf(xmtd, cdw2, 0.0, 0.0, mach)
        cd0w = cdw + cdf
        tra.update({39: rlcoff, 40: rnn, 41: rl, 42: cf, 67: cdw, 68: cdf})
        tra.update({43 + i: float(v) for i, v in enumerate(cdw2)})
        result.update({'cd0': float(cd0w), 'wave_drag': float(cdw),
                       'friction_drag': float(cdf),
                       'stale_wave_points': stale_points,
                       'a62': g[62], 'a86': g[86]})

    if body is None:
        return result
    b = {k: float(v) for k, v in body.items()}
    arg1 = (mach - .60) / .80
    result['body_cla'] = b['cla_06'] + (b['cla_14'] - b['cla_06']) * arg1
    result['body_cma'] = b['cma_06'] + (b['cma_14'] - b['cma_06']) * arg1
    cdfb = b['cf_06'] * b['wetted_area'] / sref
    cdpb = b['cd0_06'] - b['cd_base_06'] - cdfb
    if 1.0 <= mach <= 1.2:
        cdpb = cdpb * (1.0 - (mach - 1.0) / .2)
    if mach > 1.2:
        cdpb = 0.0
    dmax = b['max_diameter']
    db = b['base_diameter']
    if db < 0.3 * dmax:
        db = 0.3 * dmax
    arg = (db / dmax)**2
    cdbfig = interx(2, _T424, [mach, b['cd_base_06'] / arg * 4.0 * sref /
                               (PI * dmax**2)], [13, 3], _D424, lind=13,
                    lx2l=2, lx2u=2)
    cdbb = cdbfig * arg * PI * dmax**2 / (4.0 * sref)
    if mach < 1.0:
        cdwb = 0.0
    else:
        cdwb = interx(2, _T426, [b['length'] / dmax, mach], [11, 5], _D426,
                      lind=11, lx1u=2, lx2u=2) * PI * dmax**2 / (4.0 * sref)
    cd0b = cdfb + cdpb + cdbb + cdwb
    if mach > 1.2:
        cd0b = cd0b + (b['cd0_14'] - cd0b) * (mach - 1.2) / 0.2
    tra.update({76: cdbb, 77: cdwb, 78: cd0b, 79: cdfb, 80: cdpb,
                81: cdbfig, 73: cd0b + cd0w})
    result.update({
        'body_cd': cd0b + (alpha / RAD)**2 * b['base_area'] / sref,
        'cd0_wing_body': float(cd0b + cd0w), 'base_diameter': float(db),
    })
    return result


def calculate_trnht(a: Mapping[int, float], deltay: float,
                    sref: float) -> Dict[str, float]:
    """Translate TRNHT: TRANWG on the horizontal tail's blocks.

    The two source routines differ only in their COMMON blocks
    (``/HTDATA/``, ``/HTI/`` for ``/WINGD/``, ``/WINGI/``) and in layout.
    """
    return calculate_tranwg(a, deltay, sref)


def calculate_trsonj(mach: float, alpha_deg: Sequence[float],
                     tail: Mapping, a: Mapping[int, float], sref: float,
                     roughness: float, body: Optional[Mapping] = None,
                     stale_wave_drag: Optional[Sequence[float]] = None,
                     nf: int = 0) -> Dict[str, object]:
    """Translate TRSONJ: TRSONI on the horizontal tail's blocks.

    TRSONI with TRNHT for TRANWG, the tail's ``A``, ``HTIN`` and ``HT``
    blocks, and its own ``TRA`` (``/SBETA/`` word 244 on).  The body words
    are the same ``BODY``/``BD`` words TRSONI writes, and the tail-body drag
    goes to ``BH(1)``.  Arguments and returns are those of
    :func:`calculate_trsoni` with the tail in place of the wing.

    One statement differs.  TRSONJ's label 1090 is on the store
    ``CDW2(I)=CDW1(I)*...``, so the supersonic branch's ``GO TO 1090``
    stores its Figure 4.1.5.1-29 value; TRSONI's 1090 is a ``CONTINUE``
    after the store, which that branch skips.  The tail's wave drag is
    therefore faired through all fifteen points, and ``stale_wave_drag``
    has no effect.
    """
    return _transonic_surface(mach, alpha_deg, tail, a, sref, roughness,
                              body, stale_wave_drag, nf, True)
