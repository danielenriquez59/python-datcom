"""
HINGE: control-surface hinge-moment derivatives.

Two derivatives come out: ``CHA`` with angle of attack and ``CHD`` with
flap deflection, one value of ``CHD`` per deflection in the schedule.

Each is built the same way.  A two-dimensional section value is read from
Section 6.1.3 as a theoretical value times a correction for the real
section, adjusted for the trailing-edge angle, and divided by the
Prandtl-Glauert factor.  A nose balance, if present, multiplies it by a
further factor whose curve depends on the nose shape.  The section value is
then swept to three dimensions and an induced-camber increment from Section
6.1.6 is added, weighted across the flap span by a ``K`` factor read at both
ends.

The flap chord is first converted to a chord *normal to the quarter-chord
line*, which is what the Section 6.1.6 figures want, through the hinge-line
and leading-edge sweeps.

The routine serves both a wing-mounted and a horizontal-tail-mounted
device; the only difference is which geometry block it reads, so one
translation covers both.

Reference: datcom-legacy/datcom_2000/hinge.f
"""

import numpy as np
from typing import Dict, Optional, Sequence
import logging

from pydatcom.utils.constants import RAD, UNUSED
from pydatcom.utils.legacy_tables import tlinex
from pydatcom.utils.legacy_interp import interx

logger = logging.getLogger(__name__)

# The source's NTYPE, reached through a computed GO TO.  An index outside
# 1 to 3 makes a FORTRAN computed GO TO fall through to the next statement,
# which in both of this routine's uses is the sharp-nose branch.
ROUND_NOSE = 1
ELLIPTIC_NOSE = 2
SHARP_NOSE = 3

# Figure 6.1.1.1-39A: thickness ratio  (X1125A)
_F6111_39A_TC = np.array([
    0.0, 0.02, 0.04, 0.06, 0.08, 0.1, 0.12, 0.15,
])

# flap-chord ratio  (X2125A)
_F6111_39A_CFOCA = np.array([
    0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.4, 0.5,
])

# (cl_delta)theory  (Y1125A)
_F6111_39A = np.array([
    1.77, 2.5, 3.0, 3.46, 3.82, 4.16, 4.69, 5.14, 1.77, 2.515, 3.03,
    3.5, 3.873, 4.22, 4.78, 5.24, 1.77, 2.53, 3.06, 3.54, 3.926, 4.29,
    4.87, 5.35, 1.77, 2.545, 3.09, 3.58, 3.979, 4.35, 4.95, 5.46, 1.77,
    2.56, 3.12, 3.62, 4.032, 4.4, 5.04, 5.56, 1.77, 2.575, 3.15, 3.66,
    4.085, 4.48, 5.12, 5.69, 1.77, 2.59, 3.18, 3.7, 4.138, 4.55, 5.21,
    5.79, 1.77, 2.6, 3.22, 3.74, 4.19, 4.62, 5.33, 5.96,
]).reshape((8, 8), order='F')

# Figure 6.1.1.1-39B: cl_alpha ratio  (X1125B)
_F6111_39B_CLOCLT = np.array([
    0.7, 0.72, 0.74, 0.76, 0.78, 0.8, 0.82, 0.84, 0.86, 0.88, 0.9,
    0.92, 0.94, 0.96, 0.98, 1.0,
])

# flap-chord ratio  (X2125B)
_F6111_39B_CFOCA = np.array([
    0.05, 0.1, 0.15, 0.2, 0.25, 0.5,
])

# cl_delta over theory  (Y1125B)
_F6111_39B = np.array([
    0.356, 0.382, 0.409, 0.431, 0.452, 0.548, 0.399, 0.426, 0.452,
    0.477, 0.498, 0.583, 0.442, 0.471, 0.499, 0.523, 0.543, 0.619,
    0.485, 0.521, 0.548, 0.569, 0.589, 0.659, 0.53, 0.569, 0.594,
    0.613, 0.63, 0.693, 0.578, 0.614, 0.639, 0.657, 0.671, 0.729,
    0.619, 0.655, 0.678, 0.692, 0.709, 0.761, 0.659, 0.696, 0.713,
    0.733, 0.746, 0.793, 0.7, 0.734, 0.75, 0.765, 0.778, 0.819, 0.742,
    0.771, 0.789, 0.8, 0.81, 0.85, 0.784, 0.809, 0.824, 0.838, 0.843,
    0.875, 0.826, 0.843, 0.86, 0.865, 0.873, 0.9, 0.865, 0.885, 0.895,
    0.9, 0.903, 0.921, 0.91, 0.921, 0.928, 0.931, 0.933, 0.938, 0.951,
    0.962, 0.964, 0.966, 0.967, 0.968, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
]).reshape((6, 16), order='F')

# Figure 4.1.1.2-8A: log10(Reynolds * MAC)  (X1128A)
_F4112_8A_LOG_RF = np.array([
    6.0, 7.0, 8.0,
])

# trailing-edge angle tangent  (X2128A)
_F4112_8A_TANPHP = np.array([
    0.0, 0.02, 0.04, 0.06, 0.08, 0.1, 0.12, 0.14, 0.16, 0.18, 0.2,
])

# cl_alpha over theory  (Y1128A)
_F4112_8A = np.array([
    0.9, 0.878, 0.858, 0.836, 0.815, 0.794, 0.772, 0.75, 0.728, 0.708,
    0.685, 0.95, 0.938, 0.924, 0.907, 0.894, 0.878, 0.86, 0.842, 0.822,
    0.802, 0.78, 0.966, 0.957, 0.947, 0.936, 0.924, 0.91, 0.896, 0.88,
    0.862, 0.842, 0.822,
]).reshape((11, 3), order='F')

# Figure 6.1.3.1-11A: thickness ratio  (X1317A)
_F6131_11A_TC = np.array([
    0.0, 0.04, 0.06, 0.08, 0.1, 0.12, 0.15,
])

# flap-chord ratio  (X2317A)
_F6131_11A_CFOCA = np.array([
    0.0, 0.05, 0.1, 0.175, 0.25, 0.4,
])

# (ch_alpha)theory  (Y1317A)
_F6131_11A = np.array([
    0.0, -0.245, -0.345, -0.465, -0.565, -0.745, 0.0, -0.225, -0.325,
    -0.445, -0.54, -0.72, 0.0, -0.205, -0.305, -0.425, -0.52, -0.71,
    0.0, -0.185, -0.285, -0.405, -0.505, -0.7, 0.0, -0.17, -0.27,
    -0.385, -0.485, -0.685, 0.0, -0.15, -0.25, -0.363, -0.465, -0.67,
    0.0, -0.125, -0.225, -0.336, -0.435, -0.646,
]).reshape((6, 7), order='F')

# Figure 6.1.3.1-11B: cl_alpha ratio  (X1317B)
_F6131_11B_CLACLT = np.array([
    0.7, 0.72, 0.74, 0.76, 0.78, 0.8, 0.82, 0.84, 0.86, 0.88, 0.9,
    0.92, 0.94, 0.96, 0.98, 1.0,
])

# flap-chord ratio  (X2317B)
_F6131_11B_CFOCA = np.array([
    0.1, 0.4,
])

# ch_alpha over theory  (Y1317B)
_F6131_11B = np.array([
    -0.11, 0.13, -0.01, 0.21, 0.08, 0.3, 0.175, 0.38, 0.27, 0.46, 0.35,
    0.54, 0.43, 0.61, 0.51, 0.66, 0.58, 0.71, 0.65, 0.76, 0.71, 0.8,
    0.77, 0.84, 0.82, 0.89, 0.88, 0.93, 0.94, 0.96, 1.0, 1.0,
]).reshape((2, 16), order='F')

# Figure 6.1.3.1-12A: nose-balance ratio  (X1318A)
_F6131_12A_BALANCE = np.array([
    0.0, 0.15, 0.185, 0.3, 0.35, 0.4, 0.5,
])

# sharp nose  (Y318A1)
_F6131_12A_SHARP = np.array([
    1.0, 1.0, 1.0, 0.81, 0.7, 0.57, 0.26,
])

# elliptic nose  (Y318A2)
_F6131_12A_ELLIPTIC = np.array([
    1.0, 0.98, 0.9, 0.63, 0.51, 0.4, 0.16,
])

# round nose  (Y318A3)
_F6131_12A_ROUND = np.array([
    1.0, 0.93, 0.84, 0.54, 0.42, 0.28, 0.03,
])

# Figure 6.1.3.2-12A: thickness ratio  (X1327A)
_F6132_12A_TC = np.array([
    0.0, 0.04, 0.06, 0.08, 0.1, 0.12, 0.15,
])

# flap-chord ratio  (X2327A)
_F6132_12A_CFOCA = np.array([
    0.1, 0.15, 0.2, 0.25, 0.4,
])

# (ch_delta)theory  (Y1327A)
_F6132_12A = np.array([
    -0.883, -0.901, -0.92, -0.944, -1.01, -0.83, -0.855, -0.885,
    -0.913, -0.995, -0.8, -0.83, -0.862, -0.895, -0.984, -0.77, -0.805,
    -0.84, -0.875, -0.972, -0.735, -0.775, -0.814, -0.85, -0.958,
    -0.696, -0.74, -0.783, -0.824, -0.94, -0.639, -0.683, -0.73,
    -0.777, -0.92,
]).reshape((5, 7), order='F')

# Figure 6.1.3.2-12B: cl_alpha ratio  (X1327B)
_F6132_12B_CLACLT = np.array([
    0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0,
])

# flap-chord ratio  (X2327B)
_F6132_12B_CFOCA = np.array([
    0.1, 0.2, 0.25, 0.3, 0.35, 0.4,
])

# ch_delta over theory  (Y1327B)
_F6132_12B = np.array([
    0.646, 0.595, 0.56, 0.52, 0.47, 0.419, 0.705, 0.67, 0.65, 0.62,
    0.585, 0.545, 0.755, 0.735, 0.72, 0.704, 0.685, 0.66, 0.8, 0.788,
    0.779, 0.767, 0.755, 0.739, 0.845, 0.836, 0.83, 0.821, 0.814, 0.8,
    0.884, 0.876, 0.87, 0.868, 0.864, 0.856, 0.925, 0.919, 0.915,
    0.913, 0.91, 0.909, 0.964, 0.961, 0.96, 0.959, 0.958, 0.955, 1.0,
    1.0, 1.0, 1.0, 1.0, 1.0,
]).reshape((6, 9), order='F')

# Figure 6.1.3.2-13A: sharp nose, balance ratio  (X1328A)
_F6132_13A_BALANCE = np.array([
    0.0, 0.185, 0.5,
])

# sharp nose  (Y1328A)
_F6132_13A = np.array([
    1.0, 1.0, 0.5,
])

# Figure 6.1.3.2-13B: elliptic nose, thickness  (X1328B)
_F6132_13B_TC = np.array([
    0.09, 0.15,
])

# balance ratio  (X2328B)
_F6132_13B_BALANCE = np.array([
    0.0, 0.185, 0.3, 0.4, 0.5,
])

# elliptic nose  (Y1328B)
_F6132_13B = np.array([
    1.0, 0.86, 0.66, 0.44, 0.2, 1.0, 0.87, 0.7, 0.54, 0.36,
]).reshape((5, 2), order='F')

# Figure 6.1.3.2-13C: round nose, thickness  (X1328C)
_F6132_13C_TC = np.array([
    0.09, 0.15,
])

# balance ratio  (X2328C)
_F6132_13C_BALANCE = np.array([
    0.0, 0.175, 0.3, 0.4, 0.46,
])

# round nose  (Y1328C)
_F6132_13C = np.array([
    1.0, 0.74, 0.31, -0.1, -0.3, 1.0, 0.78, 0.47, 0.17, 0.0,
]).reshape((5, 2), order='F')

# Figure 6.1.6.1-19A: aspect ratio  (X6115A)
_F6161_19A_AR = np.array([
    2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0,
])

# delta ch_alpha over constant  (Y6115A)
_F6161_19A = np.array([
    0.0182, 0.014, 0.0108, 0.0085, 0.0068, 0.0055, 0.0046, 0.0039,
    0.0035,
])

# Figure 6.1.6.1-19B: span station  (X6115B)
_F6161_19B_ETA = np.array([
    0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.74, 0.8, 0.9, 1.0,
])

# K_alpha  (Y6115B)
_F6161_19B = np.array([
    1.0, 1.12, 1.25, 1.43, 1.65, 1.92, 2.22, 2.62, 2.8, 3.06, 3.63,
    4.26,
])

# Figure 6.1.6.1-19C: balance-chord ratio  (X16116)
_F6161_19C_CBOCF = np.array([
    0.0, 0.2, 0.3, 0.4, 0.5, 0.6,
])

# normal flap-chord ratio  (X26116)
_F6161_19C_CFOCAP = np.array([
    0.0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55,
    0.6,
])

# B2  (Y16116)
_F6161_19C = np.array([
    0.0, 0.49, 0.65, 0.8, 0.92, 1.02, 1.09, 1.16, 1.22, 1.28, 1.33,
    1.38, 1.42, 0.0, 0.44, 0.6, 0.73, 0.85, 0.93, 1.01, 1.08, 1.14,
    1.19, 1.25, 1.29, 1.34, 0.0, 0.39, 0.54, 0.65, 0.75, 0.84, 0.92,
    0.99, 1.05, 1.1, 1.16, 1.21, 1.25, 0.0, 0.32, 0.45, 0.55, 0.63,
    0.71, 0.77, 0.85, 0.92, 0.98, 1.04, 1.09, 1.15, 0.0, 0.22, 0.31,
    0.39, 0.47, 0.54, 0.6, 0.66, 0.72, 0.79, 0.86, 0.93, 1.01, 0.0,
    0.09, 0.16, 0.21, 0.27, 0.33, 0.4, 0.45, 0.52, 0.59, 0.67, 0.78,
    0.91,
]).reshape((13, 6), order='F')

# Figure 6.1.6.2-15A: normal flap-chord ratio  (X1629A)
_F6162_15A_CFOCAP = np.array([
    0.2, 0.4, 0.6,
])

# aspect ratio  (X2629A)
_F6162_15A_AR = np.array([
    2.0, 2.5, 3.0, 3.5, 4.0, 5.0, 6.0, 7.0, 8.0, 10.0,
])

# delta ch_delta over constant  (Y1629A)
_F6162_15A = np.array([
    0.039, 0.03, 0.025, 0.0214, 0.0183, 0.0142, 0.0113, 0.009, 0.0072,
    0.005, 0.035, 0.028, 0.023, 0.0195, 0.0168, 0.013, 0.01, 0.0082,
    0.0065, 0.0046, 0.0305, 0.0246, 0.0205, 0.0175, 0.0151, 0.0118,
    0.0094, 0.0075, 0.0061, 0.0043,
]).reshape((10, 3), order='F')

# Figure 6.1.6.2-15B: span station  (X1629B)
_F6162_15B_ETA = np.array([
    0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.744, 0.8, 0.9, 1.0,
])

# K_delta  (Y1629B)
_F6162_15B = np.array([
    1.0, 1.08, 1.2, 1.34, 1.52, 1.75, 2.05, 2.4, 2.6, 2.91, 3.56, 4.34,
])


def _nose_balance_ratio(chord_balance: float, chord_flap_average: float,
                        thickness_at_hinge: float):
    """The source's BRATIO, or None when no nose balance applies.

    ``BRATIO`` is ``sqrt((cb/cf)^2 - (t/2cf)^2)``.  The source skips the
    whole nose-balance block when the balance chord is absent or the
    radicand is negative, and detects that later by testing ``BRATIO``
    against UNUSED -- which works because ``INITZ2`` fills the ``FHG``
    block with UNUSED before the case runs.
    """
    if (
        chord_balance is None
        or chord_balance == UNUSED
        or chord_balance == 0.0
        or chord_flap_average == 0.0
    ):
        return None
    chord_over_flap = chord_balance / chord_flap_average
    half_thickness_over_flap = thickness_at_hinge / (2.0 * chord_flap_average)
    radicand = chord_over_flap ** 2 - half_thickness_over_flap ** 2
    if radicand < 0.0:
        return None
    return float(np.sqrt(radicand))


def _alpha_nose_factor(nose_type: int, balance_ratio: float) -> float:
    """Figure 6.1.3.1-12A, one curve per nose shape."""
    if nose_type == ROUND_NOSE:
        curve = _F6131_12A_ROUND
    elif nose_type == ELLIPTIC_NOSE:
        curve = _F6131_12A_ELLIPTIC
    else:
        curve = _F6131_12A_SHARP
    return float(interx(1, _F6131_12A_BALANCE, [balance_ratio], [7], curve,
                        lind=7, lx1u=1))


def _delta_nose_factor(nose_type: int, balance_ratio: float,
                       thickness: float) -> float:
    """Figures 6.1.3.2-13A, -13B and -13C, one per nose shape.

    The sharp-nose curve is a function of the balance ratio alone; the
    elliptic and round ones also carry thickness.
    """
    if nose_type == ELLIPTIC_NOSE:
        return float(tlinex(_F6132_13B_TC, _F6132_13B_BALANCE, _F6132_13B,
                            thickness, balance_ratio, 1, 0, 1, 1))
    if nose_type == ROUND_NOSE:
        return float(tlinex(_F6132_13C_TC, _F6132_13C_BALANCE, _F6132_13C,
                            thickness, balance_ratio, 1, 0, 1, 1))
    return float(interx(1, _F6132_13A_BALANCE, [balance_ratio], [3],
                        _F6132_13A, lind=3, lx1u=1))


def calculate_hinge(mach: float,
                    deflections: Sequence[float],
                    surface: Dict[str, float],
                    flap: Dict[str, float],
                    section_lift_slope_ratio: float,
                    reynolds_per_length: float,
                    section_cl_alpha: float,
                    section_dcl: Optional[Sequence[float]] = None,
                    section_dcl_table: Optional[Sequence[float]] = None,
                    alpha_delta_table: Optional[Sequence[float]] = None,
                    nose_type: int = SHARP_NOSE) -> Dict[str, object]:
    """Translate HINGE: hinge-moment derivatives of a trailing-edge flap.

    Args:
        mach: Free-stream Mach number, the source's ``FLC(M+2)``.  Must be
            subsonic; the routine forms ``sqrt(1-M^2)``.
        deflections: The deflection schedule ``DELTA``, the source's
            ``F(1)`` onward.  Its length is the source's ``NDELTA``.
        surface: Carrying-surface geometry, whether wing or horizontal
            tail, as ``{'tovc', 'tovco', 'semispan', 'semispan_exposed',
            'sspnop', 'aspect_ratio', 'sin_c4', 'cos_c4', 'sweep_c4_deg',
            'sweep_le_deg', 'cos_le', 'tan_te', 'mac_exposed'}``.  These
            are ``WINGIN(16)``, ``WINGIN(65)``, ``WINGIN(4)``,
            ``WINGIN(3)``, ``WINGIN(2)``, ``A(120)``, ``A(66)``, ``A(67)``,
            ``A(64)``, ``A(58)``, ``A(61)``, ``A(80)`` and ``A(16)``, or
            their ``HTIN``/``AHT`` counterparts.
        flap: Flap geometry as ``{'span_inboard', 'span_outboard',
            'chord_inboard', 'chord_outboard', 'chord_ratio',
            'chord_balance', 'thickness_at_hinge', 'tan_te_angle',
            'tan_te_angle_effective'}``, the source's ``F(14)``, ``F(15)``,
            ``F(12)``, ``F(13)``, ``FLP(61)``, ``F(59)``, ``F(60)``,
            ``F(11)`` and ``F(61)``.
        section_lift_slope_ratio: ``CLACLT``, the source's ``FLP(33)``.
        reynolds_per_length: The source's ``FLC(M+42)``.
        section_cl_alpha: ``CLASEC``, the section lift-curve slope.
        section_dcl: ``SDCL``, one lift increment per deflection.  When
            supplied the routine uses it directly; otherwise it averages
            the four-per-deflection tables below.
        section_dcl_table: ``DCL``, four values per deflection.
        alpha_delta_table: ``ALDAG``, four values per deflection.
        nose_type: ``NTYPE``.  See the module constants.

    Returns:
        Dictionary with ``cha``, the angle-of-attack derivative, and
        ``chd``, one deflection derivative per deflection, together with
        the section values and figure results behind them.

    Raises:
        ValueError: If the Mach number is not subsonic, if the flap span is
            degenerate, if a deflection is zero, or if neither ``SDCL`` nor
            the two four-per-deflection tables are supplied.

    Notes:
        The source's section-increment tables are averaged four at a time,
        ``SUM/4.`` over ``ALDAG`` and ``SUM2/4.`` over ``DCL``, with a
        running index that walks straight through both arrays.  That
        averaging is reproduced rather than reinterpreted.

        The thickness that enters every section figure is the outboard
        value whenever the flap starts outboard of the planform break --
        the source's ``IF(BIF.GE.(BO2-SSPNOP)) TC=TOVCO`` -- so a flap that
        straddles the break is treated wholly as outboard.
    """
    if not np.isfinite(mach) or mach >= 1.0:
        raise ValueError(
            f"HINGE forms sqrt(1-M^2); Mach {mach} is not subsonic",
        )

    deflections = np.atleast_1d(np.asarray(deflections, dtype=float))
    if deflections.size == 0 or np.any(deflections == 0.0):
        raise ValueError("HINGE divides by each deflection; none may be zero")

    beta = float(np.sqrt(1.0 - mach ** 2))

    span_inboard = float(flap['span_inboard'])
    span_outboard = float(flap['span_outboard'])
    chord_inboard = float(flap['chord_inboard'])
    chord_outboard = float(flap['chord_outboard'])
    span_delta = span_outboard - span_inboard
    if span_delta == 0.0:
        raise ValueError("HINGE needs a nonzero flap span")

    # Hinge-line sweep, and the flap chord normal to the quarter chord.
    chord_delta = chord_inboard - chord_outboard
    sweep_hinge = float(np.arctan(
        (span_delta * float(surface['tan_te']) + chord_delta) / span_delta))
    cos_hinge = float(np.cos(sweep_hinge))
    sweep_c4 = float(surface['sweep_c4_deg'])
    sweep_le = float(surface['sweep_le_deg'])
    cos_c4 = float(surface['cos_c4'])
    sin_c4 = float(surface['sin_c4'])
    chord_ratio = float(flap['chord_ratio'])
    normal_ratio = (cos_hinge / np.cos(sweep_c4 / RAD - sweep_hinge) *
                    np.cos((sweep_c4 - sweep_le) / RAD) /
                    float(surface['cos_le']) * chord_ratio)

    chord_balance = flap.get('chord_balance')
    balance_over_flap = (0.0 if chord_balance in (None, UNUSED) else
                         float(chord_balance) * 2.0 /
                         (chord_inboard + chord_outboard))
    tan_bl = ((chord_delta * (1.0 + balance_over_flap) +
               span_delta * float(surface['tan_te'])) / span_delta)
    balance_over_flap = (balance_over_flap *
                         (cos_c4 + float(surface['tan_te']) * sin_c4) /
                         (cos_c4 + tan_bl * sin_c4))

    # The section lift increments, either supplied or averaged in fours.
    count = deflections.size
    if section_dcl is not None:
        increments = np.atleast_1d(np.asarray(section_dcl, dtype=float))
        alpha_delta = -increments / (deflections * section_cl_alpha)
    elif section_dcl_table is not None and alpha_delta_table is not None:
        dcl = np.asarray(section_dcl_table, dtype=float)
        aldag = np.asarray(alpha_delta_table, dtype=float)
        if dcl.size < 4 * count or aldag.size < 4 * count:
            raise ValueError(
                "HINGE averages four DCL and ALDAG entries per deflection")
        increments = dcl[:4 * count].reshape(count, 4).sum(axis=1) / 4.0
        alpha_delta = aldag[:4 * count].reshape(count, 4).sum(axis=1) / 4.0
    else:
        raise ValueError(
            "HINGE needs either SDCL or both the DCL and ALDAG tables")

    # A flap starting outboard of the break takes the outboard thickness.
    thickness = float(surface['tovc'])
    if span_inboard >= (float(surface['semispan']) -
                        float(surface['sspnop'])):
        thickness = float(surface['tovco'])

    slope_ratio = float(section_lift_slope_ratio)
    tan_te_effective = float(flap['tan_te_angle_effective'])

    # ---- section derivative with angle of attack ------------------------
    cl_alpha_theory = 5.0 * thickness + 6.28
    ch_alpha_theory = float(tlinex(_F6131_11A_TC, _F6131_11A_CFOCA,
                                   _F6131_11A, thickness, chord_ratio,
                                   0, 0, 1, 1))
    ch_alpha_ratio = float(tlinex(_F6131_11B_CLACLT, _F6131_11B_CFOCA,
                                  _F6131_11B, slope_ratio, chord_ratio,
                                  1, 1, 0, 1))
    ch_alpha = ch_alpha_ratio * ch_alpha_theory / RAD
    ch_alpha_corrected = (ch_alpha + 2.0 * cl_alpha_theory *
                          (1.0 - slope_ratio) *
                          (tan_te_effective - thickness) / RAD)

    balance_ratio = _nose_balance_ratio(
        chord_balance, (chord_inboard + chord_outboard) / 2.0,
        float(flap.get('thickness_at_hinge', 0.0) or 0.0))
    alpha_nose = None
    if balance_ratio is not None:
        alpha_nose = _alpha_nose_factor(int(nose_type), balance_ratio)
        cha_mac = ch_alpha_corrected * alpha_nose / beta
    else:
        cha_mac = ch_alpha_corrected / beta

    # ---- section derivative with deflection -----------------------------
    ch_delta_ratio = float(tlinex(_F6132_12B_CLACLT, _F6132_12B_CFOCA,
                                  _F6132_12B, slope_ratio, chord_ratio,
                                  1, 1, 0, 1))
    ch_delta_theory = float(tlinex(_F6132_12A_TC, _F6132_12A_CFOCA,
                                   _F6132_12A, thickness, chord_ratio,
                                   0, 1, 1, 1))
    ch_delta = ch_delta_ratio * ch_delta_theory / RAD

    log_reynolds = float(np.log10(reynolds_per_length *
                                  float(surface['mac_exposed'])))
    cl_alpha_ratio = float(tlinex(_F4112_8A_LOG_RF, _F4112_8A_TANPHP,
                                  _F4112_8A, log_reynolds,
                                  float(flap['tan_te_angle']), 1, 0, 0, 1))
    cl_delta_ratio = float(tlinex(_F6111_39B_CLOCLT, _F6111_39B_CFOCA,
                                  _F6111_39B, cl_alpha_ratio, chord_ratio,
                                  0, 0, 0, 0))
    cl_delta_theory = float(tlinex(_F6111_39A_TC, _F6111_39A_CFOCA,
                                   _F6111_39A, thickness, chord_ratio,
                                   0, 2, 1, 2))
    ch_delta_corrected = (ch_delta + 2.0 * cl_delta_theory *
                          (1.0 - cl_delta_ratio) *
                          (tan_te_effective - thickness) / RAD)

    delta_nose = None
    if balance_ratio is not None:
        delta_nose = _delta_nose_factor(int(nose_type), balance_ratio,
                                        thickness)
        chd_mac = ch_delta_corrected * delta_nose / beta
    else:
        chd_mac = ch_delta_corrected / beta

    # ---- three-dimensional derivatives ----------------------------------
    aspect_ratio = float(surface['aspect_ratio'])
    dcha_over_k = float(interx(1, _F6161_19A_AR, [aspect_ratio], [9],
                               _F6161_19A, lind=9, lx1l=2))
    b2 = float(tlinex(_F6161_19C_CBOCF, _F6161_19C_CFOCAP, _F6161_19C,
                      balance_over_flap, normal_ratio, 0, 0, 1, 1))

    semispan = float(surface['semispan'])
    eta_inboard = span_inboard / semispan
    eta_outboard = span_outboard / semispan
    eta_span = eta_outboard - eta_inboard
    if eta_span == 0.0:
        raise ValueError("HINGE divides by the flap span fraction")

    k_alpha_in = float(interx(
        1, _F6161_19B_ETA, [eta_inboard], [12], _F6161_19B, lind=12, lx1u=2,
    ))
    k_alpha_out = float(interx(
        1, _F6161_19B_ETA, [eta_outboard], [12], _F6161_19B, lind=12, lx1u=2,
    ))
    k_alpha = (
        k_alpha_in * (1.0 - eta_inboard) - k_alpha_out * (1.0 - eta_outboard)
    ) / eta_span
    delta_cha = dcha_over_k * (section_cl_alpha * b2 * k_alpha * cos_c4)
    cha = (aspect_ratio * cos_c4 * cha_mac /
           (aspect_ratio + 2.0 * cos_c4) + delta_cha)

    k_delta_in = float(interx(
        1, _F6162_15B_ETA, [eta_inboard], [12], _F6162_15B, lind=12, lx1u=2,
    ))
    k_delta_out = float(interx(
        1, _F6162_15B_ETA, [eta_outboard], [12], _F6162_15B, lind=12, lx1u=2,
    ))
    k_delta = (
        k_delta_in * (1.0 - eta_inboard) - k_delta_out * (1.0 - eta_outboard)
    ) / eta_span
    dchd_over_k = float(tlinex(_F6162_15A_CFOCAP, _F6162_15A_AR,
                               _F6162_15A, normal_ratio, aspect_ratio,
                               1, 1, 1, 0))

    sweep_factor = b2 * k_delta * cos_c4 * cos_hinge
    rotation = cos_c4 * cos_hinge
    induced = 2.0 * cos_c4 / (aspect_ratio + 2.0 * cos_c4)
    dchd = dchd_over_k * increments * sweep_factor / deflections
    chd = rotation * (chd_mac + alpha_delta * cha_mac * induced) + dchd

    return {
        'cha': float(cha),
        'chd': chd,
        'cha_section': float(cha_mac),
        'chd_section': float(chd_mac),
        'delta_cha': float(delta_cha),
        'dchd': dchd,
        'k_alpha': float(k_alpha),
        'k_delta': float(k_delta),
        'b2': float(b2),
        'dcha_over_k': float(dcha_over_k),
        'dchd_over_k': float(dchd_over_k),
        'balance_ratio': balance_ratio,
        'alpha_nose_factor': alpha_nose,
        'delta_nose_factor': delta_nose,
        'has_nose_balance': balance_ratio is not None,
        'thickness_used': float(thickness),
        'thickness_is_outboard': thickness == float(surface['tovco']),
        'normal_chord_ratio': float(normal_ratio),
        'balance_over_flap': float(balance_over_flap),
        'sweep_hinge_deg': float(sweep_hinge * RAD),
        'beta': beta,
        'section_increments': increments,
        'alpha_delta': alpha_delta,
        'method': 'legacy_hinge',
    }
