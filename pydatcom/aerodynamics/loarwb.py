"""
Low-aspect-ratio wing-body: LOARWB and its overlay M14O16 (DATCOM Section
4.8 and the Section 5.5 lateral method).

For a blended low-aspect-ratio configuration described as a whole
(``LARWB`` namelist, the ``/POWER/`` ``LBIN`` block) rather than as a wing
plus a body:

- Normal force from Figures 4.8.1.2-12 to -14, by the leading-edge radius
  parameter or, when given, the leading-edge angle, with the nonlinear
  term to 20 degrees.
- Axial force: skin friction (Figure 4.1.5.1-27 cutoff and FIG26), base
  pressure (Figures 4.8.2.1-7, 4.8.2.2-11A/B) and the Figure 4.8.2.2-10
  increment.
- Pitching moment about the centre of pressure of Figures 4.8.3.2-6B and
  -7A plus the frontal-area term.
- Side force, yawing and rolling moment per degree of sideslip, Figures
  5.5.1.1-6 to 5.5.3.1-6.

Curves formed on the angle measured from the zero-lift line are
interpolated back to the flight angles.  M14O16 then forms CN, CA and the
two slopes.

The tables were extracted from the source DATA statements by parsing, with
``tools/fortran_data.py``, rather than by hand.

Reference: datcom-legacy/datcom_2000/loarwb.f, m14o16.f
"""

import math
from typing import Dict, Mapping, Sequence

import numpy as np

from pydatcom.utils.constants import PI, RAD, UNUSED
from pydatcom.utils.legacy_numeric import tbfunx
from pydatcom.utils.legacy_tables import tlinex
from pydatcom.utils.table_lookup import fig26

# Figure 4.1.5.1-27 intercept (the cutoff Reynolds number).
_X27M = np.array([
    0.0, 1.0, 2.0, 3.0,
])
_X27I = np.array([
    1.5778, 1.67221, 1.98509, 2.28874,
])
# Figures 4.8.1.2-12A, -13A, -14A: against the leading-edge angle.
_DELTAE = np.array([
    0.0, 10.0, 20.0, 40.0, 60.0, 80.0, 100.0, 120.0,
])
_SCNA0 = np.array([
    1.18, 1.16, 1.11, 1.0, 0.91, 0.82, 0.73, 0.65,
])
# Figures 4.8.1.2-12B, -13B, -14B: against the leading-edge radius parameter.
_RLEOB = np.array([
    0.0, 0.02, 0.04, 0.06, 0.08, 0.1, 0.12,
])
_RCNA0 = np.array([
    1.18, 1.01, 0.89, 0.82, 0.8, 0.78, 0.77,
])
_SCNA20 = np.array([
    1.18, 1.17, 1.15, 1.09, 1.03, 0.98, 0.93, 0.89,
])
_RCNA20 = np.array([
    1.18, 0.97, 0.82, 0.76, 0.73, 0.7, 0.68,
])
_ZS = np.array([
    1.0, 1.0, 1.0, 1.0, 1.0, 0.98, 0.96, 0.92,
])
_ZR = np.array([
    1.0, 0.96, 0.83, 0.6, 0.37, 0.2, 0.2,
])
# Figure 4.8.2.1-7: base pressure against the base shape parameter.
_SBOLHB = np.array([
    0.0, 0.02, 0.04, 0.06, 0.08, 0.1, 0.12,
])
_CPBOP = np.array([
    0.0, -0.08, -0.14, -0.19, -0.225, -0.245, -0.255,
])
_X210 = np.array([
    0.0, 0.1, 0.2, 0.3, 0.4, 0.5,
])
# Figure 4.8.2.2-10: the axial-force increment (TLINEX, 11 by 6).
_X110 = np.array([
    0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7,
    0.8, 0.9, 1.0,
])
_Y10 = np.array([
    0.0, 0.2, 0.31, 0.38, 0.43, 0.48, 0.0, 0.23,
    0.37, 0.46, 0.52, 0.58, 0.0, 0.27, 0.42, 0.52,
    0.6, 0.66, 0.0, 0.31, 0.49, 0.605, 0.69, 0.73,
    0.0, 0.34, 0.56, 0.68, 0.77, 0.82, 0.0, 0.39,
    0.605, 0.75, 0.84, 0.91, 0.0, 0.415, 0.67, 0.82,
    0.92, 1.0, 0.0, 0.45, 0.71, 0.89, 1.0, 1.08,
    0.0, 0.49, 0.775, 0.95, 1.07, 1.14, 0.0, 0.51,
    0.81, 1.0, 1.115, 1.19, 0.0, 0.55, 0.87, 1.06,
    1.16, 1.195,
])
# Figures 4.8.2.2-11A/B: base pressure at 20 degrees.
_B2OHS = np.array([
    2.0, 4.0, 8.0, 12.0, 16.0, 20.0, 24.0,
])
_CP2OCP = np.array([
    1.3, 1.5, 1.8, 2.0, 2.1, 2.15, 2.2,
])
_AP2OCP = np.array([
    1.15, 1.15, 1.25, 1.4, 1.47, 1.51, 1.54,
])
# Figure 4.8.3.2-6B: centre of pressure against the nose angle.
_THETA = np.array([
    0.0, 10.0, 20.0, 30.0, 40.0,
])
_XCPOCR = np.array([
    0.65, 0.615, 0.585, 0.565, 0.55,
])
_X27A = np.array([
    0.0, 0.01, 0.02, 0.04, 0.08, 0.12, 0.16, 0.21,
])
# Figure 4.8.3.2-7A: blunt-nose centre-of-pressure shift (7 by 8).
_X17A = np.array([
    0.0, 5.0, 10.0, 15.0, 20.0, 25.0, 30.0,
])
_Y7A = np.array([
    0.0, -0.056, -0.072, -0.09, -0.11, -0.122, -0.13, -0.136,
    0.0, -0.05, -0.066, -0.082, -0.102, -0.114, -0.123, -0.129,
    0.0, -0.045, -0.06, -0.076, -0.096, -0.107, -0.115, -0.122,
    0.0, -0.042, -0.056, -0.072, -0.089, -0.1, -0.109, -0.117,
    0.0, -0.039, -0.052, -0.067, -0.084, -0.095, -0.103, -0.11,
    0.0, -0.036, -0.048, -0.062, -0.0795, -0.09, -0.098, -0.105,
    0.0, -0.032, -0.044, -0.058, -0.074, -0.085, -0.092, -0.098,
])
# Figure 5.5.1.1-6: side force at zero lift.
_X5516 = np.array([
    0.0, 0.05, 0.1, 0.2, 0.25, 0.3, 0.35, 0.4,
    0.45, 0.5, 0.55, 0.6,
])
_Y5516 = np.array([
    0.0, -0.005, -0.009, -0.021, -0.037, -0.064, -0.107, -0.172,
    -0.25, -0.34, -0.44, -0.55,
])
_X1558 = np.array([
    0.0, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08,
    0.09, 0.1, 0.105, 0.11, 0.115, 0.12,
])
_X2558 = np.array([
    0.0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35,
    0.4, 0.45, 0.5,
])
_Y558 = np.array([
    0.0, 0.0, -0.03, -0.09, -0.15, -0.24, -0.35, -0.47,
    -0.61, -0.78, -0.95, 0.0, 0.0, -0.02, -0.05, -0.09,
    -0.14, -0.2, -0.26, -0.35, -0.44, -0.54, 0.0, 0.0,
    -0.001, -0.02, -0.04, -0.08, -0.1, -0.15, -0.2, -0.25,
    -0.3, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.03, 0.05, 0.08, 0.1, 0.13, 0.19, 0.22, 0.0,
    0.0, 0.0, 0.04, 0.07, 0.11, 0.18, 0.24, 0.31,
    0.4, 0.5, 0.0, 0.0, 0.02, 0.07, 0.13, 0.2,
    0.29, 0.39, 0.5, 0.63, 0.8, 0.0, 0.0, 0.03,
    0.1, 0.19, 0.29, 0.4, 0.54, 0.7, 0.96, 1.1,
    0.0, 0.0, 0.04, 0.12, 0.22, 0.35, 0.5, 0.69,
    0.9, 1.11, 1.39, 0.0, 0.0, 0.05, 0.15, 0.25,
    0.4, 0.58, 0.79, 1.02, 1.3, 1.6, 0.0, 0.0,
    0.04, 0.12, 0.22, 0.35, 0.5, 0.69, 0.9, 1.11,
    1.39, 0.0, 0.0, 0.03, 0.1, 0.19, 0.29, 0.4,
    0.54, 0.7, 0.9, 1.1, 0.0, 0.0, 0.0, 0.05,
    0.1, 0.15, 0.21, 0.29, 0.38, 0.49, 0.6, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.02, 0.03, 0.04, 0.06,
    0.08, 0.1,
])
# Figures 5.5.1.2-8 and -9: side force with lift.
_X1559 = np.array([
    0.0, 10.0, 20.0, 30.0, 40.0, 50.0, 60.0,
])
_X2559 = np.array([
    0.0, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4,
    0.45, 0.5,
])
_Y559 = np.array([
    0.0, -0.04, -0.09, -0.155, -0.24, -0.342, -0.472, -0.62,
    -0.77, -0.95, 0.0, -0.02, -0.065, -0.11, -0.17, -0.244,
    -0.338, -0.44, -0.55, -0.67, 0.0, -0.01, -0.04, -0.07,
    -0.11, -0.151, -0.2, -0.26, -0.329, -0.4, 0.0, 0.0,
    -0.02, -0.023, -0.035, -0.05, -0.06, -0.08, -0.1, -0.12,
    0.0, 0.0, 0.01, 0.022, 0.04, 0.06, 0.08, 0.11,
    0.14, 0.18, 0.0, 0.01, 0.039, 0.07, 0.11, 0.159,
    0.22, 0.29, 0.37, 0.46, 0.0, 0.028, 0.065, 0.123,
    0.198, 0.279, 0.37, 0.48, 0.61, 0.75,
])
# Figures 5.5.2.1-8A/B: rolling moment at zero lift.
_X5528 = np.array([
    3.5, 4.0, 4.8, 6.0, 8.0, 9.5, 12.0, 13.5,
    16.0, 18.0, 22.0, 25.0, 28.0, 31.0, 34.0, 39.0,
    42.0, 45.0,
])
_Y5528 = np.array([
    -2.7, -2.4, -2.0, -1.6, -1.2, -1.0, -0.8, -0.7,
    -0.58, -0.52, -0.4, -0.36, -0.3, -0.27, -0.24, -0.2,
    -0.18, -0.16,
])
_X1558B = np.array([
    5.0, 10.0, 15.0, 20.0, 25.0, 30.0,
])
_X2558B = np.array([
    0.0, 0.01, 0.02, 0.04, 0.06, 0.08, 0.1, 0.12,
    0.14, 0.16, 0.18, 0.2,
])
_Y558B = np.array([
    0.0, 0.24, 0.3, 0.403, 0.477, 0.537, 0.585, 0.628,
    0.667, 0.704, 0.74, 0.77, 0.0, 0.092, 0.13, 0.176,
    0.21, 0.236, 0.257, 0.275, 0.294, 0.314, 0.33, 0.347,
    0.0, 0.054, 0.074, 0.105, 0.125, 0.139, 0.152, 0.168,
    0.179, 0.188, 0.199, 0.204, 0.0, 0.035, 0.05, 0.07,
    0.08, 0.09, 0.101, 0.11, 0.12, 0.13, 0.138, 0.144,
    0.0, 0.024, 0.033, 0.05, 0.059, 0.066, 0.071, 0.08,
    0.086, 0.09, 0.1, 0.102, 0.0, 0.019, 0.025, 0.036,
    0.042, 0.048, 0.052, 0.06, 0.063, 0.067, 0.072, 0.076,
])
_X5522A = np.array([
    0.0, 0.005, 0.01, 0.015, 0.02, 0.03, 0.035, 0.045,
    0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.11, 0.12,
    0.13, 0.14, 0.145,
])
_Y5522A = np.array([
    1.0, 0.84, 0.76, 0.7, 0.66, 0.61, 0.61, 0.66,
    0.7, 0.81, 0.94, 1.1, 1.24, 1.34, 1.42, 1.5,
    1.54, 1.6, 1.61,
])
_X5522B = np.array([
    0.0, 10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0,
])
_Y5522B = np.array([
    1.0, 0.94, 0.83, 0.69, 0.5, 0.32, 0.19, 0.1,
])
# Figures 5.5.2.2-12A/B and -13: rolling-moment increments.
_X55213 = np.array([
    0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.65,
    0.7, 0.8, 0.9, 1.0,
])
_Y55213 = np.array([
    0.0, 0.03, 0.09, 0.12, 0.16, 0.2, 0.26, 0.3,
    0.38, 0.53, 0.8, 0.125,
])
# Figure 5.5.3.1-6: centre of the side area.
_X5536 = np.array([
    0.0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07,
    0.08, 0.09, 0.1, 0.11, 0.12, 0.13, 0.14, 0.15,
    0.16, 0.165, 0.17,
])
_Y5536 = np.array([
    1.04, 1.04, 1.04, 1.04, 1.04, 1.04, 1.04, 1.04,
    1.04, 1.04, 1.04, 1.0, 0.93, 0.81, 0.655, 0.51,
    0.32, 0.18, 0.0,
])


def _tlinex(x1, x2, flat, q1, q2, l1, l2, u1, u2) -> float:
    """TLINEX on a source ``Y(NX2,NX1)`` table given flat."""
    y = np.asarray(flat, dtype=float).reshape(len(x1), len(x2)).T
    return float(tlinex(x1, x2, y, q1, q2, l1, l2, u1, u2))


def _value(x, y, query, lower=0, upper=0) -> float:
    return float(tbfunx(x, y, query, lower, upper)[0])


def _at_flight_angles(alpha, alphap, curve, alpha0):
    """The source's pass back to flight angles: the curve itself when the
    zero-lift angle is zero, else TBFUNX on the shifted angles."""
    if alpha0 == 0.0:
        return np.array(curve, dtype=float)
    return np.array([_value(alphap, curve, a) for a in alpha])


def calculate_loarwb(alpha_deg: Sequence[float], lbin: Mapping[int, object],
                     mach: float, reynolds_per_length: float,
                     roughness: float, stale_xocrb: float = 0.0
                     ) -> Dict[str, object]:
    """Translate LOARWB: low-aspect-ratio wing-body forces and moments.

    Args:
        alpha_deg: ``FLC(23)`` onward.
        lbin: The ``LARWB`` inputs, ``LBIN(1..21)`` one-based: 1 ``ZB``,
            2 ``SREF``, 3 ``DELTEP``, 4 ``SFRONT``, 5 ``AR``, 6 ``R3LEOB``,
            7 ``DELTAL``, 8 ``L``, 9 ``SWET``, 10 ``PERBAS``, 11 ``SBASE``,
            12 ``HB``, 13 ``BB``, 14 ``BLF`` (bool), 15 ``XCG``, 16
            ``THETAD``, 17 ``ROUNDN`` (bool), 18 ``SBS``, 19 ``SBSLB``,
            20 ``XCENSB``, 21 ``XCENW``.  ``DELTEP`` or ``R3LEOB`` is
            ``UNUSED`` when not given.
        mach, reynolds_per_length: ``FLC(I+2)``, ``FLC(I+42)``.
        roughness: ``ROUGFC``.
        stale_xocrb: ``LB(119)`` before the call; see Notes.

    Returns:
        ``cl``, ``cd``, ``cn``, ``ca``, ``cm`` (``BW(21)``, ``(1)``,
        ``(61)``, ``(81)``, ``(41)`` onward), ``xcp`` (``BW(201)``),
        ``kyb``, ``knb``, ``klb`` (``BW(141)``, ``(161)``, ``(181)``) at
        the flight angles, the ``/SUPDW/`` words as ``lb`` (one-based), and
        ``xcg`` (written to ``SYNTSS(1)``).

    Notes:
        The blunt-nose shift ``XOCRB`` (``LB(119)``) is computed only for a
        round nose, and otherwise keeps whatever the block held; it is an
        input.  The normal-force branch is chosen by whether ``DELTEP`` is
        given, the lateral ones by whether ``R3LEOB`` is; with ``DELTEP``
        given the axial-force parameter ``2*R3LEOB*AR/SFOSR`` is formed
        from the ``UNUSED`` radius.  ``DX`` and ``KCCA20`` share ``LB(22)``,
        which ends as ``KCCA20``.  Kept.
    """
    v = {int(k): val for k, val in lbin.items()}
    zb = float(v[1])
    sref = float(v[2])
    deltep = float(v[3])
    sfront = float(v[4])
    ar = float(v[5])
    r3leob = float(v[6])
    deltal = float(v[7])
    length = float(v[8])
    swet = float(v[9])
    perbas = float(v[10])
    sbase = float(v[11])
    hb = float(v[12])
    bb = float(v[13])
    blf = bool(v[14])
    xcg = float(v[15])
    thetad = float(v[16])
    roundn = bool(v[17])
    sbs = float(v[18])
    sbslb = float(v[19])
    xcensb = float(v[20])
    xcenw = float(v[21])

    alpha = np.asarray(alpha_deg, dtype=float)
    lb: Dict[int, float] = {}

    alpha0 = math.atan(zb / length) * RAD
    alphap = alpha - alpha0
    cnac0 = (PI * ar / 2.0 + 2.0 * sfront / sref) * (4.0 / (4.0 + ar))
    cnc20 = 2.195 * ((ar + 0.61) / (ar + 4.0))
    if deltep == UNUSED:
        acna0 = _value(_RLEOB, _RCNA0, r3leob, 0, 2)
        acna20 = _value(_RLEOB, _RCNA20, r3leob, 0, 2)
        z = _value(_RLEOB, _ZR, r3leob, 0, 2)
        cn20 = cnc20 * ((1.0 + math.cos(deltal / RAD)) / 2.0) * acna20
    else:
        acna0 = _value(_DELTAE, _SCNA0, deltep, 0, 2)
        acna20 = _value(_DELTAE, _SCNA20, deltep, 0, 2)
        z = _value(_DELTAE, _ZS, deltep, 0, 2)
        cn20 = cnc20 * acna20
    cna0 = cnac0 * acna0
    alpapr = alphap / RAD
    cnp = (cna0 * alpapr + 8.21 * (cn20 - 0.349 * cna0) *
           (z + 2.81 * (1.0 - z) * alpapr) * alpapr**2)

    shapep = 2.0 * sbase / (PI * length * (hb + bb))
    cpbops = _value(_SBOLHB, _CPBOP, shapep, 0, 2)
    cpb0 = cpbops * perbas * sbase / (2.0 * math.sqrt(PI * sbase) * sref)
    rn = reynolds_per_length * length
    lok = 12.0 * length / roughness
    cept = _value(_X27M, _X27I, mach, 0, 0)
    rl = lok**1.0482 * 10.0**cept
    if rl < rn:
        rn = rl
    cf = fig26(rn, mach)
    cx0p = -cf * swet / sref + cpb0
    sfosr = sfront / sref
    geopar = 2.0 * r3leob * ar / sfosr
    dcxcxc = _tlinex(_X110, _X210, _Y10, geopar, sfosr, 0, 0, 2, 2)
    acx = dcxcxc * 0.349 * ((ar + 2) / (ar + 4))
    shapeb = bb**2 / (hb * sbase**.50)
    cp20o0 = _value(_B2OHS, _AP2OCP if blf else _CP2OCP, shapeb, 2, 2)
    acpb0 = cpb0 * (cp20o0 - 1.0)
    cxp = -(cx0p + acx * cnp * np.sqrt(np.abs(alphap) / 20.0) +
            acpb0 * (alphap / 20.0)**2)

    cm0 = 0.0003333 * alpha0
    xocrb = float(stale_xocrb)
    if roundn:
        bluntp = 1.0 - (4.0 * math.tan(thetad / RAD)) / ar
        xocrb = _tlinex(_X17A, _X27A, _Y7A, thetad, bluntp, 0, 0, 1, 2)
        lb[117] = bluntp
    xocrd = _value(_THETA, _XCPOCR, thetad, 0, 2)
    xocrt = 0.1020 * sfosr
    xcpoc = xocrd + xocrb + xocrt
    dx = xcg / length - xcpoc
    cm_alphap = cm0 + dx * cnp

    cn = _at_flight_angles(alpha, alphap, cnp, alpha0)
    ca = _at_flight_angles(alpha, alphap, cxp, alpha0)
    cm = _at_flight_angles(alpha, alphap, cm_alphap, alpha0)
    sin_alpha = np.sin(alpha / RAD)
    cos_alpha = np.cos(alpha / RAD)
    cl = cn * cos_alpha - ca * sin_alpha
    cd = ca * cos_alpha + cn * sin_alpha
    cn = np.where(cn == 0.0, 1.0e-20, cn)
    xcp = cm / cn

    # Lateral: Section 5.5.
    arg6 = sbs / sref
    kybno = _value(_X5516, _Y5516, arg6, 0, 2)
    xcpxc = _value(_X5536, _Y5536, sbslb / sbs, 0, 0)
    xcpp = xcpxc * xcensb
    akn1 = -.006 * (xcg - xcenw) / (RAD * bb)
    akn2 = (xcg - xcpp) / (RAD * bb)
    dklcno = _value(_X5528, _Y5528, thetad, 2, 2)
    dklcnb = 0.0
    if roundn:
        dklcnb = _tlinex(_X1558B, _X2558B, _Y558B, thetad,
                         1. - (4. * math.tan(thetad / RAD) / ar), 2, 0, 2, 2)
    klbcno = 1. / RAD * (dklcno + dklcnb)
    if r3leob == UNUSED:
        kycn20 = _tlinex(_X1559, _X2559, _Y559, deltep, arg6, 0, 0, 1, 1)
        kckcc2 = _value(_X5522B, _Y5522B, deltep, 0, 2)
    else:
        kycn20 = _tlinex(_X1558, _X2558, _Y558, r3leob, arg6, 0, 0, 0, 2)
        kckcc2 = _value(_X5522A, _Y5522A, r3leob, 0, 2)
    dkckcc = _value(_X55213, _Y55213, arg6, 0, 0)
    kcca20 = 1. / (1. + .152 / math.tan(thetad / RAD)) * dklcno
    bpart = kckcc2 + dkckcc
    kyb = (1. / RAD) * (kybno + kycn20 * cnp**2)
    knb = (kyb + .006) * akn2 + akn1
    klb = cnp / RAD * (dklcno * (1. - (alphap / 20.)**2) +
                       (bpart * kcca20 + dklcnb) * (alphap / 20.)**2)

    lb.update({1: alpha0, 22: kcca20, 23: dkckcc, 24: kckcc2, 25: kycn20,
               26: klbcno, 27: dklcnb, 28: cnac0, 29: cnc20, 30: acna0,
               31: acna20, 32: z, 33: cn20, 34: cna0, 75: shapep,
               76: cpbops, 77: dklcno, 78: 0.0, 79: xcpxc, 80: kybno,
               83: cpb0, 84: rn, 85: lok, 86: cf, 87: cx0p, 88: sfosr,
               89: geopar, 90: dcxcxc, 91: acx, 92: shapeb, 93: cp20o0,
               94: acpb0, 115: cm0, 116: xcpoc, 118: xocrd, 119: xocrb,
               120: xocrt})
    del lb[78]
    for j in range(len(alpha)):
        lb[2 + j], lb[35 + j], lb[55 + j] = alphap[j], alpapr[j], cnp[j]
        lb[95 + j], lb[121 + j] = cxp[j], cm_alphap[j]
        lb[141 + j], lb[161 + j], lb[181 + j] = kyb[j], knb[j], klb[j]
    return {
        'cl': cl, 'cd': cd, 'cn': cn, 'ca': ca, 'cm': cm, 'xcp': xcp,
        'kyb': _at_flight_angles(alpha, alphap, kyb, alpha0),
        'knb': _at_flight_angles(alpha, alphap, knb, alpha0),
        'klb': _at_flight_angles(alpha, alphap, klb, alpha0),
        'lb': lb, 'xcg': xcg, 'dx': float(dx), 'method': 'legacy_loarwb'}


def calculate_m14o16(alpha_deg: Sequence[float], lbin: Mapping[int, object],
                     mach: float, reynolds_per_length: float,
                     roughness: float, stale_xocrb: float = 0.0
                     ) -> Dict[str, object]:
    """Translate M14O16: LOARWB, then CN, CA and the slopes from CL and CD.

    ``BW(61)``/``(81)`` are rebuilt from CL and CD (recovering LOARWB's own
    CN and CA, bar the ``1e-20`` guard), and ``BW(101)``/``(121)`` onward are
    TBFUNX's derivatives of CL and CM over the flight angles.  Arguments
    are those of :func:`calculate_loarwb`.
    """
    r = calculate_loarwb(alpha_deg, lbin, mach, reynolds_per_length,
                         roughness, stale_xocrb)
    alpha = np.asarray(alpha_deg, dtype=float)
    cos_alpha = np.cos(alpha / RAD)
    sin_alpha = np.sin(alpha / RAD)
    r['cn'] = r['cl'] * cos_alpha + r['cd'] * sin_alpha
    r['ca'] = r['cd'] * cos_alpha - r['cl'] * sin_alpha
    r['cla'] = np.array([
        tbfunx(alpha, r['cl'], angle_deg, 0, 0)[1] for angle_deg in alpha
    ])
    r['cma'] = np.array([
        tbfunx(alpha, r['cm'], angle_deg, 0, 0)[1] for angle_deg in alpha
    ])
    r['method'] = 'legacy_m14o16'
    return r
