"""
Subsonic wing-body buildup: WBAERO and the routines it calls.

``WBAERO`` combines a lifting surface with the body in four steps, each a
separate legacy routine:

- ``WBDRAG``: zero-lift drag with the Figure 4.3.3.1-37 interference factor,
  plus both components' drag due to lift.
- ``WBLIFT``: the Figure 4.3.1.2-10 carryover factors K_W(B) and K_B(W),
  the -12 incidence factors, and the body-vortex lift increment.
- ``WBCM``: the Section 4.3.2 aerodynamic centre, the zero-lift moment
  (``WBCM0``'s regression where it applies), and the moment buildup at each
  angle.
- ``WBAERO`` itself: the lift- and moment-curve slopes by TBFUNX, and the
  normal and axial forces.

``EXSUBT`` (experimental data substitution) is an input-file operation and
has no numerical content here.  ``EXIT``, which WBCM calls when its ellipse
fit fails, only closes files; the translation raises instead.

The tables were extracted from the source DATA statements by parsing, with
``tools/fortran_data.py``, rather than by hand.

Reference: datcom-legacy/datcom_2000/wbaero.f, wbdrag.f, wblift.f, wbcm.f,
wbcm0.f
"""

import numpy as np
from typing import Dict, Optional, Sequence
import logging

from pydatcom.aerodynamics.cdrag import STRAIGHT_TAPERED as STRAIGHT_TAPERED_WB
from pydatcom.aerodynamics.tablec import calculate_tablec
from pydatcom.aerodynamics.tbsub import calculate_tbsub
from pydatcom.aerodynamics.tbsup import calculate_tbsup
from pydatcom.aerodynamics.tbtrn import calculate_tbtrn
from pydatcom.interactions.body_vortex import calculate_bodowg, getmax
from pydatcom.utils.constants import PI, RAD, UNUSED
from pydatcom.utils.legacy_numeric import tbfunx
from pydatcom.utils.legacy_tables import tlinex

logger = logging.getLogger(__name__)

# The source's "not available" marker for a moment, 2*UNUSED.
NOT_AVAILABLE = 2.0 * UNUSED

# WBDRAG: Figure 4.3.3.1-37, the wing-body interference factor R_WB.
# Y37 is (19,7): Reynolds number fastest, one run per Mach number.
_X137 = np.array([
    0.25, 0.4, 0.6, 0.7, 0.8, 0.85, 0.9,
])
_X237 = np.array([
    3000000.0, 5000000.0, 7000000.0, 10000000.0, 15000000.0, 20000000.0, 25000000.0, 30000000.0, 35000000.0, 40000000.0,
    45000000.0, 50000000.0, 60000000.0, 70000000.0, 80000000.0, 100000000.0, 150000000.0, 200000000.0, 700000000.0,
])
_Y37 = np.array([
    1.063, 1.07, 1.073, 1.076, 1.072, 1.066, 1.057, 1.045, 1.025, 0.9935,
    0.965, 0.9515, 0.939, 0.9335, 0.931, 0.928, 0.923, 0.9225, 0.954, 1.02,
    1.023, 1.028, 1.036, 1.05, 1.058, 1.058, 1.05, 1.032, 1.018, 1.008,
    1.001, 0.992, 0.9875, 0.9845, 0.98, 0.977, 0.9755, 0.975, 0.98, 0.984,
    0.989, 0.9965, 1.008, 1.02, 1.0325, 1.0375, 1.035, 1.0315, 1.028, 1.023,
    1.019, 1.0155, 1.015, 1.015, 1.015, 1.015, 1.015, 0.955, 0.96, 0.9647,
    0.9725, 0.983, 0.995, 1.0085, 1.013, 1.014, 1.0145, 1.015, 1.015, 1.015,
    1.015, 1.015, 1.015, 1.015, 1.015, 1.015, 0.925, 0.93, 0.934, 0.942,
    0.953, 0.9655, 0.9775, 0.9885, 0.9965, 1.0025, 1.008, 1.011, 1.0145, 1.015,
    1.015, 1.015, 1.015, 1.015, 1.015, 0.9025, 0.907, 0.912, 0.919, 0.931,
    0.9425, 0.954, 0.969, 0.982, 0.992, 0.9985, 1.0035, 1.011, 1.014, 1.015,
    1.015, 1.015, 1.015, 1.015, 0.87, 0.8715, 0.8775, 0.884, 0.896, 0.909,
    0.9225, 0.94, 0.957, 0.9725, 0.9865, 0.9935, 1.0065, 1.012, 1.015, 1.015,
    1.015, 1.015, 1.015,
])

# WBLIFT: Figures 4.3.1.2-10A/B (K_W(B), K_B(W)), -12A1/A2 (the incidence
# factors k_W(B), k_B(W)) and -12C (a wing running the body's length).
_X10A = np.array([
    0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0,
])
_Y10A = np.array([
    1.0, 1.08, 1.16, 1.26, 1.36, 1.46, 1.56, 1.67, 1.78, 1.89, 2.0,
])
_X10B = np.array([
    0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0,
])
_Y10B = np.array([
    0.0, 0.13, 0.29, 0.45, 0.62, 0.8, 1.0, 1.22, 1.45, 1.7, 2.0,
])
_X12A1 = np.array([
    0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0,
])
_Y12A1 = np.array([
    1.0, 0.97, 0.95, 0.94, 0.94, 0.94, 0.94, 0.95, 0.96, 0.98, 0.99,
])
_X12A2 = np.array([
    0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0,
])
_Y12A2 = np.array([
    0.0, 0.11, 0.21, 0.31, 0.41, 0.51, 0.6, 0.7, 0.8, 0.9, 1.0,
])
_X12C = np.array([
    0.0, 0.1, 0.15, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.83, 0.9, 1.0,
])
_Y12C = np.array([
    1.0, 1.0, 0.999, 0.99, 0.98, 0.965, 0.95, 0.933, 0.92, 0.92, 0.928, 0.95, 1.0,
])

# WBLIFT: Figure 4.3.1.4-12B/C, the wing-body CLMAX and stall-angle
# ratios.  Each is (7,5): d/b fastest, one run per A(160) value.
_XA12 = np.array([
    1.0, 2.0, 4.0, 6.0, 12.0,
])
_XB12 = np.array([
    0.0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3,
])
_Y412B = np.array([
    1.0, 0.998, 0.995, 0.989, 0.982, 0.968, 0.95,
    1.0, 0.995, 0.989, 0.978, 0.963, 0.937, 0.917,
    1.0, 0.992, 0.982, 0.965, 0.942, 0.9, 0.864,
    1.0, 0.967, 0.932, 0.911, 0.907, 0.945, 1.038,
    1.0, 0.983, 0.968, 0.962, 0.977, 1.027, 1.134,
])
_Y412C = np.array([
    1.0, 1.013, 1.017, 1.013, 1.0, 0.982, 0.956,
    1.0, 1.003, 1.001, 0.992, 0.973, 0.937, 0.894,
    1.0, 0.989, 0.973, 0.944, 0.898, 0.822, 0.745,
    1.0, 0.98, 0.943, 0.879, 0.758, 0.648, 0.613,
    1.0, 0.942, 0.845, 0.688, 0.594, 0.552, 0.532,
])

# WBCM: Figure 4.3.2.2-36B, and the carryover aerodynamic-centre curve
# its source labels only 21C.
_X38B = np.array([
    0.0, 1.0, 1000000.0,
])
_Y38B = np.array([
    0.0, 0.5, 0.5,
])
_X21C = np.array([
    0.0, 0.025, 0.05, 0.075, 0.1, 0.15, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8,
])
_Y21C = np.array([
    0.0, 0.056, 0.101, 0.13, 0.152, 0.19, 0.22, 0.266, 0.301, 0.33, 0.348, 0.365, 0.375,
])


def _div(a: float, b: float) -> float:
    """IEEE division, as the source gets it: inf or NaN rather than an
    exception, so a degenerate input propagates the way it does there."""
    with np.errstate(divide='ignore', invalid='ignore'):
        return float(np.float64(a) / np.float64(b))


def calculate_wbdrag(mach: float, reynolds_per_length: float,
                     body_length: float, wing_cd0: float,
                     body_cd_friction: float, body_cd_base: float,
                     body_cdl: Sequence[float], wing_cdl: Sequence[float],
                     experimental: Optional[Dict[str, object]] = None
                     ) -> Dict[str, object]:
    """Translate WBDRAG: subsonic wing-body drag.

    ``CD0 = (CD0_wing + CD0_body,friction) * R_WB + CD_body,base`` with R_WB
    from Figure 4.3.3.1-37 at the body-length Reynolds number, and at each
    angle ``CD = CD0 + CDL_body + CDL_wing``.

    Args:
        mach: ``B(1)``.
        reynolds_per_length: ``A(129)``, the flight Reynolds number per unit
            length.
        body_length: ``BD(1)``.
        wing_cd0: ``D(20)``.
        body_cd_friction: ``BD(59)``, body friction and pressure drag.
        body_cd_base: ``BD(60)``, body base drag.
        body_cdl: ``BD(215)`` onward, the body drag due to lift.
        wing_cdl: ``D(36)`` onward, the wing drag due to lift.
        experimental: Optional ``{'kbody', 'kwing', 'body_cd', 'wing_cd'}``.
            With either flag set, an angle at which both supplied drag
            values are available takes their plain sum instead.

    Returns:
        Dictionary with ``cd`` (``BW(1)`` onward), ``cd0`` (``WB(17)``),
        ``interference`` (``WB(18)``) and ``reynolds`` (``WB(19)``).
    """
    reynolds = reynolds_per_length * body_length
    interference = tlinex(_X137, _X237, _Y37.reshape(7, 19).T, mach,
                          reynolds, 2, 2, 2, 1)
    cd0 = (wing_cd0 + body_cd_friction) * interference + body_cd_base
    cd = [
        cd0 + body_dcl + wing_dcl
        for body_dcl, wing_dcl in zip(body_cdl, wing_cdl)
    ]
    if experimental and (experimental.get('kbody') or
                         experimental.get('kwing')):
        for angle_slot, (b, w) in enumerate(zip(experimental['body_cd'],
                                                experimental['wing_cd'])):
            if b != UNUSED and w != UNUSED:
                cd[angle_slot] = b + w
    return {
        'cd': np.array(cd), 'cd0': float(cd0),
        'interference': float(interference), 'reynolds': float(reynolds),
        'method': 'legacy_wbdrag',
    }


def calculate_wblift(alpha_deg: Sequence[float],
                     local_alpha_deg: Sequence[float],
                     body_diameter: float, incidence: float,
                     semispan: float, alpha_zero_lift: float,
                     surface_alone: Dict[str, object],
                     body_cl: Sequence[float], body_cla: float,
                     vortex: Dict[str, object],
                     stall: Optional[Dict[str, float]] = None,
                     stale: Optional[Dict[str, float]] = None
                     ) -> Dict[str, object]:
    """Translate WBLIFT: wing-body lift with carryover and body vortex.

    ``CL = CL_body + (K_W(B)+K_B(W)) * CL_wing(alpha_eff) + vortex`` with
    ``alpha_eff = alpha - alpha_0L + (k_W(B)+k_B(W))/(K_W(B)+K_B(W)) * i``,
    looked up on the wing curve by TBFUNX end modes 1 and 1.

    Args:
        alpha_deg: ``FLC(23)`` onward.
        local_alpha_deg: ``B(23)`` onward, the surface's local angles, which
            the vortex term multiplies.
        body_diameter: ``DB``, the body width at the surface.
        incidence: ``AIW``, degrees.
        semispan: ``WINGIN(4)``, the theoretical semispan.
        alpha_zero_lift: ``B(49)``.
        surface_alone: ``{'cl', 'cla'}``, ``WING(21)`` onward and
            ``WING(101)``.
        body_cl: ``BODY(21)`` onward.
        body_cla: ``BODY(101)``.
        vortex: ``{'ratio', 'ivbw', 'go2pav'}``: ``FACT(1)`` and BODOWG's
            ``FACT(2)`` and ``FACT(22)`` arrays.
        stall: ``{'a160', 'clmax', 'alpha_clmax'}`` for the Figure
            4.3.1.4-12 ratios, which apply when ``d/b <= 0.30``.
        stale: The ``WB(2)`` to ``WB(5)``, ``WB(7)`` and ``WB(8)`` a
            previous case left, keyed ``kwb``, ``kbw``, ``wb4``, ``wb5``,
            ``kwb_incidence`` and ``kbw_incidence``, and zero by default.
            See Notes.

    Returns:
        Dictionary with ``cl`` (``BW(21)`` onward), ``cla`` (``BW(101)``)
        and the ``WB`` factors: ``kwb`` (2), ``kbw`` (3), ``wb4``, ``wb5``,
        and with incidence ``kwb_incidence`` (7), ``kbw_incidence`` (8),
        ``wb9``, ``wb10``, ``wb11``; ``wb20`` to ``wb23`` for a small body.

    Notes:
        Above ``d/b = 0.8`` the source takes Figure 4.3.1.2-12C for the
        slope but skips the lookups of K_W(B) and K_B(W), and still
        multiplies the wing lift by ``WB(2)+WB(3)``.  It therefore reads
        whatever an earlier case left there, as WBCM does with the
        carryover slopes ``WB(4)`` and ``WB(5)``, also skipped.  Likewise ``WB(7)`` and
        ``WB(8)`` are left untouched with no incidence, where they are
        multiplied by zero.  Those values are taken from ``stale`` and the
        branch is flagged; with fresh zeros the source's own arithmetic
        gives ``0/0``, which propagates as NaN here as it does there.
    """
    stale = dict({'kwb': 0.0, 'kbw': 0.0, 'wb4': 0.0, 'wb5': 0.0,
                  'kwb_incidence': 0.0, 'kbw_incidence': 0.0},
                 **(stale or {}))
    wing_cla = float(surface_alone['cla'])
    ratio = body_diameter / (2.0 * semispan)
    result = {'ratio': float(ratio), 'method': 'legacy_wblift'}

    kwb, kbw = stale['kwb'], stale['kbw']
    kwb_inc, kbw_inc = stale['kwb_incidence'], stale['kbw_incidence']
    if ratio > 0.80:
        extended, _ = tbfunx(_X12C, _Y12C, ratio, 0, 0)
        cla = extended * wing_cla
        result.update({'wb4': float(stale['wb4']),
                       'wb5': float(stale['wb5']),
                       'source_defect': 'ratio_above_0.8_reads_stale_kwb_kbw'})
    else:
        kwb, _ = tbfunx(_X10A, _Y10A, ratio, 0, 0)
        kbw, _ = tbfunx(_X10B, _Y10B, ratio, 0, 0)
        result.update({'wb4': kwb * wing_cla, 'wb5': kbw * wing_cla})
        cla = body_cla + result['wb4'] + result['wb5']
        if incidence != 0.0:
            kwb_inc, _ = tbfunx(_X12A1, _Y12A1, ratio, 0, 0)
            kbw_inc, _ = tbfunx(_X12A2, _Y12A2, ratio, 0, 0)
            result.update({
                'wb9': kwb_inc * wing_cla,
                'wb10': kbw_inc * wing_cla,
                'wb11': (kwb_inc + kbw_inc) * wing_cla,
            })
    result.update({
        'kwb': float(kwb),
        'kbw': float(kbw),
        'kwb_incidence': float(kwb_inc),
        'kbw_incidence': float(kbw_inc),
        'cla': float(cla),
    })

    alpha = np.asarray(alpha_deg, dtype=float)
    lookup_angles = alpha + incidence - alpha_zero_lift
    incidence_shift = _div(kwb_inc + kbw_inc, kwb + kbw) * incidence
    cl = []
    for angle_slot, angle_deg in enumerate(alpha):
        wing_cl, _ = tbfunx(
            lookup_angles, surface_alone['cl'],
            angle_deg - alpha_zero_lift + incidence_shift, 1, 1,
        )
        cl.append(
            body_cl[angle_slot] + (kwb + kbw) * wing_cl +
            vortex['ivbw'][angle_slot] * vortex['go2pav'][angle_slot] *
            local_alpha_deg[angle_slot] * vortex['ratio'] * wing_cla
        )
    result['cl'] = np.array(cl)

    if ratio <= 0.30 and stall is not None:
        wb20 = tlinex(_XA12, _XB12, _Y412B.reshape(5, 7).T,
                      stall['a160'], ratio, 0, 0, 0, 0)
        wb21 = tlinex(_XA12, _XB12, _Y412C.reshape(5, 7).T,
                      stall['a160'], ratio, 0, 0, 0, 0)
        result.update({'wb20': float(wb20), 'wb21': float(wb21),
                       'wb22': float(wb20 * stall['clmax']),
                       'wb23': float(wb21 * stall['alpha_clmax'])})
    return result


def calculate_wbcm0(aspect_ratio: float, tan_le: float, tovc: float,
                    nose_length: float, afterbody_length: float,
                    taper_ratio: float, leading_edge_radius: float,
                    twist: float, ycm: float, cld: float,
                    reynolds: float, tr: float, wing_height: float,
                    vt: float, hd: float, body_radius_ratio: float,
                    mach: float) -> Optional[float]:
    """Translate WBCM0: the regression zero-lift moment of a wing-body.

    Returns ``None`` outside the regression's range, where the source
    leaves its caller's value untouched.  The Reynolds number is clamped
    to 8e5 to 8e6 before use, as in the source.

    Args follow the source's names: ``AR``, ``TANLE``, ``TOVC``, ``LN``,
    ``LA``, ``TAPR``, ``LER``, ``TWIST``, ``YCM``, ``CLD``, ``RN``, ``TR``,
    ``WL``, ``VT``, ``HD``, ``DB``, ``MACH``.
    """
    inside = (mach <= 2.5 and 1.6 <= aspect_ratio <= 6.0 and
              0.0 <= tan_le <= 2.74748 and 0.025 <= tovc <= 0.100 and
              2.2 <= nose_length <= 8.4 and 0.3 <= afterbody_length <= 5.6 and
              0.0 <= taper_ratio <= 1.0 and
              0.0 <= leading_edge_radius <= 0.015 and
              -9.4 <= twist <= UNUSED and 0.0 <= ycm <= 0.0263 and
              0.0 <= cld <= 0.45)
    if not inside:
        return None
    reynolds = min(max(reynolds, 8.0e5), 8.0e6)
    c = calculate_tablec(mach)['c']
    return float(c[0] + c[1] / aspect_ratio + c[2] * aspect_ratio +
                 c[3] * tan_le + c[4] * tovc + c[5] * nose_length +
                 c[6] * afterbody_length + c[7] * taper_ratio +
                 c[8] * taper_ratio**2 + c[9] * tr +
                 c[10] * leading_edge_radius + c[11] * twist + c[12] * ycm +
                 c[13] * cld + c[14] * wing_height + c[15] * vt +
                 c[16] * hd + c[17] * body_radius_ratio +
                 c[18] * reynolds / 1.0e6)


def calculate_wbcm(alpha_deg: Sequence[float],
                   local_alpha_deg: Sequence[float],
                   surface: Dict[str, float],
                   surface_alone: Dict[str, object],
                   body: Dict[str, object],
                   synthesis: Dict[str, float],
                   lift: Dict[str, object],
                   vortex: Dict[str, object],
                   flight: Dict[str, float],
                   cbarr: float, c6: float) -> Dict[str, object]:
    """Translate WBCM: wing-body pitching moment.

    Args:
        alpha_deg: ``FLC(23)`` onward.
        local_alpha_deg: ``B(23)`` onward.
        surface: ``WINGIN`` entries ``sspn`` (4), ``chrdr`` (6),
            ``twista`` (11), ``tovc`` (16), ``ler`` (62), ``ycm`` (93) and
            ``cld`` (94); ``A`` entries ``a7``, ``a10``, ``a27``, ``a38``,
            ``a44``, ``a62``, ``a80``, ``a118``, ``a120``, ``a122``,
            ``a173``; ``B`` entries ``beta`` (2), ``cm0`` (47) and
            ``alpha_zero_lift`` (49).
        surface_alone: ``cm``, ``cn``, ``ca`` (``WINGC(1)``, ``(21)``,
            ``(41)`` onward: ``WING(41)``, ``(61)``, ``(81)``), ``cla``
            (``WING(101)``) and ``cma``, the caller's ``CMA`` argument.
        body: ``cm`` (``BODY(41)`` onward), ``cla`` and ``cma``
            (``BODY(101)`` and ``BODY(121)`` onward), ``cm0`` (``CMOB``),
            ``alpha_zero_lift`` (``BAL0``), ``length`` ``BD(1)`` and
            ``max_area`` ``BD(3)``.
        synthesis: ``xcg``, ``xw``, ``zw``, ``zcg``, ``incidence`` (``ALI``),
            ``cos_incidence`` (``COSAIW``) and ``body_diameter`` (``DB``).
        lift: WBLIFT's result: ``cl``, ``kwb``, ``kbw``, ``kwb_incidence``,
            ``kbw_incidence``, ``wb4``, ``wb5``.
        vortex: ``{'ratio', 'ivbw', 'go2pav'}``, as for WBLIFT.
        flight: ``mach`` (``FLC(II+2)``), ``reynolds_per_length``
            (``FLC(II+42)``) and ``tr`` (``FLC(96)``).
        cbarr: ``CBARR``.
        c6: ``C(6)`` of the surface's ``/WHAERO/`` block.

    Returns:
        Dictionary with ``cm`` (``BW(41)`` onward), ``cm0`` (``WB(16)``),
        ``alpha_zero_lift`` (``ALOWB``), and ``WB(12)`` to ``WB(15)``.

    Raises:
        ValueError: If the ellipse fit has no real root, where the source
            prints an error, "exits" into a routine that only closes files,
            and carries on with the square root of a negative number.

    Notes:
        WBCM0 is called with ``TWIST = WINGIN(11)/RAD``, in radians, and its
        term ``C11*TWIST`` is therefore in radians, as the companion
        regression WBCDL's ``-B12*TWIST/RAD`` is.  But WBCM0 range-checks
        that radian value against ``-9.4``, a degree bound (WBCDL checks
        the degree value), so the check admits washout far beyond the
        regression's data.  The source form is kept.
    """
    a7, a10 = float(surface['a7']), float(surface['a10'])
    sspn = float(surface['sspn'])
    db = float(synthesis['body_diameter'])
    incidence = float(synthesis['incidence'])

    temp0 = 0.25 * a7 * (1.0 + float(surface['a27'])) * float(surface['a38'])
    wb15 = 0.50
    if temp0 < 1.0:
        wb15, _ = tbfunx(_X38B, _Y38B, temp0, 0, 0)
    ratio = db / (2.0 * sspn)
    brac, _ = tbfunx(_X21C, _Y21C, ratio, 0, 0)
    temp4 = 0.25 + (2.0 * sspn - db) / (2.0 * a10) * float(surface['a44']) * brac

    arg = float(surface['beta']) * a7
    if arg >= 4.0:
        wb14 = temp4
    else:
        # The ellipse through (0, WB(15)) and (4, TEMP4).
        endpoint_gap = abs(temp4 - wb15)
        linear_coeff = -2.0 * wb15
        constant_term = (wb15 * wb15 - endpoint_gap ** 2 +
                         ((endpoint_gap / 4.0) ** 2) * (arg - 4.0) ** 2)
        discriminant = linear_coeff ** 2 - 4.0 * constant_term
        if not discriminant > 0.0:
            raise ValueError("WBCM: ellipse curve fit in error")
        root = np.sqrt(discriminant)
        wb14 = ((-linear_coeff - root) / 2.0 if temp4 < wb15 else
                (-linear_coeff + root) / 2.0)
    wb13 = wb14 * a10 / cbarr

    dxcg = float(synthesis['xcg']) - (float(synthesis['xw']) + .50 * db *
                                      float(surface['a62']) *
                                      float(synthesis['cos_incidence']))
    alpha_zero, _ = tbfunx(lift["cl"], alpha_deg, 0.0, 1, 1, ordered=False)
    alpha_zero_surface = float(surface['alpha_zero_lift'])
    cm0_surface = float(surface['cm0'])
    body_cma = np.asarray(body['cma'], dtype=float)
    cmowb = (float(body['cm0']) + cm0_surface +
             body_cma[0] * (alpha_zero - float(body['alpha_zero_lift'])) +
             float(surface_alone['cma']) *
             (alpha_zero - alpha_zero_surface + incidence))

    xw = float(synthesis['xw'])
    diameter = 2.0 * np.sqrt(float(body['max_area']) / PI)
    regression = calculate_wbcm0(
        float(surface['a120']), float(surface['a38']), float(surface['tovc']),
        (xw + 0.5 * db * float(surface['a38'])) / db,
        float(body['length']) / db - (xw + float(surface['chrdr']) +
                                      0.5 * db * float(surface['a80'])) / db,
        float(surface['a118']), float(surface['ler']),
        float(surface['twista']) / RAD, float(surface['ycm']),
        float(surface['cld']),
        float(flight['reynolds_per_length']) * float(surface['a122']),
        float(flight['tr']), 0.5 + float(synthesis['zw']) / diameter, 0.0,
        0.5, 0.5 * diameter / sspn, float(flight['mach']))
    if regression is not None:
        cmowb = regression
    if abs(cmowb) <= UNUSED:
        cmowb = 0.0

    body_cla = np.asarray(body['cla'], dtype=float)
    wb4, wb5 = float(lift['wb4']), float(lift['wb5'])
    xac = _div(-body_cma[0] * cbarr, body_cla[0]) + dxcg
    anum = xac * body_cla[0] / cbarr + c6 * wb4 * a10 / cbarr + wb13 * wb5
    wb12 = _div(anum, body_cla[0] + wb4 + wb5)

    cm_surface = np.asarray(surface_alone['cm'], dtype=float)
    cn_surface = np.asarray(surface_alone['cn'], dtype=float)
    ca_surface = np.asarray(surface_alone['ca'], dtype=float)
    cla_surface = float(surface_alone['cla'])
    dxcpbw = float(surface['a173']) / cbarr - wb13
    lever = (float(synthesis['zw']) - float(synthesis['zcg'])) / cbarr
    kwb, kbw = float(lift['kwb']), float(lift['kbw'])
    kwb_inc = float(lift['kwb_incidence'])
    kbw_inc = float(lift['kbw_incidence'])
    cm = []
    for angle_slot in range(len(alpha_deg)):
        dxcp_wing = (0.0 if cn_surface[angle_slot] == 0.0 else
                     ((cm_surface[angle_slot] - cm0_surface) /
                      cn_surface[angle_slot]))
        vortex_cn = (vortex['ivbw'][angle_slot] *
                     vortex['go2pav'][angle_slot] *
                     vortex['ratio'] * local_alpha_deg[angle_slot] *
                     cla_surface)
        cn_basic = cn_surface[angle_slot] - cla_surface * incidence
        value = (body['cm'][angle_slot] + cm0_surface +
                 cn_basic * kwb * dxcp_wing +
                 cla_surface * incidence * kwb_inc * dxcp_wing +
                 cn_basic * kbw * dxcpbw +
                 cla_surface * incidence * kbw_inc * dxcpbw +
                 vortex_cn * dxcp_wing + ca_surface[angle_slot] * lever)
        cm.append(NOT_AVAILABLE if cm_surface[angle_slot] == NOT_AVAILABLE
                  else value)
    return {
        'cm': np.array(cm), 'cm0': float(cmowb),
        'cm0_regression': regression is not None,
        'alpha_zero_lift': float(alpha_zero),
        'wb12': float(wb12), 'wb13': float(wb13), 'wb14': float(wb14),
        'wb15': float(wb15), 'method': 'legacy_wbcm',
    }


def calculate_wbaero(alpha_deg: Sequence[float],
                     surface: Dict[str, float],
                     surface_alone: Dict[str, object],
                     body: Dict[str, object],
                     synthesis: Dict[str, float],
                     flight: Dict[str, float],
                     cbarr: float,
                     experimental: Optional[Dict[str, object]] = None,
                     stale: Optional[Dict[str, float]] = None,
                     moment_cutoff: bool = True
                     ) -> Dict[str, object]:
    """Translate WBAERO's surface-body pass: the complete buildup.

    The same routines run for the wing (``BW``) and, with the tail's own
    blocks, for the horizontal tail (``BH``); see
    :func:`tail_body_inputs` for how the tail pass differs.

    Args:
        alpha_deg: ``FLC(23)`` onward.
        surface: Everything :func:`calculate_wbcm` reads from ``surface``,
            plus ``sspne`` (``WINGIN(3)``), ``a129``, ``a160``, ``clmax``
            and ``alpha_clmax`` (``B(44)``, ``B(43)``), ``local_alpha``
            (``B(23)`` onward), ``x_quarter_chord`` (``BD(83)``), ``cd0``
            and ``cdl`` (``D(20)``, ``D(36)`` onward) and ``c6``.
        surface_alone: ``cd``, ``cl``, ``cm``, ``cn``, ``ca`` (``WING(1)``,
            ``(21)``, ``(41)``, ``(61)``, ``(81)`` onward), ``cla`` and
            ``cma`` (``WING(101)``, ``WING(121)``).
        body: ``cd``, ``cl``, ``cm``, ``cla``, ``cma`` (``BODY(1)``,
            ``(21)``, ``(41)``, ``(101)``, ``(121)`` onward), ``cm0``
            ``BD(62)``, ``alpha_zero_lift`` ``BD(81)``, ``cd_friction``
            ``BD(59)``, ``cd_base`` ``BD(60)``, ``cdl`` ``BD(215)``
            onward, and the station arrays ``x`` and ``s``, from which
            GETMAX gives ``BD(3)`` and ``BD(1)`` is the last station.
        synthesis: ``xcg``, ``xw``, ``zw``, ``zcg``, ``aliw`` (``BD(77)``).
        flight: ``mach``, ``reynolds_per_length`` and ``tr``.
        cbarr: ``CBARR``.
        experimental: See :func:`calculate_wbdrag`.
        stale: See :func:`calculate_wblift`.
        moment_cutoff: The wing pass stops the CMa slope at the first
            unavailable moment and marks later ones unavailable; the tail
            pass has no such cutoff and takes the slope over every angle.

    Returns:
        Dictionary with ``cd``, ``cl``, ``cm``, ``cn``, ``ca``, ``cla`` and
        ``cma`` (``BW(1)``, ``(21)``, ``(41)``, ``(61)``, ``(81)``,
        ``(101)``, ``(121)`` onward), and the component results under
        ``drag``, ``lift`` and ``moment``.

    Notes:
        ``CMA`` is the TBFUNX slope over the angles before the first
        unavailable moment only, and is itself marked unavailable past
        them.  ``CN`` and ``CA`` here replace the ones WBLIFT and WBCM form
        with the local angle; theirs are never read.
    """
    alpha = np.asarray(alpha_deg, dtype=float)
    local = np.asarray(surface['local_alpha'], dtype=float)
    sspn, sspne = float(surface['sspn']), float(surface['sspne'])
    body_diameter = 2.0 * (sspn - sspne)
    incidence = float(synthesis['aliw'])

    # BODOWG, on the body's own angles BD(255) onward.
    _, max_area, _ = getmax(body['x'], body['s'])
    radius = np.sqrt(max_area / PI)
    body_alpha = alpha + float(body['alpha_zero_lift'])
    bodowg = [
        calculate_bodowg(angle_deg, float(surface['x_quarter_chord']),
                         radius, sspn, float(surface['a27']))
        for angle_deg in body_alpha
    ]
    vortex = {'ratio': (sspn - sspne) / sspn,
              'ivbw': [b['ivbw'] for b in bodowg],
              'go2pav': [b['go2pav'] for b in bodowg]}

    drag = calculate_wbdrag(
        float(flight['mach']), float(surface['a129']), float(body['x'][-1]),
        float(surface['cd0']), float(body['cd_friction']),
        float(body['cd_base']), body['cdl'], surface['cdl'],
        dict(experimental, body_cd=body['cd'], wing_cd=surface_alone['cd'])
        if experimental else None)
    lift = calculate_wblift(
        alpha, local, body_diameter, incidence, sspn,
        float(surface['alpha_zero_lift']), surface_alone, body['cl'],
        float(body['cla'][0]), vortex,
        {'a160': float(surface['a160']), 'clmax': float(surface['clmax']),
         'alpha_clmax': float(surface['alpha_clmax'])}, stale)
    moment = calculate_wbcm(
        alpha, local, surface, surface_alone,
        dict(body, length=float(body['x'][-1]), max_area=max_area),
        dict(synthesis, incidence=incidence,
             cos_incidence=np.cos(incidence / RAD),
             body_diameter=body_diameter),
        lift, vortex, flight, cbarr, float(surface['c6']))

    cd, cl, cm = drag['cd'], lift['cl'], moment['cm']
    available = len(alpha)
    for angle_index, value in enumerate(cm):
        if moment_cutoff and value == NOT_AVAILABLE:
            available = angle_index
            break
    cla = np.array([
        tbfunx(alpha, cl, angle_deg, 0, 0)[1] for angle_deg in alpha
    ])
    cma = np.array([
        tbfunx(alpha[:available], cm[:available], angle_deg, 0, 0)[1]
        if angle_index < available else NOT_AVAILABLE
        for angle_index, angle_deg in enumerate(alpha)
    ])
    cos_alpha = np.cos(alpha / RAD)
    sin_alpha = np.sin(alpha / RAD)
    return {
        'cd': cd, 'cl': cl, 'cm': cm,
        'cn': cl * cos_alpha + cd * sin_alpha,
        'ca': cd * cos_alpha - cl * sin_alpha,
        'cla': cla, 'cma': cma,
        'drag': drag, 'lift': lift, 'moment': moment, 'vortex': vortex,
        'method': 'legacy_wbaero',
    }


def tail_body_inputs(tail: Dict[str, float], tail_alone: Dict[str, object],
                     synthesis: Dict[str, float], bd63: float,
                     cbarr: float) -> Dict[str, object]:
    """The arguments WBAERO's horizontal-tail pass hands the shared routines.

    The tail pass differs from the wing pass only in what it passes: the
    tail's ``AHT``/``BHT``/``HTIN`` blocks, ``XH``, ``ZH`` and ``ALIH`` for
    ``XW``, ``ZW`` and ``ALIW``, a quarter-chord station ``AHT(161)+XH``,
    and for WBCM's ``CMA`` the product ``(BD(63)/CBARR)*HT(101)`` in place
    of the surface's own moment slope.

    Args:
        tail: The tail's ``surface`` entries, as for the wing pass, with
            ``a161`` for the quarter-chord station.
        tail_alone: The tail-alone curves; ``cla`` is ``HT(101)``.
        synthesis: ``xcg``, ``xh``, ``zh``, ``zcg``, ``alih``.
        bd63: ``BD(63)``, the tail moment arm.  See Notes.
        cbarr: ``CBARR``.

    Returns:
        ``{'surface', 'surface_alone', 'synthesis'}`` for
        :func:`calculate_wbaero`, called with ``moment_cutoff=False``.

    Notes:
        ``BD(63)`` is written by WBTAIL, in overlay M10O12, which the main
        program runs *after* WBAERO's overlay M07O07 within each Mach
        iteration.  The tail pass therefore reads the previous Mach's arm,
        and zero on the first.  It enters only the tail-body zero-lift
        moment ``HB(16)``, not the moment curve.  The caller supplies the
        value the source would read.
    """
    surface = dict(
        tail,
        x_quarter_chord=float(tail['a161']) + float(synthesis['xh']),
    )
    tail_cla = float(tail_alone['cla'])
    alone = dict(tail_alone, cma=(bd63 / cbarr) * tail_cla)
    synthesis_for_tail = {
        'xcg': synthesis['xcg'],
        'xw': synthesis['xh'],
        'zw': synthesis['zh'],
        'zcg': synthesis['zcg'],
        'aliw': synthesis['alih'],
    }
    return {
        'surface': surface,
        'surface_alone': alone,
        'synthesis': synthesis_for_tail,
    }


def calculate_body_vertical(alpha_deg: Sequence[float],
                            body: Dict[str, Sequence[float]],
                            vertical_cd0: float) -> Dict[str, np.ndarray]:
    """Translate WBAERO's body-vertical pass, the BV set.

    The body's curves with the vertical tail and ventral fin zero-lift drag
    ``DVT(20)+DVF(20)`` added to its drag; the lift and moment slopes are
    retaken by TBFUNX from the second angle, the first keeping the body's.

    Args:
        alpha_deg: ``FLC(23)`` onward.
        body: ``cd``, ``cl``, ``cm``, ``cla``, ``cma`` (``BODY(1)``, ``(21)``,
            ``(41)``, ``(101)``, ``(121)`` onward).
        vertical_cd0: ``DVT(20)+DVF(20)``.
    """
    alpha = np.asarray(alpha_deg, dtype=float)
    out = {k: np.array(body[k], dtype=float)
           for k in ('cd', 'cl', 'cm', 'cla', 'cma')}
    out['cd'] = out['cd'] + vertical_cd0
    cos_alpha = np.cos(alpha / RAD)
    sin_alpha = np.sin(alpha / RAD)
    out['cn'] = out['cl'] * cos_alpha + out['cd'] * sin_alpha
    out['ca'] = out['cd'] * cos_alpha - out['cl'] * sin_alpha
    for angle_index in range(1, len(alpha)):
        out['cla'][angle_index] = tbfunx(alpha, out['cl'],
                                         alpha[angle_index], 0, 0)[1]
        out['cma'][angle_index] = tbfunx(alpha, out['cm'],
                                         alpha[angle_index], 0, 0)[1]
    return out


# TABLES: the Mach grid of the WBCDL regression, and the angle past which
# each speed range has no data.
_TABLES_MACH = np.array([0.00, 0.60, 0.70, 0.80, 0.90, 0.95, 1.00, 1.10,
                         1.20, 1.30, 1.40, 1.50, 2.00, 2.50])


def _angle_limit(mach: float) -> int:
    if mach <= 0.9:
        return 18
    if mach < 1.0:
        return 11
    if mach < 1.1:
        return 12
    return 15


def calculate_tables(mach: float, alpha_deg: float) -> Optional[np.ndarray]:
    """Translate TABLES: the sixteen WBCDL coefficients at one condition.

    Bilinear in angle (one-degree rows from TBSUB, TBTRN or TBSUP by Mach
    index) and Mach.  Returns ``None`` where the source sets ``NDM``: past
    the speed range's angle limit (18, 11, 12 or 15 degrees) or outside
    the Mach grid.  The source then leaves its coefficient array as it was;
    WBCDL overwrites the result it forms from them, so nothing stale
    reaches the output.
    """
    alp = abs(float(alpha_deg))
    iam = _angle_limit(mach)
    if alp > iam:
        return None
    ia = int(alp)
    if ia == iam:
        ia = iam - 1
    base = float(ia)
    mach_index_high = None
    for mach_bracket in range(2, 15):
        if (_TABLES_MACH[mach_bracket - 2] <= mach <=
                _TABLES_MACH[mach_bracket - 1]):
            mach_index_high = mach_bracket
            break
    if mach_index_high is None:
        return None
    mach_index_low = mach_index_high - 1
    coeff_grid = np.zeros((2, 2, 16))
    for mach_row, mach_index in enumerate((mach_index_low, mach_index_high)):
        for angle_col, angle_index in enumerate((ia + 1, ia + 2)):
            if mach_index <= 4:
                block = calculate_tbsub(mach_index, angle_index)
            elif mach_index <= 10:
                block = calculate_tbtrn(mach_index, angle_index)
            else:
                block = calculate_tbsup(mach_index, angle_index)
            coeff_grid[mach_row, angle_col] = block['coefficients']
    at_low_mach = (coeff_grid[0, 0] +
                   (coeff_grid[0, 1] - coeff_grid[0, 0]) * (alp - base))
    at_high_mach = (coeff_grid[1, 0] +
                    (coeff_grid[1, 1] - coeff_grid[1, 0]) * (alp - base))
    mach_low = _TABLES_MACH[mach_index_low - 1]
    mach_high = _TABLES_MACH[mach_index_high - 1]
    return at_low_mach + (at_high_mach - at_low_mach) * (
        (mach - mach_low) / (mach_high - mach_low))


def calculate_wbcdl(aspect_ratio: float, tan_le: float, tovc: float,
                    nose_length: float, afterbody_length: float,
                    taper_ratio: float, leading_edge_radius: float,
                    twist_deg: float, ycm: float, cld: float,
                    reynolds: float, tr: float, mach: float,
                    alpha_deg: Sequence[float]) -> Optional[np.ndarray]:
    """Translate WBCDL: the regression drag due to lift of a wing-body.

    Returns ``None`` (the source's ``NA``) outside the regression's range,
    else ``CDL`` at each angle, ``UNUSED`` where TABLES has no data.  The
    source also sets ``NA`` when the first angle has no data; its caller
    then discards the whole curve, which ``calculate_wbcd`` reproduces.
    The Reynolds number is clamped to 8e5 to 8e6.
    """
    inside = (1.6 <= aspect_ratio <= 6.0 and 0.0 <= tan_le <= 2.74748 and
              0.025 <= tovc <= 0.100 and 2.2 <= nose_length <= 8.4 and
              0.3 <= afterbody_length <= 5.6 and 0.0 <= taper_ratio <= 1.0
              and 0.0 <= leading_edge_radius <= 0.015 and
              -9.4 <= twist_deg <= UNUSED and 0.0 <= ycm <= 0.0263 and
              0.0 <= cld <= 0.45)
    if not inside:
        return None
    reynolds = min(max(reynolds, 8.0e5), 8.0e6)
    cdl = []
    for angle_deg in alpha_deg:
        b = calculate_tables(mach, angle_deg)
        if b is None:
            cdl.append(UNUSED)
            continue
        cdl.append(b[0] + b[1] / aspect_ratio + b[2] * aspect_ratio +
                   b[3] * np.sqrt(tan_le) + b[4] * tovc + b[5] * nose_length +
                   b[6] * afterbody_length + b[7] * taper_ratio +
                   b[8] * taper_ratio**2 + b[9] * taper_ratio**3 +
                   b[10] * tr + b[11] * leading_edge_radius -
                   b[12] * twist_deg / RAD + b[13] * ycm + b[14] * cld +
                   b[15] * reynolds / 1.0e6)
    return np.array(cdl, dtype=float)


def surface_regression_drag(alpha_deg: Sequence[float],
                            surface: Dict[str, float],
                            body: Dict[str, float], x_surface: float,
                            flight: Dict[str, float]):
    """WBCDL as WBCD and TRANCD call it for one surface.

    Forms the nose and afterbody lengths from the surface's position and
    calls WBCDL.  Returns ``(cdl, ln, la)``, or ``None`` when the surface is
    not straight tapered, the regression is out of range, or its first
    angle has no data (the source's ``NA``).  Arguments are those of
    :func:`calculate_wbcd`.
    """
    if float(surface['type']) != STRAIGHT_TAPERED_WB:
        return None
    db = 2.0 * (float(surface['sspn']) - float(surface['sspne']))
    ln = (x_surface + 0.5 * db * float(surface['a38'])) / db
    la = (x_surface + float(surface['chrdr']) +
          0.5 * db * float(surface['a56'])) / db
    if float(body['x_max_area']) < ln * db:
        ln = ln * db / float(body['max_diameter'])
    if float(body['x_max_area']) > la * db:
        la = la * db / float(body['max_diameter'])
    la = float(body['length']) / db - la
    alpha = np.asarray(alpha_deg, dtype=float)
    cdl = calculate_wbcdl(
        float(surface['a120']), float(surface['a38']), float(surface['tovc']),
        ln, la, float(surface['a118']), float(surface['ler']),
        float(surface['twista']), float(surface['ycm']),
        float(surface['cld']),
        float(flight['reynolds_per_length']) * float(surface['a122']),
        float(flight['tr']), float(flight['mach']), alpha)
    if cdl is None or cdl[0] == UNUSED:
        return None
    return cdl, ln, la


def calculate_wbcd(alpha_deg: Sequence[float], surface: Dict[str, float],
                   combination: Dict[str, Sequence[float]],
                   body: Dict[str, float], x_surface: float,
                   flight: Dict[str, float]) -> Optional[Dict[str, object]]:
    """Translate one half of WBCD: the regression drag of a surface-body.

    The wing and horizontal-tail halves are the same code on their own
    blocks.  Only a straight tapered surface is treated; otherwise, or when
    WBCDL is out of range, the source leaves the buildup drag and this
    returns ``None``.

    Args:
        alpha_deg: ``FLC(23)`` onward.
        surface: ``type`` (``WINGIN(15)``), ``sspn``, ``sspne``, ``chrdr``
            (4, 3, 6), ``tovc`` (16), ``ler`` (62), ``twista`` (11),
            ``ycm`` (93), ``cld`` (94), and ``a38``, ``a56``, ``a118``,
            ``a120``, ``a122``.
        combination: The surface-body ``cd0`` (``WB(17)`` or ``HB(17)``) and
            ``cl`` (``BW(21)`` or ``BH(21)`` onward).
        body: ``length`` ``BD(1)``, ``x_max_area`` ``BD(2)`` and
            ``max_diameter`` ``BD(85)``.
        x_surface: ``XW`` or ``XH``.
        flight: ``mach`` (``FLC(I+2)``), ``reynolds_per_length``
            (``FLC(I+42)``) and ``tr`` (``FLC(96)``).

    Returns:
        ``{'cd', 'cn', 'ca', 'cdl'}``, or ``None`` when not applicable.
        Angles without regression data are marked ``-UNUSED`` in ``cd``,
        ``cn`` and ``ca``, as the source marks them.

    Notes:
        When the body's maximum-area station lies ahead of the surface's
        leading or trailing edge, the source renormalises the nose or
        afterbody station by the maximum diameter ``BD(85)`` instead of
        ``DB``, and then forms the afterbody length as ``BD(1)/DB - LA``,
        mixing the two scales in one subtraction.  Kept.
    """
    regression = surface_regression_drag(alpha_deg, surface, body, x_surface,
                                         flight)
    if regression is None:
        return None
    cdl, ln, la = regression
    alpha = np.asarray(alpha_deg, dtype=float)
    cl = np.asarray(combination['cl'], dtype=float)
    cd = float(combination['cd0']) + cdl
    cos_alpha = np.cos(alpha / RAD)
    sin_alpha = np.sin(alpha / RAD)
    cn = cl * cos_alpha + cd * sin_alpha
    ca = cd * cos_alpha - cl * sin_alpha
    missing = cdl == UNUSED
    cd, cn, ca = (np.where(missing, -UNUSED, v) for v in (cd, cn, ca))
    return {'cd': cd, 'cn': cn, 'ca': ca, 'cdl': cdl,
            'nose_length': float(ln), 'afterbody_length': float(la),
            'method': 'legacy_wbcd'}
