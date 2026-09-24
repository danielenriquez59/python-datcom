"""
Subsonic wing-body-tail buildup: overlay M10O12 with WGEOTL and WBTAIL.

``M10O12`` runs after the wing-body pass (WBAERO) and the downwash (DWASH):

- ``WGEOTL``: the Figure 4.4.1-71 vortex position and the ALI interference
  factor at the tail, which the canard method uses; for a canard it also
  replaces DWASH's vortex span.
- ``WBTAIL``: the tail's lift, moment and drag added to the wing-body at
  each angle, through the tail carryover factors K_H(B) and K_B(H), the
  downwash and dynamic-pressure ratio, and the body vortex on the tail.
- ``M10O12`` itself: normal and axial force, the lift- and moment-curve
  slopes for the second angle onward, and the body-wing-vertical and
  body-wing-tail-vertical sets.

``TWASH`` above 2.5 (``WINGIN(101)``) selects the canard method, in which
the forward surface's downwash does not act on the wing and a vortex-lift
increment is added instead.

The tables were extracted from the source DATA statements by parsing, with
``tools/fortran_data.py``, rather than by hand.  WBTAIL's five carryover
tables are the same digitisations as WBLIFT's, and are shared.

Reference: datcom-legacy/datcom_2000/m10o12.f, wgeotl.f, wbtail.f
"""

import math
import numpy as np
from typing import Dict, Optional, Sequence
import logging

from pydatcom.aerodynamics.wing_body import (
    NOT_AVAILABLE, _X10A, _Y10A, _X10B, _Y10B, _X12A1, _Y12A1, _X12A2,
    _Y12A2,
)
from pydatcom.interactions.body_vortex import ali, calculate_bodowg, getmax
from pydatcom.utils.constants import PI, RAD, UNUSED
from pydatcom.utils.legacy_numeric import tbfunx
from pydatcom.utils.legacy_tables import tlin3x

logger = logging.getLogger(__name__)

# WGEOTL: Figure 4.4.1-71A-C, the vortex lateral position.  Y44160 is
# (9,5,3): beta*A fastest, then A*tan(LE sweep), then taper ratio.
_X41602 = np.array([
    0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0,
])
_X41601 = np.array([
    -4.0, -2.0, 0.0, 2.0, 4.0,
])
_X41603 = np.array([
    0.0, 0.5, 1.0,
])
_Y44160 = np.array([
    0.785, 0.6014, 0.601, 0.5992, 0.5964, 0.5932, 0.5898, 0.5864, 0.583,
    0.785, 0.6735, 0.6633, 0.6526, 0.6428, 0.634, 0.6263, 0.6194, 0.6133,
    0.785, 0.7532, 0.7254, 0.7039, 0.6868, 0.6728, 0.6612, 0.6513, 0.6427,
    0.785, 0.7732, 0.754, 0.7348, 0.7176, 0.7027, 0.6898, 0.6786, 0.6687,
    0.785, 0.7722, 0.7613, 0.7485, 0.7353, 0.7226, 0.7109, 0.7002, 0.6904,
    0.785, 0.661, 0.6656, 0.6702, 0.6747, 0.679, 0.6831, 0.6868, 0.6903,
    0.785, 0.7326, 0.7264, 0.7241, 0.724, 0.7251, 0.7267, 0.7283, 0.7299,
    0.785, 0.7795, 0.7754, 0.7739, 0.7734, 0.7732, 0.7731, 0.773, 0.7728,
    0.785, 0.8137, 0.8204, 0.8234, 0.824, 0.8232, 0.8217, 0.8199, 0.818,
    0.785, 0.8922, 0.8876, 0.884, 0.8803, 0.8764, 0.8724, 0.8683, 0.8643,
    0.785, 0.6787, 0.6882, 0.6992, 0.7106, 0.7222, 0.7334, 0.7441, 0.7542,
    0.785, 0.747, 0.7464, 0.7524, 0.7615, 0.7718, 0.7822, 0.7923, 0.8017,
    0.785, 0.7876, 0.7965, 0.8078, 0.8194, 0.8304, 0.8404, 0.8496, 0.8579,
    0.785, 0.8495, 0.8665, 0.8809, 0.8926, 0.9022, 0.9102, 0.917, 0.9229,
    0.785, 0.9866, 0.984, 0.9848, 0.9869, 0.9894, 0.9917, 0.9938, 0.9956,
])

_FIG_71 = _Y44160.reshape(3, 5, 9).transpose(2, 1, 0)

# WBTAIL's tail geometry takes the tail incidence as zero.
_SINAI, _COSAI, _TANAI = 0.0, 1.0, 0.0


def _f(value) -> np.float64:
    """A float64, so a zero dynamic-pressure ratio divides as it does in
    the source rather than raising."""
    return np.float64(value)


def calculate_wgeotl(alpha_deg: Sequence[float], wing: Dict[str, float],
                     tail: Dict[str, float], aliw: float) -> Dict[str, object]:
    """Translate WGEOTL: the wing vortex position and its effect on the tail.

    Figure 4.4.1-71 gives the vortex semispan fraction ``YT``; the vortex
    sits at ``BVTO2 = SSPNE*YT + (SSPN - SSPNE)``, and ALI evaluates its
    interference on the tail at each angle.

    Args:
        alpha_deg: ``FLC(23)`` onward.
        wing: ``a120`` (theoretical aspect ratio), ``a38`` (tan LE),
            ``a118`` (theoretical taper), ``a12`` and ``a24`` (INFTGM's
            tail height and arm), ``a80`` (tan TE), ``beta`` ``B(2)``,
            ``sspn`` and ``sspne``, and ``twash`` ``WINGIN(101)``.
        tail: ``sspn`` and ``sspne`` (``HTIN(4)``, ``HTIN(3)``) and ``a27``.
        aliw: Wing incidence, degrees.

    Returns:
        Dictionary with ``ali`` (``FACT(42)`` onward), ``yt``, ``bvto2``,
        and ``vortex_span`` (``2*BVTO2``, which replaces ``FACT(82)``
        onward for a canard) or ``None``.
    """
    ata = float(wing['a120']) * float(wing['a38'])
    ba = float(wing['a120']) * float(wing['beta'])
    yt = tlin3x(_X41601, _X41602, _X41603, _FIG_71, ata, ba,
                float(wing['a118']), 0, 0, 2, 0, 1, 2)
    sspn, sspne = float(wing['sspn']), float(wing['sspne'])
    bvto2 = sspne * yt + (sspn - sspne)
    radius = float(tail['sspn']) - float(tail['sspne'])
    factors = []
    for a in alpha_deg:
        z0 = (float(wing['a12']) - (float(wing['a24']) - float(wing['a80']) *
                                    bvto2) * math.tan((a + aliw) / RAD))
        factors.append(ali(z0, bvto2, float(tail['sspn']), radius,
                           float(tail['a27'])))
    canard = float(wing.get('twash', 0.0)) > 2.5
    return {
        'ali': np.array(factors), 'yt': float(yt), 'bvto2': float(bvto2),
        'vortex_span': 2.0 * bvto2 if canard else None,
        'method': 'legacy_wgeotl',
    }


def calculate_wbtail(alpha_deg: Sequence[float],
                     tail: Dict[str, object],
                     tail_alone: Dict[str, Sequence[float]],
                     wing_body: Dict[str, object],
                     downwash: Dict[str, Sequence[float]],
                     wing: Dict[str, object],
                     synthesis: Dict[str, float],
                     body: Dict[str, object],
                     sref: float, cbarr: float,
                     vertical_cd0: float = 0.0) -> Dict[str, object]:
    """Translate WBTAIL: the wing-body-tail buildup at each angle.

    Args:
        alpha_deg: ``FLC(23)`` onward.
        tail: ``/HTI/`` ``sspn`` and ``sspne``; ``AHT`` entries ``a3``,
            ``a7``, ``a10``, ``a27``, ``a62``, ``a161``; ``c6`` (``CHT(6)``);
            ``BHT`` entries ``local_alpha`` (23 onward), ``cd0`` (46),
            ``cm0`` (47) and ``alpha_zero_lift`` (49).
        tail_alone: ``cd``, ``cl`` and ``cla``: ``HT(1)``, ``(21)`` and
            ``(101)`` onward.
        wing_body: ``cd``, ``cl``, ``cm`` (``BW(1)``, ``(21)``, ``(41)``
            onward), ``cla`` and ``cma`` (``BW(101)``, ``BW(121)``),
            ``cd0`` (``WB(17)``) and ``kwb`` (``WB(2)``).
        downwash: ``qoqi``, ``angle``, ``gradient`` (``DWASH(1)``, ``(21)``,
            ``(41)`` onward); for a canard also ``ali`` and
            ``vortex_span`` (``FACT(42)`` and ``FACT(82)`` onward).
        wing: ``sspn`` and ``sspne`` (``WINGIN(4)``, ``(3)``), ``twash``
            (``WINGIN(101)``); for a canard ``cla`` (``WING(101)``),
            ``local_alpha`` (``B(23)`` onward) and ``alpha_zero_lift``.
        synthesis: ``xcg``, ``xh``, ``zh``, ``zcg``, ``alih``.
        body: ``x`` and ``s`` (for BODOWG's GETMAX) and ``bd70``.
        sref, cbarr: ``SREF``, ``CBARR``.
        vertical_cd0: ``DVT(20)+DVF(20)``.

    Returns:
        Dictionary with the body-wing-tail ``cd``, ``cl``, ``cm``, ``cla``
        and ``cma`` (``BWH(1)``, ``(21)``, ``(41)``, ``(101)``, ``(121)``
        onward), ``cd_with_vertical`` (``BWHV(1)`` onward), the tail
        geometry written to ``BD``, the ``WBT`` quantities, and for a canard
        ``FACT(102)`` and ``FACT(122)`` onward.

    Notes:
        WBTAIL hands ``WBT(108) = SSPN-SSPNE`` to BODOWG as the argument
        BODOWG overwrites with the maximum body radius ``sqrt(BD(3)/pi)``,
        and only reads ``WBT(108)`` afterwards.  So the tail's body-vortex
        lift scales with the maximum body radius, not with the body radius
        at the tail as the wing's does.  This is what executes, and it is
        kept; ``body_radius`` is returned.

        On the canard path the first-loop lookups of the tail's CL, CLa and
        CD use ``BHT(J+22)`` on a grid of ``BHT(J+22)-BHT(49)``, without the
        zero-lift shift the other path applies, so the canard is read at
        its angle plus its zero-lift angle.  Kept.
    """
    canard = float(wing.get('twash', 0.0)) > 2.5
    alpha = np.asarray(alpha_deg, dtype=float)
    n = len(alpha)
    sspn, sspne = float(tail['sspn']), float(tail['sspne'])
    a62, a10 = float(tail['a62']), float(tail['a10'])
    xcg, xh = float(synthesis['xcg']), float(synthesis['xh'])
    zh, zcg = float(synthesis['zh']), float(synthesis['zcg'])

    bd58 = (sspn - sspne) * a62 * _COSAI
    bd84 = xcg - (xh + bd58)
    bd762 = zh - bd58 * _TANAI
    bd761 = float(tail['c6']) * a10
    bd64 = bd762 - bd761 * _SINAI - zcg
    bd63 = xcg - (bd761 + bd58 + xh)
    bd31 = (bd63 + bd64 * _TANAI) * _COSAI
    bd30 = bd64 / _COSAI - (bd63 + bd64 * _TANAI) * _SINAI
    bd8 = xh - float(tail['a161']) * _COSAI

    local = np.asarray(tail['local_alpha'], dtype=float)
    alpha_zero = float(tail['alpha_zero_lift'])
    ang = local - alpha_zero
    q = np.asarray(downwash['qoqi'], dtype=float)
    eps = np.asarray(downwash['angle'], dtype=float)
    deda = np.asarray(downwash['gradient'], dtype=float)
    ht_cl = np.asarray(tail_alone['cl'], dtype=float)
    ht_cd = np.asarray(tail_alone['cd'], dtype=float)
    ht_cla = float(np.asarray(tail_alone['cla'], dtype=float)[0])
    bw_cd = np.asarray(wing_body['cd'], dtype=float)
    bw_cl = np.asarray(wing_body['cl'], dtype=float)
    bw_cm = np.asarray(wing_body['cm'], dtype=float)
    bw_cla, bw_cma = float(wing_body['cla']), float(wing_body['cma'])

    clh, cdh, bd94 = [], [], []
    for j in range(n):
        alpat = local[j] if canard else local[j] - eps[j] - alpha_zero
        clh.append(tbfunx(ang, ht_cl, alpat, 1, 1)[0])
        cdh.append(tbfunx(ang, ht_cd, alpat, 1, 1)[0])
        bd94.append(tbfunx(alpha, bw_cd, alpha[j], 0, 0)[1])

    ratio = (sspn - sspne) / sspn
    akhbi, _ = tbfunx(_X12A1, _Y12A1, ratio, 0, 0)
    akbhi, _ = tbfunx(_X12A2, _Y12A2, ratio, 0, 0)
    wbt1, _ = tbfunx(_X10A, _Y10A, ratio, 0, 0)
    wbt2, _ = tbfunx(_X10B, _Y10B, ratio, 0, 0)
    wbt3, wbt4 = wbt1 * ht_cla, wbt2 * ht_cla
    wbt109 = xh + float(tail['a161'])

    # BODOWG, on the free-stream angles; it leaves the maximum body radius
    # in WBT(108).
    _, max_area, _ = getmax(body['x'], body['s'])
    radius = math.sqrt(max_area / PI)
    vortex = [calculate_bodowg(a, wbt109, radius, sspn, float(tail['a27']))
              for a in alpha]

    incidence_shift = (akhbi + akbhi) / (wbt1 + wbt2) * float(synthesis['alih'])
    bd70 = float(body['bd70'])
    out = {k: np.zeros(n) for k in ('cd', 'cd_with_vertical', 'cl', 'cm',
                                    'cla', 'cma', 'tail_lift',
                                    'vortex_lift', 'wbt87')}
    fact101, fact121, wbt25s = np.zeros(n), np.zeros(n), np.zeros(n)
    with np.errstate(divide='ignore', invalid='ignore'):
        for j in range(n):
            qj = _f(q[j])
            alpat = local[j] if canard else local[j] - eps[j]
            cla = bw_cla + (wbt3 + wbt4) * (1.0 - deda[j]) * qj
            wbt25 = 0.0
            if canard:
                anum = (float(wing['cla']) * ht_cla * qj *
                        float(wing_body['kwb']) * downwash['ali'][j] * sspne)
                aden = 2.0 * PI * float(tail['a7']) * (
                    downwash['vortex_span'][j] / 2.0 - float(wing['sspn']) +
                    float(wing['sspne']))
                wbt25 = anum / aden * RAD * sref / float(tail['a3'])
                cla = bw_cla + (wbt3 + wbt4) * qj + wbt25
                fact121[j] = -wbt25 / ((wbt3 + wbt4) * qj)
                wbt25s[j] = wbt25
            out['cla'][j] = cla

            alpef = alpha[j] - alpha_zero + incidence_shift
            alpaht = local[j] - alpha_zero
            if not canard:
                alpef -= eps[j]
                alpaht -= eps[j]
            vortex_lift = (vortex[j]['ivbw'] * vortex[j]['go2pav'] * radius *
                           qj * ht_cla * alpaht / sspn)
            alpa = alpha[j] - eps[j]
            tail_lift = (tbfunx(ang, ht_cl, alpef, 1, 1)[0] * qj *
                         (wbt1 + wbt2))
            if canard:
                cl = bw_cl[j] + tail_lift + vortex_lift + wbt25 * alpha[j]
                if alpaht != 0.0:
                    fact101[j] = ((-wbt25 * (float(wing['local_alpha'][j]) -
                                             float(wing['alpha_zero_lift']))) /
                                  (ht_cla * (wbt1 + wbt2) * qj))
                alpa = alpha[j]
            else:
                cl = bw_cl[j] + tail_lift + vortex_lift
            out['cl'][j], out['tail_lift'][j] = cl, tail_lift
            out['vortex_lift'][j] = vortex_lift

            sa, ca = math.sin(alpa / RAD), math.cos(alpa / RAD)
            sf, cf = math.sin(alpha[j] / RAD), math.cos(alpha[j] / RAD)
            dxacwb = bw_cma / bw_cla
            apart = dxacwb * ((-bw_cl[j] / RAD + bd94[j]) * sf +
                              (bw_cla + bw_cd[j] / RAD) * cf)
            bpart = (bd70 / cbarr) * ((bw_cla + bw_cd[j] / RAD) * sf +
                                      (bw_cl[j] / RAD - bd94[j]) * cf)
            cdht, dcdda = tbfunx(local, ht_cd, alpat, 0, 0)
            clht = (cl - bw_cl[j]) / qj
            dclda = wbt3 + wbt4
            if canard:
                dclda = dclda + wbt25 / qj
            cpart = (-clht / RAD + dcdda) * sa
            dpart = (clht / RAD - dcdda) * ca
            epart = qj * (1.0 - deda[j])
            fpart = (dclda + cdht / RAD) * sa
            gpart = (dclda + cdht / RAD) * ca
            wbt87 = ((bd63 / cbarr) * (cpart + gpart) * epart -
                     (bd64 / cbarr) * (fpart + dpart) * epart)
            if canard:
                wbt87 = wbt87 / (1.0 - deda[j])
            out['wbt87'][j] = wbt87
            out['cma'][j] = (NOT_AVAILABLE if bw_cm[j] == NOT_AVAILABLE
                             else apart - bpart + wbt87)

            if canard:
                dclht = cl - bw_cl[j]
            else:
                dclht = (((cl - bw_cl[j]) / qj +
                          cdh[j] * math.sin(eps[j] / RAD)) /
                         math.cos(eps[j] / RAD)) * qj
            cdhq = cdh[j] * qj
            cm = (bw_cm[j] + (bd63 / cbarr) * (dclht * ca + cdhq * sa) +
                  (bd64 / cbarr) * (cdhq * ca - dclht * sa) +
                  float(tail['cm0']) * qj)
            out['cm'][j] = NOT_AVAILABLE if bw_cm[j] == NOT_AVAILABLE else cm

            downwash_angle = 0.0 if canard else eps[j]
            with_vertical = (bw_cd[j] + vertical_cd0 +
                             (cdh[j] * math.cos(downwash_angle / RAD) +
                              clh[j] * math.sin(downwash_angle / RAD)) * qj)
            out['cd_with_vertical'][j] = with_vertical
            out['cd'][j] = with_vertical - vertical_cd0

    out.update({
        'canard': canard,
        'geometry': {'bd8': bd8, 'bd30': bd30, 'bd31': bd31, 'bd58': bd58,
                     'bd63': bd63, 'bd64': bd64, 'bd84': bd84,
                     'bd761': bd761, 'bd762': bd762},
        'cd_slope': np.array(bd94),
        'khb': float(wbt1), 'kbh': float(wbt2),
        'khb_incidence': float(akhbi), 'kbh_incidence': float(akbhi),
        'wbt3': float(wbt3), 'wbt4': float(wbt4),
        'wbt66': float(vertical_cd0),
        'wbt67': float(wing_body['cd0']) + float(tail['cd0']) + vertical_cd0,
        'body_radius': radius,
        'ivbw': np.array([v['ivbw'] for v in vortex]),
        'go2pav': np.array([v['go2pav'] for v in vortex]),
        'method': 'legacy_wbtail',
    })
    if canard:
        out.update({'fact101': fact101, 'fact121': fact121,
                    'canard_lift_slope': wbt25s})
    return out


def calculate_m10o12(alpha_deg: Sequence[float],
                     wing_body: Dict[str, Sequence[float]],
                     vertical_cd: float = 0.0, ventral_cd: float = 0.0,
                     tail: Optional[Dict[str, object]] = None
                     ) -> Dict[str, Dict[str, np.ndarray]]:
    """Translate M10O12's closing pass: the BWH, BWV and BWHV sets.

    Args:
        alpha_deg: ``FLC(23)`` onward.
        wing_body: The BW set: ``cd``, ``cl``, ``cm``, ``cla``, ``cma``.
        vertical_cd, ventral_cd: ``VT(1)`` and ``VF(1)``, which the source
            adds to the drag at every angle.
        tail: :func:`calculate_wbtail`'s result, or ``None`` without a
            horizontal tail (``HTPL`` false).

    Returns:
        ``{'bwv': ..., 'bwh': ..., 'bwhv': ...}``, each with ``cd``, ``cl``,
        ``cm``, ``cn``, ``ca``, ``cla``, ``cma``; the last two only with a
        tail.

    Notes:
        The lift and moment slopes are recomputed by TBFUNX from the second
        angle on only; at the first angle WBTAIL's analytic ``BWH(101)``
        and ``BWH(121)`` stand.  The moment slope uses only the angles
        before the first unavailable wing-body moment.
    """
    alpha = np.asarray(alpha_deg, dtype=float)
    ca, sa = np.cos(alpha / RAD), np.sin(alpha / RAD)
    bw = {k: np.asarray(wing_body[k], dtype=float)
          for k in ('cd', 'cl', 'cm', 'cla', 'cma')}
    available = len(alpha)
    for j, value in enumerate(bw['cm']):
        if value == NOT_AVAILABLE:
            available = j
            break

    bwv = {'cd': bw['cd'] + vertical_cd + ventral_cd, 'cl': bw['cl'].copy(),
           'cm': bw['cm'].copy(), 'cla': bw['cla'].copy(),
           'cma': bw['cma'].copy()}
    bwv['cn'] = bwv['cl'] * ca + bwv['cd'] * sa
    bwv['ca'] = bwv['cd'] * ca - bwv['cl'] * sa
    result = {'bwv': bwv}

    if tail is not None:
        bwh = {k: np.array(tail[k], dtype=float)
               for k in ('cd', 'cl', 'cm', 'cla', 'cma')}
        bwh['cn'] = bwh['cl'] * ca + bwh['cd'] * sa
        bwh['ca'] = bwh['cd'] * ca - bwh['cl'] * sa
        for j in range(1, len(alpha)):
            bwh['cla'][j] = tbfunx(alpha, bwh['cl'], alpha[j], 0, 0)[1]
            bwh['cma'][j] = (
                tbfunx(alpha[:available], bwh['cm'][:available], alpha[j],
                       0, 0)[1] if j < available else NOT_AVAILABLE)
        bwhv = {'cd': np.array(tail['cd_with_vertical'], dtype=float),
                'cl': bwh['cl'].copy(), 'cm': bwh['cm'].copy(),
                'cla': bwh['cla'].copy(), 'cma': bwh['cma'].copy()}
        bwhv['cn'] = bwhv['cl'] * ca + bwhv['cd'] * sa
        bwhv['ca'] = bwhv['cd'] * ca - bwhv['cl'] * sa
        result.update({'bwh': bwh, 'bwhv': bwhv})

    missing = bw['cd'] == -UNUSED
    for key in result:
        for component in ('cd', 'cn', 'ca'):
            result[key][component][missing] = -UNUSED
    return result
