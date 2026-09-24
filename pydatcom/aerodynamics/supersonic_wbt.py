"""
Supersonic wing-body-tail: SUPWBT.

The supersonic counterpart of the subsonic WBTAIL buildup: the tail-body
carryover K_B(T) by INTKBW's integration or Figure 4.3.1.2-10 and K_T(B)
from Figure 4.3.1.2-10, the incidence carryover of Figure 4.3.1.2-12A/B,
the body-vortex lift on the tail (BODOWG), the canard's wing-vortex
interference on the aft surface (Figure 4.4.1-80 and ALI), and the
wing-body-tail lift, drag, force and moment curves with the tail seen
through the supersonic downwash (SDWASH's ``QOQINF``, ``EPSLON``,
``DEDALP`` and the local Mach number ``HMACH``).

The tables were extracted from the source DATA statements by parsing, with
``tools/fortran_data.py``.  Figure 4.3.1.2-10 and 4.3.1.2-12A are WBTRAN's
and SUPWB's, pinned equal by test; SUPWBT's Figure 4.3.1.2-12B keeps the
chart's end abscissae where SUPWB's moves them inward.

Reference: datcom-legacy/datcom_2000/supwbt.f
"""

import math
from typing import Dict, Mapping

import numpy as np

from pydatcom.aerodynamics.supersonic_wing_body import _D4312A, _T4312A
from pydatcom.aerodynamics.transonic_buildup import (_DKBW10, _DKWB10,
                                                     _TFIG10)
from pydatcom.interactions.body_vortex import ali, calculate_bodowg, getmax
from pydatcom.interactions.carryover import intkbw
from pydatcom.utils.constants import PI, RAD, UNUSED
from pydatcom.utils.legacy_interp import interx
from pydatcom.utils.legacy_numeric import tbfunx

# Figure 4.3.1.2-12B: K_W(B) at incidence.  SUPWBT keeps the chart's
# own end abscissae 0 and 1, where SUPWB/SUPHB use 0.015 and 0.975.
_T4312B = np.array([
    0.0, 0.1, 0.2, 0.3, 0.4, 0.6, 0.8, 1.0,
])
_D4312B = np.array([
    1.0, 0.975, 0.956, 0.947, 0.941, 0.95, 0.978, 1.0,
])
# Figure 4.4.1-80A-C: the lateral position of the wing tip vortex,
# y_v/(b'/2), against beta*A (9), A*tan(LE) (5, padded to 9) and taper (3).
_T4467 = np.array([
    0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0,
    -4.0, -2.0, 0.0, 2.0, 4.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.5, 1.0,
])
_D4467 = np.array([
    0.785, 0.7, 0.625, 0.565, 0.515, 0.49, 0.475, 0.47, 0.47,
    0.785, 0.7, 0.625, 0.565, 0.515, 0.49, 0.475, 0.47, 0.47,
    0.785, 0.71, 0.645, 0.6, 0.57, 0.555, 0.54, 0.54, 0.54,
    0.785, 0.735, 0.695, 0.67, 0.645, 0.63, 0.62, 0.61, 0.6,
    0.785, 0.75, 0.73, 0.71, 0.7, 0.69, 0.685, 0.68, 0.675,
    0.785, 0.68, 0.6, 0.54, 0.5, 0.48, 0.475, 0.485, 0.505,
    0.785, 0.71, 0.645, 0.6, 0.58, 0.565, 0.56, 0.565, 0.57,
    0.785, 0.735, 0.7, 0.675, 0.66, 0.65, 0.64, 0.635, 0.63,
    0.785, 0.8, 0.81, 0.81, 0.8, 0.79, 0.77, 0.74, 0.71,
    0.785, 0.81, 0.83, 0.845, 0.855, 0.855, 0.85, 0.835, 0.81,
    0.785, 0.685, 0.605, 0.55, 0.51, 0.505, 0.53, 0.58, 0.645,
    0.785, 0.71, 0.67, 0.64, 0.635, 0.65, 0.67, 0.72, 0.77,
    0.785, 0.775, 0.78, 0.79, 0.805, 0.82, 0.845, 0.875, 0.905,
    0.785, 0.815, 0.84, 0.865, 0.885, 0.905, 0.915, 0.925, 0.93,
    0.785, 0.83, 0.87, 0.9, 0.92, 0.935, 0.94, 0.94, 0.935,
])


def _floats(values, n):
    return [float(v) for v in list(values)[:n]]


def calculate_supwbt(data: Mapping[str, object],
                     user_downwash: bool = False) -> Dict[str, object]:
    """Translate SUPWBT: supersonic wing-body-tail aerodynamics.

    Args:
        data: The COMMON words the routine reads, by name:

            - ``alpha`` (``FLC(23..)``), ``mach`` (``FLC(I+2)``), ``sref``,
              ``cbarr``.
            - ``position``: ``xcg``, ``aliw``, ``zcg``, ``xh``, ``zh``,
              ``alih`` (``/SYNTSS/``).
            - ``dwash``: ``q`` (``QOQINF``, ``DWASH(1..)``), ``eps``
              (``EPSLON``, ``DWASH(21..)``), ``dedalp`` (``DWASH(41..)``),
              ``hmach`` (``DWA(189..)``) and ``jdetch`` (``DWA(237)``).
            - ``tail``: ``span``, ``spans`` (``HTIN(4)``, ``(3)``), ``a``
              (the ``AHT`` words 3, 7, 10, 27, 34, 38, 62, 161), ``cd``,
              ``cl``, ``cla`` (``HT(1..)``, ``(21..)``, ``(101..)``),
              ``alpha`` (``STG(33..)``, the tail's angle table), ``cd0``
              (``STG(80)``), ``xac`` and ``xac136`` (``STG(134)``,
              ``(136)``).
            - ``wing``: ``span``, ``spans``, ``type`` (``WINGIN(4)``,
              ``(3)``, ``(101)``), ``a`` (the ``A`` words 7, 10, 12, 24,
              27, 38, 80, 120), ``cla`` (``WING(101)``), ``xac``,
              ``xac136`` (``SLG(134)``, ``(136)``), ``kwb``, ``kkwb``,
              ``kkbw``, ``cd0`` (``SWB(35)``, ``(2)``, ``(37)``, ``(4)``).
            - ``wing_body``: ``cd``, ``cl`` (``BW(1..)``, ``(21..)``),
              ``cla``, ``cma`` (``BW(101)``, ``(121)``).
            - ``body``: ``x``, ``s`` (``/BODYI/``), ``dn``, ``d1``
              (``SBD(4)``, ``(5)``), ``bd68``, ``bd78``.
            - ``stale``: ``STP`` words the routine reads without setting,
              or sets only on some paths: ``cd0v`` (1), ``cd0vf`` (156),
              ``kkbw`` (131), ``kkwb`` (132).

        user_downwash: ``KDEODA`` and ``KQOQIN`` both nonzero (downwash
            and dynamic pressure supplied as experimental data): the tail
            Mach number is then the free stream's.  Otherwise the routine
            returns at once unless SDWASH reached a detachment angle
            (``JDETCH``), and a positive ``JDETCH`` shortens the angle
            loop to that many angles.

    Returns:
        ``None`` when the routine returns without computing; otherwise the
        words set: ``nalpha``, ``kbw``, ``kwb``, ``kkbw``, ``kkwb``,
        ``clahb``, ``clabh``, ``dd``, ``trino``, ``rkbw``, ``yt``,
        ``rcreo2`` (the body radius BODOWG leaves), ``dxacwb``,
        ``cd0wbt``, ``cd0wbv``, the ``STP`` curves ``cmah``, ``cltb``,
        ``cdawb``, ``ivwh``, ``deltat``, ``gamma``, ``ivbh``, the ``BWH``
        curves ``cd``, ``cl``, ``cm``, ``cn``, ``ca``, ``cla``, ``cma``,
        the ``BWHV`` curves ``cdv``, ``cnv``, ``cav``, the canard's
        ``fact102`` and ``fact122`` (``FACT(J+101)``, ``FACT(J+121)``),
        and ``bd`` (the ``BD`` words 2, 3, 58, 63, 64, 84, 761, 762).

    Notes:
        Kept as executed:

        - The canard path reads ``EPSLON(J+20)``, which is ``DEDALP(J)``
          by the ``DWASH`` layout.
        - The wing-vortex interference ``ALI`` is handed ``SPAN-SPANS``
          as the body radius before BODOWG, later in the routine,
          overwrites that word with the radius from the maximum area.
        - The canard's vortex height ``A(24)-A(80)*BVTO2`` subtracts a
          length times a tangent from a length.
    """
    alpha = [float(v) for v in data['alpha']]
    mach = float(data['mach'])
    sref, cbarr = float(data['sref']), float(data['cbarr'])
    pos = {k: float(v) for k, v in data['position'].items()}
    dw = data['dwash']
    tl, wg = data['tail'], data['wing']
    wb, body = data['wing_body'], data['body']
    st = {k: float(v) for k, v in data['stale'].items()}
    aht = {int(k): float(v) for k, v in tl['a'].items()}
    aw = {int(k): float(v) for k, v in wg['a'].items()}

    canard = float(wg['type']) > 2.5
    nalph = len(alpha)
    nalpha = nalph
    na1 = nalph
    for j in range(nalph):
        if float(tl['cl'][j]) == UNUSED and float(wb['cl'][j]) == UNUSED:
            na1 = j
    tanle = aht[62] if aht[62] != 0.0 else .00001
    jdetch = int(dw['jdetch'])
    if user_downwash:
        hmach = [mach] * nalpha
    else:
        if jdetch == 0:
            return None
        if jdetch > 0:
            nalpha = jdetch
        hmach = _floats(dw['hmach'], nalpha)
    span, spans = float(tl['span']), float(tl['spans'])
    beta = math.sqrt(mach**2 - 1.)
    rlb = float(body['x'][-1])
    dd = 2.0 * (span - spans)
    nalpha = min(nalpha, na1)
    arstar, crstar, tapexp = aht[7], aht[10], aht[27]
    clah = [float(v) for v in tl['cla']]
    r: Dict[str, object] = {'dd': dd}

    if tapexp == 0.0:
        integrate = beta * arstar > 1.
    else:
        trino = beta * arstar * (1.0 + tapexp) * (1. + tanle / beta)
        r['trino'] = trino
        integrate = trino > 4.
    if integrate:
        dx = rlb - pos['xh'] - dd * aht[38] - crstar
        rkbw = 0.0 if dx <= -crstar else intkbw(mach, aht[34], crstar, dd,
                                                 dx)[0]
        kbw = rkbw / (RAD * beta * (sref / aht[3]) * clah[0] *
                      (tapexp + 1.) * (2. * span / dd - 1.))
        r['rkbw'] = rkbw
    else:
        kbw = interx(1, _TFIG10, [dd / (2. * span)], [11], _DKBW10, lind=11)
    ratio = (span - spans) / span
    kwb = interx(1, _TFIG10, [ratio], [11], _DKWB10, lind=11)
    clahb, clabh = clah[0] * kwb, clah[0] * kbw

    q = _floats(dw['q'], nalpha)
    eps = _floats(dw['eps'], nalpha)
    dedalp = _floats(dw['dedalp'], nalpha)
    clawb = float(wb['cla'])
    rcreo2 = span - spans
    ivwh, deltat, fact122 = [], [], []
    if canard:
        dbwo2 = float(wg['span']) - float(wg['spans'])
        yt = interx(3, _T4467, [beta * aw[7], aw[120] * aw[38], aw[27]],
                    [9, 5, 3], _D4467, lind=9, lx1u=1, lx2u=1, lx3u=1)
        bvto2 = float(wg['spans']) * yt + dbwo2
        r['yt'] = yt
        for j in range(nalpha):
            z0 = aw[12] - (aw[24] - aw[80] * bvto2) * math.tan(
                (alpha[j] + pos['aliw']) / RAD)
            ivwh.append(ali(z0, bvto2, spans, rcreo2, tapexp))
            anum = (float(wg['cla']) * clah[0] * q[j] * float(wg['kwb']) *
                    spans * ivwh[j])
            aden = 2.0 * PI * arstar * (bvto2 - dbwo2)
            deltat.append((anum / aden) * RAD * sref / aht[3])
            fact122.append(-deltat[j] / ((clahb + clabh) * q[j]))
        clawbt = [clawb + (clahb + clabh) * q[j] + deltat[j]
                  for j in range(nalpha)]
    else:
        clawbt = [clawb + (clahb + clabh) * (1.0 - dedalp[j]) * q[j]
                  for j in range(nalpha)]

    alih, aliw = pos['alih'], pos['aliw']
    kkbw, kkwb = st['kkbw'], st['kkwb']
    if not (alih == 0.0 or alih == UNUSED):
        kkbw = interx(1, _T4312A, [ratio], [11], _D4312A, lind=11)
        kkwb = interx(1, _T4312B, [ratio], [8], _D4312B, lind=8)
    cli = clah[0] * alih
    xcbr4h = pos['xh'] + aht[161]
    xmax, smax, _ = getmax(body['x'], body['s'])
    rcreo2 = math.sqrt(smax / PI)
    vortex = [calculate_bodowg(a, xcbr4h, rcreo2, span, tapexp)
              for a in alpha[:nalpha]]
    ivbh = [v['ivbw'] for v in vortex]
    gamma = [v['go2pav'] for v in vortex]
    arg1 = kwb + kbw
    arg2 = (kkwb + kkbw) * cli
    arg3 = clah[0] * (dd / (2. * span))
    sinai, cosai = math.sin(alih / RAD), math.cos(alih / RAD)
    tanai = math.tan(alih / RAD)
    bd58 = (span - spans) * tanle * cosai
    bd84 = pos['xcg'] - (pos['xh'] + bd58)
    bd762 = pos['zh'] - bd58 * tanai
    bd761 = float(tl['xac']) * crstar
    bd64 = bd762 - bd761 * sinai - pos['zcg']
    bd63 = bd84 - bd761 * cosai

    fl = np.asarray(alpha[:nalpha])
    stga = np.asarray(_floats(tl['alpha'], nalpha))
    cdwb = _floats(wb['cd'], nalpha)
    clwb = _floats(wb['cl'], nalpha)
    htcl = _floats(tl['cl'], nalpha)
    htcla = clah[:nalpha]
    htcd = _floats(tl['cd'], nalpha)
    cmawb = float(wb['cma'])
    xacw, xach = float(wg['xac']), float(tl['xac'])
    cd0wbt = float(wg['cd0']) + float(tl['cd0'])
    cd0wbv = cd0wbt + st['cd0v'] + st['cd0vf']
    cdawb, cltb, clwbt, fact102, cmah, cmawbt, cmwbt = ([] for _ in
                                                         range(7))
    cdwbt, cdwbtv, clh = [], [], []
    dxacwb = cmawb / clawb
    for j in range(nalpha):
        cdawb.append(tbfunx(fl, cdwb, alpha[j], 0, 0)[1])
        alpaht = alpha[j] - eps[j] + alih
        alpat = (stga[j] * RAD - eps[j]) / RAD
        if canard:
            alpat = stga[j]
        clh.append(tbfunx(stga, htcl, alpat, 0, 0)[0])
        claha = tbfunx(stga, htcla, alpat, 0, 0)[0]
        detcl = arg3 * alpaht * ivbh[j] * gamma[j] * q[j]
        dcloal = arg3 * ivbh[j] * gamma[j] * q[j]
        alpa = alpha[j] - eps[j]
        mratio = ((mach**2 - 1.0) / (hmach[j]**2 - 1.0))**0.50
        cltb.append((arg2 + (clh[j] - cli) * arg1) * q[j] * mratio)
        if canard:
            clwbt.append(clwb[j] + cltb[j] + detcl + deltat[j] * alpha[j])
            f = 0.0
            if alpaht != 0.0:
                f = ((dcloal * dedalp[j] - deltat[j] * alpha[j]) /
                     (dcloal + claha * arg1 * q[j] * mratio))
            fact102.append(f)
            alpa = alpha[j]
        else:
            clwbt.append(clwb[j] + cltb[j] + detcl)
        sa, ca = math.sin(alpa / RAD), math.cos(alpa / RAD)
        sf, cf = math.sin(alpha[j] / RAD), math.cos(alpha[j] / RAD)
        apart = dxacwb * ((-clwb[j] / RAD + cdawb[j]) * sf +
                          (clawb + cdwb[j] / RAD) * cf)
        zac = float(body['bd68']) - pos['zcg'] - float(body['bd78']) * \
            xacw * aw[10]
        bpart = (zac / cbarr) * ((clawb + cdwb[j] / RAD) * sf +
                                 (clwb[j] / RAD - cdawb[j]) * cf)
        cdht, dcdda = tbfunx(stga, htcd, alpat, 0, 0)
        clht = (clwbt[j] - clwb[j]) / q[j]
        dclda = clahb + clabh
        if canard:
            dclda = dclda + deltat[j] / q[j]
        cpart = (-clht + dcdda) * sa / RAD
        dpart = (clht - dcdda) * ca / RAD
        epart = q[j] * (1.0 - dedalp[j])
        fpart = (dclda + cdht / RAD) * sa
        gpart = (dclda + cdht / RAD) * ca
        cm_h = ((bd63 / cbarr) * (cpart + gpart) * epart -
                (bd64 / cbarr) * (fpart + dpart) * epart)
        if canard:
            cm_h = cm_h / (1.0 - dedalp[j])
        cmah.append(cm_h)
        cmawbt.append(apart - bpart + cm_h)
        cmwbt.append(cmawbt[j] * alpha[j] +
                     (float(wg['xac136']) / RAD) *
                     (float(wg['kkwb']) + float(wg['kkbw'])) * aliw +
                     (float(tl['xac136']) / RAD) * (kkbw + kkwb) * alih *
                     q[j])
        cdwbt.append(cdwb[j] + (cdht * math.cos(eps[j] / RAD) +
                                clh[j] * math.sin(eps[j] / RAD)) * q[j])
        cdwbtv.append(cdwbt[j] + st['cd0v'] + st['cd0vf'])
    cn, ca_, cnv, cav = [], [], [], []
    for j in range(nalpha):
        c, s = math.cos(alpha[j] / RAD), math.sin(alpha[j] / RAD)
        cn.append(clwbt[j] * c + cdwbt[j] * s)
        ca_.append(cdwbt[j] * c - clwbt[j] * s)
        cnv.append(clwbt[j] * c + cdwbtv[j] * s)
        cav.append(cdwbtv[j] * c - clwbt[j] * s)
    r.update({
        'nalpha': nalpha, 'kbw': float(kbw), 'kwb': float(kwb),
        'kkbw': float(kkbw), 'kkwb': float(kkwb), 'clahb': float(clahb),
        'clabh': float(clabh), 'rcreo2': rcreo2, 'dxacwb': dxacwb,
        'cd0wbt': cd0wbt, 'cd0wbv': cd0wbv, 'cmah': cmah, 'cltb': cltb,
        'cdawb': cdawb, 'ivwh': ivwh, 'deltat': deltat, 'gamma': gamma,
        'ivbh': ivbh, 'cd': cdwbt, 'cl': clwbt, 'cm': cmwbt, 'cn': cn,
        'ca': ca_, 'cla': clawbt, 'cma': cmawbt, 'cdv': cdwbtv,
        'cnv': cnv, 'cav': cav, 'fact102': fact102, 'fact122': fact122,
        'bd': {2: xmax, 3: smax, 58: bd58, 63: bd63, 64: bd64, 84: bd84,
               761: bd761, 762: bd762},
        'method': 'legacy_supwbt'})
    return r
