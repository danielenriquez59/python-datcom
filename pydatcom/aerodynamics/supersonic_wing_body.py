"""
Supersonic wing-body and tail-body: SUPWB, SUPHB and their overlay M20O24.

The supersonic counterpart of WBTRAN plus the curve buildup: the body
carryover K_B(W) by INTKBW's integration or Figure 4.3.1.2-10, K_W(B),
the incidence carryover and body-vortex lift for a straight surface, the
aerodynamic centre from the nose, surface and carryover centres, and the
surface-body lift, drag, force and moment curves.  SUPHB is SUPWB on the
tail's blocks, with the differences listed on :func:`calculate_supwb`.

The tables were extracted from the source DATA statements by parsing, with
``tools/fortran_data.py``, rather than by hand; the Figure 4.3.1.2-10 and
4.3.2.2-37A/B tables are WBTRAN's, pinned equal by test.

Reference: datcom-legacy/datcom_2000/supwb.f, suphb.f, m20o24.f
"""

import math
from typing import Dict, Mapping, Optional, Sequence

import numpy as np

from pydatcom.aerodynamics.cdrag import STRAIGHT_TAPERED
from pydatcom.aerodynamics.transonic_buildup import (_D4337A, _D4337B,
                                                     _DKBW10, _DKWB10,
                                                     _T4337A, _T4337B,
                                                     _TFIG10)
from pydatcom.interactions.body_vortex import calculate_bodowg, getmax
from pydatcom.interactions.carryover import intkbw
from pydatcom.utils.constants import PI, RAD, UNUSED
from pydatcom.utils.legacy_interp import interx
from pydatcom.utils.legacy_numeric import tbfunx

# Figure 4.2.1.2-18A/B: the nose centre of pressure, ahead of and past the Mach line (pointed and blunt noses).
_T4218A = np.array([
    0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 0.0, 0.0,
    0.0, 0.0, 0.4, 0.8, 1.2, 1.6, 2.0, 3.0,
    4.0, 5.0,
])
_DL218A = np.array([
    0.543, 0.542, 0.541, 0.54, 0.534, 0.526, 0.4, 0.409,
    0.418, 0.43, 0.441, 0.45, 0.305, 0.328, 0.35, 0.369,
    0.387, 0.4, 0.238, 0.265, 0.295, 0.318, 0.339, 0.356,
    0.198, 0.221, 0.246, 0.274, 0.298, 0.32, 0.16, 0.185,
    0.21, 0.239, 0.262, 0.288, 0.065, 0.095, 0.122, 0.15,
    0.177, 0.21, 0.0, 0.005, 0.035, 0.062, 0.089, 0.13,
    0.0, 0.0, 0.0, 0.0, 0.002, 0.05,
])
_DR218A = np.array([
    0.445, 0.464, 0.485, 0.5, 0.518, 0.526, 0.448, 0.455,
    0.46, 0.46, 0.459, 0.45, 0.46, 0.449, 0.438, 0.424,
    0.412, 0.4, 0.45, 0.43, 0.412, 0.394, 0.375, 0.356,
    0.432, 0.41, 0.388, 0.365, 0.343, 0.32, 0.42, 0.394,
    0.369, 0.34, 0.314, 0.288, 0.388, 0.354, 0.322, 0.278,
    0.244, 0.21, 0.357, 0.314, 0.273, 0.216, 0.171, 0.13,
    0.325, 0.274, 0.225, 0.154, 0.1, 0.05,
])
_T4218B = np.array([
    0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 0.0, 0.0,
    0.5, 1.0, 2.0, 3.0, 4.0, 5.0,
])
_DL218B = np.array([
    0.665, 0.665, 0.665, 0.665, 0.665, 0.665, 0.425, 0.492,
    0.539, 0.55, 0.55, 0.55, 0.33, 0.37, 0.405, 0.438,
    0.459, 0.47, 0.184, 0.215, 0.25, 0.284, 0.318, 0.35,
    0.06, 0.097, 0.133, 0.17, 0.206, 0.24, 0.0, 0.0,
    0.044, 0.083, 0.127, 0.17, 0.0, 0.0, 0.0, 0.02,
    0.063, 0.105,
])
_DR218B = np.array([
    0.665, 0.665, 0.665, 0.665, 0.665, 0.665, 0.48, 0.5,
    0.519, 0.536, 0.546, 0.55, 0.338, 0.388, 0.43, 0.458,
    0.471, 0.47, 0.338, 0.372, 0.394, 0.394, 0.375, 0.35,
    0.41, 0.375, 0.341, 0.308, 0.272, 0.24, 0.377, 0.338,
    0.294, 0.251, 0.211, 0.17, 0.342, 0.3, 0.246, 0.194,
    0.146, 0.1,
])
# Figure 4.3.1.2-12A/B: the incidence carryover, K_B(W) and K_W(B).
_T4312A = np.array([
    0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7,
    0.8, 0.9, 1.0,
])
_D4312A = np.array([
    0.0, 0.11, 0.21, 0.31, 0.41, 0.51, 0.6, 0.7,
    0.8, 0.9, 1.0,
])
_T4312B = np.array([
    0.015, 0.1, 0.2, 0.3, 0.4, 0.6, 0.8, 0.975,
])
_D4312B = np.array([
    1.0, 0.975, 0.956, 0.947, 0.941, 0.95, 0.978, 1.0,
])


def calculate_supwb(alpha_deg: Sequence[float], mach: float,
                    surface: Mapping[str, float], a: Mapping[int, float],
                    position: Mapping[str, float], body: Mapping[str, object],
                    aero: Mapping[str, object], sref: float, cbarr: float,
                    tail: bool = False,
                    moment_block: Optional[Mapping[str, object]] = None,
                    stale: Optional[Mapping[str, object]] = None
                    ) -> Dict[str, object]:
    """Translate SUPWB (and SUPHB): supersonic surface-body aerodynamics.

    K_B(W) by INTKBW's integration above a similarity parameter of four
    (untapered: above ``beta*A = 1``) or Figure 4.3.1.2-10 below it; K_W(B)
    from Figure 4.3.1.2-10; for a straight surface the incidence carryover
    of Figure 4.3.1.2-12A/B and BODOWG's vortex lift; the nose centre of
    pressure (Figure 4.2.1.2-18A/B) and the carryover centre (Figure
    4.3.2.2-37A/B) for the aerodynamic centre and CMa; then CD, CN, CA and
    the moment curve.

    Args:
        alpha_deg, mach: ``FLC(23..)``, ``FLC(I+2)``.
        surface: ``spans``, ``span``, ``cr``, ``type`` (``WINGIN(3)``,
            ``(4)``, ``(6)``, ``(15)``).
        a: The surface's ``A`` block: 3, 7, 10, 27, 34, 38, 62, 161, 173.
        position: ``x`` (``XW``, or ``XH`` for SUPHB), ``incidence``
            (``ALIW``/``ALIH``), and ``zw``, ``zcg`` (``SYNA(3)``, ``(5)``).
        body: ``x``, ``s`` (the ``/BODYI/`` stations, for BODOWG), ``rln``,
            ``bnose``, ``dn``, ``d1`` (``SBD(4)``, ``(5)``), ``alpha0``
            (``BD(81)``), ``cla`` (``SBD(18)``), ``cd0`` (``SBD(124)``),
            ``cl``, ``cd``, ``cm`` (``BODY(21..)``, ``(1..)``, ``(41..)``).
        aero: The surface's ``cla`` (``WING(101)``), ``cl`` (``(21..)``),
            ``cn``, ``ca`` (``(61..)``, ``(81..)``), ``cd0`` (``SLG(80)``),
            ``cdl`` (``SLG(53..)``) and ``xac`` (``SLG(134)``, root chords).
        sref, cbarr: ``/OPTION/``.
        tail: SUPHB's afterbody length ``DD*AHT(38)`` (SUPWB halves it) and
            moment loop (below).
        moment_block: For SUPHB, what its moment loop reads from the
            *wing's* blocks: ``kwb``, ``kkwb``, ``kbw``, ``kkbw``,
            ``xacbw``, ``dd``, ``ivbw``, ``gamma`` (``SWB``), ``xacw``
            (``SLG(134)``) and ``incidence`` (``SYNA(4)``).
        stale: Words the routine leaves when it does not set them:
            ``kkwb``, ``kkbw`` (a straight surface at zero or unused
            incidence), and for a surface that is not straight ``cl``
            (``BW(21..)``), ``ivbw`` and ``gamma``.

    Returns:
        ``kbw``, ``kwb``, ``kkbw``, ``kkwb``, ``clawb``, ``clabw``, ``cla``
        (``BW(101)``), ``cli``, ``ivbw``, ``gamma``, ``xcpln``, ``xacn``,
        ``xaca``, ``xacbw``, ``xacw`` (rescaled to ``CBARR``), ``xac``,
        ``cma`` (``BW(121)``), ``cd0``, the curves ``cl``, ``cd``, ``cn``,
        ``ca``, ``cm``, ``alpha_body`` (``BD(255..)``), and the other
        ``SWB`` words set (``dd``, ``beta``, ``rkbw``, ``trino``, ``fa``,
        ``fn``, ``rlap``, ``delxw``).

    Notes:
        SUPHB's moment loop reads the wing's carryover factors, centres,
        vortex factors, span and incidence (its ``/SUPWBB/`` and
        ``/SUPWH/`` declarations put the wing's block first), and its
        vortex term omits the ``/RAD`` of SUPWB's.  SUPWB's own vortex
        moment divides by ``RAD`` where its lift term does not.  BODOWG
        is handed the exposed-root offset and assigns it the body radius.
        All kept.
    """
    st = dict(stale or {})
    g = {int(k): float(v) for k, v in a.items()}
    alpha = np.asarray(alpha_deg, dtype=float)
    nalpha = len(alpha)
    span, spans = float(surface['span']), float(surface['spans'])
    cr = float(surface['cr'])
    arstar, crstar, tapexp = g[7], g[10], g[27]
    claw = float(aero['cla'])
    xw, aliw = float(position['x']), float(position['incidence'])
    xs = np.asarray(body['x'], dtype=float)
    rlb = float(xs[-1])
    dcyl = (float(body['dn']) + float(body['d1'])) / 2.
    beta = math.sqrt(mach**2 - 1.)
    dd = 2.0 * (span - spans)
    tanle = g[62] if g[62] != 0.0 else .00001
    r: Dict[str, object] = {'dd': dd, 'beta': beta}
    if tapexp == 0.0:
        integrate = beta * arstar > 1.
    else:
        trino = beta * arstar * (1.0 + tapexp) * (1. + tanle / beta)
        r['trino'] = trino
        integrate = trino > 4.
    if integrate:
        dx = rlb - xw - (dd if tail else dd / 2.) * g[38] - crstar
        rkbw = 0.0 if dx <= -crstar else intkbw(mach, g[34], crstar, dd,
                                                 dx)[0]
        kbw = rkbw / (RAD * beta * (sref / g[3]) * claw * (tapexp + 1.) *
                      (2. * span / dd - 1.))
        r.update({'rkbw': rkbw, 'dx': dx})
    else:
        kbw = interx(1, _TFIG10, [dd / (2. * span)], [11], _DKBW10, lind=11)
    albo = 0.0 if float(body['alpha0']) == UNUSED else float(body['alpha0'])
    alphab = alpha + albo
    ratio = (span - spans) / span
    kwb = interx(1, _TFIG10, [ratio], [11], _DKWB10, lind=11)
    clab = float(body['cla'])
    clawb, clabw = claw * kwb, claw * kbw
    cla = clabw + clawb + clab
    kkbw = float(st.get('kkbw', 0.0))
    kkwb = float(st.get('kkwb', 0.0))
    cl = [float(v) for v in st.get('cl', [0.0] * nalpha)]
    ivbw = [float(v) for v in st.get('ivbw', [0.0] * nalpha)]
    gamma = [float(v) for v in st.get('gamma', [0.0] * nalpha)]
    cli = None
    if float(surface['type']) == STRAIGHT_TAPERED:
        if not (aliw == 0.0 or aliw == UNUSED):
            kkbw = interx(1, _T4312A, [ratio], [11], _D4312A, lind=11)
            kkwb = interx(1, _T4312B, [ratio], [8], _D4312B, lind=8)
        cli = claw * aliw
        _, max_area, _ = getmax(body['x'], body['s'])
        radius = math.sqrt(max_area / PI)
        vortex = [calculate_bodowg(ab, xw + g[161], radius, span, tapexp)
                  for ab in alphab]
        ivbw = [v['ivbw'] for v in vortex]
        gamma = [v['go2pav'] for v in vortex]
        arg3 = claw * (dd / (2. * span))
        clw, clb = aero['cl'], body['cl']
        cl = [
            float(clb[angle_slot]) +
            (kwb + kbw) * (float(clw[angle_slot]) - cli) +
            (kkwb + kkbw) * cli +
            arg3 * alphab[angle_slot] * ivbw[angle_slot] * gamma[angle_slot]
            for angle_slot in range(nalpha)
        ]

    delxw = (span - spans) * tanle * math.cos(aliw / RAD)
    rln = float(body['rln'])
    rlap = xw + delxw - rln
    arg1 = rln
    if rlap < 0.:
        arg1 = rln + rlap
        rlap = 0.0
    fa, fn = rlap / dcyl, arg1 / dcyl
    var1, var2 = beta / fn, fa / fn
    if float(body['bnose']) == 1.:
        t, dl, dr, n2 = _T4218B, _DL218B, _DR218B, 7
    else:
        t, dl, dr, n2 = _T4218A, _DL218A, _DR218A, 9
    if var1 > 1.:
        xcpln = interx(2, t, [1. / var1, var2], [6, n2], dr, lind=n2)
    else:
        xcpln = interx(2, t, [var1, var2], [6, n2], dl, lind=n2)
    xacn = (xcpln - 1.) * (arg1 + rlap) / cbarr
    var = [beta * dd / crstar, beta / tanle]
    if (xw + cr) / rlb > 1.:
        xaca = interx(2, _T4337B, var, [10, 2], _D4337B, lind=10)
    else:
        xaca = interx(2, _T4337A, var, [8, 3], _D4337A, lind=8, lx1u=1)
    xacbw = xaca * crstar / cbarr
    xacw = float(aero['xac']) * crstar / cbarr
    xac = (xacn * clab + xacw * clawb + xacbw * clabw) / (clab + clawb +
                                                          clabw)
    cma = (g[173] / cbarr - xac) * cla
    cd0w = float(aero['cd0'])

    mb = dict(moment_block) if tail else {
        'kwb': kwb, 'kkwb': kkwb, 'kbw': kbw, 'kkbw': kkbw, 'xacbw': xacbw,
        'dd': dd, 'ivbw': ivbw, 'gamma': gamma, 'xacw': xacw,
        'incidence': aliw}
    zarm = (float(position['zw']) - float(position['zcg'])) / cbarr
    dxcpwb = g[173] / cbarr - float(mb['xacw'])
    dxcpbw = g[173] / cbarr - float(mb['xacbw'])
    cd, cn, ca, cm = [], [], [], []
    for angle_slot in range(nalpha):
        cosa, sina = (math.cos(alpha[angle_slot] / RAD),
                      math.sin(alpha[angle_slot] / RAD))
        cdj = (cd0w + float(aero['cdl'][angle_slot]) +
               float(body['cd'][angle_slot]))
        dcnv = (float(mb['ivbw'][angle_slot]) *
                float(mb['gamma'][angle_slot]) *
                float(mb['dd']) / (2 * span) * alphab[angle_slot] * claw)
        if not tail:
            dcnv = dcnv / RAD
        cnw = float(aero['cn'][angle_slot])
        caw = float(aero['ca'][angle_slot])
        inc = float(mb['incidence'])
        cm.append(float(body['cm'][angle_slot]) +
                  cnw * float(mb['kwb']) * dxcpwb +
                  claw * inc * float(mb['kkwb']) * dxcpwb / RAD +
                  cnw * float(mb['kbw']) * dxcpbw +
                  claw * inc * float(mb['kkbw']) * dxcpbw / RAD +
                  dcnv * (dxcpbw if tail else dxcpwb) + caw * zarm)
        cd.append(cdj)
        cn.append(cl[angle_slot] * cosa + cdj * sina)
        ca.append(cdj * cosa - cl[angle_slot] * sina)
    r.update({
        'kbw': float(kbw), 'kwb': float(kwb), 'kkbw': float(kkbw),
        'kkwb': float(kkwb), 'clawb': float(clawb), 'clabw': float(clabw),
        'cla': float(cla), 'cli': cli, 'ivbw': ivbw, 'gamma': gamma,
        'xcpln': float(xcpln), 'xacn': float(xacn), 'xaca': float(xaca),
        'xacbw': float(xacbw), 'xacw': float(xacw), 'xac': float(xac),
        'cma': float(cma), 'cd0': cd0w + float(body['cd0']), 'fa': fa,
        'fn': fn, 'rlap': rlap, 'delxw': delxw, 'alpha_body': alphab,
        'cl': cl, 'cd': cd, 'cn': cn, 'ca': ca, 'cm': cm,
        'method': 'legacy_suphb' if tail else 'legacy_supwb'})
    return r


def m20o24_slopes(alpha_deg: Sequence[float], result: Mapping[str, object]
                  ) -> Dict[str, list]:
    """M20O24's pass over SUPWB's or SUPHB's curves.

    CN and CA are re-formed from CL and CD (as SUPWB forms them), the
    first angle keeps the routine's own CLa and CMa, and the others take
    TBFUNX's slope of CL for CLa and ``-UNUSED`` for CMa.
    """
    alpha = np.asarray(alpha_deg, dtype=float)
    cl = np.asarray(result['cl'], dtype=float)
    cd = np.asarray(result['cd'], dtype=float)
    ca_, sa_ = np.cos(alpha / RAD), np.sin(alpha / RAD)
    cla = [float(result['cla'])] + [float(tbfunx(alpha, cl, x, 0, 0)[1])
                                    for x in alpha[1:]]
    cma = [float(result['cma'])] + [-UNUSED] * (len(alpha) - 1)
    return {'cn': list(cl * ca_ + cd * sa_), 'ca': list(cd * ca_ - cl * sa_),
            'cla': cla, 'cma': cma}
