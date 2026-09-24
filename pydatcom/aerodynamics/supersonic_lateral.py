"""
Supersonic lateral stability of the wing and wing-body (SUPLAT) and of the
horizontal tail and tail-body (SUPLAH).

The wing alone's CY_beta and Cn_beta for a rectangular planform with
supersonic edges or a delta with subsonic ones (Figures 5.1.1.1-6 and
7.1.1.1-8), its Cl_beta from the dihedral and the roll damping of Figure
7.1.2.2-25; the body's side force, the body yawing moment of Figure
5.2.3.1-8 with the Reynolds-number factor of Figure 5.2.3.1-9, and the
wing-height rolling moment; and in SUPLAT the horizontal tail's side-
force increment on the body (Figure 5.3.1.1-25OO).  SUPLAH is the same
buildup on the tail's blocks.

The tables were extracted from the source DATA statements by parsing, with
``tools/fortran_data.py``.

Reference: datcom-legacy/datcom_2000/suplat.f, suplah.f
"""

import math
from typing import Dict, Mapping

import numpy as np

from pydatcom.utils.constants import PI, RAD, UNUSED
from pydatcom.utils.legacy_interp import interx
from pydatcom.utils.legacy_numeric import tbfunx, trapz
from pydatcom.utils.legacy_tables import tlinex

# Figure 5.2.3.1-8A/B/C: the body yawing-moment factor K_N.
_X158A = np.array([
    20.0, 14.0, 10.0, 8.0, 7.0, 6.0, 5.0, 4.0, 3.0, 2.5,
])
_X258A = np.array([
    0.2, 0.8,
])
_Y58A = np.array([
    0.1, 1.88, 0.4, 2.21, 0.74, 2.6, 0.98, 2.8, 1.3, 3.13,
    1.61, 3.5, 2.0, 3.88, 2.5, 4.4, 2.99, 5.0, 3.45, 5.4,
])
_X158B = np.array([
    0.8, 1.0, 1.2, 1.4, 1.6,
])
_X258B = np.array([
    0.0, 3.0, 6.0,
])
_Y58B = np.array([
    0.0, 2.35, 4.68, 0.0, 3.0, 6.0, 0.0, 3.6, 7.25, 0.0,
    4.18, 8.5, 0.0, 4.79, 9.5,
])
_X158C = np.array([
    0.5, 0.6, 0.8, 1.0, 2.0,
])
_X258C = np.array([
    0.0, 6.0,
])
_Y58C = np.array([
    -0.00048, 0.00251, -0.00048, 0.0035, -0.00048, 0.00477, -0.00048, 0.00559, -0.00048, 0.00641,
])
# Figure 7.1.1.1-8 (E) and 5.1.1.1-6 (1/Q): elliptic-integral factors.
_X71118 = np.array([
    0.0, 0.05, 0.1, 0.165, 0.25, 0.35, 0.8, 1.0,
])
_Y71118 = np.array([
    1.0, 0.995, 0.985, 0.966, 0.94, 0.9, 0.705, 0.631,
])
_X51116 = np.array([
    0.0, 0.1, 0.2, 0.45, 0.55, 0.65, 0.7, 0.75,
    0.8, 0.85, 0.875, 0.9, 0.95, 0.974, 0.99, 1.0,
])
_Y51116 = np.array([
    1.0, 1.028, 1.078, 1.241, 1.284, 1.3, 1.287, 1.26,
    1.207, 1.125, 1.05, 0.97, 0.75, 0.55, 0.35, 0.0,
])
# Figure 5.3.1.1-25OO: the horizontal tail's apparent-mass factor.
_XA25OO = np.array([
    0.0, 0.2, 0.4, 0.6, 0.8, 1.0,
])
_YA25OO = np.array([
    0.0, 0.01, 0.04, 0.12, 0.22, 0.35,
])
_XB25OO = np.array([
    0.0, 0.1, 0.2, 0.3, 0.4, 0.6, 0.8, 1.0,
])
_YB25OO = np.array([
    1.275, 1.24, 1.18, 1.08, 0.95, 0.69, 0.51, 0.35,
])
# Figure 7.1.2.2-25A-E: roll damping, against beta*A (17), A*tan(c/2)
# (7, padded to 17) and taper (5).
_X12225 = np.array([
    0.0, 1.0, 2.0, 2.45, 3.0, 3.31, 4.0, 4.7, 5.0, 5.3, 5.75, 6.0, 6.7, 7.0, 8.0, 9.0, 10.0,
    0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.25, 0.5, 0.75, 1.0,
])
_Y12225 = np.array([
    -0.098, -0.1025, -0.1067, -0.1088, -0.0954, -0.087, -0.075, -0.065, -0.0615, -0.0585, -0.054, -0.052, -0.047, -0.0455, -0.0408, -0.0363, -0.0338,
    -0.098, -0.099, -0.1008, -0.1015, -0.099, -0.09, -0.077, -0.0666, -0.063, -0.06, -0.055, -0.0535, -0.048, -0.0465, -0.0414, -0.0368, -0.0335,
    -0.098, -0.0958, -0.0925, -0.0908, -0.0888, -0.087, -0.085, -0.071, -0.0665, -0.064, -0.0578, -0.0555, -0.05, -0.048, -0.0424, -0.0376, -0.034,
    -0.085, -0.082, -0.0814, -0.0806, -0.079, -0.078, -0.0764, -0.074, -0.073, -0.068, -0.062, -0.06, -0.053, -0.051, -0.044, -0.0387, -0.0345,
    -0.0748, -0.074, -0.073, -0.073, -0.073, -0.0729, -0.0713, -0.0696, -0.0686, -0.0678, -0.066, -0.065, -0.056, -0.0531, -0.046, -0.04, -0.035,
    -0.0665, -0.0668, -0.067, -0.067, -0.0672, -0.067, -0.0663, -0.065, -0.0647, -0.064, -0.063, -0.0627, -0.061, -0.06, -0.0497, -0.0425, -0.037,
    -0.0575, -0.0585, -0.06, -0.0605, -0.0614, -0.0616, -0.0626, -0.062, -0.0615, -0.061, -0.06, -0.0599, -0.0582, -0.058, -0.0558, -0.0454, -0.04,
    -0.0985, -0.116, -0.135, -0.1285, -0.1185, -0.113, -0.098, -0.0865, -0.081, -0.077, -0.072, -0.0695, -0.0635, -0.061, -0.054, -0.049, -0.0445,
    -0.0985, -0.1115, -0.125, -0.123, -0.1165, -0.11, -0.101, -0.0895, -0.0833, -0.08, -0.073, -0.071, -0.0645, -0.062, -0.055, -0.0495, -0.045,
    -0.0985, -0.1058, -0.113, -0.112, -0.1108, -0.11, -0.101, -0.0905, -0.085, -0.0815, -0.075, -0.0728, -0.0655, -0.063, -0.056, -0.05, -0.0455,
    -0.092, -0.097, -0.1025, -0.1028, -0.102, -0.101, -0.098, -0.0925, -0.088, -0.0845, -0.079, -0.076, -0.069, -0.0655, -0.058, -0.0516, -0.0465,
    -0.084, -0.0871, -0.0905, -0.092, -0.0936, -0.0931, -0.0913, -0.089, -0.088, -0.087, -0.082, -0.0795, -0.0725, -0.069, -0.061, -0.0542, -0.049,
    -0.0765, -0.0776, -0.0788, -0.081, -0.083, -0.084, -0.086, -0.0843, -0.0835, -0.083, -0.0814, -0.0806, -0.076, -0.0731, -0.0644, -0.057, -0.0509,
    -0.069, -0.0696, -0.07, -0.072, -0.0739, -0.0746, -0.077, -0.0786, -0.0789, -0.0783, -0.0773, -0.077, -0.075, -0.0735, -0.0684, -0.0606, -0.0536,
    -0.098, -0.119, -0.1375, -0.131, -0.123, -0.118, -0.1055, -0.095, -0.09, -0.086, -0.0815, -0.0777, -0.0712, -0.0682, -0.0619, -0.056, -0.0512,
    -0.098, -0.114, -0.131, -0.127, -0.121, -0.1175, -0.1055, -0.0977, -0.0938, -0.086, -0.0815, -0.0782, -0.0721, -0.069, -0.062, -0.056, -0.0512,
    -0.098, -0.1088, -0.119, -0.1182, -0.1155, -0.1133, -0.1055, -0.0977, -0.094, -0.089, -0.0842, -0.081, -0.073, -0.0692, -0.062, -0.056, -0.0512,
    -0.095, -0.1011, -0.1075, -0.1091, -0.1081, -0.1078, -0.104, -0.0977, -0.095, -0.0905, -0.086, -0.0831, -0.076, -0.073, -0.0645, -0.0576, -0.052,
    -0.0875, -0.091, -0.094, -0.0958, -0.0982, -0.0996, -0.0985, -0.0977, -0.095, -0.0915, -0.0879, -0.0855, -0.0791, -0.076, -0.0675, -0.0605, -0.0548,
    -0.0816, -0.0822, -0.0826, -0.0843, -0.0865, -0.0875, -0.0905, -0.0906, -0.0899, -0.0886, -0.0879, -0.0855, -0.0806, -0.0778, -0.07, -0.063, -0.056,
    -0.0738, -0.073, -0.0725, -0.074, -0.076, -0.077, -0.08, -0.0827, -0.084, -0.085, -0.0845, -0.0841, -0.083, -0.081, -0.073, -0.0657, -0.059,
    -0.0989, -0.12, -0.1415, -0.1335, -0.124, -0.119, -0.108, -0.099, -0.094, -0.09, -0.0861, -0.083, -0.0775, -0.074, -0.066, -0.06, -0.0545,
    -0.0989, -0.1145, -0.1311, -0.1265, -0.12, -0.117, -0.108, -0.099, -0.094, -0.09, -0.0861, -0.083, -0.0775, -0.074, -0.066, -0.06, -0.0545,
    -0.0989, -0.108, -0.1171, -0.121, -0.117, -0.1144, -0.107, -0.1, -0.096, -0.093, -0.0882, -0.085, -0.0795, -0.075, -0.067, -0.0605, -0.055,
    -0.0959, -0.0999, -0.1038, -0.1056, -0.1085, -0.111, -0.1054, -0.1, -0.097, -0.0935, -0.09, -0.0865, -0.0816, -0.0775, -0.069, -0.0623, -0.0565,
    -0.0899, -0.0906, -0.0915, -0.096, -0.098, -0.099, -0.101, -0.1, -0.0972, -0.095, -0.092, -0.089, -0.0841, -0.0806, -0.0725, -0.065, -0.0585,
    -0.083, -0.082, -0.081, -0.0845, -0.088, -0.0898, -0.093, -0.094, -0.0935, -0.0934, -0.092, -0.09, -0.086, -0.0825, -0.0745, -0.067, -0.06,
    -0.0765, -0.0746, -0.073, -0.076, -0.0798, -0.0815, -0.0846, -0.086, -0.0867, -0.0875, -0.0875, -0.0865, -0.0841, -0.0825, -0.0755, -0.0686, -0.062,
    -0.098, -0.1135, -0.1283, -0.13, -0.1245, -0.12, -0.109, -0.1, -0.096, -0.0915, -0.0881, -0.086, -0.0794, -0.076, -0.0681, -0.062, -0.057,
    -0.098, -0.112, -0.1255, -0.1265, -0.123, -0.119, -0.109, -0.1, -0.096, -0.0915, -0.0881, -0.086, -0.0794, -0.076, -0.0681, -0.062, -0.057,
    -0.098, -0.1051, -0.1125, -0.118, -0.1178, -0.116, -0.109, -0.102, -0.0985, -0.0935, -0.09, -0.0868, -0.0815, -0.0775, -0.0695, -0.063, -0.0578,
    -0.1, -0.091, -0.104, -0.1083, -0.114, -0.1128, -0.1075, -0.102, -0.0985, -0.095, -0.092, -0.089, -0.0831, -0.0795, -0.072, -0.065, -0.059,
    -0.09, -0.0912, -0.092, -0.0955, -0.0994, -0.101, -0.1045, -0.102, -0.098, -0.095, -0.093, -0.09, -0.085, -0.0811, -0.073, -0.066, -0.0606,
    -0.0832, -0.0823, -0.0815, -0.085, -0.0889, -0.091, -0.094, -0.0956, -0.096, -0.096, -0.093, -0.0905, -0.086, -0.083, -0.075, -0.0682, -0.0625,
    -0.076, -0.0743, -0.0725, -0.0751, -0.0787, -0.081, -0.084, -0.087, -0.088, -0.089, -0.0895, -0.09, -0.0867, -0.084, -0.0767, -0.0705, -0.065,
])


def _tl(x1, x2, flat, q1, q2, l1, l2, u1, u2):
    grid = np.asarray(flat, dtype=float).reshape(len(x1), len(x2)).T
    return float(tlinex(x1, x2, grid, q1, q2, l1, l2, u1, u2))


def calculate_suplat(data: Mapping[str, object], tail: bool = False
                     ) -> Dict[str, object]:
    """Translate SUPLAT (or SUPLAH, ``tail=True``): supersonic sideslip
    derivatives of the wing (tail) and the wing-body (tail-body).

    Args:
        data: By name: ``alpha`` (``FLC(23..)``), ``mach``, ``rl``
            (``FLC(I+42)``), ``i`` (the Mach index), ``sref``, ``cbarr``,
            ``blref``, ``bo``, ``htpl``, ``syna`` (``/SYNTSS/`` 1-8),
            ``win`` (the surface's ``WINGIN``/``HTIN`` words 3, 4, 6, 12,
            13, 14, 15 (1.0 straight)), ``a`` (its ``A`` words 10, 59, 61,
            62, 74, 118, 120), ``clpcty``, ``cnaw`` (its ``SLG(3)``,
            ``(7)``), ``cn`` (its ``WING(61..)``), ``clab`` (``SBD(18)``),
            ``bd`` (``BD(66)`` for the wing, ``BD(89)`` for the tail),
            ``body`` (``/BODYI/``: ``x``, ``r``), and for SUPLAT's tail
            section ``htin`` (3, 4, and ``SHB``, ``SEXT``, ``RLPH`` at the
            Mach index), ``aht`` (30, 62), and ``stale``: the surface's
            ``cybw``, ``cnbw``, ``clbw`` curves and ``bwi`` words where
            the routine does not set them.

    Returns:
        ``sla`` (the 23 ``SBETA`` words by index), ``win`` (the dihedral
        words UNUSED set to 0), ``cybw``, ``cnbw``, ``clbw`` (the
        surface's ``(141..)``, ``(161..)``, ``(181..)``), ``bwi`` (the
        body block's 141, 161, 181..), ``rlb`` (``BD(1)``), and for SUPLAT
        ``bwh``, ``bwv``, ``bwhv`` words and ``lf``.

    Notes:
        Kept as executed: SUPLAH forms ``(TANLE/4)**1.3333`` without
        SUPLAT's ``ABS``, which is NaN for a swept-forward leading edge,
        and takes the delta branch at ``beta/tan(LE) <= 1`` where SUPLAT
        uses ``< 0.998``.
    """
    st = dict(data.get('stale', {}))
    w = {int(k): float(v) for k, v in data['win'].items()}
    g = {int(k): float(v) for k, v in data['a'].items()}
    syna = {int(k): float(v) for k, v in data['syna'].items()}
    alpha = [float(v) for v in data['alpha']]
    na = len(alpha)
    i = int(data['i'])
    sr, crbar, blref = (float(data['sref']), float(data['cbarr']),
                        float(data['blref']))
    mach = float(data['mach'])
    beta = math.sqrt(mach**2 - 1.)
    sla = {1: mach, 2: beta}
    out: Dict[str, object] = {}
    if not tail:
        out['lf'] = 0
    scale = 2. * w[4] / blref
    bx = [float(v) for v in data['body']['x']]
    br = [float(v) for v in data['body']['r']]
    rlb = bx[-1]
    rm, rb = mach**2, beta**2
    for k in (13, 12, 14):
        if w[k] == UNUSED:
            w[k] = 0.0
    diheq = (w[13] * (w[3] - w[12]) + w[14] * w[12]) / w[3]
    sla[4] = diheq
    cybw = list(st.get('cybw', [0.0] * na))
    cnbw = list(st.get('cnbw', [0.0] * na))
    clbw = list(st.get('clbw', [0.0] * na))
    swepe, tapr, ar, tanle = g[59], g[118], g[120], g[62]
    xw = syna[6] if tail else syna[2]
    xcg, cr, span, spano = syna[1], w[6], w[4], w[3]
    straight = w[15] == 1.0
    if straight:
        if swepe == 0 and tapr == 1 and ar * beta >= 1.0:
            x = xw + 0.5 * cr - xcg
            sla[3] = x
            arg1 = 1. / (PI * ar**2 * rb)
            arg2 = 4. * rm / 3. + 8. * rm * x / crbar
            arg3 = PI * ar * (1. - rb) * (3. + rb) / (3. * beta**3)
            for j in range(na):
                alp = (alpha[j] / RAD)**2
                cybw[j] = -alp * 8.0 * rm / (RAD * PI * rb * ar) - \
                    .0001 * abs(diheq)
                cnbw[j] = scale * alp * arg1 * (arg2 - arg3) / RAD
        else:
            arg = beta / tanle if tanle != 0.0 else beta / .00001
            delta_branch = (tapr == 0.0 and
                            (arg <= 1.0 if tail else arg < .998))
            if delta_branch:
                # At beta/tan(LE) = 1 the chart gives 1/Q = 0 and the
                # source divides by it: IEEE semantics, as compiled.
                qbc = np.float64(interx(1, _X51116, [arg], [16], _Y51116,
                                        lind=16))
                ebc = interx(1, _X71118, [arg], [8], _Y71118, lind=8)
                sla[5], sla[6] = float(qbc), ebc
                with np.errstate(all='ignore'):
                    arg1 = PI * ar * rm / (4. * RAD * qbc)
                x = xw + .6666 * span * tanle - xcg
                sla[3] = x
                with np.errstate(all='ignore'):
                    arg2 = PI / 3. * (ebc + (ar**2 / 16. + x / crbar) * rm /
                                      qbc) / RAD
                    for j in range(na):
                        alp = (alpha[j] / RAD)**2
                        cybw[j] = float(-alp * arg1 - .0001 * abs(diheq))
                        cnbw[j] = float(scale * alp * arg2)
        clptoa = interx(3, _X12225, [beta * ar, ar * g[74], tapr],
                        [17, 7, 5], _Y12225, lind=17, lx1u=2, lx2u=2)
        clp = clptoa * ar * float(data['clpcty'])
        clbd = 2. * diheq * (1. + 2. * tapr) * clp / (RAD**2. *
                                                     (1. + 3. * tapr))
        sla.update({7: clptoa, 8: clp, 9: clbd})
        arg1 = (1. + tapr * (1. + swepe)) * (1. + swepe / 2.) * tanle / beta
        with np.errstate(invalid='ignore'):
            tpow = (float(np.float64(tanle / 4.) ** 1.3333) if tail else
                    abs(tanle / 4.)**1.3333)
        arg2 = rm * g[61]**2 / ar + tpow
        cnw = [float(v) for v in data['cn']]
        cnaw = float(data['cnaw'])
        clbw = [scale * (clbd - 0.061 * cnw[j] * cnaw * arg1 * arg2 / RAD)
                for j in range(na)]
    out.update({'cybw': cybw, 'cnbw': cnbw, 'clbw': clbw, 'win': w})
    if not data['bo']:
        out['sla'] = sla
        return out
    bwi = {int(k): float(v) for k, v in st.get('bwi', {}).items()}
    spanin = span - spano
    inc = syna[8] if tail else syna[4]
    height = syna[7] if tail else syna[3]
    zw = -height + (.25 * g[10] + float(data['bd'])) * math.sin(inc / RAD) - \
        spanin * math.tan(w[13] / RAD)
    arg = zw / spanin
    rki = 1.0 - .85 * arg if arg < 0.0 else 1.0 + .49 * arg
    clab = float(data['clab'])
    bwi[141] = -rki * clab - 0.0001 * abs(diheq)
    rnn = float(data['rl']) * rlb
    rkrl = 1. + math.log(1.e-6 * rnn) / 4.86
    rh1 = 2. * tbfunx(bx, br, rlb * .25, 0, 0)[0]
    rh2 = 2. * tbfunx(bx, br, rlb * .75, 0, 0)[0]
    sbs = 2. * float(trapz(br, bx)[0])
    ydumy = _tl(_X158A, _X258A, _Y58A, rlb**2 / sbs, xcg / rlb, 2, 1, 2, 1)
    ydumy2 = _tl(_X158B, _X258B, _Y58B, math.sqrt(rh1 / rh2), ydumy,
                 2, 0, 2, 1)
    rkn = _tl(_X158C, _X258C, _Y58C, 1.0, ydumy2, 2, 0, 2, 1)
    bwi[161] = -rkn * rkrl * sbs * rlb / (sr * blref)
    sla.update({10: zw, 11: rki, 12: rnn, 13: rkrl, 14: rh1, 15: rh2,
                16: sbs, 17: rkn})
    bwv, bwh, bwhv = {}, {}, {}
    if straight:
        zwp = zw + cr * math.sin((syna[8] if tail else syna[4]) / RAD) / 4.0
        arg1 = math.sqrt(ar)
        dbody = 2.0 * spanin
        clbzw = 0.6 * arg1 * zwp * dbody / (RAD * span**2)
        dclb = -0.0005 * arg1 * (dbody / (2. * span))**2 * diheq
        sla.update({18: zwp, 19: clbzw, 20: dclb})
        for j in range(na):
            bwi[181 + j] = clbw[j] + (clbzw + dclb) * scale
            if not tail:
                bwv[181 + j] = bwi[181 + j]
        if not tail:
            bwv[141], bwv[161] = bwi[141], bwi[161]
    out.update({'sla': sla, 'bwi': bwi, 'rlb': rlb})
    if tail:
        return out
    h = {int(k): float(v) for k, v in data['htin'].items()}
    aht = {int(k): float(v) for k, v in data['aht'].items()}
    zh = syna[7] - ((h[4] - h[3]) * aht[62] + aht[30]) * \
        math.sin(syna[8] / RAD)
    htpl = bool(data['htpl'])
    dcyhwb = 0.0
    rh = h[4] - h[3] if htpl else float(st.get('rh', 0.0))
    if htpl and not (abs(zh) > rh or zh / rh == 0.0):
        arg = rh / h[4]
        if arg < 1.0:
            rkhbhl = interx(1, _XB25OO, [arg], [8], _YB25OO, lind=8)
        else:
            rkhbhl = interx(1, _XA25OO, [1. / arg], [6], _YA25OO, lind=6)
        rkhb = rkhbhl * (1. - (1. - (zh / rh)**2)**0.5)
        dcyhwb = -rkhb * clab * h[114 + i] / h[134 + i]
        sla.update({21: rkhbhl, 22: rkhb})
    elif htpl and zh / rh != 0.0:
        out['lf'] = 1
    sla[23] = dcyhwb
    if htpl:
        bwh[141] = bwi[141] + dcyhwb
        bwh[161] = bwi[161] - dcyhwb * h[94 + i] / blref
        bwhv[141], bwhv[161] = bwh[141], bwh[161]
        if straight:
            for j in range(na):
                bwh[181 + j] = bwi[181 + j]
                bwv[181 + j] = bwi[181 + j]
                bwhv[181 + j] = bwh[181 + j]
    out.update({'bwh': bwh, 'bwv': bwv, 'bwhv': bwhv, 'zh': zh})
    return out


def m23o27_words(nalpha: int, body: Mapping[int, float],
                 vt: Mapping[int, float], vf: Mapping[int, float]
                 ) -> Dict[str, Dict[int, float]]:
    """M23O27's pass after SUPLAT, SUPLAH, SUPLAV and SUPLAF: the angle-
    dependent side-force and yawing-moment words past the first angle are
    set to -UNUSED in every block, and the body-vertical block ``BV``
    sums the body and both vertical panels.

    Args:
        nalpha: ``FLC(2)+0.5``.
        body, vt, vf: The ``BODY``, ``VT`` and ``VF`` words 141, 161 and
            181.. .

    Returns:
        The words set, by block: ``surface`` (the value written to
        ``WING``, ``HT``, ``VT``, ``VF``, ``BW``, ``BWH``, ``BWV``,
        ``BWHV`` and ``BH`` words 142.., 162..) and ``bv``.
    """
    fill = {}
    for j in range(2, nalpha + 1):
        fill[j + 140] = -UNUSED
        fill[j + 160] = -UNUSED
    b = {int(k): float(v) for k, v in body.items()}
    t = {int(k): float(v) for k, v in vt.items()}
    f = {int(k): float(v) for k, v in vf.items()}
    bv = {k: b[k] + t[k] + f[k] for k in (141, 161, 181)}
    for j in range(2, nalpha + 1):
        bv[j + 140] = -UNUSED
        bv[j + 160] = -UNUSED
        bv[j + 180] = b[j + 180] + t[j + 180] + f[j + 180]
    return {'surface': fill, 'bv': bv}
