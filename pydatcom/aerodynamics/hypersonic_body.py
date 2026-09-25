"""
Hypersonic body: HYPBOD and its overlay M26O32.

HYPBOD finds the nose, centre-section and tail diameters from the body
stations, the small-angle normal-force and moment slopes of each section
(Figures 4.2.1.1-24 and 4.2.2.1-27 for an ogive nose, the cone equations
for a cone nose, the Figure 4.2.1.1-26 and 4.2.2.1-25A equations for the
flares), and the Newtonian normal force, axial force and moment at each
angle by NASA TN D-176's integration along the body.  M26O32 adds the
slopes of the lift and moment curves and the stability-axis words.

:func:`hypersonic._hypbod_coefficients` is the angle loop alone, driven by
a state dictionary; this module is the whole routine on its COMMON words.

The tables were extracted from the source DATA statements by parsing, with
``tools/fortran_data.py``.  Figure 4.2.3.1-68 (``T468``/``D468``) and the
Figure 4.1.5.1-27 grid are declared but never read.

Reference: datcom-legacy/datcom_2000/hypbod.f, m26o32.f
"""

import math
from typing import Dict, Mapping, Sequence

import numpy as np

from pydatcom.utils.constants import PI, RAD
from pydatcom.utils.legacy_interp import interx
from pydatcom.utils.legacy_numeric import tbfunx, trapz

# Figure 4.2.1.1-24: the nose normal-force slope, against l_N/d_N (8,
# padded to 9) and d_s/d_N (9).
_T425 = np.array([
    1.0, 1.5, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 0.0,
    0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8,
])
_D425 = np.array([
    1.6, 1.766, 1.85, 1.93, 1.96, 1.976, 1.987, 1.99,
    1.6, 1.764, 1.85, 1.921, 1.95, 1.968, 1.975, 1.979,
    1.592, 1.758, 1.84, 1.904, 1.93, 1.942, 1.95, 1.95,
    1.581, 1.74, 1.81, 1.864, 1.889, 1.9, 1.901, 1.9,
    1.56, 1.704, 1.761, 1.809, 1.822, 1.83, 1.835, 1.838,
    1.53, 1.65, 1.693, 1.729, 1.739, 1.743, 1.75, 1.75,
    1.488, 1.565, 1.604, 1.629, 1.634, 1.64, 1.64, 1.642,
    1.419, 1.463, 1.489, 1.506, 1.51, 1.513, 1.516, 1.515,
    1.32, 1.332, 1.342, 1.354, 1.36, 1.365, 1.369, 1.37,
])
# Figure 4.2.2.1-27: the nose centre of pressure, x_cp/l_N, against
# l_N/d_N (9) and d_s/d_N (9).
_T422 = np.array([
    1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0,
    0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8,
])
_D422 = np.array([
    0.496, 0.47, 0.468, 0.468, 0.467, 0.467, 0.467, 0.467, 0.467,
    0.486, 0.457, 0.453, 0.451, 0.45, 0.448, 0.447, 0.446, 0.445,
    0.475, 0.441, 0.436, 0.433, 0.43, 0.428, 0.426, 0.425, 0.424,
    0.462, 0.426, 0.419, 0.414, 0.41, 0.407, 0.405, 0.403, 0.401,
    0.45, 0.408, 0.4, 0.394, 0.389, 0.385, 0.382, 0.38, 0.377,
    0.438, 0.39, 0.379, 0.371, 0.365, 0.36, 0.356, 0.353, 0.349,
    0.423, 0.367, 0.353, 0.344, 0.336, 0.329, 0.324, 0.32, 0.317,
    0.405, 0.337, 0.317, 0.306, 0.297, 0.291, 0.286, 0.282, 0.278,
    0.381, 0.296, 0.274, 0.26, 0.249, 0.242, 0.236, 0.231, 0.227,
])


def calculate_hypbod(alpha_deg: Sequence[float], mach: float,
                     body: Mapping[str, object], xcg: float, sref: float,
                     cbar: float, stale: Mapping[str, float]
                     ) -> Dict[str, object]:
    """Translate HYPBOD: hypersonic body aerodynamics.

    Args:
        alpha_deg: ``FLC(23..)``; mach: ``FLC(NMN+2)``.
        body: The ``/BODYI/`` words ``x``, ``r`` (stations and radii),
            ``rln``, ``rla`` (nose and centre-section lengths), ``ds``
            (nose bluntness diameter) and ``bnose`` (1 for a cone nose).
        xcg: ``XCG``; sref, cbar: ``/OPTION/``.
        stale: ``SBD`` words read without being set on some paths:
            ``thetaa`` (130, no centre section) and ``thetat`` (135, no
            tail, stored into ``SBD(14)``).

    Returns:
        The ``SBD`` words set, keyed by index (1-6, 10, 14, 16, 18, 102,
        110, 125-140, and the curves ``theta`` 141.., ``lx`` 161..,
        ``intgcn`` 181.., ``intgcm`` 201..), and the ``BODY`` curves
        ``cd``, ``cl``, ``cm``, ``cn``, ``ca``.

    Notes:
        Kept as executed: the centre-section moment ``CMAA`` omits the
        ``/RAD`` that the nose and tail moments carry, so it is per
        radian where they are per degree; and the tail area ``ATT`` is
        formed on the centre-section diameter ``D1``.
    """
    x = np.asarray(body['x'], dtype=float)
    r = np.asarray(body['r'], dtype=float)
    nx = len(x)
    rln, rla = float(body['rln']), float(body['rla'])
    ds = float(body['ds'])
    rlbp = rln + rla
    rlb = float(x[-1])
    rlbt = rlb - rlbp
    tail = abs(rlbt / rlb) >= 0.01
    if not tail:
        rlbt = 0.0

    def diameter(station):
        return 2. * interx(1, x, [station], [nx], r, lind=nx)

    dn = diameter(rln)
    d1 = diameter(rlbp) if rla != 0.0 else dn
    d2 = diameter(rlb) if tail else d1
    fn = rln / dn
    thetan = math.atan(.5 * (dn - ds) / (dn * fn - .5 * ds))
    if float(body['bnose']) != 1.:
        var = [rln / dn, ds / dn]
        cnanf = interx(2, _T425, var, [8, 9], _D425, lind=9)
        xcpln = interx(2, _T422, var, [9, 9], _D422, lind=9)
    else:
        c2 = math.cos(thetan)**2
        cnanf = 2. * c2 * (1. - .5 * (ds / dn)**2 * c2)
        arg = ds * math.cos(thetan) / dn
        sint = math.sin(thetan)
        xcpln = (2. / (3. * c2)) * (
            ((1. - 1.5 * (1. - sint) * arg) + .25 * (2. - 3. * sint) *
             arg**3) / ((1. - .5 * arg**2) * (1. - arg / (1. + sint))))
    ann = PI * dn**2 / 4.
    cnan = cnanf * ann / (RAD * sref)
    cman = cnan * (xcg / rln - xcpln) * rln / cbar
    aaa = PI * d1**2 / 4.
    thetaa = float(stale['thetaa'])
    cnaaf = cmaaf = 0.0
    if rla != 0.0:
        ratioa = dn / d1
        thetaa = math.atan(.5 * (d1 - dn) / rla)
        if ratioa < 1.0:
            cnaaf = 2. * math.cos(thetaa)**2 * (1. - ratioa**2)
            cmaaf = -((2. * (1. - ratioa**3)) -
                      3. * ratioa * math.cos(thetaa)**2 *
                      (1. - ratioa**2)) / (3. * math.tan(thetaa))
    cnaa = cnaaf * aaa / (RAD * sref)
    cmaa = cmaaf * aaa * d1 / (cbar * sref) + cnaa * (xcg - rln) / cbar
    att = PI * d1**2 / 4.
    sbd: Dict[int, object] = {}
    thetat = float(stale['thetat'])
    cnatf = cmatf = 0.0
    flare = False
    if tail:
        ratiot = d2 / d1
        thetat = math.atan(.5 * (d2 - d1) / rlbt)
        flare = ratiot > 1.0
        if flare:
            sbd[16] = sbd[102] = thetat
            cnatf = 2. * math.cos(thetat)**2 * (1. - ratiot**2)
            cmatf = -((2. * (1. - ratiot**3)) -
                      3. * ratiot * math.cos(thetat)**2 *
                      (1. - ratiot**2)) / (3. * math.tan(thetat))
    if not flare:
        sbd[14] = thetat
    cnat = cnatf * att / (RAD * sref)
    cmat = cmatf * att * d2 / (cbar * RAD * sref) + \
        cnat * (xcg - rln - rla) / cbar

    theta = [RAD * math.atan(tbfunx(x, r, xn, -1, -1)[1]) for xn in x[:-1]]
    theta.append(theta[-1])
    lx = [xcg - xn for xn in x]
    k = 1.833 * (1. - 0.4545 / mach**2)
    factor = k * rlb / sref
    indep = x / rlb
    intgcn, intgcm, cn, cm, caf, cl, cd = ([] for _ in range(7))
    for angle_deg in alpha_deg:
        a = abs(angle_deg / RAD)
        sa, ca, ta = math.sin(a), math.cos(a), math.tan(a)
        dep, dep1, dep2 = [], [], []
        for station in range(nx):
            th = theta[station] / RAD
            tn, cnn, sn = math.tan(th), math.cos(th), math.sin(th)
            if abs(angle_deg / RAD) > abs(th):
                phe = math.acos(tn / ta)
            else:
                phe = 0.0 if th > 0.0 else PI
            sp, cp = math.sin(phe), math.cos(phe)
            ktheta = ((2. / 3.) * (cnn * sa)**2 * sp * (cp**2 + 2.) +
                      4. * sn * cnn * ca * sa *
                      (PI / 2. - .5 * sp * cp - phe / 2.) +
                      2. * (sn * ca)**2 * sp)
            kaf = (2. * (ca * sn)**2 * tn * (PI - phe) +
                   4. * ca * sa * sp * sn**2 +
                   cnn * sn * sa**2 * (PI - phe - sp * cp))
            dep.append(ktheta * r[station])
            dep1.append(dep[station] * lx[station])
            dep2.append(kaf * r[station])
        intgcn.append(float(trapz(dep, indep)[0]))
        intgcm.append(float(trapz(dep1, indep)[0]))
        intgca = float(trapz(dep2, indep)[0])
        sign = -1.0 if angle_deg < 0.0 else 1.0
        cn.append(sign * factor * intgcn[-1])
        cm.append(sign * factor * intgcm[-1] / cbar)
        caf.append(factor * intgca)
    for angle_index, angle_deg in enumerate(alpha_deg):
        arg = angle_deg / RAD
        cl.append(cn[angle_index] * math.cos(arg) -
                  caf[angle_index] * math.sin(arg))
        cd.append(caf[angle_index] * math.cos(arg) +
                  cn[angle_index] * math.sin(arg))
    sbd.update({
        1: rlbp, 2: rlb, 3: rlbt, 4: dn, 5: d1, 6: d2, 10: fn,
        18: cnan + cnaa + cnat, 110: cman + cmaa + cmat, 125: cnanf,
        126: xcpln, 127: thetan, 128: cnan, 129: cman, 130: thetaa,
        131: cnaaf, 132: cmaaf, 133: cnaa, 134: cmaa, 135: thetat,
        136: cnatf, 137: cmatf, 138: cnat, 139: cmat, 140: k})
    return {'sbd': sbd, 'theta': theta, 'lx': lx, 'intgcn': intgcn,
            'intgcm': intgcm, 'cn': cn, 'cm': cm, 'ca': caf, 'cl': cl,
            'cd': cd, 'method': 'legacy_hypbod'}


def m26o32_slopes(alpha_deg: Sequence[float], cl: Sequence[float],
                  cm: Sequence[float], cbarr: float, blref: float
                  ) -> Dict[str, list]:
    """M26O32's pass after HYPBOD: ``BODY(101..)`` and ``(121..)`` take
    TBFUNX's slopes of the lift and moment curves, ``(141..)`` and
    ``(161..)`` their negatives (the moment scaled by ``CBARR/BLREF``),
    and ``(181..)`` zero."""
    alpha = np.asarray(alpha_deg, dtype=float)
    cla = [tbfunx(alpha, cl, a, 0, 0)[1] for a in alpha]
    cma = [tbfunx(alpha, cm, a, 0, 0)[1] for a in alpha]
    bfact = cbarr / blref
    return {'cla': cla, 'cma': cma, 'b141': [-v for v in cla],
            'b161': [-bfact * v for v in cma], 'b181': [0.0] * len(cla)}
