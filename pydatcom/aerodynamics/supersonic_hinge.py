"""
Supersonic hinge moments: SSHING.

SSHING forms the supersonic hinge-moment derivatives of a trailing-edge
control (Section 6.1.6): the hinge-line geometry, the Mach-cone regions
at the control's tip and root edges (PTCP's pressures over up to eight
regions), the flat-sided and biconvex ``Ch_alpha`` and ``Ch_delta`` with
DFLCON's control derivatives, and the area moment about the hinge.  For a
translating flap it forms the chord-extension factor and lift-curve slope
instead.  SSSYM follows in the source and is translated separately.

Reference: datcom-legacy/datcom_2000/sshing.f
"""

import math
from typing import Dict, List, Mapping

from pydatcom.aerodynamics.ptcp import calculate_ptcp
from pydatcom.aerodynamics.supersonic_control import calculate_dflcon
from pydatcom.utils.constants import PI, RAD
from pydatcom.utils.math_utils import arcsin

_ITRANS = (2, 3, 4, 7, 8)


def _alpha(num, den):
    if den == 0.:
        return PI / 2.
    v = math.atan(num / den)
    return PI + v if v < 0. else v


def calculate_sshing(data: Mapping[str, object]) -> Dict[str, object]:
    """Translate SSHING: supersonic hinge-moment derivatives.

    Args:
        data: By name: ``htpl``; ``a``, ``aht`` (the wing's and tail's
            ``A`` words 10, 25, 28, 59, 62, 77, 80, 86, 104); ``win``,
            ``htin`` (``WINGIN``/``HTIN`` 1-6, 18, 66); ``claw``,
            ``clah`` (``WING(101)``, ``HT(101)``); ``f`` (``/FLAPIN/``
            11-17, 39.., 49..); ``mach``; ``sref``; ``spr`` (``SPR(1..59)``
            before the call); and ``stale``: the saved locals ``lamhl``,
            ``k``, ``mu``, ``xt``, ``xr``, ``ci``, ``ct``,
            ``tovca`` and ``kaseno``, read when the routine does not set
            them first.

    Returns:
        ``returned`` (an edge or the hinge line is subsonic: nothing after
        the hinge-line block is set), ``spr`` (the 59 words as left),
        ``clalds`` (``WING(241..)``, translating flaps), ``chabc``
        (``WING(251)``), ``chdelr`` (``WING(261)``), and ``kaseno``.

    Notes:
        Kept as executed: the hinge-moment sum ``SPAMT+SPAMR`` is taken
        after the loop, when ``SPAMT`` holds the root's value, so the tip's
        area moment is dropped and the root's counted twice; with an
        unswept leading edge the root's region moments ``SPR(45..52)``
        repeat the tip's; the tail branch takes ``A(77)``, the wing's
        trailing-edge sweep, and PTCP reads the wing's leading- and
        trailing-edge tangents; a flap type other than 1 skips the
        hinge-line block and reads the previous call's geometry.
    """
    htpl = bool(data['htpl'])
    a = {int(k): float(v) for k, v in data['a'].items()}
    at = {int(k): float(v) for k, v in data['aht'].items()}
    f = {int(k): float(v) for k, v in data['f'].items()}
    spr = [0.0] + [float(v) for v in data['spr']]
    st = {k: float(v) for k, v in data['stale'].items()}
    ct = st['ct']
    if not htpl:
        w = {int(k): float(v) for k, v in data['win'].items()}
        g = a
        sai100 = a[77]
        ct = w[1]
        claw = float(data['claw'])
    else:
        w = {int(k): float(v) for k, v in data['htin'].items()}
        g = at
        sai100 = a[77]
        claw = float(data['clah'])
    sspne, sai000, chrdre, chrdtp = w[3], g[59], g[10], w[1]
    tanle, tante, bo2, cr = g[62], g[80], w[4], w[6]
    tapri, bstro2, tanteo, tanleo = g[25], w[2], g[104], g[86]
    cb, tapro, tovc, tovco = w[5], g[28], w[18], w[66]
    sr = float(data['sref'])
    mach = float(data['mach'])
    aloci, aloco, cfi, cfo = f[14], f[15], f[12], f[13]
    beta = math.sqrt(mach**2 - 1.)
    spr[1] = beta
    ndelta = int(f[16] + .5)
    iftype = int(f[17] + .5)
    lamhl, k_, mu = st['lamhl'], st['k'], st['mu']
    xt, xr, ci = st['xt'], st['xr'], st['ci']
    co = spr[36]                    # CO is SPR(36), not a saved local
    out: Dict[str, object] = {'returned': False}
    if iftype == 1:
        lamhl = math.atan(((aloco - aloci) * math.tan(sai100) + cfi - cfo) /
                          (aloco - aloci))
        spr[4] = lamhl * RAD
        spr[5] = 2 * math.atan(f[11] / math.cos(lamhl)) * RAD
        spr[14] = math.tan(lamhl)
        spr[32] = tanle / beta
        spr[33] = spr[14] / beta
        spr[34] = tante / beta
        if spr[32] >= 1. or spr[33] >= 1. or spr[34] >= 1.:
            out.update({'returned': True, 'spr': spr[1:], 'clalds': None,
                        'kaseno': None,
                        'locals': dict(st, lamhl=lamhl, ct=ct)})
            return out
        spr[35] = cfo / cfi
        co = cr + aloco * (tante - tanle)
        ci = cr + aloci * (tante - tanle)
        spr[36] = co
        spr[31] = ci
        if aloci > bo2 - bstro2:
            spr[31] = ct + (bo2 - aloci) * (tanleo - tanteo)
        mu = arcsin(1. / mach)
        spr[2] = 2. / beta
        spr[3] = (2.4 * mach**4 - 4. * beta**2) / (2. * beta**4)
        spr[6] = 1. - (spr[3] / spr[2]) * (spr[5] / RAD)
        k_ = math.tan(sai000 - lamhl) * math.tan(sai000 - sai100)
        xt = co - cfo + (sspne - aloco) * (spr[14] - tanle)
        xr = ci + aloci * (tanle - spr[14]) - cfi
    transl = iftype in _ITRANS
    clalds = None
    if transl:
        deln4 = 0.25 * (aloco - aloci) / bo2
        eta = [aloci / bo2]
        arg2 = (tante - tanle) * bo2
        arg3, arg4, arg5 = bo2 * cr, cr, tapri
        if not aloci < bo2 - bstro2:
            arg2 = (tanteo - tanleo) * bo2
            arg3, arg4, arg5 = bo2 * cr, cb, tapro
        chrd = [arg4 + eta[0] * arg2]
        spr[31] = chrd[0]
        swf = []
        for m in range(1, 5):
            eta.append(eta[m - 1] + deln4)
            chrd.append(arg4 + eta[m] * arg2)
            swf.append(arg3 * (2. - (1. - arg5) * (eta[m - 1] + eta[m])))
        clalds = []
        for i in range(ndelta):
            cpi, cpo = f[39 + i], f[49 + i]
            cp = [cpi] + [cpi - (cpi - cpo) / (4. * deln4) *
                          (eta[m] - eta[0]) for m in range(1, 5)]
            cfactr = [((cp[m] / chrd[m] + cp[m - 1] / chrd[m - 1]) / 2. -
                       1.) * swf[m - 1] / sr for m in range(1, 5)]
            cfact = (cfactr[0] + cfactr[1] + cfactr[2] + cfactr[3]) / 4.
            clalds.append(claw * (1. + cfact) + claw)
        out.update({'spr': spr[1:], 'clalds': clalds, 'kaseno': None,
                    'locals': dict(st, lamhl=lamhl, k=k_, mu=mu, xt=xt,
                                   xr=xr, co=co, ci=ci, ct=ct)})
        return out

    d = chrdre - ((aloci + aloco) / 2.) * (tanle - tante)
    c = (cfi + cfo) / 2.
    b = d - c
    tovca = st['tovca']
    if aloci < bo2 - bstro2:
        tovca = tovc
    if aloci > bo2 - bstro2:
        tovca = tovco
    gg = math.sin(0.5 * PI - sai000) * b / math.sin(0.5 * PI + sai000 -
                                                     lamhl)
    e = math.sin(0.5 * PI - sai100) * c / math.sin(0.5 * PI + sai100 -
                                                    lamhl)
    chrdpr = gg + e
    tovcp = tovca * d / chrdpr
    xhoc = gg / chrdpr
    inbord = not aloco / bo2 >= .99
    spr[7] = 2. * ((cfi + cfo) / 2.) * (aloco - aloci)
    if spr[35] == 1.:
        aforlf = ((2. * (aloco - aloci))**2) / spr[7] * beta
        taperd = False
    else:
        aforlf = spr[35]
        taperd = True
    dfl = calculate_dflcon(spr[33], spr[34], aforlf, beta, inbord, taperd)
    if dfl is not None:
        spr[12], spr[13] = dfl['cld'], dfl['clld']
        spr[17], spr[18] = dfl['cmd'], dfl['chd']
    tanle_w, tante_w = a[62], a[80]
    kaseno = int(st['kaseno'])
    pam = [0.0] * 8
    spamt = 0.0
    spamr = 0.0
    for i in (1, 2):
        y = float(i)
        if i == 1:
            v1, v2, v3, v4, v5, v6, v7 = (aloci, ci, cfi, aloco, co, cfo,
                                          sspne)
        else:
            v1, v2, v3, v4, v5, v6, v7 = (aloco, co, cfo, aloci, ci, cfi,
                                          0.)
        xyz = v7 - v1 if i == 1 else -(v7 - v1)
        zyx = v7 - v4 if i == 1 else -(v7 - v4)
        al1 = _alpha(xyz, v1 * tanle + v2 - v3 - v7 * tanle)
        al2 = _alpha(xyz, v1 * tanle + v2 - v7 * tanle)
        al3 = _alpha(zyx, v4 * tanle + v5 - v6 - v7 * tanle)
        al4 = _alpha(zyx, v4 * tanle + v5 - v7 * tanle)
        if mu > al1:
            kaseno = 1
        if mu <= al1 and mu > al2 and mu > al3:
            kaseno = 2
        if mu <= al2 and mu > al3:
            kaseno = 3
        if mu <= al3 and mu > al2:
            kaseno = 4
        if mu <= al2 and mu < al3 and mu > al4:
            kaseno = 5
        if mu <= al4:
            kaseno = 6
        skip = False
        if i == 2 and tanle == 0.0:
            spamt = 0.
            skip = True
        if not skip:
            if i == 2:
                x1, x4, x5, x6 = chrdre, xr, aloco, aloci
                x7, x8 = 1. - spr[34], 1. - spr[33]
                x9, x10 = -spr[34], -spr[33]
            else:
                x1, x4, x5, x6 = chrdtp, xt, sspne - aloci, sspne - aloco
                x7, x8 = 1. + spr[34], 1. + spr[33]
                x9, x10 = spr[34], spr[33]
            pam = [0.0] * 8
            spamt = 0.

            def ptcp(r, region):
                p = calculate_ptcp(r, region, y, tanle_w, beta, tante_w,
                                   spr[14])
                return p['pressure_ratio'], p['centre_of_pressure']

            def areas14():
                r1 = x1 / (beta * x5) - x7
                p1, t1 = ptcp(r1, 1)
                pam[0] = r1 * p1 * beta * x5**2 * (
                    2. * beta * x5 * (1. + x10 * t1) / t1 - 3. * x4)
                r4 = 1. - beta * x5 * x7 / x1
                p4, t4 = ptcp(r4, 4)
                pam[3] = r4 * p4 * x1**2 * (1. / (beta * x7)) * (
                    3. * x4 - 2. * x1 * (1. + x10 * t4) / (1. + x9 * t4))
                return pam[0] + pam[3]

            if kaseno != 6:
                r3 = 1. - (beta * x6 * x7) / x1
                p3, t3 = ptcp(r3, 3)
                pam[2] = r3 * p3 * x1**2 * (1. / (beta * x7)) * (
                    2. * x1 * (1. + x10 * t3) / (1. + x9 * t3) - 3. * x4)
                r7 = x1 / (beta * x6) - x7
                p7, t7 = ptcp(r7, 7)
                pam[6] = r7 * p7 * beta * x6**2 * (
                    3. * x4 - 2. * beta * x6 * (1. + x10 * t7) / t7)
                spamt = pam[2] + pam[6]
                if kaseno == 5:
                    pass
                elif kaseno >= 4:
                    spamt += areas14()
                else:
                    r5 = 1. - beta * x6 * x8 / x4
                    p5, t5 = ptcp(r5, 5)
                    pam[4] = r5 * p5 * x4**3 * (1. / (beta * x8))
                    r8 = x4 / (beta * x6) - x8
                    p8, t8 = ptcp(r8, 8)
                    pam[7] = r8 * p8 * beta * x6**2 * (
                        2. * beta * x6 * (1. + x10 * t8) / t8 - 3. * x4)
                    spamt += pam[4] + pam[7]
                    if kaseno < 3:
                        spamt += areas14()
                        if kaseno < 2:
                            r2 = x4 / (beta * x5) - x8
                            p2, t2 = ptcp(r2, 2)
                            pam[1] = r2 * p2 * beta * x5**2 * (
                                3. * x4 - 2. * beta * x5 * (1. + x10 * t2) /
                                t2)
                            r6 = 1. - beta * x5 * x8 / x4
                            p6, t6 = ptcp(r6, 6)
                            pam[5] = -r6 * p6 * x4**3 * (1. / (beta * x8))
                            spamt += pam[1] + pam[5]
        if i == 1:
            spr[37:45] = pam
        else:
            spamr = spamt
            spr[45:53] = pam
    spr[53] = ((-2. / (RAD * beta * math.sqrt(1. - spr[32]**2))) *
               (1. - (spamt + spamr) /
                ((aloco - aloci) * (cfi**2 + cfi * cfo + cfo**2))))
    spr[54] = spr[6] * spr[53]
    spr[55] = (cfi * (aloco - aloci) / 6.) * math.cos(lamhl) * (
        (cfo / cfi)**2 + cfo / cfi + 1.)
    c1, c2 = spr[2], spr[3]
    chabc = spr[53] * (1. - (.666667 * c2 * tovcp /
                             (c1 * (1 + k_) * math.cos(sai000 - lamhl))) *
                       (2. * (1. + 2. * xhoc) - k_ * (1. - xhoc)**2))
    spr[56] = spr[6] * spr[18]
    chdelr = (1. - 1.33333 * (c2 / c1) * tovcp * (1. + 2. * xhoc)) * spr[18]
    out.update({'spr': spr[1:], 'clalds': None, 'chabc': chabc,
                'chdelr': chdelr, 'kaseno': kaseno,
                'locals': dict(st, lamhl=lamhl, k=k_, mu=mu, xt=xt, xr=xr,
                               co=co, ci=ci, ct=ct, tovca=tovca,
                               kaseno=kaseno)})
    return out
