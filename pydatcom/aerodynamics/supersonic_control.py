"""
Supersonic trailing-edge control derivatives: DFLCON.

DFLCON evaluates the linear-theory closed forms of Section 6.1.6 for a
trailing-edge control on a supersonic wing: the lift, rolling-moment,
pitching-moment and hinge-moment effectiveness of a flap bounded by Mach
lines from its inboard and outboard edges.  The pressures ``P`` and their
centres ``TCP`` over each region of the control (the regions 1, 1A, 1B, 1C
about the inboard edge; 2 about the outboard edge for a control inboard of
the tip; 3 about the tip for a tip control) are combined with the region
areas and moments ``SL``, ``SLBX``, ``SLBY``.  The tapered forms take
``D``, the trailing-edge sweep parameter, and a second pass mirrors the
planform (``A``, ``D`` negated, ``FLD`` inverted) for the outboard edge.

Reference: datcom-legacy/datcom_2000/dflcon.f
"""

import math
from typing import Dict, Optional

from pydatcom.utils.constants import PI, RAD
from pydatcom.utils.legacy_numeric import arccos


def _sqra(x):
    return math.sqrt(abs(x))


def calculate_dflcon(aa: float, dd: float, aforlf: float, beta: float,
                     inbord: bool, taperd: bool
                     ) -> Optional[Dict[str, float]]:
    """Translate DFLCON: supersonic control effectiveness derivatives.

    Args:
        aa: ``AA``, the leading-edge (hinge-line) sweep parameter.
        dd: ``DD``, the trailing-edge sweep parameter (tapered only).
        aforlf: ``AFORLF``, the control's span-to-chord Mach parameter.
        beta: ``BETA``.
        inbord: The control lies inboard of the tip.
        taperd: The tapered forms.

    Returns:
        ``cld``, ``clld``, ``cmd``, ``chd`` (``SPR(12)``, ``(13)``,
        ``(17)``, ``(18)``), or ``None`` for a tapered control with
        ``|D| > 1``, where the routine returns without setting them.
    """
    p = [0.5, 0.5]
    p3 = 0.5
    a, d = float(aa), float(dd)
    two3rd = 2. / 3.
    a2 = a**2
    opa2 = 1. + a2
    oma2 = 1. - a2
    tcp0od = 4. / (RAD * beta * _sqra(oma2))
    a4 = a2**2
    sq3 = _sqra(oma2)
    oma4 = 1. - a4
    pa, pb, pc = [0.0, 0.0], [0.0, 0.0], [0.0, 0.0]
    tcp, tcpa, rcpb, tcpc = [0.0, 0.0], [0.0, 0.0], [0.0, 0.0], [0.0, 0.0]
    v: Dict[str, float] = {}
    if taperd:
        fld = aforlf
        if abs(d) > 1.0:
            return None
        if fld == 0.0:
            fld = 1.e-04
        d2 = d**2
        omd2 = 1. - d2
        sq1 = _sqra(oma2 * omd2)
        sq2 = _sqra(oma2 / omd2)
        opd2 = 1. + d2
        savfld = fld
        tip_done = False
        for k in range(2):
            apd, opa, opd = a + d, 1. + a, 1. + d
            amd, oma, omd = a - d, 1. - a, 1. - d
            fld2 = fld**2
            ad = a * d
            opad, omad = 1. + ad, 1. - ad
            omfld = 1. - fld
            amd2 = amd**2
            if amd < 0.0:
                sign, sign2 = 1., -1.
            else:
                sign, sign2 = -1., 1.
            g1 = 1. / (fld * omd - oma)
            g2 = fld * omad
            g3 = fld * amd
            g4 = fld / PI * arccos((oma2 - g2) / g3)
            g5 = sq2 / PI * arccos((omad - fld * omd2) / amd)
            g6 = 2. * fld * omad - fld2 * omd2 - oma2
            g7 = omfld * sq3 / PI * math.log(
                abs((a - fld * d + sign * _sqra(g6)) / omfld))
            g8 = arccos(d) / PI
            g9 = arccos(a) / PI
            g10 = 1. + 3. * (ad - d2) - ad * d2
            g11 = opd / omd * (oma2 * omd2 - 2. * d * oma**2)
            g115 = .5 / (fld * omd - oma)
            g12 = g115 / opd
            g15 = _sqra(oma2 * g6) / PI
            g16 = 1. / (opd * amd2)
            g17 = sq2 * (1. - g8)
            p[k] = .5 * (sq1 - oma * opd) / amd
            pa[k] = omd * g1 * (g4 - g5)
            pb[k] = (amd * g4 + g7) * g1
            pc[k] = omd / amd * (g17 + g9 - oma / omd)
            tcp[k] = .25 / (p[k] * amd2) * (sq2 * g10 - g11)
            tcpa[k] = g12 / (pa[k] * amd) * (
                omd2 * (2. * opad - fld * opd2) * g4 - g5 * g10 +
                sign * opd2 * g15)
            rcpb[k] = g115 / pb[k] * (amd * (2. * a - fld * apd) * g4 /
                                      omfld + sign * g15 + a * g7)
            tcpc[k] = .5 * g16 / pc[k] * (
                g10 * g17 + amd * opd2 * sq3 / PI - g11 +
                omd2 * (1. + 2. * ad - d2) * g9)
            if k == 1:
                a, d = -a, -d
                continue
            if not inbord:
                t2 = 2. + a + d
                sq4 = _sqra(opa * opd)
                t1 = 1.0 / (opd - fld * opa)
                t3 = arccos((t2 - 2. * fld * opa) / amd) / PI
                t4 = fld * arccos((2. * opd - fld * t2) / (fld * amd)) / PI
                t45 = omfld * (fld * opa - opd)
                sq5 = 2. * _sqra(opa * t45) / PI
                t5 = (2. - amd) * omd2 + 2. * (ad * opd + d * opa)
                p3 = (opa - sq4) / amd
                p3a = t1 * (opd * t3 - sq4 * t4)
                p3b = -t1 * (amd * t3 + sign * sq5)
                tcp3 = .25 * g16 / p3 * (sq4 * t5 - 2. * (
                    oma2 * omd2 + 2. * d * (opa**2)))
                tcp3a = -.25 * t1 / amd * (
                    2. * opd / fld * (2. * fld * (1. + ad) - opd2) * t3 +
                    sign2 * opd2 * sq5 - t4 * _sqra(opa / opd) * t5) / p3a
                rcp3b = 1. / (6. * p3b * t45) * (
                    3. * amd * (2. * a * fld - (a + d)) * t3 +
                    sign2 * (2. * fld * (1. - 2. * a) - (2. - 3. * a - d)) *
                    sq5)
            opfld = 1. + fld
            omfld2 = 1. - fld2
            sl1 = 2. * amd / (omfld2 * omd2)
            sl2 = fld2 * sl1
            sl1a = fld * omd - oma
            sl1b = sl1a / (opfld * amd)
            sl1a = sl1a / (omfld2 * omd)
            sl3 = .5 * omd2 * sl2 / opd
            sl0 = (oma / omd - fld2 * opa / opd) / omfld2
            h1 = two3rd * amd / (omfld * opd2)
            slby1 = sl1 * h1 * (tcp[0] - d)
            slby1a = sl1a * h1 * (tcpa[0] - d)
            slby1b = two3rd * sl1b
            if not inbord:
                slby3 = sl3 * (1. - h1 * fld * (tcp3 + d))
            sly0 = (opa * omfld - amd2 * amd / ((omfld * omd)**2) -
                    (opd - fld * opa)**3 / ((omfld * opd)**2)) / \
                (3. * amd * opfld)
            fld3 = fld2 * fld
            omfld3 = 1. - fld3
            rat = omfld2 / omfld3
            rat1 = rat / opd2
            rat2 = amd / (omfld3 * opd2)
            rat3 = 1. / (omfld3 * opd * opd2)
            slx0 = .5 * ((oma / omd)**2 - fld3 * (opa / opd)**2) / omfld3
            slbx1 = sl1 * rat1 * (opad - amd * tcp[0])
            slbx1a = sl1a * rat1 * (opad - amd * tcpa[0])
            slbx1b = sl1b * rat * (omfld * (rcpb[0] - a)) / amd
            slbx1c = rat2 * (opad - amd * tcpc[0]) / omd
            if not inbord:
                slbx3 = fld3 * rat2 * (opad + amd * tcp3) / opd
                slbx3a = rat3 * fld2 * (fld * opa - opd) * (opad +
                                                             amd * tcp3a)
                slbx3b = omfld**2 * (fld * opa - opd) * (rcp3b + a) / \
                    (omfld3 * (amd**2))
                tip_done = True
                break
            fld = 1. / fld
            a, d = -a, -d
        if not tip_done:
            fld = savfld
            fld2 = fld**2
            omfld = 1. - fld
            opad = 1. + a * d
            amd = a - d
            opd2 = 1. + d**2
            opa = 1. + a
            opd = 1. + d
            slby2 = sl2 * (1. - h1 * fld * (tcp[1] + d))
            slbx2 = sl2 * fld * rat * (opad + amd * tcp[1]) / opd2
            tmp = fld * opa - opd
            slbx2a = fld2 * tmp * (opad + amd * tcpa[1]) * rat3
            slbx2b = omfld**2 * tmp * (rcpb[1] + a) / (omfld3 * (amd**2))
            slbx2c = fld3 * amd * (opad + amd * tcpc[1]) * rat3
    else:
        fld = aforlf
        fld2 = fld**2
        for k in range(2):
            opa, oma = 1. + a, 1. - a
            g1 = 1. + 2. * a * fld - oma2 * fld2
            g2 = oma2 * fld - a
            g3 = arccos(g2)
            g4 = 1. / ((1. - oma * fld) * PI)
            g5 = fld * sq3
            g7 = _sqra(oma2 * g1)
            g8 = g4 / opa
            g9 = math.log(abs(1. + a * fld - _sqra(g1)) / fld)
            g10 = arccos(a)
            g11 = 1. + a / PI * g10 - sq3 / PI
            g12 = g8 / (4. * oma2)
            g13 = 7. * a2 - 2. * a4 + 1.
            g14 = 2. * a + opa2 * fld
            g15 = 2. * fld * oma2**2
            g16 = a * (7. - a2)
            g17 = 1. / (opa * oma2)
            pa[k] = g8 * (g7 - g2 * g3)
            pb[k] = g4 * (g3 + g5 * g9)
            pc[k] = g11 / opa
            tcp[k] = (8. * a + opa2) / (4. * oma2)
            tcpa[k] = g12 / pa[k] * ((g13 - g14 * g15) * g3 +
                                     (g16 + oma4 * fld) * g7)
            rcpb[k] = g4 / (2. * pb[k]) * ((1. + 2. * a * fld) * g3 / fld -
                                           g7 + a * g5 * g9)
            tcpc[k] = g17 / (4. * pc[k]) * (8. * a + opa2 + g13 / PI * g10 -
                                            g16 / PI * sq3)
            if k == 1:
                a = -a
                oma, opa = 1. - a, 1. + a
                sl2 = sl1
                slby2 = sl2 * (1. - (1. - 4. * a) / (6. * fld * oma2))
                slbx2 = two3rd * sl2
                slbx2a = (1. - fld * opa) / (3. * fld * opa)
                slbx2b = fld * (rcpb[1] + a) * (1. - fld * opa) / 3.
                slbx2c = 1. / (3. * fld * opa)
                break
            if not inbord:
                t1 = 1. - opa * fld
                t2 = 2. * fld * opa - 1.
                t3 = 5. * a2 + 8. * a - 3.
                t4 = opa * t1
                t45 = _sqra(t4 * fld)
                t5 = arccos(t2)
                p3a = 0.5 / (t1 * PI) * (2. * t45 - t2 * t5)
                p3b = 1. / (PI * t1) * (t5 - 2. * t45)
                tcp3a = -1. / (t4 * 16. * p3a * PI) * (
                    (t3 + 8. * fld * opa**2 * (fld * opa2 - 2. * a)) * t5 +
                    2. * (t3 - 2. * fld * opa * opa2) * t45)
                rcp3b = 1. / (6. * p3b * fld * t1 * PI) * (
                    3. * (1. - 2. * a * fld) * t5 -
                    2. * (1. + 2. * fld * (1. - 2. * a)) * t45)
                sl3 = 1. / (2. * fld * opa)
                slby3 = sl3 * (1. - (5. / (12. * fld * opa)))
                slbx3 = two3rd * sl3
                slbx3a = (1. - fld * opa) / (3. * fld * opa)
                slbx3b = fld * (rcp3b + a) * (1. - fld * opa) / 3.
            sl1 = 1. / (fld * oma2)
            sl1a = (1. - fld * oma) / (2. * fld * oma)
            sl1b = (1. - fld * oma) / 2.
            slby1 = sl1 * (1. + 4. * a) / (6. * fld * oma2)
            slby1a = sl1a * two3rd * (tcpa[0] - a) / (fld * opa2)
            slby1b = two3rd * sl1b
            slbx1 = two3rd * sl1
            slbx1b = two3rd * sl1b * fld * (rcpb[0] - a)
            slbx1a = two3rd * sl1a
            slbx1c = 1. / (3. * fld * oma)
            if not inbord:
                break
            a = -a
        sl0 = (fld * oma2 - 1.) / (fld * oma2)
        sly0 = (3. * fld * oma * oma2 * (fld * opa - 1.) - 4. * a) / \
            (6. * fld2 * oma2**2)
        slx0 = (3. * fld * oma2 - 4.) / (6. * fld * oma2)

    psl1 = p[0] * sl1
    pslx1 = p[0] * slbx1
    psly1 = p[0] * slby1
    pslx1c = pc[0] * slbx1c
    pslx1a = pa[0] * slbx1a
    pslx1b = pb[0] * slbx1b
    if not taperd:
        fld = 1. / fld
    testr1, testr2 = 1. - a, 1. + a
    if taperd:
        testr1 = (1. - a) / (1. - d)
        testr2 = (1. + d) / (1. + a)

    def first_test():
        if taperd and not fld < 1.:
            return not fld >= testr1
        return not fld <= testr1

    def second_test():
        if taperd and not fld < 1.:
            return not fld >= testr2
        return not testr2 >= fld

    if inbord:
        psl2 = p[1] * sl2
        psly2 = p[1] * slby2
        pslx2 = p[1] * slbx2
        pslx2c = pc[1] * slbx2c
        pslx2a = pa[1] * slbx2a
        pslx2b = pb[1] * slbx2b
        cld = sl0 + psl1 + psl2
        cmd = -(slx0 + pslx1 + pslx2)
        clld = sly0 + psly1 + psly2
        chd = -(slx0 + pslx1c + pslx2c)
        if first_test():
            chd = chd - (pslx1a - pslx1b)
        if second_test():
            chd = chd - (pslx2a - pslx2b)
    else:
        psl1a = pa[0] * sl1a
        psl1b = pb[0] * sl1b
        psl3 = p3 * sl3
        pslx3 = p3 * slbx3
        psly1a = pa[0] * slby1a
        psly1b = pb[0] * slby1b
        psly3 = p3 * slby3
        pslx3a = p3a * slbx3a
        pslx3b = p3b * slbx3b
        cld = sl0 + psl1 + psl3
        cmd = -(slx0 + pslx1 + pslx3)
        clld = sly0 + psly1 + psly3
        chd = -(slx0 + pslx1c + pslx3)
        if first_test():
            cld = cld + psl1a - psl1b
            cmd = cmd - (pslx1a - pslx1b)
            clld = clld + psly1a - psly1b
            chd = chd - (pslx1a - pslx1b)
        if second_test():
            chd = chd - (pslx3a - pslx3b)
    return {'cld': tcp0od * cld, 'clld': tcp0od * clld,
            'cmd': tcp0od * cmd, 'chd': tcp0od * chd,
            'method': 'legacy_dflcon'}
