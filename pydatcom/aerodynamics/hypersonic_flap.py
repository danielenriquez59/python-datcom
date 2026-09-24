"""
Hypersonic flap increments: HYPFLP and its overlay M42O52.

HYPFLP forms the normal-force, axial-force and hinge-moment increments of
a flat-plate trailing-edge flap at hypersonic speed (Section 6.3.1): the
local flow ahead of the hinge from HYPROP, the attached-flow pressure on
the flap (FIG68's oblique shock, averaged with the isentropic compression
below Mach 6 for a laminar layer), the incipient-separation pressure, and,
where the flap separates the boundary layer, the separation point found by
SIMUL2, the plateau pressure and the interaction lengths of Figures
6.3.1-59 to -70.

Reference: datcom-legacy/datcom_2000/hypflp.f, m42o52.f
"""

import math
from typing import Dict, Mapping, Optional

from pydatcom.aerodynamics.hyprop import calculate_hyprop
from pydatcom.utils.constants import RAD
from pydatcom.utils.legacy_numeric import simul2
from pydatcom.utils.table_lookup import fig68


def calculate_hypflp(data: Mapping[str, object],
                     hyprop=calculate_hyprop) -> Dict[str, object]:
    """Translate HYPFLP: hypersonic flap force and moment increments.

    Args:
        data: By name: ``alpha`` (``FLC(23..)``), ``mach``, ``rl``
            (``FLC(I+42)``, the unit Reynolds number), ``f`` (``/FLAPIN/``:
            1 altitude, 2 hinge-line station ``XHL``, 3 wall temperature
            ``TWOTI``, 4 flap chord ``CF``, 5.. deflections, 16 count),
            ``laminar`` (``F(15)``), ``sref``, ``cbar``, ``stale``: the
            saved local ``phe`` and the previous ``cpi2``.
        hyprop: The local-flow routine, HYPROP by default.

    Returns:
        ``dcn`` (``BODY``), ``dca`` (``WING``), ``dcmcn`` (``HT``),
        ``dcmca`` (``VT``), keyed ``10(J-1)+N`` as the source, the ``HYP``
        curves ``paopi``, ``taoti``, ``malp``, ``raori``, and ``stale``.

    Notes:
        An attached-flow deflection reads ``CPIP``, the plateau pressure
        left by the last separated one, but with ``D1 = D2 = D3 = LFI = 0``
        it cancels from all four increments.  Kept as executed: the
        turbulent branch forms the single-shock pressure before
        setting a detached shock's angle to 90 degrees, where the laminar
        branch sets it first; separation at the hinge line's leading edge
        skips ``CPI2``, so the previous deflection's value is used; and
        the dead-air angle there is ``ATAN(D2*SINDF)/(D1+D2*COSDF)``, the
        arctangent of the numerator alone.
    """
    f = {int(k): float(v) for k, v in data['f'].items()}
    lamnr = bool(data['laminar'])
    alt, xhl, twoti, cf = f[1], f[2], f[3], f[4]
    mach, rl = float(data['mach']), float(data['rl'])
    sr, cbar = float(data['sref']), float(data['cbar'])
    ndelta = int(f[16] + .5)
    st = dict(data.get('stale', {}))
    phe = float(st.get('phe', 0.0))
    cpi2 = float(st.get('cpi2', 0.0))
    cpip = 0.0
    out = {k: {} for k in ('dcn', 'dca', 'dcmcn', 'dcmca')}
    hyp = {k: {} for k in ('paopi', 'taoti', 'malp', 'raori')}
    d1 = d2 = d3 = lfi = 0.0
    rlstrl = 0.0
    for angle_index, alpha in enumerate(data['alpha']):
        alpha = float(alpha)
        if alpha < 0 or alpha > 20.0:
            continue
        h = hyprop(alt, mach, alpha)
        paopi, taoti = float(h['pressure']), float(h['temperature'])
        malp, raori = float(h['mach']), float(h['density'])
        hyp['paopi'][angle_index], hyp['taoti'][angle_index] = paopi, taoti
        hyp['malp'][angle_index], hyp['raori'][angle_index] = malp, raori
        cpia = (paopi - 1.) / (.7 * mach**2)
        rlhl = xhl * raori * rl
        m2 = malp**2
        for deflection_index in range(1, ndelta + 1):
            delta = f[deflection_index + 4]
            sindf, cosdf = math.sin(delta / RAD), math.cos(delta / RAD)
            if lamnr:
                cpinc = 2.03 * (malp**2 - 1.)**(-0.306) / rlhl**0.25
                theta, ier = fig68(malp, delta)
                if ier == 2:
                    theta = 90.0
                p2pass = (7. * m2 * math.sin(theta / RAD)**2 - 1.) / 6.
                if delta == 0.0:
                    p2pass = 1.0
                if malp > 6.0:
                    p2pa = p2pass
                else:
                    arg = 0.2 * m2 * sindf / math.sqrt(m2 - 1.)
                    p2pa = 0.5 * (p2pass + (1. + arg)**7)
            else:
                cpinc = 2.2 / rlhl**0.1
                theta, ier = fig68(malp, delta)
                p2pass = (7. * m2 * math.sin(theta / RAD)**2 - 1.) / 6.
                if ier == 2:
                    theta = 90.0
                if delta == 0.0:
                    p2pass = 1.0
                p2pa = p2pass
            cpa2 = (p2pa - 1.) / (.7 * malp**2)
            separated = cpa2 >= cpinc
            skip_cpi2 = False
            x0 = 0.0
            if not separated:
                lfi = d1 = d2 = d3 = 0.0
            else:
                x8, x0 = xhl / 8.0, xhl / 20.
                arg = raori * rl
                a1, a2, a3 = [], [], []
                for probe in range(1, 9):
                    x0 = x0 + x8
                    rlx0 = arg * x0
                    twota = twoti / taoti
                    if not lamnr:
                        cpap = 1.91 * (malp**2 - 1.)**(-0.309) / rlx0**0.1
                        arg4 = 0.7 * cpap * malp**2 + 1.
                        d1od0 = 1.1e6 * (malp**(-1.67) * (arg4 - 1.))**8.55
                        if probe == 1:
                            rlstrl = (0.28 + 0.5 * twota + 0.22 *
                                      (1. + 0.2 * malp**2 * 0.72**0.3333))
                        rlstr = rlstrl * arg
                        d0 = 0.154 / rlstr**0.142857 * x0**0.85714
                    else:
                        cpap = 1.56 * (malp**2 - 1.)**(-0.262) / rlx0**0.25
                        arg4 = 0.7 * cpap * malp**2 + 1.
                        d1od0 = 5.69e5 * malp**(-4.1) * (arg4 - 1.)**3.5
                        if probe == 1:
                            rlstrl = (0.28 + 0.5 * twota + 0.22 *
                                      (1. + 0.2 * malp**2 *
                                       0.72**0.3333))**(-1.67)
                        rlstr = rlstrl * arg
                        d0 = 5.2 * math.sqrt(x0) / math.sqrt(rlstr)
                    d1 = d1od0 * d0
                    a1.append(x0 + d1)
                    a2.append(xhl)
                    a3.append(x0)
                x0, x0pd1 = simul2(a3, a2, a1)
                if x0 == -1000.:
                    x0 = 0.0
                if x0 == 0.0:
                    skip_cpi2 = True
                else:
                    d1 = x0pd1 - x0
                    rlx0 = arg * x0
                    arg1 = (malp / mach)**2 * paopi
                    rlstr = rlstrl * arg
                    if not lamnr:
                        d0 = 0.154 / rlstr**0.142857 * x0**0.85714
                        cpap = 1.91 * (malp**2 - 1.)**(-0.309) / rlx0**0.1
                        cpip = cpap * arg1 + cpia
                        arg4 = 0.7 * cpap * malp**2 + 1.
                        lfiod0 = 1.84e4 * ((arg4 - 1.) / malp**1.325)**8.4
                        if malp >= 3.9:
                            lfiod0 = d1 / d0
                        d1od0 = 1.1e6 * (malp**(-1.67) * (arg4 - 1.))**8.55
                    else:
                        d0 = 5.2 * math.sqrt(x0) / math.sqrt(rlstr)
                        cpap = 1.56 * (m2 - 1.)**(-0.262) / rlx0**0.25
                        cpip = cpap * arg1 + cpia
                        pppa = 0.7 * cpap * m2 + 1.
                        lfiod0 = 2.47e5 * malp**(-4.2) * (pppa - 1.)**3.45
                        d1od0 = 5.69e5 * malp**(-4.1) * (pppa - 1.)**3.5
                    lfi = lfiod0 * d0
                    d1 = d1od0 * d0
                    if not lamnr:
                        d3 = 0.0
                    else:
                        pppa = 0.7 * cpap * malp**2 + 1.
                        a = math.sqrt((6.0 * pppa + 1.0) / (7. * malp**2))
                        th = math.atan(a / math.sqrt(1. - a**2))
                        a = m2 * math.sin(2. * th) - 2. / math.tan(th)
                        a1_ = 10. + m2 * (7. + 5. * math.cos(2. * th))
                        phe = math.atan(5. * a / a1_) * RAD
            if not skip_cpi2:
                cpi2 = cpa2 * (malp / mach)**2 * paopi + cpia
            d3 = 0.0
            if cpa2 >= cpinc:
                cfod1 = cf / d1
                madf = malp * delta / RAD
                if lamnr:
                    if madf >= 5.:
                        d2 = .344 * d1 * math.sqrt(cfod1)
                        if cfod1 >= 1.:
                            d2 = d1 * 0.344
                        if cfod1 <= 0.25:
                            d2 = d1 * 0.172
                    else:
                        d2 = d1 * (0.545 - 0.0403 * madf) * math.sqrt(cfod1)
                        if cfod1 >= 1.:
                            d2 = d1 * (0.545 - 0.04 * madf)
                        if cfod1 <= 0.25:
                            d2 = d1 * (.273 - 0.02 * madf)
                else:
                    if madf >= 2.4:
                        d2 = d1 * 0.37 * math.sqrt(cfod1)
                        if cfod1 >= 1.0:
                            d2 = d1 * 0.37
                        if cfod1 <= 0.25:
                            d2 = d1 * 0.37 * 0.5
                    else:
                        d2 = (1.16 - 0.33 * madf) * math.sqrt(cfod1) * d1
                        if cfod1 >= 1.:
                            d2 = d1 * (1.16 - 0.33 * madf)
                        if cfod1 <= 0.25:
                            d2 = d1 * (0.58 - 0.165 * madf)
                if x0 == 0.0:
                    d1 = xhl
                    lfi = 0.0
                    phe = math.atan(d2 * sindf) / (d1 + d2 * cosdf) * RAD
                    theta, ier = fig68(malp, phe)
                    if ier == 2:
                        theta = 90.0
                    pppa = (7. * m2 * math.sin(theta / RAD)**2 - 1.) / 6.
                    if phe == 0.0:
                        pppa = 1.0
                    cpap = (pppa - 1.) / (0.7 * m2)
                    cpip = cpap * (malp / mach)**2 * paopi + cpia
                if lamnr:
                    arg1 = math.tan(delta / RAD) - math.tan(phe / RAD)
                    d3 = d1 * ((math.tan(phe / RAD) / arg1) /
                               math.cos(delta / RAD))
            key = angle_index * 10 + (deflection_index - 1)
            dp, dq = cpip - cpia, cpi2 - cpip
            out['dcn'][key] = (dp * (d1 - lfi / 2. + cf * cosdf) +
                               dq * cosdf * (cf - 0.5 * (d2 + d3))) / sr
            out['dca'][key] = sindf / sr * (cf * dp +
                                            dq * (cf - 0.5 * (d2 + d3)))
            out['dcmcn'][key] = (
                dp * (lfi**2 / 6. + .5 * (-lfi * d1 + d1**2 -
                                          (cf * cosdf)**2)) -
                dq * cosdf**2 * (cf**2 / 2. - (d2**2 + d2 * d3 + d3**2) /
                                 6.)) / (sr * cbar)
            out['dcmca'][key] = (-math.sin(delta / RAD)**2 / (sr * cbar)) * (
                dp * cf**2 / 2.0 +
                dq * (cf**2 / 2.0 - (d2**2 + d2 * d3 + d3**2) / 6.))
    out.update(hyp)
    out['stale'] = {'phe': phe, 'cpi2': cpi2}
    out['method'] = 'legacy_hypflp'
    return out
