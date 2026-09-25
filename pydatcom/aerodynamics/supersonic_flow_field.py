"""
Supersonic wing flow field at the tail: SDWASH and its overlay M21O25.

SDWASH forms the downwash gradient and angle at the horizontal tail and
the dynamic-pressure ratio there.  For a wing with subsonic leading and
trailing edges it builds the downwash from the vortex-sheet charts of
Figure 4.4.1-76 (SDWA-SDWE), with the vortex core dropped by Figure
4.4.1-74 (SDDVC), averaged over the tail's root and tip and integrated in
angle of attack; a supersonic edge instead takes the simple
``1.62 CL / (pi A)`` estimate.  The pressure ratio comes from the wake
(Section 4.4.1's viscous wake) when the tail is inside it and from
DPRESR's wave calculation when it is not; a tail Mach number at or below 1
marks the angle where the method stops (``JDETCH``).

Reference: datcom-legacy/datcom_2000/sdwash.f, m21o25.f
"""

import math
from typing import Dict, Mapping

import numpy as np

from pydatcom.aerodynamics.supersonic_downwash import (calculate_sdwa,
                                                       calculate_sdwb,
                                                       calculate_sdwc,
                                                       calculate_sdwd,
                                                       calculate_sdwe)
from pydatcom.aerodynamics.tail_pressure import calculate_dpresr
from pydatcom.aerodynamics.vortex_core import calculate_sddvc
from pydatcom.utils.constants import PI, RAD
from pydatcom.utils.legacy_numeric import tbfunx, trapz

_TAPER = np.array([0.00, 0.25, 0.50, 1.00])


def _icase(swepte: float, swepc2: float) -> int:
    if swepte == 0.0:
        return 1
    if swepte > 0.0:
        return 4
    return 2 if swepc2 <= 0.0 else 3


def calculate_sdwash(data: Mapping[str, object],
                     user_epsilon: bool = False) -> Dict[str, object]:
    """Translate SDWASH: downwash and dynamic pressure at the tail.

    Args:
        data: The COMMON words read, by name:

            - ``alpha`` (``FLC(23..)``), ``mach`` (``FLC(NZ+2)``),
              ``sref``, ``nf`` (``/OVERLY/``'s message flag).
            - ``position``: ``xw``, ``zw``, ``aliw``, ``xh``, ``zh``.
            - ``a``: the wing's ``A`` words 3, 11, 12, 16, 24, 58, 70, 76,
              118, 120 (INFTGM's 11, 12, 16 and 24 as it leaves them).
            - ``win``: ``WINGIN`` 1, 4, 6 and the slopes 95-100.
            - ``tail``: ``tanle`` (``AHT(62)``), ``span``, ``spandi``,
              ``dihei``, ``diheo`` (``HTIN(4)``, ``(12)``, ``(13)``,
              ``(14)``).
            - ``wing``: ``cl`` (``WING(21..)``), ``cla`` (``WING(101..)``),
              ``cd0`` (``SLG(80)``).
            - ``dwangl``: ``DWANGL`` as supplied (read only with
              ``user_epsilon``), ``qqinfy``: ``QQINFY`` before the call
              (kept where DPRESR leaves it unset).

        user_epsilon: ``KEPSLN``, the downwash angles supplied as data.

    Returns:
        ``nalpha`` (the angle count as the routine leaves ``/OVERLY/``),
        ``nf``, ``jdetch``, ``beta``, and when the vortex-sheet method runs
        ``x``, ``y``, ``z``, ``dhb``, ``zeff``, ``depx``, ``depavg``; always
        ``alpha`` (with incidence), ``clanl``, the ``/IDWASH/`` curves
        ``qqinfy``, ``dwangl``, ``depda``, and the ``/SUPDW/`` words ``m``,
        ``zc``, ``zwakec``, ``delqo`` and ``dpresr`` (the last DPRESR's
        ``DWA(231..236)``, by name).
    """
    a = {int(k): float(v) for k, v in data['a'].items()}
    w = {int(k): float(v) for k, v in data['win'].items()}
    pos = {k: float(v) for k, v in data['position'].items()}
    tl = {k: float(v) for k, v in data['tail'].items()}
    mach, sref = float(data['mach']), float(data['sref'])
    nalpha = len(data['alpha'])
    nf = int(data['nf'])
    alpha = [float(v) + pos['aliw'] for v in data['alpha']]
    clw = [float(v) for v in data['wing']['cl']]
    cla = [float(v) for v in data['wing']['cla']]
    claw = cla[0]
    cdow = float(data['wing']['cd0'])
    jdetch = -1
    beta = math.sqrt(mach**2 - 1.0)
    abeta = a[120] * beta
    visdw = (math.tan(a[58] / RAD) / beta > 1.0 or
             math.tan(a[76] / RAD) / beta > 1.0)
    r: Dict[str, object] = {'beta': beta, 'alpha': alpha, 'clanl': cla}
    qqinfy = [float(v) for v in data['qqinfy']]
    dwangl = [float(v) for v in data['dwangl']]
    depda = [0.0] * nalpha
    alpha_arr = np.asarray(alpha)

    if not (visdw or user_epsilon):
        span, spanh = w[4], tl['span']
        x = [(pos['xh'] - pos['xw']) / (span * beta),
             (pos['xh'] - pos['xw'] + tl['tanle'] * spanh) / (span * beta)]
        y = [0.0, spanh / span]
        z1 = pos['zh'] - pos['zw']
        z2 = z1 + math.tan(tl['dihei'] / RAD) * (spanh - tl['spandi']) + \
            math.tan(tl['diheo'] / RAD) * tl['spandi']
        z = [z1 / span, z2 / span]
        swepte, swepc2, tapr = a[76], a[70], a[118]
        icase = _icase(swepte, swepc2)
        swepr = math.atan((w[6] - w[1]) / (2.0 * span)) * RAD
        dhb = list(calculate_sddvc(x, abeta, tapr, icase, swepte,
                                   swepr)['dhb'])
        zeff, depx, depavg = [], [], []
        for angle_slot in range(nalpha):
            na = min(max(int(abs(2.0 * alpha[angle_slot]) + 1.5), 2), 21)
            xna = float(na - 1)
            alp, depa = [], []
            for substep in range(na):
                if nf >= 0:
                    nf = -1
                if substep == na - 1 and nf == -1:
                    nf = 0
                alp.append(alpha[angle_slot] * substep / xna)
                sdw = [0.0, 0.0]
                for station in range(2):
                    ze = (z[station] + dhb[station] * beta *
                          alp[substep] / RAD)
                    if substep == na - 1:
                        zeff.append(ze)
                    ze = abs(ze)
                    if icase in (1, 3):
                        dep = [calculate_sdwa(x[station], y[station], ze,
                                              abeta)['sdw'],
                               *calculate_sdwc(x[station], y[station], ze,
                                               abeta)['sdw'],
                               calculate_sdwb(x[station], y[station], ze,
                                              abeta)['sdw']]
                        sdw[station] = tbfunx(_TAPER, dep, tapr, 0, 0)[0]
                        sdw3 = sdw[station]
                    if icase in (2, 3):
                        dep = [*calculate_sdwd(x[station], y[station], ze,
                                               abeta)['sdw'],
                               calculate_sdwb(x[station], y[station], ze,
                                              abeta)['sdw']]
                        sdw[station] = tbfunx(_TAPER, dep, tapr, 0, 0)[0]
                        if icase == 3:
                            sdw[station] = (sdw[station] +
                                            swepc2 * (sdw3 - sdw[station]) /
                                            swepr)
                    if icase == 4:
                        sdw[station] = calculate_sdwe(
                            x[station], y[station], ze, abeta, tapr)['sdw']
                    if substep == na - 1:
                        depx.append(sdw[station])
                clanl_at_substep = tbfunx(alpha_arr, cla, alp[substep], 0, 0)[0]
                depa.append((sdw[0] + sdw[1]) / 2.0 *
                            clanl_at_substep / claw)
            depavg.append(depa[-1] * claw / clanl_at_substep)
            depda[angle_slot] = depa[-1]
            dwangl[angle_slot] = float(trapz(depa, alp)[0])
        r.update({'x': x, 'y': y, 'z': z, 'dhb': dhb, 'zeff': zeff,
                  'depx': depx, 'depavg': depavg, 'icase': icase})

    m = [0.0] * nalpha
    zc = [0.0] * nalpha
    gamma = a[11]
    zwakec = delqo = None
    dpresr = None
    stopped = False
    for angle_slot in range(nalpha):
        if not user_epsilon and visdw:
            dwangl[angle_slot] = (1.62 * clw[angle_slot] / (PI * a[120]) *
                                  RAD * sref / a[3])
        arg = (dwangl[angle_slot] - alpha[angle_slot]) / RAD + gamma
        arg4 = a[24] * math.cos(arg) / (math.cos(gamma) * a[16])
        arg1 = cdow * (arg4 + 0.15) * sref / a[3]
        zwakec = 0.68 * math.sqrt(arg1)
        zwaket = zwakec * a[16]
        delqo = 2.42 * math.sqrt(cdow * sref / a[3]) / (arg4 + 0.3)
        zc[angle_slot] = arg4 * math.tan(arg)
        test = abs(zc[angle_slot] / zwakec)
        if test <= 1.:
            qqinfy[angle_slot] = 1. - delqo * math.cos(PI * test / 2.0)**2
            m[angle_slot] = mach
            continue
        dpresr = calculate_dpresr(
            zc[angle_slot] * a[16], zwaket, alpha[angle_slot],
            dwangl[angle_slot] / RAD, mach, w[6], a[12], a[24],
            [w[95 + word] for word in range(6)])
        if dpresr['qqinfy'] is not None:
            qqinfy[angle_slot] = dpresr['qqinfy']
        m[angle_slot] = dpresr['mj']
        if m[angle_slot] <= 1.0:
            jdetch = angle_slot
            stopped = True
            break
    r.update({'qqinfy': qqinfy, 'm': m, 'zc': zc, 'zwakec': zwakec,
              'delqo': delqo, 'dpresr': dpresr, 'jdetch': jdetch, 'nf': nf})
    if stopped and not visdw:
        r.update({'nalpha': nalpha, 'dwangl': dwangl, 'depda': depda})
        return r
    if visdw and not user_epsilon and jdetch > 0:
        nalpha = jdetch
    if jdetch == 0:
        nalpha = 0
    if visdw or user_epsilon:
        for angle_slot in range(nalpha):
            depda[angle_slot] = tbfunx(
                alpha_arr[:nalpha], dwangl[:nalpha], alpha[angle_slot], 0, 0)[1]
    r.update({'nalpha': nalpha, 'dwangl': dwangl, 'depda': depda})
    return r
