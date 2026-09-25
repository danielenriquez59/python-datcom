"""
Wing pitching moment due to flaps: FLAPCM and its overlay M37O45.

FLAPCM spreads the flap's lift over the span with GDELTA's loading (the
difference between the outboard- and inboard-edge flap spans), places
each strip's centre of pressure (Figure 6.1.5.1-67A/B for the carry-over
outside the flap, Figure 6.1.2.1-35B for a plain flap, the section moment
of slotted and Fowler flaps, or Figure 6.1.2.1-36 for leading-edge
devices), and integrates the moment increment across the span for each
deflection.

The arrays mirror the source's 1-based COMMON layout, and the routine's
saved state (the span stations ``ET``, which it edits, and ``KOUNT``,
``INDEX``, ``ARG`` and the ``/SUPWH/`` words it reuses) is explicit.

Reference: datcom-legacy/datcom_2000/flapcm.f, m37o45.f
"""

import math
from typing import Dict, List, Mapping

import numpy as np

from pydatcom.aerodynamics.flap_loading import calculate_gdelta
from pydatcom.utils.constants import RAD, UNUSED
from pydatcom.utils.legacy_interp import interx
from pydatcom.utils.legacy_numeric import trapz
from pydatcom.utils.legacy_tables import tlinex

_ETAG = [0.924, .707, .383, 0.0]
ET_DATA = [0.0, .1423, .2817, .4153, .5407, .6549, .7557, .8412, .9097,
           .9595, .9898, 1.0, .5, .75]
_X5126A = [-.02, 0., .02, .04, .06, .08, .10, .12, .14, .16, .18, .20, .22]
_Y5126A = [1., 1., .955, .845, .665, .485, .33, .21, .12, .06, .018, 0., 0.]
_X2126B = [0., .1, .2, .3, .4, .5, .6, .7, .8, .9, 1.]
_X1126B = [0., 20., 20.01, 60.]
_Y5126B = [.745, .698, .65, .6, .55, .5, .45, .4, .35, .3, .25,
           .745, .698, .65, .6, .55, .5, .45, .4, .35, .3, .25,
           .75, .704, .675, .63, .6, .593, .59, .592, .598, .6, .6,
           .75, .704, .675, .63, .6, .593, .59, .592, .598, .6, .6]
_X1215B = [.1, .2, .25, .3, .5]
_X2215B = [0., 10., 20., 30., 40., 50., 70.]
_Y1215B = [0., -.05, -.105, -.14, -.163, -.175, -.180,
           0., -.086, -.16, -.20, -.219, -.23, -.24,
           0., -.11, -.195, -.245, -.27, -.28, -.29,
           0., -.09, -.165, -.22, -.24, -.26, -.26,
           0., -.078, -.145, -.2, -.235, -.26, -.28]
_X12136 = [0., .05, .1, .2, .3, .35, .4, .45, .5]
_Y12136 = [0., -.00025, -.00065, -.00165, -.0026, -.0031, -.00345,
           -.0036, -.00375]


def _tl(x1, x2, flat, q1, q2):
    grid = np.asarray(flat, dtype=float).reshape(len(x1), len(x2)).T
    return float(tlinex(x1, x2, grid, q1, q2, 0, 0, 0, 0))


def calculate_flapcm(data: Mapping[str, object]) -> Dict[str, object]:
    """Translate FLAPCM: the pitching-moment increment of a flap.

    Args:
        data: By name: ``htpl``, ``asyfp``, ``mach`` (``FLC(II+2)``),
            ``ii`` (the Mach index, for the section lift slope),
            ``surface`` (the wing's or tail's words: ``clasec``
            (``WINGIN(II+20)``), ``tanc4`` (``A(68)``), ``swstr``
            (``A(3)``), ``bsto2``, ``bo2`` (``WINGIN(3)``, ``(4)``),
            ``tante``, ``tanle`` (``A(80)``, ``(62)``), ``x`` (``XW`` or
            ``XH``), ``cr`` (``A(10)``), ``tapexp`` (``A(27)``),
            ``arstar`` (``A(7)``), ``alpo`` (``A(134)``), ``cmo``
            (``WINGIN(61)``)), ``tail`` (for ``asyfp``: GDELTA's tail
            words), ``f`` (``/FLAPIN/`` 1-17, 19.., 29.., 39.., 49..),
            ``flp`` (``/POWR/``'s ``FLP`` words 1-5, 28-32, 60, 110..149,
            150..189), ``xcg``, ``sref``, ``cbarr``, ``fcm`` (the 282
            ``/SUPWH/`` words before the call), ``tcd`` (the 58 ``TCD``
            words), and ``state``: ``et`` (the 14 saved span stations,
            ``ET_DATA`` on the first call) and ``kinbd``, ``koutbd``.

    Returns:
        ``delcm`` (``WING(211..)``), ``fcm``, ``tcd``, ``flp`` (with
        ``ETA(1)``, ``ETA(5)`` as nudged off a station) and ``state`` as
        the next call will find it.

    Notes:
        Kept as executed: ``KOUNT`` and ``INDEX`` are set once before the
        deflection loop, so from the second deflection on the stations
        inboard of the flap are referenced to the outboard edge; the
        Figure 6.1.5.1-67A factors ``KK`` are computed only where UNUSED
        and then reused; ``DXCP`` keeps its previous value where the
        strip's lift is zero; and the routine overwrites one of its
        standard span stations ``ET`` with the flap's inboard edge, which
        later calls inherit.
    """
    s = {k: float(v) for k, v in data['surface'].items()}
    f = [0.0] + [float(v) for v in data['f']]
    flp = [0.0] + [float(v) for v in data['flp']]
    fcm = [0.0] + [float(v) for v in data['fcm']]
    tcd = [0.0] + [float(v) for v in data['tcd']]
    st = data['state']
    et = [0.0] + [float(v) for v in st['et']]
    kinbd, koutbd = int(st['kinbd']), int(st['koutbd'])
    clasec, tanc4, swstr = s['clasec'], s['tanc4'], s['swstr']
    bsto2, bo2, tante, tanle = s['bsto2'], s['bo2'], s['tante'], s['tanle']
    xw, cr, tapexp, arstar = s['x'], s['cr'], s['tapexp'], s['arstar']
    alpo, cmo = s['alpo'], s['cmo']
    xcg, sref, cbarr = (float(data['xcg']), float(data['sref']),
                        float(data['cbarr']))
    mach = float(data['mach'])
    beta = math.sqrt(1. - mach**2)
    ndelta = int(f[16] + 0.5)
    fcm[1] = math.atan(tanc4 / beta) * RAD
    sweepb = fcm[1]
    fcm[6] = cavg = swstr / (2. * bsto2)
    iftype = int(f[17] + .5)
    arg1 = (tante - tanle) * bsto2
    boc = []
    for quarter_index in range(4):
        cc = cr + _ETAG[quarter_index] * arg1
        boc.append(2. * beta * bsto2 / cc)
    for eta_slot in range(1, 13):
        if flp[1] == et[eta_slot]:
            flp[1] = flp[1] + .0001
        if flp[5] == et[eta_slot]:
            flp[5] = flp[5] - .0001
    eta1, eta5 = flp[1], flp[5]
    gd = calculate_gdelta(eta1, eta5, boc, sweepb, bool(data['asyfp']),
                          data.get('tail'))
    if data['asyfp']:
        boc = gd['boch']
        tcd[43:47] = gd['gdh']
    else:
        tcd[1:15], tcd[15:29], tcd[29:43] = gd['gd2'], gd['gd3'], gd['gd1']
    fcm[2:6] = boc
    gdi, gdo = [0.0] + tcd[1:15], [0.0] + tcd[15:29]
    et[13], et[14] = eta1, eta5
    alpdel = [0.0] * 11
    if f[19] != UNUSED:
        for deflection_index in range(1, ndelta + 1):
            alpdel[deflection_index] = (-f[18 + deflection_index] /
                                        (f[deflection_index] * clasec))
    else:
        nn = 0
        for deflection_index in range(1, ndelta + 1):
            total = 0.0
            for _ in range(4):
                nn += 1
                total += flp[149 + nn]
            alpdel[deflection_index] = total / 4.
    etak = [0.0] + fcm[7:21]
    gdinbd = [0.0] + fcm[35:49]
    gdoutb = [0.0] + fcm[49:63]
    kout = 1
    done = False
    for eta_slot in range(1, 12):
        etak[eta_slot], gdinbd[eta_slot], gdoutb[eta_slot] = (et[eta_slot],
                                                              gdi[eta_slot],
                                                              gdo[eta_slot])
        if not (eta1 > et[eta_slot] and eta1 < et[eta_slot + 1]):
            continue
        kinbd = eta_slot + 1
        etak[kinbd], gdinbd[kinbd], gdoutb[kinbd] = eta1, gdi[13], gdo[13]
        jj, shift = eta_slot + 2, 1
        while True:
            for strip_slot in range(jj, 15):
                etak[strip_slot], gdinbd[strip_slot], gdoutb[strip_slot] = (
                    et[strip_slot - shift], gdi[strip_slot - shift],
                    gdo[strip_slot - shift])
            if kout == 2:
                done = True
                break
            inboard_strip = eta_slot + 1
            et[inboard_strip] = eta1
            for outboard_slot in range(inboard_strip, 12):
                if eta5 > et[outboard_slot] and eta5 < et[outboard_slot + 1]:
                    koutbd = outboard_slot + 2
                    etak[koutbd] = eta5
                    kout = 2
                    gdinbd[koutbd], gdoutb[koutbd] = gdi[14], gdo[14]
                    jj, shift = outboard_slot + 3, 2
                    break
            else:
                break
        if done:
            break
    deln4 = flp[60]
    arg = (f[12] - f[13]) / (4. * deln4)
    arg7 = bsto2 * tanle
    ck = [0.0] * 15
    cfoc = [0.0] * 15
    xle = [0.0] * 15
    deltgd = [0.0] * 15
    for strip_slot in range(1, 15):
        ck[strip_slot] = cr + etak[strip_slot] * arg1
        if kinbd <= strip_slot <= koutbd:
            cfoc[strip_slot] = (
                f[12] - arg * (etak[strip_slot] - eta1)) / ck[strip_slot]
        xle[strip_slot] = xw + (bo2 - bsto2) * tanle + etak[strip_slot] * arg7
        deltgd[strip_slot] = gdoutb[strip_slot] - gdinbd[strip_slot]
    arg2 = (1. - tapexp) / (1. + tapexp)
    kount = 1
    reference_strip = kinbd
    nn = 0
    kk = [0.0] + fcm[101:115]
    cloald = [0.0] + fcm[21:35]
    dxcp = [0.0] + fcm[143:283]
    cl = [0.0] * 15
    xcp = [0.0] * 15
    delcmf = [0.0] * 15
    swepb = [0.0] * 15
    deltp = [0.0] * 15
    delcm = [float(v) for v in data.get('delcm', [0.0] * 10)]
    for deflection_index in range(1, ndelta + 1):
        delta = f[deflection_index]
        argz = abs(delta)
        for kx in (kinbd, koutbd):
            cl[kx] = (-4. * bsto2 * deltgd[kx] * alpdel[deflection_index] *
                      delta / (ck[kx] * RAD))
        dcmf = [0.0] * 15
        for strip_index in range(1, 14):
            ll = strip_index
            cl[strip_index] = (-4. * bsto2 * deltgd[strip_index] *
                               alpdel[deflection_index] * delta /
                               (ck[strip_index] * RAD))
            if deflection_index == 1:
                cloald[strip_index] = (-4. * bsto2 * deltgd[strip_index] /
                                       (ck[strip_index] * RAD))
            if cl[strip_index] != 0.0:
                if not (etak[strip_index] <= eta1 - 0.2 or
                        etak[strip_index] > eta5 + 0.2):
                    arg = cfoc[kinbd]
                if not (etak[strip_index] < eta1 or etak[strip_index] > eta5):
                    kount = 3
                    arg = cfoc[strip_index]
                if not etak[strip_index] <= eta5:
                    kount = 2
                    arg = cfoc[koutbd]
                if iftype in (1, 5):
                    xcpbi = _tl(_X1126B, _X2126B, _Y5126B, argz, arg)
                else:
                    xcpbi = 0.54
                arg3 = 4.0 * (xcpbi - 0.25) * arg2 / arstar
                swepbi = math.atan(tanc4 - arg3) * RAD
                deltpi = math.atan(math.tan(argz / RAD) /
                                   math.cos(swepbi / RAD)) * RAD
                if kount == 3:
                    deltp[strip_index], swepb[strip_index] = deltpi, swepbi
                    kp = 1
                else:
                    if kount == 2:
                        reference_strip = koutbd
                    if kk[strip_index] == UNUSED:
                        kk[strip_index] = interx(
                            1, _X5126A,
                            [abs(etak[strip_index]
                                 - etak[reference_strip])],
                            [13], _Y5126A, lind=13)
                    swepb[strip_index], deltp[strip_index] = swepbi, deltpi
                    if kount == 2:
                        ll = koutbd
                    kp = 2
                cossb2 = math.cos(swepb[strip_index] / RAD)**2
                if f[29] != UNUSED:
                    delcmf[strip_index] = f[28 + deflection_index]
                    xc = 0.25 + abs((kk[strip_index] * delcmf[strip_index] /
                                     cl[reference_strip] * cossb2)
                                    if kp == 2 else
                                    (delcmf[strip_index] * cossb2 /
                                     cl[strip_index]))
                elif iftype <= 1:
                    delcmf[strip_index] = _tl(_X1215B, _X2215B, _Y1215B,
                                                cfoc[ll], deltp[strip_index])
                    xc = 0.25 + abs((kk[strip_index] * delcmf[strip_index] /
                                     cl[reference_strip] * cossb2)
                                    if kp == 2 else
                                    (delcmf[strip_index] * cossb2 /
                                     cl[strip_index]))
                else:
                    xrefoc, xcpocp = .25, .44
                    if iftype == 5:
                        xcpocp = 0.5 - 0.25 * cfoc[ll]
                    cpoc = (-((f[48 + deflection_index] -
                               f[38 + deflection_index]) / (eta1 - eta5) *
                              (etak[strip_index] - eta5) -
                              f[48 + deflection_index]) / ck[strip_index])
                    if etak[strip_index] < eta1 or etak[strip_index] > eta5:
                        cpoc = 1.0
                    ind = 4 * (deflection_index - 1)
                    scl = (flp[110 + ind] + flp[111 + ind] + flp[112 + ind] +
                           flp[113 + ind]) / 4.
                    if iftype < 5:
                        delcmf[strip_index] = scl * (xrefoc - xcpocp * cpoc)
                    else:
                        if iftype == 6:
                            cpoc = 1.0
                        cmdlep = interx(1, _X12136, [cfoc[ll] / cpoc], [9],
                                        _Y12136, lind=9)
                        delcmf[strip_index] = (cmdlep * cpoc**2 * delta +
                                               (xrefoc + cpoc - 1.) * scl +
                                               .75 * clasec * alpo * cpoc *
                                               (cpoc - 1.) +
                                               (cpoc**2 - 1.) * cmo)
                    xc = 0.25 - ((kk[strip_index] * delcmf[strip_index] /
                                  cl[reference_strip] * cossb2)
                                 if kp == 2 else
                                 (delcmf[strip_index] * cossb2 /
                                  cl[strip_index]))
                xcp[strip_index] = xle[strip_index] + xc * ck[strip_index]
            nn += 1
            if cl[strip_index] != 0.0:
                dxcp[nn] = (xcg - xcp[strip_index]) / cbarr
            dcmf[strip_index] = (cl[strip_index] * dxcp[nn] * ck[strip_index] *
                                 swstr / (cavg * sref))
        dcmf[14] = 0.0
        nn += 1
        dxcp[nn] = dxcp[nn - 1]
        cloald[14] = 0.0
        delcm[deflection_index - 1] = float(trapz(dcmf[1:15], etak[1:15])[0])
    fcm[7:21] = etak[1:15]
    fcm[21:35] = cloald[1:15]
    fcm[35:49] = gdinbd[1:15]
    fcm[49:63] = gdoutb[1:15]
    fcm[63:73] = alpdel[1:11]
    fcm[73:87] = ck[1:15]
    fcm[87:101] = deltgd[1:15]
    fcm[101:115] = kk[1:15]
    fcm[115:129] = xle[1:15]
    fcm[129:143] = cfoc[1:15]
    fcm[143:283] = dxcp[1:141]
    return {'delcm': delcm, 'fcm': fcm[1:], 'tcd': tcd[1:],
            'flp': flp[1:], 'method': 'legacy_flapcm',
            'state': {'et': et[1:], 'kinbd': kinbd, 'koutbd': koutbd}}


def m37o45(asyfp: bool, control_type: float, mach: float, aht68: float,
           run_flapcm, tail: Mapping[str, float]):
    """M37O45: FLAPCM, or for an all-moving tail (``ASYFP`` and control
    type 5) GDELTA's tail loadings at the sweep ``atan(AHT(68)/beta)``.
    The GDELTA call passes the overlay's unset locals for the flap span,
    which the all-moving-tail path does not read."""
    if asyfp and control_type == 5.:
        beta = math.sqrt(1. - mach**2)
        sb = math.atan(aht68 / beta) * RAD
        return calculate_gdelta(0.0, 0.0, [1., 1., 1., 1.], sb, True, tail)
    return run_flapcm()
