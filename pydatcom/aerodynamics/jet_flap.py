"""
Jet-flap increments: JETFP and its overlay M55O67.

JETFP adds the lift, lift-curve slope, maximum lift and pitching moment of
a blown flap (Section 6.1.4 and 6.1.5.1): a pure jet flap, an internally
blown flap (IBF), an externally blown flap (EBF) or a combination, over
the flap deflections of ``/FLAPIN/``.  Section lift comes from Figure
6.1.1.1-49 at the section and three-dimensional momentum coefficients;
the wing increment from the IBF equation (also used for the pure jet
flap) or, for the EBF, the jet-turning factor of Figure 6.1.4.1-18; the
blown lift-curve slope from Figures 6.1.4.2-9 and 6.1.4.1-15; the EBF's
maximum-lift increment from Figure 6.1.4.3-12; and the moment from a
five-strip span breakdown with the centres of Figures 6.1.2.1-37 and
6.1.5.1-68.

The tables were extracted from the source DATA statements by parsing, with
``tools/fortran_data.py``.

Reference: datcom-legacy/datcom_2000/jetfp.f, m55o67.f
"""

import math
from typing import Dict, Mapping

import numpy as np

from pydatcom.utils.constants import DEG, PI, RAD, UNUSED
from pydatcom.utils.legacy_interp import interx
from pydatcom.utils.legacy_numeric import tbfunx

# Figure 6.1.1.1-49: section lift, against the momentum coefficient
# (9) and flap-chord ratio (7), padded to 18.
_X149 = np.array([
    0.0, 0.5, 1.0, 2.0, 3.0, 4.0, 6.0, 8.0, 10.0,
    0.0, 0.05, 0.1, 0.2, 0.3, 0.5, 1.0, 0.0, 0.0,
])
_Y149 = np.array([
    0.0, 2.8, 4.0, 6.0, 7.9, 9.7, 12.8, 15.6, 18.4,
    1.9, 3.5, 4.8, 6.8, 8.5, 10.2, 13.2, 16.1, 19.0,
    2.3, 3.9, 5.2, 7.2, 8.9, 10.6, 13.6, 16.5, 19.4,
    3.5, 5.0, 6.0, 8.0, 9.6, 11.3, 14.2, 17.3, 20.0,
    4.2, 5.5, 6.5, 8.5, 10.2, 11.9, 14.6, 17.8, 20.5,
    5.1, 6.2, 7.5, 9.1, 10.3, 12.6, 15.6, 18.6, 21.1,
    6.3, 7.5, 8.5, 10.5, 12.1, 13.6, 16.5, 19.2, 22.1,
])
# Figure 6.1.4.1-18: the EBF turning factor (12 by 4, padded to 24).
_X418 = np.array([
    0.0, 0.5, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0,
    0.1, 0.2, 0.3, 0.4, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
])
_Y418 = np.array([
    0.0, 1.7, 2.8, 4.8, 6.5, 8.1, 9.75, 11.1, 12.6, 14.05, 15.5, 16.9,
    0.0, 1.5, 2.6, 4.5, 6.2, 7.85, 9.3, 10.9, 12.25, 13.8, 15.1, 16.5,
    0.0, 1.45, 2.45, 4.3, 6.05, 7.75, 9.1, 10.75, 12.05, 13.55, 14.98, 16.4,
    0.0, 1.45, 2.45, 4.25, 5.95, 7.55, 9.0, 10.55, 11.95, 13.25, 14.8, 16.05,
])
# Figure 6.1.4.2-9: K(A, CJ') (11 by 10, padded to 22).
_X409 = np.array([
    0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0,
    1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 0.0,
])
_Y409 = np.array([
    1.0, 1.52, 2.02, 2.49, 2.9, 3.35, 3.75, 4.15, 4.55, 4.95, 5.35,
    1.0, 1.35, 1.65, 1.95, 2.21, 2.5, 2.75, 3.03, 3.28, 3.56, 3.8,
    1.0, 1.3, 1.56, 1.8, 2.03, 2.24, 2.46, 2.7, 2.9, 3.13, 3.35,
    1.0, 1.3, 1.53, 1.75, 1.95, 2.16, 2.35, 2.55, 2.75, 2.95, 3.15,
    1.0, 1.3, 1.53, 1.73, 1.93, 2.12, 2.3, 2.5, 2.67, 2.84, 3.07,
    1.0, 1.3, 1.53, 1.73, 1.92, 2.1, 2.28, 2.46, 2.65, 2.82, 3.03,
    1.0, 1.3, 1.53, 1.73, 1.92, 2.1, 2.28, 2.45, 2.63, 2.8, 2.98,
    1.0, 1.3, 1.53, 1.73, 1.92, 2.1, 2.28, 2.45, 2.63, 2.8, 2.98,
    1.0, 1.3, 1.53, 1.73, 1.92, 2.1, 2.28, 2.45, 2.63, 2.8, 2.98,
    1.0, 1.3, 1.53, 1.73, 1.92, 2.1, 2.28, 2.45, 2.63, 2.8, 2.98,
])
# Figure 6.1.4.1-15: the span factor K_b (11 by 3, padded to 22).
_X415 = np.array([
    0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0,
    0.0, 0.5, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
])
_Y415 = np.array([
    0.0, 0.17, 0.3, 0.43, 0.57, 0.67, 0.78, 0.86, 0.93, 0.98, 1.0,
    0.0, 0.14, 0.28, 0.4, 0.52, 0.52, 0.73, 0.82, 0.91, 0.97, 1.0,
    0.0, 0.13, 0.26, 0.38, 0.39, 0.5, 0.7, 0.8, 0.89, 0.96, 1.0,
])
# Figure 6.1.4.3-12: the maximum-lift increment.
_X412 = np.array([
    0.0, 0.4, 0.8, 1.2, 1.6, 2.0, 2.4,
])
_Y412 = np.array([
    0.0, 1.6, 2.65, 3.65, 4.55, 5.4, 6.2,
])
# Figure 6.1.2.1-37: the 2-D centre of lift (11 by 8, padded to 22).
_X237 = np.array([
    0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0,
    0.0, 0.05, 0.1, 0.2, 0.3, 0.5, 0.7, 1.0, 0.0, 0.0, 0.0,
])
_Y237 = np.array([
    0.5, 0.575, 0.61, 0.63, 0.655, 0.67, 0.685, 0.695, 0.705, 0.71, 0.715,
    0.49, 0.55, 0.58, 0.605, 0.62, 0.645, 0.66, 0.67, 0.675, 0.68, 0.69,
    0.48, 0.525, 0.56, 0.585, 0.6, 0.61, 0.62, 0.635, 0.64, 0.65, 0.655,
    0.445, 0.485, 0.51, 0.54, 0.555, 0.565, 0.575, 0.585, 0.59, 0.595, 0.599,
    0.42, 0.45, 0.475, 0.49, 0.5, 0.505, 0.51, 0.515, 0.52, 0.525, 0.53,
    0.38, 0.385, 0.39, 0.399, 0.402, 0.405, 0.41, 0.415, 0.42, 0.422, 0.426,
    0.325, 0.32, 0.315, 0.315, 0.315, 0.315, 0.315, 0.315, 0.315, 0.315, 0.315,
    0.25, 0.23, 0.215, 0.21, 0.202, 0.198, 0.195, 0.19, 0.185, 0.178, 0.175,
])
# Figure 6.1.5.1-68: the 3-D to 2-D centre fraction (11 by 11).
_X568 = np.array([
    0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0,
    0.0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5,
])
_Y568 = np.array([
    1.0, 1.098, 1.18, 1.255, 1.32, 1.38, 1.43, 1.48, 1.52, 1.575, 1.605,
    1.0, 1.095, 1.175, 1.245, 1.31, 1.375, 1.415, 1.46, 1.5, 1.55, 1.58,
    1.0, 1.09, 1.17, 1.235, 1.3, 1.365, 1.4, 1.445, 1.48, 1.52, 1.57,
    1.0, 1.085, 1.16, 1.225, 1.29, 1.345, 1.385, 1.43, 1.46, 1.5, 1.54,
    1.0, 1.08, 1.15, 1.215, 1.28, 1.325, 1.375, 1.4, 1.445, 1.48, 1.505,
    1.0, 1.075, 1.14, 1.205, 1.26, 1.305, 1.355, 1.385, 1.42, 1.455, 1.48,
    1.0, 1.07, 1.13, 1.195, 1.24, 1.29, 1.33, 1.36, 1.39, 1.425, 1.46,
    1.0, 1.065, 1.12, 1.18, 1.22, 1.265, 1.3, 1.345, 1.37, 1.39, 1.42,
    1.0, 1.06, 1.11, 1.16, 1.2, 1.24, 1.28, 1.3, 1.33, 1.36, 1.38,
    1.0, 1.055, 1.1, 1.145, 1.18, 1.21, 1.245, 1.27, 1.29, 1.315, 1.34,
    1.0, 1.05, 1.09, 1.12, 1.16, 1.185, 1.205, 1.23, 1.25, 1.275, 1.29,
])


def _fig(table, dep, lengths, lind, v1, v2):
    return interx(2, table, [v1, v2], lengths, dep, lind=lind,
                  lx1l=1, lx2l=1, lx1u=1, lx2u=1)


def _f149(v1, v2):
    return _fig(_X149, _Y149, [9, 7], 9, v1, v2)


def _f415(v1, v2):
    return _fig(_X415, _Y415, [11, 3], 11, v1, v2)


def _f237(v1, v2):
    return _fig(_X237, _Y237, [11, 8], 11, v1, v2)


def _f568(v1, v2):
    return _fig(_X568, _Y568, [11, 11], 11, v1, v2)


def calculate_jetfp(data: Mapping[str, object]) -> Dict[str, object]:
    """Translate JETFP: jet-flap lift, slope, maximum lift and moment.

    Args:
        data: By name: ``f`` (the ``/FLAPIN/`` words by index: 1.. the flap
            deflections, 12-15 the flap chords and spans, 16 count, 17
            type, 39.., 49.. the extended chords, 63 ``CMU``, 64.. the jet
            deflections, 74 the jet-flap type, 75.. the jet turning
            angles), ``jet`` (``aietlj``, ``jevloc``, ``jealoc``,
            ``jeangl``, ``jelloc``, ``jerad``), ``a`` (``A`` 4, 38, 118,
            120), ``win`` (``WINGIN`` 1, 3, 4, 6, 11, 16), ``claub``
            (``WING(101..)``), ``position`` (``xcg``, ``xw``, ``zw``,
            ``zcg``), ``sref``, ``cbarr``, and ``stale``: ``etat`` (the
            saved turning efficiency) and the ``WING`` curves ``deccl``,
            ``delcm``, ``dclmax``, ``clab`` as they stood.

    Returns:
        ``deccl`` (``WING(201..)``), ``delcm`` (``(211..)``), ``dclmax``
        (``(221..)``), ``clab`` (``(241..)``), ``jeangl`` (12 degrees
        when UNUSED), ``etat``, ``skipped`` (no momentum coefficient or
        jet-flap type) and ``terminated`` (an EBF whose jet misses the
        flap: the source prints a message and stops).

    Notes:
        Kept as executed: the combination jet flap (type 4) forms its
        section increments but never stores a wing increment; the pure
        jet flap and IBF read the turning efficiency ``ETAT`` that only
        the EBF sets.
    """
    f = {int(k): float(v) for k, v in data['f'].items()}
    j = {k: float(v) for k, v in data['jet'].items()}
    a = {int(k): float(v) for k, v in data['a'].items()}
    w = {int(k): float(v) for k, v in data['win'].items()}
    p = {k: float(v) for k, v in data['position'].items()}
    st = data['stale']
    out = {key: [float(v) for v in st[key]]
           for key in ('deccl', 'delcm', 'dclmax', 'clab')}
    out.update({'jeangl': j['jeangl'], 'etat': float(st['etat']),
                'skipped': False, 'terminated': False,
                'method': 'legacy_jetfp'})
    cmu = f[63]
    if cmu == UNUSED or f[74] == UNUSED:
        out['skipped'] = True
        return out
    sref, cbarr = float(data['sref']), float(data['cbarr'])
    sw, tanle, tapr, ar = a[4], a[38], a[118], a[120]
    chrdtp, sspn, chrdr = w[1], w[4], w[6]
    twista, toc = w[11], w[16]
    chrdfi, chrdfo, spanfi, spanfo = f[12], f[13], f[14], f[15]
    cf = (chrdfi + chrdfo) / 2.
    if j['jeangl'] == UNUSED:
        j['jeangl'] = 12.0
    out['jeangl'] = j['jeangl']
    tanjet = math.tan(j['jeangl'] * DEG)
    c = chrdr + (spanfi + spanfo) * (chrdtp - chrdr) / (2.0 * sspn)
    jetflp = int(f[74] + .5)
    ndelta = int(f[16] + .5)
    ftype = int(f[17] + .5)
    etat = out['etat']
    claub = [float(v) for v in data['claub']]
    for n in range(ndelta):
        delflp, deljet, effjet = f[1 + n], f[64 + n], f[75 + n]
        cprmei, cprmeo = f[39 + n], f[49 + n]
        cprime = (cprmei + cprmeo) / 2.
        plain = ftype in (5, 1, 6)
        if plain:
            cprime = c
        tcprm = toc * c / cprime * 0.80
        st_ = sw if plain else sw + 2.0 * (cprime - c) * (spanfo - spanfi)
        at = ar * sw / st_
        yi, yo = spanfi, spanfo
        if jetflp == 3:
            cosdf = math.cos(delflp * DEG) - 1.
            bk = (cprmei - cprmeo) / (spanfi - spanfo)
            fk = (chrdfi - chrdfo) / (spanfi - spanfo)
            yi = (j['jelloc'] - j['jerad'] - tanjet *
                  (p['xw'] - j['jealoc'] + cprmeo - bk * spanfo -
                   fk * cosdf * spanfo + chrdfo)) / \
                (1. + tanjet * (tanle - bk - fk * cosdf))
            yo = 2. * j['jelloc'] - yi
            if yi < spanfi:
                yi = spanfi
            if yo > spanfo:
                yo = spanfo
        sj = (yo - yi) * chrdr * (2. - (yi + yo) * (1. - tapr) / sspn)
        mu = cmu * c / cprime
        if jetflp not in (2, 3):
            cldj = _f149(mu, 0.0)
            delcl = ((1. + tcprm) * deljet * (cldj - mu) +
                     mu * deljet) * cprime / (c * RAD)
            delcl1 = delcl if jetflp == 4 else None
        if jetflp != 1:
            cldf = _f149(mu, cf / cprime)
            delcl = ((1. + tcprm) * delflp * (cldf - mu) +
                     mu * delflp) * cprime / (c * RAD)
            if jetflp == 4:
                delcl += delcl1
        cj = 2. / sref * c * cmu * (spanfo - spanfi)
        cjprm = cj * sref / sj
        if jetflp not in (2, 3):
            cldjt = _f149(cjprm, 0.0)
            delclt = ((1. + tcprm) * deljet * (cldjt - cjprm) +
                      cjprm * deljet) * cprime / (c * RAD)
            delcl2 = delclt if jetflp == 4 else None
        if jetflp != 1:
            cldft = _f149(cjprm, cf / cprime)
            delclt = ((1. + tcprm) * delflp * (cldft - cjprm) +
                      cjprm * delflp) * cprime / (c * RAD)
            if jetflp == 4:
                delclt += delcl2
        if jetflp in (1, 2):
            out['deccl'][n] = (delclt * sj / sref * (at + 2. * cjprm / PI) /
                               (at + 2. + .604 * math.sqrt(cjprm) +
                                .876 * cjprm))
        if jetflp == 3:
            claprm = _f149(cjprm, 1.)
            pido4 = _fig(_X418, _Y418, [12, 4], 12, cjprm, cf / cprime)
            out['deccl'][n] = (pido4 * effjet * sj / (sref * RAD) *
                               (PI * at + 2. * cjprm) /
                               (PI * at + claprm + 2.01 * cjprm))
        if jetflp == 4:
            continue
        term2 = cj * (math.cos((effjet + deljet) * DEG) - 1.) / RAD
        atcjk = _fig(_X409, _Y409, [11, 10], 11, cjprm, at)
        bk = _f415(spanfo / sspn, 1.) - _f415(spanfi / sspn, 1.)
        term1k = (atcjk - 1.) * bk + 1.
        out['clab'][n] = term1k * claub[n] + term2
        if jetflp == 3:
            etat = 0.0
            if ftype == 4:
                etat = 1. - 2. * effjet / 300.
                cons = etat * cj * math.sin(effjet * DEG)
                out['dclmax'][n] = tbfunx(_X412, _Y412, cons, 1, 1)[0]
        if ftype > 5:
            continue
        if jetflp == 3:
            eh = (j['jerad'] + tanjet * (cprime - cf + p['xw'] +
                                         j['jelloc'] * tanle)) * \
                math.cos(j['aietlj'] * DEG)
            if eh < j['jevloc'] - p['zw']:
                out['terminated'] = True
                break
        eta1, eta2, eta3, eta4 = (spanfi / sspn, yi / sspn, yo / sspn,
                                  spanfo / sspn)
        at = ar * c / cprime
        cmuprm = cj * sref / sj * c / cprime
        ak2 = at / (2. + at) * cprime / c
        ak3 = ((at + 2. * cmuprm / PI) /
               (at + 2. + .604 * math.sqrt(cmuprm) + .876 * cmuprm) *
               cprime / c)
        ak1, ak4, ak5 = 0., ak2, 0.
        bk10 = _f415(eta1, tapr)
        bk22 = _f415(eta2, tapr)
        bk21 = bk22 - bk10
        bk33 = _f415(eta3, tapr)
        bk32 = bk33 - bk22
        bk44 = _f415(eta4, tapr)
        bk43 = bk44 - bk33
        bk54 = 1. - bk44
        if tapr != 1.:
            xmoc = (1. + tapr) * ar * tanle / (4. * (1. - tapr))
        else:
            xmoc = (p['xw'] - p['xcg']) / chrdr
        xf2f = _f237(cmuprm, cf / cprime)
        x3d2f = _f568(1. / at, cf / cprime)
        xf2j = _f237(cmuprm, 0.)
        x3d2j = _f568(1. / at, 0.)
        al = twista
        xficdd = xf2f * x3d2f
        xjicdd = xf2j * x3d2j
        xacdd = _f237(cmuprm, 1.0)
        cla1 = _f149(0., 1.)
        cla2 = cla4 = cla5 = cla1
        cla3 = _f149(cmuprm, 1.)
        delc41 = al * ak1 * cla1 / RAD
        delc42 = al * ak2 * cla2 / RAD
        delc43 = al * ak3 * cla3 / RAD
        delc44 = al * ak4 * cla4 / RAD
        delc45 = al * ak5 * cla5 / RAD
        x23 = xmoc - xacdd * cprime / c
        delcm4 = -cj * etat * sw / sj * xmoc * al / RAD
        dcmda1 = delc41 * xmoc
        dcmda2 = delc42 * xmoc
        dcmda3 = delc43 * x23 + delcm4
        dcmda4 = delc44 * xmoc
        dcmda5 = delc45 * xmoc
        x5 = xmoc - xficdd * cprime / c
        cldfi = _f149(cmuprm, cf / cprime)
        delc5 = delflp * ak3 * cldfi / RAD
        dcmdf = delc5 * x5
        cldji = _f149(cmuprm, 0.)
        dj = effjet - delflp if effjet > delflp else 0.0
        if jetflp < 5 and jetflp != 3:
            dj = deljet
        delc6 = dj * ak3 * cldji / RAD
        x6 = xmoc - xjicdd * cprime / c
        dcmdj = delc6 * x6
        cmm = (dcmdj * bk32 + dcmdf * (bk44 - bk10) + dcmda1 * bk10 +
               dcmda2 * bk21 + dcmda3 * bk32 + dcmda4 * bk43 +
               dcmda5 * bk54)
        cl1 = bk10 * delc41
        cl2 = bk21 * delc42
        cl3 = bk32 * (delc43 + delc5 + delc6)
        cl4 = bk43 * delc44
        cl5 = bk54 * delc45
        dxocb = (p['xw'] + xmoc * chrdr - p['xcg']) / cbarr
        out['delcm'][n] = (cmm + etat * cj * (p['zcg'] - p['zw']) / c +
                           dxocb * (-cl1 - cl2 - cl4 - cl5 - cl3 +
                                    cj * etat * sw / sj * al * sj /
                                    (RAD * sw)))
    out['etat'] = etat
    return out
