"""
Jet power effects: FG6115, JETPWE and their overlay M30O36.

JETPWE adds the jet's lift, pitching moment and tail downwash to the
subsonic buildup: the thrust vector's lift, the inlet normal force
(through the wing upwash of Figure 4.4.1-73), the offset thrust line's
moment, and the jet-induced downwash at the tail from Figure 4.6.1-28 and
either Figure 4.6.1-30 (FG6115, a subsonic jet near the tail), Figure
4.6.1-31 (a subsonic jet far upstream) or Figure 4.6.1-32 (a supersonic
jet, moved to its fully expanded orifice).  M30O36 adds the body-axis
forces and the lift and moment slopes.

The tables were extracted from the source DATA statements by parsing, with
``tools/fortran_data.py``.

Reference: datcom-legacy/datcom_2000/fg6115.f, jetpwe.f, m30o36.f
"""

import math
from typing import Dict, Mapping

import numpy as np

from pydatcom.utils.constants import DEG, RAD, UNUSED
from pydatcom.utils.legacy_numeric import tbfunx, tlinvs
from pydatcom.utils.legacy_tables import tlinex

# FG6115: Figure 4.6.1-30A-C, the jet downwash increment, as TLINVS
# tables: the increment (X1), x/r_j (X2) and z/r_j stored X2-fastest.
_D6115A = np.array([
    0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.1, 0.12, 0.15,
])
_X6115A = np.array([
    -4.5, -4.0, 0.0, 4.0, 8.0, 12.0, 16.0, 20.0, 24.0, 28.0,
    30.0,
])
_Z6115A = np.array([
    2.0, 6.0, 10.8, 13.4, 15.2, 15.9, 15.9, 15.4, 14.8, 13.7, 13.2,
    0.0, 3.0, 7.6, 9.5, 10.9, 11.6, 11.5, 10.9, 10.0, 9.0, 8.6,
    0.0, 2.0, 6.2, 7.8, 8.8, 9.2, 9.0, 8.5, 8.0, 7.2, 6.6,
    0.0, 2.0, 4.5, 6.5, 7.5, 8.0, 7.8, 7.4, 6.8, 6.0, 5.5,
    0.0, 2.0, 4.1, 6.0, 6.8, 7.2, 7.1, 6.6, 5.9, 5.1, 4.9,
    0.0, 0.0, 3.3, 5.2, 6.2, 6.5, 6.3, 5.8, 5.2, 5.4, 4.0,
    0.0, 0.0, 2.6, 4.5, 5.7, 6.0, 5.9, 5.5, 5.25, 3.5, 3.2,
    0.0, 0.0, 2.0, 3.5, 4.7, 5.2, 5.0, 4.2, 3.4, 2.3, 2.3,
    0.0, 0.0, 1.0, 3.0, 4.0, 4.3, 4.0, 3.1, 2.2, 2.25, 2.3,
    0.0, 0.0, 1.0, 2.3, 3.0, 3.1, 2.8, 2.1, 2.2, 2.25, 2.3,
])
_D6115B = np.array([
    0.1, 0.15, 0.2, 0.3, 0.4, 0.6, 0.8, 1.0, 1.5, 2.0,
])
_X6115B = np.array([
    -6.0, -4.0, -2.0, 0.0, 4.0, 8.0, 12.0, 14.0, 16.0, 20.0,
    24.0, 28.0, 30.0,
])
_Z6115B = np.array([
    6.0, 13.0, 15.5, 17.0, 19.3, 20.9, 21.8, 22.0, 22.0, 21.5, 20.7, 20.0, 19.6,
    0.0, 5.5, 10.0, 12.0, 14.7, 16.0, 16.48, 16.5, 16.4, 16.0, 15.3, 14.2, 13.8,
    0.0, 0.0, 6.9, 7.5, 12.0, 13.2, 13.6, 13.5, 13.4, 12.9, 12.0, 11.0, 10.5,
    0.0, 0.0, 3.1, 6.5, 9.1, 10.0, 10.2, 10.1, 9.95, 9.4, 8.3, 7.3, 7.0,
    0.0, 0.0, 1.0, 4.0, 7.3, 8.0, 8.0, 7.9, 7.7, 7.0, 6.3, 5.6, 5.25,
    0.0, 0.0, 1.0, 3.2, 5.0, 5.8, 5.7, 5.5, 5.25, 4.8, 4.0, 3.8, 3.9,
    0.0, 0.0, 0.0, 2.0, 4.0, 4.6, 4.5, 4.2, 4.0, 3.5, 3.4, 3.8, 3.9,
    0.0, 0.0, 0.0, 2.0, 3.5, 3.75, 3.4, 3.2, 3.0, 2.9, 3.4, 3.8, 3.9,
    0.0, 0.0, 0.0, 1.5, 2.6, 2.8, 2.6, 2.5, 2.5, 2.9, 3.4, 3.8, 3.9,
    0.0, 0.0, 0.0, 1.4, 1.8, 2.0, 2.3, 2.4, 2.5, 2.9, 3.4, 3.8, 3.9,
])
_D6115C = np.array([
    0.25, 0.5, 0.75, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0,
    6.0,
])
_X6115C = np.array([
    -8.0, -4.0, -2.0, 0.0, 4.0, 8.0, 14.0, 20.0, 24.0, 28.0,
    30.0,
])
_Z6115C = np.array([
    14.0, 20.5, 22.7, 24.5, 27.0, 28.0, 28.0, 28.0, 28.0, 28.0, 28.0,
    0.0, 6.0, 11.5, 14.0, 16.4, 17.9, 19.0, 19.0, 18.6, 18.0, 17.6,
    0.0, 0.0, 6.0, 7.6, 12.1, 13.3, 14.0, 13.7, 13.0, 12.2, 11.8,
    0.0, 0.0, 4.0, 7.3, 10.0, 11.0, 11.0, 10.5, 10.0, 9.4, 9.0,
    0.0, 0.0, 1.0, 4.5, 7.2, 8.1, 8.2, 7.7, 7.0, 6.5, 6.0,
    0.0, 0.0, 1.5, 3.4, 4.8, 6.5, 6.3, 5.8, 5.4, 4.8, 5.0,
    0.0, 0.0, 1.5, 3.5, 4.9, 5.5, 5.2, 4.8, 4.3, 4.5, 5.0,
    0.0, 0.0, 0.0, 2.5, 4.2, 5.0, 4.3, 3.9, 4.3, 4.5, 5.0,
    0.0, 0.0, 0.0, 1.9, 3.5, 3.6, 3.3, 3.9, 4.3, 4.5, 5.0,
    0.0, 0.0, 0.0, 1.0, 2.9, 3.1, 3.0, 3.9, 4.3, 4.5, 5.0,
    0.0, 0.0, 0.0, 1.0, 2.5, 2.5, 3.0, 3.9, 4.3, 4.5, 5.0,
])
# JETPWE: Figure 4.4.1-73, the upwash gradient at the plane of
# symmetry, against aspect ratio (X1) and chord distance (X2).
_X14161 = np.array([
    4.0, 6.0, 9.0, 12.0,
])
_X24161 = np.array([
    0.25, 0.4, 0.5, 0.6, 0.7, 0.8, 1.0, 1.2, 1.6, 2.0,
])
_Y44161 = np.array([
    1.08, 0.545, 0.4, 0.31, 0.24, 0.2, 0.13, 0.1, 0.06, 0.04,
    1.18, 0.68, 0.52, 0.4, 0.32, 0.27, 0.19, 0.15, 0.1, 0.08,
    1.3, 0.81, 0.62, 0.49, 0.4, 0.34, 0.25, 0.2, 0.13, 0.12,
    1.4, 0.88, 0.7, 0.56, 0.445, 0.39, 0.3, 0.24, 0.165, 0.14,
])
# Figure 4.6.1-28: the mean jet downwash over the tail, against
# spanwise (X1) and vertical (X2) position.
_X16114 = np.array([
    0.0, 0.3, 0.5, 0.7, 0.9, 1.1,
])
_X26114 = np.array([
    0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4,
])
_Y46114 = np.array([
    0.0, 0.54, 0.765, 0.86, 0.915, 0.95, 0.969, 0.988,
    0.0, 0.51, 0.721, 0.825, 0.89, 0.923, 0.95, 0.97,
    0.0, 0.45, 0.675, 0.79, 0.86, 0.9, 0.922, 0.945,
    0.0, 0.395, 0.595, 0.725, 0.81, 0.861, 0.892, 0.914,
    0.0, 0.3, 0.522, 0.66, 0.748, 0.81, 0.859, 0.89,
    0.0, 0.2, 0.415, 0.57, 0.682, 0.758, 0.81, 0.85,
])
# Figure 4.6.1-29: equivalent against actual velocity ratio, with the
# temperature ratio (X1, descending) and velocity ratio (X2).
_X16117 = np.array([
    0.46, 0.42, 0.38, 0.34, 0.32, 0.3, 0.28, 0.26, 0.24, 0.22, 0.2, 0.16, 0.12,
])
_X26117 = np.array([
    1.0, 6.0, 14.0, 24.0,
])
_Y46117 = np.array([
    1.0, 4.25, 9.66, 16.6,
    1.0, 4.1, 9.21, 15.7,
    1.0, 3.9, 8.8, 14.9,
    1.0, 3.85, 8.35, 14.1,
    1.0, 3.6, 8.1, 13.7,
    1.0, 3.55, 7.89, 13.3,
    1.0, 3.4, 7.6, 12.9,
    1.0, 3.3, 7.4, 12.4,
    1.0, 3.2, 7.1, 12.0,
    1.0, 3.1, 6.8, 11.4,
    1.0, 3.0, 6.5, 11.0,
    1.0, 2.72, 5.9, 9.85,
    1.0, 2.42, 5.15, 8.55,
])
# Figure 4.6.1-31: z/x times the downwash increment, far downstream.
_X6118A = np.array([
    0.0, 0.4, 0.8, 1.2, 1.6, 2.0, 2.4, 2.6,
])
_Y6118A = np.array([
    0.0, 0.34, 0.63, 0.88, 1.1, 1.31, 1.51, 1.61,
])
# Figure 4.6.1-32A/B: the equivalent orifice radius and its
# downstream displacement for a supersonic jet.
_X6118B = np.array([
    3.0, 4.0, 6.0, 8.0, 10.0, 12.0, 14.0, 16.0,
])
_Y6118B = np.array([
    1.08, 1.3, 1.7, 2.02, 2.23, 2.41, 2.55, 2.7,
])
_X6119A = np.array([
    1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6,
])
_Y6119A = np.array([
    0.0, 1.0, 1.6, 2.07, 2.4, 2.65, 2.86, 3.04, 3.16,
])


def _grid(flat, n1, n2):
    """A source ``Y(N2,N1)`` table as TLINEX's ``(n2, n1)`` array."""
    return np.asarray(flat, dtype=float).reshape(n1, n2).T


_FG_A = (_D6115A, _X6115A, _grid(_Z6115A, 10, 11))
_FG_B = (_D6115B, _X6115B, _grid(_Z6115B, 10, 13))
_FG_C = (_D6115C, _X6115C, _grid(_Z6115C, 11, 11))


def fg6115(xeporj: float, zjporj: float, vjpovi: float) -> float:
    """Translate FG6115: Figure 4.6.1-30, the downwash increment of a
    subsonic jet in a subsonic stream.

    Part A up to an equivalent velocity ratio of 2, part C from 8, part B
    within 0.01 of 4, and linear blends of A and B or B and C between; each
    part is inverted by TLINVS at ``x/r_j`` for the ``z/r_j`` given.
    """
    def part(table):
        d, x, z = table
        return tlinvs(d, x, z, xeporj, zjporj)

    if vjpovi <= 2.0:
        return part(_FG_A)
    if vjpovi >= 8.0:
        return part(_FG_C)
    vtest = vjpovi - 4.0
    if abs(vtest) < 1.e-2:
        return part(_FG_B)
    if vtest < 0.0:
        da, db = part(_FG_A), part(_FG_B)
        return da + (db - da) * (vjpovi - 2.) / (4. - 2.)
    db, dc = part(_FG_B), part(_FG_C)
    return db + (dc - db) * (vjpovi - 4.) / (8. - 4.)


def _tmpink(t: float) -> float:
    return 273. + .556 * (t - 492.)


def calculate_jetpwe(alpha_deg, mach: float, jet: Mapping[str, float],
                     win: Mapping[int, float], a: Mapping[int, float],
                     position: Mapping[str, float], htpl: bool,
                     tail: Mapping[str, object], sref: float, cbarr: float,
                     stale: Mapping[str, float]) -> Dict[str, object]:
    """Translate JETPWE: jet power effects on lift and pitching moment.

    Args:
        alpha_deg: ``FLC(23..)``; mach: ``MACH(M)``.
        jet: The ``/POWER/`` words by name: ``aietlj``, ``nengsj``,
            ``thstcj``, ``jialoc``, ``jevloc``, ``jealoc``, ``jinlta``,
            ``jeangl``, ``jevelo``, ``ambtmp``, ``jestmp``, ``jelloc``,
            ``jetotp``, ``ambstp``, ``jerad``.
        win: ``WINGIN`` 1, 2, 4, 5, 6.
        a: ``A`` words 38, 86, 120, 134.
        position: ``xcg``, ``xw``, ``aliw``, ``zcg``, ``xh``, ``zh``,
            ``alih``.
        htpl: The tail is present.
        tail: ``span`` (``HTIN(4)``), ``cla`` (``HT(101)``), ``xbarr``
            (``AHT(161)``), ``q`` (``DWASH(1..)``, the dynamic-pressure
            ratio).
        sref, cbarr: ``/OPTION/``.
        stale: Words read before being set: ``cosaih`` (``COSAIH``, a
            saved local, when there is no tail) and ``comp1``
            (``/POWR/``'s ``COMP1``, likewise).

    Returns:
        The ``/POWR/`` words set, by name (``deuda``, the tail words
        ``xep`` ... ``xhpc``, ``zbart``, ``xl``, ``dlh``, ``comp1``, the
        last angle's ``atp``, ``epslon``, ``atj``, and the curves ``dclt``,
        ``dclnj``, ``dclhe``, ``dcmnj``, ``dcme``, ``dcmt1``), the
        ``/IPOWER/`` curves ``dcdpow``, ``dclpow``, ``dcmpow``, ``kase``,
        and ``win1``, ``WINGIN(1)`` as the routine leaves it.

    Notes:
        Kept as executed: for a jet on the inboard panel of a cranked wing
        the routine sets ``WINGIN(1)`` to the break chord and keeps it (the
        ``SAVE`` meant to restore the tip chord is taken after the
        overwrite); with no tail the moment arm ``DLH`` uses the previous
        call's tail-incidence cosine.
    """
    j = {k: float(v) for k, v in jet.items()}
    w = {int(k): float(v) for k, v in win.items()}
    g = {int(k): float(v) for k, v in a.items()}
    p = {k: float(v) for k, v in position.items()}
    out: Dict[str, object] = {}
    tpcne = j['nengsj'] * j['thstcj']
    cosaiw = math.cos(DEG * p['aliw'])
    jelloc = j['jelloc']
    win1 = w[1]
    cranked = j['nengsj'] != 1. and w[2] != UNUSED
    if cranked and jelloc > w[4] - w[2]:
        crj = w[5] - (w[5] - w[1]) * ((jelloc - (w[4] - w[2])) / w[2])
        xbarj = (p['xw'] + g[38] * (w[4] - w[2]) + g[86] *
                 (jelloc - (w[4] - w[2])) - j['jialoc']) * cosaiw / crj
    else:
        if cranked:
            win1 = w[5]
        crj = w[6] - (w[6] - win1) * jelloc / w[4]
        xbarj = (p['xw'] + g[38] * jelloc + crj / 4.0 - j['jialoc']) * \
            cosaiw / crj
    out['win1'] = win1
    deuda = tlinex(_X14161, _X24161, _grid(_Y44161, 4, 10), g[120], xbarj,
                   2, 2, 2, 2)
    ainne2 = 2.0 * j['jinlta'] * j['nengsj'] / sref
    sinait = math.sin(DEG * j['aietlj'])
    cosait = math.cos(DEG * j['aietlj'])
    cosaih = float(stale['cosaih'])
    comp1 = float(stale['comp1'])
    kase = None
    if htpl:
        cosaih = math.cos(DEG * p['alih'])
        xep = p['xh'] + float(tail['xbarr']) * cosaih - j['jealoc']
        tanthj = math.tan(DEG * j['jeangl'])
        zt = j['jevloc'] + j['jerad'] * sinait / tanthj
        zjp = xep * sinait + (p['zh'] - zt) * cosait
        xep = xep / cosait
        xjp = 4.6 * j['jerad']
        xhp = xjp + xep
        ain = 49.0 * math.sqrt(j['ambtmp'])
        vin = mach * ain
        vjovin = j['jevelo'] / vin
        xeporj = xep / j['jerad']
        bo2h = float(tail['span'])
        zjpobh = zjp / (2.0 * bo2h)
        ytob2h = jelloc / bo2h
        fig28 = _grid(_Y46114, 6, 8)
        debode = tlinex(_X16114, _X26114, fig28, ytob2h, zjpobh, 0, 0, 2, 2)
        out.update({'xep': xep, 'zjp': zjp, 'xjp': xjp, 'xhp': xhp,
                    'ain': ain, 'vin': vin, 'zjpobh': zjpobh,
                    'ytob2h': ytob2h})
        kase = 3 if j['jevelo'] >= ain else (2 if xeporj > 16.0 else 1)
        if kase == 2:
            zjpxhp = zjp / xhp
            srtpco = sref * j['thstcj'] / xhp**2
            zjdexh = tbfunx(_X6118A, _Y6118A, srtpco, 0, 1)[0]
            de = zjdexh / zjpxhp
            out.update({'zjpxhp': zjpxhp, 'srtpco': srtpco,
                        'zjdexh': zjdexh})
        else:
            zjporj = None
            if kase == 3:
                pteopi = j['jetotp'] / j['ambstp']
                rjporj = tbfunx(_X6118B, _Y6118B, pteopi, 2, 2)[0]
                rjp = j['jerad'] * rjporj
                dxporj = tbfunx(_X6119A, _Y6119A, rjporj, 0, 2)[0]
                dxp = dxporj * j['jerad']
                xepc = xep + dxp
                xhpc = xepc + rjp / tanthj
                ztp = j['jevloc'] - dxp * sinait
                zjp = xhpc * sinait + cosait * (p['zh'] - ztp)
                zjpobh = zjp / (2.0 * bo2h)
                debode = tlinex(_X16114, _X26114, fig28, ytob2h, zjpobh,
                                0, 0, 2, 2)
                zjporj = zjp / rjp
                xeporj = xepc / rjp
                out.update({'pteopi': pteopi, 'rjporj': rjporj, 'rjp': rjp,
                            'dxporj': dxporj, 'dxp': dxp, 'xepc': xepc,
                            'xhpc': xhpc, 'ztp': ztp, 'zjp': zjp,
                            'zjpobh': zjpobh})
            tinotj = _tmpink(j['ambtmp']) / _tmpink(j['jestmp'])
            vjpovi = tlinex(_X16117, _X26117, _grid(_Y46117, 13, 4),
                            tinotj, vjovin, 2, 0, 2, 1)
            if kase == 1:
                zjporj = zjp / j['jerad']
            de = fg6115(xeporj, zjporj, vjpovi)
            out.update({'tinotj': tinotj, 'vjpovi': vjpovi,
                        'zjporj': zjporj})
        comp1 = float(tail['cla']) * de * debode
        out.update({'de': de, 'debode': debode})
    zbart = cosait * (j['jevloc'] - p['zcg']) + sinait * \
        (j['jealoc'] - p['xcg'])
    dcmt1 = -j['thstcj'] * zbart * j['nengsj'] / cbarr
    xl = j['jialoc'] - p['xcg']
    dlh = p['xh'] + float(tail['xbarr']) * cosaih - p['xcg']
    curves = {k: [] for k in ('dclt', 'dclnj', 'dclhe', 'dclpow', 'dcmnj',
                              'dcme', 'dcmpow', 'dcdpow')}
    for jj, alpha in enumerate(alpha_deg):
        atp = alpha + j['aietlj']
        curves['dclt'].append(tpcne * math.sin(DEG * atp))
        epslon = deuda * (p['aliw'] + alpha - g[134])
        atj = alpha + j['aietlj'] + epslon
        curves['dclnj'].append(ainne2 * math.sin(DEG * atj))
        dclhe = comp1 * float(tail['q'][jj]) if htpl else 0.0
        curves['dclhe'].append(dclhe)
        curves['dclpow'].append(curves['dclnj'][-1] + curves['dclt'][-1] +
                                dclhe)
        curves['dcmnj'].append(-curves['dclnj'][-1] * xl / cbarr)
        curves['dcme'].append(-dclhe * dlh / cbarr)
        curves['dcmpow'].append(dcmt1 + curves['dcmnj'][-1] +
                                curves['dcme'][-1])
        curves['dcdpow'].append(0.0)
    out.update(curves)
    out.update({'deuda': deuda, 'comp1': comp1, 'zbart': zbart,
                'dcmt1': dcmt1, 'xl': xl, 'dlh': dlh, 'atp': atp,
                'epslon': epslon, 'atj': atj, 'kase': kase,
                'xbarj': xbarj, 'cosaih': cosaih,
                'method': 'legacy_jetpwe'})
    return out


def m30o36_block(alpha_deg, dcd, dcl, dcm) -> Dict[str, list]:
    """M30O36's pass over ``/IPOWER/``: CN and CA from the lift and drag
    increments (``POWER(61..)``, ``(81..)``) and TBFUNX's slopes of lift
    and moment (``(101..)``, ``(121..)``)."""
    alpha = np.asarray(alpha_deg, dtype=float)
    cn, ca, cla, cma = [], [], [], []
    for a_, d_, l_ in zip(alpha, dcd, dcl):
        s, c = math.sin(a_ / RAD), math.cos(a_ / RAD)
        cn.append(l_ * c + d_ * s)
        ca.append(d_ * c - l_ * s)
        cla.append(tbfunx(alpha, dcl, a_, 0, 0)[1])
        cma.append(tbfunx(alpha, dcm, a_, 0, 0)[1])
    return {'cn': cn, 'ca': ca, 'cla': cla, 'cma': cma}
