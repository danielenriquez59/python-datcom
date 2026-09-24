"""
Vertical-panel area shadowed by Mach lines: PTINT1 and VTAREA.

At supersonic speed only the part of a vertical panel outside the Mach
cones from the wing's (or horizontal tail's) leading and trailing edges
feels that surface's sidewash.  VTAREA lays out the panel's inboard (and,
if not straight tapered, outboard) quadrilateral, rotates it into each
surface's chord plane (PTINT1), finds where the Mach line from each edge
crosses the panel's sides, and takes the areas cut off (AREA1) to give the
panel area shadowed by the wing-body (``VTIN(94+I)``) and by the tail
(``VTIN(134+I)``), and what remains (``VTIN(114+I)``).

Reference: datcom-legacy/datcom_2000/ptint1.f, vtarea.f
"""

import math
from typing import Dict, List, Mapping, Sequence

from pydatcom.utils.constants import PI, RAD
from pydatcom.utils.math_utils import area1


def ptint1(xp: Sequence[float], yp: Sequence[float], a: Sequence[float],
           amu: float, j: int, ncon: int, vertup: bool,
           avt: Mapping[int, float], xi: List[float], yi: List[float]
           ) -> Dict[str, object]:
    """Translate PTINT1: the Mach line's crossings of a panel's sides.

    Args:
        xp, yp: The panel's five corners (closed), in body axes.
        a: ``A(3)``: the surface edge's x, z and incidence (degrees).
        amu: The Mach angle, radians.
        j: 1 for the leading edge, 2 for the trailing edge.
        ncon: 0 for the inboard panel, 1 for the outboard.
        vertup: The panel is above the body.
        avt: The *vertical tail's* ``AVT`` words 59, 77, 83, 101 (the
            sweeps PTINT1 reads from ``/VTDATA/``).
        xi, yi: The crossings, updated in place; a side of zero length
            keeps what the list held.

    Returns:
        ``x``, ``y`` (the rotated corners), ``nsum`` (the sum of the side
        numbers crossed), ``effect`` and ``k``, the loop index at exit
        (5 when the loop runs out).
    """
    cos_incidence, sin_incidence = (math.cos(-a[2] / RAD),
                                    math.sin(-a[2] / RAD))
    x = [(xp[n] - a[0]) * cos_incidence + (yp[n] - a[1]) * sin_incidence
         for n in range(5)]
    y = [(yp[n] - a[1]) * cos_incidence - (xp[n] - a[0]) * sin_incidence
         for n in range(5)]
    half_pi = PI / 2.0
    inorot = 24 if ncon == 1 else 0
    le = avt[59 + inorot] + a[2] / RAD
    te = avt[77 + inorot] + a[2] / RAD
    le_up = avt[59 + inorot] - a[2] / RAD
    te_up = avt[77 + inorot] - a[2] / RAD
    icon = nsum = index = 0
    k = 1
    exit_now = False
    while k <= 4:
        i = k - 1
        xdif, ydif = x[k] - x[i], y[k] - y[i]
        crossed = False
        if not (xdif == 0.0 and ydif == 0.0):
            tan_mach_angle = math.sin(amu) / math.cos(amu)
            xi[i] = (y[i] * xdif - x[i] * ydif) / (xdif * tan_mach_angle - ydif)
            yi[i] = xi[i] * tan_mach_angle
            if k == 2:
                crossed = x[1] <= xi[1] <= x[2]
            elif k == 4:
                crossed = x[4] <= xi[3] <= x[3]
            elif not vertup:
                if k == 1:
                    if yi[0] <= y[0]:
                        if yi[0] >= y[1]:
                            crossed = True
                        elif not half_pi + amu < le:
                            index += 1 if j != 2 else 0
                            exit_now = True
                    elif half_pi + amu < le:
                        index += 1 if j != 2 else 0
                        exit_now = True
                else:
                    if y[2] <= yi[2] <= y[3]:
                        crossed = True
                    elif yi[2] < y[2] and half_pi + amu > te:
                        pass
                    else:
                        index += 1 if j == 2 else 0
                        exit_now = True
            else:
                if k == 1:
                    if yi[0] >= y[0]:
                        if yi[0] <= y[1]:
                            crossed = True
                        elif not half_pi - amu < le_up:
                            index += 1 if j != 2 else 0
                            exit_now = True
                    elif half_pi - amu < le_up:
                        index += 1 if j != 2 else 0
                        exit_now = True
                else:
                    if y[3] <= yi[2] <= y[2]:
                        crossed = True
                    elif yi[2] > y[2] and half_pi - amu > te_up:
                        pass
                    else:
                        index += 1 if j == 2 else 0
                        exit_now = True
        if exit_now:
            break
        if crossed:
            icon += 1
            nsum += k
            index += 1
            if icon == 2:
                break
        k += 1
    return {'x': x, 'y': y, 'nsum': nsum, 'effect': index > 0, 'k': k}


def _shadow(nsum, k, x, y, xi, yi):
    """Labels 1080-1170: the corners cut off, AREA1's area and FLIP."""
    flip = False
    if nsum in (1, 2, 3) or nsum > 7:
        xx = [xi[0], x[1], xi[1], 0.0]
        yy = [yi[0], y[1], yi[1], 0.0]
    elif nsum == 4:
        xx = [xi[0], x[1], x[2], xi[2]]
        yy = [yi[0], y[1], y[2], yi[2]]
        flip = abs(xi[0]) > abs(xi[2])
    elif nsum == 5 and k != 4:
        xx = [xi[1], x[2], xi[2], 0.0]
        yy = [yi[1], y[2], yi[2], 0.0]
        flip = True
    elif nsum == 5:
        xx = [x[0], xi[0], xi[3], 0.0]
        yy = [y[0], yi[0], yi[3], 0.0]
    elif nsum == 6:
        xx = [x[0], x[1], xi[1], xi[3]]
        yy = [y[0], y[1], yi[1], yi[3]]
    else:
        xx = [xi[2], x[3], xi[3], 0.0]
        yy = [yi[2], y[3], yi[3], 0.0]
        flip = True
    return area1(xx, yy, nsum), flip


def calculate_vtarea(vtin: Mapping[int, float], avt: Mapping[int, float],
                     vertup: bool, xv: float, zv: float, mach: float,
                     wing: Mapping[str, object], tail: Mapping[str, object],
                     syna: Mapping[int, float], htpl: bool,
                     vt_common: Mapping[int, float], mach_index: int,
                     stale_vtin: Mapping[int, float]
                     ) -> Dict[str, object]:
    """Translate VTAREA: the vertical panel's area in the Mach shadows.

    Args:
        vtin: The panel's ``VTIN`` (or ``VFIN``) words 1-5 and 15 (15 is
            1.0 for straight tapered).
        avt: The panel's ``A`` words 1, 2, 3, 10, 21, 62, 86.
        vertup: The panel is above the body; xv, zv: its position.
        mach: ``FLC(I+2)``.
        wing: ``span``, ``spans`` (``WINGIN(4)``, ``(3)``), ``a62``,
            ``a10`` (``A(62)``, ``A(10)``).
        tail: ``span``, ``spans`` (``HTIN(4)``, ``(3)``), ``a10``,
            ``a16``, ``a30``, ``a62``.
        syna: ``/SYNTSS/`` words 2, 3, 4, 6, 7, 8 (``XW``, ``ZW``,
            ``ALIW``, ``XH``, ``ZH``, ``ALIH``).
        htpl: The horizontal tail is present.
        vt_common: The vertical tail's ``AVT`` words 59, 77, 83, 101,
            which PTINT1 reads from ``/VTDATA/`` whichever panel this is.
        mach_index: ``I``.
        stale_vtin: The panel's ``VTIN(134+I)`` before the call (read when
            there is no tail).

    Returns:
        ``svwb`` (``VTIN(94+I)``), ``svhb`` (``VTIN(134+I)``), ``svb``
        (``VTIN(114+I)``), and ``syna6``, ``syna7`` as the routine leaves
        ``XH`` and ``ZH``.

    Notes:
        Kept as executed: PTINT1 reads the vertical tail's sweeps for the
        ventral fin too; the outboard panel's heights are not mirrored for
        a panel below the body; and turning the tail to the wing's
        incidence rewrites ``XH`` and ``ZH`` in ``/SYNTSS/`` for the rest
        of the run.
    """
    v = {int(k): float(val) for k, val in vtin.items()}
    g = {int(k): float(val) for k, val in avt.items()}
    s = {int(k): float(val) for k, val in syna.items()}
    xi, yi = [0.0] * 4, [0.0] * 4
    out = {'svhb': float(stale_vtin[134 + mach_index])}
    htail = angchg = False
    rad_w = s[4] / RAD
    while True:
        ncon = 0
        sv = [0.0, 0.0]
        yp1 = v[4] - v[3] + zv if vertup else v[4] - v[3] - zv
        xp1 = xv + (v[4] - v[3]) * g[62]
        xp2 = xv + g[21] * g[62]
        xp = [xp1, xp2, xp2 + v[5], xp1 + g[10], xp1]
        yp = [yp1, g[21], g[21], yp1, yp1]
        if not vertup:
            yp = [-y_ for y_ in yp]
        while True:
            amuu = math.atan(1. / math.sqrt(mach**2 - 1.))
            area = [0.0, 0.0]
            done = False
            for j in (1, 2):
                if not htail:
                    a1 = s[2] + (float(wing['span']) - float(wing['spans'])) \
                        * float(wing['a62']) * math.cos(rad_w)
                    a2 = s[3] - (a1 - s[2]) * math.sin(rad_w) / \
                        math.cos(rad_w)
                    if j == 2:
                        a1 += float(wing['a10']) * math.cos(rad_w)
                        a2 -= float(wing['a10']) * math.sin(rad_w)
                else:
                    if not (s[4] == s[8] or angchg):
                        hacle = (float(tail['span']) - float(tail['spans'])) \
                            * float(tail['a62']) + float(tail['a30']) - \
                            float(tail['a16']) / 4.
                        xhac = s[6] + hacle * math.cos(s[8] / RAD)
                        zhac = s[7] - hacle * math.sin(s[8] / RAD)
                        s[6] = xhac - hacle * math.cos(rad_w)
                        s[7] = zhac + hacle * math.sin(rad_w)
                        angchg = True
                    a1 = s[6] + (float(tail['span']) - float(tail['spans'])) \
                        * float(tail['a62']) * math.cos(rad_w)
                    a2 = s[7] - (a1 - s[6]) * math.sin(rad_w) / \
                        math.cos(rad_w)
                    if j == 2:
                        a1 += float(tail['a10']) * math.cos(rad_w)
                        a2 -= float(tail['a10']) * math.sin(rad_w)
                p = ptint1(xp, yp, [a1, a2, s[4]], amuu, j, ncon, vertup,
                           vt_common, xi, yi)
                if not p['effect']:
                    if j == 1:
                        done = True
                        break
                    area[1] = g[1 + ncon] if p['k'] == 1 else 0.
                    sv[ncon] = g[1 + ncon] - (area[0] + area[1])
                    continue
                flip = False
                if p['nsum'] == 0:
                    area[j - 1] = 0.
                    flip = j != 1
                else:
                    area[j - 1], flip = _shadow(p['nsum'], p['k'], p['x'],
                                                p['y'], xi, yi)
                if j == 1:
                    if flip:
                        area[0] = g[1 + ncon] - area[0]
                    continue
                if not flip:
                    area[1] = g[1 + ncon] - area[1]
                sv[ncon] = g[1 + ncon] - (area[0] + area[1])
            if done or ncon == 1 or v[15] == 1.0:
                break
            ncon = 1
            xp1, yp1 = xp[1], yp[1]
            xp2 = xp[1] + v[2] * g[86]
            xp = [xp1, xp2, xp2 + v[1], xp1 + v[5], xp1]
            yp = [yp1, v[4], v[4], yp1, yp1]
        if htail:
            out['svhb'] = sv[0] + sv[1]
            break
        out['svwb'] = sv[0] + sv[1]
        if not htpl:
            break
        htail = True
    out['svb'] = g[3] - (out['svwb'] + out['svhb'])
    out.update({'syna6': s[6], 'syna7': s[7], 'method': 'legacy_vtarea'})
    return out
