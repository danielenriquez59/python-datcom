"""
Body area shadowed by the horizontal tail's Mach lines: PTINT2 and BDAREA.

BDAREA rotates the body's side profile (the stations and radii of
``/BODYI/``, closed across the base) into the horizontal tail's chord
plane, finds where the Mach lines from the tail's leading and trailing
edges cross the upper and lower profiles (PTINT2), and sums the triangles
between them (AREA2) to give the body area in the tail's zone of influence
and its centroid.

Reference: datcom-legacy/datcom_2000/ptint2.f, bdarea.f
"""

import math
from typing import Dict, Mapping, Sequence

from pydatcom.utils.constants import RAD
from pydatcom.utils.math_utils import area2


def ptint2(xb: Sequence[float], rb: Sequence[float], a: Sequence[float],
           amu: float, xi: list, yi: list) -> Dict[str, object]:
    """Translate PTINT2: the Mach line's crossings of the body profile.

    Arrays use the source's 1-based layout: ``x[1..nx+1]`` is the upper
    profile (closed down the base), ``x[22..nx+22]`` the lower.

    Args:
        xb, rb: The body stations and radii.
        a: The tail edge's x, z and incidence (degrees).
        amu: The Mach angle, radians.
        xi, yi: 43-word crossing arrays, updated in place.

    Returns:
        ``x``, ``y`` (43-word rotated profiles), ``indxui``, ``indxli``
        (the upper and lower segments crossed), ``inxuie``, ``inxlie``
        (nonzero where the line passes behind the base) and ``abort``.
    """
    nx = len(xb)
    tan_mach_angle = math.sin(amu) / math.cos(amu)
    neg_tan_mach = -tan_mach_angle
    dx = -0.0001 * xb[-1]
    xp = list(xb) + [xb[-1]]
    yp = list(rb) + [-rb[-1]]
    cos_incidence, sin_incidence = (math.cos(-a[2] / RAD),
                                    math.sin(-a[2] / RAD))
    x = [0.0] * 43
    y = [0.0] * 43
    for station in range(1, nx + 2):
        px, py = xp[station - 1], yp[station - 1]
        x[station] = (px - a[0]) * cos_incidence + (py - a[1]) * sin_incidence
        y[station] = (py - a[1]) * cos_incidence - (px - a[0]) * sin_incidence
        x[station + 21] = ((px - a[0]) * cos_incidence +
                             (-py - a[1]) * sin_incidence)
        y[station + 21] = ((-py - a[1]) * cos_incidence -
                           (px - a[0]) * sin_incidence)
    nx22m1 = nx + 21
    if dx < x[nx] < 0.0:
        x[nx] = 0.0
    if dx < x[nx22m1] < 0.0:
        x[nx22m1] = 0.0
    indexu = next((j for j in range(1, nx + 1) if 0.0 <= x[j]), nx + 1) - 1
    indexl = next((j for j in range(1, nx + 1) if 0.0 <= x[j + 21]),
                  nx + 1) + 20
    intersection = {'x': x, 'y': y, 'inxuie': 0, 'inxlie': 0, 'abort': True}

    def run(first, last, slope, sign, below, tag):
        for seg in range(first, last + 1):
            xdif, ydif = x[seg + 1] - x[seg], y[seg + 1] - y[seg]
            if xdif == 0.0 and ydif == 0.0:
                continue
            xi[seg] = (y[seg] * xdif - x[seg] * ydif) / (xdif * slope - ydif)
            yi[seg] = sign * xi[seg] * tan_mach_angle
            if seg == first:
                if x[last] < 0.0:
                    return None
                yin = (ydif / xdif) * (-x[seg]) + y[seg]
                if (yin < 0.0) if below else (yin > 0.0):
                    return None
            if x[seg] <= xi[seg] <= x[seg + 1]:
                return seg
            if seg == last - 1 and xi[seg] > x[seg + 1]:
                intersection[tag] = last - 1
        return last

    ku = run(indexu, nx, tan_mach_angle, 1.0, True, 'inxuie')
    if ku is None:
        return intersection
    kl = run(indexl, nx22m1, neg_tan_mach, -1.0, False, 'inxlie')
    if kl is None:
        return intersection
    intersection.update({'indxui': ku, 'indxli': kl, 'abort': False})
    return intersection


def calculate_bdarea(xb: Sequence[float], rb: Sequence[float], mach: float,
                     syna: Mapping[int, float], htin: Mapping[int, float],
                     aht: Mapping[int, float]) -> Dict[str, object]:
    """Translate BDAREA: the body area in the horizontal tail's Mach zone.

    Args:
        xb, rb: The ``/BODYI/`` stations and radii.
        mach: ``FLC(I+2)``.
        syna: ``/SYNTSS/`` words 1, 6, 7, 8 (``XCG``, ``XH``, ``ZH``,
            ``ALIH``).
        htin: ``HTIN`` 3, 4; aht: ``AHT`` 10, 62.

    Returns:
        ``abort`` (a Mach line misses the body: the routine returns without
        setting anything), else ``sb`` (``HTIN(114+I)``, the area between
        the lines, less what lies behind the base), ``s`` (``HTIN(134+I)``,
        including it) and ``xbar`` (``HTIN(94+I)``, its centroid aft of
        the centre of gravity).
    """
    syna_local = {int(k): float(v) for k, v in syna.items()}
    htin_local = {int(k): float(v) for k, v in htin.items()}
    aht_local = {int(k): float(v) for k, v in aht.items()}
    nx = len(xb)
    amuu = math.atan(1. / math.sqrt(mach**2 - 1.))
    rad8 = syna_local[8] / RAD
    a1 = (syna_local[6] + (htin_local[4] - htin_local[3]) * aht_local[62] *
          math.cos(rad8))
    a2 = (syna_local[7] - (a1 - syna_local[6]) * math.sin(rad8) /
          math.cos(rad8))
    xli, yli = [0.0] * 43, [0.0] * 43
    xti, yti = [0.0] * 43, [0.0] * 43
    le = ptint2(xb, rb, [a1, a2, syna_local[8]], amuu, xli, yli)
    if le['abort']:
        return {'abort': True}
    a1 += aht_local[10] * math.cos(rad8)
    a2 -= aht_local[10] * math.sin(rad8)
    te = ptint2(xb, rb, [a1, a2, syna_local[8]], amuu, xti, yti)
    if te['abort']:
        return {'abort': True}
    xl, yl = le['x'], le['y']
    c10 = aht_local[10]
    parts = []
    for profile_pass in (1, 2):
        if profile_pass == 1:
            ili, ilie = le['indxui'], le['inxuie']
            iti, itie = te['indxui'], te['inxuie']
            last = nx
        else:
            ili, ilie = le['indxli'], le['inxlie']
            iti, itie = te['indxli'], te['inxlie']
            last = nx + 21
        area = areab = saxb = sayb = 0.0

        def tri(x, y, inum):
            return area2(x, y, inum)

        if itie == 0:
            if ili == iti:
                ar, ax, ay = tri([0., xli[ili], xti[iti] + c10, c10],
                                 [0., yli[ili], yti[iti], 0.], 2)
                area = areab = ar
                saxb, sayb = ax, ay
            else:
                ar, ax, ay = tri([0., xli[ili], xl[ili + 1], c10],
                                 [0., yli[ili], yl[ili + 1], 0.], 2)
                area, saxb, sayb = ar, ax, ay
                ii = 1
                while True:
                    ii += 1
                    if xl[ili + ii] > xti[iti] + c10:
                        break
                    ar, ax, ay = tri([xl[ili + ii - 1], xl[ili + ii], c10],
                                     [yl[ili + ii - 1], yl[ili + ii], 0.], 3)
                    area += ar
                    saxb += ax
                    sayb += ay
                ar, ax, ay = tri([xl[ili + ii - 1], xti[iti] + c10, c10],
                                 [yl[ili + ii - 1], yti[iti], 0.], 3)
                area += ar
                areab = area
                saxb += ax
                sayb += ay
        elif ilie == 0:
            ar, ax, ay = tri([0., xli[ili], xl[ili + 1], c10],
                             [0., yli[ili], yl[ili + 1], 0.], 2)
            area, saxb, sayb = ar, ax, ay
            ii = 1
            while ili + ii != last:
                ar, ax, ay = tri([xl[ili + ii], xl[ili + ii + 1], c10],
                                 [yl[ili + ii], yl[ili + ii + 1], 0.], 3)
                area += ar
                saxb += ax
                sayb += ay
                ii += 1
            x3 = [xl[ili + ii], xti[itie] + c10, c10]
            y3 = [yl[ili + ii], yti[itie], 0.]
            ar, ax, ay = tri(x3, y3, 3)
            area += ar
            saxb += ax
            sayb += ay
            extare, ax, ay = tri([x3[0], x3[1], xti[iti] + c10],
                                 [y3[0], y3[1], yti[iti]], 3)
            areab = area - extare
            saxb -= ax
            sayb -= ay
        else:
            ar, _, _ = tri([0., xli[ilie], xti[itie] + c10, c10],
                           [0., yli[ilie], yti[itie], 0.], 2)
            area += ar
            ar, ax, ay = tri([0., xli[ili], xti[iti] + c10, c10],
                             [0., yli[ili], yti[iti], 0.], 2)
            areab += ar
            saxb += ax
            sayb += ay
        parts.append((area, areab, saxb, sayb))
    (au, aub, sxu, syu), (al, alb, sxl, syl) = parts
    sb = aub + alb
    xcb = (sxu + sxl) / sb
    ycb = (syu + syl) / sb
    xbarb = (xcb * math.cos(-rad8) - ycb * math.sin(-rad8) + syna_local[6] +
             (htin_local[4] - htin_local[3]) * aht_local[62] * math.cos(rad8))
    return {'abort': False, 'sb': sb, 's': au + al,
            'xbar': xbarb - syna_local[1],
            'method': 'legacy_bdarea'}
