"""
Supersonic dynamic pressure at the horizontal tail: DPRESR.

The inviscid part of SDWASH's dynamic-pressure ratio.  The wing's surface
turns the flow at six chord stations (the slopes ``WINGIN(95..100)``):
an expansion at the leading edge takes the flow through Prandtl-Meyer
turns, a compression through an oblique shock (FIG68) and then shocks
until the first turn that expands, after which the rest are Prandtl-Meyer.
The waves from those stations cut the survey plane at the heights ``Z``;
a last turn at the trailing edge returns the flow to the free-stream
direction.  The tail height ``ZJ`` then reads the pressure ratio and Mach
number off the six stations, or, inside the last wave, blends linearly to
the trailing-edge value at the wake.

Reference: datcom-legacy/datcom_2000/dpresr.f
"""

import math
from typing import Dict, Optional, Sequence

from pydatcom.utils.constants import RAD, UNUSED
from pydatcom.utils.legacy_numeric import mach2, tbfunx
from pydatcom.utils.table_lookup import fig68

_KN = (0.0, 0.2, 0.4, 0.6, 0.8, 1.0)


def _qpt(xm: float) -> float:
    """The statement function QPT: q over total pressure."""
    return 0.7 * xm**2 / (1.0 + 0.2 * xm**2)**3.5


def _nu(xm: float) -> float:
    """The Prandtl-Meyer angle as DPRESR writes it, degrees."""
    arg2 = math.sqrt(xm**2 - 1.)
    return (2.4495 * math.atan(.40825 * arg2) - math.atan(arg2)) * RAD


def _shock(m1: float, theta: float, turn: float):
    """The pressure ratio and Mach number behind one oblique shock."""
    arg = math.sin(theta / RAD)**2
    ang = abs((theta + turn) / RAD)
    qq = (m1**2 * arg + 5.0) / (6.0 * m1**2 * math.sin(ang)**2)
    m = ((m1**2 * arg + 5.) / (7. * m1**2 * arg - 1.))**.5 / math.sin(ang)
    return qq, m


def calculate_dpresr(zj: float, zwake: float, alpha_deg: float,
                     dwangl: float, mach: float, cr: float, a12: float,
                     rl2: float, slope: Sequence[float],
                     max_detached: int = 100) -> Dict[str, object]:
    """Translate DPRESR: the dynamic-pressure ratio and Mach number at the
    horizontal tail.

    Args:
        zj: ``ZJ``, the tail height, positive above the wing.
        zwake: ``ZWAKE``, the wake height.
        alpha_deg: ``ALPHJ``.
        dwangl: ``DWANGL``, the downwash angle in radians.
        mach: ``DWA(1)``.
        cr: ``WINGIN(6)``; a12, rl2: ``A(12)``, ``A(24)``.
        slope: ``WINGIN(95..100)``, the surface slope at the six chord
            stations 0, 0.2, ... 1.0, degrees.
        max_detached: A guard on the detached-shock retry, which the
            source repeats without limit.

    Returns:
        ``qqinfy`` (``None`` where the source returns without setting it:
        a compression whose shock leaves the flow subsonic), ``mj``, the
        ``DWA`` words set (``dle`` 231, ``deltaz`` 232, ``xsur`` 233,
        ``theta1`` 234 for a leading-edge compression, ``delte`` 235 and
        ``thete`` 236 for a trailing-edge compression), and the station
        arrays ``z``, ``qq``, ``m``.

    Notes:
        Kept as executed: ``KNUINF``, the free-stream Prandtl-Meyer angle
        of a leading-edge expansion, is implicitly INTEGER, so the angle is
        truncated to whole degrees before the turn is added; a descending ``Z`` (the expansion above the
        wing) is handed to TBFUNX, which then never interpolates and
        returns an end value; the trailing-edge turn is ``-alpha+slope``
        for a compression and ``alpha+slope`` for an expansion whichever
        side the tail is on.
    """
    zupper = zj > 0.0
    alpha = alpha_deg / RAD
    xsur = (a12 - (cr + rl2) * math.tan(alpha)) * math.sin(alpha) + \
        (rl2 + cr) / math.cos(alpha)
    deltaz = cr * math.sin(alpha) + (xsur - cr * math.cos(alpha)) * \
        math.tan(dwangl)
    dle = (alpha * RAD if zupper else -alpha * RAD) - slope[0]
    out: Dict[str, object] = {'dle': dle, 'deltaz': deltaz, 'xsur': xsur}
    x = [k * cr for k in _KN]
    z = []
    m = [0.0] * 7
    qq = [0.0] * 7
    qqpt = [0.0] * 7
    knu = [0.0] * 7
    dknu = [0.0] * 7

    def expand_from(start_knot):
        for knot in range(start_knot, 6):
            dknu[knot] = slope[knot - 1] - slope[knot]
            knu[knot] = knu[knot - 1] + dknu[knot]
            m[knot] = mach2(knu[knot])[0]
            qqpt[knot] = _qpt(m[knot])
            qq[knot] = qqpt[knot] * qq[knot - 1] / qqpt[knot - 1]

    if dle >= 0.0:
        u = math.atan(1. / math.sqrt(mach**2 - 1.))
        ang = u + alpha
        arg1 = xsur * math.sin(alpha) / math.tan(ang)
        arg2 = xsur * math.cos(alpha)
        arg3 = math.sin(alpha) + math.cos(alpha) / math.tan(ang)
        for xn in x:
            zn = (-xn + arg2 - arg1) / arg3
            z.append(zn + deltaz if zupper else -abs(zn) - deltaz)
        dknu[0] = dle
        # KNUINF is implicitly INTEGER: the free-stream angle is truncated.
        knu[0] = int(_nu(mach)) + dknu[0]
        m[0] = mach2(knu[0])[0]
        qqpt[0] = _qpt(m[0])
        qq[0] = qqpt[0] / _qpt(mach)
        expand_from(1)
    else:
        delta = abs(dle)
        for _ in range(max_detached):
            theta1, ier = fig68(mach, delta)
            if ier != 2:
                break
            delta = theta1
        else:
            raise RuntimeError("DPRESR: the shock stays detached")
        out['theta1'] = theta1
        ang = alpha - theta1 / RAD
        arg1 = xsur * math.cos(alpha) - xsur * math.sin(alpha) / math.tan(ang)
        arg2 = math.sin(alpha) + math.cos(alpha) / math.tan(ang)
        for xn in x:
            zn = (-xn + arg1) / arg2
            z.append(abs(zn) + deltaz if zupper else zn + deltaz)
        dknu[0] = dle
        qq[0], m[0] = _shock(mach, theta1, dle)
        if m[0] < 1.0:
            out.update({'qqinfy': None, 'mj': m[0], 'z': z})
            return out
        knu[0] = _nu(m[0])
        for shock_knot in range(1, 6):
            dknu[shock_knot] = slope[shock_knot - 1] - slope[shock_knot]
            theta = fig68(m[shock_knot - 1], dknu[shock_knot])[0]
            if dknu[shock_knot] > 0.0:
                qqpt[shock_knot - 1] = _qpt(m[shock_knot - 1])
                knu[shock_knot - 1] = _nu(m[shock_knot - 1])
                expand_from(shock_knot)
                break
            ratio, m[shock_knot] = _shock(
                m[shock_knot - 1], theta, dknu[shock_knot])
            qq[shock_knot] = ratio * qq[shock_knot - 1]

    if m[5] <= mach:
        dknu[6] = alpha * RAD + slope[5]
        qqpt[5] = _qpt(m[5])
        knu[5] = _nu(m[5])
        knu[6] = knu[5] + dknu[6]
        m[6] = mach2(knu[6])[0]
        qqpt[6] = _qpt(m[6])
        qq[6] = qqpt[6] * qq[5] / qqpt[5]
    else:
        dknu[6] = -alpha * RAD + slope[5]
        delte = abs(dknu[6])
        thete, ier = fig68(m[5], delte)
        if ier == 2:
            dknu[6] = -thete
            delte = thete
            thete, ier = fig68(m[5], delte)
        out.update({'delte': delte, 'thete': thete})
        ratio, m[6] = _shock(m[5], thete, dknu[6])
        qq[6] = ratio * qq[5]

    if zj / z[5] <= 1.0:
        qqinfy = 1.0 + (qq[6] - 1.0) * (zj - zwake) / (z[5] - zwake)
        mj = mach
    else:
        qqinfy = tbfunx(z, qq[:6], zj, 0, 0, ordered=False)[0]
        mj = tbfunx(z, m[:6], zj, 0, 0, ordered=False)[0]
    out.update({'qqinfy': qqinfy, 'mj': mj, 'z': z, 'qq': qq, 'm': m})
    return out
