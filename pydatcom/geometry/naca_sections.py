"""
NACA section coordinates as the legacy program forms them: COORD1, COORD4,
COORD5, COORD6, CORD4M, CORD5M and XYCORD.

Each routine fills the upper and lower surface coordinates, the mean line
and the half-thickness at the chord stations ``X`` from the digits of the
designation (``/IBW/``'s ``I, J, K, II, JJ, KK, III, JJJ``): the 1-series,
four- and five-digit (standard and modified) and 6-series sections, and
XYCORD for a section given by its own ordinates.  The modified thickness
forms solve two small linear systems with SLEQ, as the source does.

``geometry/airfoil.py`` has an independent NACA generator; these are the
source's own forms, with its quirks, for the translated pipeline.

Reference: datcom-legacy/datcom_2000/coord1.f, coord4.f, coord5.f,
coord6.f, cord4m.f, cord5m.f, xycord.f
"""

import math
from typing import Dict, List, Mapping, Optional, Sequence

import numpy as np

from pydatcom.utils.constants import PI, RAD
from pydatcom.utils.legacy_numeric import arccos, sleq, tbfunx


def _pow(x: float, y: float) -> float:
    """``X**Y`` for a real exponent: NaN for a negative base, as compiled."""
    return float(np.float64(x) ** np.float64(y)) if x < 0 else x ** y


def _digits(d: Mapping[str, int]):
    return [float(d.get(k, 0)) for k in ('i', 'j', 'k', 'ii', 'jj', 'kk',
                                         'iii', 'jjj')]


def _four_digit_thickness(t: float, x: float) -> float:
    return 5. * t * (.2969 * math.sqrt(x) - .126 * x - .3516 * x**2 +
                     .2843 * x**3 - .1015 * x**4)


def _modified_thickness(zt: float, t: float, d0: float, d1: float,
                        a0: float):
    """The modified thickness form's coefficients ``D2, D3, A1, A2, A3``,
    from two SLEQ solutions."""
    xn, _ = sleq([[-2. * (1. - zt), -3. * ((1. - zt)**2)],
                  [(1. - zt)**2, (1. - zt)**3]],
                 [d1, t / 2. - d0 - d1 * (1. - zt)])
    d2, d3 = float(xn[0]), float(xn[1])
    xm, _ = sleq([[0.0, 2.0, 6. * zt], [1.0, 2. * zt, 3. * zt**2],
                  [zt, zt**2, zt**3]],
                 [2. * d2 + 6. * d3 * (1. - zt) + a0 / (4. * zt**1.5),
                  -a0 / (2. * zt**.5), -a0 * zt**.5 + t / 2.])
    return d2, d3, float(xm[0]), float(xm[1]), float(xm[2])


def _yt_modified(x, zt, t, d0, d1, d2, d3, a0, a1, a2, a3, yt):
    if x == zt:
        yt = t / 2.
    if x < zt:
        yt = a0 * x**.5 + a1 * x + a2 * x**2 + a3 * x**3
    if x > zt:
        yt = d0 + d1 * (1. - x) + d2 * (1. - x)**2 + d3 * (1. - x)**3
    return yt


def _five_digit_mean(zp: float):
    """ZM from the cubic of the five-digit mean line."""
    a = 6. * zp - 3.
    b = -2. + 6. * zp - 3. * zp**2
    g = b * b / 4. + a * a * a / 27.
    if g < 0.:
        phi = arccos((-b / 2.) / _pow(-a**3 / 27., .5))
        return 1. + 2. * _pow(-a / 3., .5) * math.cos(phi / 3. + 4.18879)
    d = _pow(-b / 2. + g**.5, 1. / 3.)
    e = _pow(-b / 2. - g**.5, 1. / 3.)
    return d + e + 1.


def _surfaces(coords, m, x, yc, yt, alpha):
    coords['xu'][m] = x - yt * math.sin(alpha)
    coords['yun'][m] = yc + yt * math.cos(alpha)
    coords['xl'][m] = x + yt * math.sin(alpha)
    coords['yln'][m] = yc - yt * math.cos(alpha)
    coords['cam'][m] = 0.0 if yc < 1.e-05 else yc
    coords['thn'][m] = yt


def _start(x, prev):
    n = len(x)
    coords = {k: list(prev[k]) if prev and k in prev else [0.0] * n
              for k in ('xu', 'xl', 'yun', 'yln', 'thn', 'cam')}
    return coords


def _close(coords, n, lower_first=False):
    """The shared closing assignments; COORD1 writes ``XL(1)`` twice where
    the others write ``XU(1)``."""
    coords['thn'][0] = coords['thn'][n - 1] = 0.0
    coords['cam'][0] = coords['cam'][n - 1] = 0.0
    coords['xu'][n - 1] = coords['xl'][n - 1] = 1.
    coords['yun'][n - 1] = coords['yln'][n - 1] = 0.0
    if lower_first:
        coords['xl'][0] = 0.0
    else:
        coords['xu'][0] = 0.0
    coords['yln'][0] = 0.0


def coord4(digits: Mapping[str, int], x: Sequence[float],
           prev: Optional[Mapping[str, list]] = None) -> Dict[str, object]:
    """Translate COORD4: NACA four-digit coordinates."""
    ai, aj, ak, aii, ajj, akk, _, _ = _digits(digits)
    zm, zp = ai * .01, aj * .1
    t = ak * .1 + aii * .01 + ajj * .001 + akk * .0001
    coords = _start(x, prev)
    yc = alpha = 0.0
    for m, xm in enumerate(x):
        yt = _four_digit_thickness(t, xm)
        if xm == zp:
            yc, alpha = zm, 0.0
        if xm < zp:
            yc = (2. * zp * xm - xm**2) * zm / zp**2
            alpha = math.atan((2. * zm / (zp**2)) * (zp - xm))
        if xm > zp:
            yc = (zm / ((1. - zp)**2)) * (1. - 2. * zp + 2. * zp * xm - xm**2)
            alpha = math.atan((2. * zm / ((1. - zp)**2)) * (zp - xm))
        _surfaces(coords, m, xm, yc, yt, alpha)
    _close(coords, len(x))
    coords.update({'rho': 1.1019 * t**2, 't': t, 'zm': zm, 'zp': zp})
    return coords


def coord5(digits: Mapping[str, int], x: Sequence[float],
           prev: Optional[Mapping[str, list]] = None,
           stale: Optional[Mapping[str, float]] = None) -> Dict[str, object]:
    """Translate COORD5: NACA five-digit coordinates, standard (``K=0``)
    or reflexed.

    ``stale`` holds ``yc`` and ``alpha`` as the previous call left them:
    a station exactly at ``ZM`` (and not at ``ZP``) sets neither."""
    ai, aj, ak, aii, ajj, akk, aiii, _ = _digits(digits)
    t = aii * .1 + ajj * .01 + akk * .001 + aiii * .0001
    zp = aj * .1 / 2.
    zm = _five_digit_mean(zp)
    xk = (6. * ai * .01) / (zp**3 - 3. * zm * zp**2 + zm**2 * (3. - zm) * zp)
    stale_state = dict(stale or {})
    yc, alpha = float(stale_state.get('yc', 0.0)), float(stale_state.get('alpha', 0.0))
    coords = _start(x, prev)
    for m, xm in enumerate(x):
        yt = _four_digit_thickness(t, xm)
        if ak == 0.:
            if xm < zm:
                yc = (1. / 6.) * xk * (xm**3 - 3. * zm * xm**2 +
                                       zm**2 * (3. - zm) * xm)
                alpha = math.atan((1. / 6.) * xk * (3. * xm**2 - 6. * zm * xm +
                                                    zm**2 * (3. - zm)))
            if xm == zp:
                yc, alpha = ai * .01, 0.0
            if xm > zm:
                yc = (1. / 6.) * xk * zm**3 * (1. - xm)
                alpha = math.atan(-(1. / 6.) * xk * zm**3)
        else:
            rk = (3. * ((zm - zp)**2) - zm**3) / ((1. - zm)**3)
            xk = (6. * ai * .01) / ((zp - zm)**3 - rk * (1. - zm)**3 * zp -
                                    (zm**3) * zp + zm**3)
            if xm < zm:
                yc = (1. / 6.) * xk * ((xm - zm)**3 - rk * xm * (1. - zm)**3 -
                                       zm**3 * xm + zm**3)
                alpha = math.atan((1. / 6.) * xk * (3. * (xm - zm)**2 -
                                                    rk * (1. - zm)**3 -
                                                    zm**3))
            if xm == zp:
                yc, alpha = ai * .01, 0.0
            if xm > zm:
                yc = (1. / 6.) * xk * (rk * (xm - zm)**3 -
                                       rk * xm * (1. - zm)**3 -
                                       xm * zm**3 + zm**3)
                alpha = math.atan((1. / 6.) * xk * (3. * rk * (xm - zm)**2 -
                                                    rk * (1. - zm)**3 -
                                                    zm**3))
        _surfaces(coords, m, xm, yc, yt, alpha)
    _close(coords, len(x))
    coords.update({'rho': 1.1019 * t**2, 't': t, 'zm': zm, 'zp': zp,
                   'stale': {'yc': yc, 'alpha': alpha}})
    return coords


def _modified_d1(zt, t, stale_d1):
    for lo, hi, f in ((.18, .22, 1.), (.28, .32, 1.17), (.38, .42, 1.575),
                      (.48, .52, 2.325), (.58, .62, 3.5)):
        if lo < zt < hi:
            return f * t
    return stale_d1


def cord4m(digits: Mapping[str, int], x: Sequence[float],
           prev: Optional[Mapping[str, list]] = None,
           stale: Optional[Mapping[str, float]] = None) -> Dict[str, object]:
    """Translate CORD4M: NACA four-digit modified coordinates.

    ``stale['d1']`` is the trailing-edge slope the previous call left,
    used when the maximum-thickness station is not one of 0.2 ... 0.6."""
    ai, aj, ak, aii, _, akk, aiii, _ = _digits(digits)
    stale_state = dict(stale or {})
    zm, zp, zt = ai * .01, aj * .1, aiii * .1
    t = ak * .1 + aii * .01
    d1 = _modified_d1(zt, t, float(stale_state.get('d1', 0.0)))
    d0 = .01 * t
    a0 = math.sqrt(2. * 1.1019 * ((t * akk / 6.)**2))
    d2, d3, a1, a2, a3 = _modified_thickness(zt, t, d0, d1, a0)
    coords = _start(x, prev)
    yc = alpha = yt = 0.0
    for m, xm in enumerate(x):
        if xm == zp:
            yc, alpha = zm, 0.0
        if xm < zp:
            yc = (zm / zp**2) * (2. * zp * xm - xm**2)
            alpha = math.atan((2. * zm / zp**2) * (zp - xm))
        if xm > zp:
            yc = (zm / ((1. - zp)**2)) * (1. - 2. * zp + 2. * zp * xm - xm**2)
            alpha = math.atan((2. * zm / ((1. - zp)**2)) * (zp - xm))
        yt = _yt_modified(xm, zt, t, d0, d1, d2, d3, a0, a1, a2, a3, yt)
        _surfaces(coords, m, xm, yc, yt, alpha)
    _close(coords, len(x))
    coords.update({'rho': .5 * a0**2, 't': t, 'zm': zm, 'zp': zp,
                   'stale': {'d1': d1}})
    return coords


def cord5m(digits: Mapping[str, int], x: Sequence[float],
           prev: Optional[Mapping[str, list]] = None,
           stale: Optional[Mapping[str, float]] = None) -> Dict[str, object]:
    """Translate CORD5M: NACA five-digit modified coordinates.

    ``stale`` holds ``d1``, ``yc`` and ``alpha`` as the previous call left
    them (see :func:`cord4m` and :func:`coord5`)."""
    ai, aj, ak, aii, ajj, _, aiii, ajjj = _digits(digits)
    stale_state = dict(stale or {})
    t = aii * .1 + ajj * .01
    zp = aj * .1 / 2.
    zt = ajjj * .1
    zm = _five_digit_mean(zp)
    xk = (6. * ai * .01) / (zp**3 - 3. * zm * zp**2 + zm**2 * (3. - zm) * zp)
    d1 = _modified_d1(zt, t, float(stale_state.get('d1', 0.0)))
    d0 = .01 * t
    a0 = math.sqrt(2. * 1.1019 * ((t * aiii / 6.)**2))
    d2, d3, a1, a2, a3 = _modified_thickness(zt, t, d0, d1, a0)
    coords = _start(x, prev)
    yc, alpha = float(stale_state.get('yc', 0.0)), float(stale_state.get('alpha', 0.0))
    yt = 0.0
    for m, xm in enumerate(x):
        if ak == 0.:
            if xm == zt:
                yt = t / 2.
            if xm == zp:
                yc, alpha = ai * .01, 0.0
            if xm < zm:
                yc = (1. / 6.) * xk * (xm**3 - 3. * zm * xm**2 +
                                       zm**2 * (3. - zm) * xm)
                alpha = math.atan((1. / 6.) * xk * (3. * xm**2 - 6. * zm * xm +
                                                    zm**2 * (3. - zm)))
            if xm < zt:
                yt = a0 * xm**.5 + a1 * xm + a2 * xm**2 + a3 * xm**3
            if xm > zt:
                yt = d0 + d1 * (1. - xm) + d2 * (1. - xm)**2 + \
                    d3 * (1. - xm)**3
            if xm > zm:
                yc = (1. / 6.) * xk * zm**3 * (1. - xm)
                alpha = math.atan(-(1. / 6.) * xk * zm**3)
        else:
            rk = (3. * ((zm - zp)**2) - zm**3) / ((1. - zm)**3)
            xk = (6. * ai * .01) / ((zp - zm)**3 - rk * (1. - zm)**3 * zp -
                                    zm**3 * zp + zm**3)
            if xm < zt:
                yt = a0 * xm**.5 + a1 * xm + a2 * xm**2 + a3 * xm**3
            if xm < zm:
                yc = (1. / 6.) * xk * ((xm - zm)**3 - rk * xm * (1. - zm)**3 -
                                       zm**3 * xm + zm**3)
                alpha = math.atan((1. / 6.) * xk * (3. * (xm - zm)**2 -
                                                    rk * (1. - zm)**3 -
                                                    zm**3))
            if xm == zt:
                yt = t / 2.
            if xm == zp:
                yc, alpha = ai * .01, 0.0
            if xm > zt:
                yt = d0 + d1 * (1. - xm) + d2 * (1. - xm)**2 + \
                    d3 * (1. - xm)**3
            if xm > zm:
                yc = (1. / 6.) * xk * (rk * (xm - zm)**3 -
                                       rk * xm * (1. - zm)**3 -
                                       xm * zm**3 + zm**3)
                alpha = math.atan((1. / 6.) * xk * (3. * rk * (xm - zm)**2 -
                                                    rk * (1. - zm)**3 -
                                                    zm**3))
        _surfaces(coords, m, xm, yc, yt, alpha)
    _close(coords, len(x))
    coords.update({'rho': .5 * a0**2, 't': t, 'zm': zm, 'zp': zp,
                   'stale': {'d1': d1, 'yc': yc, 'alpha': alpha}})
    return coords


def coord1(digits: Mapping[str, int], x: Sequence[float],
           prev: Optional[Mapping[str, list]] = None,
           stale: Optional[Mapping[str, float]] = None) -> Dict[str, object]:
    """Translate COORD1: NACA 1-series coordinates.

    Only the interior stations are computed; the closing assignments write
    ``XL(1)`` twice where the other routines write ``XU(1)``, so ``XU(1)``
    and ``YUN(1)`` keep what ``prev`` held.  ``stale`` holds ``d1`` and
    ``sm``, left unset for a series digit other than 6, 8 or 9."""
    ai, aj, ak, aii, ajj, akk, aiii, ajjj = _digits(digits)
    stale_state = dict(stale or {})
    j = int(aj)
    zt = aj * .1 - .1
    t = ajj * .1 + akk * .01 + aiii * .001 + ajjj * .0001
    if j == 6:
        zt = aj * .1 - .2
    d0 = 0.0
    d1, sm = float(stale_state.get('d1', 0.0)), float(stale_state.get('sm', 0.0))
    if j == 6:
        d1, sm = 2.157 * t, 4.
    if j == 8:
        d1, sm = 3.6833 * t, 3.
    if j == 9:
        d1, sm = 5.5283 * t, 3.
    a0 = math.sqrt(2. * 1.1019 * ((t * sm / 6.)**2))
    d2, d3, a1, a2, a3 = _modified_thickness(zt, t, d0, d1, a0)
    cl = aii * .1
    coords = _start(x, prev)
    n = len(x)
    yt = float(stale_state.get('yt', 0.0))
    for m in range(1, n - 1):
        xm = x[m]
        yc = -(cl / (4. * PI)) * ((1. - xm) * math.log(1. - xm) +
                                  xm * math.log(xm))
        alpha = math.atan((-cl / (4. * PI)) * (math.log(xm) -
                                               math.log(1. - xm)))
        yt = _yt_modified(xm, zt, t, d0, d1, d2, d3, a0, a1, a2, a3, yt)
        _surfaces(coords, m, xm, yc, yt, alpha)
    _close(coords, n, lower_first=True)
    coords.update({'rho': .5 * a0**2, 't': t, 'alphai': 0.0,
                   'alphao': -RAD * cl / (2. * PI), 'aii': cl,
                   'stale': {'d1': d1, 'sm': sm, 'yt': yt}})
    return coords


def coord6(digits: Mapping[str, int], x: Sequence[float],
           prev: Optional[Mapping[str, list]] = None,
           stale: Optional[Mapping[str, float]] = None) -> Dict[str, object]:
    """Translate COORD6: NACA 6-series coordinates, with the ``a=1`` or
    general ``a`` mean line and, for a series with a subscript digit, the
    straight trailing edge from the first station past 0.8 chord.

    ``stale`` holds ``no``, ``sxu``, ``sxl``, ``syu``, ``syl``, ``smu`` and
    ``sml``: the trailing-edge line when the first interior station
    already lies past 0.8 chord."""
    ai, aj, ak, aii, ajj, akk, aiii, ajjj = _digits(digits)
    stale_state = dict(stale or {})
    subscript_digit = int(aii)
    series_digit = int(aj)
    t = akk * .1 + aiii * .01
    if series_digit == 3:
        zt, sm, r0, sub = .35, -.6116, .46, 1.
    elif series_digit == 4:
        zt, sm, r0, sub = .40, -.6888, .523, 1.04
    elif series_digit == 5:
        zt, sm, r0, sub = .4, -.8833, .65, 1.17
    else:
        zt, sm, r0, sub = .45, -1.268, .873, None
    d1 = (sm * (t - .06) + r0) * t
    if sub is not None and subscript_digit > 0:
        d1 = sub * t
    d0 = 0.0
    rle = .01 * (68.682 * t**2 + .0182 * t + .0014)
    a0 = math.sqrt(2. * rle)
    d2, d3, a1, a2, a3 = _modified_thickness(zt, t, d0, d1, a0)
    za = ajjj * .1
    if ajjj < 1.:
        za = 1.
    mean_denominator = mean_g = mean_h = None
    if za != 1.:
        mean_denominator = 1. - za
        mean_g = (-1. / mean_denominator) * (
            (za**2) * (.5 * math.log(za) - .25) + .25)
        mean_h = ((1. / mean_denominator) *
                  ((.5 * mean_denominator**2) * math.log(mean_denominator) -
                   .25 * mean_denominator**2) + mean_g)
    cl = ajj * .1
    coords = _start(x, prev)
    n = len(x)
    no = int(stale_state.get('no', 1))
    sxu, sxl = float(stale_state.get('sxu', 0.)), float(stale_state.get('sxl', 0.))
    syu, syl = float(stale_state.get('syu', 0.)), float(stale_state.get('syl', 0.))
    smu, sml = float(stale_state.get('smu', 0.)), float(stale_state.get('sml', 0.))
    yt = float(stale_state.get('yt', 0.0))
    for m in range(1, n - 1):
        xm = x[m]
        if za != 1.:
            complement_x = 1. - xm
            za_minus_x = za - xm
            if za_minus_x == 0.0:
                za_minus_x = 1.0e-10
            yc = (cl / (2. * PI * (za + 1.))) * (
                (1. / mean_denominator) * (
                    (.5 * za_minus_x**2) * math.log(abs(za_minus_x)) -
                    (.5 * complement_x**2) * math.log(complement_x) +
                    .25 * complement_x**2 - .25 * za_minus_x**2) -
                xm * math.log(xm) + mean_g - xm * mean_h)
            alpha = math.atan((cl / (2. * PI * (1. + za))) * (
                (1. / mean_denominator) * (
                    -za_minus_x * math.log(abs(za_minus_x)) +
                    complement_x * math.log(complement_x)) -
                math.log(xm) - 1. - mean_h))
        else:
            yc = -(cl / (4. * PI)) * ((1. - xm) * math.log(1. - xm) +
                                      xm * math.log(xm))
            alpha = math.atan((-cl / (4. * PI)) * (math.log(xm) -
                                                   math.log(1. - xm)))
        yt = _yt_modified(xm, zt, t, d0, d1, d2, d3, a0, a1, a2, a3, yt)
        coords['xu'][m] = xm - yt * math.sin(alpha)
        coords['yun'][m] = yc + yt * math.cos(alpha)
        coords['xl'][m] = xm + yt * math.sin(alpha)
        coords['yln'][m] = yc - yt * math.cos(alpha)
        if coords['xu'][m] >= .80 and subscript_digit > 0:
            if no == 1:
                sxu, sxl = coords['xu'][m], coords['xl'][m]
                syu, syl = coords['yun'][m], coords['yln'][m]
                smu = -syu / (1. - sxu)
                sml = -syl / (1. - sxl)
                no = 2
            else:
                xu_, xl_ = coords['xu'][m] - sxu, coords['xl'][m] - sxl
                coords['yun'][m] = smu * xu_ + syu
                coords['yln'][m] = sml * xl_ + syl
                coords['xu'][m] = xu_ + sxu
                coords['xl'][m] = xl_ + sxl
        else:
            no = 1
        coords['cam'][m] = 0.0 if yc < 1.e-05 else yc
        coords['thn'][m] = yt
    _close(coords, n)
    coords['xl'][0] = 0.0
    coords['yun'][0] = 0.0
    alphao = -RAD * cl / (2. * PI)
    alphai = 0.0
    if za != 1.:
        alphai = alphao * mean_h / (za + 1.)
        alphao = alphao + alphai
    coords.update({'rho': rle, 't': t, 'alphai': alphai, 'alphao': alphao,
                   'ajj': cl,
                   'stale': {'no': no, 'sxu': sxu, 'sxl': sxl, 'syu': syu,
                             'syl': syl, 'smu': smu, 'sml': sml, 'yt': yt}})
    return coords


def xycord(x: Sequence[float], yu: Sequence[float], yl: Sequence[float],
           ival: int = 0, thn: Optional[Sequence[float]] = None,
           cam: Optional[Sequence[float]] = None) -> Dict[str, object]:
    """Translate XYCORD: the mean line and thickness from given ordinates.

    With ``ival`` 0 the half-thickness and mean line are formed from the
    upper and lower ordinates; otherwise ``thn`` and ``cam`` are taken as
    given.  The surfaces are rebuilt about the mean line's TBFUNX slope.
    (The source's optional print, and its doubling and halving of the
    thickness around it, leave the arrays unchanged.)"""
    n = len(x)
    if ival == 0:
        thn = [0.5 * (a - b) for a, b in zip(yu, yl)]
        cam = [0.5 * (a + b) for a, b in zip(yu, yl)]
    else:
        thn, cam = list(thn), list(cam)
    thn[0] = thn[n - 1] = 0.0
    cam[0] = cam[n - 1] = 0.0
    xu, xl, yun, yln = [0.0] * n, [0.0] * n, [0.0] * n, [0.0] * n
    for station in range(n):
        dydx = tbfunx(x, cam, x[station], 0, 0)[1]
        slope_angle = math.atan(dydx)
        sin_slope = math.sin(slope_angle)
        cos_slope = math.cos(slope_angle)
        xu[station] = x[station] - thn[station] * sin_slope
        xl[station] = x[station] + thn[station] * sin_slope
        yun[station] = cam[station] + thn[station] * cos_slope
        yln[station] = cam[station] - thn[station] * cos_slope
    xu[0] = xl[0] = 0.0
    xu[n - 1] = xl[n - 1] = 1.0
    yun[0] = yun[n - 1] = yln[0] = yln[n - 1] = 0.0
    return {'xu': xu, 'xl': xl, 'yun': yun, 'yln': yln, 'thn': thn,
            'cam': cam, 'method': 'legacy_xycord'}


def cordsp(digits: Mapping[str, int], x: Sequence[float],
           surface_in: Mapping[int, float]) -> Dict[str, object]:
    """Translate CORDSP: supersonic section coordinates from the NACA card.

    ``I`` selects the section (1 double wedge, 2 biconvex, 3 hexagonal;
    any other value falls through to the biconvex, as the computed GO TO
    does), ``J K II`` the maximum-thickness station, ``JJ KK III`` the
    thickness and ``JJJ KKK LLL`` the hexagonal flat.  All 60 stations of
    ``X`` are filled, whatever the station count.

    Args:
        digits: ``i``, ``j``, ``k``, ``ii``, ``jj``, ``kk``, ``iii``,
            ``jjj``, ``kkk``, ``lll``.
        x: The 60 words of ``X``.
        surface_in: The named surface's input words 15, 16, 18, 62, 63,
            70 and 71, whose UNUSED entries CORDSP fills.

    Returns:
        ``xu``, ``xl``, ``yu``, ``yl`` (60 words), ``toc``, ``xt``, ``xf``,
        ``ksharp``, ``rho`` (0), and ``surface_in`` as left.

    Notes:
        Kept as executed: the planform test compares the REAL ``WGIN(15)``
        with an INTEGER Hollerith constant, a numeric comparison that never
        finds them equal, so ``WGIN(63)`` is zeroed when UNUSED even for a
        straight wing.
    """
    card_digits = {k: float(digits.get(k, 0)) for k in (
        'i', 'j', 'k', 'ii', 'jj', 'kk', 'iii', 'jjj', 'kkk', 'lll')}
    xt = (100. * card_digits['j'] + 10. * card_digits['k'] +
          card_digits['ii']) / 1000.
    toc = (100. * card_digits['jj'] + 10. * card_digits['kk'] +
           card_digits['iii']) / 1000.
    xf = (100. * card_digits['jjj'] + 10. * card_digits['kkk'] +
          card_digits['lll']) / 1000.
    kind = int(card_digits['i'])
    yu = []
    if kind == 1:
        ksharp = (1. / xt) / (1. - xt)
        for xi in x:
            yu.append(xi * toc / (2. * xt) if not xi > xt else
                      toc / 2. - (xi - xt) * toc / (2. * (1. - xt)))
    elif kind == 3:
        ksharp = (1. - xf) / (xt * (1. - xt - xf))
        for xi in x:
            if xi >= xt + xf:
                yu.append(toc / 2. - (xi - xt - xf) * toc /
                          (2. * (1. - xt - xf)))
            elif xi >= xt:
                yu.append(toc / 2.)
            else:
                yu.append(xi * toc / (2. * xt))
    else:
        xt = 0.50
        ksharp = 16. / 3.
        rc = (toc**2 + 1.) / (4. * toc)
        yu = [toc / 2. - rc + math.sqrt(rc**2 - (xi - .5)**2) for xi in x]
    unused = 1.0e-30
    surface_words = {int(k): float(v) for k, v in surface_in.items()}
    for word, value in ((16, toc), (18, xt), (70, toc), (71, ksharp),
                        (62, 0.)):
        if surface_words[word] == unused:
            surface_words[word] = value
    if surface_words[63] == unused:
        surface_words[63] = 0.0
    return {'xu': list(x), 'xl': list(x), 'yu': yu, 'yl': [-v for v in yu],
            'toc': toc, 'xt': xt, 'xf': xf, 'ksharp': ksharp, 'rho': 0.,
            'surface_in': surface_words, 'method': 'legacy_cordsp'}


_DIGIT = {'1': 1, '2': 2, '3': 3, '4': 4, '5': 5, '6': 6, '7': 7, '8': 8,
          '9': 9, '0': 0}
_SERIES = {'1': 5, '4': 1, '5': 3, '6': 6, 'S': 7}


def decode(card: str, stale_na: int = 0) -> Dict[str, object]:
    """Translate DECODE: the NACA card's designation and the chord grid.

    Column 8 names the family (1, 4, 5, 6 or S) and the digits start in
    column 10.  For the four- and five-digit and 1-series cards each
    column from 10 takes the next of ten digit slots: a digit fills it,
    a blank or point gives it back, a hyphen consumes it (and, except for
    the 1-series, marks the modified section); the 6-series and
    supersonic cards follow their own rules for ``A``, ``=``, ``,``, ``-``
    and ``.``.

    Args:
        card: The 80-column card image.
        stale_na: ``NA`` as the caller held it, kept when column 8 names
            no family.

    Returns:
        ``na`` (1 four-digit, 2 modified four-digit, 3 five-digit, 4
        modified five-digit, 5 1-series, 6 6-series, 7 supersonic),
        ``digits`` (``IUM(1..10)``: ``I, J, K, II, JJ, KK, III, JJJ, KKK,
        LLL``), ``x`` (the chord stations) and ``l`` (their count).
    """
    padded_card = ' ' + card.ljust(80)[:80]
    na = _SERIES.get(padded_card[8], stale_na)
    digit_slots = [0] * 13
    if na == 6:
        slot_index = 0
        saw_equals = saw_decimal = False
        for column_index in range(10, 81):
            column_char = padded_card[column_index]
            slot_index += 1
            if slot_index > 8:
                break
            if column_char in _DIGIT:
                digit_slots[slot_index] = _DIGIT[column_char]
            if slot_index == 3 and column_char == 'A':
                slot_index += 1
            if slot_index == 4 and column_char == 'A':
                digit_slots[slot_index] = 1
            if slot_index > 4 and column_char == 'A':
                slot_index -= 1
            if column_char == ' ':
                slot_index -= 1
            if column_char == ',':
                slot_index -= 1
            if column_char == '-' and slot_index == 3:
                slot_index += 1
            if column_char == '=':
                saw_equals = True
            if saw_equals and not saw_decimal:
                slot_index -= 1
            if column_char == ' ' and saw_equals and not saw_decimal:
                slot_index += 1
            if column_char == '.':
                saw_decimal = True
    elif na == 7:
        slot_index, digits_in_group = 0, 2
        for column_index in range(10, 81):
            column_char = padded_card[column_index]
            if slot_index == 10:
                continue
            action = None
            if column_char == ' ' and column_index == 80:
                action = 1090
            elif column_char == ' ':
                continue
            else:
                if column_char in _DIGIT:
                    digits_in_group += 1
                    if digits_in_group > 3:
                        continue
                    slot_index += 1
                    digit_slots[slot_index] = _DIGIT[column_char]
                    action = 1090 if column_index == 80 else None
                elif column_char == '-' and slot_index > 0:
                    action = 1090
                elif column_char == '.' and slot_index > 0:
                    action = 1100
                elif column_index == 80:
                    action = 1090
            if action == 1090:
                if digits_in_group == 0:
                    slot_index += 3
                if digits_in_group == 1:
                    digit_slots[slot_index + 1] = digit_slots[slot_index]
                    digit_slots[slot_index] = 0
                    slot_index += 2
                if digits_in_group == 2:
                    slot_index += 1
                if slot_index == 1:
                    digits_in_group = 0
                if digits_in_group > 0:
                    digits_in_group = 0
            elif action == 1100:
                if digits_in_group == 0:
                    slot_index += 2
                    digits_in_group = 2
                if digits_in_group == 1:
                    digit_slots[slot_index + 1] = digit_slots[slot_index]
                    digit_slots[slot_index] = 0
                    slot_index += 1
                    digits_in_group = 2
    else:
        slot_index = 0
        for column_index in range(10, 81):
            column_char = padded_card[column_index]
            slot_index += 1
            if slot_index > 10:
                break
            if column_char in _DIGIT:
                digit_slots[slot_index] = _DIGIT[column_char]
            if column_char == ' ':
                slot_index -= 1
            if column_char == '-' and na != 5:
                na += 1
            if column_char == '.':
                slot_index -= 1
    chord_stations = [0.0]
    delx = 0.00100
    station_count = 1
    for index in range(2, 61):
        station_count = index
        chord_stations.append(chord_stations[-1] + delx)
        if chord_stations[-1] >= .01:
            delx = .01
        if chord_stations[-1] > .29:
            delx = .05
        if chord_stations[-1] > .79:
            delx = .02
        if chord_stations[-1] >= 1.0:
            break
    chord_stations[station_count - 1] = 1.0
    names = ('i', 'j', 'k', 'ii', 'jj', 'kk', 'iii', 'jjj', 'kkk', 'lll')
    return {'na': na, 'digits': dict(zip(names, digit_slots[1:11])),
            'x': chord_stations, 'l': station_count, 'method': 'legacy_decode'}


def airfol(card: str, surface_in: Mapping[int, float],
           stale_na: int = 0) -> Dict[str, object]:
    """Translate AIRFOL: decode the NACA card and compute the section.

    Dispatches to COORD4, CORD4M, COORD5, CORD5M, COORD1, COORD6 or, for a
    supersonic card, CORDSP then XYCORD (a biconvex card fixing its
    maximum-thickness digits at 0.5).  ``surface_in`` holds the named
    surface's input words CORDSP fills."""
    decoded = decode(card, stale_na)
    x, digits = decoded['x'], dict(decoded['digits'])
    routine = {1: coord4, 2: cord4m, 3: coord5, 4: cord5m, 5: coord1,
               6: coord6}.get(decoded['na'])
    if routine is not None:
        section = routine(digits, x)
    else:
        if digits['i'] != 2:
            pass
        else:
            digits.update({'j': 5, 'k': 0, 'ii': 0})
        grid = x + [x[-1]] * (60 - len(x))
        cordsp_result = cordsp(digits, grid, surface_in)
        section = xycord(x, cordsp_result['yu'][:len(x)],
                         cordsp_result['yl'][:len(x)])
        section['cordsp'] = cordsp_result
    section.update({'na': decoded['na'], 'digits': digits, 'x': x,
                    'l': decoded['l']})
    return section
