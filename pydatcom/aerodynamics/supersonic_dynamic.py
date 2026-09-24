"""
Supersonic wing acceleration derivatives: CALCA.

CALCA evaluates the linear-theory closed forms of Sections 7.1.1.3 and
7.1.2.3 for a supersonic wing with a subsonic leading edge: the pitching
and plunging acceleration derivatives (CL_alpha-dot, Cm_alpha-dot) and the
pitch-rate terms (CL_q, Cm_q), through the elliptic-integral factors E, G,
K and E' (the four INTERX lookups) and the planform parameters
``K6 = beta/tan(LE)``, ``W`` and ``SN``.

The long expressions were converted mechanically from the source text
(continuations joined, names lowercased) and are evaluated on NumPy
scalars, so a negative square root or a fractional power of a negative
number gives NaN as the compiled FORTRAN does rather than raising.

Reference: datcom-legacy/datcom_2000/calca.f
"""

from typing import Dict, Mapping, Optional

import numpy as np

from pydatcom.utils.constants import PI as _PI, RAD as _RAD
from pydatcom.utils.legacy_interp import interx
from pydatcom.utils.math_utils import arcsin as _arcsin

_XEGIN = [0., .025, .05, .1, .15, .2, .3, .4, .5, .6, .7, .8, .9, 1.0]
_YEGIN = [1., 1., .996, .984, .969, .951, .912, .868, .827, .782, .743,
          .706, .669, .637]
_XGM = [0., .05, .1, .2, .3, .4, .5, .6, .7, .8, .9, 1.]
_YGM = [1., .984, .96, .883, .802, .728, .657, .595, .540, .495, .459,
        .425]
_XKGIN = [0., .025, .05, .1, .15, .2, .25, .3, .35, .4, .45, .5, .6, .7,
          .8, .9, 1.]
_YKGIN = [.14, .198, .229, .272, .305, .334, .358, .381, .401, .424, .444,
          .464, .501, .537, .571, .605, .637]
_XEDM = [0., .025, .05, .1, .15, .2, .3, .4, .5, .6, .7, .8, .9, 1.]
_YEDM = [1., 1., .997, .985, .969, .952, .913, .87, .827, .782, .743,
         .707, .669, .637]


def _look(x, y, q):
    return np.float64(interx(1, x, [float(q)], [len(x)], y, lind=len(x),
                             lx1l=1, lx1u=1))


def arcsin(value):
    return np.float64(_arcsin(float(value)))


def calculate_calca(mach: float, a: Mapping[int, float],
                    wingin: Mapping[int, float]) -> Optional[Dict[str, object]]:
    """Translate CALCA: supersonic wing acceleration and pitch-rate
    derivatives.

    Args:
        mach: ``FLC(IM+2)``.
        a: The wing's ``A`` words 3, 7, 9, 27 and 58 (``A(58)`` the
            leading-edge sweep in degrees).
        wingin: ``WINGIN`` 1, 3, 4.

    Returns:
        ``dyn``: ``DYN(16)`` (CL_alpha-dot), ``(22)`` (CL_q), ``(26)``
        (Cm_q) and ``(27)`` (Cm_alpha-dot), with the intermediate factors;
        ``None`` where the source returns without setting them (``A2`` or
        ``A3`` negative, or ``SN`` outside [0, 1)).
    """
    pi, rad = np.float64(_PI), np.float64(_RAD)
    g = {int(k): np.float64(v) for k, v in a.items()}
    win = {int(k): np.float64(v) for k, v in wingin.items()}
    mach = np.float64(mach)
    dyn = {}
    with np.errstate(all='ignore'):
        beta = np.sqrt(mach**2-1.0)
        b2 = beta
        pio2 = pi/2.
        d = (win[4]-win[3])*2.
        b = 2.0*win[4]
        db = d/b
        k6 = beta/np.tan(g[58]/rad)
        cr = g[9]
        s = (cr+win[1])*win[4]
        se = g[3]
        are = g[7]
        lame = g[27]
        gg = beta/(np.tan(g[58]/rad)+(4.0/are)*((lame-1.0)/(lame+1.0)))
        egin = _look(_XEGIN, _YEGIN, gg)
        gm = _look(_XGM, _YGM, k6)
        kgin = _look(_XKGIN, _YKGIN, gg)
        edm = _look(_XEDM, _YEDM, k6)
        epm = 1.0/edm
        epg = 1.0/egin
        kpg = 1.0/kgin
        w = (4.0*k6)/(are*b2*(1.0+lame))
        sn = 1.0-(1.0-lame)*w
        a1 = (1.0+k6)
        a2 = (1.0+sn)
        a3 = (a1*a2+w*(k6-1.0))*(w+sn-1.0)
        a4 = 1.0-k6
        a5 = a1*sn+k6*w
        a6 = a2+w
        a7 = a2-w
        a8 = a1-w
        a9 = sn+k6
        if a2 < 0.0 or a3 < 0.0 or sn >= 1.0 or sn < 0.0:
            return None
        a4sq = a4**2
        a1sr = np.sqrt(a1)
        a2sr = np.sqrt(a2)
        a3sr = np.sqrt(a3)
        a4sr = np.sqrt(a4)
        a31p5 = a3*a3sr
        a11p5 = a1*a1sr
        a3a1sr = a3sr*a1sr
        a42p5 = a4sq*a4sr
        a21p5 = a2*a2sr
        a41p5 = a4*a4sr
        a43pr = a4sq*a4*a4sr
        sn2 = sn**2
        sn3 = sn*sn2
        sn2m1 = sn2-1.
        sn2m1s = np.sqrt(-sn2m1)
        cm1 = 3.0*(1.0-pio2/kpg)
        cm2 = w*cm1
        cm3 = (epg-gg**2*kpg)/(kpg-epg)
        sn4 = sn2**2
        sn5 = sn4*sn
        cm4 = 4.0*sn5-12.0*sn3+23.0*sn
        a10 = arcsin((a1*(sn2m1)+w*(1.0+sn*k6))/(w*a9))-arcsin(sn)
        a11 = arcsin((w-1.0-k6*(w+sn))/a9)+pio2
        a12 = (sn+w)*(sn-k6)-a9+2.0*(w-1.0)
        a13 = arcsin(a12/(a6*a9))+pio2
        oo = ooo = o4 = np.float64(0.)
        if not (gg >= 1. or gg <= 0.):
            oo = ((8.0*k6*w**2)/gg)*(cm1+(sn/w)*(2.0-(epg/kpg)-cm3)+(sn2/(2.0*w**2))*(1.0-((pi*(1.0-gg**2))/(4.0*(kpg-epg)))))
            ooo = (2.0*k6/(3.0*gg))*(cm2+sn*(1.0-cm3))
            o4 = (k6/gg)*(sn*(1.0-epg/kpg)+cm2)
        q = 3.0*gm*((w**4/(sn2m1)**3)*((a8**2*(2.0*sn**5+12.0*sn**3+sn)+a5*(a8*(-27.0*sn**2-3.0)+a5*cm4))/(w**3*a9**3)*(a3a1sr)+a10*(3.0*(2.0*sn**2+3.0)/(sn2m1s))-cm4)+(a1sr/a4**3)*((3.0*(2.0*k6**2+3.0)/a4sr)*a11+((a8**2*(-2.0*k6**5-12.0*k6**3-k6)+a5*(a8*(-27.0*k6**2-3.0)+a5*(-4.0*k6**5+12.0*k6**3-23.0*k6)))/(a9**3*a1**2))*a3sr)-oo)
        c = ((8.0*a4*a9*(3.0+27.0*k6+37.0*k6**2+13.0*k6**3)+(12.0+14.0*k6+21.0*k6**2+13.0*k6**3)*(-3.0*a4*a8+5.0*k6*a9))/(3.0*a4**2*a9**4))*a3**1.5+((15.0*(8.0+6.0*k6+3.0*k6**2+3.0*k6**3))/(2.0*a4**3))*((((w-1.0)-k6*(w+sn))/(-a9**2))*a3**.5-(1.0/a4**.5)*a11)-((((4.0-2.0*sn+9.0*sn**2)*(3.0+2.0*k6)+5.0*k6**2*(sn-2.0))*(-6.0*a2*a8-5.0*a9*a7)-2.0*a2*a9*(8.0*w*(3.0+2.0*k6)*(9.0*sn-1.0)+20.0*k6*(k6*w+(2.0-sn)*a1)))/(6.0*a2**2*a9**4))*a3**1.5
        f = (((4.0-2.0*sn+9.0*sn**2)*(3.0+2.0*k6)+5.0*k6**2*(sn-2.0))*(5.0*a7**2+4.0*w*a2)+2.0*a2*a7*(8.0*w*(3.0+2.0*k6)*(9.0*sn-1.0)+20.0*k6*(k6*w+(2.0-sn)*a1))+16.0*w*a2**2*(9.0*w*(3.0+2.0*k6)-5.0*k6*a1))*((a12/(-a9**2))*a3**.5-((a6**2)/(2.0*a2**.5))*a13)
        kcmq = (-2.0*k6**2)/(are*b2**2*(3.0*w**2-3.0*w*(1.0-sn)+(1.0-sn)**2)**2)*(q+(8.0/(5.0*pi*a1**1.5))*(c-(1.0/(16.0*a2**3))*f))
        o = ((w*((2.0*sn**2-1.0)*(sn**2-1.0)*a1-w*(2.0*sn**2*(1.0-sn*k6)+(5.0*sn*k6+1.0))))/(3.0*(sn**2-1.0)**2*a9**2))*(a3*a1)**.5-((w**3*(2.0*sn**3-5.0*sn))/(3.0*(sn**2-1.0)**2))-(w**3/((sn**2-1.0)**2*(1.0-sn**2)**.5))*a10+(a1**.5/a4**2.5)*a11-(((-2.0*k6**3*sn+5.0*sn*k6+2.0*k6**2+1.0)-w*(2.0*k6**2-1.0)*(k6-1.0))/(3.0*(k6-1.0)**2*a9**2))*(a3*a1)**.5-ooo
        p = (((a9*(3.0-k6)*a1))/(3.0*a2*a4*a9**3))*a3**1.5+((3.0+k6)/a4**2)*((((w-1.0)-k6*(w+sn))/(-2.0*a9**2))*a3**.5-(1.0/(2.0*a4**.5))*a11)-((a7*(3.0*sn+2.0*sn*k6+k6**2)+2.0*a2*(w*(3.0+2.0*k6)-k6*a1))/a2**2)*((a12/(-8.0*a9**2))*a3**.5-(a6**2/(16.0*a2**.5))*a13)
        kclq = ((4.0*k6)/(b2*(3.0*w**2-3.0*w*(1.0-sn)+(1.0-sn)**2)))*(3.0*gm*o+(8.0/(pi*a1**1.5))*p)
        cma = ((2.0*k6)/(b2*(3.0*w**2-3.0*w*(1.0-sn)+(1.0-sn)**2)))*((1.0/epm)*(((w*sn*(k6-1.0)**2*(3.0*w*sn-w*k6*(sn**2-4.0)-sn*(sn**2-1.0)*a1)+k6*(sn**2-1.0)**2*(4.0*sn+3.0*k6-sn*k6**2+k6*w*a4))/((sn**2-1.0)**2*(k6-1.0)**2*a9**2))*(a3*a1)**.5+((w**3*(sn**2+2.0))/((sn**2-1.0)**2*(1.0-sn**2)**.5))*a10+((w**3*(sn**3-4.0*sn))/((sn**2-1.0)**2))-((a1**.5*(k6**2+2.0))/a4**2.5)*a11+o4)+(8.0/(pi*a1**.5))*(((-(a3**1.5))/(a9**2*a4*a2))-((a7*(2.0-sn)-2.0*w*a2)/a2)*((a12/(-8.0*a2*a9**2))*a3**.5-((a6**2)/(16.0*a2**1.5))*a13)-((2.0*k6+1.0)/a4)*((((w-1.0)-k6*(w+sn))/(-2.0*a4*a9**2))*a3**.5-(1.0/(2.0*a4**1.5))*a11)))
        clda2 = ((4.0*k6)/(b2**3*(3.0*w**2-3.0*w*(1.0-sn)+(1.0-sn)**2)))*(edm*(((2.0*w*(k6**2-1.0)*(sn**2-1.0)+(1.0+sn*k6)*(w**2*(k6-1.0)+a1*(sn**2-1.0)))/((sn**2-1.0)*(k6-1.0)*a9**2))*(a3*a1)**.5+(w**3/((sn**2-1.0)*(1.0-sn**2)**.5))*a10-((sn*w**3)/(sn**2-1.0))+(a1/a4)**1.5*a11)+(8.0/(pi*a1**.5))*((((k6-1.0)*a6*(2.0*a2*a8-a7*a9)+4.0*a2*((1.0+sn*k6)+w*(k6-1.0)))/(4.0*a2*a4*a9**2))*a3**.5+(a6**3/(8.0*a2**1.5))*a13-(1.0/a4**1.5)*a11))
        cmda2a = ((6.0*k6**2)/(are*b2**4*(3.0*w**2-3.0*w*(1.0-sn)+(1.0-sn)**2)**2))
        cmda2b = (3.0*edm*((1.0/(3.0*(sn**2-1.0)**3*(k6-1.0)**3*a1**1.5*a9**3))*(a8**2*(w*(k6-1.0)**3*a1**2*(-2.0*sn**5+sn**3+sn)-(sn**2-1.0)**3*(2.0*k6**5-k6**3-k6))+a5*a8*(w*(k6-1.0)**3*a1**2*(3.0*sn**4-3.0)-(sn**2-1.0)**3*(3.0*k6**4-3.0))+a5**2*(w*(k6-1.0)**3*a1**2*(2.0*sn**5-7.0*sn**3+5.0*sn)+(sn**2-1.0)**3*(2.0*k6**5-7.0*k6**3+5.0*k6)))*a3**.5+(a1**1.5/a4**2.5)*a11-((w**4*(2.0*sn**5-7.0*sn**3+5.0*sn))/(3.0*(sn**2-1.0)**3))-(w**4/(1.0-sn**2)**2.5)*a10))
        cmda2c = ((8.0/(5.0*pi*a1**.5))*((((5.0*sn*((sn+w)*a4+a1)+10.0*k6)*a3**1.5)/(a2*a4*a9**3))+((5.0*(2.0+k6))/a4**2)*(((w*a4-(1.0+sn*k6))/(-a9**2))*a3**.5-(1.0/a4**.5)*a11)-((5.0*(2.0*(3.0*w+2.0)*(w+sn**2)-(2.0+sn-10.0*sn*w-3.0*sn**3-3.0*sn*w**2)))/(8.0*a2**2))*(((-2.0*a2*a8+a9*a7)/(-a9**2))*a3**.5-(a6**2/(2.0*a2**.5))*a13)))
        cmda2 = cmda2a*(cmda2b+cmda2c)
        cmda2 = -cmda2
        cmda3a = (2.0*mach**2/(b2**2))
        cmda3b = ((3.0*k6**2)/(are*b2**2*(3.0*w**2-3.0*w*(1.0-sn)+(1.0-sn)**2)**2))*((3.0/epm)*((1.0/(3.0*a9**3*(sn**2-1.0)**3*(k6-1.0)**3*a1**1.5))*(a8**2*(w*(k6-1.0)**3*a1**2*(11.0*sn**3+4.0*sn**5)+(sn**2-1.0)**3*(11.0*k6**3+4.0*k6**5))-a5*a8*(w*(k6-1.0)**3*a1**2*(27.0*sn**2+3.0*sn**4)-(sn**2-1.0)**3*(27.0*k6**2+3.0*k6**4))+a5**2*(w*(k6-1.0)**3*a1**2*(2.0*sn**5-5.0*sn**3+18.0*sn)-(sn**2-1.0)**3*(5.0*k6**3-18.0*k6-2.0*k6**5)))*a3**.5-((w**4*(2.0*sn**5-5.0*sn**3+18.0*sn))/(3.0*(sn**2-1.0)**3))+((w**4*(3.0*sn**2+2.0))/((sn**2-1.0)**3*(1.0-sn**2)**.5))*a10+(((3.0*k6**2+2.0)*a1**.5)/a4**3.5)*a11))
        cmda3c = ((4.0/(5.0*pi*a1**.5))*(((a4**2*a9*(5.0*a7*(8.0-4.0*sn+3.0*sn**2)-16.0*w*a2*(2.0-3.0*sn))+6.0*a2*a4*a8*(a4*(8.0-4.0*sn+3.0*sn**2)-a2*(8.0+4.0*k6+3.0*k6**2))+2.0*a2**2*a9*(16.0+64.0*k6+4.0*k6**2-9.0*k6**3))/(3.0*a2**2*a4**2*a9**4))*a3**1.5-((5.0*(4.0*w*a2*(3.0*sn**2+4.0*w)+a7**2*(8.0-4.0*sn+3.0*sn**2)))/(8.0*a2**3))*(((-2.0*a2*a8+a7*a9)/(-a9**2))*a3**.5-(a6**2/(2.0*a2**.5))*a13)+((5.0*(4.0+4.0*k6+7.0*k6**2))/a4**3)*(((w*a4-(1.0+sn*k6))/(-a9**2))*a3**.5-(1.0/a4**.5)*a11)))
        cmda3 = cmda3a*(cmda3b+cmda3c)
        dyn[16] = clda2*b2**2
        dyn[27] = cmda2*b2**2
        dyn[22] = kclq+2.0*cma
        dyn[26] = kcmq+cmda3*b2**2/mach**2
    return {'dyn': {k: float(v) for k, v in dyn.items()},
            'egin': float(egin), 'gm': float(gm), 'kgin': float(kgin),
            'edm': float(edm), 'k6': float(k6), 'w': float(w),
            'sn': float(sn), 'gg': float(gg), 'method': 'legacy_calca'}
