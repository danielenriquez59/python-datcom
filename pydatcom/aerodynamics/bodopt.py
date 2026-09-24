"""
BODOPT and overlay M04O04: the asymmetric (cambered) body.

The main program takes this path instead of M06O06 when the body's upper
and lower contours ``ZU``/``ZL`` are given and Jorgensen's method was not
asked for.  BODOPT builds the same friction and base drag as BODYRT, then:

- the camber line from the contour midpoints, less the chord line joining
  its ends, and from it the zero-lift angle and zero-lift moment by thin
  airfoil integrals over the body length;
- a potential normal force ``KP*sin*cos`` and a viscous one ``KV*sin**2``
  past an onset angle, with ``KP`` and ``KV`` from the body's planform
  aspect ratio (Figures BA-1 and BA-2) and the onset from Figure
  4.2.1.2-37;
- moments about the centre of pressure of each, located from the
  planform.

The tables were extracted from the source DATA statements by parsing, with
``tools/fortran_data.py``, rather than by hand.

Reference: datcom-legacy/datcom_2000/bodopt.f, m04o04.f
"""

import numpy as np
from typing import Dict, Sequence
import logging

from pydatcom.aerodynamics.bodyrt import body_slope_pass
from pydatcom.interactions.body_vortex import getmax
from pydatcom.utils.constants import PI, RAD, UNUSED
from pydatcom.utils.legacy_interp import eqspce
from pydatcom.utils.legacy_numeric import tbfunx, trapz
from pydatcom.utils.legacy_tables import tlinex
from pydatcom.utils.table_lookup import fig26

logger = logging.getLogger(__name__)

_XBA1 = np.array([
    0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0,
])
_YBA1 = np.array([
    0.0, 0.72, 1.33, 1.83, 2.225, 2.575, 2.88, 3.125, 3.37,
])
_XBA2 = np.array([
    0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0,
])
_YBA2 = np.array([
    3.14, 3.14, 3.15, 3.19, 3.21, 3.626, 3.31, 3.9, 3.47,
])
_X1BA3 = np.array([
    1.0, 2.0, 3.0, 4.0, 100.0,
])
_X2BA3 = np.array([
    0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 100.0,
])
_YBA3 = np.array([
    34.85, 31.2, 27.55, 23.9, 20.25, 16.6, 12.95, 11.25, 10.3, 9.75, 9.4, 0.0,
    30.65, 25.9, 21.15, 16.4, 11.65, 9.1, 7.6, 6.7, 6.1, 5.7, 5.5, 0.0,
    26.7, 20.0, 13.3, 9.1, 7.5, 4.9, 4.0, 3.5, 3.1, 2.9, 2.75, 0.0,
    22.7, 16.0, 9.3, 7.2, 5.0, 3.95, 3.0, 2.4, 2.05, 1.85, 1.7, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
])
_X27M = np.array([
    0.0, 1.0, 2.0, 3.0,
])
_X27I = np.array([
    1.5778, 1.67221, 1.98509, 2.28874,
])

_FIG_BA3 = _YBA3.reshape(5, 12).T


def calculate_bodopt(x: Sequence[float], s: Sequence[float],
                     p: Sequence[float], r: Sequence[float],
                     zu: Sequence[float], zl: Sequence[float],
                     alpha_deg: Sequence[float], mach: float,
                     reynolds_per_length: float, sref: float, lref: float,
                     xcg: float, roughness: float = 1.6e-4
                     ) -> Dict[str, object]:
    """Translate BODOPT: asymmetric-body lift, drag and moment.

    Args:
        x, s, p, r: ``/BODYI/`` stations, areas, perimeters, half-widths.
        zu, zl: Upper and lower contour ordinates.
        alpha_deg: ``FLC(23)`` onward.
        mach: ``B(1)``.
        reynolds_per_length: ``FLC(M+42)``.
        sref: ``SREF``.
        lref: ``LREF``, the second word of ``/OPTION/``, i.e. ``CBARR``.
        xcg: ``XCG``.
        roughness: ``ROUGFC``.

    Returns:
        Dictionary with ``cd``, ``cl``, ``cm``, ``cn``, ``ca`` at each
        angle (the ``BODY`` slots), ``cla`` and ``cma`` (``BODY(101)``,
        ``BODY(121)``), ``alpha_zero_lift`` ``BD(81)``, ``cm0`` ``BD(62)``,
        ``cd_zero_lift`` ``BD(61)``, ``cdl`` ``BD(215)`` onward, and the
        planform quantities behind them.

    Notes:
        The fineness ratio and the width-to-depth ratio for Figure
        4.2.1.2-37 are taken at the *last* station, ``S(NX)``, although
        the source's comment reads "FIND SMAX ... AT DSDX=ZERO".  A body
        that closes to a point therefore always takes the defaults
        ``RATIO=3``, ``FR=10``.  Kept.

        The potential force works at ``alpha - alpha_0`` but the viscous
        force at the schedule angle itself, as the source writes it.

        Figure BA-2 (``KV``) reads ``3.19, 3.21, 3.626, 3.31, 3.9, 3.47``,
        an otherwise smooth rise broken by two entries that look like
        transposed digits of 3.26 and 3.39.  Chart data; kept.
    """
    x = np.asarray(x, dtype=float)
    s = np.asarray(s, dtype=float)
    p = np.asarray(p, dtype=float)
    r = np.asarray(r, dtype=float)
    zu = np.asarray(zu, dtype=float)
    zl = np.asarray(zl, dtype=float)
    nx = len(x)
    length = x[-1]

    # --- Zero-lift drag (same buildup as BODYRT) ---
    base = s[-1]
    _, max_area, _ = getmax(x, s)
    if base <= 0.3 * max_area:
        base = 0.30 * max_area

    cept, _ = tbfunx(_X27M, _X27I, mach, 0, 0)
    cutoff = (12.0 * length / roughness) ** 1.0482 * 10.0 ** cept
    reynolds = reynolds_per_length * length
    cf = fig26(min(reynolds, cutoff), mach)

    resampled = eqspce(x, r, p, s, 20)
    wetted = trapz(resampled['pe'], resampled['xe'], 1)[0]
    base_diameter = np.sqrt(base * 4.0 / PI)
    max_diameter = np.sqrt(max_area * 4.0 / PI)
    fineness = length / max_diameter
    friction = (
        cf * (1.0 + 60.0 / fineness ** 3 + 0.0025 * fineness) * wetted / max_area
    )
    base_drag = (
        0.029 * (base_diameter / max_diameter) ** 3
        / np.sqrt(friction) * max_area / sref
    )
    friction = friction * max_area / sref
    cd0 = friction + base_drag

    # --- Camber line: contour midpoints minus the chord through the ends ---
    z0 = 0.5 * (zu[0] + zl[0])
    znx = 0.5 * (zu[-1] + zl[-1])
    dzdx = (znx - z0) / (x[-1] - x[0])
    depth = zu - zl
    depth = np.where((depth == 0.0) & (np.arange(nx) != 0), UNUSED, depth)
    zprime = zl + 0.5 * depth - (x * dzdx + z0)
    zpol = zprime / length
    xol = x / length
    grid = eqspce(x, xol, zpol, s, 20)
    xole, zpole = grid['re'], grid['pe']
    alpha_integrand = np.zeros(20)
    cm_integrand = np.zeros(20)
    inner = slice(1, 19)
    root = np.sqrt(xole[inner] - xole[inner] ** 2)
    alpha_integrand[inner] = zpole[inner] / ((1.0 - xole[inner]) * root)
    cm_integrand[inner] = zpole[inner] * ((1.0 - 2.0 * xole[inner]) / root)
    alzero = trapz(alpha_integrand[:19], xole[:19], 1)[0]
    alpha_zero = -alzero / PI * RAD

    # --- Planform, BA-1/BA-2 factors, onset angle (Figure BA-3) ---
    planform = eqspce(x, r, xol, s, 20)
    xe, ye = planform['xe'], planform['re']
    xole = planform['pe']
    sp = 2.0 * trapz(ye, xe, 1)[0]
    cm0 = (
        2.0 * trapz(cm_integrand, xole, 1)[0] * (sp / sref) * (length / lref)
    )

    _, ymax, imax = getmax(x, r)
    arb = 4.0 * ymax * ymax / sp
    y95 = 0.96 * ymax
    im = imax - 1
    x95 = x[im] + (x[imax] - x[im]) * (y95 - r[im]) / (r[imax] - r[im])
    akp, _ = tbfunx(_XBA1, _YBA1, arb, 0, 2)
    akv, _ = tbfunx(_XBA2, _YBA2, arb, 0, 2)
    smax = s[-1]
    deff = np.sqrt(4.0 * smax / PI)
    if deff == 0.0:
        ratio, fr = 3.0, 10.0
    else:
        fr = length / deff
        ratio = 2.0 * r[-1] / depth[-1]

    tx95, ty95 = x.copy(), r.copy()
    n95 = None
    for i in range(nx):
        if x95 > x[i]:
            continue
        n95 = i + 1
        tx95[i], ty95[i] = x95, y95
        break
    if n95 is None:
        raise ValueError("BODOPT: the 0.96*YMAX station lies beyond the body")
    short = eqspce(tx95[:n95], ty95[:n95], np.ones(n95), np.ones(n95), 20)
    a95max = 2.0 * trapz(short['re'], short['xe'], 1)[0]
    alphv = tlinex(_X1BA3, _X2BA3, _FIG_BA3, ratio, fr, 2, 2, 2, 2)
    expon = (2.0 * y95 * x95) / a95max - 1.0
    xolcnp = 2.0 * expon / (2.0 * expon + 1.0)
    xovlcg = xcg / length
    delxol = xovlcg - xolcnp * (x95 / length)
    cumulative = trapz(ye, xe, 0)
    strips = np.diff(cumulative)
    moment = np.sum(2.0 * strips * (xe[:19] + (xe[1:] - xe[:19]) / 2.0))
    xovlca = (moment / sp) / length

    # --- Forces at each schedule angle ---
    alpha = np.asarray(alpha_deg, dtype=float)
    incidence_rad = (alpha - alpha_zero) / RAD
    sin_incidence = np.sin(incidence_rad)
    cos_incidence = np.cos(incidence_rad)

    cnp = akp * sin_incidence * cos_incidence
    cnv = np.where(
        alpha <= alphv,
        0.0,
        akv * np.sin((alpha - alphv) / RAD) ** 2,
    )
    cn = (cnp + cnv) * (sp / sref)
    cd = cn * sin_incidence + cd0

    moment_scale = (sp / sref) * (length / lref)
    cm = (
        cm0
        + cnp * delxol * moment_scale
        + cnv * (xovlcg - xovlca) * moment_scale
    )

    alpha_rad = alpha / RAD
    cl = (cn - cd * np.sin(alpha_rad)) / np.cos(alpha_rad)
    axial = cd * np.cos(alpha_rad) - cl * np.sin(alpha_rad)

    return {
        'cd': cd,
        'cl': cl,
        'cm': cm,
        'cn': cn,
        'ca': axial,
        'cla': float(akp * sp / (RAD * sref)),
        'cma': float(akp * delxol * moment_scale / RAD),
        'alpha_zero_lift': float(alpha_zero),
        'cm0': float(cm0),
        'cd_zero_lift': float(cd0),
        'cd_friction': float(friction),
        'cd_base': float(base_drag),
        'cdl': cn * sin_incidence,
        'planform_area': float(sp),
        'aspect_ratio': float(arb),
        'kp': float(akp),
        'kv': float(akv),
        'onset_angle': float(alphv),
        'method': 'legacy_bodopt',
    }


def calculate_m04o04(x, s, p, r, zu, zl, alpha_deg, mach,
                     reynolds_per_length, sref, cbar, blref, xcg,
                     roughness: float = 1.6e-4,
                     experimental: bool = False) -> Dict[str, object]:
    """Translate overlay M04O04: BODOPT and the body slope pass.

    The slope pass is M06O06's (see
    :func:`pydatcom.aerodynamics.bodyrt.body_slope_pass`).
    """
    opt = calculate_bodopt(
        x, s, p, r, zu, zl, alpha_deg, mach,
        reynolds_per_length, sref, cbar, xcg, roughness,
    )
    result = body_slope_pass(
        alpha_deg, opt['cd'], opt['cl'], opt['cm'],
        opt['cla'], opt['cma'], cbar, blref, experimental,
    )
    result.update({'bodopt': opt, 'method': 'legacy_m04o04'})
    return result
