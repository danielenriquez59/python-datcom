"""
Supersonic vertical-tail and ventral-fin sideslip: MASRAT, SUPLAV, SUPLAF.

MASRAT is Datcom method 3 for the apparent-mass ratio of a vertical panel
mounted on a body with a horizontal surface (Figures 5.3.1.1-25B to -25P,
R1/b = 1 only): a mid-height surface, a surface tangent on the panel's
side, or tangent on the opposite side, and a three-point interpolation in
height between those for any other position.  SUPLAV and SUPLAF use it to
form the panel's carryover factor on the wing-body and on the
wing-body-tail, and the resulting side-force, yawing-moment and rolling-
moment increments.  SUPLAF is SUPLAV on the ventral fin's blocks, with the
height measured the other way.

The tables were extracted from the source DATA statements by parsing, with
``tools/fortran_data.py``.

Reference: datcom-legacy/datcom_2000/masrat.f, suplav.f, suplaf.f
"""

import math
from typing import Dict, Mapping, Sequence

import numpy as np

from pydatcom.utils.constants import RAD
from pydatcom.utils.legacy_numeric import tbfunx

_R225U = [0.0, 0.2, 0.4, 0.6, 0.8]
# Figure 5.3.1.1-25B, C, E, F (mid-height surface).
_R1125B = [0.0, 0.2, 0.4, 0.6, 0.8, 1.0]
_F1125B = [1.0, 1.44, 1.97, 2.58, 3.25, 4.0]
_R1125C = [0.0, 0.1, 0.2, 0.4, 0.525, 0.61, 0.675, 0.74, 1.0]
_F1125C = [0.39, 0.96, 1.3, 1.9, 2.28, 2.56, 2.77, 3.0, 3.9]
_R1125E = [0.0, 0.23, 0.33, 0.41, 0.53, 0.68, 0.7, 0.745, 1.0]
_F1125E = [0.51, 1.23, 1.6, 1.83, 2.25, 2.75, 2.85, 3.0, 4.0]
_R1125F = [0.0, 0.22, 0.32, 0.4, 0.52, 0.66, 0.7, 0.735, 1.0]
_F1125F = [0.505, 1.19, 1.5, 1.79, 2.2, 2.68, 2.84, 3.0, 4.0]
# Figure 5.3.1.1-25M (opposite side).
_R1125M = [0.0, 1.0]
_F1125M = [0.5, 5.75]
# Figure 5.3.1.1-25H, I, J, K (same side).
_R1125H = [0.0, 0.1, 0.15, 0.2, 0.25, 0.3, 0.5, 1.0]
_F1125H = [0.5, 0.84, 0.93, 0.99, 1.01, 1.02, 1.06, 1.06]
_R1125I = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 1.0]
_F1125I = [0.5, 0.8, 1.0, 1.1, 1.18, 1.2, 1.22, 1.23]
_R1125J = [0.0, 0.2, 0.3, 0.4, 0.5, 0.7, 1.0]
_F1125J = [0.5, 1.02, 1.2, 1.33, 1.4, 1.45, 1.45]
_R1125K = [0.0, 0.3, 0.4, 0.5, 0.7, 1.0]
_F1125K = [0.5, 1.26, 1.43, 1.55, 1.66, 1.7]
_ZK = [-1.0, 0.0, 1.0]
# Figure 5.3.1.1-25A: K_V(B).
_XAMF = [0.0, 0.40, 1.0]
_YAMF = [0.5, 1.79, 4.0]


def _pline(x1, x2, y1, y2, xv):
    return y1 + (xv - x1) * (y2 - y1) / (x2 - x1)


def _lines(r, knee, y0, yk, y1):
    """Two straight lines through (0, y0), (knee, yk), (1, y1), clamped
    outside 0..1, as MASRAT writes Figures 25D, 25L, 25N, 25O and 25P."""
    if r < 0.0:
        return y0
    if r > 1.0:
        return y1
    if r <= knee:
        return _pline(0., knee, y0, yk, r)
    return _pline(knee, 1.0, yk, y1, r)


def _look(x, y, q):
    return tbfunx(x, y, q, 0, 0)[0]


def masrat(r2obo2: float, r1oba: float, zhor1: float) -> float:
    """Translate MASRAT: the apparent-mass ratio K by Datcom method 3.

    Args:
        r2obo2: ``R2OBO2``, the body radius over the surface semispan.
        r1oba: ``R1OBA``, the body radius over the panel span.
        zhor1: ``ZHOR1``, the surface height over the body radius: 0 mid,
            1 tangent on the panel's side, -1 on the opposite side.
    """
    def mid():
        an = [_look(_R1125B, _F1125B, r1oba),
              _look(_R1125C, _F1125C, r1oba),
              _lines(r1oba, .74, .52, 3.0, 4.),
              _look(_R1125E, _F1125E, r1oba),
              _look(_R1125F, _F1125F, r1oba)]
        return _look(_R225U, an, r2obo2)

    def same():
        an = [1.0, _look(_R1125H, _F1125H, r1oba),
              _look(_R1125I, _F1125I, r1oba),
              _look(_R1125J, _F1125J, r1oba),
              _look(_R1125K, _F1125K, r1oba)]
        return _look(_R225U, an, r2obo2)

    def opposite():
        an = [_lines(r1oba, .45, 1.03, 3., 6.),
              _look(_R1125M, _F1125M, r1oba),
              _lines(r1oba, .56, .5, 3.0, 5.2),
              _lines(r1oba, .63, .5, 3., 4.7),
              _lines(r1oba, .67, .5, 3., 4.5)]
        return _look(_R225U, an, r2obo2)

    if zhor1 == 0.0:
        return mid()
    if zhor1 == 1.0:
        return same()
    if zhor1 == -1.0:
        return opposite()
    return _look(_ZK, [opposite(), mid(), same()], zhor1)


def calculate_suplav(alpha_deg: Sequence[float], ventral: bool,
                     data: Mapping[str, object]) -> Dict[str, object]:
    """Translate SUPLAV (or SUPLAF, ``ventral=True``): the supersonic
    sideslip increments of a vertical panel.

    Args:
        alpha_deg: ``FLC(23..)``.
        ventral: SUPLAF: the ventral fin's blocks, the height sign
            reversed, and the ``DCVWHB`` slip below.
        data: The words read, by name: ``vtin`` (the panel's ``VTIN``, or
            ``VFIN``, words 3 and 4), ``sv`` (``SVWB``, ``SVB``, ``SVHB``,
            its ``VTIN(95..)``, ``(115..)``, ``(135..)`` words at the
            Mach index), ``avt`` (its ``A`` words 3, 30, 31, 62), ``cnav``
            (``SLA(31)``), ``span`` (``WINGIN(4)``), ``straight``
            (``WINGIN(15)`` is STRA), ``htin`` (3, 4), ``aht`` (30, 62),
            ``position`` (``xcg``, ``zw``, ``zcg``, ``zh``, ``alih``,
            ``xv``, ``zv``, ``vertup``), ``htpl``, ``blref``, and
            ``bwv``, ``bwhv``: the ``BWV``/``BWHV`` words 141, 161 and
            181.. the routine adds to.

    Returns:
        The ``SLA`` words (``rkvwb``, ``rkvb``, ``rkpvwb``, ``dcybv``,
        ``rkvhb``, ``zp``, ``rlp``), ``dcvwhb``, the panel block's
        ``vt141``, ``vt161`` and ``vt181`` (SUPLAF leaves ``vt181`` unset
        without a tail), and the updated ``bwv`` and ``bwhv`` words.

    Notes:
        Kept as executed: in SUPLAF ``DCVWHB`` shares ``VF(141)``, which is
        set to the tail-off ``DCYBV`` before the tail-on sums read it, so
        the ventral fin's wing-body-tail side force and moments repeat its
        wing-body ones.
    """
    p = data['position']
    vtin = {int(k): float(v) for k, v in data['vtin'].items()}
    avt = {int(k): float(v) for k, v in data['avt'].items()}
    htin = {int(k): float(v) for k, v in data['htin'].items()}
    aht = {int(k): float(v) for k, v in data['aht'].items()}
    sv = {k: float(v) for k, v in data['sv'].items()}
    span, blref = float(data['span']), float(data['blref'])
    cnav = float(data['cnav'])
    vertup = bool(p['vertup'])
    zh = float(p['zh']) - ((htin[4] - htin[3]) * aht[62] + aht[30]) * \
        math.sin(float(p['alih']) / RAD)
    svstar = avt[3]
    scale = 2.0 * span / blref
    r1 = vtin[4] - vtin[3]
    r1oba = r1 / vtin[4]

    def flip(z):
        return -z if (vertup if ventral else not vertup) else z

    arg1 = masrat(r1 / span, r1oba, flip(-float(p['zw']) / r1))
    rkvb = tbfunx(_XAMF, _YAMF, r1oba, 0, 1)[0]
    rkvwb = arg1
    rkpvwb = (rkvb * sv['svb'] + rkvwb * sv['svwb']) / svstar
    dcybv = -rkpvwb * cnav
    r = {'rkvwb': rkvwb, 'rkvb': rkvb, 'rkpvwb': rkpvwb, 'dcybv': dcybv}
    htpl = bool(data['htpl'])
    dcvwhb = None
    if htpl:
        rkvhb = masrat(r1 / htin[4], r1oba, flip(zh / r1))
        dcvwhb = dcybv - rkvhb * sv['svhb'] / svstar * cnav
        r['rkvhb'] = rkvhb
    bw = 2. * span
    delx = avt[62] * (vtin[4] - vtin[3])
    rlp = float(p['xv']) - float(p['xcg']) + delx + avt[30]
    zp = -float(p['zcg']) + vtin[4] - vtin[3] + avt[31] + float(p['zv'])
    vt141 = dcybv
    vt161 = -dcybv * rlp / blref
    bwv = {int(k): float(v) for k, v in data['bwv'].items()}
    bwhv = {int(k): float(v) for k, v in data['bwhv'].items()}
    bwv[141] += vt141
    bwv[161] += vt161
    vt181 = None
    if ventral and htpl:
        dcvwhb = vt141
    if data['straight']:
        dz = [(zp * math.cos(a / RAD) - rlp * math.sin(a / RAD)) / bw
              for a in alpha_deg]
        vt181 = None if ventral else [dcybv * d * scale for d in dz]
        for j, d in enumerate(dz):
            bwv[181 + j] += dcybv * d * scale
        if htpl:
            if not ventral:
                vt141 = dcvwhb
                vt161 = -dcvwhb * rlp / blref
            bwhv[141] += dcvwhb
            bwhv[161] -= dcvwhb * rlp / blref
            vt181 = [dcvwhb * d * scale for d in dz]
            for j, v in enumerate(vt181):
                bwhv[181 + j] += v
    r.update({'zp': zp, 'rlp': rlp, 'dcvwhb': dcvwhb, 'vt141': vt141,
              'vt161': vt161, 'vt181': vt181, 'bwv': bwv, 'bwhv': bwhv,
              'method': 'legacy_suplaf' if ventral else 'legacy_suplav'})
    return r
