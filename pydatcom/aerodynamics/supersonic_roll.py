"""
Supersonic roll and yaw of controls: SPRYAW.

SPRYAW forms the supersonic rolling effectiveness of a plain trailing-edge
flap (from DFLCON's derivatives) and its yawing moment (Figure 6.2.2.1-13),
the spoiler's rolling and yawing moments (Figures 6.2.2.1-14 and
6.2.1.1-30), and the rolling moment of a differentially deflected
all-moving horizontal tail (Figure 4.3.1.2-12A1/A2 with SUPWBT's body
vortex terms).

The tables were extracted from the source DATA statements by parsing, with
``tools/fortran_data.py``.

Reference: datcom-legacy/datcom_2000/spryaw.f
"""

import math
from typing import Dict, Mapping

from pydatcom.aerodynamics.supersonic_control import calculate_dflcon
from pydatcom.utils.constants import RAD
from pydatcom.utils.legacy_numeric import tbfunx
from pydatcom.utils.legacy_tables import tlinex_flat

_X113A1 = [0., 5., 10., 15.]
_X113A2 = [0., .01]
_Y2113A = [-.02, -.02, -.02, -.0176, -.02, -.01068, -.02, -.0004]
_X113B1 = [0., 5., 10., 15.]
_X113B2 = [0., .01]
_Y2113B = [0., 0., 0., 1.2, 0., 5.15, 0., 11.2]
_X113C1 = [0., 2., 4., 6.]
_X113C2 = [0., .06]
_Y2113C = [0., 0., 0., 1.65, 0., 3.38, 0., 7.65]
_X113E1 = [0., 45.]
_X113E2 = [0., .8]
_Y2113E = [0., .00795, 0., .00474]
_X114F1 = [1.]
_X114F2 = [3.2, 2.6, 2.2, 2., 1.8, 1.7, 1.6, 1.4]
_Y2114F = [2.15, 1.9, 1.85, 2., 2.25, 2.5, 2.75, 4.]
_X11301 = [0., 45.]
_X11302 = [0., 1.9]
_Y21130 = [0., .02, 0., .0184]
_X114A1 = [0., .5, 1.]
_X114A2 = [.4, 1.]
_Y2114A = [8.7, 5.5, 6.65, 5.5, 5.5, 5.5]
_X114E1 = [1.]
_X114E2 = [3.2, 2.8, 2.4, 2.2, 2.0, 1.8, 1.6]
_Y2114E = [.35, .75, 1.35, 1.75, 2.25, 2.75, 3.5]
_X12A1 = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
_Y12A1 = [1.0, 0.97, 0.95, 0.94, 0.94, 0.94, 0.94, 0.95, 0.96, 0.98, 0.99]
_X12A2 = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
_Y12A2 = [0.0, 0.11, 0.21, 0.31, 0.41, 0.51, 0.6, 0.7, 0.8, 0.9, 1.0]


def calculate_spryaw(data: Mapping[str, object]) -> Dict[str, object]:
    """Translate SPRYAW: supersonic control roll and yaw.

    Args:
        data: By name: ``alpha`` (``FLC(23..)``), ``mach``, ``f``
            (``/FLAPIN/`` 11-18, 19.., 29.., 39..), ``win`` (``WINGIN``
            2, 3, 4), ``a`` (``A`` 3, 26, 28, 62, 70, 76, 77), ``htin``
            (3, 4), ``sh`` (``AHT(3)``), ``cnahs`` (``HT(101)``),
            ``cladeg`` (``WING(101)``), ``ivbh``, ``gamma`` (``STP(133..)``,
            ``(111..)``), ``sref``, ``blref``, ``spr`` (``SPR(1..59)``),
            and ``sae025``, the local SPRYAW reads without setting.

    Returns:
        ``returned`` (``STYPE`` 3, or a subsonic hinge line or trailing
        edge), ``spr`` as left, ``clrlal`` (``HT(211..)``, flaps),
        ``cnywal`` (``BODY(201..)``, 10 per angle), ``clrlsp``, ``cnywsp``
        (``HT(211..)``, ``(221..)``, spoilers), ``clrlht``
        (``WING(201..)``, all-moving tail).

    Notes:
        Kept as executed: the control is taken as inboard when its outer
        edge is at the tip (the test is reversed from SSHING's); ``TANTE``
        is ``A(77)``, the trailing-edge sweep in radians, divided by beta
        as if a tangent; the spoiler's Figure 6.2.1.1-30 lookup reads
        ``SAE025``, which is never set.
    """
    f = {int(k): float(v) for k, v in data['f'].items()}
    w = {int(k): float(v) for k, v in data['win'].items()}
    a = {int(k): float(v) for k, v in data['a'].items()}
    h = {int(k): float(v) for k, v in data['htin'].items()}
    spr = [0.0] + [float(v) for v in data['spr']]
    out: Dict[str, object] = {'returned': False}
    stype = f[18]
    if stype == 3.:
        out.update({'returned': True, 'spr': spr[1:]})
        return out
    sspn, sspne, bstro2 = w[4], w[3], w[2]
    aloci, aloco, cfi, cfo = f[14], f[15], f[12], f[13]
    inbord = not aloco / sspne < 0.99
    mach = float(data['mach'])
    ndelta = int(f[16] + .5)
    beta = math.sqrt(mach**2 - 1.)
    spr[1] = beta
    sr, blref = float(data['sref']), float(data['blref'])
    alpha = [float(v) for v in data['alpha']]
    deltal = [f[19 + slot] for slot in range(ndelta)]
    deltar = [f[29 + slot] for slot in range(ndelta)]
    if stype == 5.:
        sspnh, sspneh = h[4], h[3]
        rads = sspnh - sspneh
        ratio = rads / sspnh
        yhs = .4 * sspneh + rads
        khb = tbfunx(_X12A1, _Y12A1, ratio, 0, 0)[0]
        kbh = tbfunx(_X12A2, _Y12A2, ratio, 0, 0)[0]
        spr[9], spr[10], spr[11] = khb, kbh, yhs
        clrlht = {}
        for angle_index in range(len(alpha)):
            for deflection_index in range(ndelta):
                clrls = (.006108 * (float(data['ivbh'][angle_index]) *
                                    float(data['gamma'][angle_index]) *
                                    (rads / sspneh) + (kbh + khb)) *
                         float(data['cnahs']) * yhs * float(data['sh']) /
                         (blref * sr))
                clrlht[10 * angle_index + deflection_index] = (
                    clrls * (deltal[deflection_index] -
                             deltar[deflection_index]))
        out.update({'spr': spr[1:], 'clrlht': clrlht})
        return out
    spr[2] = 2. / beta
    spr[3] = (2.4 * mach**4 - 4. * beta**2) / (2. * beta**4)
    trtoe = a[28] if aloci >= sspn - bstro2 else a[26]
    lamhl = math.atan(((aloco - aloci) * math.tan(a[76] / RAD) + cfi - cfo) /
                      (aloco - aloci))
    spr[4] = lamhl
    spr[14] = math.tan(lamhl)
    spr[5] = 2. * math.atan(f[11] / math.cos(lamhl)) * RAD
    spr[6] = 1. - (spr[3] / spr[2]) * (spr[5] / RAD)
    spr[7] = 2. * (((cfo + cfi) / 2.) * (aloco - aloci))
    amgcln = 1.
    spnspo = (aloco - aloci) / sspne
    geom = (aloci + aloco) / sspne
    bef = 2. * (aloco - aloci)
    spr[32] = a[62] / beta
    spr[33] = spr[14] / beta
    spr[34] = a[77] / beta
    sw = a[3]
    if stype == 4.:
        if spr[33] >= 1. or spr[34] >= 1.:
            out.update({'returned': True, 'spr': spr[1:]})
            return out
        trtofl = cfo / cfi
        if trtofl == 1.:
            aforlf = ((2. * (aloco - aloci))**2) / spr[7] * beta
            taperd = False
        else:
            aforlf, taperd = trtofl, True
        dfl = calculate_dflcon(spr[33], spr[34], aforlf, beta, inbord,
                               taperd)
        if dfl is not None:
            spr[12], spr[13] = dfl['cld'], dfl['clld']
            spr[17], spr[18] = dfl['cmd'], dfl['chd']
        bcld1, bcld2 = spr[12], spr[13]
        spr[8] = spr[6] * (spr[7] / sr) * (bcld1 / beta) * .5 * (
            (aloci / blref) + (bef / (2. * blref)) * (bcld1 / bcld2))
        clrlal = [spr[8] * (deltal[deflection_index] -
                             deltar[deflection_index]) / 2.
                  for deflection_index in range(ndelta)]
        scale = sw * 2. * sspne / (sr * blref)
        cldg = abs(float(data['cladeg']))
        clrl = abs(spr[8])
        cnywal = {}
        for angle_index, al in enumerate(alpha):
            dfg = tlinex_flat(_X113C1, _X113C2, _Y2113C, abs(al), cldg, 0, 0,
                              2, 1)
            for deflection_index in range(ndelta):
                tempo = [0.0, 0.0]
                for side in range(2):
                    delflp = abs(deltal[deflection_index] if side == 0 else
                                 deltar[deflection_index])
                    cld1 = tlinex_flat(_X113A1, _X113A2, _Y2113A, delflp, clrl,
                                       0, 0, 2, 1)
                    abc = tlinex_flat(_X113B1, _X113B2, _Y2113B, delflp, clrl,
                                      0, 1, 2, 1)
                    cld2 = (10. - abc + cld1 * 500.) / 500.
                    subax1 = dfg + cld2 * (2.5 + dfg * 1.5) / (-.005)
                    subax2 = geom * subax1 * .1 * .5
                    tempo[side] = tlinex_flat(_X113E1, _X113E2, _Y2113E, a[70],
                                              subax2, 0, 0, 2, 1)
                cnywal[10 * angle_index + deflection_index] = (
                    (tempo[1] - tempo[0]) * scale)
        out.update({'spr': spr[1:], 'clrlal': clrlal, 'cnywal': cnywal})
        return out
    sae025 = float(data['sae025'])
    mach_lookup = mach
    clrlsp, cnywsp = [], []
    for deflection_index in range(ndelta):
        delsoc = f[39 + deflection_index]
        vert = tlinex_flat(_X114A1, _X114A2, _Y2114A, trtoe, geom, 0, 1, 0, 1)
        nugeom = geom * delsoc * spnspo * 10.
        trnsgf = vert * nugeom * .5
        arbit = tlinex_flat(_X114F1, _X114F2, _Y2114F, amgcln, mach_lookup, 0,
                            1, 0, 2)
        finaly = trnsgf * (.2 + arbit * .2)
        taktim = tlinex_flat(_X11301, _X11302, _Y21130, sae025, finaly, 0, 0,
                             1, 1)
        clrlsp.append(taktim * sw * 2. * sspne / (blref * sr))
    vert = tlinex_flat(_X114A1, _X114A2, _Y2114A, trtoe, geom, 0, 1, 0, 1)
    arbit = tlinex_flat(_X114E1, _X114E2, _Y2114E, amgcln, mach_lookup, 0, 1,
                        0, 1)
    for deflection_index in range(ndelta):
        vert2 = (f[39 + deflection_index] * spnspo) / 1.2
        nugeom = vert2 * (geom + .1) * 12. - .1
        trnsgf = vert * (.05 + .5 * nugeom)
        pastim = .2 * trnsgf * (1. + arbit) * .005
        cnywsp.append(pastim * 2. * sspne * sw / (blref * sr))
    out.update({'spr': spr[1:], 'clrlsp': clrlsp, 'cnywsp': cnywsp})
    return out
