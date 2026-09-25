"""
Flap spanwise loading: AGENR and GDELTA.

GDELTA solves the four-point spanwise loading of a flapped wing (the
DeYoung method behind Datcom's G/delta charts): AGENR forms the 4x4
influence matrix from the span-to-chord ratios at four stations and the
sweep, SIMUL4 solves it for four flap spans, and the four solutions are
spread to nine span stations with the method's interpolation constants
(``CN``) and read at the flap's inboard and outboard edges.  For an
all-moving horizontal tail (``ASYFP``) it returns the four loadings at
the root for the tail's own planform instead.

The tables were extracted from the source DATA statements by parsing, with
``tools/fortran_data.py``.

Reference: datcom-legacy/datcom_2000/agenr.f, gdelta.f
"""

import math
from typing import Dict, Mapping, Sequence

import numpy as np

from pydatcom.utils.constants import DEG
from pydatcom.utils.legacy_numeric import simul4, tbfunx
from pydatcom.utils.legacy_tables import tlinex

_C1 = [-0.07612, 1e-05, 0.21677, 0.5412, -0.29289, -0.21677, 1e-05,
       0.32443, -0.61732, -0.5412, -0.32443, 1e-05, -1.0, -0.92388,
       -0.70711, -0.38267]
_C2 = [1.92388, 1.84776, 1.63099, 1.30656, 1.70711, 1.63099, 1.41422,
       1.08979, 1.38268, 1.30655, 1.08979, 0.76536, 1.0, 0.92388, 0.70711,
       0.38267]
_F1 = [5.2262, 1.036, 0.0, 0.11208]
_F2 = [1.91433, 2.8284, 0.91418, 0.0]
_F3 = [0.0, 1.1944, 2.1634, 1.5772]
_F4 = [0.14645, 0.0, 0.85357, 2.0]
_C3 = [0.92388, 0.70711, 0.38268, 1e-05]

_X1 = [0.0, 0.195, 0.556, 0.831, 1.0]
_SPU = [0.0, 0.1423, 0.2817, 0.4153, 0.5407, 0.6549, 0.7557, 0.8412,
        0.9097, 0.9595, 0.9898, 1.0, 0.5, 0.75]
_SD = [0.0, 0.383, 0.707, 0.924, 1.0]
_SPT = [0.0, 0.195, 0.383, 0.556, 0.707, 0.831, 0.924, 0.981, 1.0]
_EQ = [-0.017, 0.029, -0.014, 0.988, 0.032, -0.003, 0.994, 0.976, -0.041,
       1.021, 0.955, 1.04, 1.0, 1.0, 1.0, 1.0]
_CN = [0.05775384, -0.1585488, 0.61908, 0.5162876, -0.1332188, 0.5903109,
       0.6376293, -0.09812725, 0.528724, 0.6979618, -0.2243549,
       0.05177605, 0.8184079, -0.2989922, 0.1550801, 0.03851272,
       0.04523051, -0.1196368, 0.4375114, 0.6320983, -0.1533014,
       0.6545052, 0.6621265, -0.1765276, 0.5301887, 0.6743488, -0.2033035,
       0.08116567, 0.8319939, -0.2928609, 0.142265, 0.06120639,
       0.03238602, -0.07837254, 0.4212096, 0.6365698, -0.1345339,
       0.5254995, 0.7812855, -0.2178882, 0.5871445, 0.6832385, -0.3022919,
       0.1264219, 0.7802445, -0.2512727, 0.1793873, 0.08073133]


def agenr(boak: Sequence[float], sb: float) -> list:
    """Translate AGENR: the 16 influence coefficients for GDELTA.

    Args:
        boak: ``BOAK(4)``, span over chord at the four stations.
        sb: ``SB``, the sweep, degrees.

    Returns:
        ``A(16)``, equation ``i``'s coefficient on unknown ``m`` at
        ``4(i-1)+m`` (zero-based ``4i+m``), as SIMUL4 reads it.
    """
    tansb = math.tan(DEG * sb)
    al = [0.0] * 16
    boco16 = [0.0] * 4
    for coeff_index in range(16):
        boch_block = coeff_index // 4
        boch = float(boak[boch_block])
        boco16[boch_block] = -boch / 16.
        rcplbc = 1.0 / boch
        boc2 = boch**2
        boctn = boch * tansb
        denom = 1. + _C3[boch_block] * 2.0 * boctn
        tmp = (1. + _C1[coeff_index] * boctn)**2
        al[coeff_index] = (
            rcplbc / _C1[coeff_index] *
            (math.sqrt(tmp + boc2 * _C1[coeff_index]**2) - 1.)
            - rcplbc / _C2[coeff_index] *
            (math.sqrt(tmp + boc2 * _C2[coeff_index]**2) / denom - 1.0)
            - 2. * tansb *
            math.sqrt((1. + _C3[boch_block] * boctn)**2 +
                      boc2 * _C3[boch_block]**2) / denom)
    a = []
    for equation_index in range(4):
        coeff_base = 4 * equation_index
        b = boco16[equation_index]
        f = [-2. * _F1[equation_index], -2. * _F2[equation_index],
             -2. * _F3[equation_index], -2. * _F4[equation_index]]
        f[equation_index] = -f[equation_index]
        a.append(f[0] + b * (2.6131 * al[coeff_base] + 2. * (
            -.70711 * al[coeff_base + 1] - .76537 * al[coeff_base + 2] +
            .20711 * al[coeff_base + 3])))
        a.append(f[1] + b * (-1.4142 * al[coeff_base] + 2. * (
            1.8478 * al[coeff_base + 1] - .50000 * al[coeff_base + 2] -
            .76537 * al[coeff_base + 3])))
        a.append(f[2] + b * (1.0824 * al[coeff_base] + 2. * (
            -1.2071 * al[coeff_base + 1] + 1.8478 * al[coeff_base + 2] -
            .70711 * al[coeff_base + 3])))
        a.append(f[3] + b * (-.5 * al[coeff_base] +
                             1.0824 * al[coeff_base + 1] -
                             1.4142 * al[coeff_base + 2] +
                             2.6131 * al[coeff_base + 3]))
    return a


def calculate_gdelta(efi: float, efo: float, boch: Sequence[float],
                     sb: float, asyfp: bool = False,
                     tail: Mapping[str, float] = None) -> Dict[str, object]:
    """Translate GDELTA: the flap spanwise loading coefficient G/delta.

    Args:
        efi, efo: ``EFI``, ``EFO``, the flap's inboard and outboard span
            stations (fractions of the semispan).
        boch: ``BOCH(4)``, span over chord at the four stations (replaced
            for ``asyfp``).
        sb: ``SB``, the sweep, degrees.
        asyfp: The all-moving tail path.
        tail: For ``asyfp``: ``tante``, ``tanle`` (``AHT(80)``, ``(62)``),
            ``bsto2`` (``HTIN(3)``) and ``crh`` (``AHT(10)``).

    Returns:
        For a flap: ``gd1`` (the full-span loading at the 14 stations
        ``SPU``, the last two the flap edges), ``gd2`` and ``gd3`` (the
        loading of flaps from the root to ``EFI`` and to ``EFO``), ``fgc``
        (the four nine-point curves) and ``gd`` (the SIMUL4 solutions as
        rearranged); for ``asyfp``, ``gdh`` (``TCD(43..46)``) and the
        ``boch`` it formed.
    """
    boch = [float(v) for v in boch]
    if asyfp:
        t = {k: float(v) for k, v in tail.items()}
        arg1 = (t['tante'] - t['tanle']) * t['bsto2']
        boch = [2. * t['bsto2'] / (t['crh'] + _SD[3 - span_station] * arg1)
                for span_station in range(4)]
    a = agenr(boch, sb)
    gd = []
    for block_start in range(0, 16, 4):
        gd += simul4(a, _EQ[block_start:block_start + 4])
    gd.append(0.0)
    if asyfp:
        return {'gdh': [gd[15], gd[14], gd[13], gd[12]], 'boch': boch,
                'gd': gd, 'method': 'legacy_gdelta'}
    gi = []
    for row in range(3):
        for col in range(4):
            gi.append(sum(_CN[4 * col + 16 * row + coeff] * gd[4 * row + coeff]
                          for coeff in range(4)))
    fgc = [0.0] * 36
    for curve in range(4):
        for point in range(4):
            fgc[4 * curve + point] = gd[4 * curve + 3 - point]
    gd[:16] = fgc[:16]
    for l in (2, 4, 6, 8):
        fgc[l + 26] = tbfunx(_SD, gd[12:17], _SPT[l - 1], 0, 0)[0]
    ll = 0
    j = 1
    for l in range(16):
        fgc[ll] = gd[l]
        if l < 12:
            fgc[ll + 1] = gi[l]
        ll += 2
        j += 1
        if j >= 5:
            fgc[ll] = 0.
            ll += 1
            j = 1
    for i in (8, 17, 26, 35):
        fgc[i - 1] = .25 * fgc[i - 2]
    spu = list(_SPU)
    spu[12], spu[13] = efi, efo
    curves = [fgc[9 * curve:9 * curve + 9] for curve in range(4)]
    zrx = np.zeros((12, 5))
    for curve in range(4):
        zrx[:, curve + 1] = [
            tbfunx(_SPT, curves[curve], spu[station], 0, 0)[0]
            for station in range(12)]
    gd1, gd2, gd3 = [], [], []
    for station in range(14):
        gd1.append(tbfunx(_SPT, curves[3], spu[station], 0, 0)[0])
        gd2.append(float(tlinex(_X1, spu[:12], zrx, efi, spu[station])))
        gd3.append(float(tlinex(_X1, spu[:12], zrx, efo, spu[station])))
    return {'gd1': gd1, 'gd2': gd2, 'gd3': gd3, 'fgc': fgc, 'gd': gd,
            'boch': boch, 'method': 'legacy_gdelta'}
