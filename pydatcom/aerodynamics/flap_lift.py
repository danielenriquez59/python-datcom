"""
Flap lift increments and maximum lift: LIFTFP.

LIFTFP forms, for each flap deflection, the section lift increment over
four spanwise strips of the flap (plain, split, slotted, double-slotted,
Fowler, leading-edge and slat/Krueger devices, Figures 6.1.1.1-39 to -50),
the three-dimensional wing lift increment with the span factors of
Figure 6.1.4.1-15 and the flap-chord factors of -14, the lift-curve slope
with the chord extension of a translating flap, and the maximum-lift
increment (Figures 6.1.1.3-12 and -13).

The ``/FLAPIN/``, ``/POWR/`` (``FLP``), ``/SUPWH/`` and ``WING`` blocks
are mirrored as flat 1-based arrays with the source's EQUIVALENCE
offsets, because the double-slotted path indexes one word below some of
its arrays.

The tables were extracted from the source DATA statements by parsing, with
``tools/fortran_data.py``.

Reference: datcom-legacy/datcom_2000/liftfp.f
"""

import math
from typing import Dict, List, Mapping

import numpy as np

from pydatcom.utils.constants import RAD, UNUSED
from pydatcom.utils.legacy_interp import interx
from pydatcom.utils.legacy_tables import tlinex_flat

_X2128A = np.array([
    0.0, 0.02, 0.04, 0.06, 0.08, 0.1, 0.12, 0.14, 0.16, 0.18,
    0.2,
])
_X1128A = np.array([
    6.0, 7.0, 8.0,
])
_Y1128A = np.array([
    0.9, 0.878, 0.858, 0.836, 0.815, 0.794, 0.772, 0.75, 0.728, 0.708,
    0.685, 0.95, 0.938, 0.924, 0.907, 0.894, 0.878, 0.86, 0.842, 0.822,
    0.802, 0.78, 0.966, 0.957, 0.947, 0.936, 0.924, 0.91, 0.896, 0.88,
    0.862, 0.842, 0.822,
])
_X2125A = np.array([
    0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.4, 0.5,
])
_X1125A = np.array([
    0.0, 0.02, 0.04, 0.06, 0.08, 0.1, 0.12, 0.15,
])
_Y1125A = np.array([
    1.77, 2.5, 3.0, 3.46, 3.82, 4.16, 4.69, 5.14, 1.77, 2.515,
    3.03, 3.5, 3.873, 4.22, 4.78, 5.24, 1.77, 2.53, 3.06, 3.54,
    3.926, 4.29, 4.87, 5.35, 1.77, 2.545, 3.09, 3.58, 3.979, 4.35,
    4.95, 5.46, 1.77, 2.56, 3.12, 3.62, 4.032, 4.4, 5.04, 5.56,
    1.77, 2.575, 3.15, 3.66, 4.085, 4.48, 5.12, 5.69, 1.77, 2.59,
    3.18, 3.7, 4.138, 4.55, 5.21, 5.79, 1.77, 2.6, 3.22, 3.74,
    4.19, 4.62, 5.33, 5.96,
])
_X2125B = np.array([
    0.05, 0.1, 0.15, 0.2, 0.25, 0.5,
])
_X1125B = np.array([
    0.7, 0.72, 0.74, 0.76, 0.78, 0.8, 0.82, 0.84, 0.86, 0.88,
    0.9, 0.92, 0.94, 0.96, 0.98, 1.0,
])
_Y1125B = np.array([
    0.356, 0.382, 0.409, 0.431, 0.452, 0.548, 0.399, 0.426, 0.452, 0.477,
    0.498, 0.583, 0.442, 0.471, 0.499, 0.523, 0.543, 0.619, 0.485, 0.521,
    0.548, 0.569, 0.589, 0.659, 0.53, 0.569, 0.594, 0.613, 0.63, 0.693,
    0.578, 0.614, 0.639, 0.657, 0.671, 0.729, 0.619, 0.655, 0.678, 0.692,
    0.709, 0.761, 0.659, 0.696, 0.713, 0.733, 0.746, 0.793, 0.7, 0.734,
    0.75, 0.765, 0.778, 0.819, 0.742, 0.771, 0.789, 0.8, 0.81, 0.85,
    0.784, 0.809, 0.824, 0.838, 0.843, 0.875, 0.826, 0.843, 0.86, 0.865,
    0.873, 0.9, 0.865, 0.885, 0.895, 0.9, 0.903, 0.921, 0.91, 0.921,
    0.928, 0.931, 0.933, 0.938, 0.951, 0.962, 0.964, 0.966, 0.967, 0.968,
    1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
])
_X21126 = np.array([
    0.0, 10.0, 12.0, 14.0, 16.0, 18.0, 20.0, 23.0, 27.0, 30.0,
    35.0, 40.0, 50.0, 60.0,
])
_X11126 = np.array([
    0.1, 0.15, 0.2, 0.25, 0.3, 0.4, 0.5,
])
_Y11126 = np.array([
    1.0, 1.0, 0.994, 0.989, 0.97, 0.938, 0.9, 0.829, 0.755, 0.722,
    0.672, 0.641, 0.596, 0.562, 1.0, 1.0, 0.994, 0.989, 0.97, 0.937,
    0.89, 0.809, 0.737, 0.698, 0.65, 0.618, 0.569, 0.531, 1.0, 1.0,
    0.994, 0.989, 0.968, 0.936, 0.87, 0.783, 0.71, 0.673, 0.63, 0.595,
    0.542, 0.5, 1.0, 1.0, 0.994, 0.989, 0.965, 0.935, 0.85, 0.74,
    0.677, 0.644, 0.6, 0.569, 0.518, 0.48, 1.0, 1.0, 0.994, 0.989,
    0.963, 0.905, 0.8, 0.7, 0.643, 0.61, 0.57, 0.541, 0.496, 0.461,
    1.0, 1.0, 0.993, 0.969, 0.924, 0.86, 0.75, 0.656, 0.606, 0.579,
    0.54, 0.513, 0.471, 0.44, 1.0, 1.0, 0.981, 0.943, 0.88, 0.79,
    0.695, 0.625, 0.571, 0.542, 0.512, 0.49, 0.45, 0.423,
])
_X21419 = np.array([
    0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.85,
    0.9, 1.0,
])
_X11419 = np.array([
    0.0, 0.5, 1.0,
])
_Y61419 = np.array([
    0.0, 0.16, 0.305, 0.44, 0.56, 0.67, 0.772, 0.86, 0.93, 0.96,
    0.981, 1.0, 0.0, 0.14, 0.27, 0.4, 0.515, 0.63, 0.735, 0.83,
    0.912, 0.947, 0.972, 1.0, 0.0, 0.125, 0.255, 0.37, 0.49, 0.6,
    0.705, 0.8, 0.885, 0.921, 0.955, 1.0,
])
_Y11127 = np.array([
    -0.373, -0.373, -0.368, -0.345, -0.32, -0.296, -0.268, -0.23, -0.19, -0.165,
    -0.151, -0.142, -0.14, -0.45, -0.446, -0.43, -0.412, -0.39, -0.364, -0.328,
    -0.29, -0.24, -0.207, -0.188, -0.172, -0.167, -0.512, -0.51, -0.491, -0.48,
    -0.458, -0.43, -0.398, -0.356, -0.3, -0.247, -0.229, -0.2, -0.19, -0.55,
    -0.55, -0.53, -0.52, -0.5, -0.472, -0.44, -0.396, -0.34, -0.28, -0.25,
    -0.224, -0.21, -0.6, -0.6, -0.585, -0.57, -0.55, -0.521, -0.495, -0.455,
    -0.4, -0.33, -0.29, -0.25, -0.23,
])
_X21127 = np.array([
    0.0, 10.0, 20.0, 25.0, 30.0, 35.0, 40.0, 45.0, 50.0, 55.0,
    60.0, 70.0, 80.0,
])
_X11127 = np.array([
    0.15, 0.2, 0.25, 0.3, 0.4,
])
_X21418 = np.array([
    0.0, 0.25, 0.5, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 8.0,
])
_X11418 = np.array([
    -1.0, -0.9, -0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1,
])
_Y61418 = np.array([
    1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
    1.11, 1.075, 1.06, 1.049, 1.03, 1.019, 1.011, 1.01, 1.009, 1.007,
    1.2, 1.15, 1.125, 1.095, 1.06, 1.04, 1.03, 1.02, 1.019, 1.015,
    1.365, 1.245, 1.2, 1.145, 1.095, 1.065, 1.05, 1.04, 1.032, 1.025,
    1.5, 1.35, 1.285, 1.205, 1.132, 1.095, 1.072, 1.06, 1.05, 1.039,
    1.7, 1.5, 1.4, 1.29, 1.185, 1.135, 1.1, 1.085, 1.071, 1.055,
    2.0, 1.63, 1.53, 1.39, 1.25, 1.185, 1.145, 1.118, 1.1, 1.079,
    2.4, 1.9, 1.7, 1.52, 1.35, 1.26, 1.2, 1.165, 1.142, 1.11,
    2.73, 2.43, 2.1, 1.73, 1.49, 1.375, 1.3, 1.24, 1.2, 1.16,
    3.08, 2.85, 2.55, 2.18, 1.79, 1.58, 1.46, 1.385, 1.325, 1.245,
])
_X137AD = np.array([
    0.0, 2.0, 4.0, 5.0, 6.0, 8.0, 9.0, 10.0, 11.0, 12.0,
    14.0, 15.0, 16.0, 17.0, 18.0, 19.0, 20.0,
])
_Y137AD = np.array([
    1.0, 1.0, 0.979, 0.95, 0.92, 0.82, 0.8, 0.82, 0.85, 0.91,
    1.09, 1.19, 1.31, 1.43, 1.51, 1.57, 1.6,
])
_X137B2 = np.array([
    0.0, 30.0,
])
_Y137B2 = np.array([
    0.0, 1.2,
])
_X137B1 = np.array([
    0.0, 2.0, 4.0, 6.0, 8.0, 10.0, 12.0, 14.0, 16.0, 20.0,
    24.0, 28.0, 30.0,
])
_Y137B1 = np.array([
    0.0, 0.2, 0.34, 0.47, 0.57, 0.65, 0.72, 0.78, 0.83, 0.92,
    0.99, 1.04, 1.06,
])
_X137AB = np.array([
    0.0, 2.0, 5.0, 7.0, 9.0, 11.0, 13.0, 15.0, 16.0, 17.0,
    18.0, 19.0,
])
_Y137AB = np.array([
    1.0, 1.0, 1.04, 1.09, 1.17, 1.29, 1.45, 1.64, 1.73, 1.77,
    1.8, 1.82,
])
_X137AC = np.array([
    0.0, 5.0, 6.0, 8.0, 10.0, 12.0, 14.0, 16.0, 17.0, 18.0,
    19.0, 20.0,
])
_Y137AC = np.array([
    1.0, 1.0, 1.02, 1.08, 1.17, 1.3, 1.47, 1.67, 1.71, 1.73,
    1.715, 1.68,
])
_X1138A = np.array([
    0.0, 5.0, 10.0, 15.0, 20.0, 25.0, 30.0, 35.0, 40.0, 45.0,
    50.0, 60.0,
])
_Y138A1 = np.array([
    0.4, 0.5, 0.61, 0.71, 0.79, 0.87, 0.94, 0.98, 1.0, 1.0,
    1.0, 1.0,
])
_Y138A2 = np.array([
    0.18, 0.33, 0.47, 0.59, 0.7, 0.79, 0.87, 0.93, 0.97, 1.0,
    1.0, 1.0,
])
_Y138A3 = np.array([
    0.18, 0.32, 0.44, 0.56, 0.66, 0.76, 0.84, 0.9, 0.95, 0.99,
    1.0, 1.0,
])
_Y138A4 = np.array([
    0.0, 0.17, 0.33, 0.46, 0.57, 0.67, 0.76, 0.83, 0.87, 0.92,
    0.95, 1.0,
])
_X1138B = np.array([
    0.0, 0.2, 0.3, 0.4, 0.45, 0.55, 0.6, 0.8, 1.0,
])
_Y138B1 = np.array([
    0.0, 0.26, 0.39, 0.5, 0.57, 0.66, 0.7, 0.87, 1.0,
])
_Y138B2 = np.array([
    0.0, 0.11, 0.23, 0.4, 0.52, 0.66, 0.7, 0.87, 1.0,
])
_X1418A = np.array([
    0.0, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8,
    0.9, 1.0,
])
_Y1418A = np.array([
    0.0, -0.26, -0.4, -0.55, -0.66, -0.75, -0.82, -0.88, -0.93, -0.97,
    -0.99, -1.0,
])
_X61142 = np.array([
    0.0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45,
    0.5,
])
_Y61142 = np.array([
    0.0, 0.03, 0.042, 0.052, 0.06, 0.068, 0.072, 0.078, 0.082, 0.086,
    0.09,
])
_X6143A = np.array([
    15.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 0.1, 0.25, 0.4,
    0.0, 0.0, 0.0, 0.0,
])
_Y6143A = np.array([
    0.66, 0.66, 0.645, 0.505, 0.5, 0.4, 0.35, 0.73, 0.73, 0.7,
    0.645, 0.55, 0.43, 0.38, 0.78, 0.78, 0.745, 0.68, 0.585, 0.475,
    0.4,
])
_X6143B = np.array([
    0.0, 15.0, 20.0, 30.0, 40.0, 60.0, 90.0, 0.0, 10.0, 20.0,
    0.0, 0.0, 0.0, 0.0,
])
_Y6143B = np.array([
    1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.945,
    0.89, 0.885, 0.875, 0.865, 1.0, 1.0, 0.9, 0.79, 0.75, 0.725,
    0.72,
])
_X11147 = np.array([
    5.0, 10.0, 15.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 75.0,
    0.1, 0.2, 0.3, 0.4, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
])
_Y11147 = np.array([
    -0.36, -0.295, -0.24, -0.22, -0.195, -0.18, -0.17, -0.155, -0.14, -0.13,
    -0.435, -0.37, -0.32, -0.29, -0.26, -0.245, -0.23, -0.215, -0.195, -0.185,
    -0.57, -0.48, -0.415, -0.375, -0.33, -0.305, -0.28, -0.27, -0.255, -0.25,
    -0.6, -0.545, -0.48, -0.44, -0.378, -0.35, -0.325, -0.3, -0.275, -0.26,
])
_X11150 = np.array([
    0.0, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5,
])
_Y11150 = np.array([
    0.0, -0.0005, -0.0016, -0.0044, -0.0082, -0.0108, -0.02,
])


class _View:
    """A source array EQUIVALENCEd into a COMMON block: index ``k`` is
    word ``offset + k`` of the flat 1-based block, so an index of 0 reads
    and writes the word before, as the compiled program does."""

    def __init__(self, block: List[float], offset: int):
        self.block, self.offset = block, offset

    def __getitem__(self, k):
        return self.block[self.offset + k]

    def __setitem__(self, k, value):
        self.block[self.offset + k] = value


def _ix1(x, y, q, l1=0, u1=0):
    return float(interx(1, x, [q], [len(x)], y, lind=len(x), lx1l=l1,
                        lx1u=u1))


def _ix2(table, dep, q1, q2, n1, n2, lind):
    return float(interx(2, table, [q1, q2], [n1, n2], dep, lind=lind,
                        lx1l=1, lx2l=1, lx1u=1, lx2u=1))


_TRANSLATING = (2, 3, 4, 7, 8)


def calculate_liftfp(data: Mapping[str, object]) -> Dict[str, object]:
    """Translate LIFTFP: flap lift, lift-curve slope and maximum lift.

    Args:
        data: By name: ``htpl``, ``transn``, ``m`` (the Mach index),
            ``mach``, ``rl`` (``FLC(M+42)``), ``sref``, ``surface`` (the
            wing's or tail's words: ``cbarex`` (``A(16)``), ``bo2``
            (``WINGIN(3)``), ``tante``, ``tanle``, ``cr``, ``tapri``,
            ``bstro2``, ``tanteo``, ``tanleo``, ``cb``, ``tapro``, ``aw``,
            ``btheo``, ``tovc``, ``tovco``, ``cosc4``, ``clasec``
            (``WINGIN(M+20)``, or ``WINGIN(69)/.8`` transonic), ``cla``
            (``WING(101)`` or ``TRA(70)``)), ``f`` (the 116 ``/FLAPIN/``
            words), ``flp`` (the 189 ``FLP`` words), ``fcm282`` (the
            ``/SUPWH/`` words 282-287), ``wing`` (the ``WING`` words
            201-250 as they stood), and ``state``: the saved local ``N``.

    Returns:
        ``f``, ``flp`` (the blocks as left, 1-based lists with a leading
        pad), ``fcm282``, ``wing`` (201-250), ``state``.

    Notes:
        Kept as executed: the double-slotted flap enters the strip
        computation ahead of the strip index update, so it uses and fills
        the previous strip's words (at the first strip, the fourth strip's,
        since the geometry loop leaves the index at 4); once any deflection gives
        its section lift (``SDCL``) every later one takes that path too;
        and with experimental data on a translating flap the chord factors
        are formed on a stale strip index.
    """
    s = {k: float(v) for k, v in data['surface'].items()}
    f = [0.0] + [float(v) for v in data['f']]
    flp = [0.0] + [float(v) for v in data['flp']]
    fcm = {282 + word_offset: float(v)
           for word_offset, v in enumerate(data['fcm282'])}
    wing = {200 + word: float(v)
            for word, v in enumerate(data['wing'], 1)}
    delta = _View(f, 0)
    sdcl, cpi, cpo = _View(f, 18), _View(f, 38), _View(f, 48)
    capi, capo, df2 = _View(f, 84), _View(f, 94), _View(f, 104)
    eta, chrd, cf = _View(flp, 0), _View(flp, 5), _View(flp, 10)
    aldavg, dkb, swf = _View(flp, 15), _View(flp, 19), _View(flp, 23)
    cp, cldoct, cldthy = _View(flp, 27), _View(flp, 33), _View(flp, 38)
    cfoc, adcads, cfact = _View(flp, 61), _View(flp, 66), _View(flp, 70)
    dsclmx, rk2, alphad = _View(flp, 80), _View(flp, 90), _View(flp, 104)
    delcla, aldag = _View(flp, 109), _View(flp, 149)

    class _FcmView:
        def __getitem__(self, k):
            return fcm[282 + k]

        def __setitem__(self, k, v):
            fcm[282 + k] = v
    delcl = _FcmView()
    cpocf = [0.0] * 5
    # CFACTR is a saved local: with experimental data on a translating
    # flap only one strip's entry is written, the rest keep the last call's.
    cfactr = [0.0] + [float(v) for v in
                      data.get('state', {}).get('cfactr', [0.0] * 4)]
    dclk = [0.0] * 5
    rkb = [0.0] * 6
    cf2 = [0.0] * 6
    cf2oc = [0.0] * 6
    tanphe, cfi, cfo = f[11], f[12], f[13]
    bif, bof, cf2i, cf2o = f[14], f[15], f[115], f[116]
    transn = bool(data['transn'])
    clasec, claw = s['clasec'], s['cla']
    iftype = int(f[17] + 0.5)
    ndelta = int(f[16] + 0.5)
    rf = float(data['rl'])
    mach = 0.6 if transn else float(data['mach'])
    sr, aw, btheo = float(data['sref']), s['aw'], s['btheo']
    for deflection_index in range(1, ndelta + 1):
        if delta[deflection_index] == 0.0:
            delta[deflection_index] = 0.01
    transl = iftype in _TRANSLATING
    deln4 = 0.25 * (bof - bif) / btheo
    flp[60] = deln4
    eta[1] = bif / btheo
    cf[1] = cfi
    arg1 = (cfi - cfo) / (4. * deln4)
    arg2 = (s['tante'] - s['tanle']) * btheo
    arg3 = (bof - bif) * s['cr'] / 4.0
    arg4, arg5, tc = s['cr'], s['tapri'], s['tovc']
    if not bif < btheo - s['bstro2']:
        arg2 = (s['tanteo'] - s['tanleo']) * s['bo2']
        arg3 = (bof - bif) * s['cb'] / 4.0
        arg4, arg5, tc = s['cb'], s['tapro'], s['tovco']
    rkb[1] = tlinex_flat(_X11419, _X21419, _Y61419, arg5, eta[1])
    chrd[1] = arg4 + eta[1] * arg2
    cfoc[1] = cf[1] / chrd[1]
    alphad[1] = _ix1(_X1418A, _Y1418A, cfoc[1])
    for strip_index in range(2, 6):
        inboard = strip_index - 1
        eta[strip_index] = eta[inboard] + deln4
        cf[strip_index] = cfi - arg1 * (eta[strip_index] - eta[1])
        rkb[strip_index] = tlinex_flat(_X11419, _X21419, _Y61419, arg5,
                                       eta[strip_index])
        dkb[inboard] = rkb[strip_index] - rkb[inboard]
        chrd[strip_index] = arg4 + eta[strip_index] * arg2
        cfoc[strip_index] = cf[strip_index] / chrd[strip_index]
        alphad[strip_index] = _ix1(_X1418A, _Y1418A, cfoc[strip_index])
        aldavg[inboard] = 0.50 * (alphad[strip_index] + alphad[inboard])
        swf[inboard] = (arg3 * (2. - (1. - arg5) *
                                (eta[inboard] + eta[strip_index])))
    n = 4                      # the geometry loop's N (source strip count)
    expdcl = False
    nn = 0
    arg1 = math.log10(rf * s['cbarex'])
    if iftype == 1:
        flp[33] = tlinex_flat(_X1128A, _X2128A, _Y1128A, arg1, tanphe,
                              1, 0, 0, 1)
    adcad = {}
    for deflection_index in range(1, ndelta + 1):
        arg8 = delta[deflection_index] * clasec
        if sdcl[deflection_index] != UNUSED:
            expdcl = True
        if not (expdcl and not transl):
            cp[1] = cpi[deflection_index]
            for strip_index in range(1, 6):
                if strip_index != 1 and transl:
                    arg1 = (cpi[deflection_index] - cpo[deflection_index]) / (
                        4. * deln4)
                    cp[strip_index] = (cpi[deflection_index] -
                                       arg1 * (eta[strip_index] - eta[1]))
                label = 1150 if expdcl else None
                if label is None:
                    argz = abs(delta[deflection_index])
                    kind = iftype if 1 <= iftype <= 8 else 1
                    if kind == 1:
                        if deflection_index <= 1:
                            cldoct[strip_index] = tlinex_flat(
                                 _X1125B, _X2125B, _Y1125B, flp[33],
                                cfoc[strip_index])
                            cldthy[strip_index] = tlinex_flat(
                                 _X1125A, _X2125A, _Y1125A, tc,
                                cfoc[strip_index])
                        kfprm = tlinex_flat(_X11126, _X21126, _Y11126,
                                            cfoc[strip_index], argz)
                        delcl[strip_index] = (delta[deflection_index] *
                                              cldoct[strip_index] *
                                              cldthy[strip_index] * kfprm /
                                              RAD)
                        label = 1140
                    elif kind in (2, 3):
                        alphad[strip_index] = tlinex_flat(
                            _X11127, _X21127, _Y11127, cfoc[strip_index],
                            argz)
                        delcl[strip_index] = (-clasec * alphad[strip_index] *
                                              delta[deflection_index])
                        label = 1140
                    elif kind == 4:
                        label = 1150
                    elif kind == 5:
                        alfad = _ix2(_X11147, _Y11147, delta[deflection_index],
                                     cfoc[strip_index], 10, 4, 10)
                        delcl[strip_index] = (-clasec * alfad *
                                              delta[deflection_index])
                        label = 1140
                    else:
                        cldk = _ix1(_X11150, _Y11150, cfoc[strip_index], 0, 1)
                        if kind == 6:
                            delcl[strip_index] = cldk * delta[deflection_index]
                        else:
                            delcl[strip_index] = (cldk * delta[deflection_index] *
                                                  cp[strip_index] /
                                                  chrd[strip_index])
                        label = 1140
                while label is not None:
                    if label == 1140:
                        n = strip_index - 1
                        if n == 0:
                            label = None
                            continue
                        nn += 1
                        delcla[nn] = (delcl[strip_index] + delcl[n]) / 2.
                        aldag[nn] = -delcla[nn] / arg8
                        if iftype == 4 or not transl:
                            label = 1160
                        else:
                            label = 1150
                    elif label == 1150:
                        cpocf[n] = (cp[strip_index] / chrd[strip_index] +
                                    cp[n] / chrd[n]) / 2.
                        cfactr[n] = (cpocf[n] - 1.) * swf[n] / sr
                        if expdcl:
                            label = None
                            continue
                        if iftype != 4:
                            label = 1160
                            continue
                        aarg1 = (cf2i - cf2o) / (4. * deln4)
                        aarg2 = (capi[deflection_index] -
                                 capo[deflection_index]) / (4. * deln4)
                        capr = (capi[deflection_index] -
                                aarg2 * (eta[strip_index] - eta[1]))
                        cf2[strip_index] = (cf2i - aarg1 *
                                            (eta[strip_index] - eta[1]))
                        cf2oc[strip_index] = (cf2[strip_index] /
                                              chrd[strip_index])
                        phi1 = delta[deflection_index] + math.atan(tanphe) * RAD
                        phi2 = phi1 + df2[deflection_index]

                        def clamp(v):
                            return min(max(v, .10), .40)
                        atea1 = _ix2(_X6143A, _Y6143A, phi1,
                                     clamp(cfoc[strip_index]), 7, 3, 7)
                        cldf1 = _ix1(_X61142, _Y61142, cfoc[strip_index], 1, 1)
                        atea2 = _ix2(_X6143A, _Y6143A, phi2,
                                     clamp(cf2oc[strip_index]), 7, 3, 7)
                        cldf2 = _ix1(_X61142, _Y61142, cfoc[strip_index], 1, 1)
                        if cfoc[strip_index] / cf2oc[strip_index] <= 0.60:
                            delcl[n] = (atea1 * cldf1 * delta[deflection_index] *
                                        (1. + cfoc[strip_index]) + atea2 *
                                        cldf2 * (delta[deflection_index] +
                                                 df2[deflection_index]) *
                                        cp[strip_index] / chrd[strip_index])
                        else:
                            ateat = _ix2(_X6143B, _Y6143B,
                                           df2[deflection_index],
                                           delta[deflection_index], 7, 3, 7)
                            delcl[n] = (atea1 * cldf1 * delta[deflection_index] *
                                        capr / chrd[strip_index] + atea2 *
                                        ateat * cldf2 * df2[deflection_index] *
                                        (1. + (cp[strip_index] - capr) /
                                         chrd[strip_index]))
                        label = 1140
                    elif label == 1160:
                        adcads[n] = tlinex_flat(_X11418, _X21418, _Y61418,
                                                aldavg[n], aw, 0, 0, 0, 1)
                        dclk[n] = delcla[nn] * adcads[n] * dkb[n] * claw / \
                            clasec
                        label = None
            if not expdcl:
                wing[200 + deflection_index] = (dclk[1] + dclk[2] + dclk[3] +
                                                dclk[4])
            if expdcl or transl:
                cfact[deflection_index] = (cfactr[1] + cfactr[2] +
                                           cfactr[3] + cfactr[4]) / 4.0
                wing[240 + deflection_index] = (cfact[deflection_index] * claw +
                                                claw)
            if not expdcl:
                continue
        adcad[deflection_index] = tlinex_flat(
             _X11418, _X21418, _Y61418,
            sdcl[deflection_index] / (clasec * delta[deflection_index]), aw, 0,
            0, 0, 1)
        wing[200 + deflection_index] = (sdcl[deflection_index] *
                                        adcad[deflection_index] *
                                        (rkb[5] - rkb[1]) * claw / clasec)
    if iftype < 6:
        for deflection_index in range(1, ndelta + 1):
            if deflection_index == 1:
                trofs = chrd[5] / chrd[1]
                cbarfs = 2. / 3. * chrd[1] * (1. + trofs * (1. + trofs)) / \
                    (1. + trofs)
                etafs = (1. + 2. * trofs) / (3. * (1. + trofs))
                flp[61] = (cfi - etafs * (cfi - cfo)) / cbarfs
                tcp = tc * 100.
                if iftype == 2 or iftype < 1:
                    flp[102] = _ix1(_X137AC, _Y137AC, tcp)
                elif iftype in (3, 4):
                    flp[102] = _ix1(_X137AB, _Y137AB, tcp)
                else:
                    flp[102] = _ix1(_X137AD, _Y137AD, tcp)
                if iftype in (3, 4):
                    flp[101] = _ix1(_X137B2, _Y137B2, flp[61] * 100., 0, 1)
                else:
                    flp[101] = _ix1(_X137B1, _Y137B1, flp[61] * 100., 0, 1)
            ad = abs(delta[deflection_index])
            table = {2: _Y138A2, 3: _Y138A1, 4: _Y138A3}.get(iftype,
                                                             _Y138A4)
            rk2[deflection_index] = _ix1(_X1138A, table, ad)
            if iftype in (1, 5):
                flp[103] = 1.
            elif iftype == 4:
                flp[103] = _ix1(_X1138B, _Y138B2, ad / 50.)
            else:
                ref = 45. if iftype == 2 else 40.
                flp[103] = _ix1(_X1138B, _Y138B1, ad / ref)
            dsclmx[deflection_index] = (flp[101] * rk2[deflection_index] *
                                        flp[103] * flp[102])
            swft = swf[1] + swf[2] + swf[3] + swf[4]
            flp[104] = (1. - 0.08 * s['cosc4']**2) * s['cosc4']**0.75
            wing[220 + deflection_index] = (dsclmx[deflection_index] * swft *
                                            flp[104] / sr)
    return {'f': f, 'flp': flp,
            'fcm282': [fcm[282 + word] for word in range(6)],
            'wing': [wing[200 + word] for word in range(1, 51)],
            'state': {'cfactr': cfactr[1:5]}, 'method': 'legacy_liftfp'}
