"""
Propeller power effects: PRPWEF and its overlay M13O15 (DATCOM Section 4.6).

The routine works through the propeller's normal force and the slipstream:
the normal-force derivative (Figures 4.6.1-25A/B) and the downwash behind
the disk (Figure 4.6.1-26), the wing area immersed in the slipstream and
its lift through the slipstream factor (Figure 4.6.1-27), the upwash ahead
of the wing (Figure 4.4.1-61), the tail's slipstream downwash and dynamic
pressure (Figures 4.6.3-14 to -16), the thrust and normal-force moments,
the immersed-area zero-lift moment, and the drag of the immersed surfaces
with the change in drag due to lift (Figures 4.6.4-13A/B).  M13O15 then
forms CN, CA and the slopes of the increments.

The tables were extracted from the source DATA statements by parsing, with
``tools/fortran_data.py``, rather than by hand.

Reference: datcom-legacy/datcom_2000/prpwef.f, m13o15.f
"""

import math
from typing import Dict, Mapping, Optional

import numpy as np

from pydatcom.utils.constants import DEG, PI, RAD, UNUSED
from pydatcom.utils.legacy_numeric import angles, tbfunx, zerang
from pydatcom.utils.legacy_tables import tlin3x, tlinex

# Figure 4.6.1-25B: the propeller inflow factor F.
_X6111A = np.array([
    0.0, 1.0, 2.0, 3.0, 4.0, 6.0, 8.0, 10.0,
    14.0, 19.0, 22.0,
])
_F6111A = np.array([
    1.0, 1.55, 1.94, 2.2, 2.4, 2.75, 3.05, 3.3,
    3.75, 4.25, 4.54,
])
# Figure 4.6.1-25A: the normal-force derivative at K_N = 80.7, by blade count and blade angle (and the dual-rotation curve).
_X1S11B = np.array([
    2.0, 3.0, 4.0, 6.0,
])
_X2S11B = np.array([
    15.0, 20.0, 25.0, 30.0, 35.0, 40.0, 60.0,
])
_YSR11B = np.array([
    0.08, 0.1, 0.115, 0.126, 0.14, 0.15, 0.192, 0.11,
    0.139, 0.16, 0.182, 0.2, 0.216, 0.275, 0.136, 0.172,
    0.2, 0.226, 0.25, 0.272, 0.35, 0.196, 0.237, 0.275,
    0.315, 0.35, 0.382, 0.5,
])
_X6111B = np.array([
    15.0, 20.0, 30.0, 40.0, 50.0, 60.0,
])
_Y6111B = np.array([
    0.25, 0.295, 0.366, 0.421, 0.468, 0.51,
])
# Figure 4.6.1-27: the slipstream lift factor (four chained parts).
_X26112 = np.array([
    2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0,
    10.0,
])
_X16112 = np.array([
    1.0, 1.5, 2.0,
])
_Y46112 = np.array([
    3.8, 2.3, 1.5, 0.9, 0.65, 0.35, 0.2, 0.1,
    0.0, 9.0, 6.6, 5.4, 4.5, 3.9, 3.6, 3.4,
    3.1, 3.0, 13.0, 10.3, 9.0, 8.05, 7.45, 6.9,
    6.55, 6.25, 6.0,
])
_XB6112 = np.array([
    0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5,
    4.0,
])
_DM2 = np.array([
    0.0, 1.3, 1.95, 2.4, 2.75, 3.1, 3.35, 3.56,
    3.75,
])
_DC2 = np.array([
    0.0, 14.0,
])
_DC1 = np.array([
    1.0, 2.0, 3.0, 4.0, 5.0,
])
_DC3 = np.array([
    0.3, 1.0, 0.6, 2.0, 0.9, 3.0, 1.2, 4.0,
    1.5, 5.0,
])
_D4 = np.array([
    0.0, 5.0,
])
_D3 = np.array([
    100.0, 4.0, 3.0, 2.0, 1.5, 1.0, 0.8, 0.6,
    0.4, 0.2, 0.0,
])
_AK6112 = np.array([
    0.0, 1.0, 0.21, 1.21, 0.32, 1.32, 0.46, 1.46,
    0.54, 1.54, 0.61, 1.61, 0.68, 1.68, 0.72, 1.72,
    0.8, 1.8, 0.89, 1.89, 1.0, 2.0,
])
# Figure 4.4.1-61: the upwash gradient ahead of the wing.
_X14161 = np.array([
    4.0, 6.0, 9.0, 12.0,
])
_X24161 = np.array([
    0.25, 0.4, 0.5, 0.6, 0.7, 0.8, 1.0, 1.2,
    1.6, 2.0,
])
_Y44161 = np.array([
    1.08, 0.545, 0.4, 0.31, 0.24, 0.2, 0.13, 0.1,
    0.06, 0.04, 1.18, 0.68, 0.52, 0.4, 0.32, 0.27,
    0.19, 0.15, 0.1, 0.08, 1.3, 0.81, 0.62, 0.49,
    0.4, 0.34, 0.25, 0.2, 0.13, 0.12, 1.4, 0.88,
    0.7, 0.56, 0.445, 0.39, 0.3, 0.24, 0.165, 0.14,
])
# Figure 4.6.1-26: the downwash gradient behind the propeller (C1, C2).
_X46113 = np.array([
    0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5,
    4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0,
    12.0, 19.0, 20.0,
])
_C16113 = np.array([
    0.0, 0.25, 0.36, 0.44, 0.49, 0.54, 0.575, 0.61,
    0.64, 0.675, 0.708, 0.73, 0.75, 0.76, 0.78, 0.79,
    0.798, 0.84,
])
_C26113 = np.array([
    0.25, 0.24, 0.225, 0.205, 0.19, 0.175, 0.165, 0.16,
    0.15, 0.14, 0.13, 0.12, 0.11, 0.1, 0.099, 0.098,
    0.097, 0.078, 0.075,
])
# Figure 4.6.3-16: the dynamic-pressure change at the tail.
_X14637 = np.array([
    0.0, 0.2, 0.4, 0.6, 0.8, 1.0,
])
_X24637 = np.array([
    0.0, 0.2, 0.4, 0.8, 1.2, 1.8,
])
_Y4637 = np.array([
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.95,
    1.25, 2.0, 2.75, 3.8, 0.0, 1.99, 2.8, 4.1,
    5.5, 7.7, 0.0, 2.8, 4.25, 6.3, 8.3, 11.3,
    0.0, 3.7, 5.75, 8.5, 11.05, 15.0, 0.0, 4.61,
    6.8, 10.25, 13.75, 18.8,
])
_Z24637 = np.array([
    0.0, 20.0,
])
_XT4637 = np.array([
    0.0, 0.2, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9,
    1.0, 1.1, 1.2, 1.3, 1.4,
])
_YF4637 = np.array([
    -0.1, 1.9, -0.1, 1.8, -0.1, 1.72, -0.1, 1.63,
    -0.1, 1.46, -0.1, 1.31, -0.1, 1.175, -0.1, 1.025,
    -0.1, 0.85, -0.1, 0.63, -0.1, 0.43, -0.1, 0.22,
    -0.1, 0.03,
])
# Figures 4.6.4-13A/B: the drag factors.
_X4648A = np.array([
    0.0, 5.0, 10.0, 15.0, 20.0, 25.0, 30.0, 35.0,
    40.0, 45.0, 48.5,
])
_Y4648A = np.array([
    4.0, 2.84, 2.2, 1.75, 1.45, 1.25, 1.08, 0.98,
    0.87, 0.82, 0.76,
])
_X1648B = np.array([
    0.1, 0.2, 0.3, 0.4, 0.5,
])
_X2648B = np.array([
    0.0, 2.5, 5.0, 7.5, 10.0, 12.5, 15.0, 20.0,
    25.0, 30.0, 35.0, 40.0, 45.0, 50.0,
])
_Y4648B = np.array([
    0.0, 0.22, 0.32, 0.39, 0.44, 0.49, 0.53, 0.59,
    0.64, 0.68, 0.71, 0.73, 0.75, 0.76, 0.3, 0.45,
    0.542, 0.6, 0.655, 0.7, 0.73, 0.776, 0.81, 0.835,
    0.85, 0.87, 0.88, 0.89, 0.55, 0.66, 0.73, 0.78,
    0.81, 0.83, 0.85, 0.88, 0.9, 0.915, 0.93, 0.94,
    0.942, 0.95, 0.81, 0.84, 0.87, 0.89, 0.905, 0.92,
    0.93, 0.95, 0.96, 0.966, 0.971, 0.975, 0.98, 0.98,
    1.0, 0.99, 0.978, 0.97, 0.967, 0.968, 0.972, 0.982,
    0.989, 0.995, 1.0, 1.0, 1.0, 1.0,
])
# Figure 4.1.4.1-5 (the 16-sweep version): twist effect on CM0.
_X31412 = np.array([
    0.0, 0.5, 1.0,
])
_X11412 = np.array([
    10.0, 8.0, 6.0, 3.5, 1.5,
])
_X21412 = np.array([
    0.0, 2.5, 5.0, 7.5, 10.0, 12.5, 15.0, 20.0,
    25.0, 30.0, 35.0, 40.0, 45.0, 50.0, 55.0, 60.0,
])
_Y41412 = np.array([
    -0.0, -0.0005, -0.0011, -0.0016, -0.0022, -0.0028, -0.0033, -0.0045,
    -0.0057, -0.0069, -0.0081, -0.0094, -0.0107, -0.012, -0.0133, -0.0147,
    -0.0, -0.0004, -0.0008, -0.0012, -0.0016, -0.002, -0.0025, -0.0033,
    -0.0041, -0.005, -0.0059, -0.0068, -0.0078, -0.0088, -0.0099, -0.0111,
    -0.0, -0.0002, -0.0005, -0.0008, -0.001, -0.0013, -0.0015, -0.0021,
    -0.0027, -0.0033, -0.0039, -0.0045, -0.0052, -0.0059, -0.0066, -0.0074,
    -0.0, -0.0001, -0.0002, -0.0003, -0.0004, -0.0005, -0.0007, -0.0009,
    -0.0011, -0.0014, -0.0017, -0.002, -0.0024, -0.0027, -0.0032, -0.0036,
    -0.0, -0.0, -0.0, -0.0001, -0.0001, -0.0001, -0.0001, -0.0002,
    -0.0003, -0.0003, -0.0004, -0.0005, -0.0006, -0.0007, -0.0008, -0.0009,
    -0.0, -0.0009, -0.0018, -0.0027, -0.0037, -0.0046, -0.0056, -0.0076,
    -0.0097, -0.0118, -0.014, -0.0162, -0.0185, -0.0208, -0.0232, -0.0257,
    -0.0, -0.0006, -0.0013, -0.0019, -0.0026, -0.0033, -0.004, -0.0054,
    -0.0069, -0.0084, -0.0099, -0.0116, -0.0134, -0.0153, -0.0172, -0.0193,
    -0.0, -0.0004, -0.0008, -0.0012, -0.0016, -0.0021, -0.0025, -0.0034,
    -0.0043, -0.0053, -0.0063, -0.0074, -0.0086, -0.0099, -0.0113, -0.0128,
    -0.0, -0.0002, -0.0003, -0.0005, -0.0007, -0.0008, -0.001, -0.0014,
    -0.0018, -0.0022, -0.0026, -0.0031, -0.0037, -0.0043, -0.005, -0.0058,
    -0.0, -0.0, -0.0001, -0.0001, -0.0001, -0.0002, -0.0002, -0.0003,
    -0.0004, -0.0005, -0.0006, -0.0007, -0.0008, -0.001, -0.0011, -0.0013,
    -0.0, -0.001, -0.0019, -0.0029, -0.0039, -0.0049, -0.0059, -0.008,
    -0.0101, -0.0122, -0.0143, -0.0166, -0.019, -0.0214, -0.024, -0.0266,
    -0.0, -0.0007, -0.0013, -0.002, -0.0027, -0.0034, -0.0041, -0.0056,
    -0.0071, -0.0086, -0.0101, -0.0118, -0.0136, -0.0155, -0.0176, -0.0197,
    -0.0, -0.0004, -0.0008, -0.0012, -0.0017, -0.0021, -0.0026, -0.0035,
    -0.0044, -0.0054, -0.0065, -0.0076, -0.0089, -0.0102, -0.0117, -0.0132,
    -0.0, -0.0002, -0.0003, -0.0005, -0.0007, -0.0008, -0.001, -0.0014,
    -0.0018, -0.0022, -0.0026, -0.0031, -0.0037, -0.0043, -0.005, -0.0057,
    -0.0, -0.0, -0.0001, -0.0002, -0.0002, -0.0003, -0.0003, -0.0004,
    -0.0005, -0.0006, -0.0007, -0.0008, -0.0009, -0.001, -0.0012, -0.0013,
])
# Figure 4.6.3-14: slipstream downwash at the tail, one engine.
_X1638A = np.array([
    14.0, 12.0, 10.0, 8.0, 6.0, 4.0, 2.0, 0.0,
])
_X2638A = np.array([
    0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4,
])
_Y4638A = np.array([
    0.0, 1.5, 2.17, 2.68, 3.2, 3.64, 4.06, 4.47,
    0.0, 1.27, 1.93, 2.4, 2.8, 3.12, 3.4, 3.68,
    0.0, 1.03, 1.55, 1.96, 2.28, 2.55, 2.79, 3.0,
    0.0, 0.8, 1.2, 1.5, 1.7, 1.9, 2.06, 2.2,
    0.0, 0.56, 0.9, 1.11, 1.3, 1.45, 1.58, 1.68,
    0.0, 0.4, 0.65, 0.8, 0.94, 1.08, 1.11, 1.18,
    0.0, 0.3, 0.45, 0.52, 0.61, 0.65, 0.7, 0.72,
    0.0, 0.15, 0.21, 0.26, 0.31, 0.32, 0.32, 0.31,
])
_X1638B = np.array([
    0.0, 0.6, 1.0, 1.2,
])
_X2638B = np.array([
    0.0, 4.8,
])
_Y4638B = np.array([
    0.0, 9.85, 0.0, 8.7, 0.0, 7.5, 0.0, 5.9,
])
# Figure 4.6.3-15: the same for several engines.
_X1639A = np.array([
    10.0, 8.0, 6.0, 5.0, 4.0, 3.0, 2.0, 1.0,
])
_X2639A = np.array([
    0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4,
])
_Y4639A = np.array([
    0.0, 1.5, 2.25, 2.77, 3.1, 3.37, 3.6, 3.8,
    0.0, 1.4, 2.1, 2.57, 2.87, 3.12, 3.31, 3.48,
    0.0, 1.17, 1.8, 2.2, 2.47, 2.68, 2.87, 3.01,
    0.0, 1.08, 1.6, 1.95, 2.2, 2.4, 2.54, 2.66,
    0.0, 0.9, 1.32, 1.6, 1.82, 2.0, 2.13, 2.25,
    0.0, 0.7, 1.05, 1.28, 1.48, 1.6, 1.7, 1.78,
    0.0, 0.5, 0.78, 0.92, 1.05, 1.13, 1.2, 1.25,
    0.0, 0.28, 0.43, 0.52, 0.6, 0.61, 0.62, 0.62,
])
_X4639B = np.array([
    0.0, 4.0,
])
_Y4639B = np.array([
    0.0, 4.0,
])
_X1639C = np.array([
    0.0, 0.4, 0.7, 1.0, 1.2,
])
_X2639C = np.array([
    0.0, 4.0,
])
_Y4639C = np.array([
    0.0, 10.0, 0.0, 9.6, 0.0, 8.75, 0.0, 7.45,
    0.0, 5.9,
])
_PISQRD = 9.8696044

# /POWR/ as PRPWEF lays it out: (name, first word, length).
POWR_LAYOUT = [
    ('dclt', 1, 20), ('xbarp', 21, 1), ('deuda', 22, 1), ('dclnp', 23, 20),
    ('dclq', 43, 20), ('dclaw', 63, 20), ('dclhq', 83, 20),
    ('dcmnp', 103, 20), ('dcmq', 123, 1), ('dcml', 124, 20),
    ('dcmhq', 144, 20), ('dcmhe', 164, 20),
] + [(name, 184 + k, 1) for k, name in enumerate(
    'sinapx prprd2 cti bstio2 sstri bst0i2 ctih sst0i sratio cnap80 cnap '
    'c1 c2 depdap srtpco f combo1 combo cosaih siosrh sih dcd0s cd0pow '
    'rpnob aak ebroep dcmt astari trpsti xbrsrr alphat alphap ep sinap zs '
    'bio2 cosat sinat si tri cbarli sweepa trpsi scapi trs0i cbsr0i cossw '
    'atovca cm0in cm02 cm0ova cm0te0 cm0i bs1 bs2 bs3 ak1 delalp dxhmac '
    'zheff zhorp dqhoqi zht zhtorp xcp dlh cnp clp ebar clww cdlrat '
    'cdlpow epowr ytemp step1'.split())] + [('dclhe', 259, 20),
                                            ('argcs', 279, 7)]

# PRPWEF assigns these scalars; the rest of /POWR/ it leaves alone.
_ASSIGNED = set(
    'xbarp deuda dcmq prprd2 cnap80 cnap c1 c2 depdap srtpco f combo1 '
    'combo dcd0s cd0pow rpnob aak ebroep dcmt cm0in cm02 cm0ova cm0te0 '
    'cm0i xbrsrr cossw'.split())


def _tlinex(x1, x2, flat, q1, q2, l1, l2, u1, u2) -> float:
    """TLINEX on a source ``Y(NX2,NX1)`` table given flat."""
    y = np.asarray(flat, dtype=float).reshape(len(x1), len(x2)).T
    return float(tlinex(x1, x2, y, q1, q2, l1, l2, u1, u2))


def _value(x, y, query, lower=0, upper=0) -> float:
    return float(tbfunx(x, y, query, lower, upper)[0])


def _sqrt(x: float) -> float:
    """SQRT as the source runs it: NaN, not an error, below zero."""
    return math.sqrt(x) if x >= 0.0 else math.nan


def _div(dividend: float, divisor: float) -> float:
    """IEEE division as the source runs it: 0/0 is NaN, x/0 infinite."""
    if divisor == 0.0:
        return math.nan if dividend == 0.0 or math.isnan(dividend) else math.copysign(
            math.inf, dividend) * math.copysign(1.0, divisor)
    return dividend / divisor


def _funcb(tri: float) -> float:
    return 2. * (1. + tri * (1. + tri)) / (3. * (1. + tri))


def calculate_prpwef(inputs: Mapping[str, object],
                     stale: Optional[Mapping[str, float]] = None
                     ) -> Dict[str, object]:
    """Translate PRPWEF: propeller power effects, Section 4.6.

    The slipstream and normal-force effects of running propellers on lift,
    drag and pitching moment: the propeller normal force (Figures
    4.6.1-25A/B), the downwash behind the disk (Figure 4.6.1-26), the
    immersed wing area and its lift through the slipstream factor (Figure
    4.6.1-27), the upwash ahead of the wing (Figure 4.4.1-61), the tail's
    slipstream downwash and dynamic pressure (Figures 4.6.3-14 to -16),
    thrust and normal-force moments, the immersed-area zero-lift moment
    (with the Figure 4.1.4.1-5 twist term), and the drag of the immersed
    surfaces and the drag-due-to-lift change (Figures 4.6.4-13A/B).

    Args:
        inputs: ``alpha`` (``FLC(23..)``, as many as ``FLC(2)``), ``sref``,
            ``cbarr``, ``htpl``;
            ``power``: ``aietlp``, ``nengsp``, ``thstcp``, ``phaloc``,
            ``phvloc``, ``prprad``, ``kn`` (``ENGFCT``, ``UNUSED`` to be
            formed from the blade widths ``bwapr3``, ``bwapr6``,
            ``bwapr9``), ``nopbpe``, ``bapr75``, ``yp``, ``crot``;
            ``wing``: ``ct``, ``bst0o2``, ``bsto2``, ``bo2``, ``cb``,
            ``cr`` (``WINGIN(1..6)``), ``twista``, ``cmo``, ``cmot``,
            ``cl`` (``WING(21..)``), ``cla`` (``WING(101)``), ``cd0``
            (``B(46)``), ``cdl`` (``D(36..)``), ``cf`` (``D(10)+D(11)``),
            ``d12`` (``D(12)``);
            ``a``: the wing's ``A`` block (1, 2, 3, 10, 16, 23, 26, 32, 33,
            62, 67, 69, 86, 91, 106, 112, 120, 121, 134, 161);
            ``position``: ``xcg``, ``xw``, ``zw``, ``aliw``, ``zcg``,
            ``xh``, ``zh``, ``alih``;
            ``tail``: ``bo2``, ``cr``, ``ct`` (``HTIN(4)``, ``(6)``,
            ``(1)``), ``area`` and ``xbarr`` (``AHT(3)``, ``(161)``),
            ``cf`` (``DHT(10)+DHT(11)``), ``cl`` (``WBT(110..)``);
            ``vertical``: ``area`` (``AVT(3)``), ``cf``
            (``DVT(10)+DVT(11)``);
            ``body``: ``cf`` (``BD(92)``), ``wetted_area`` (``BD(93)``);
            ``dwash``: ``q_ratio`` (``DWASH(1..)``), ``epsilon``
            (``DWASH(21..)``).
        stale: ``/POWR/`` words PRPWEF reads without always setting:
            ``sih`` (without a tail) and ``ytemp``.

    Returns:
        ``cdpow``, ``dclpon``, ``dcm`` (``POWER(1)``, ``(21)``, ``(41)``
        onward), the ``/POWR/`` words by the names of :data:`POWR_LAYOUT`
        (scalars as the last angle leaves them), ``kn``, ``cosaiw``
        (``BD(79)``), and ``nalpha`` (``FLC(2)+.5``, which the source
        stores in a local ``NALPHA``: PRPWEF does not declare ``/OVERLY/``).

    Notes:
        For a single or centreline propeller the quarter-chord sweep is
        taken from ``A(69)`` or from ``ANGLES``' ``ARGCS(2)``, both
        radians, and used as degrees (``COS(DEG*SWEEPA)`` and the Figure
        4.1.4.1-5 grid).  Only the two-engine path reads a degree word.
        The slipstream height is clamped with ``ZS = PRPRAD`` where the
        test is ``ZS-ZW > PRPRAD``; the tail arm ``TN`` places the wing's
        MAC with the tail's incidence cosine.  All kept.
    """
    power_in = {k: v for k, v in inputs['power'].items()}
    wing_in = inputs['wing']
    a_block = {int(k): float(v) for k, v in inputs['a'].items()}
    pos = {k: float(v) for k, v in inputs['position'].items()}
    tail_in, vertical_in, body_in = (inputs['tail'], inputs['vertical'],
                                     inputs['body'])
    dwash_in = inputs['dwash']
    srw, cbarr, htpl = float(inputs['sref']), float(inputs['cbarr']), \
        bool(inputs['htpl'])
    alpha_schedule = [float(x) for x in inputs['alpha']]
    nalpha = len(alpha_schedule)
    stale_words = dict(stale or {})
    result: Dict[str, object] = {}
    nengsp, prprad, yp = float(power_in['nengsp']), float(power_in['prprad']), \
        float(power_in['yp'])
    phaloc, phvloc = float(power_in['phaloc']), float(power_in['phvloc'])
    thstcp, aietlp = float(power_in['thstcp']), float(power_in['aietlp'])
    ct, bst0o2, bsto2, bo2 = (float(wing_in[k]) for k in ('ct', 'bst0o2', 'bsto2',
                                                          'bo2'))
    cb, cr = float(wing_in['cb']), float(wing_in['cr'])
    xw, zw, aliw = pos['xw'], pos['zw'], pos['aliw']
    xh, zh, alih, xcg = pos['xh'], pos['zh'], pos['alih'], pos['xcg']
    ar, crstr, alpha0, xbarrw = a_block[120], a_block[10], a_block[134], a_block[161]
    argcs = zerang()

    # The propeller's chord station and the wing's upwash.
    if nengsp != 1.0 and bst0o2 != UNUSED and yp > bo2 - bst0o2:
        crp = cb - (cb - ct) * ((yp - (bo2 - bst0o2)) / bst0o2)
        xbarp = (xw + a_block[62] * (bo2 - bst0o2) + a_block[86] *
                 (yp - (bo2 - bst0o2)) - phaloc) * math.cos(DEG * aliw)
    else:
        ct_ = cb if (nengsp != 1.0 and bst0o2 != UNUSED) else ct
        crp = cr - (cr - ct_) * yp / (bo2 - bst0o2)
        xbarp = (xw + a_block[62] * yp + crp / 4.0 - phaloc) * math.cos(DEG * aliw)
    deuda = -1.0
    if xbarp / crp >= 0.25:
        deuda = _tlinex(_X14161, _X24161, _Y44161, ar, xbarp, 2, 2, 2, 2)
    cosaiw = math.cos(DEG * aliw)
    prprd2 = prprad**2
    srtpco = srw * thstcp / (8.0 * prprd2)
    kn = float(power_in['kn'])
    if kn == UNUSED:
        kn = float(power_in['nopbpe']) * (262.0 * float(power_in['bwapr3']) +
                                   262. * float(power_in['bwapr6']) +
                                   135. * float(power_in['bwapr9'])) / prprad
    if bool(power_in['crot']) and float(power_in['nopbpe']) >= 6:
        cnap80 = _value(_X6111B, _Y6111B, float(power_in['bapr75']), 2, 2)
    else:
        cnap80 = _tlinex(_X1S11B, _X2S11B, _YSR11B, float(power_in['nopbpe']),
                         float(power_in['bapr75']), 0, 2, 2, 2)
    cnap = cnap80 * (1. + .8 * (kn / 80.7 - 1.))
    c1 = _value(_X46113[:18], _C16113, srtpco, 0, 1)
    c2 = _value(_X46113, _C26113, srtpco, 0, 1)
    depdap = c1 + c2 * cnap
    f = _value(_X6111A, _F6111A, srtpco, 0, 2)
    combo1 = PI * nengsp * f * cnap * prprd2 / (RAD * srw)
    combo = nengsp * thstcp
    alphap = aietlp + deuda * (aliw - alpha0)
    alphap = alphap - depdap * alphap
    sinap = math.sin(DEG * alphap)
    zs = phvloc + xbarp * math.tan(DEG * alphap)
    powr_words = {'xbarp': xbarp, 'deuda': deuda, 'prprd2': prprd2,
         'srtpco': srtpco, 'cnap80': cnap80, 'cnap': cnap, 'c1': c1,
         'c2': c2, 'depdap': depdap, 'f': f, 'combo1': combo1,
         'combo': combo, 'alphap': alphap, 'sinap': sinap, 'zs': zs,
         'sih': float(stale_words.get('sih', 0.0)),
         'ytemp': float(stale_words.get('ytemp', 0.0))}

    # The tail's immersed area.
    bo2h, crh, cth = (float(tail_in[k]) for k in ('bo2', 'cr', 'ct'))
    srh, xbarrh = float(tail_in['area']), float(tail_in['xbarr'])
    if htpl:
        cosaih = math.cos(DEG * alih)
        tn = xh + xbarrh * cosaih - (xw + xbarrw * cosaih)
        zheff = zs - zh + tn * math.tan(DEG * alphap)
        powr_words.update({'cosaih': cosaih, 'zheff': zheff})
        if prprd2 < zheff**2:
            sih = siosrh = 0.0
        elif nengsp == 2.0:
            if yp - prprad >= bo2h:
                sih = siosrh = 0.0
            else:
                sih = (bo2h - yp + prprad) / bo2h * srh
                siosrh = sih / srh
        else:
            siosrh, sih = 1.0, srh
            if prprad < bo2h:
                ctih = crh - prprad * (crh - cth) / bo2h
                sih = _sqrt(prprd2 - zheff**2) * (crh + ctih)
                siosrh = sih / srh
                powr_words['ctih'] = ctih
        powr_words.update({'sih': sih, 'siosrh': siosrh,
                  'ytemp': _tlinex(_X14637, _X24637, _Y4637, siosrh,
                                   srtpco, 0, 0, 0, 1)})

    # The wing's immersed area, its sweep and its zero-lift moment.
    bio2 = _sqrt(prprd2 - (zs - (zw - xbarrw * sinap))**2)
    powr_words['bio2'] = bio2
    if nengsp == 2.0:
        sstri = bio2 * crp * 4.0
        astari = _div(4.0 * bio2**2, .50 * sstri)
        trpsi = 1.0
        sweepa = a_block[112] if yp > bo2 - bst0o2 else a_block[106]
        cbarli = crp
        si = sstri
        dcd0s = srtpco * 8.0 / (PI * srw) * (sstri * float(wing_in['cf']))
        if yp > bo2 - bst0o2:
            dcd0s = srtpco * 8.0 / (PI * srw) * (sstri * float(wing_in['d12']))
    else:
        if bio2 > bo2 - bst0o2:
            bst0i2 = bst0o2 - (bo2 - bio2)
            cti = cb - bst0i2 * (cb - ct) / bst0o2
            sst0i = (cb + cti) * bst0i2
            sstri = a_block[1] + sst0i
            scapi = (cr + cb) * (bo2 - bst0o2)
            si = scapi + sst0i
            trs0i = cti / cb
            cbsr0i = cb * _funcb(trs0i)
            cbarli = (scapi * a_block[121] + sst0i * cbsr0i) / si
            argcs[3] = (sstri * a_block[67] + sst0i * a_block[91]) / sstri
            argcs = angles(4, argcs)
            sweepa = argcs[1]
            trpsi = a_block[26] * trs0i
            bstio2 = a_block[23] + bst0i2
            powr_words.update({'bst0i2': bst0i2, 'sst0i': sst0i, 'scapi': scapi,
                      'trs0i': trs0i, 'cbsr0i': cbsr0i})
        else:
            cti = cr - bio2 * ((cr - cb) / (bo2 - bst0o2))
            bstio2 = bio2 - bo2 + bsto2
            sstri = (crstr + cti) * bstio2
            si = (cr + cti) * bio2
            tri = cti / cr
            cbarli = cr * _funcb(tri)
            sweepa = a_block[69]
            trpsi = cti / crstr
            powr_words['tri'] = tri
        astari = 4.0 * (bstio2**2) / sstri
        powr_words.update({'cti': cti, 'bstio2': bstio2})
        dcd0s = (srtpco * 8.0 / (PI * srw)) * (
            float(wing_in['cf']) * sstri + float(tail_in['cf']) * powr_words['sih'] +
            .50 * float(vertical_in['cf']) * float(vertical_in['area']) +
            float(body_in['cf']) * float(body_in['wetted_area']))
    cd0pow = float(wing_in['cd0']) + dcd0s
    rpnob = .5 * nengsp * prprad / bo2
    aak = _value(_X4648A, _Y4648A, srtpco, 0, 2)
    ebroep = _tlinex(_X1648B, _X2648B, _Y4648B, rpnob, srtpco, 2, 0, 2, 2)
    dcmt = thstcp * (pos['zcg'] - phvloc) * nengsp / cbarr
    cossw = math.cos(DEG * sweepa)
    cm0in, cm02 = float(wing_in['cmo']), float(wing_in['cmot'])
    cm0ova = cm0in
    if not (abs(cm0in) < 1.e-10 or abs(cm02) < 1.e-10):
        cm0ova = 0.5 * (cm0in + cm02)
    cm0te0 = astari * cossw**2 * cm0ova / (astari + 2.0 * cossw)
    twista = float(wing_in['twista'])
    cm0i = cm0te0
    if not twista < 1.e-10:
        y = np.asarray(_Y41412).reshape(3, 5, 16).transpose(2, 1, 0)
        cm0i = cm0te0 + float(tlin3x(_X11412, _X21412, _X31412, y, astari,
                                     sweepa, trpsi, 2, 0, 0, 2, 2,
                                     0)) * twista
    dcmq = (srtpco * si / srw * cbarli / cbarr * cm0i) * 8.0 / PI
    xbrsrr = a_block[16] / 4. + (a_block[1] * a_block[32] * a_block[62] + a_block[2] *
                           (a_block[23] * a_block[62] + (a_block[33] - a_block[23]) * a_block[86])) / a_block[3]
    powr_words.update({'sstri': sstri, 'astari': astari, 'trpsi': trpsi,
              'sweepa': sweepa, 'cbarli': cbarli, 'si': si, 'dcd0s': dcd0s,
              'cd0pow': cd0pow, 'rpnob': rpnob, 'aak': aak,
              'ebroep': ebroep, 'dcmt': dcmt, 'cossw': cossw,
              'cm0in': cm0in, 'cm02': cm02, 'cm0ova': cm0ova,
              'cm0te0': cm0te0, 'cm0i': cm0i, 'dcmq': dcmq,
              'xbrsrr': xbrsrr})

    # The angle loop.
    arrays = {k: [0.0] * nalpha for k in (
        'dclt', 'dclnp', 'dclq', 'dclaw', 'dclhq', 'dcmnp', 'dcml', 'dcmhq',
        'dcmhe', 'dclhe', 'cdpow', 'dclpon', 'dcm')}
    cl_w = [float(x) for x in wing_in['cl']]
    cdl = [float(x) for x in wing_in['cdl']]
    cla = float(wing_in['cla'])
    dlh = float(stale_words.get('dlh', 0.0))
    for j in range(nalpha):
        alphat = alpha_schedule[j] + aietlp
        cosat, sinat = math.cos(DEG * alphat), math.sin(DEG * alphat)
        alphap = alphat + deuda * (aliw + alpha_schedule[j] - alpha0)
        cnp = cnap * alphap / RAD * PI * prprad**2 / srw
        ebar = ebroep * depdap * alphap
        arrays['dclnp'][j] = combo1 * alphap * cosat
        ep = depdap * alphap
        alphap = alphap - ep
        sinap = math.sin(DEG * alphap)
        zs = phvloc + xbarp * math.tan(DEG * alphap)
        if zs - zw > prprad:
            zs = prprad
        bio2 = _sqrt(prprd2 - (zs - zw)**2)
        if deuda == -1.0:
            zs = phvloc
        immersed = True
        if nengsp == 2.0:
            sstri = 4.0 * bio2 * crp
            astari = _div(4.0 * bio2**2, .50 * sstri)
        elif bio2 > bo2 - bst0o2:
            bst0i2 = bst0o2 - (bo2 - bio2)
            cti = cb - bst0i2 * (cb - ct) / bst0o2
            sst0i = (cb + cti) * bst0i2
            sstri = a_block[1] + sst0i
            bstio2 = a_block[23] + bst0i2
            astari = 4.0 * (bstio2**2) / sstri
            powr_words.update({'bst0i2': bst0i2, 'cti': cti, 'sst0i': sst0i,
                      'bstio2': bstio2})
        elif bio2 <= bo2 - bsto2:
            immersed = False
        else:
            cti = cr - bio2 * ((cr - cb) / (bo2 - bst0o2))
            bstio2 = bio2 - bo2 + bsto2
            sstri = (crstr + cti) * bstio2
            astari = 4.0 * (bstio2**2) / sstri
            powr_words.update({'cti': cti, 'bstio2': bstio2})
        if immersed:
            bs1 = _tlinex(_X16112, _X26112, _Y46112, astari, ar, 0, 0, 2, 0)
            bs2 = _value(_XB6112, _DM2, srtpco, 0, 2)
            bs3 = _tlinex(_DC1, _DC2, _DC3, bs2, bs1, 0, 0, 0, 0)
            ak1 = _tlinex(_D3, _D4, _AK6112, srtpco, bs3, 0, 0, 0, 0)
            powr_words.update({'bs1': bs1, 'bs2': bs2, 'bs3': bs3, 'ak1': ak1})
            immersed = deuda != -1.0
        if immersed:
            arrays['dclq'][j] = 8.0 * powr_words['ak1'] * srtpco * sstri / (
                PI * srw) * cl_w[j]
            delalp = -ep / (1.0 + deuda)
            arrays['dclaw'][j] = cla * delalp * sstri * powr_words['ak1'] / srw * (
                1.0 + 8.0 * srtpco / PI)
            powr_words['delalp'] = delalp
        arrays['dclt'][j] = combo * sinat
        if htpl:
            powr_words['dxhmac'] = xh + xbarrh * math.cos(DEG * alih) - phaloc
            zht = zh - phvloc + powr_words['dxhmac'] * math.tan(DEG * aietlp)
            zhtorp = zht / prprad
            eps = float(dwash_in['epsilon'][j])
            if nengsp > 1.0:
                step1 = _tlinex(_X1639A, _X2639A, _Y4639A, eps, srtpco,
                                2, 2, 2, 2)
                flpsup = _value(_X4639B, _Y4639B, step1, 0, 1)
                depowr = _tlinex(_X1639C, _X2639C, _Y4639C, zhtorp, flpsup,
                                 0, 0, 2, 1)
            else:
                step1 = _tlinex(_X1638A, _X2638A, _Y4638A, eps, srtpco,
                                0, 0, 2, 2)
                depowr = _tlinex(_X1638B, _X2638B, _Y4638B, zhtorp, step1,
                                 0, 1, 2, 1)
            clh, clalph = tbfunx(alpha_schedule, tail_in['cl'], alpha_schedule[j], 1, 1)
            cosaih = math.cos(DEG * alih)
            tn = xh + xbarrh * cosaih - (xw + xbarrw * cosaih)
            epowr = eps + depowr
            zheff = zs - zh + tn * math.tan(DEG * (alphap - epowr))
            zhorp = zheff / prprad
            dqhoqi = _tlinex(_XT4637, _Z24637, _YF4637, abs(zhorp),
                             powr_words['ytemp'], 0, 0, 2, 1)
            arrays['dclhq'][j] = dqhoqi * clh
            arrays['dclhe'][j] = -clalph * depowr * (
                float(dwash_in['q_ratio'][j]) + dqhoqi)
            dlh = xh + xbarrh * cosaih - xcg
            arrays['dcmhe'][j] = -dlh * arrays['dclhe'][j] / cbarr
            powr_words.update({'zht': zht, 'zhtorp': zhtorp, 'step1': step1,
                      'epowr': epowr, 'zheff': zheff, 'zhorp': zhorp,
                      'dqhoqi': dqhoqi})
        arrays['dclpon'][j] = (arrays['dclt'][j] + arrays['dclnp'][j] +
                               arrays['dclq'][j] + arrays['dclaw'][j] +
                               arrays['dclhq'][j] + arrays['dclhe'][j])
        arrays['dcmnp'][j] = arrays['dclnp'][j] * (xcg - phaloc) / (
            cbarr * cosat)
        xcp = (xw + xbrsrr * cosaiw) - xcg
        arrays['dcml'][j] = -(arrays['dclq'][j] + arrays['dclaw'][j]) * \
            xcp / cbarr
        arrays['dcmhq'][j] = -dlh * arrays['dclhq'][j] / cbarr
        arrays['dcm'][j] = (dcmt + arrays['dcmnp'][j] + dcmq +
                            arrays['dcml'][j] + arrays['dcmhq'][j] +
                            arrays['dcmhe'][j])
        clp = (thstcp * sinat + cnp * cosat) * nengsp
        if cl_w[j] == 0.0:
            cdlrat = 0.0
        else:
            clww = cl_w[j] + arrays['dclq'][j]
            cdlrat = ((clww / cl_w[j])**2 *
                      (1. + _PISQRD * ar * ebar / (180. * clww)) +
                      aak * (bo2 / prprad * clp / cl_w[j])**2)
            powr_words['clww'] = clww
        cdlpow = cdlrat * cdl[j]
        arrays['cdpow'][j] = dcd0s + cdlpow - cdl[j]
        powr_words.update({'alphat': alphat, 'cosat': cosat, 'sinat': sinat,
                  'alphap': alphap, 'cnp': cnp, 'ebar': ebar, 'ep': ep,
                  'sinap': sinap, 'zs': zs, 'bio2': bio2, 'sstri': sstri,
                  'astari': astari, 'xcp': xcp, 'dlh': dlh, 'clp': clp,
                  'cdlrat': cdlrat, 'cdlpow': cdlpow})
    result.update(arrays)
    result.update({'powr': powr_words, 'kn': kn, 'cosaiw': cosaiw, 'nalpha': nalpha,
                   'argcs': argcs, 'method': 'legacy_prpwef'})
    return result


def calculate_m13o15(inputs: Mapping[str, object],
                     stale: Optional[Mapping[str, float]] = None
                     ) -> Dict[str, object]:
    """Translate M13O15: PRPWEF, then the power increments' CN, CA and
    slopes (``POWER(61)``, ``(81)``, ``(101)``, ``(121)`` onward) from its
    lift, drag and moment increments."""
    r = calculate_prpwef(inputs, stale)
    alpha = np.asarray(inputs['alpha'], dtype=float)
    ca_, sa_ = np.cos(alpha / RAD), np.sin(alpha / RAD)
    cl, cd = np.asarray(r['dclpon']), np.asarray(r['cdpow'])
    r['cn'] = cl * ca_ + cd * sa_
    r['ca'] = cd * ca_ - cl * sa_
    r['cla'] = np.array([tbfunx(alpha, cl, x, 0, 0)[1] for x in alpha])
    r['cma'] = np.array([tbfunx(alpha, r['dcm'], x, 0, 0)[1]
                         for x in alpha])
    r['method'] = 'legacy_m13o15'
    return r
