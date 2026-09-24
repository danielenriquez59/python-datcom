"""
Asymmetric flap, spoiler and all-movable-tail rolling and yawing moments:
LATFLP.

LATFLP forms, for each pair of left and right deflections:

* plain flaps (``STYPE`` 4): the section lift of four spanwise strips
  (Figures 6.1.1.1-39 and -40), the span-loading derivative of Figure
  6.2.1.1-23 and the rolling moment, with the yawing moment from
  Figure 6.2.2.1-9.  A flap that stops short of the tip is built by
  superposition: a panel from the flap's inner edge to the tip less a
  panel from its outer edge to the tip;
* spoilers (``STYPE`` 1-3): the rolling moment from the spoiler lift of
  Figure 6.1.1.1-51B and the effective span of 6.2.1.1-26A, the slot
  factor of -26B, and the yawing moment of Figures 6.2.2.1-10 to -12;
* a differentially deflected horizontal tail (``STYPE`` 5): the rolling
  moment of Section 6.2.1.2.

The ``/FLAPIN/``, ``FLA`` (``/POWR/`` 60-104) and ``HT`` words are
mirrored as flat 1-based arrays with the source's EQUIVALENCE offsets.

The tables were extracted from the source DATA statements by parsing, with
``tools/fortran_data.py``.

Reference: datcom-legacy/datcom_2000/latflp.f
"""

import math
from typing import Dict, Mapping

import numpy as np

from pydatcom.utils.constants import PI, RAD, UNUSED
from pydatcom.utils.legacy_interp import interx
from pydatcom.utils.legacy_tables import tlin3x, tlin4x, tlinex

_X2128A = np.array([
    0.0, 0.02, 0.04, 0.06, 0.08, 0.1, 0.12, 0.14, 0.16, 0.18,
    0.2,
])
_X1128A = np.array([
    1.0, 10.0, 100.0,
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
_X21123 = np.array([
    0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9,
    1.0,
])
_X11123 = np.array([
    -40.0, 0.0, 40.0, 60.0,
])
_X31123 = np.array([
    2.0, 4.0, 8.0,
])
_X41123 = np.array([
    0.0, 0.5, 1.0,
])
_Y1 = np.array([
    0.0, 0.005, 0.014, 0.03, 0.055, 0.087, 0.126, 0.163, 0.197, 0.222,
    0.245, 0.0, 0.005, 0.017, 0.036, 0.063, 0.1, 0.144, 0.185, 0.224,
    0.262, 0.282, 0.0, 0.005, 0.017, 0.036, 0.064, 0.103, 0.149, 0.192,
    0.23, 0.264, 0.282, 0.0, 0.005, 0.017, 0.036, 0.062, 0.098, 0.136,
    0.178, 0.215, 0.243, 0.259,
])
_Y2 = np.array([
    0.0, 0.007, 0.027, 0.054, 0.097, 0.145, 0.197, 0.243, 0.292, 0.333,
    0.363, 0.0, 0.008, 0.029, 0.064, 0.113, 0.173, 0.235, 0.3, 0.363,
    0.414, 0.445, 0.0, 0.008, 0.029, 0.064, 0.119, 0.18, 0.245, 0.313,
    0.374, 0.417, 0.443, 0.0, 0.008, 0.029, 0.063, 0.108, 0.163, 0.22,
    0.273, 0.319, 0.35, 0.363,
])
_Y3 = np.array([
    0.0, 0.006, 0.032, 0.08, 0.145, 0.212, 0.275, 0.336, 0.394, 0.446,
    0.485, 0.0, 0.012, 0.043, 0.093, 0.165, 0.255, 0.354, 0.445, 0.525,
    0.587, 0.627, 0.0, 0.013, 0.045, 0.106, 0.189, 0.28, 0.375, 0.465,
    0.542, 0.592, 0.612, 0.0, 0.013, 0.045, 0.106, 0.187, 0.258, 0.327,
    0.39, 0.443, 0.475, 0.485,
])
_Y4 = np.array([
    0.0, 0.006, 0.017, 0.038, 0.066, 0.1, 0.142, 0.185, 0.232, 0.273,
    0.293, 0.0, 0.007, 0.023, 0.048, 0.079, 0.114, 0.154, 0.2, 0.249,
    0.284, 0.311, 0.0, 0.007, 0.023, 0.048, 0.079, 0.114, 0.154, 0.198,
    0.24, 0.276, 0.304, 0.0, 0.007, 0.023, 0.043, 0.072, 0.106, 0.148,
    0.19, 0.232, 0.264, 0.286,
])
_Y5 = np.array([
    0.0, 0.006, 0.025, 0.06, 0.1, 0.151, 0.214, 0.294, 0.37, 0.43,
    0.479, 0.0, 0.009, 0.033, 0.075, 0.125, 0.185, 0.259, 0.348, 0.427,
    0.491, 0.538, 0.0, 0.009, 0.033, 0.075, 0.125, 0.185, 0.259, 0.34,
    0.413, 0.472, 0.515, 0.0, 0.008, 0.026, 0.065, 0.113, 0.17, 0.231,
    0.294, 0.35, 0.398, 0.432,
])
_Y6 = np.array([
    0.0, 0.007, 0.034, 0.072, 0.125, 0.187, 0.27, 0.371, 0.48, 0.58,
    0.662, 0.0, 0.014, 0.052, 0.108, 0.18, 0.275, 0.39, 0.51, 0.62,
    0.718, 0.8, 0.0, 0.014, 0.052, 0.113, 0.19, 0.28, 0.385, 0.485,
    0.578, 0.658, 0.729, 0.0, 0.014, 0.052, 0.103, 0.163, 0.234, 0.308,
    0.379, 0.445, 0.504, 0.554,
])
_Y7 = np.array([
    0.0, 0.003, 0.015, 0.037, 0.065, 0.103, 0.147, 0.19, 0.237, 0.28,
    0.306, 0.0, 0.004, 0.02, 0.045, 0.079, 0.118, 0.159, 0.205, 0.252,
    0.291, 0.315, 0.0, 0.004, 0.02, 0.045, 0.079, 0.114, 0.156, 0.2,
    0.246, 0.286, 0.31, 0.0, 0.004, 0.02, 0.045, 0.073, 0.108, 0.147,
    0.187, 0.23, 0.268, 0.291,
])
_Y8 = np.array([
    0.0, 0.008, 0.033, 0.063, 0.108, 0.162, 0.225, 0.303, 0.389, 0.453,
    0.5, 0.0, 0.008, 0.033, 0.072, 0.122, 0.187, 0.26, 0.348, 0.437,
    0.504, 0.558, 0.0, 0.008, 0.033, 0.072, 0.122, 0.187, 0.258, 0.335,
    0.413, 0.474, 0.522, 0.0, 0.008, 0.033, 0.07, 0.11, 0.165, 0.229,
    0.29, 0.347, 0.4, 0.445,
])
_Y9 = np.array([
    0.0, 0.01, 0.035, 0.07, 0.118, 0.192, 0.291, 0.405, 0.52, 0.626,
    0.721, 0.0, 0.014, 0.049, 0.104, 0.179, 0.275, 0.393, 0.52, 0.653,
    0.768, 0.85, 0.0, 0.014, 0.05, 0.109, 0.186, 0.284, 0.393, 0.505,
    0.602, 0.69, 0.763, 0.0, 0.014, 0.04, 0.094, 0.155, 0.23, 0.306,
    0.38, 0.458, 0.519, 0.578,
])
_X22219 = np.array([
    0.0, 0.2, 0.4, 0.6, 0.7, 0.8, 0.9,
])
_X12219 = np.array([
    3.0, 4.0, 6.0, 8.0,
])
_X32219 = np.array([
    0.25, 0.5, 0.75, 1.0,
])
_Y62219 = np.array([
    -0.285, -0.278, -0.275, -0.282, -0.289, -0.3, -0.315, -0.234, -0.223, -0.217,
    -0.22, -0.223, -0.226, -0.23, -0.16, -0.152, -0.146, -0.146, -0.15, -0.162,
    -0.175, -0.12, -0.116, -0.11, -0.106, -0.112, -0.125, -0.14, -0.338, -0.336,
    -0.338, -0.34, -0.347, -0.355, -0.364, -0.25, -0.25, -0.252, -0.258, -0.264,
    -0.27, -0.28, -0.17, -0.17, -0.172, -0.182, -0.19, -0.202, -0.219, -0.13,
    -0.131, -0.132, -0.136, -0.145, -0.158, -0.175, -0.327, -0.327, -0.33, -0.34,
    -0.362, -0.39, -0.43, -0.261, -0.261, -0.265, -0.28, -0.294, -0.307, -0.323,
    -0.179, -0.179, -0.186, -0.2, -0.212, -0.225, -0.24, -0.138, -0.138, -0.14,
    -0.159, -0.172, -0.184, -0.195, -0.361, -0.361, -0.361, -0.361, -0.361, -0.361,
    -0.361, -0.262, -0.267, -0.279, -0.295, -0.305, -0.315, -0.325, -0.182, -0.19,
    -0.2, -0.219, -0.232, -0.248, -0.27, -0.145, -0.15, -0.16, -0.175, -0.187,
    -0.202, -0.223,
])
_X2110A = np.array([
    0.0, 0.3, 0.5, 0.7, 0.9, 1.0,
])
_X1110A = np.array([
    0.4, 0.6, 0.8, 0.9, 1.0,
])
_Y2110A = np.array([
    0.0, 0.2, 0.24, 0.32, 0.44, 0.49, 0.0, 0.5, 0.8, 1.2,
    1.53, 1.7, 0.0, 0.91, 1.45, 2.0, 2.42, 2.68, 0.0, 1.38,
    2.1, 2.75, 3.3, 3.55, 0.0, 1.68, 2.6, 3.5, 4.2, 4.4,
])
_X2110B = np.array([
    0.0, 1.0, 2.0, 3.0, 4.0, 5.0,
])
_X1110B = np.array([
    2.0, 3.0, 6.0, 10.0,
])
_Y2110B = np.array([
    0.0, 0.5, 1.0, 1.35, 1.65, 1.85, 0.0, 0.7, 1.35, 1.85,
    2.23, 2.48, 0.0, 0.98, 1.75, 2.45, 3.05, 3.44, 0.0, 1.0,
    1.95, 2.95, 3.9, 4.88,
])
_X2110C = np.array([
    0.0, 5.0,
])
_X1110C = np.array([
    0.4, 1.0,
])
_Y2110C = np.array([
    0.0, 5.05, 0.0, 8.4,
])
_X2110D = np.array([
    0.0, 5.0,
])
_X1110D = np.array([
    0.5, 0.6, 0.7, 0.8,
])
_Y2110D = np.array([
    0.0, 0.31, 0.0, 0.268, 0.0, 0.2, 0.0, 0.146,
])
_X211A0 = np.array([
    0.0, 0.1, 0.2, 0.3, 0.4,
])
_X111A0 = np.array([
    0.4, 0.6, 0.8, 1.0,
])
_Y111A0 = np.array([
    0.0, 0.55, 1.0, 1.35, 1.5, 0.0, 0.79, 1.45, 2.0, 2.5,
    0.0, 1.0, 2.0, 2.8, 3.55, 0.0, 1.6, 3.0, 4.25, 5.3,
])
_X211A1 = np.array([
    0.4, 0.5, 0.6,
])
_X111A1 = np.array([
    0.6, 0.8, 1.0,
])
_Y111A1 = np.array([
    2.5, 2.75, 2.8, 3.55, 4.2, 4.4, 5.3, 6.0, 6.55,
])
_X211A2 = np.array([
    0.6, 0.7, 0.8,
])
_X111A2 = np.array([
    0.8, 1.0,
])
_Y111A2 = np.array([
    4.4, 4.5, 4.3, 6.55, 6.86, 7.0,
])
_X111A3 = np.array([
    0.8, 0.9, 1.0,
])
_Y111A3 = np.array([
    7.0, 6.8, 6.28,
])
_X1111B = np.array([
    2.0, 4.0, 6.0,
])
_X2111B = np.array([
    0.0, 7.0,
])
_Y2111B = np.array([
    0.0, 2.85, 0.0, 3.5, 0.0, 3.75,
])
_X1111C = np.array([
    0.0, 0.4, 0.6, 1.0,
])
_X2111C = np.array([
    0.0, 5.0,
])
_Y2111C = np.array([
    0.0, 7.45, 0.0, 8.6, 0.0, 10.0, 0.0, 16.5,
])
_X1111D = np.array([
    20.0, 45.0, 60.0,
])
_X2111D = np.array([
    0.0, 8.0,
])
_Y2111D = np.array([
    0.0, 0.153, 0.0, 0.172, 0.0, 0.193,
])
_X2111E = np.array([
    0.0, 0.1,
])
_Y2111E = np.array([
    0.0, 0.15,
])
_X22112 = np.array([
    0.175, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1,
    1.15, 1.25,
])
_Y22112 = np.array([
    2.3, 2.33, 2.34, 2.32, 2.29, 2.22, 2.09, 1.81, 1.53, 1.38,
    1.37, 1.49,
])
_X1126A = np.array([
    0.0, 2.5, 5.0, 7.5, 10.0, 12.5, 15.0, 17.5, 20.0, 25.0,
    30.0, 35.0, 40.0, 45.0,
])
_Y1126A = np.array([
    0.0, 9.0, 13.0, 16.0, 18.5, 20.6, 22.4, 23.7, 25.0, 27.3,
    29.0, 30.4, 31.2, 31.6,
])
_X1132B = np.array([
    0.5, 0.6, 0.7, 0.8,
])
_X2132B = np.array([
    0.06, 0.07, 0.08, 0.09, 0.1, 0.11, 0.12, 0.14, 0.16, 0.18,
    0.19,
])
_Y1132B = np.array([
    0.012, 0.046, 0.068, 0.086, 0.102, 0.116, 0.126, 0.143, 0.158, 0.167,
    0.17, 0.046, 0.064, 0.08, 0.095, 0.108, 0.121, 0.13, 0.149, 0.164,
    0.175, 0.18, 0.054, 0.072, 0.088, 0.102, 0.116, 0.128, 0.137, 0.154,
    0.17, 0.181, 0.186, 0.08, 0.093, 0.105, 0.117, 0.126, 0.136, 0.144,
    0.16, 0.175, 0.189, 0.196,
])
_X6226B = np.array([
    0.68, 0.77, 1.0, 1.17, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0,
])
_Y6226B = np.array([
    0.7, 1.0, 1.36, 1.5, 1.65, 1.73, 1.72, 1.64, 1.5, 1.33,
])
_X62122 = np.array([
    0.165, 0.2, 0.25, 0.3, 0.35, 0.4,
])
_Y62122 = np.array([
    0.83, 0.818, 0.779, 0.738, 0.68, 0.615,
])

_Y21123 = np.concatenate([_Y1, _Y2, _Y3, _Y4, _Y5, _Y6, _Y7, _Y8, _Y9])


def _tl(x1, x2, flat, q1, q2, l1=0, l2=0, u1=0, u2=0):
    grid = np.asarray(flat, dtype=float).reshape(len(x1), len(x2)).T
    return float(tlinex(x1, x2, grid, q1, q2, l1, l2, u1, u2))


def _ix1(x, y, q, l1=0, u1=0):
    return float(interx(1, x, [q], [len(x)], y, lind=len(x), lx1l=l1,
                        lx1u=u1))


def _span_loading(sweep, eta, bak, taper):
    """Figure 6.2.1.1-23, ``Y(11,4,3,3)`` in source order."""
    y = _Y21123.reshape((len(_X21123), len(_X11123), len(_X31123),
                         len(_X41123)), order='F')
    return float(tlin4x(_X11123, _X21123, _X31123, _X41123, y, sweep, eta,
                        bak, taper))


def _yaw_factor(aw, eta, taper):
    """Figure 6.2.2.1-9, ``Y(7,4,4)`` in source order."""
    y = _Y62219.reshape((len(_X22219), len(_X12219), len(_X32219)),
                        order='F')
    return float(tlin3x(_X12219, _X22219, _X32219, y, aw, eta, taper))


def calculate_latflp(data: Mapping[str, object]) -> Dict[str, object]:
    """Translate LATFLP: rolling and yawing moments of asymmetric controls.

    Args:
        data: By name: ``transn``, ``nalpha``, ``mach`` (``FLC(M+2)``),
            ``rl`` (``FLC(M+42)``), ``alpha`` (``FLC(23..)``), ``sref``,
            ``blref``, ``surface`` (the wing's words: ``bo2``
            (``WINGIN(4)``), ``cr``, ``tovc``, ``clasec``
            (``WINGIN(M+20)``, or ``WINGIN(69)/.8`` transonic), ``area``
            (``A(4)``), ``cbarex``, ``sweple``, ``tanle``, ``cosc4``,
            ``tanc4``, ``swepte``, ``coste``, ``tante``, ``taprw``,
            ``aw``, ``cla`` (``WING(101)``), ``clw`` (``WING(21..)``)),
            ``htail`` (``bo2h`` (``HTIN(4)``), ``bo2hst`` (``HTIN(3)``),
            ``sh`` (``AHT(3)``), ``clah`` (``HT(101)``)), ``gd``
            (``TCD(43..46)``), ``dedalp``, ``rivbh``, ``gamvr`` (per
            angle of attack), ``f`` (the 69 ``/FLAPIN/`` words), ``fla``
            (the 45 ``FLA`` words), ``ht201`` (``HT`` 201-230), ``clrol``
            (``WING`` 201-400) and ``cn`` (``BODY`` 201-400), as they stood.

    Returns:
        ``f``, ``fla`` (1-based lists with a leading pad), ``ht201``,
        ``clrol`` and ``cn``.

    Notes:
        Kept as executed: the spoiler yawing moment of Figure
        6.2.2.1-10A is read at half the rolling-moment derivative where
        the chart's abscissa is the spoiler span (the span is formed only
        after the lookup), and the plain flap's outboard loading stores
        ``BCLOKO`` as zero.
    """
    s = {k: v for k, v in data['surface'].items()}
    f = [0.0] + [float(v) for v in data['f']]
    fla = [0.0] + [float(v) for v in data['fla']]
    ht = {200 + k: float(v) for k, v in enumerate(data['ht201'], 1)}
    clrol = [0.0] + [float(v) for v in data['clrol']]
    cn = [0.0] + [float(v) for v in data['cn']]
    nalpha = int(data['nalpha'])
    alpha = [0.0] + [float(v) for v in data['alpha']]
    clw = [0.0] + [float(v) for v in s['clw']]
    sref, blref = float(data['sref']), float(data['blref'])
    # EQUIVALENCE offsets into /FLAPIN/.
    deltal = lambda deflection_index: f[18 + deflection_index]  # noqa: E731
    deltar = lambda deflection_index: f[28 + deflection_index]  # noqa: E731
    stype, xsprme = f[18], f[59]
    ndelta = int(f[16] + 0.5)
    rf = float(data['rl'])
    transn = bool(data['transn'])
    mach = 0.6 if transn else float(data['mach'])
    beta = math.sqrt(1. - mach ** 2)
    clasec = float(s['clasec'])
    bo2, aw, taprw = float(s['bo2']), float(s['aw']), float(s['taprw'])
    scale = (2 * bo2 * float(s['area'])) / (blref * sref)

    if stype == 5.:
        _tail_roll(data, f, clrol, nalpha, sref, blref, aw)
    else:
        arg1 = float(s['tanc4']) / beta
        sweepb = math.atan(arg1) * RAD
        fla[1] = sweepb
        kc = clasec * beta / (2. * PI) * RAD
        deln4 = 0.25 * (f[15] - f[14]) / bo2
        eta = [0.0] * 6
        eta[1] = f[14] / bo2
        eta[5] = f[15] / bo2
        if stype == 4.:
            _plain_flap(s, f, fla, ht, cn, clw, alpha, nalpha, ndelta,
                        deltal, deltar, beta, kc, sweepb, deln4, eta, rf,
                        clasec, scale)
        else:
            _spoiler(s, f, fla, ht, ndelta, stype, xsprme, beta, kc, deln4,
                     eta, mach, scale)

    for deflection_index in range(1, ndelta + 1):
        ht[200 + deflection_index] = (deltal(deflection_index) -
                                      deltar(deflection_index))
    return {'f': f, 'fla': fla,
            'ht201': [ht[200 + word_index] for word_index in range(1, 31)],
            'clrol': clrol[1:], 'cn': cn[1:]}


def _plain_flap(s, f, fla, ht, cn, clw, alpha, nalpha, ndelta, deltal,
                deltar, beta, kc, sweepb, deln4, eta, rf, clasec, scale):
    """Labels 1000-1150: the plain flap, with the tip superposition."""
    bo2, aw, taprw = float(s['bo2']), float(s['aw']), float(s['taprw'])
    cr, tovc, cla = float(s['cr']), float(s['tovc']), float(s['cla'])
    tante, tanle = float(s['tante']), float(s['tanle'])
    cf, chrd, cfoc = [0.0] * 6, [0.0] * 6, [0.0] * 6
    ans, cldpm = [0.0] * 6, [0.0] * 5
    cldthy, cldoct, dcll, dclr = [0.0] * 6, [0.0] * 6, [0.0] * 6, [0.0] * 6
    aldl, aldr = [0.0] * 5, [0.0] * 5
    cldlt, cldrt, clrolt = [0.0] * 11, [0.0] * 11, [0.0] * 11
    cntemp = {1: [0.0] * 201, 2: [0.0] * 201}
    tipcal = False
    level = 1
    saved = cft = None
    arg1 = (f[12] - f[13]) / (4. * deln4)
    arg2 = (tante - tanle) * bo2
    while True:
        # Label 1000: strip geometry and loading.
        for strip_index in range(1, 6):
            if strip_index != 1:
                eta[strip_index] = eta[strip_index - 1] + deln4
            cf[strip_index] = f[12] - arg1 * (eta[strip_index] - eta[1])
            chrd[strip_index] = cr + eta[strip_index] * arg2
            cfoc[strip_index] = cf[strip_index] / chrd[strip_index]
        arg1 = beta * aw / kc
        arg2 = eta[1]
        for span_station in range(1, 6):
            ans[span_station] = _span_loading(sweepb, arg2, arg1, taprw)
            arg2 = arg2 + deln4
        for span_station in range(1, 5):
            cldpm[span_station] = ((ans[span_station + 1] - ans[span_station]) *
                                   kc / beta)
        fla[2] = ans[1]
        fla[3] = ans[5]
        fla[4] = fla[2]
        fla[3] = fla[2] - fla[2]
        fla[5] = fla[4] * kc / beta
        arg1 = math.log10(rf * float(s['cbarex']))
        cloclt = _tl(_X1128A, _X2128A, _Y1128A, arg1, f[11])
        for deflection_index in range(1, ndelta + 1):
            arg = abs(deltal(deflection_index))
            arg1 = abs(deltar(deflection_index))
            for strip_index in range(1, 6):
                cldthy[strip_index] = _tl(_X1125A, _X2125A, _Y1125A, tovc,
                                          cfoc[strip_index])
                cldoct[strip_index] = _tl(_X1125B, _X2125B, _Y1125B, cloclt,
                                          cfoc[strip_index])
                kprml = _tl(_X11126, _X21126, _Y11126, cfoc[strip_index], arg)
                arg2 = cldoct[strip_index] * cldthy[strip_index]
                dcll[strip_index] = arg2 * kprml
                kprmr = _tl(_X11126, _X21126, _Y11126, cfoc[strip_index], arg1)
                dclr[strip_index] = arg2 * kprmr
                if strip_index == 1:
                    continue
                pair = strip_index - 1
                aldl[pair] = abs(((dcll[strip_index] + dcll[pair]) / 2.0) /
                                 clasec) / RAD
                aldr[pair] = abs(((dclr[strip_index] + dclr[pair]) / 2.0) /
                                 clasec) / RAD
            cldl = (cldpm[1] * aldl[1] + cldpm[2] * aldl[2] +
                    cldpm[3] * aldl[3] + cldpm[4] * aldl[4])
            cldr = (cldpm[1] * aldr[1] + cldpm[2] * aldr[2] +
                    cldpm[3] * aldr[3] + cldpm[4] * aldr[4])
            roll = (cldl * deltal(deflection_index) / RAD -
                    cldr * deltar(deflection_index) / RAD) / 2. * scale
            if tipcal:
                cldlt[deflection_index], cldrt[deflection_index], \
                    clrolt[deflection_index] = cldl, cldr, roll
            else:
                fla[5 + deflection_index], fla[15 + deflection_index], \
                    ht[210 + deflection_index] = cldl, cldr, roll
        if eta[5] < .98:
            # The flap stops short of the tip: a panel from its inner edge
            # to the tip, then one from its outer edge (label 1070).
            tipcal = True
            saved = (f[12], f[13], f[15], f[14], fla[2], fla[3], fla[4],
                     fla[5])
            level = 1
            deln4 = 0.25 * (bo2 - f[14]) / bo2
            eta[1] = f[14] / bo2
            eta[5] = 1.0
            cft = f[12] - (f[12] - f[13]) * (bo2 - f[14]) / (f[15] - f[14])
            arg1 = (f[12] - cft) / (4. * deln4)
            arg2 = (tante - tanle) * bo2
            continue
        # Label 1100.
        kyaw = _yaw_factor(aw, eta[1], taprw)
        fla[45] = kyaw
        if not tipcal:
            nn = 0
            for angle_slot in range(1, nalpha + 1):
                for deflection_index in range(1, ndelta + 1):
                    nn += 1
                    cn[nn] = kyaw * clw[angle_slot] * ht[210 + deflection_index]
                    if abs(clw[angle_slot]) == UNUSED:
                        cn[nn] = (kyaw * alpha[angle_slot] * cla *
                                  ht[210 + deflection_index])
            return
        nn = 0
        for angle_slot in range(1, nalpha + 1):
            for deflection_index in range(1, ndelta + 1):
                nn += 1
                cntemp[level][nn] = kyaw * clw[angle_slot] * clrolt[
                    deflection_index]
                if abs(clw[angle_slot]) == UNUSED:
                    cntemp[level][nn] = (kyaw * alpha[angle_slot] * cla *
                                         clrolt[deflection_index])
        level += 1
        if level == 2:
            f[12] = f[13]
            deln4 = 0.25 * (bo2 - f[15]) / bo2
            eta[1] = f[15] / bo2
            arg1 = (f[12] - cft) / (4. * deln4)
            arg2 = (tante - tanle) * bo2
            for deflection_index in range(1, ndelta + 1):
                ht[210 + deflection_index] = clrolt[deflection_index]
            continue
        nn = 0
        for angle_slot in range(1, nalpha + 1):
            for deflection_index in range(1, ndelta + 1):
                nn += 1
                cn[nn] = cntemp[1][nn] - cntemp[2][nn]
        (f[12], f[13], f[15], f[14], fla[2], fla[3], fla[4],
         fla[5]) = saved
        for deflection_index in range(1, ndelta + 1):
            ht[210 + deflection_index] = (ht[210 + deflection_index] -
                                          clrolt[deflection_index])
        return


def _spoiler(s, f, fla, ht, ndelta, stype, xsprme, beta, kc, deln4, eta,
             mach, scale):
    """Labels 1160-1260: spoilers."""
    aw, taprw = float(s['aw']), float(s['taprw'])
    swepte, sweple = float(s['swepte']), float(s['sweple'])
    xsoc = lambda deflection_index: f[48 + deflection_index]  # noqa: E731
    dsoc = lambda deflection_index: f[38 + deflection_index]  # noqa: E731
    hsoc = lambda deflection_index: f[59 + deflection_index]  # noqa: E731
    arg1 = 0.75 - (1. - xsprme)
    arg2 = (1. - taprw) / (1. + taprw)
    tansi = float(s['tanc4']) - 4.0 * arg1 / aw * arg2
    sbacki = math.atan(tansi) * RAD
    fla[36] = sbacki
    thetai = _ix1(_X1126A, _Y1126A, sbacki)
    fla[37] = thetai
    arg1 = 4. * (1. - xsprme) / (aw * (1. + taprw))
    arg2 = 1. - (1. - taprw) * eta[1]
    arg3 = (float(s['coste']) * math.sin(thetai / RAD) /
            math.cos((swepte + thetai) / RAD))
    fla[39] = arg1 * arg2 * arg3
    arg2 = 1. - (1. - taprw) * eta[5]
    fla[38] = arg1 * arg2 * arg3
    fla[40] = eta[1] + fla[39]
    fla[41] = eta[5] + fla[38]
    arg4 = beta * aw / kc
    fla[42] = _span_loading(sbacki, fla[40], arg4, taprw)
    fla[43] = _span_loading(sbacki, fla[41], arg4, taprw)
    fla[5] = kc * (fla[43] - fla[42]) / beta
    arg1 = fla[5] / 2.0
    for deflection_index in range(1, ndelta + 1):
        deltas = _tl(_X1132B, _X2132B, _Y1132B, xsoc(deflection_index),
                     hsoc(deflection_index))
        ht[210 + deflection_index] = arg1 * deltas * scale
    if stype == 3.:
        for deflection_index in range(1, ndelta + 1):
            # DDOC is /FLAPIN/ 1-10.
            fla[25 + deflection_index] = _ix1(_X6226B, _Y6226B,
                                               dsoc(deflection_index) /
                                               f[deflection_index])
            ht[210 + deflection_index] = (fla[25 + deflection_index] *
                                          ht[210 + deflection_index])
    if abs(sweple - swepte) <= 4.0:
        # ARG1 still holds CLDPRM/2 here; the span 4*DELN4 is formed after.
        dumya = _tl(_X1110A, _X2110A, _Y2110A, eta[5], arg1)
        dumyb = _tl(_X1110B, _X2110B, _Y2110B, aw, dumya)
        dumyc = _tl(_X1110C, _X2110C, _Y2110C, taprw, dumyb)
        for deflection_index in range(1, ndelta + 1):
            cnods = _tl(_X1110D, _X2110D, _Y2110D, xsoc(deflection_index),
                        dumyc, u2=1)
            ht[220 + deflection_index] = (cnods * dsoc(deflection_index) *
                                          scale)
    else:
        bs = 4. * deln4
        if bs <= 0.4:
            dumya = _tl(_X111A0, _X211A0, _Y111A0, bs, eta[5])
        if 0.4 < bs <= 0.6:
            dumya = _tl(_X111A1, _X211A1, _Y111A1, bs, eta[5])
        if 0.6 < bs <= 0.8:
            dumya = _tl(_X111A2, _X211A2, _Y111A2, bs, eta[5])
        if bs > 0.8:
            dumya = _ix1(_X111A3, _Y111A3, bs)
        dumyb = _tl(_X1111B, _X2111B, _Y2111B, aw, dumya)
        dumyc = _tl(_X1111C, _X2111C, _Y2111C, taprw, dumyb)
        dumyd = _tl(_X1111D, _X2111D, _Y2111D, sweple, dumyc)
        cnodsb = (_ix1(_X2111E, _Y2111E, dumyd, u1=1) if stype == 1.
                  else dumyd)
        for deflection_index in range(1, ndelta + 1):
            ht[220 + deflection_index] = (cnodsb * dsoc(deflection_index) *
                                          scale)
    if stype == 3.0:
        kssd = _ix1(_X22112, _Y22112, mach * float(s['cosc4']))
        for deflection_index in range(1, ndelta + 1):
            ht[220 + deflection_index] = ht[220 + deflection_index] * kssd


def _tail_roll(data, f, clrol, nalpha, sref, blref, aw):
    """Label 1270: a differentially deflected horizontal tail."""
    h = data['htail']
    gd1, gd2, gd3, gd4 = (float(v) for v in data['gd'])
    dedalp = [0.0] + [float(v) for v in data['dedalp']]
    rivbh = [0.0] + [float(v) for v in data['rivbh']]
    gamvr = [0.0] + [float(v) for v in data['gamvr']]
    bo2h, bo2hst, shst = float(h['bo2h']), float(h['bo2hst']), float(h['sh'])
    arg1 = 0.352 * gd1 + 0.503 * gd2 + 0.344 * gd3 + 0.041 * gd4
    arg2 = 0.383 * gd1 + 0.707 * gd2 + 0.924 * gd3 + 0.500 * gd4
    etacp = arg1 / arg2
    dho2 = bo2h - bo2hst
    eqhoq = _ix1(_X62122, _Y62122, dho2 / bo2h, l1=2, u1=2)
    yh = etacp * bo2hst + dho2
    clahs = float(h['clah']) * sref / shst
    a2 = dho2 / bo2hst
    arg1 = eqhoq * yh * shst * clahs / (blref * sref)
    a1 = PI * aw / RAD
    ndelta = int(f[16] + 0.5)
    for angle_slot in range(1, nalpha + 1):
        b1 = 1. - a1 * dedalp[angle_slot]
        b2 = rivbh[angle_slot] * gamvr[angle_slot] * a2
        cldh = 0.5 * (b1 + b2) * arg1
        for deflection_index in range(1, ndelta + 1):
            nn = 10 * (angle_slot - 1) + deflection_index
            clrol[nn] = (cldh * f[18 + deflection_index] -
                         cldh * f[28 + deflection_index])


def m52o64(data: Mapping[str, object]) -> Dict[str, object]:
    """M52O64: the overlay that sets its number (52) and runs LATFLP."""
    return calculate_latflp(data)
