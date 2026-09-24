"""
Transonic buildup: the wing-body and wing-body-tail routines of overlays 24,
25 and 35 that follow TRSONI and TRSONJ.

- ``TRANCD``: the surface-body drag due to lift, WBCDL's regression on top
  of the transonic zero-lift drag.
- ``TRACM0``: the surface-body zero-lift moment, WBCM0's regression.
- ``TRANCM``/``TRHTCM`` with ``WBTRAN``/``HBTRAN`` and ``WBCM1``: the
  aerodynamic centre faired across the transonic range and the wing-body
  moment slope.
- ``TRAWBT`` (through ``WBTRA``): the wing-body-tail lift and moment slopes.

Each wing routine and its tail twin (TRANCM/TRHTCM, WBTRAN/HBTRAN) is the
same code on different COMMON blocks, so one translation serves both.

The tables were extracted from the source DATA statements by parsing, with
``tools/fortran_data.py``, rather than by hand.

Reference: datcom-legacy/datcom_2000/trancd.f, tracm0.f, trancm.f,
trhtcm.f, wbtran.f, hbtran.f, wbcm1.f, trawbt.f, wbtra.f, m24o30.f,
m25o31.f
"""

import math
from typing import Dict, Mapping, Optional, Sequence

import numpy as np

from pydatcom.aerodynamics.cdrag import STRAIGHT_TAPERED
from pydatcom.aerodynamics.wing_body import (_X21C, _X38B, _Y21C, _Y38B,
                                             calculate_wbcm0,
                                             surface_regression_drag)
from pydatcom.interactions.body_vortex import calculate_bodowg, getmax
from pydatcom.interactions.carryover import (_FIG_431212A_KKBW,
                                             _FIG_431212A_KKWB,
                                             _FIG_431212A_RATIO)
from pydatcom.utils.constants import PI, RAD, UNUSED
from pydatcom.utils.legacy_interp import interx
from pydatcom.utils.legacy_numeric import tbfunx
from pydatcom.utils.legacy_tables import tlin4x, tlinex
from pydatcom.utils.tranac import tranac

# TRANCM: the Mach numbers of the Figure 4.1.4.2-26 lookups and their differences.
_XM = np.array([
    0.5, 0.7, 0.6, 1.3, 1.5, 1.4,
])
# Figure 4.1.4.2-26A-F: abscissa, A*tan(LE) and taper grids (INTERX, six each).
_T422AF = np.array([
    0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.0, 2.0,
    3.0, 4.0, 5.0, 6.0, 0.0, 0.2, 0.25, 0.33,
    0.5, 1.0,
])
_SUBAF1 = np.array([
    0.25, 0.245, 0.24, 0.235, 0.23, 0.225, 0.335, 0.335,
    0.335, 0.335, 0.335, 0.335, 0.42, 0.43, 0.435, 0.445,
    0.45, 0.455, 0.5, 0.515, 0.53, 0.54, 0.55, 0.56,
    0.58, 0.6, 0.63, 0.645, 0.66, 0.68, 0.68, 0.695,
    0.72, 0.74, 0.76, 0.78, 0.285, 0.275, 0.27, 0.265,
    0.26, 0.255, 0.4, 0.41, 0.415, 0.415, 0.415, 0.41,
    0.51, 0.53, 0.535, 0.54, 0.545, 0.55, 0.64, 0.65,
    0.66, 0.675, 0.685, 0.69, 0.75, 0.765, 0.78, 0.785,
    0.8, 0.815, 0.87, 0.88, 0.895, 0.905, 0.92, 0.93,
    0.3, 0.295, 0.285, 0.28, 0.275, 0.265, 0.42, 0.42,
    0.425, 0.425, 0.425, 0.43, 0.545, 0.55, 0.56, 0.565,
    0.575, 0.58, 0.67, 0.68, 0.69, 0.7, 0.71, 0.72,
    0.795, 0.805, 0.815, 0.83, 0.84, 0.85, 0.925, 0.945,
    0.96, 0.965, 0.975, 0.98, 0.325, 0.32, 0.315, 0.305,
    0.3, 0.29, 0.46, 0.46, 0.46, 0.46, 0.455, 0.455,
    0.595, 0.6, 0.6, 0.6, 0.61, 0.62, 0.735, 0.74,
    0.75, 0.76, 0.765, 0.775, 0.885, 0.89, 0.895, 0.9,
    0.91, 0.925, 1.045, 1.05, 1.05, 1.06, 1.065, 1.075,
    0.355, 0.35, 0.345, 0.34, 0.33, 0.32, 0.53, 0.53,
    0.525, 0.525, 0.52, 0.52, 0.7, 0.7, 0.7, 0.705,
    0.71, 0.71, 0.88, 0.88, 0.885, 0.89, 0.89, 0.895,
    1.04, 1.045, 1.05, 1.055, 1.06, 1.065, 1.2, 1.205,
    1.21, 1.215, 1.225, 1.23, 0.51, 0.49, 0.48, 0.47,
    0.46, 0.45, 0.75, 0.75, 0.75, 0.75, 0.745, 0.74,
    1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.25, 1.25,
    1.25, 1.25, 1.25, 1.25, 1.5, 1.5, 1.5, 1.49,
    1.49, 1.49, 1.74, 1.74, 1.74, 1.73, 1.73, 1.73,
])
_SUBAF2 = np.array([
    0.165, 0.18, 0.2, 0.21, 0.22, 0.225, 0.335, 0.335,
    0.335, 0.335, 0.335, 0.335, 0.5, 0.48, 0.465, 0.46,
    0.46, 0.455, 0.67, 0.625, 0.595, 0.58, 0.575, 0.56,
    0.83, 0.75, 0.73, 0.705, 0.695, 0.68, 0.99, 0.86,
    0.835, 0.81, 0.795, 0.78, 0.2, 0.215, 0.23, 0.24,
    0.25, 0.255, 0.4, 0.4, 0.4, 0.405, 0.41, 0.41,
    0.6, 0.58, 0.565, 0.56, 0.555, 0.55, 0.795, 0.76,
    0.735, 0.715, 0.7, 0.69, 0.97, 0.91, 0.87, 0.84,
    0.825, 0.815, 1.15, 1.05, 1.0, 0.965, 0.94, 0.93,
    0.23, 0.24, 0.245, 0.25, 0.26, 0.265, 0.415, 0.42,
    0.425, 0.425, 0.43, 0.43, 0.63, 0.615, 0.6, 0.59,
    0.585, 0.58, 0.83, 0.785, 0.76, 0.74, 0.73, 0.72,
    1.03, 0.95, 0.905, 0.88, 0.865, 0.85, 1.25, 1.09,
    1.05, 1.015, 0.99, 0.98, 0.22, 0.24, 0.25, 0.265,
    0.28, 0.29, 0.44, 0.445, 0.45, 0.45, 0.455, 0.455,
    0.67, 0.655, 0.64, 0.63, 0.625, 0.62, 0.88, 0.83,
    0.805, 0.79, 0.78, 0.775, 1.07, 1.0, 0.96, 0.94,
    0.935, 0.925, 1.27, 1.17, 1.12, 1.1, 1.085, 1.075,
    0.25, 0.27, 0.295, 0.31, 0.315, 0.32, 0.5, 0.505,
    0.51, 0.515, 0.52, 0.52, 0.75, 0.74, 0.73, 0.72,
    0.715, 0.71, 0.98, 0.94, 0.915, 0.9, 0.9, 0.895,
    1.19, 1.12, 1.09, 1.08, 1.07, 1.065, 1.38, 1.3,
    1.27, 1.25, 1.24, 1.23, 0.34, 0.38, 0.41, 0.43,
    0.44, 0.45, 0.68, 0.7, 0.72, 0.73, 0.73, 0.74,
    0.95, 0.98, 1.0, 1.0, 1.0, 1.0, 1.2, 1.23,
    1.25, 1.25, 1.25, 1.25, 1.44, 1.47, 1.48, 1.48,
    1.49, 1.49, 1.68, 1.71, 1.71, 1.72, 1.72, 1.73,
])
_SUPAF1 = np.array([
    0.165, 0.21, 0.25, 0.29, 0.31, 0.345, 0.335, 0.365,
    0.39, 0.415, 0.445, 0.47, 0.5, 0.54, 0.56, 0.56,
    0.56, 0.56, 0.67, 0.67, 0.67, 0.67, 0.67, 0.67,
    0.83, 0.775, 0.775, 0.775, 0.775, 0.775, 0.99, 0.93,
    0.895, 0.895, 0.895, 0.895, 0.2, 0.23, 0.28, 0.305,
    0.335, 0.36, 0.4, 0.445, 0.485, 0.5, 0.52, 0.53,
    0.6, 0.63, 0.65, 0.66, 0.665, 0.665, 0.795, 0.8,
    0.8, 0.805, 0.81, 0.815, 0.97, 0.965, 0.955, 0.955,
    0.955, 0.955, 1.15, 1.135, 1.12, 1.1, 1.1, 1.105,
    0.23, 0.275, 0.3, 0.33, 0.35, 0.37, 0.415, 0.47,
    0.5, 0.53, 0.545, 0.55, 0.63, 0.67, 0.68, 0.685,
    0.69, 0.69, 0.83, 0.835, 0.835, 0.84, 0.845, 0.85,
    1.03, 1.015, 1.005, 1.0, 1.005, 1.01, 1.25, 1.225,
    1.2, 1.17, 1.165, 1.16, 0.22, 0.28, 0.315, 0.345,
    0.375, 0.39, 0.44, 0.5, 0.535, 0.56, 0.57, 0.58,
    0.67, 0.7, 0.72, 0.725, 0.74, 0.74, 0.88, 0.885,
    0.895, 0.9, 0.9, 0.9, 1.07, 1.07, 1.075, 1.075,
    1.08, 1.08, 1.27, 1.26, 1.26, 1.255, 1.255, 1.255,
    0.25, 0.3, 0.33, 0.38, 0.415, 0.445, 0.5, 0.56,
    0.6, 0.62, 0.635, 0.64, 0.75, 0.78, 0.8, 0.82,
    0.82, 0.825, 0.98, 0.99, 1.0, 1.02, 1.02, 1.02,
    1.19, 1.2, 1.2, 1.21, 1.22, 1.225, 1.38, 1.39,
    1.4, 1.41, 1.42, 1.42, 0.34, 0.38, 0.41, 0.46,
    0.5, 0.54, 0.68, 0.7, 0.73, 0.77, 0.79, 0.84,
    0.95, 0.99, 1.01, 1.05, 1.08, 1.12, 1.2, 1.24,
    1.29, 1.33, 1.37, 1.42, 1.44, 1.5, 1.55, 1.61,
    1.67, 1.72, 1.68, 1.76, 1.82, 1.89, 1.95, 2.02,
])
_SUPAF2 = np.array([
    0.415, 0.41, 0.4, 0.385, 0.37, 0.345, 0.5, 0.5,
    0.495, 0.485, 0.48, 0.47, 0.585, 0.58, 0.58, 0.575,
    0.57, 0.56, 0.67, 0.67, 0.67, 0.67, 0.67, 0.67,
    0.75, 0.75, 0.755, 0.76, 0.765, 0.775, 0.83, 0.84,
    0.845, 0.855, 0.87, 0.895, 0.46, 0.455, 0.445, 0.42,
    0.39, 0.36, 0.575, 0.575, 0.57, 0.56, 0.545, 0.53,
    0.695, 0.695, 0.69, 0.685, 0.68, 0.665, 0.8, 0.805,
    0.805, 0.81, 0.815, 0.815, 0.92, 0.93, 0.935, 0.945,
    0.97, 0.955, 1.04, 1.045, 1.05, 1.075, 1.11, 1.105,
    0.475, 0.465, 0.45, 0.43, 0.4, 0.37, 0.6, 0.6,
    0.595, 0.585, 0.575, 0.55, 0.725, 0.73, 0.73, 0.725,
    0.715, 0.69, 0.85, 0.85, 0.855, 0.865, 0.87, 0.85,
    0.97, 0.975, 0.98, 1.0, 1.02, 1.01, 1.11, 1.11,
    1.11, 1.13, 1.18, 1.16, 0.5, 0.49, 0.47, 0.45,
    0.425, 0.39, 0.64, 0.635, 0.63, 0.62, 0.6, 0.58,
    0.77, 0.775, 0.78, 0.775, 0.765, 0.74, 0.92, 0.915,
    0.92, 0.93, 0.935, 0.9, 1.05, 1.055, 1.06, 1.08,
    1.105, 1.08, 1.195, 1.2, 1.205, 1.225, 1.265, 1.255,
    0.55, 0.535, 0.525, 0.5, 0.475, 0.445, 0.72, 0.715,
    0.71, 0.69, 0.67, 0.64, 0.89, 0.89, 0.89, 0.885,
    0.87, 0.825, 1.06, 1.05, 1.05, 1.06, 1.065, 1.02,
    1.215, 1.215, 1.22, 1.245, 1.27, 1.225, 1.38, 1.38,
    1.395, 1.42, 1.47, 1.42, 0.76, 0.73, 0.7, 0.65,
    0.6, 0.54, 1.0, 1.0, 0.97, 0.93, 0.89, 0.84,
    1.24, 1.23, 1.23, 1.22, 1.19, 1.12, 1.5, 1.48,
    1.48, 1.49, 1.47, 1.42, 1.75, 1.72, 1.73, 1.76,
    1.78, 1.72, 2.0, 1.97, 1.98, 2.02, 2.07, 2.02,
])
# Figure 4.1.4.2-30A-D: A*t^(1/3) (7), A*tan(LE) (9), transonic similarity (4) and taper (3) grids, padded as in the source.
_T425AD = np.array([
    0.0, 0.75, 1.0, 1.25, 2.0, 3.0, 6.0, 0.0,
    0.0, 0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0,
    10.0, 15.0, -2.0000001, -1.0, 0.0, 1.0000001, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.2, 0.5,
])
_D425AD = np.array([
    0.0, 0.1, 0.125, 0.14, 0.16, 0.173, 0.21, 0.175,
    0.23, 0.26, 0.264, 0.28, 0.3, 0.37, 0.34, 0.375,
    0.39, 0.392, 0.4, 0.412, 0.45, 0.496, 0.504, 0.509,
    0.51, 0.518, 0.53, 0.565, 0.655, 0.625, 0.623, 0.623,
    0.625, 0.627, 0.63, 0.81, 0.735, 0.735, 0.732, 0.729,
    0.725, 0.705, 0.97, 0.855, 0.845, 0.835, 0.835, 0.825,
    0.805, 1.61, 1.28, 1.26, 1.245, 1.215, 1.17, 1.035,
    2.405, 1.8, 1.725, 1.705, 1.635, 1.555, 1.28, 0.0,
    0.09, 0.115, 0.14, 0.168, 0.2, 0.29, 0.173, 0.245,
    0.27, 0.288, 0.305, 0.32, 0.365, 0.34, 0.39, 0.405,
    0.425, 0.435, 0.445, 0.49, 0.5, 0.525, 0.53, 0.533,
    0.54, 0.553, 0.585, 0.664, 0.65, 0.652, 0.654, 0.66,
    0.667, 0.685, 0.815, 0.78, 0.776, 0.776, 0.778, 0.78,
    0.785, 0.97, 0.9, 0.88, 0.88, 0.88, 0.88, 0.88,
    1.61, 1.4, 1.325, 1.27, 1.28, 1.29, 1.31, 2.41,
    1.98, 1.84, 1.72, 1.745, 1.765, 1.84, 0.0, 0.11,
    0.15, 0.199, 0.199, 0.199, 0.199, 0.18, 0.26, 0.285,
    0.31, 0.321, 0.321, 0.321, 0.35, 0.405, 0.427, 0.44,
    0.445, 0.445, 0.445, 0.5, 0.533, 0.543, 0.55, 0.55,
    0.55, 0.55, 0.665, 0.663, 0.662, 0.661, 0.66, 0.66,
    0.66, 0.818, 0.793, 0.788, 0.78, 0.78, 0.78, 0.78,
    0.975, 0.923, 0.905, 0.89, 0.89, 0.89, 0.89, 1.615,
    1.43, 1.38, 1.32, 1.32, 1.32, 1.32, 2.415, 2.08,
    1.96, 1.845, 1.845, 1.845, 1.845, 0.0, 0.12, 0.165,
    0.2, 0.2, 0.2, 0.2, 0.175, 0.275, 0.3, 0.324,
    0.324, 0.324, 0.324, 0.335, 0.417, 0.445, 0.445, 0.445,
    0.445, 0.445, 0.5, 0.548, 0.548, 0.548, 0.548, 0.548,
    0.548, 0.665, 0.665, 0.665, 0.665, 0.665, 0.665, 0.665,
    0.815, 0.775, 0.775, 0.775, 0.775, 0.775, 0.775, 0.985,
    0.88, 0.88, 0.88, 0.88, 0.88, 0.88, 1.63, 1.28,
    1.28, 1.28, 1.28, 1.28, 1.28, 2.435, 1.73, 1.73,
    1.73, 1.73, 1.73, 1.73, 0.0, 0.12, 0.135, 0.15,
    0.168, 0.186, 0.24, 0.2, 0.27, 0.28, 0.284, 0.296,
    0.315, 0.375, 0.4, 0.425, 0.426, 0.43, 0.44, 0.45,
    0.483, 0.6, 0.575, 0.575, 0.575, 0.576, 0.577, 0.589,
    0.795, 0.73, 0.718, 0.719, 0.72, 0.721, 0.725, 0.98,
    0.86, 0.837, 0.835, 0.83, 0.825, 0.805, 1.175, 1.02,
    0.98, 0.975, 0.965, 0.955, 0.91, 1.955, 1.59, 1.47,
    1.405, 1.4, 1.39, 1.37, 2.93, 2.36, 2.16, 1.96,
    1.94, 1.91, 1.83, 0.0, 0.085, 0.12, 0.145, 0.186,
    0.215, 0.295, 0.205, 0.27, 0.3, 0.325, 0.348, 0.368,
    0.43, 0.4, 0.449, 0.47, 0.485, 0.495, 0.51, 0.53,
    0.6, 0.61, 0.614, 0.618, 0.632, 0.649, 0.7, 0.8,
    0.778, 0.77, 0.773, 0.78, 0.794, 0.83, 0.978, 0.935,
    0.915, 0.916, 0.919, 0.919, 0.919, 1.165, 1.085, 1.06,
    1.06, 1.055, 1.05, 1.045, 1.94, 1.75, 1.695, 1.645,
    1.62, 1.585, 1.48, 2.905, 2.505, 2.485, 2.45, 2.38,
    2.275, 1.965, 0.0, 0.125, 0.175, 0.21, 0.215, 0.215,
    0.215, 0.2, 0.29, 0.32, 0.355, 0.36, 0.36, 0.36,
    0.4, 0.465, 0.485, 0.505, 0.51, 0.51, 0.51, 0.6,
    0.63, 0.64, 0.648, 0.65, 0.65, 0.65, 0.8, 0.795,
    0.794, 0.793, 0.793, 0.793, 0.793, 0.98, 0.94, 0.925,
    0.91, 0.91, 0.91, 0.91, 1.175, 1.07, 1.04, 1.04,
    1.04, 1.04, 1.04, 1.96, 1.61, 1.49, 1.49, 1.49,
    1.49, 1.49, 2.94, 2.24, 2.0, 2.0, 2.0, 2.0,
    2.0, 0.0, 0.15, 0.2, 0.237, 0.237, 0.237, 0.237,
    0.2, 0.32, 0.37, 0.382, 0.382, 0.382, 0.382, 0.4,
    0.49, 0.523, 0.523, 0.523, 0.523, 0.523, 0.6, 0.64,
    0.65, 0.65, 0.65, 0.65, 0.65, 0.8, 0.782, 0.78,
    0.78, 0.78, 0.78, 0.78, 0.98, 0.895, 0.888, 0.888,
    0.888, 0.888, 0.888, 1.175, 1.04, 1.0, 1.0, 1.0,
    1.0, 1.0, 1.96, 1.54, 1.42, 1.42, 1.42, 1.42,
    1.42, 2.94, 1.91, 1.91, 1.91, 1.91, 1.91, 1.91,
    0.0, 0.12, 0.15, 0.163, 0.178, 0.195, 0.255, 0.26,
    0.335, 0.35, 0.351, 0.365, 0.385, 0.43, 0.5, 0.53,
    0.535, 0.54, 0.55, 0.565, 0.61, 0.753, 0.738, 0.738,
    0.738, 0.735, 0.735, 0.725, 0.975, 0.912, 0.908, 0.904,
    0.894, 0.884, 0.845, 1.22, 1.14, 1.11, 1.08, 1.06,
    1.04, 0.96, 1.46, 1.33, 1.28, 1.24, 1.21, 1.18,
    1.09, 2.42, 2.07, 1.96, 1.84, 1.785, 1.725, 1.52,
    3.6, 3.18, 3.04, 2.9, 2.48, 2.365, 2.0, 0.0,
    0.095, 0.14, 0.173, 0.23, 0.278, 0.42, 0.26, 0.34,
    0.37, 0.394, 0.42, 0.455, 0.55, 0.51, 0.56, 0.585,
    0.593, 0.61, 0.634, 0.705, 0.75, 0.775, 0.78, 0.786,
    0.798, 0.81, 0.85, 0.98, 0.955, 0.95, 0.95, 0.948,
    0.945, 0.935, 1.18, 1.105, 1.097, 1.094, 1.08, 1.065,
    1.015, 1.385, 1.27, 1.24, 1.23, 1.22, 1.205, 1.17,
    2.15, 1.86, 1.76, 1.65, 1.64, 1.62, 1.56, 3.05,
    2.5, 2.31, 2.14, 2.115, 2.085, 2.0, 0.0, 0.15,
    0.2, 0.26, 0.275, 0.275, 0.275, 0.26, 0.37, 0.405,
    0.445, 0.455, 0.455, 0.455, 0.51, 0.583, 0.605, 0.63,
    0.634, 0.634, 0.634, 0.755, 0.788, 0.8, 0.8, 0.8,
    0.8, 0.8, 0.98, 0.96, 0.95, 0.943, 0.943, 0.943,
    0.943, 1.175, 1.11, 1.093, 1.064, 1.064, 1.064, 1.064,
    1.36, 1.24, 1.21, 1.163, 1.163, 1.163, 1.163, 2.07,
    1.75, 1.645, 1.54, 1.54, 1.54, 1.54, 2.91, 2.32,
    2.125, 1.94, 1.94, 1.94, 1.94, 0.0, 0.165, 0.245,
    0.28, 0.28, 0.28, 0.28, 0.26, 0.405, 0.475, 0.475,
    0.475, 0.475, 0.475, 0.51, 0.61, 0.65, 0.65, 0.65,
    0.65, 0.65, 0.755, 0.8, 0.82, 0.82, 0.82, 0.82,
    0.82, 0.98, 0.961, 0.96, 0.96, 0.96, 0.96, 0.96,
    1.18, 1.08, 1.08, 1.08, 1.08, 1.08, 1.08, 1.37,
    1.25, 1.25, 1.25, 1.25, 1.25, 1.25, 2.09, 1.78,
    1.78, 1.78, 1.78, 1.78, 1.78, 2.91, 2.4, 2.4,
    2.4, 2.4, 2.4, 2.4,
])
# Figure 4.1.4.2-33: the thickness increment.
_T428 = np.array([
    0.0, 6.6, 10.0, 14.0, 0.0, 0.0, 0.0, 0.0,
    1.7, 2.0, 2.5, 3.0, 4.0, 6.0,
])
_D428 = np.array([
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, -0.1, -0.22, 0.0, 0.0, -0.21, -0.45,
    0.0, 0.0, -0.305, -0.665, 0.0, 0.0, -0.475, -1.04,
    0.0, 0.0, -0.65, -1.41,
])
# WBTRAN: Figure 4.3.1.2-10, K_W(B) and K_B(W) against d/b.
_TFIG10 = np.array([
    0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7,
    0.8, 0.9, 1.0,
])
_DKWB10 = np.array([
    1.0, 1.08, 1.16, 1.26, 1.36, 1.46, 1.56, 1.67,
    1.78, 1.89, 2.0,
])
_DKBW10 = np.array([
    0.0, 0.13, 0.29, 0.45, 0.62, 0.8, 1.0, 1.22,
    1.45, 1.7, 2.0,
])
# Figure 4.3.1.2-11A (afterbody) and -11B (no afterbody): K_B(W) at supersonic speed.
_T4311A = np.array([
    0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4,
    1.6, 1.8, 2.0, 2.4, 2.8, 3.2, 4.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.1, 0.2, 0.3, 0.4,
    0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.2, 1.6,
    2.0, 3.0, 4.0, 8.0, 10.0, 999999.0,
])
_D4311A = np.array([
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.75,
    0.58, 0.47, 0.38, 0.33, 0.32, 0.3, 0.3, 0.29,
    0.28, 0.27, 0.25, 0.23, 0.21, 0.19, 1.45, 1.2,
    1.0, 0.87, 0.78, 0.75, 0.7, 0.68, 0.64, 0.6,
    0.58, 0.55, 0.5, 0.48, 0.4, 1.95, 1.7, 1.48,
    1.33, 1.2, 1.12, 1.07, 1.0, 0.97, 0.92, 0.88,
    0.82, 0.75, 0.72, 0.64, 2.5, 2.2, 1.93, 1.74,
    1.58, 1.5, 1.42, 1.38, 1.3, 1.25, 1.2, 1.1,
    1.0, 0.93, 0.8, 2.95, 2.6, 2.25, 2.05, 1.89,
    1.8, 1.7, 1.65, 1.57, 1.5, 1.48, 1.36, 1.27,
    1.2, 1.03, 3.3, 2.92, 2.6, 2.35, 2.18, 2.04,
    1.95, 1.86, 1.8, 1.73, 1.7, 1.59, 1.49, 1.4,
    1.21, 3.7, 3.32, 2.99, 2.7, 2.5, 2.33, 2.23,
    2.12, 2.03, 1.97, 1.9, 1.8, 1.7, 1.6, 1.42,
    4.1, 3.64, 3.25, 2.97, 2.77, 2.6, 2.45, 2.32,
    2.23, 2.15, 2.1, 1.95, 1.85, 1.76, 1.62, 4.3,
    3.83, 3.49, 3.2, 2.99, 2.8, 2.65, 2.5, 2.42,
    2.33, 2.29, 2.14, 2.0, 1.89, 1.74, 4.6, 4.08,
    3.69, 3.35, 3.14, 2.95, 2.8, 2.69, 2.6, 2.49,
    2.4, 2.25, 2.16, 2.02, 1.85, 5.0, 4.37, 3.98,
    3.64, 3.39, 3.19, 3.03, 2.89, 2.8, 2.7, 2.61,
    2.43, 2.33, 2.21, 2.0, 5.45, 4.75, 4.29, 3.95,
    3.7, 3.5, 3.34, 3.2, 3.1, 3.0, 2.91, 2.75,
    2.6, 2.48, 2.28, 5.9, 5.1, 4.63, 4.25, 4.0,
    3.8, 3.65, 3.5, 3.33, 3.23, 3.14, 2.95, 2.8,
    2.67, 2.4, 6.2, 5.5, 5.0, 4.65, 4.35, 4.15,
    3.98, 3.8, 3.67, 3.53, 3.4, 3.24, 3.1, 2.92,
    2.68, 6.6, 5.75, 5.25, 4.9, 4.6, 4.4, 4.22,
    4.05, 3.91, 3.78, 3.68, 3.45, 3.26, 3.1, 2.83,
    7.05, 6.15, 5.6, 5.2, 4.95, 4.7, 4.49, 4.3,
    4.13, 4.0, 3.85, 3.68, 3.5, 3.34, 3.1, 7.4,
    6.35, 5.8, 5.4, 5.12, 4.85, 4.67, 4.49, 4.33,
    4.2, 4.07, 3.85, 3.65, 3.5, 3.24, 8.0, 6.73,
    6.2, 5.8, 5.5, 5.25, 5.0, 4.82, 4.67, 4.5,
    4.35, 4.12, 3.9, 3.8, 3.59,
])
_T4311B = np.array([
    0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4,
    1.6, 1.8, 2.0, 2.4, 2.8, 3.2, 4.0, 0.0,
    0.2, 0.4, 0.6, 0.8, 1.0, 2.0, 4.0, 999999.0,
])
_D4311B = np.array([
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.2,
    1.0, 0.75, 0.58, 0.45, 0.37, 0.3, 0.28, 0.24,
    0.2, 0.18, 0.15, 0.13, 0.11, 0.1, 2.4, 1.85,
    1.42, 1.1, 0.92, 0.77, 0.66, 0.6, 0.51, 0.48,
    0.41, 0.31, 0.26, 0.22, 0.22, 3.5, 2.6, 2.0,
    1.58, 1.28, 1.07, 0.91, 0.78, 0.7, 0.6, 0.52,
    0.42, 0.39, 0.33, 0.3, 4.3, 3.1, 2.35, 1.8,
    1.45, 1.2, 1.04, 0.92, 0.83, 0.75, 0.68, 0.55,
    0.48, 0.4, 0.38, 5.0, 3.65, 2.78, 2.22, 1.82,
    1.5, 1.3, 1.15, 1.0, 0.92, 0.82, 0.69, 0.58,
    0.5, 0.42, 5.75, 4.55, 3.66, 2.9, 2.3, 1.92,
    1.62, 1.42, 1.25, 1.1, 0.97, 0.82, 0.7, 0.63,
    0.52, 6.7, 5.25, 4.18, 3.32, 2.63, 2.2, 1.88,
    1.62, 1.45, 1.28, 1.15, 1.0, 0.85, 0.75, 0.62,
    7.6, 6.2, 4.91, 3.95, 3.18, 2.6, 2.22, 1.91,
    1.7, 1.5, 1.35, 1.12, 0.97, 0.85, 0.7,
])
# Figure 4.3.2.2-37A/B: the carryover aerodynamic centre at Mach 1.4.
_T4337A = np.array([
    0.0, 0.4, 0.8, 1.2, 1.6, 2.0, 2.4, 2.8,
    0.1, 1.0, 999999.0,
])
_D4337A = np.array([
    0.5, 0.72, 0.9, 1.08, 1.24, 1.39, 1.53, 1.68,
    0.5, 0.72, 0.91, 1.09, 1.25, 1.41, 1.57, 1.72,
    0.5, 0.73, 0.92, 1.11, 1.27, 1.43, 1.59, 1.74,
])
_T4337B = np.array([
    0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.8,
    1.0, 2.8, 0.2, 999999.0,
])
_D4337B = np.array([
    0.5, 0.56, 0.595, 0.62, 0.64, 0.65, 0.66, 0.669,
    0.669, 0.671, 0.5, 0.54, 0.578, 0.6, 0.62, 0.638,
    0.649, 0.66, 0.669, 0.671,
])
# TRAWBT: Figure 4.4.1-66, effective aspect ratio and span of the wing wake (three parts).
_X155A1 = np.array([
    60.0, 45.0, 30.0, 0.0,
])
_X155A2 = np.array([
    0.4, 0.55, 0.725, 1.0,
])
_Y4155A = np.array([
    5.0, 4.0, 3.0, 1.25, 5.0, 5.0, 4.2, 2.75,
    5.0, 5.0, 5.0, 4.0, 5.0, 5.0, 5.0, 5.0,
])
_X155B1 = np.array([
    0.0, 0.5, 1.0,
])
_X155B2 = np.array([
    5.0, 0.3,
])
_Y4155B = np.array([
    1.0, 0.3, 1.0, 0.425, 1.0, 0.525,
])
_X155C1 = np.array([
    0.0, 1.0,
])
_X155C2 = np.array([
    0.41, 0.5, 0.58, 0.66, 0.75, 0.9, 1.0,
])
_Y4155C = np.array([
    0.56, 0.66, 0.725, 0.8, 0.86, 0.94, 0.99, 0.58,
    0.58, 0.58, 0.675, 0.74, 0.9, 0.99,
])
# Figure 4.4.1-67: the downwash gradient in two parts.
_X156A1 = np.array([
    1.5, 2.0, 3.0, 5.0, 7.0, 10.0,
])
_X156A2 = np.array([
    0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9,
    1.0, 1.2, 1.4, 1.8, 2.2,
])
_Y156A = np.array([
    3.4, 2.88, 2.45, 2.1, 1.8, 1.6, 1.4, 1.28,
    1.12, 0.9, 0.76, 0.58, 0.45, 3.2, 2.8, 2.35,
    1.93, 1.65, 1.46, 1.26, 1.12, 1.0, 0.8, 0.64,
    0.46, 0.35, 2.75, 2.33, 1.95, 1.6, 1.34, 1.2,
    1.0, 0.88, 0.75, 0.6, 0.47, 0.32, 0.26, 2.25,
    1.85, 1.46, 1.2, 0.98, 0.8, 0.7, 0.6, 0.55,
    0.4, 0.35, 0.25, 0.2, 1.95, 1.54, 1.23, 1.0,
    0.8, 0.68, 0.58, 0.5, 0.42, 0.33, 0.26, 0.18,
    0.18, 1.55, 1.2, 0.9, 0.74, 0.6, 0.5, 0.4,
    0.35, 0.3, 0.22, 0.18, 0.1, 0.1,
])
_X156B1 = np.array([
    0.0, 30.0, 45.0, 60.0,
])
_X156B2 = np.array([
    1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 6.0, 8.0,
    10.0,
])
_Y4156B = np.array([
    1.25, 0.92, 0.8, 0.72, 0.64, 0.52, 0.4, 0.33,
    0.3, 1.0, 0.88, 0.76, 0.68, 0.6, 0.49, 0.36,
    0.296, 0.26, 0.9, 0.76, 0.66, 0.572, 0.52, 0.432,
    0.32, 0.252, 0.2, 0.8, 0.64, 0.52, 0.464, 0.4,
    0.33, 0.24, 0.18, 0.16,
])
# Figure 4.4.1-68B: the tail-height correction.
_X157B1 = np.array([
    0.0, 0.2, 0.6, 1.0,
])
_X157B2 = np.array([
    0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.833,
    1.0, 1.1, 1.2, 1.3, 1.4, 1.5,
])
_Y4157B = np.array([
    1.0, 1.03, 1.068, 1.11, 1.176, 1.26, 1.36, 1.34,
    1.2, 1.11, 1.0, 0.88, 0.74, 0.6, 0.96, 0.96,
    0.965, 0.976, 0.98, 1.016, 1.036, 1.04, 0.946, 0.88,
    0.8, 0.69, 0.58, 0.46, 0.74, 0.73, 0.72, 0.71,
    0.69, 0.67, 0.644, 0.638, 0.582, 0.56, 0.52, 0.48,
    0.44, 0.4, 0.5, 0.5, 0.5, 0.49, 0.476, 0.46,
    0.448, 0.44, 0.4, 0.372, 0.352, 0.32, 0.3, 0.27,
])
# TLIN4X's Y(NX2,NX1,NX3,NX4): A*t^(1/3) fastest, then A*tan(LE),
# similarity and taper.
_Y425AD = _D425AD.reshape(7, 9, 4, 3, order='F')


def _tlinex(x1, x2, flat, q1, q2, l1, l2, u1, u2) -> float:
    """TLINEX on a source ``Y(NX2,NX1)`` table given flat."""
    y = np.asarray(flat, dtype=float).reshape(len(x1), len(x2)).T
    return float(tlinex(x1, x2, y, q1, q2, l1, l2, u1, u2))


def _stra(value) -> bool:
    return float(value) == STRAIGHT_TAPERED


def calculate_trancd(alpha_deg: Sequence[float], surface: Dict[str, float],
                     cd0: float, body: Dict[str, float], x_surface: float,
                     flight: Dict[str, float]) -> Optional[np.ndarray]:
    """Translate one half of TRANCD: transonic surface-body drag.

    The wing half runs when the wing is straight tapered and both the
    wing and the body are present; the tail half when the tail is straight
    tapered and present (the source does not test for the body there).
    The caller applies those gates and the ``NF < 0`` skip.

    Args:
        cd0: The surface-body zero-lift drag, ``BW(1)`` or ``BH(1)``, as
            TRSONI or TRSONJ left it.
        Others: as :func:`wing_body.calculate_wbcd`.

    Returns:
        ``CD0 + CDL`` at each angle, ``-UNUSED`` where TABLES has no data;
        ``None`` where the source leaves the drag unchanged.
    """
    if not _stra(surface['type']):
        return None
    regression = surface_regression_drag(alpha_deg, surface, body, x_surface,
                                         flight)
    if regression is None:
        return None
    cdl = regression[0]
    return np.where(cdl == UNUSED, -UNUSED, float(cd0) + cdl)


def calculate_tracm0(surface: Dict[str, float], x_surface: float,
                     wing_height: float, body_length: float,
                     max_diameter: float, mach: float, tr: float,
                     stale_cm0: float = 0.0) -> Dict[str, object]:
    """Translate one half of TRACM0: the transonic surface-body CM0.

    WBCM0's regression with the surface's geometry.  Differences from the
    subsonic call in WBCM, all kept:

    - ``RN = FLC(I+2)*A(122)`` multiplies the MAC by the *Mach number*
      (WBCM uses ``FLC(II+42)``, the unit Reynolds number), so the product
      is far below WBCM0's lower clamp and the regression always sees
      ``RN = 8e5``.
    - ``WL = 0.5*ZW/DBM`` where WBCM forms ``0.5+ZW/DBM``, ``HD`` is 0
      where WBCM uses 0.5, and the diameter is ``SBD(120)``.
    - The tail half reuses the *wing's* height ``ZW``.

    Args:
        surface: ``sspn``, ``sspne``, ``chrdr``, ``tovc``, ``ler``,
            ``twista``, ``ycm``, ``cld`` and ``a38``, ``a80``, ``a118``,
            ``a120``, ``a122``.
        x_surface: ``XW`` or ``XH``.
        wing_height: ``ZW``, for both halves.
        body_length: ``BD(1)``.
        max_diameter: ``SBD(120)``.
        mach: ``FLC(I+2)``.
        tr: ``FLC(96)``.
        stale_cm0: ``BW(41)`` or ``BH(41)`` before the call, which the
            source keeps when the regression is out of range.

    Returns:
        ``cm0`` (``BW(41)``/``BH(41)``, and ``TRA(74)``/``TRAH(74)``),
        ``computed`` and the regression inputs formed.
    """
    surf = {k: float(v) for k, v in surface.items()}
    exposed_diameter = 2.0 * (surf['sspn'] - surf['sspne'])
    xln = (x_surface + 0.5 * exposed_diameter * surf['a38']) / exposed_diameter
    xla = (x_surface + surf['chrdr'] + 0.5 * exposed_diameter * surf['a80'])
    xla = xla / exposed_diameter
    xla = body_length / exposed_diameter - xla
    rn = mach * surf['a122']
    wl = 0.5 * wing_height / max_diameter
    dob = 0.5 * max_diameter / surf['sspn']
    cm0 = calculate_wbcm0(surf['a120'], surf['a38'], surf['tovc'], xln, xla,
                          surf['a118'], surf['ler'], surf['twista'], surf['ycm'],
                          surf['cld'], rn, tr, wl, 0.0, 0.0, dob, mach)
    return {'cm0': float(stale_cm0 if cm0 is None else cm0),
            'computed': cm0 is not None, 'nose_length': xln,
            'afterbody_length': xla, 'reynolds': rn, 'wing_height': wl,
            'body_ratio': dob}


def calculate_wbcm1(a: Mapping[int, float], sspn: float, bd87: float,
                    cbarr: float, stale_wb13: float = 0.0,
                    stale_wb14: float = 0.0) -> Dict[str, object]:
    """Translate WBCM1: the carryover aerodynamic centre at Mach 0.6.

    WBCM's Section 4.3.2 method at ``beta = 0.8``: Figure 4.3.2.2-36B at
    ``beta*A = 0``, the Figure 4.3.2.2-35 value at ``beta*A = 4``, and the
    ellipse between them.

    Args:
        a: The surface's ``A`` block: 7, 10, 27, 38, 44.
        sspn: ``WINGIN(4)``.
        bd87: ``BD(87)``, the body diameter at the surface.
        cbarr: ``CBARR``.
        stale_wb13, stale_wb14: ``WB(13)``, ``WB(14)`` before the call.
            When the ellipse has no real root the source prints a message
            and returns with them unchanged.

    Returns:
        ``wb13`` ((x_ac/c)_B(W) on the reference chord), ``wb14`` (on the
        root chord), ``wb15``, and ``ellipse_failed``.
    """
    a_block = {int(k): float(v) for k, v in a.items()}
    aspect_taper_term = 0.25 * a_block[7] * (1.0 + a_block[27]) * a_block[38]
    wb15 = 0.50
    if aspect_taper_term < 1.0:
        wb15, _ = tbfunx(_X38B, _Y38B, aspect_taper_term, 0, 0)
    brac, _ = tbfunx(_X21C, _Y21C, bd87 / (2.0 * sspn), 0, 0)
    ac_offset = (0.25 + (2.0 * sspn - bd87) / (2.0 * a_block[10]) *
                 a_block[44] * brac)
    beta_ar = 0.80 * a_block[7]
    result = {'wb15': float(wb15), 'temp4': float(ac_offset),
              'ellipse_failed': False}
    if beta_ar >= 4.0:
        wb14 = ac_offset
    else:
        ellipse_y1 = wb15
        ellipse_radius = abs(ac_offset - wb15)
        quad_b = -2.0 * ellipse_y1
        quad_c = (ellipse_y1 * ellipse_y1 - ellipse_radius * ellipse_radius +
                  ((ellipse_radius / 4.0)**2) * (beta_ar - 4.0)**2)
        discriminant = quad_b * quad_b - 4.0 * quad_c
        if not discriminant > 0.0:
            result.update({'wb13': float(stale_wb13),
                           'wb14': float(stale_wb14),
                           'ellipse_failed': True})
            return result
        ymax = (-quad_b + math.sqrt(discriminant)) / 2.0
        ymin = (-quad_b - math.sqrt(discriminant)) / 2.0
        wb14 = ymin if ac_offset < ellipse_y1 else ymax
    result.update({'wb14': float(wb14),
                   'wb13': float(wb14 * a_block[10] / cbarr)})
    return result


def calculate_wbtran(mach: float, alpha_deg: Sequence[float],
                     surface: Dict[str, float], a: Mapping[int, float],
                     x_surface: float, body_length: float, alpha0_body: float,
                     cla_surface: float, cla_body: float, sref: float,
                     cbarr: float) -> Dict[str, object]:
    """Translate WBTRAN (and HBTRAN): the Mach 1.4 surface-body carryover.

    HBTRAN is the same code on the tail's blocks.  K_B(W) from Figure
    4.3.1.2-11A/B where the similarity parameter exceeds four (or, untapered,
    where ``beta*A > 1``), else from Figure 4.3.1.2-10; K_W(B) from Figure
    4.3.1.2-10; and the carryover aerodynamic centre of Figure 4.3.2.2-37A/B,
    both read at ``beta = 0.98`` whatever the Mach number.

    Args:
        mach: ``FLC(I+2)``.
        alpha_deg: ``FLC(23)`` onward.
        surface: ``sspn`` (``WINGIN(4)``), ``sspne`` (3), ``chrdr`` (6).
        a: The surface's ``A`` block: 3, 7, 10, 27, 62.
        x_surface: ``XW`` (``XH`` for HBTRAN).
        body_length: ``XCOOR(NX)``, the last body station.
        alpha0_body: ``BD(81)``; ``UNUSED`` counts as zero.
        cla_surface: ``WING(101)`` (``HT(101)``).
        cla_body: ``BODY(101)``.
        sref, cbarr: ``/OPTION/``.

    Returns:
        ``kbw``, ``kwb``, ``clawb`` (``TRA(71)``), ``clabw`` (``TRA(72)``),
        ``cla`` (``BW(101)``), ``xaca`` (``SWB(39)``), ``xacbw``
        (``SWB(8)``), ``dd`` (``SWB(5)``), ``alpha_body`` (``BD(255)``
        onward), and where reached ``rkbw`` (``SWB(32)``) and ``trino``
        (``SWB(60)``).
    """
    a_block = {int(k): float(v) for k, v in a.items()}
    span, spans, cr = (float(surface['sspn']), float(surface['sspne']),
                       float(surface['chrdr']))
    dd = 2.0 * (span - spans)
    tan_le = a_block[62] if a_block[62] != 0.0 else 0.00001
    taper, arstar, crstar = a_block[27], a_block[7], a_block[10]
    result: Dict[str, object] = {'dd': dd}
    supersonic_kbw = False
    if mach == 1.0:
        beta = 0.0000001
    else:
        beta = math.sqrt(abs(mach**2 - 1.0))
        if taper == 0.0:
            supersonic_kbw = beta * arstar > 1.0
        else:
            trino = beta * arstar * (1.0 + taper) * (1.0 + tan_le / beta)
            result['trino'] = trino
            supersonic_kbw = trino > 4.0
    if supersonic_kbw:
        fig311_args = [beta * dd / crstar, beta / tan_le]
        if (x_surface + cr) / body_length <= 1.0:
            rkbw = interx(2, _T4311A, fig311_args, [15, 19], _D4311A, lind=19,
                          lx1l=2, lx1u=1)
        else:
            rkbw = interx(2, _T4311B, fig311_args, [15, 9], _D4311B, lind=15,
                          lx1l=2, lx1u=1)
        kbw = rkbw / (RAD * beta * (sref / a_block[3]) * cla_surface *
                      (taper + 1.0) * (2.0 * span / dd - 1.0))
        result['rkbw'] = rkbw
    else:
        kbw = interx(1, _TFIG10, [dd / (2.0 * span)], [11], _DKBW10,
                     lind=11)
    alpha0_offset = 0.0 if alpha0_body == UNUSED else alpha0_body
    kwb = interx(1, _TFIG10, [(span - spans) / span], [11], _DKWB10,
                 lind=11)
    clawb, clabw = cla_surface * kwb, cla_surface * kbw
    fig337_args = [0.98 * dd / crstar, 0.98 / tan_le]
    if (x_surface + cr) / body_length <= 1.0:
        xaca = interx(2, _T4337A, fig337_args, [8, 3], _D4337A, lind=8, lx1u=1)
    else:
        xaca = interx(2, _T4337B, fig337_args, [10, 2], _D4337B, lind=10)
    result.update({
        'kbw': float(kbw), 'kwb': float(kwb), 'clawb': float(clawb),
        'clabw': float(clabw), 'cla': float(clabw + clawb + cla_body),
        'xaca': float(xaca), 'xacbw': float(xaca * crstar / cbarr),
        'alpha_body': np.asarray(alpha_deg, dtype=float) + alpha0_offset,
        'method': 'legacy_wbtran'})
    return result


def _fig26_af(tan_over_beta: float, ar_tanle: float, taper: float,
              supersonic: bool) -> float:
    """Figure 4.1.4.2-26A-F at one Mach, as TRANCM's two loops read it."""
    below, above = (_SUPAF1, _SUPAF2) if supersonic else (_SUBAF1, _SUBAF2)
    table = below if tan_over_beta <= 1.0 else above
    abscissa = tan_over_beta if tan_over_beta <= 1.0 else 1.0 / tan_over_beta
    return interx(3, _T422AF, [abscissa, ar_tanle, taper], [6, 6, 6], table,
                  lind=6, lx2l=2, lx2u=2)


def calculate_trancm(mach: float, mfb: float, tovc: float,
                     a: Mapping[int, float], cla: float, cbarr: float,
                     xcg: float, x_surface: float, wgpl: bool,
                     body: Optional[Dict[str, object]] = None
                     ) -> Dict[str, object]:
    """Translate TRANCM (and TRHTCM): transonic aerodynamic centre and CMa.

    TRHTCM is the same code on the tail's blocks; its ``SYNTSS`` puts
    ``XH`` where TRANCM reads ``XW``.  The aerodynamic centre is faired by
    TRANAC through six Mach numbers: Figure 4.1.4.2-26A-F at 0.6 and 1.4
    (with slopes from 0.5/0.7 and 1.3/1.5) and Figure 4.1.4.2-30A-D at the
    four transonic similarity points ``sqrt(1 + v*t^(2/3))``,
    ``v = -2, -1, 0, 1``.  Above 7% thickness it is refaired through eight
    points with the Figure 4.1.4.2-33 increment at the fifth.

    Args:
        mach: ``TRA(4)``.
        mfb: ``TRA(6)``, the force-break Mach number.
        tovc: ``WINGIN(16)``.
        a: The surface's ``A`` block: 7, 10, 27, 62, 73, and 173 (the CG
            arm ``DXCG``, replaced by ``XCG - XW`` when the wing is present
            without a body); with a body also 3 and 44.
        cla: ``WING(101)`` (``HT(101)``).
        cbarr: ``CBARR``.
        xcg, x_surface: ``XCG``, and ``XW`` (``XH`` for TRHTCM).
        wgpl: ``WGPL``, whether the *wing* is present, for both routines.
        body: With a body (``BO``): ``alpha_deg``, ``surface``,
            ``body_length``, ``alpha0_body``, ``cla_body`` and ``sref`` as
            :func:`calculate_wbtran` takes them, ``cma_body``
            (``BODY(121)``), and for WBCM1 ``bd87`` and optionally
            ``stale_wb13``/``stale_wb14``.  For TRHTCM, ``xacbw_14`` is
            the wing's ``SWB(8)``; see Notes.

    Returns:
        ``xac`` (``TRA(105)``, on the root chord), ``xacw`` (``TRA(95)``),
        ``cma`` (``WING(121)``), ``xmv``/``xacv`` (``TRA(83..94)``), and
        ``delxac`` (``TRA(96)``) and ``xacp`` (``TRA(97..104)``) when
        refaired.  With a body: ``xacbw`` (``TRA(106)``), ``xacwb``
        (``TRA(107)``), ``cma_wing_body`` (``BW(121)``), and the
        ``wbtran`` and ``wbcm1`` results.

    Notes:
        The carryover centre is interpolated from Mach 0.6 to 1.4 with
        ``ABS(XACBW4-XACBW6)/0.80``, so it always moves aft with Mach
        number whichever end is further aft.  Kept.

        TRHTCM declares ``/SUPWBB/ SWB(61),SS(61)`` and reads ``SWB(8)`` as
        its Mach 1.4 anchor, but HBTRAN declares ``SSS(61),SWB(61)`` and
        writes the tail's value into the second block.  The tail therefore
        interpolates towards the *wing's* WBTRAN centre (stale when the
        wing pass did not run).  Kept: pass the wing's ``SWB(8)`` as
        ``body['xacbw_14']``; without it the tail's own value is used.
    """
    a_block = {int(k): float(v) for k, v in a.items()}
    aspect_star, tan_le, taper, crstar = (a_block[7], a_block[62],
                                          a_block[27], a_block[10])
    dxcg = a_block[173]
    if wgpl and body is None:
        dxcg = xcg - x_surface
    ar_times_tovc_cbrt = aspect_star * tovc**0.3333
    ar_times_tan_le = aspect_star * tan_le
    xmv = [0.60, 0.0, 0.0, 0.0, 0.0, 1.40]
    xacv = [0.0] * 6
    fig425_x1, fig425_x2 = _T425AD[9:18], _T425AD[0:7]
    fig425_x3, fig425_x4 = _T425AD[18:22], _T425AD[27:30]
    for i, vbar in enumerate((-2.0, -1.0, 0.0, 1.0), start=1):
        xmv[i] = math.sqrt(1.0 + vbar * tovc**.6666)
        xacv[i] = float(tlin4x(fig425_x1, fig425_x2, fig425_x3, fig425_x4,
                               _Y425AD, ar_times_tan_le, ar_times_tovc_cbrt,
                               vbar, taper, 0, 0, 0, 0, 2, 1, 0, 1))
    with np.errstate(divide='ignore'):
        subsonic_tan_beta = [
            float(np.float64(tan_le) / math.sqrt(1.0 - m**2))
            for m in _XM[:3]]
        supersonic_beta_tan = [
            float(math.sqrt(m**2 - 1.0) / np.float64(tan_le))
            for m in _XM[3:]]
    fig26_sub = [_fig26_af(t, ar_times_tan_le, taper, False)
                 for t in subsonic_tan_beta]
    dxac1, xacv[0] = (fig26_sub[1] - fig26_sub[0]) / 0.2, fig26_sub[2]
    fig26_sup = [_fig26_af(t, ar_times_tan_le, taper, True)
                 for t in supersonic_beta_tan]
    dxac2, xacv[5] = (fig26_sup[1] - fig26_sup[0]) / 0.2, fig26_sup[2]
    xac = tranac(xmv, xacv, dxac1, dxac2, mach)['value']
    result: Dict[str, object] = {
        'xmv': [float(x) for x in xmv], 'xacv': [float(x) for x in xacv],
        'dxac1': float(dxac1), 'dxac2': float(dxac2), 'dxcg': float(dxcg)}
    if tovc > 0.07:
        delxac = interx(2, _T428, [tovc * 100.0, aspect_star * a_block[73]**2],
                        [4, 7], _D428, lind=7, lx2l=-1, lx1u=1)
        zmt = [0.60, (mfb + .6) / 2.0, mfb, mfb + .03, mfb + .07, mfb + .14,
               xmv[4], 1.4]
        if zmt[5] + 0.01 >= zmt[6]:
            zmt[6] = (zmt[5] + 1.4) / 2.0
        xacp = [tranac(xmv, xacv, dxac1, dxac2, z)['value'] for z in zmt]
        xac = tranac(zmt, xacp, dxac1, dxac2, mach, delxac)['value']
        result.update({'delxac': float(delxac), 'zmt': zmt,
                       'xacp': [float(x) for x in xacp]})
    xacw = xac * crstar / cbarr
    result.update({'xac': float(xac), 'xacw': float(xacw),
                   'cma': float((dxcg / cbarr - xacw) * cla),
                   'method': 'legacy_trancm'})
    if body is None:
        return result
    body_inputs = dict(body)
    wbtran = calculate_wbtran(mach, body_inputs['alpha_deg'],
                              body_inputs['surface'], a_block,
                              x_surface, body_inputs['body_length'],
                              body_inputs['alpha0_body'],
                              cla, body_inputs['cla_body'],
                              body_inputs['sref'], cbarr)
    wbcm1 = calculate_wbcm1(
        a_block, float(body_inputs['surface']['sspn']),
        body_inputs['bd87'], cbarr,
        body_inputs.get('stale_wb13', 0.0),
        body_inputs.get('stale_wb14', 0.0))
    xacbw_m14 = float(body_inputs.get('xacbw_14', wbtran['xacbw']))
    xacbw_m06 = wbcm1['wb13']
    xacbw = xacbw_m06 + abs(xacbw_m14 - xacbw_m06) / 0.80 * (mach - 0.60)
    clab, cmab = float(body_inputs['cla_body']), float(body_inputs['cma_body'])
    cnob = (-cmab / clab * cbarr + dxcg) * clab / cbarr
    xacwb_num = (cnob + xacw * wbtran['clawb'] + xacbw * wbtran['clabw'])
    xacwb_den = clab + wbtran['clawb'] + wbtran['clabw']
    xacwb = xacwb_num / xacwb_den
    result.update({'xacbw': float(xacbw), 'xacwb': float(xacwb),
                   'cma_wing_body': float((dxcg / cbarr - xacwb) *
                                          wbtran['cla']),
                   'wbtran': wbtran, 'wbcm1': wbcm1})
    return result


def calculate_trawbt(wing: Dict[str, float], a: Mapping[int, float],
                     b48: float, tail: Dict[str, float],
                     aht: Mapping[int, float], position: Dict[str, float],
                     cla_wing: float, cd0_wing: float, cla_wing_body: float,
                     clawb_tail: float, clabw_tail: float, cd0_tail: float,
                     xac_wing: float, sref: float, cbarr: float,
                     downwash_gradient: Optional[float] = None,
                     dynamic_pressure_ratio: Optional[float] = None
                     ) -> Optional[Dict[str, float]]:
    """Translate TRAWBT: the transonic wing-body-tail lift and moment slopes.

    The Section 4.4.1 downwash gradient (Figures 4.4.1-66, -67 and -68B, the
    same method as DWASH's) and the wake dynamic-pressure loss, applied to
    the tail-body carryover lift.  Returns ``None`` when the tail is not
    straight tapered or the wing semispan is under 1.5 tail semispans.

    Args:
        wing: ``sspn`` (``WINGIN(4)``), ``span_break`` (12), and the
            inboard and outboard dihedrals ``dihedral_in`` (13) and
            ``dihedral_out`` (14), degrees.
        a: The wing's ``A`` block: 3, 11, 12, 16, 24, 40, 118, 120, 126,
            127.
        b48: ``B(48)``, the wing's Mach-zero lift-curve slope.
        tail: ``type`` (``HTIN(15)``) and ``sspn`` (``HTIN(4)``).
        aht: The tail's ``A`` block: 161.
        position: ``xcg``, ``xw``, ``xh``, ``aliw`` from ``/SYNTSS/``.
        cla_wing, cd0_wing: ``WING(101)``, ``WING(1)``.
        cla_wing_body: ``BW(101)``.
        clawb_tail, clabw_tail: ``TRAH(71)``, ``TRAH(72)`` from HBTRAN.
        cd0_tail: ``HT(1)``.
        xac_wing: ``TRA(105)`` from TRANCM.
        downwash_gradient, dynamic_pressure_ratio: ``DWASH(41)`` and
            ``DWASH(1)`` when supplied as experimental data (``KDEODA``,
            ``KQOQIN``); otherwise computed.

    Returns:
        ``deda`` (``DWASH(41)``), ``q_ratio`` (``DWASH(1)``), ``cla``
        (``BWHV(101)``), ``cma`` (``BWHV(121)``, also ``BWH(121)``),
        ``cd0_tail`` (``TRAH(108)``), and the intermediate ``TRA`` words.

    Notes:
        The moment slope subtracts ``TRA(105)``, TRANCM's aerodynamic centre
        as a fraction of the *root chord*, from the dimensional ``XCG-XW``.
        Kept.
    """
    if not _stra(tail['type']) or float(wing['sspn']) < 1.5 * float(tail['sspn']):
        return None
    a_block = {int(k): float(v) for k, v in a.items()}
    sspn = float(wing['sspn'])
    dihedral_blend = ((float(position['aliw']) - a_block[126]) /
                      (a_block[127] - a_block[126]))
    tip_to_ob = a_block[24] / sspn
    sweep_c4, taper_ratio = a_block[40], a_block[118]
    fig155a = _tlinex(_X155A1, _X155A2, _Y4155A, sweep_c4, dihedral_blend,
                      0, 1, 2, 0)
    aeefoa = _tlinex(_X155B1, _X155B2, _Y4155B, taper_ratio, fig155a,
                     0, 0, 0, 0)
    beffob = _tlinex(_X155C1, _X155C2, _Y4155C, taper_ratio, aeefoa,
                     0, 1, 0, 0)
    aeff = a_block[120] * aeefoa
    an1 = _tlinex(_X156A1, _X156A2, _Y156A, aeff, tip_to_ob, 2, 1, 2, 2)
    an2 = _tlinex(_X156B1, _X156B2, _Y4156B, sweep_c4, aeff, 0, 2, 2, 2)
    downwash_y = .8 - (-an1 * an2 + an1) / 5.0
    fans = .8 - downwash_y + an2
    beff = 2.0 * sspn * beffob
    span_outboard = sspn - float(wing['span_break'])
    dihedral_height = beff * math.tan(abs(float(wing['dihedral_in']) / RAD))
    if beff / 2.0 > span_outboard:
        dihedral_height = (span_outboard * math.tan(float(wing['dihedral_in']) / RAD) +
                           (beff / 2.0 - span_outboard) *
                           math.tan(float(wing['dihedral_out']) / RAD))
    vortex_axis = a_block[12] - dihedral_height / 2.0
    debode = _tlinex(_X157B1, _X157B2, _Y4157B, abs(2.0 * vortex_axis / beff),
                     2.0 * float(tail['sspn']) / beff, 0, 2, 0, 2)
    deda = (debode * fans * cla_wing / b48 if downwash_gradient is None
            else float(downwash_gradient))
    span_ratio = a_block[24] / a_block[16]
    zwc = 0.68 * math.sqrt(cd0_wing * sref / a_block[3] * (span_ratio + 0.15))
    zc = span_ratio * math.tan(a_block[11])
    dqoq = 2.42 * math.sqrt(cd0_wing * sref / a_block[3]) / (span_ratio + 0.3)
    wake_angle = PI * zc / (2.0 * zwc)
    qoq = (1.0 - dqoq * math.cos(wake_angle)**2
           if dynamic_pressure_ratio is None
           else float(dynamic_pressure_ratio))
    tail_cla_contrib = (clawb_tail + clabw_tail) * (1.0 - deda) * qoq
    cma = (-(float(position['xcg']) - float(position['xw']) - xac_wing) /
           cbarr * cla_wing_body -
           (float(position['xcg']) - float(position['xh']) -
            float(aht[161])) / cbarr * tail_cla_contrib)
    return {'deda': float(deda), 'q_ratio': float(qoq),
            'cla': float(cla_wing_body + tail_cla_contrib), 'cma': float(cma),
            'cd0_tail': float(qoq * cd0_tail), 'fans': float(fans),
            'debode': float(debode), 'aeff': float(aeff),
            'beff': float(beff), 'zwc': float(zwc), 'zc': float(zc),
            'dqoq': float(dqoq), 'dj': float(wake_angle),
            'method': 'legacy_trawbt'}


def m24o30_body_words(cla_body: float, cma_body: float, cbarr: float,
                      blref: float) -> Dict[str, float]:
    """M24O30's body words: the normal-force slope and the yawing slope.

    ``BODY(141) = -BODY(101)`` and ``BODY(161) = -BODY(121)*CBARR/BLREF``;
    the overlay otherwise only marks the body's per-angle words unused and
    sequences TRSONI, TRSONJ and TRANCD.
    """
    return {'body_141': -float(cla_body),
            'body_161': -float(cma_body) * cbarr / blref}


def calculate_clbclc(data: Sequence[float], nalpha: int) -> float:
    """Translate CLBCLC: the first usable ``CLB/CL`` ratio of an IOM block.

    ``DATA(J+180)/DATA(J+20)`` at the first angle where the lift is not
    zero or unused and the ratio's numerator is not unused; ``UNUSED`` when
    no angle qualifies.  ``data`` is the block one-based, ``data[1]`` being
    ``DATA(1)``; a sequence is taken as zero-based.
    """
    get = (data.__getitem__ if isinstance(data, Mapping)
           else (lambda k: data[k - 1]))
    for angle_slot in range(1, nalpha + 1):
        cl, clb = float(get(angle_slot + 20)), float(get(angle_slot + 180))
        if abs(cl) > UNUSED and abs(clb) != UNUSED:
            return clb / cl
    return UNUSED


def calculate_wbclb(alpha_deg: Sequence[float], body_alpha_deg: Sequence[float],
                    surface: Dict[str, float], body_x: Sequence[float],
                    body_s: Sequence[float], x_quarter_chord: float,
                    taper_ratio: float, incidence: float,
                    cl_surface: Sequence[float], cla_surface: float,
                    cl_body: Sequence[float], cla_body: float,
                    cd_surface: Sequence[float], cd_body: Sequence[float],
                    kwb: float, kbw: float, cla_wing_body: float,
                    clb_cl: Dict[str, float], mach: float, mfb: float,
                    cl_wing_body: Sequence[float],
                    cd_wing_body: Sequence[float],
                    clb_wing_body: Sequence[float]) -> Dict[str, object]:
    """Translate WBCLB: second-level transonic surface-body lift and CLB.

    Fills the surface-body lift where it is still unused: the carryover
    buildup ``CLB + (K_W(B)+K_B(W))*(CL_W - CL_i)`` (or the linear
    ``(CLa_B + (K_W(B)+K_B(W))*CLa_W)*alpha`` where a component curve is
    missing), the Figure 4.3.1.2-12A incidence carryover, and at
    ``r/s >= 1/3`` BODOWG's vortex term.  The drag is summed where unused,
    and ``CLB`` is ``CL`` times ``CLB/CL^2`` interpolated from the
    force-break Mach number to 1.4, times ``CL_a,WB^2``.

    Args:
        alpha_deg: ``FLC(23)`` onward.
        body_alpha_deg: ``BD(255)`` onward (WBTRAN's body angles).
        surface: ``sspn`` (``WINGIN(4)``).
        body_x, body_s: The body stations and areas (``/BODYI/``), from
            which BODOWG takes the maximum area.
        x_quarter_chord: ``A(161)+XW``.
        taper_ratio: ``A(27)``.
        incidence: ``ALIW``.
        cl_surface, cla_surface: ``WING(21)`` onward, ``WING(101)``.
        cl_body, cla_body: ``BODY(21)`` onward, ``BODY(101)``.
        cd_surface, cd_body: ``WING(1)``, ``BODY(1)`` onward.
        kwb, kbw: ``SWB(35)``, ``SWB(11)``.
        cla_wing_body: ``BW(101)``.
        clb_cl: ``clb_14``, ``cna_14`` (``SEC(6)``, ``SEC(9)``), and
            ``clb_mfb``, ``cla_mfb`` (``SEC(5)``, ``TRA(12)``).
        mach, mfb: ``FLC(I+2)``, ``TRA(6)``.
        cl_wing_body, cd_wing_body, clb_wing_body: ``BW(21)``, ``BW(1)``,
            ``BW(181)`` onward as they stand; unused entries are filled.

    Returns:
        ``cl``, ``cd``, ``clb`` (the filled curves), ``clbcl``
        (``SEC(21)``), ``kkwb``/``kkbw`` (``SWB(2)``/``(37)``), ``ratio``
        (``FACT(1)``), ``ivbw``/``go2pav`` (``FACT(2..)``/``(22..)``).

    Notes:
        WBCLB passes its local ``YB = SSPN-SSPNE`` as BODOWG's radius
        argument, and BODOWG assigns that argument ``SQRT(BD(3)/PI)`` from
        the body's maximum area.  ``FACT(1) = YB/SSPN`` is formed after the
        call, so the ratio used for the incidence carryover and the vortex
        gate is ``sqrt(Smax/pi)/SSPN``, not the exposed-root offset.
    """
    sspn = float(surface['sspn'])
    _, max_area, _ = getmax(body_x, body_s)
    yb = math.sqrt(max_area / PI)
    bodowg = [calculate_bodowg(a, x_quarter_chord, yb, sspn, taper_ratio)
              for a in body_alpha_deg]
    ivbw = [b['ivbw'] for b in bodowg]
    go2pav = [b['go2pav'] for b in bodowg]
    ratio = yb / sspn
    kkwb, _ = tbfunx(_FIG_431212A_RATIO, _FIG_431212A_KKWB, ratio, 0, 0)
    kkbw, _ = tbfunx(_FIG_431212A_RATIO, _FIG_431212A_KKBW, ratio, 0, 0)
    cl_incidence = cla_surface * incidence
    clwb = [float(v) for v in cl_wing_body]
    cdwb = [float(v) for v in cd_wing_body]
    clbb = [float(v) for v in clb_wing_body]
    has_drag_sum = abs(cdwb[1]) != UNUSED
    for angle_index, alpha in enumerate(alpha_deg):
        if abs(clwb[angle_index]) == UNUSED:
            value = (cla_body + (kwb + kbw) * cla_surface) * alpha
            if abs(cl_body[angle_index]) != UNUSED and \
                    abs(cl_surface[angle_index]) != UNUSED:
                value = (cl_body[angle_index] + (kwb + kbw) *
                         (cl_surface[angle_index] - cl_incidence))
            value += (kkwb + kkbw) * cl_incidence
            if ratio >= 1.0 / 3.0:
                value += (ratio * ivbw[angle_index] * go2pav[angle_index] *
                          cla_surface * (alpha - incidence))
            clwb[angle_index] = value
        if not has_drag_sum and abs(cd_surface[angle_index]) != UNUSED and \
                abs(cd_body[angle_index]) != UNUSED:
            cdwb[angle_index] = (cd_body[angle_index] +
                                 cd_surface[angle_index])
    clb_inputs = {k: float(v) for k, v in clb_cl.items()}
    clbcl_at_mfb = clb_inputs['clb_mfb'] / clb_inputs['cla_mfb']**2
    clbcl = (((clb_inputs['clb_14'] / clb_inputs['cna_14']**2 - clbcl_at_mfb) *
              (mach - mfb) / (1.4 - mfb) + clbcl_at_mfb) * cla_wing_body**2)
    for angle_index in range(len(alpha_deg)):
        if clwb[angle_index] != UNUSED and clbb[angle_index] == UNUSED:
            clbb[angle_index] = clbcl * clwb[angle_index]
    return {'cl': clwb, 'cd': cdwb, 'clb': clbb, 'clbcl': float(clbcl),
            'kkwb': float(kkwb), 'kkbw': float(kkbw), 'ratio': float(ratio),
            'ivbw': ivbw, 'go2pav': go2pav, 'method': 'legacy_wbclb'}



def _setup2_mach(state: Dict[str, object], mach: float, beta: float,
                 i: int) -> None:
    state['mach'] = mach
    state['b'] = [mach, beta]
    state['bht'] = [mach, beta]
    state['wingin'][i + 20] = state['wingin'][69] / beta
    state['htin'][i + 20] = state['htin'][69] / beta


def setup2_step(nf: int, state: Dict[str, object]) -> int:
    """Translate SETUP2: the second-level method's Mach schedule.

    Each call runs steps ``N = -NF`` onward, stopping after the first whose
    configuration gate is open (or after the last), and returns the new
    ``NF``.  The steps set the flight Mach number for the next pass and
    harvest the ``CLB/CL`` ratios (CLBCLC) of the pass before:

    1. Mach 0.6, subsonic; keeps the flight Mach and both force-break
       Mach numbers.  Gate: wing or tail.
    2. Mach 0.7; ``SEC(1)``, ``SEC(3)`` from the wing and tail, ``SEC(11)``
       from ``WBT(67)``.  Gate: body, wing and tail.
    3. The wing's force-break Mach (at most 0.95); ``SEC(12)``.  Gate:
       body and wing.
    4. The tail's force-break Mach (at most 0.95); ``SEC(5)`` from BW.
       Gate: body and tail.
    5. Mach 1.4, supersonic; ``SEC(7)`` from BH.  Gate: wing or tail.
    6. Mach 1.1; ``SEC(2,4,6,8)``, the two carryover slopes and
       ``SEC(14)``.  Gate: body, wing and tail.
    7. Restores the flight Mach and the angle count, marks ``DONE``, and
       sets the entries of absent components ``UNUSED``.

    ``state`` is updated in place.  Keys: ``i`` (``I``, default 1),
    ``mach`` (``FLC(I+2)``), ``flc2`` (``FLC(2)``), ``nalpha``; the flags
    ``bo``, ``wgpl``, ``htpl``, ``subson``, ``transn``, ``supers`` and
    ``done``; ``sec``, ``SEC(1..23)`` as a one-based dict, whose 17-19 hold
    the saved flight Mach and the two force-break Mach numbers; ``tra6``,
    ``trah6``; ``b`` and ``bht`` (``B(1..2)``, ``BHT(1..2)``); ``wingin``
    and ``htin`` (one-based dicts: 68 and 69 read, ``I+20`` and ``I+40``
    written); ``wbt67`` and ``stp155`` (``WBT(67)``, ``WBT(155)``); the IOM
    blocks ``wing``, ``ht``, ``bw``, ``bh`` for CLBCLC; ``bw101``,
    ``bh101``.
    """
    sec = state['sec']
    i = int(state.get('i', 1))
    while True:
        step_num = -nf
        if step_num == 1:
            sec[17] = state['mach']
            state['subson'], state['transn'] = True, False
            sec[18], sec[19] = state['tra6'], state['trah6']
            _setup2_mach(state, 0.6, 0.8, i)
            state['wingin'][i + 40] = state['wingin'][68]
            state['htin'][i + 40] = state['htin'][68]
            proceed = state['wgpl'] or state['htpl']
        elif step_num == 2:
            sec[11] = state['wbt67']
            sec[1] = calculate_clbclc(state['wing'], state['nalpha'])
            sec[3] = calculate_clbclc(state['ht'], state['nalpha'])
            _setup2_mach(state, 0.7, 0.71414284, i)
            proceed = state['bo'] and state['wgpl'] and state['htpl']
        elif step_num == 3:
            sec[12] = state['wbt67']
            sec[18] = min(sec[18], 0.95)
            _setup2_mach(state, sec[18], math.sqrt(1.0 - sec[18]**2), i)
            proceed = state['bo'] and state['wgpl']
        elif step_num == 4:
            sec[5] = calculate_clbclc(state['bw'], state['nalpha'])
            sec[19] = min(sec[19], 0.95)
            _setup2_mach(state, sec[19], math.sqrt(1.0 - sec[19]**2), i)
            proceed = state['bo'] and state['htpl']
        elif step_num == 5:
            sec[7] = calculate_clbclc(state['bh'], state['nalpha'])
            state['subson'], state['supers'] = False, True
            state['mach'] = 1.4
            proceed = state['wgpl'] or state['htpl']
        elif step_num == 6:
            for k, block in ((2, 'wing'), (4, 'ht'), (6, 'bw'), (8, 'bh')):
                sec[k] = calculate_clbclc(state[block], state['nalpha'])
            sec[9], sec[10], sec[14] = state['bw101'], state['bh101'], state['stp155']
            state['mach'] = 1.1
            proceed = state['bo'] and state['wgpl'] and state['htpl']
        elif step_num == 7:
            sec[13] = state['stp155']
            state['mach'] = sec[17]
            state['nalpha'] = int(state['flc2'] + 0.5)
            state['done'] = True
            state['supers'], state['transn'] = False, True
            absent = []
            if not state['wgpl']:
                absent += [1, 2]
            if not state['htpl']:
                absent += [3, 4]
            if not (state['bo'] and state['wgpl']):
                absent += [5, 6, 9]
            if not (state['bo'] and state['htpl']):
                absent += [7, 8, 10]
            if not (state['bo'] and state['wgpl'] and state['htpl']):
                absent += [11, 12, 13, 14]
            for k in absent + list(range(17, 24)):
                sec[k] = UNUSED
            proceed = True
        else:
            raise ValueError(f"SETUP2 has no step {step_num}")
        nf -= 1
        if proceed:
            return nf


def calculate_supcm0(surface: Dict[str, float], x_surface: float,
                     wing_height: float, body_length: float,
                     max_diameter: float, mach: float, tr: float,
                     stale_cm0: float = 0.0) -> Dict[str, object]:
    """Translate SUPCM0: the supersonic surface-body CM0.

    Line for line TRACM0 (with its Mach-number Reynolds number and the
    wing's ``ZW`` for the tail), storing ``TRA(73)``/``TRAH(73)`` where
    TRACM0 stores word 74.  Arguments and returns are those of
    :func:`calculate_tracm0`.
    """
    return calculate_tracm0(surface, x_surface, wing_height, body_length,
                            max_diameter, mach, tr, stale_cm0)
