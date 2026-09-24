"""
Supersonic wing and horizontal-tail lift, pitching moment and drag: SUPLNG
and SUPLTG, with their overlays M27O33 and M22O26.

SUPLNG forms the wing's supersonic normal-force slope and its nonlinear
normal force, the aerodynamic centre and pitching-moment slope, and the
drag, normal and axial force at each angle of attack:

* a straight-tapered wing: the theoretical slope of Figures 4.1.3.2-56A
  to G with the leading-edge correction of -60A/B; the nonlinear normal
  force of 4.1.3.3-59 (subsonic leading edge) or -60 (supersonic leading
  edge, with the shock-detachment transition of -60B and -61A);
* a cranked or double-delta wing: the basic wing, glove and trailing-edge
  extension components of Figures 4.1.3.2-56, -61, -62 and -63, and the
  centre of pressure blended from the inboard and outboard panels, whose
  slopes come from re-entering the straight-wing computation.

SUPLTG is the same computation on the horizontal tail's blocks, with its
own skin-friction, wave and drag-due-to-lift terms.

The aerodynamic centre comes from Figures 4.1.4.2-26A to F (FWDXAC for a
forward-swept edge).  The ``/SUPWH/`` words (``SLG``), the wing's
``WINGD`` words and the ``WING`` result words 1-121 are mirrored as flat
1-based arrays.

The tables were extracted from the source DATA statements by parsing, with
``tools/fortran_data.py``.

Reference: datcom-legacy/datcom_2000/suplng.f, supltg.f, m27o33.f,
m22o26.f
"""

import math
from typing import Dict, Mapping

import numpy as np

from pydatcom.aerodynamics.fwdxac import calculate_fwdxac
from pydatcom.utils.constants import PI, RAD, UNUSED
from pydatcom.utils.legacy_interp import interx
from pydatcom.utils.legacy_numeric import tbfunx
from pydatcom.utils.table_lookup import angdet, fig26, fig60b

_A1350 = np.array([
    0.0, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0001, 0.0,
    0.41, 0.82, 1.24, 2.12, 3.18, 6.95, 16.1,
])
_DA50 = np.array([
    1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.05,
    1.05, 1.05, 1.05, 0.985, 0.945, 0.915, 0.9, 0.9, 1.04, 1.04,
    1.04, 1.04, 0.965, 0.908, 0.87, 0.85, 0.84, 1.12, 1.12, 1.12,
    1.015, 0.94, 0.88, 0.838, 0.81, 0.796, 1.11, 1.11, 1.11, 1.0,
    0.903, 0.84, 0.795, 0.765, 0.75, 1.08, 1.08, 1.08, 0.954, 0.865,
    0.8, 0.75, 0.72, 0.7, 1.2, 1.2, 1.043, 0.907, 0.817, 0.75,
    0.707, 0.675, 0.66, 1.14, 1.14, 0.975, 0.857, 0.772, 0.717, 0.675,
    0.65, 0.632,
])
_B1350 = np.array([
    0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9,
    1.0001, 0.0, 4.0, 8.0, 12.0, 20.0, 30.0, 50.0, 70.0,
])
_DB50 = np.array([
    1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
    1.0, 1.15, 1.15, 1.15, 1.15, 1.095, 1.04, 0.99, 0.96, 0.935,
    0.915, 0.9, 1.12, 1.12, 1.12, 1.12, 1.05, 0.985, 0.94, 0.905,
    0.88, 0.857, 0.84, 1.15, 1.15, 1.15, 1.08, 1.005, 0.945, 0.902,
    0.87, 0.842, 0.82, 0.796, 1.22, 1.14, 1.05, 0.98, 0.93, 0.89,
    0.853, 0.823, 0.795, 0.77, 0.75, 1.13, 1.05, 0.98, 0.925, 0.88,
    0.845, 0.81, 0.782, 0.752, 0.73, 0.7, 1.02, 0.942, 0.895, 0.855,
    0.82, 0.79, 0.76, 0.735, 0.71, 0.685, 0.66, 1.0, 0.92, 0.87,
    0.825, 0.79, 0.755, 0.728, 0.7, 0.678, 0.655, 0.632,
])
_G13246 = np.array([
    0.0, 0.2, 0.4, 0.45, 0.5, 0.6, 0.8, 0.9, 1.0,
])
_DG3246 = np.array([
    1.61, 1.58, 1.55, 1.57, 1.62, 1.75, 1.94, 2.0, 2.0,
])
_T13251 = np.array([
    0.015, 0.03, 0.05, 0.075, 0.11, 0.16, 0.23, 0.45, 0.9, 1.2,
    1.7, 2.0,
])
_D13251 = np.array([
    0.963, 0.94, 0.92, 0.9, 0.88, 0.86, 0.84, 0.8, 0.76, 0.74,
    0.72, 0.709,
])
_T13246 = np.array([
    0.0, 0.1, 0.2, 0.3, 0.33, 0.4, 0.5, 0.6, 0.7, 0.8,
    0.9, 1.0, 1.111, 1.25, 1.429, 1.667, 2.0, 2.5, 2.941, 4.167,
    7.143, 14.286, 30.0, 0.25, 0.5, 1.0, 2.0, 3.0, 4.0, 5.0,
    6.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.2, 0.25, 0.3333,
    0.5, 1.0,
])
_DUMY1 = np.array([
    0.39, 0.39, 0.39, 0.39, 0.39, 0.39, 0.39, 0.39, 0.39, 0.39,
    0.39, 0.41, 0.44, 0.5, 0.58, 0.65, 0.8, 1.0, 1.15, 1.55,
    2.4, 3.92, 4.0, 0.77, 0.77, 0.78, 0.78, 0.79, 0.79, 0.8,
    0.8, 0.8, 0.81, 0.82, 0.84, 0.92, 1.0, 1.18, 1.32, 1.6,
    1.9, 2.18, 2.8, 3.84, 3.92, 4.0, 1.55, 1.56, 1.57, 1.57,
    1.59, 1.59, 1.6, 1.63, 1.66, 1.68, 1.7, 1.75, 1.88, 2.09,
    2.3, 2.6, 2.94, 3.35, 3.7, 3.8, 3.98, 3.98, 4.0, 3.15,
    3.15, 3.15, 3.15, 3.15, 3.15, 3.17, 3.19, 3.23, 3.27, 3.33,
    3.4, 3.46, 3.54, 3.6, 3.7, 3.75, 3.8, 3.86, 3.91, 4.0,
    4.0, 4.0, 4.71, 4.74, 4.83, 5.09, 5.25, 5.05, 4.8, 4.55,
    4.3, 4.09, 3.9, 3.72, 3.78, 3.8, 3.83, 3.88, 3.9, 3.94,
    3.92, 3.96, 4.0, 4.0, 4.0, 6.29, 6.2, 5.99, 5.72, 5.61,
    5.42, 5.18, 4.9, 4.64, 4.42, 4.2, 4.0, 4.0, 4.0, 4.0,
    4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 6.3, 6.34,
    6.39, 6.08, 5.97, 5.78, 5.5, 5.2, 4.95, 4.7, 4.48, 4.25,
    4.2, 4.18, 4.14, 4.11, 4.07, 4.05, 4.04, 4.02, 4.0, 4.0,
    4.0, 6.32, 6.4, 6.4, 6.34, 6.32, 6.13, 5.82, 5.51, 5.25,
    4.99, 4.73, 4.5, 4.4, 4.32, 4.27, 4.2, 4.13, 4.11, 4.07,
    4.04, 4.0, 4.0, 4.0,
])
_DUMY2 = np.array([
    0.41, 0.41, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4,
    0.4, 0.41, 0.41, 0.48, 0.56, 0.67, 0.84, 1.08, 1.29, 1.83,
    2.8, 3.8, 4.0, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8,
    0.8, 0.8, 0.8, 0.83, 0.85, 0.98, 1.09, 1.27, 1.48, 1.8,
    2.25, 2.53, 3.22, 3.79, 3.92, 4.0, 1.57, 1.59, 1.6, 1.6,
    1.6, 1.61, 1.63, 1.65, 1.69, 1.73, 1.77, 1.8, 2.0, 2.26,
    2.57, 2.9, 3.27, 3.53, 3.65, 3.8, 3.91, 3.97, 4.0, 3.17,
    3.17, 3.22, 3.4, 3.5, 3.6, 3.7, 3.67, 3.58, 3.49, 3.37,
    3.23, 3.38, 3.52, 3.62, 3.71, 3.8, 3.88, 3.9, 3.95, 3.97,
    3.99, 4.0, 4.72, 5.0, 5.1, 5.02, 4.97, 4.84, 4.64, 4.45,
    4.24, 4.05, 3.89, 3.7, 3.84, 3.91, 3.97, 3.99, 4.0, 4.0,
    4.0, 4.0, 3.99, 4.0, 4.0, 5.57, 5.6, 5.62, 5.61, 5.6,
    5.42, 5.15, 4.91, 4.7, 4.48, 4.25, 4.08, 4.13, 4.18, 4.16,
    4.11, 4.09, 4.05, 4.02, 4.01, 4.0, 4.0, 4.0, 5.73, 5.77,
    5.79, 5.79, 5.77, 5.75, 5.6, 5.33, 5.08, 4.83, 4.59, 4.4,
    4.42, 4.4, 4.32, 4.23, 4.17, 4.08, 4.05, 4.03, 4.01, 4.0,
    4.0, 5.83, 5.88, 5.91, 5.93, 5.93, 5.91, 5.88, 5.68, 5.4,
    5.14, 4.9, 4.7, 4.68, 4.59, 4.46, 4.33, 4.23, 4.13, 4.1,
    4.06, 4.02, 4.01, 4.0,
])
_DUMY3 = np.array([
    0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4,
    0.4, 0.42, 0.43, 0.5, 0.6, 0.71, 0.9, 1.12, 1.35, 1.97,
    2.95, 3.8, 4.0, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8,
    0.81, 0.82, 0.83, 0.83, 0.88, 0.95, 1.06, 1.23, 1.46, 1.83,
    2.3, 2.65, 3.32, 3.76, 3.92, 4.0, 1.59, 1.59, 1.59, 1.59,
    1.59, 1.6, 1.62, 1.65, 1.68, 1.72, 1.78, 1.82, 2.0, 2.29,
    2.65, 3.0, 3.27, 3.5, 3.6, 3.78, 3.9, 3.98, 4.0, 3.14,
    3.2, 3.41, 3.62, 3.68, 3.75, 3.7, 3.62, 3.51, 3.42, 3.3,
    3.18, 3.34, 3.5, 3.6, 3.7, 3.8, 3.85, 3.88, 3.93, 3.98,
    4.0, 4.0, 4.7, 4.85, 5.02, 4.98, 4.91, 4.8, 4.6, 4.4,
    4.2, 4.0, 3.82, 3.68, 3.8, 3.9, 3.98, 4.0, 4.0, 4.0,
    4.0, 4.0, 4.0, 4.0, 4.0, 5.4, 5.44, 5.45, 5.45, 5.42,
    5.35, 5.12, 4.9, 4.69, 4.48, 4.22, 4.02, 4.12, 4.18, 4.17,
    4.15, 4.1, 4.05, 4.04, 4.01, 4.0, 4.0, 4.0, 5.61, 5.63,
    5.64, 5.62, 5.61, 5.6, 5.54, 5.3, 5.02, 4.8, 4.6, 4.39,
    4.42, 4.44, 4.38, 4.26, 4.16, 4.1, 4.07, 4.03, 4.01, 4.0,
    4.0, 5.72, 5.75, 5.76, 5.78, 5.77, 5.76, 5.72, 5.69, 5.4,
    5.15, 4.9, 4.68, 4.68, 4.6, 4.49, 4.33, 4.22, 4.15, 4.1,
    4.06, 4.02, 4.01, 4.0,
])
_DUMY4 = np.array([
    0.41, 0.41, 0.41, 0.41, 0.41, 0.41, 0.41, 0.41, 0.41, 0.41,
    0.41, 0.41, 0.43, 0.5, 0.6, 0.71, 0.89, 1.11, 1.35, 2.0,
    3.0, 3.72, 4.0, 0.82, 0.82, 0.81, 0.8, 0.8, 0.81, 0.82,
    0.83, 0.84, 0.85, 0.87, 0.89, 0.98, 1.1, 1.3, 1.5, 1.82,
    2.23, 2.68, 3.3, 3.7, 3.91, 4.0, 1.6, 1.59, 1.58, 1.58,
    1.58, 1.59, 1.6, 1.62, 1.66, 1.73, 1.81, 1.92, 2.2, 2.45,
    2.7, 2.98, 3.22, 3.45, 3.59, 3.75, 3.7, 3.91, 4.0, 3.13,
    3.18, 3.32, 3.6, 3.64, 3.72, 3.7, 3.64, 3.54, 3.45, 3.3,
    3.14, 3.32, 3.48, 3.6, 3.7, 3.79, 3.87, 3.9, 3.95, 3.99,
    4.0, 4.0, 4.71, 4.76, 4.8, 4.83, 4.84, 4.7, 4.53, 4.34,
    4.18, 4.0, 3.82, 3.63, 3.78, 3.9, 3.98, 4.0, 4.0, 4.0,
    4.0, 4.0, 4.0, 4.0, 4.0, 5.2, 5.22, 5.22, 5.2, 5.2,
    5.15, 5.1, 4.88, 4.63, 4.42, 4.21, 4.0, 4.15, 4.21, 4.21,
    4.17, 4.14, 4.07, 4.05, 4.03, 4.01, 4.0, 4.0, 5.45, 5.47,
    5.47, 5.45, 5.45, 5.41, 5.36, 5.3, 5.02, 4.81, 4.58, 4.36,
    4.45, 4.47, 4.41, 4.3, 4.19, 4.12, 4.08, 4.02, 4.0, 4.0,
    4.0, 5.58, 5.59, 5.59, 5.59, 5.59, 5.58, 5.57, 5.53, 5.41,
    5.14, 4.91, 4.65, 4.7, 4.65, 4.53, 4.4, 4.27, 4.16, 4.12,
    4.05, 4.02, 4.01, 4.0,
])
_DUMY5 = np.array([
    0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4,
    0.4, 0.4, 0.44, 0.52, 0.63, 0.78, 0.98, 1.21, 1.4, 2.07,
    3.07, 3.6, 4.0, 0.8, 0.8, 0.79, 0.79, 0.79, 0.79, 0.8,
    0.81, 0.83, 0.85, 0.88, 0.9, 1.01, 1.17, 1.38, 1.62, 2.0,
    2.5, 2.74, 3.22, 3.66, 3.87, 4.0, 1.58, 1.59, 1.62, 1.71,
    1.76, 1.77, 1.8, 1.85, 1.88, 1.91, 1.92, 1.93, 2.13, 2.4,
    2.6, 2.88, 3.12, 3.39, 3.51, 3.71, 3.87, 3.99, 4.0, 3.15,
    3.21, 3.38, 3.61, 3.7, 3.67, 3.59, 3.5, 3.4, 3.29, 3.17,
    3.03, 3.18, 3.37, 3.51, 3.67, 3.78, 3.85, 3.88, 3.93, 3.98,
    4.0, 4.0, 4.42, 4.42, 4.41, 4.4, 4.39, 4.38, 4.35, 4.25,
    4.08, 3.9, 3.72, 3.57, 3.7, 3.86, 3.93, 3.99, 4.0, 4.0,
    4.0, 4.0, 4.0, 4.0, 4.0, 4.88, 4.85, 4.81, 4.8, 4.79,
    4.77, 4.73, 4.68, 4.53, 4.35, 4.15, 3.99, 4.08, 4.2, 4.23,
    4.19, 4.1, 4.07, 4.05, 4.01, 4.0, 4.0, 4.0, 5.08, 5.09,
    5.08, 5.08, 5.07, 5.06, 5.03, 4.99, 4.92, 4.74, 4.53, 4.3,
    4.43, 4.48, 4.43, 4.32, 4.2, 4.12, 4.08, 4.03, 4.01, 4.0,
    4.0, 5.19, 5.21, 5.22, 5.23, 5.23, 5.23, 5.23, 5.21, 5.19,
    5.09, 4.88, 4.6, 4.72, 4.69, 4.58, 4.41, 4.28, 4.18, 4.12,
    4.07, 4.02, 4.01, 4.0,
])
_DUMY6 = np.array([
    0.4, 0.4, 0.4, 0.4, 0.4, 0.41, 0.44, 0.49, 0.51, 0.54,
    0.59, 0.61, 0.69, 0.75, 0.81, 0.9, 1.01, 1.26, 1.53, 2.1,
    2.9, 3.42, 4.0, 0.81, 0.81, 0.81, 0.83, 0.84, 0.86, 0.88,
    0.91, 0.99, 1.02, 1.1, 1.19, 1.3, 1.42, 1.58, 1.75, 2.0,
    2.34, 2.6, 3.0, 3.46, 3.71, 4.0, 1.6, 1.58, 1.59, 1.6,
    1.6, 1.61, 1.67, 1.7, 1.79, 1.86, 1.92, 2.0, 2.13, 2.32,
    2.51, 2.75, 2.97, 3.2, 3.35, 3.58, 3.8, 3.9, 4.0, 3.13,
    3.1, 3.08, 3.05, 3.04, 3.03, 3.01, 2.99, 2.95, 2.93, 2.91,
    2.89, 3.02, 3.2, 3.37, 3.52, 3.67, 3.77, 3.82, 3.88, 3.92,
    3.97, 4.0, 3.79, 3.84, 3.86, 3.8, 3.79, 3.72, 3.63, 3.58,
    3.5, 3.47, 3.43, 3.41, 3.57, 3.75, 3.89, 3.98, 4.0, 3.98,
    3.94, 3.93, 3.92, 3.99, 4.0, 4.12, 4.2, 4.19, 4.09, 4.06,
    4.03, 3.99, 3.95, 3.92, 3.9, 3.88, 3.87, 4.01, 4.1, 4.21,
    4.21, 4.12, 4.05, 4.0, 3.96, 3.99, 4.0, 4.0, 4.39, 4.45,
    4.44, 4.39, 4.35, 4.31, 4.27, 4.25, 4.25, 4.25, 4.25, 4.25,
    4.37, 4.5, 4.49, 4.39, 4.23, 4.1, 4.03, 3.99, 4.0, 4.0,
    4.0, 4.6, 4.68, 4.65, 4.59, 4.57, 4.54, 4.5, 4.45, 4.48,
    4.53, 4.59, 4.61, 4.71, 4.75, 4.69, 4.5, 4.31, 4.18, 4.08,
    4.01, 4.0, 4.0, 4.0,
])
_A13359 = np.array([
    0.0, 10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0, 0.0,
    0.0, 0.7, 0.8, 0.9, 1.0, 1.1, 1.25, 1.5, 2.0, 2.5,
    3.0, 3.25,
])
_D13359 = np.array([
    0.0, 3.5, 4.95, 5.85, 6.4, 6.83, 7.25, 7.57, 7.85, 0.0,
    3.0, 4.19, 4.95, 5.58, 6.04, 6.49, 6.8, 7.0, 0.0, 2.34,
    3.4, 4.0, 4.49, 4.91, 5.28, 5.63, 5.92, 0.0, 1.86, 2.78,
    3.2, 3.53, 3.9, 4.25, 4.6, 4.85, 0.0, 1.42, 2.2, 2.68,
    3.05, 3.34, 3.62, 3.87, 4.06, 0.0, 1.1, 1.62, 1.94, 2.19,
    2.4, 2.64, 2.8, 2.95, 0.0, 0.92, 1.35, 1.65, 1.84, 2.0,
    2.17, 2.28, 2.35, 0.0, 0.73, 1.1, 1.35, 1.54, 1.67, 1.76,
    1.83, 1.85, 0.0, 0.62, 0.94, 1.0, 1.1, 1.12, 1.12, 1.1,
    1.09, 0.0, 0.54, 0.72, 0.68, 0.42, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.45, 0.63, 0.51, 0.23, 0.0, 0.0, 0.0, 0.0,
])
_BL3359 = np.array([
    0.5, 0.6, 0.7, 0.8, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.6, 0.8, 1.0, 1.2, 1.6, 2.0, 2.4, 2.6, 2.8, 3.0,
])
_DL3359 = np.array([
    2.03, 1.95, 1.9, 1.86, 1.78, 1.94, 1.87, 1.82, 1.77, 1.7,
    1.82, 1.75, 1.69, 1.64, 1.57, 1.71, 1.63, 1.58, 1.54, 1.46,
    1.51, 1.44, 1.37, 1.32, 1.24, 1.32, 1.24, 1.16, 1.1, 0.99,
    0.98, 0.89, 0.82, 0.75, 0.65, 0.58, 0.52, 0.47, 0.42, 0.37,
    0.18, 0.15, 0.12, 0.13, 0.17, 0.0, 0.0, 0.0, 0.0, 0.0,
])
_BR3359 = np.array([
    0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9,
    1.0, 0.6, 0.8, 1.0, 1.2, 1.6, 2.0, 2.4, 2.6, 2.8,
    3.0,
])
_DR3359 = np.array([
    1.69, 1.69, 1.69, 1.69, 1.69, 1.69, 1.69, 1.7, 1.72, 1.75,
    1.77, 1.69, 1.66, 1.64, 1.62, 1.6, 1.6, 1.61, 1.63, 1.65,
    1.67, 1.69, 1.69, 1.64, 1.59, 1.56, 1.52, 1.51, 1.51, 1.51,
    1.52, 1.54, 1.57, 1.69, 1.61, 1.55, 1.49, 1.44, 1.41, 1.39,
    1.4, 1.41, 1.43, 1.46, 1.69, 1.56, 1.45, 1.36, 1.28, 1.22,
    1.18, 1.17, 1.18, 1.2, 1.23, 1.69, 1.52, 1.34, 1.2, 1.08,
    0.99, 0.94, 0.91, 0.91, 0.94, 0.98, 1.69, 1.4, 1.14, 0.94,
    0.8, 0.7, 0.65, 0.62, 0.61, 0.63, 0.65, 1.69, 1.32, 1.08,
    0.89, 0.74, 0.64, 0.54, 0.47, 0.42, 0.38, 0.37, 1.69, 1.32,
    1.04, 0.84, 0.67, 0.54, 0.45, 0.34, 0.27, 0.2, 0.17, 1.69,
    1.31, 1.02, 0.79, 0.61, 0.46, 0.33, 0.21, 0.12, 0.03, 0.0,
])
_TL360A = np.array([
    0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.26,
    2.4, 2.6, 2.8, 3.0,
])
_DL360A = np.array([
    1.6, 1.35, 1.12, 0.94, 0.79, 0.68, 1.6, 1.35, 1.12, 0.94,
    0.79, 0.68, 1.6, 1.35, 1.12, 0.94, 0.79, 0.68, 1.6, 1.35,
    1.12, 0.94, 0.79, 0.68, 1.6, 1.35, 1.12, 0.94, 0.79, 0.68,
    1.6, 1.35, 1.12, 0.94, 0.79, 0.68, 1.6, 1.35, 1.12, 0.94,
    0.79, 0.68, 1.6, 1.35, 1.12, 0.94, 0.79, 0.68, 1.48, 1.25,
    1.03, 0.85, 0.7, 0.6, 0.84, 0.63, 0.5, 0.42, 0.37, 0.37,
    0.34, 0.25, 0.17, 0.13, 0.11, 0.16, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0,
])
_TR360A = np.array([
    0.0, 0.2, 0.22, 0.24, 0.27, 0.3, 0.32, 0.375, 0.4, 0.471,
    0.6, 0.8, 1.0, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2,
    2.26, 2.4, 2.6, 2.8, 3.0,
])
_DUMY7 = np.array([
    1.68, 1.46, 1.45, 1.4, 1.31, 1.23, 1.16, 1.05, 1.0, 0.88,
    0.72, 0.64, 0.68, 1.68, 1.41, 1.4, 1.4, 1.31, 1.23, 1.16,
    1.05, 1.0, 0.88, 0.72, 0.64, 0.68, 1.68, 1.37, 1.35, 1.33,
    1.31, 1.23, 1.16, 1.05, 1.0, 0.88, 0.72, 0.64, 0.68, 1.68,
    1.3, 1.26, 1.24, 1.21, 1.18, 1.16, 1.05, 1.0, 0.88, 0.72,
    0.64, 0.68, 1.68, 1.24, 1.2, 1.16, 1.14, 1.11, 1.1, 1.05,
    1.0, 0.88, 0.72, 0.64, 0.68,
])
_DUMY8 = np.array([
    1.68, 1.16, 1.11, 1.1, 1.07, 1.04, 1.02, 0.96, 0.94, 0.88,
    0.72, 0.64, 0.68, 1.68, 1.11, 1.07, 1.03, 1.0, 0.97, 0.95,
    0.89, 0.87, 0.81, 0.72, 0.64, 0.68, 1.68, 1.1, 1.03, 1.02,
    0.98, 0.95, 0.93, 0.87, 0.85, 0.79, 0.72, 0.64, 0.68, 1.68,
    1.07, 1.0, 0.98, 0.95, 0.9, 0.87, 0.81, 0.79, 0.72, 0.63,
    0.56, 0.6, 1.68, 1.04, 0.98, 0.94, 0.9, 0.86, 0.83, 0.76,
    0.74, 0.66, 0.54, 0.42, 0.37, 1.68, 1.03, 0.97, 0.92, 0.86,
    0.77, 0.76, 0.69, 0.65, 0.57, 0.42, 0.25, 0.16, 1.68, 1.02,
    0.96, 0.9, 0.84, 0.775, 0.75, 0.67, 0.63, 0.52, 0.34, 0.1,
    0.0,
])
_R13252 = np.array([
    0.0, 0.5, 0.94, 1.0, 1.05, 1.5, 2.0, 2.5, 3.0, 3.5,
    4.0, 20.0,
])
_DRND52 = np.array([
    1.043, 0.97, 0.893, 0.89, 0.894, 0.987, 1.005, 1.02, 1.025, 1.03,
    1.034, 1.034,
])
_S13252 = np.array([
    0.0, 0.5, 0.94, 1.0, 1.05, 1.5, 2.0, 2.5, 3.0, 3.5,
    4.0, 20.0,
])
_DSHP52 = np.array([
    1.08, 1.013, 0.94, 0.942, 0.948, 1.005, 1.03, 1.04, 1.048, 1.05,
    1.055, 1.055,
])
_T13253 = np.array([
    0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8,
    2.0, 2.4, 2.8, 3.2, -1.0, -0.8, -0.6, -0.4, -0.2, 0.0,
    0.2, 0.4, 0.6, 0.8,
])
_D13253 = np.array([
    1.57, 1.69, 1.752, 1.78, 1.774, 1.7, 1.472, 1.288, 1.15, 1.034,
    0.938, 0.795, 0.69, 0.605, 1.57, 1.708, 1.76, 1.778, 1.727, 1.56,
    1.35, 1.18, 1.047, 0.947, 0.86, 0.72, 0.617, 0.543, 1.57, 1.708,
    1.77, 1.737, 1.575, 1.418, 1.233, 1.07, 0.94, 0.842, 0.768, 0.646,
    0.558, 0.49, 1.57, 1.708, 1.763, 1.618, 1.442, 1.278, 1.104, 0.955,
    0.845, 0.758, 0.68, 0.57, 0.49, 0.428, 1.57, 1.708, 1.548, 1.41,
    1.275, 1.133, 0.98, 0.84, 0.738, 0.66, 0.594, 0.495, 0.425, 0.37,
    1.57, 1.49, 1.367, 1.245, 1.12, 0.995, 0.845, 0.725, 0.628, 0.555,
    0.5, 0.423, 0.363, 0.312, 1.254, 1.273, 1.169, 1.062, 0.958, 0.848,
    0.7, 0.59, 0.518, 0.465, 0.415, 0.345, 0.29, 0.253, 0.94, 0.943,
    0.948, 0.875, 0.788, 0.695, 0.55, 0.465, 0.404, 0.36, 0.325, 0.268,
    0.227, 0.196, 0.62, 0.623, 0.63, 0.648, 0.597, 0.537, 0.403, 0.335,
    0.287, 0.25, 0.22, 0.18, 0.15, 0.13, 0.31, 0.315, 0.327, 0.347,
    0.391, 0.35, 0.23, 0.177, 0.147, 0.129, 0.112, 0.09, 0.073, 0.062,
])
_T61A = np.array([
    1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 6.0, 7.0,
])
_SLOPE = np.array([
    59.375, 41.667, 33.0, 27.778, 23.611, 16.667, 12.5, 8.75, 6.0,
])
_T422AF = np.array([
    0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.0, 2.0, 3.0, 4.0,
    5.0, 6.0, 0.0, 0.2, 0.25, 0.33, 0.5, 1.0,
])
_DUMYA = np.array([
    0.165, 0.21, 0.25, 0.29, 0.31, 0.345, 0.335, 0.365, 0.39, 0.415,
    0.445, 0.47, 0.5, 0.54, 0.56, 0.56, 0.56, 0.56, 0.67, 0.67,
    0.67, 0.67, 0.67, 0.67, 0.83, 0.775, 0.775, 0.775, 0.775, 0.775,
    0.99, 0.93, 0.895, 0.895, 0.895, 0.895, 0.2, 0.23, 0.28, 0.305,
    0.335, 0.36, 0.4, 0.445, 0.485, 0.5, 0.52, 0.53, 0.6, 0.63,
    0.65, 0.66, 0.665, 0.665, 0.795, 0.8, 0.8, 0.805, 0.81, 0.815,
    0.97, 0.965, 0.955, 0.955, 0.955, 0.955, 1.15, 1.135, 1.12, 1.1,
    1.1, 1.105, 0.23, 0.275, 0.3, 0.33, 0.35, 0.37, 0.415, 0.47,
    0.5, 0.53, 0.545, 0.55, 0.63, 0.67, 0.68, 0.685, 0.69, 0.69,
    0.83, 0.835, 0.835, 0.84, 0.845, 0.85, 1.03, 1.015, 1.005, 1.0,
    1.005, 1.01, 1.25, 1.225, 1.2, 1.17, 1.165, 1.16,
])
_DUMYB = np.array([
    0.22, 0.28, 0.315, 0.345, 0.375, 0.39, 0.44, 0.5, 0.535, 0.56,
    0.57, 0.58, 0.67, 0.7, 0.72, 0.725, 0.74, 0.74, 0.88, 0.885,
    0.895, 0.9, 0.9, 0.9, 1.07, 1.07, 1.075, 1.075, 1.08, 1.08,
    1.27, 1.26, 1.26, 1.255, 1.255, 1.255, 0.25, 0.3, 0.33, 0.38,
    0.415, 0.445, 0.5, 0.56, 0.6, 0.62, 0.635, 0.64, 0.75, 0.78,
    0.8, 0.82, 0.82, 0.825, 0.98, 0.99, 1.0, 1.02, 1.02, 1.02,
    1.19, 1.2, 1.2, 1.21, 1.22, 1.225, 1.38, 1.39, 1.4, 1.41,
    1.42, 1.42, 0.34, 0.38, 0.41, 0.46, 0.5, 0.54, 0.68, 0.7,
    0.73, 0.77, 0.79, 0.84, 0.95, 0.99, 1.01, 1.05, 1.08, 1.12,
    1.2, 1.24, 1.29, 1.33, 1.37, 1.42, 1.44, 1.5, 1.55, 1.61,
    1.67, 1.72, 1.68, 1.76, 1.82, 1.89, 1.95, 2.02,
])
_DUMYC = np.array([
    0.415, 0.41, 0.4, 0.385, 0.37, 0.345, 0.5, 0.5, 0.495, 0.485,
    0.48, 0.47, 0.585, 0.58, 0.58, 0.575, 0.57, 0.56, 0.67, 0.67,
    0.67, 0.67, 0.67, 0.67, 0.75, 0.75, 0.755, 0.76, 0.765, 0.775,
    0.83, 0.84, 0.845, 0.855, 0.87, 0.895, 0.46, 0.455, 0.445, 0.42,
    0.39, 0.36, 0.575, 0.575, 0.57, 0.56, 0.545, 0.53, 0.695, 0.695,
    0.69, 0.685, 0.68, 0.665, 0.8, 0.805, 0.805, 0.81, 0.815, 0.815,
    0.92, 0.93, 0.935, 0.945, 0.97, 0.955, 1.04, 1.045, 1.05, 1.075,
    1.11, 1.105, 0.475, 0.465, 0.45, 0.43, 0.4, 0.37, 0.6, 0.6,
    0.595, 0.585, 0.575, 0.55, 0.725, 0.73, 0.73, 0.725, 0.715, 0.69,
    0.85, 0.85, 0.855, 0.865, 0.87, 0.85, 0.97, 0.975, 0.98, 1.0,
    1.02, 1.01, 1.11, 1.11, 1.11, 1.13, 1.18, 1.16,
])
_DUMYD = np.array([
    0.5, 0.49, 0.47, 0.45, 0.425, 0.39, 0.64, 0.635, 0.63, 0.62,
    0.6, 0.58, 0.77, 0.775, 0.78, 0.775, 0.765, 0.74, 0.92, 0.915,
    0.92, 0.93, 0.935, 0.9, 1.05, 1.055, 1.06, 1.08, 1.105, 1.08,
    1.195, 1.2, 1.205, 1.225, 1.265, 1.255, 0.55, 0.535, 0.525, 0.5,
    0.475, 0.445, 0.72, 0.715, 0.71, 0.69, 0.67, 0.64, 0.89, 0.89,
    0.89, 0.885, 0.87, 0.825, 1.06, 1.05, 1.05, 1.06, 1.065, 1.02,
    1.215, 1.215, 1.22, 1.245, 1.27, 1.225, 1.38, 1.38, 1.395, 1.42,
    1.47, 1.42, 0.76, 0.73, 0.7, 0.65, 0.6, 0.54, 1.0, 1.0,
    0.97, 0.93, 0.89, 0.84, 1.24, 1.23, 1.23, 1.22, 1.19, 1.12,
    1.5, 1.48, 1.48, 1.49, 1.47, 1.42, 1.75, 1.72, 1.73, 1.76,
    1.78, 1.72, 2.0, 1.97, 1.98, 2.02, 2.07, 2.02,
])

_T15258 = np.array([
    0.0, 0.28, 0.4, 0.5, 1.87,
])
_DS5258 = np.array([
    0.54, 0.54, 0.559, 0.6, 1.95,
])
_DR5258 = np.array([
    0.54, 0.54, 0.593, 0.68, 2.0,
])
_X27M = np.array([
    0.0, 1.0, 2.0, 3.0,
])
_X27I = np.array([
    1.5778, 1.67221, 1.98509, 2.28874,
])
_D13246 = np.concatenate([_DUMY1, _DUMY2, _DUMY3, _DUMY4, _DUMY5, _DUMY6])
_SUBAF = np.concatenate([_DUMYA, _DUMYB])
_SUPAF = np.concatenate([_DUMYC, _DUMYD])
_DR360A = np.concatenate([_DUMY7, _DUMY8])


def _ix(table, dep, var, lengths, lind, lower=(0, 0, 0), upper=(0, 0, 0)):
    n = len(lengths)
    return float(interx(n, table, list(var[1:n + 1]), list(lengths), dep,
                        lind=lind, lx1l=lower[0], lx2l=lower[1],
                        lx3l=lower[2], lx1u=upper[0], lx2u=upper[1],
                        lx3u=upper[2]))


def calculate_suplng(data: Mapping[str, object]) -> Dict[str, object]:
    """Translate SUPLNG: supersonic wing lift, moment and drag.

    Args:
        data: By name: ``midx`` (the Mach index), ``nalpha``, ``mach``
            (``FLC(MIDX+2)``), ``alpha`` (``FLC(23..)``), ``sw``,
            ``cbarr``, ``straight`` (``WINGIN(15)`` is ``STRA``),
            ``syna`` (the 19 ``/SYNTSS/`` words), ``a`` (the 195 ``WINGD``
            words), ``wingin`` (the ``/WINGI/`` words, 100 as the source
            declares them), ``slg`` (the 141 ``/SUPWH/`` words) and
            ``wing`` (``WING`` 1-121), as they stood; ``state`` may carry
            the saved local ``alphap``.

    Returns:
        ``a``, ``slg``, ``wing`` (1-based lists with a leading pad),
        ``detach`` (``SLG(94)``, a LOGICAL) and ``state``.

    Notes:
        Kept as executed: the shock-detachment search steps a degree
        counter against the largest angle of attack in radians, so it
        tries one degree only; the cranked wing's tip correction
        (label 1350) is read with the first axis length left at 12 from
        the leading-edge lookup; the pitch pass re-enters the
        normal-force table with ``1/BOVERT`` left by Figure
        4.1.3.2-60B for a supersonic leading edge; without a glove the
        Figure 4.1.3.2-61 abscissa divides by the stale ``CLEGLV``; and
        a cranked wing's lift ``CL`` is never formed, so its drag uses
        the words as they stood.
    """
    return _buildup(data, tail=False)


def calculate_supltg(data: Mapping[str, object]) -> Dict[str, object]:
    """Translate SUPLTG: supersonic horizontal-tail lift, moment and drag.

    SUPLTG is SUPLNG for the horizontal tail, on the tail's blocks
    (``/HTDATA/``, ``/HTI/``, ``HT`` and the second half of ``/SUPWH/``,
    ``STG``), with its own drag: the skin friction of Figure 4.1.5.1-26
    with the cutoff Reynolds number of -27, the wave drag, and the
    drag-due-to-lift factor of Figure 4.1.5.2-58.  A cranked tail also
    gets a linear normal force below 0.175 radians.

    Args:
        data: As :func:`calculate_suplng`, with ``htin`` (the 154 ``/HTI/``
            words) in place of ``wingin``, ``slg`` holding ``STG``, ``a``
            the tail's ``/HTDATA/`` words and ``wing`` the ``HT`` words;
            ``syna`` supplies the tail incidence (``SYNA(8)``) and position
            (``SYNA(6)``); also ``rl`` (``FLC(MIDX+42)``) and ``ruff``.
            ``state`` may carry the saved local ``rach``.

    Notes:
        Kept as executed, beyond SUPLNG's: a cranked tail never re-forms
        the drag-due-to-lift factor, so it uses the wave-drag ``ARG`` left
        over from the wave drag, and it keeps the ``P`` and ``DRAGC`` words
        as they stood; with zero roughness the skin friction is read at the
        Mach number the previous call left in ``RACH``; and a cranked tail
        whose inboard and outboard chords are equal takes its inboard
        friction coefficient from the word as it stood.
    """
    return _buildup(data, tail=True)


def _buildup(data, tail):
    """SUPLNG, or with ``tail`` SUPLTG."""
    a = [0.0] + [float(v) for v in data['a']]
    win = [0.0] + [float(v) for v in data['htin' if tail else 'wingin']]
    slg = [0.0] + [float(v) for v in data['slg']]
    wing = [0.0] + [float(v) for v in data['wing']]
    syna = [0.0] + [float(v) for v in data['syna']]
    alpha = [0.0] + [float(v) for v in data['alpha']]
    nalpha, midx = int(data['nalpha']), int(data['midx'])
    sw, cbarr = float(data['sw']), float(data['cbarr'])
    alphap = float(data.get('state', {}).get('alphap', 0.0))
    rach = float(data.get('state', {}).get('rach', 0.0))
    alphai = syna[8] if tail else syna[4]
    detach = bool(data.get('detach', False))
    mach = float(data['mach'])
    deltay = win[17]
    if a[86] == 0.0:
        a[86] = .00001
    if a[62] == 0.0:
        a[62] = .00001
    beta = math.sqrt(mach ** 2 - 1.)
    slg[1] = beta
    for angle_slot in range(1, nalpha + 1):
        slg[32 + angle_slot] = (alpha[angle_slot] + alphai) / RAD
    alphaj = lambda angle_slot: slg[32 + angle_slot]  # noqa: E731
    vab = [0.0] * 4

    def cn_ratio(ca, bovert, suple):
        """Labels 1030-1040: (CNA)/(CNA)T, Figure 4.1.3.2-60A or B."""
        slg[8] = deltay / ca
        if not suple:
            vab[1], vab[2] = bovert, slg[8]
            cncnt = _ix(_A1350, _DA50, vab, (9, 8), 9, upper=(0, 2, 0))
        else:
            arg = deltay / (5.85 * ca)
            slg[9] = math.atan(arg) * RAD
            vab[1], vab[2] = 1. / bovert, slg[9]
            cncnt = _ix(_B1350, _DB50, vab, (11, 8), 11, upper=(0, 2, 0))
        return 1. if cncnt > 1. else cncnt

    def basic_slope(ta, suple, lengths=(23, 8, 6)):
        """Label 1320: Figures 4.1.3.2-56A to F."""
        bcna = _ix(_T13246, _D13246, vab, lengths, 23, (0, 2, 0),
                   (0, 2, 0))
        slg[4] = bcna
        return bcna, (bcna / beta if suple else bcna / ta)

    def leading_edge(bovert):
        """Figure 4.1.3.2-62: (CLE), sharp or round."""
        vab[1] = bovert
        if win[71] == UNUSED:
            return _ix(_R13252, _DRND52, vab, (12,), 12)
        return _ix(_S13252, _DSHP52, vab, (12,), 12)

    def xac(bovert):
        """Labels 1430-1455: Figures 4.1.4.2-26A to F."""
        if bovert > 1.:
            vab[1] = 1. / bovert
            x = _ix(_T422AF, _SUPAF, vab, (6, 6, 6), 6, (0, 2, 0),
                    (0, 2, 0))
        else:
            vab[1] = bovert
            x = _ix(_T422AF, _SUBAF, vab, (6, 6, 6), 6, (0, 2, 0),
                    (0, 2, 0))
        if bovert < 0.:
            vab[1] = 1. / bovert
            x = calculate_fwdxac(vab[2], vab[3], vab[1], mach)['xac']
        if a[58] < 1.E-10 and win[midx + 71] != UNUSED:
            x = win[midx + 71] * a[122 if tail else 16] / a[10]
        return x

    def cnaaa_60a(cncnt, arg2):
        """Label 1280: Figure 4.1.3.3-60A."""
        arg1 = cncnt / arg2
        vab[2] = slg[11]
        if arg1 > 1.:
            vab[1] = 1. / arg1
            return _ix(_TL360A, _DL360A, vab, (6, 12), 12)
        vab[1] = arg1
        return _ix(_TR360A, _DR360A, vab, (13, 12), 13, (2, 0, 0))

    def normal_force(angle_slot, cna):
        aj = alphaj(angle_slot)
        cn = (cna * math.sin(2. * aj) / 2. +
              slg[12 + angle_slot] * math.sin(aj) * abs(math.sin(aj)))
        wing[60 + angle_slot] = cn * a[3] / sw
        wing[20 + angle_slot] = wing[60 + angle_slot] * math.cos(aj)

    if data['straight']:
        slg[2] = bovert = beta / a[62]
        suple = bovert > 1.0
        cncnt = cn_ratio(a[61], bovert, suple)
        slg[3] = cncnt
        if not (a[25] == 1.0 and a[58] == 0.0):
            vaa = [0.0, bovert, a[7] * a[62], a[27]]
            slg[4] = _ix(_T13246, _D13246, vaa, (23, 8, 6), 23, (0, 2, 0),
                         (0, 2, 0))
            slg[5] = slg[4] / beta if suple else slg[4] / a[62]
        elif a[7] * beta > 1.0:
            slg[4] = 4. - 2. * (1. / (a[7] * beta))
            slg[5] = slg[4] / beta
        else:
            vab[1] = a[7] * beta
            slg[5] = _ix(_G13246, _DG3246, vab, (9,), 9) * a[7]
        cna = slg[5] * cncnt
        slg[7] = cna * a[3] / sw
        slg[10] = tle192 = a[62] / 1.92
        if tle192 > 1.:
            arg = deltay / (5.85 * a[61])
            slg[9] = math.atan(arg) * RAD
            vab[1], vab[2] = slg[9], cna * tle192
            slg[12] = _ix(_A13359, _D13359, vab, (9, 11), 11, (0, 2, 0),
                          (2, 2, 0))
            slg[11] = cna * (tle192 + slg[12] * (tle192 - 1.))
        else:
            slg[11] = cna
        if not suple:
            arg = slg[11] * bovert
            for angle_slot in range(1, nalpha + 1):
                if alphaj(angle_slot) == 0.0:
                    slg[12 + angle_slot] = 0.0
                else:
                    arg1 = cncnt / (beta * abs(math.tan(alphaj(angle_slot))))
                    if arg1 > 1.:
                        vab[1], vab[2] = 1. / arg1, arg
                        slg[12 + angle_slot] = _ix(_BL3359, _DL3359, vab, (5, 10),
                                          10, (2, 2, 0))
                    else:
                        vab[1], vab[2] = arg1, arg
                        slg[12 + angle_slot] = _ix(_BR3359, _DR3359, vab, (11, 10),
                                          11, (2, 2, 0))
                normal_force(angle_slot, cna)
        else:
            alphap, detach = _supersonic_edge(
                slg, a, mach, beta, cna, cncnt, nalpha, alphaj, alphap,
                cnaaa_60a, normal_force)
    else:
        spanin = win[3] - win[2]
        slg[118] = win[5] + spanin * (a[86] - a[104])
        slg[119] = (slg[118] + win[1]) * win[3]
        slg[120] = 4. * win[3] ** 2 / slg[119]
        slg[121] = win[1] / slg[118]
        slg[2] = bovert = beta / a[86]
        suple = bovert > 1.
        vab[1], vab[2], vab[3] = bovert, slg[120] * a[86], slg[121]
        bcna, cnthry = basic_slope(a[86], suple)
        slg[5] = cnthry
        cle = leading_edge(bovert)
        slg[77] = cnthry
        slg[122] = cle
        if a[86] >= a[62]:
            slg[74] = spanin ** 2 * a[86]
            slg[73] = 4. * spanin ** 2 / slg[74]
            vab[2], vab[3] = slg[73] * a[86], 0.0
            # LGB(1) is still 12 from the leading-edge lookup.
            bcna, _ = basic_slope(a[86], suple, (12, 8, 6))
            cnt2 = bcna / a[86] if vab[1] > 1. else bcna / beta
            slg[131] = (cnthry * slg[119] / sw - cnt2 * slg[74] / sw) * \
                slg[122]
        else:
            slg[131] = slg[77] * slg[119] / sw * slg[122]
        if abs(a[58] - a[82]) < 4.:
            slg[130] = 0.0
        else:
            slg[123] = a[62] * spanin
            slg[124] = slg[123] * spanin
            slg[125] = 4. * spanin ** 2 / slg[124]
            vab[1], vab[2], vab[3] = beta / a[62], slg[125] * a[62], 0.0
            slg[2] = bovert = vab[1]
            suple = bovert > 1.
            bcna, cnthry = basic_slope(a[62], suple)
            slg[5] = cnthry
            cle = leading_edge(bovert)
            slg[132] = cle
            slg[130] = cnthry * slg[132] * slg[124] / sw
        if abs(a[76] - a[100]) < 4.:
            slg[129] = 0.0
        else:
            slg[126] = 2. * spanin
            vab[1], vab[2] = beta / a[62], a[80] / a[62]
            slg[127] = _ix(_T13253, _D13253, vab, (14, 10), 14, (2, 2, 0),
                           (2, 2, 0))
            vab[2] = a[104] / a[62]
            slg[128] = _ix(_T13253, _D13253, vab, (14, 10), 14, (2, 2, 0),
                           (2, 2, 0))
            slg[129] = (slg[127] - slg[128]) * slg[126] ** 2 / sw
        vab[1] = slg[130] / (slg[132] * beta) * sw / a[4]
        slg[133] = _ix(_T13251, _D13251, vab, (12,), 12, (2, 0, 0),
                       (2, 0, 0))
        slg[7] = slg[133] * (slg[131] + slg[130] + slg[129])
        if tail:
            for angle_slot in range(1, nalpha + 1):
                if alphaj(angle_slot) <= 0.175:
                    wing[60 + angle_slot] = slg[7] * math.sin(2. * alphaj(angle_slot)) / 2.
                    wing[20 + angle_slot] = wing[60 + angle_slot] * math.cos(alphaj(angle_slot))

    # Label 1420: pitching moment.
    dxw = (win[4] - win[3]) * a[62] * math.cos(alphai / RAD)
    a[173] = syna[1] - (syna[6 if tail else 2] + dxw)
    if data['straight']:
        vab[2], vab[3] = a[7] * a[62], a[27]
        slg[134] = xac(slg[2])
    else:
        # Labels 1500-1530: the inboard and outboard panels.
        panels = ((a[5], a[61], a[62], a[26]),
                  (a[168], a[85], a[86], a[169]))
        for panel_index, (arp, ca, ta, tapr) in enumerate(panels, 1):
            slg[2] = bovert = beta / ta
            vab[1], vab[2], vab[3] = bovert, arp * ta, tapr
            suple = bovert > 1.
            slg[134] = xacc = xac(bovert)
            if panel_index == 1:
                slg[76] = xacc
            else:
                slg[78] = (xacc * a[166] / a[10] -
                           a[23] / (a[10] * 2.) * a[86] +
                           a[23] * a[62] / a[10])
            cncnt = cn_ratio(ca, bovert, suple)
            vab[2] = arp * ta
            _, cnthry = basic_slope(ta, suple)
            slg[5] = cnthry
            if panel_index == 1:
                slg[137], slg[139] = cncnt, cnthry
                slg[92] = cncnt * cnthry
            else:
                slg[138], slg[140] = cncnt, cnthry
                slg[91] = cncnt * cnthry
        arg1 = slg[92] * a[1] * slg[76] + slg[91] * a[167] * slg[78]
        arg2 = slg[92] * a[1] + slg[91] * a[167]
        slg[134] = arg1 / arg2
    slg[135] = (a[173] / a[10] - slg[134]) * a[10] / cbarr
    slg[136] = slg[135] * slg[7]
    for angle_slot in range(1, nalpha + 1):
        if alphaj(angle_slot) <= 0.175:
            wing[40 + angle_slot] = slg[136] * alphaj(angle_slot)

    # Label 1540: drag, normal and axial force.
    if tail:
        arg, rach = _tail_drag(data, a, win, slg, mach, beta, sw, rach)
    else:
        arg = (1. + slg[82]) / (PI * a[7] * slg[82])
    for angle_slot in range(1, nalpha + 1):
        cl = wing[20 + angle_slot]
        slg[52 + angle_slot] = slg[81] * arg * sw / a[3] * cl ** 2
        wing[angle_slot] = slg[80] + slg[52 + angle_slot]
        cosa = math.cos(alpha[angle_slot] / RAD)
        sina = math.sin(alpha[angle_slot] / RAD)
        if tail or cl != UNUSED:
            wing[60 + angle_slot] = cl * cosa + wing[angle_slot] * sina
            wing[80 + angle_slot] = wing[angle_slot] * cosa - cl * sina
    if tail or slg[136] != UNUSED:
        wing[121] = slg[136] / RAD
    if tail or slg[7] != UNUSED:
        wing[101] = slg[7] / RAD
    if a[86] == 0.00001:
        a[86] = 0.0
    if a[62] == 0.00001:
        a[62] = 0.0
    return {'a': a, 'slg': slg, 'wing': wing, 'detach': detach,
            'state': {'alphap': alphap, 'rach': rach}}


def _tail_drag(data, a, win, slg, mach, beta, sr, rach):
    """SUPLTG labels 1550-1670: skin friction, wave drag and the
    drag-due-to-lift factor.  Returns ``ARG`` and ``RACH``."""
    straight = bool(data['straight'])
    rnfs, ruff = float(data['rl']), float(data['ruff'])
    spans = win[3]
    cbar = a[15]
    while True:
        slg[90] = cbar * rnfs
        if ruff != 0.0:
            arg = 12. * cbar / ruff
            rach = mach
            if rach > 3.0:
                rach = 3.0
            cept = tbfunx(_X27M, _X27I, rach, 0, 0)[0]
            slg[89] = arg ** 1.0482 * 10.0 ** cept
            if slg[89] < slg[90]:
                slg[90] = slg[89]
        slg[88] = fig26(slg[90], rach)
        if straight:
            slg[87] = slg[88] * a[3] / sr * 2.
            break
        if cbar == a[17]:
            slg[85], slg[83] = slg[90], slg[88]
            slg[87] = (slg[84] * a[1] + slg[83] * a[2]) / sr * 2.
            break
        slg[86], slg[84] = slg[90], slg[88]
        cbar = a[17]
    if straight:
        ta, s, ca = a[62], a[3], a[61]
        lerbw = win[62] * ((win[6] + win[1]) / 2.)
    else:
        ta, s, ca = a[86], slg[119], a[85]
        lerbw = win[63] * ((win[6] + win[1] + 2. * win[5]) / 4.)
    slg[2] = bovert = beta / ta
    tceff, ksharp = win[70], win[71]
    if ksharp != UNUSED:
        arg = ksharp * tceff ** 2 * s / sr
        slg[79] = arg / beta
        if bovert < 1.:
            slg[79] = arg / ta
    else:
        arg1 = 1.28 * mach ** 3 * ca ** 6 / (1. + mach ** 3 * ca ** 3)
        arg2 = 2 * lerbw * (2. * spans) / (sr * ca)
        cdle = arg1 * arg2
        arg = 16. * tceff ** 2 * s / (3. * sr)
        slg[79] = cdle + arg / beta if bovert >= 1. else cdle + arg / ta
    slg[80] = slg[87] + slg[79]
    if straight:
        rlw = win[1] + a[18] * a[10]
        slg[82] = a[3] / (rlw * 2. * spans)
        arg = (1. + slg[82]) / (PI * a[7] * slg[82])
        var = [0.0, beta * spans / rlw]
        dep = _DR5258 if ksharp == UNUSED else _DS5258
        slg[81] = _ix(_T15258, dep, var, (5,), 5, upper=(1, 0, 0))
    return arg, rach


def _supersonic_edge(slg, a, mach, beta, cna, cncnt, nalpha, alphaj, alphap,
                     cnaaa_60a, normal_force):
    """Labels 1160-1300: the straight wing with a supersonic leading
    edge, with the shock attached or detached at zero incidence."""
    slg[93] = rmach = mach * math.sqrt(1. - a[60] ** 2)
    detach = True
    if rmach >= 1.:
        detach = slg[9] / RAD > angdet(rmach)
    if detach:
        for angle_slot in range(1, nalpha + 1):
            arg2 = beta * abs(math.tan(alphaj(angle_slot)))
            slg[12 + angle_slot] = (0.0 if alphaj(angle_slot) == 0.0 else
                           cnaaa_60a(cncnt, arg2))
            normal_force(angle_slot, cna)
        return alphap, detach
    # The angle of attack at which the shock starts to detach.  The
    # counter is in degrees, the limit in radians.
    angle = 1.
    angmax = abs(alphaj(1))
    if angmax < abs(alphaj(nalpha)):
        angmax = abs(alphaj(nalpha))
    found = False
    while True:
        deg = angle / RAD
        rrmach = mach * (1. - (a[60] * math.cos(deg)) ** 2) ** .5
        if rrmach >= 1. and not (slg[9] / RAD + math.tan(deg) / a[61] <
                                 angdet(rrmach)):
            slg[115] = deg
            found = True
            break
        angle = angle + 1.
        if angle > angmax:
            break
    if found:
        btana = beta * math.tan(slg[115])
        slg[116] = cnaast = fig60b(beta, btana)
        vab = [0.0, beta]
        delcn = float(interx(1, _T61A, vab[1:], [9], _SLOPE, lind=9,
                             lx1u=2))
        slg[117] = delcn * cnaast / RAD
        alphap = slg[115] + slg[117]
        slg[13] = cnaaa_60a(cncnt, beta * abs(math.tan(alphap)))
        slg[75] = slg[13]
    else:
        slg[115] = 90. / RAD
    detang, cnaast, cnaaap = slg[115], slg[116], slg[75]
    for angle_slot in range(1, nalpha + 1):
        aj = abs(alphaj(angle_slot))
        if aj < detang:
            slg[12 + angle_slot] = fig60b(beta, beta * abs(math.tan(alphaj(angle_slot))))
        else:
            ratio = (aj - detang) / (alphap - detang)
            arg2 = beta * abs(math.tan(alphaj(angle_slot)))
            if aj < alphap:
                slg[12 + angle_slot] = cnaast - ratio * (cnaast - cnaaap)
            elif alphaj(angle_slot) == 0.0:
                slg[12 + angle_slot] = 0.0
            else:
                slg[12 + angle_slot] = cnaaa_60a(cncnt, arg2)
        normal_force(angle_slot, cna)
    return alphap, detach


def m27o33(data: Mapping[str, object]) -> Dict[str, object]:
    """M27O33: SUPLNG, then the lift-curve slope at each angle of attack.

    ``WING(J+100)`` becomes TBFUNX's slope of ``CL`` against ``ALPHA`` for
    ``J`` from 2, and ``WING(J+120)`` is marked ``-UNUSED``.  For a wing
    that is not straight-tapered the per-angle words (``J``, ``J+20``, ...,
    ``J+100`` and ``J+180``) are marked ``-UNUSED`` as well, including the
    slopes just formed, and ``WING(81)`` takes ``WING(1)``.  ``wing``
    should carry the 200 words the overlay touches.
    """
    result = calculate_suplng(data)
    wing = result['wing']
    nalpha = int(data['nalpha'])
    alpha = [float(v) for v in data['alpha'][:nalpha]]
    straight = bool(data['straight'])
    for angle_slot in range(2, nalpha + 1):
        cl = wing[21:21 + nalpha]
        wing[angle_slot + 100] = tbfunx(alpha, cl, alpha[angle_slot - 1])[1]
        wing[angle_slot + 120] = -UNUSED
        if not straight:
            for base in (0, 20, 40, 60, 80, 100, 180):
                wing[angle_slot + base] = -UNUSED
    if not straight:
        wing[81] = wing[1]
    return result


def m22o26(data: Mapping[str, object]) -> Dict[str, object]:
    """M22O26: SUPLTG, then ``HT(J+100)`` the TBFUNX slope of ``CL``
    against ``ALPHA`` and ``HT(J+120)`` marked ``-UNUSED``, for ``J`` from
    2.  (The overlay's ``KEPSLN`` argument is not read by SUPLTG.)"""
    result = calculate_supltg(data)
    ht = result['wing']
    nalpha = int(data['nalpha'])
    alpha = [float(v) for v in data['alpha'][:nalpha]]
    for angle_slot in range(2, nalpha + 1):
        ht[angle_slot + 100] = tbfunx(alpha, ht[21:21 + nalpha],
                                      alpha[angle_slot - 1])[1]
        ht[angle_slot + 120] = -UNUSED
    return result
