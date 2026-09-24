"""
CMALPH: lifting-surface zero-lift moment, aerodynamic centre, and moment
at every angle of attack.

Three parts, in the source's order:

- **Zero-lift moment** ``B(47)``: the section CMO (the mean of root and tip
  when both are given) through the Section 4.1.4.1 sweep and aspect-ratio
  factor and the Figure 4.1.4.1-6 Mach correction, plus the twist term
  from Figure 4.1.4.1-5, modified for forward sweep.
- **Aerodynamic centre** ``C(6)``, root chords aft of the root leading edge:
  the MAC quarter chord (or the input section centre) for a
  high-aspect-ratio straight surface; Figure 4.1.4.2-26A/B, or FWDXAC when
  swept forward, for a low-aspect-ratio one; and an area- and
  slope-weighted blend of the two panels for double-delta, cranked and
  curved planforms.  ``CMa`` follows as ``(x_cg - x_ac)*CLa``.
- **Moment curve** ``WING(41)`` onward: linear in CN above the aspect ratio
  ``6/A(124)``; below it the Section 4.1.4.3 nonlinear method, which moves
  the centre of pressure with angle through Figures 4.1.4.3-21 to -24 up
  to a reference angle and fairs linearly to the planform centroid at 90
  degrees.

The tables were extracted from the source DATA statements by parsing, with
``tools/fortran_data.py``, rather than by hand.

Reference: datcom-legacy/datcom_2000/cmalph.f
"""

import math
import numpy as np
from typing import Dict, Sequence
import logging

from pydatcom.aerodynamics.cdrag import STRAIGHT_TAPERED
from pydatcom.aerodynamics.fwdxac import calculate_fwdxac
from pydatcom.utils.constants import PI, RAD, UNUSED

# The source's "not available" marker, 2*UNUSED.
NOT_AVAILABLE = 2.0 * UNUSED
from pydatcom.utils.legacy_numeric import tbfunx
from pydatcom.utils.legacy_tables import tlinex, tlin3x

logger = logging.getLogger(__name__)

# Figure 4.1.4.1-6: the Mach correction to the zero-lift moment.
_XCMOM = np.array([
    0.0, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55,
    0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9,
])
_YCMOM = np.array([
    1.0, 1.0, 1.005, 1.017, 1.031, 1.05, 1.072, 1.101,
    1.132, 1.162, 1.197, 1.237, 1.287, 1.355, 1.445,
])

_Y41412 = np.array([
    0.0116, 0.01, 0.0085, 0.0072, 0.006, 0.0048, 0.0036, 0.0024, 0.0012, 0.0, -0.0012,
    -0.0024, -0.0037, -0.005, -0.0063, -0.0076, -0.009, -0.0104, -0.0118, -0.0131, -0.0145, -0.0159,
    0.0081, 0.007, 0.0061, 0.0051, 0.0042, 0.0034, 0.0025, 0.0017, 0.0008, 0.0, -0.0009,
    -0.0017, -0.0026, -0.0035, -0.0045, -0.0055, -0.0065, -0.0076, -0.0086, -0.0097, -0.0108, -0.0119,
    0.0052, 0.0045, 0.0039, 0.0033, 0.0027, 0.0022, 0.0016, 0.0011, 0.0005, 0.0, -0.0005,
    -0.0011, -0.0017, -0.0023, -0.0029, -0.0035, -0.0042, -0.0049, -0.0057, -0.0063, -0.007, -0.0078,
    0.0022, 0.0019, 0.0016, 0.0014, 0.0011, 0.0009, 0.0007, 0.0004, 0.0002, 0.0, -0.0002,
    -0.0004, -0.0007, -0.0009, -0.0012, -0.0014, -0.0017, -0.0021, -0.0024, -0.0027, -0.0032, -0.0036,
    0.0005, 0.0004, 0.0004, 0.0003, 0.0002, 0.0002, 0.0001, 0.0001, 0.0, 0.0, -0.0,
    -0.0001, -0.0001, -0.0002, -0.0002, -0.0003, -0.0004, -0.0004, -0.0005, -0.0006, -0.0007, -0.0008,
    0.0204, 0.0177, 0.0152, 0.0128, 0.0105, 0.0083, 0.0061, 0.0041, 0.002, 0.0, -0.002,
    -0.004, -0.006, -0.008, -0.01, -0.0121, -0.0141, -0.0162, -0.0182, -0.0208, -0.0232, -0.0257,
    0.0144, 0.0125, 0.0107, 0.009, 0.0073, 0.0058, 0.0043, 0.0028, 0.0014, 0.0, -0.0014,
    -0.0028, -0.0042, -0.0056, -0.0071, -0.0086, -0.0102, -0.0117, -0.0133, -0.0153, -0.0172, -0.0193,
    0.0091, 0.0078, 0.0066, 0.0055, 0.0045, 0.0035, 0.0026, 0.0017, 0.0009, 0.0, -0.0009,
    -0.0017, -0.0026, -0.0035, -0.0044, -0.0054, -0.0064, -0.0075, -0.0087, -0.0099, -0.0113, -0.0128,
    0.0036, 0.003, 0.0026, 0.0021, 0.0017, 0.0013, 0.001, 0.0007, 0.0003, 0.0, -0.0003,
    -0.0007, -0.001, -0.0013, -0.0017, -0.0021, -0.0025, -0.003, -0.0036, -0.0043, -0.005, -0.0058,
    0.0007, 0.0006, 0.0005, 0.0004, 0.0003, 0.0003, 0.0002, 0.0001, 0.0001, 0.0, -0.0001,
    -0.0001, -0.0002, -0.0003, -0.0003, -0.0004, -0.0005, -0.0006, -0.0007, -0.0008, -0.0009, -0.001,
    0.0218, 0.0189, 0.0162, 0.0135, 0.0111, 0.0087, 0.0064, 0.0042, 0.0021, 0.0, -0.002,
    -0.0041, -0.0061, -0.0081, -0.0101, -0.0121, -0.0141, -0.0162, -0.0182, -0.0202, -0.0222, -0.0243,
    0.0154, 0.0133, 0.0113, 0.0095, 0.0077, 0.006, 0.0044, 0.0029, 0.0014, 0.0, -0.0014,
    -0.0028, -0.0043, -0.0057, -0.0072, -0.0087, -0.0102, -0.0117, -0.0133, -0.0148, -0.0163, -0.0179,
    0.0097, 0.0083, 0.007, 0.0058, 0.0047, 0.0037, 0.0027, 0.0018, 0.0009, 0.0, -0.0009,
    -0.0018, -0.0027, -0.0036, -0.0045, -0.0055, -0.0065, -0.0075, -0.0086, -0.0098, -0.0108, -0.0119,
    0.0038, 0.0032, 0.0027, 0.0022, 0.0018, 0.0014, 0.001, 0.0007, 0.0003, 0.0, -0.0003,
    -0.0007, -0.001, -0.0014, -0.0018, -0.0022, -0.0026, -0.0031, -0.0036, -0.0043, -0.005, -0.0057,
    0.0008, 0.0006, 0.0005, 0.0004, 0.0003, 0.0003, 0.0002, 0.0001, 0.0001, 0.0, -0.0001,
    -0.0001, -0.0002, -0.0003, -0.0003, -0.0004, -0.0005, -0.0006, -0.0007, -0.0008, -0.0009, -0.001,
])
# Figure 4.1.4.1-5, modified for forward sweep: CMO per degree of twist.
# Y41412 is Y415A-C end to end, (22,5,3): quarter-chord sweep fastest, then
# aspect ratio (descending), then taper ratio.
_X31412 = np.array([
    0.0, 0.5, 1.0,
])
_X11412 = np.array([
    10.0, 8.0, 6.0, 3.5, 1.5,
])
_X21412 = np.array([
    -45.0, -40.0, -35.0, -30.0, -25.0, -20.0, -15.0, -10.0, -5.0, 0.0, 5.0,
    10.0, 15.0, 20.0, 25.0, 30.0, 35.0, 40.0, 45.0, 50.0, 55.0, 60.0,
])

# Figure 4.1.4.2-26A/B: aerodynamic centre of aft-swept planforms, each
# (6,6,6): the abscissa fastest, then A*tan(LE sweep) (descending), then
# taper.  26B's abscissa, beta/tan, also descends.
_X322A = np.array([
    0.0, 0.2, 0.25, 0.333, 0.5, 1.0,
])
_X122A = np.array([
    6.0, 5.0, 4.0, 3.0, 2.0, 1.0,
])
_X222A = np.array([
    0.0, 0.2, 0.4, 0.6, 0.8, 1.0,
])
_Y22A = np.array([
    0.68, 0.7, 0.72, 0.74, 0.76, 0.78,
    0.58, 0.6, 0.62, 0.64, 0.66, 0.68,
    0.51, 0.52, 0.53, 0.54, 0.55, 0.56,
    0.42, 0.425, 0.435, 0.44, 0.445, 0.445,
    0.335, 0.335, 0.335, 0.335, 0.335, 0.335,
    0.25, 0.245, 0.24, 0.235, 0.23, 0.225,
    0.87, 0.88, 0.89, 0.905, 0.915, 0.925,
    0.75, 0.76, 0.775, 0.785, 0.8, 0.81,
    0.64, 0.65, 0.66, 0.675, 0.685, 0.695,
    0.52, 0.525, 0.535, 0.54, 0.5425, 0.55,
    0.4, 0.405, 0.405, 0.405, 0.405, 0.405,
    0.28, 0.275, 0.265, 0.26, 0.255, 0.25,
    0.93, 0.94, 0.95, 0.96, 0.97, 0.98,
    0.8, 0.81, 0.82, 0.83, 0.84, 0.85,
    0.67, 0.68, 0.69, 0.7, 0.71, 0.72,
    0.54, 0.55, 0.56, 0.57, 0.575, 0.58,
    0.41, 0.41, 0.42, 0.42, 0.428, 0.43,
    0.3, 0.29, 0.28, 0.28, 0.275, 0.27,
    1.04, 1.04, 1.05, 1.05, 1.06, 1.08,
    0.89, 0.89, 0.9, 0.91, 0.91, 0.92,
    0.73, 0.74, 0.75, 0.758, 0.76, 0.77,
    0.59, 0.6, 0.6, 0.6, 0.61, 0.618,
    0.46, 0.46, 0.46, 0.46, 0.45, 0.45,
    0.32, 0.318, 0.31, 0.308, 0.3, 0.29,
    1.2, 1.2, 1.2, 1.21, 1.22, 1.23,
    1.04, 1.04, 1.05, 1.05, 1.06, 1.06,
    0.88, 0.88, 0.88, 0.882, 0.89, 0.9,
    0.7, 0.7, 0.7, 0.7, 0.7, 0.705,
    0.53, 0.53, 0.52, 0.52, 0.52, 0.52,
    0.36, 0.35, 0.34, 0.335, 0.33, 0.32,
    1.74, 1.74, 1.74, 1.74, 1.73, 1.73,
    1.5, 1.5, 1.5, 1.5, 1.5, 1.5,
    1.22, 1.22, 1.22, 1.22, 1.22, 1.22,
    1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
    0.76, 0.76, 0.76, 0.76, 0.74, 0.74,
    0.52, 0.5, 0.483, 0.48, 0.463, 0.46,
])
_X322B = np.array([
    0.0, 0.2, 0.25, 0.333, 0.5, 1.0,
])
_X122B = np.array([
    6.0, 5.0, 4.0, 3.0, 2.0, 1.0,
])
_X222B = np.array([
    1.0, 0.8, 0.6, 0.4, 0.2, 0.0,
])
_Y22B = np.array([
    0.78, 0.79, 0.805, 0.83, 0.87, 0.99,
    0.68, 0.695, 0.71, 0.725, 0.76, 0.83,
    0.56, 0.57, 0.58, 0.595, 0.62, 0.66,
    0.455, 0.46, 0.46, 0.46, 0.47, 0.5,
    0.335, 0.335, 0.335, 0.335, 0.335, 0.335,
    0.225, 0.22, 0.21, 0.2, 0.185, 0.17,
    0.925, 0.94, 0.96, 1.0, 1.05, 1.15,
    0.81, 0.82, 0.84, 0.87, 0.91, 0.97,
    0.695, 0.7, 0.71, 0.73, 0.76, 0.8,
    0.55, 0.555, 0.56, 0.57, 0.58, 0.6,
    0.405, 0.402, 0.4, 0.4, 0.4, 0.4,
    0.25, 0.245, 0.24, 0.23, 0.22, 0.2,
    0.98, 0.99, 1.01, 1.05, 1.1, 1.25,
    0.85, 0.86, 0.88, 0.9, 0.95, 1.03,
    0.72, 0.73, 0.74, 0.76, 0.79, 0.82,
    0.58, 0.59, 0.595, 0.6, 0.61, 0.63,
    0.43, 0.43, 0.43, 0.43, 0.42, 0.418,
    0.27, 0.26, 0.255, 0.25, 0.24, 0.23,
    1.08, 1.09, 1.1, 1.12, 1.16, 1.26,
    0.92, 0.93, 0.94, 0.96, 1.0, 1.07,
    0.77, 0.78, 0.79, 0.8, 0.83, 0.88,
    0.618, 0.62, 0.63, 0.64, 0.65, 0.67,
    0.45, 0.45, 0.45, 0.45, 0.45, 0.44,
    0.29, 0.28, 0.27, 0.26, 0.24, 0.22,
    1.23, 1.23, 1.24, 1.26, 1.3, 1.38,
    1.06, 1.07, 1.08, 1.09, 1.12, 1.19,
    0.9, 0.9, 0.9, 0.91, 0.94, 0.98,
    0.71, 0.71, 0.715, 0.72, 0.73, 0.75,
    0.52, 0.52, 0.52, 0.51, 0.505, 0.5,
    0.32, 0.31, 0.305, 0.3, 0.28, 0.25,
    1.72, 1.72, 1.72, 1.72, 1.7, 1.66,
    1.48, 1.48, 1.48, 1.48, 1.46, 1.42,
    1.24, 1.24, 1.24, 1.24, 1.22, 1.2,
    1.0, 1.0, 1.0, 1.0, 0.98, 0.94,
    0.74, 0.723, 0.72, 0.72, 0.7, 0.66,
    0.46, 0.44, 0.42, 0.4, 0.38, 0.32,
])

# Figures 4.1.4.3-21A (the leading-edge term) and -21B (C3).
_X211A = np.array([
    0.15, 0.3, 0.42, 0.5, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0,
])
_Y11A = np.array([
    0.7, 0.62, 0.58, 0.56, 0.55, 0.56, 0.575, 0.59, 0.6, 0.61, 0.62, 0.625,
])
_X211B = np.array([
    0.0, 0.2, 0.4, 0.5, 0.6, 0.8, 1.0,
])
_Y11B = np.array([
    0.0, 0.45, 0.85, 1.05, 1.25, 1.75, 2.25,
])

# Figure 4.1.4.3-22A: (xcp/cr) against (C3+1)*A*tan; Y11C is (15,4).
_X111C = np.array([
    1.0, 2.0, 4.0, 6.0,
])
_X211C = np.array([
    0.0, 1.0, 2.0, 2.5, 3.0, 3.5, 4.0, 5.0, 6.0, 7.0, 8.0, 10.0, 12.0, 15.0, 17.0,
])
_Y11C = np.array([
    -0.17, -0.153, -0.13, -0.112, -0.07, 0.0, 0.046, 0.115, 0.171, 0.221, 0.268, 0.37, 0.48, 0.664, 0.81,
    -0.18, -0.172, -0.157, -0.14, -0.108, -0.058, 0.0, 0.073, 0.14, 0.2, 0.257, 0.368, 0.48, 0.664, 0.81,
    -0.22, -0.21, -0.19, -0.174, -0.15, -0.11, -0.047, 0.056, 0.126, 0.187, 0.244, 0.356, 0.47, 0.653, 0.8,
    -0.24, -0.228, -0.21, -0.195, -0.173, -0.132, -0.065, 0.04, 0.11, 0.173, 0.232, 0.338, 0.45, 0.632, 0.78,
])

# Figure 4.1.4.3-22B: the stability index; Y12A is (6,11).
_X112A = np.array([
    0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0,
])
_X212A = np.array([
    0.0, 0.2, 0.4, 0.6, 0.8, 1.0,
])
_Y12A = np.array([
    5.0, 4.7, 4.15, 3.5, 2.9, 2.3,
    4.35, 3.75, 3.15, 2.5, 1.9, 1.3,
    3.35, 2.7, 2.15, 1.5, 0.999, 0.43,
    2.3, 1.7, 1.1, 0.5, -0.1, -0.7,
    1.35, 0.7, 0.1, -0.5, -1.1, -1.65,
    0.35, -0.25, -0.9, -1.5, -2.1, -2.7,
    -0.6, -1.2, -1.8, -2.45, -3.05, -3.65,
    -1.65, -2.2, -2.85, -3.5, -4.1, -4.65,
    -2.65, -3.25, -3.85, -4.45, -5.0, -5.0,
    -3.65, -4.3, -4.85, -5.0, -5.0, -5.0,
    -4.6, -5.0, -5.0, -5.0, -5.0, -5.0,
])

# Figure 4.1.4.3-23A/B: the centre-of-pressure increment by stability
# index; Y12B1 is (25,6), Y12B2 (21,3).
_X112B1 = np.array([
    0.2, 0.3, 0.4, 0.6, 0.8, 1.0,
])
_X212B1 = np.array([
    4.0, 3.5, 3.0, 2.5, 2.0, 1.75, 1.5, 1.25, 1.0, 0.75, 0.5, 0.25, 0.0,
    -0.25, -0.5, -0.75, -1.0, -1.25, -1.5, -1.75, -2.0, -2.5, -3.0, -3.5, -4.0,
])
_Y12B1 = np.array([
    0.0, -0.08, -0.148, -0.203, -0.244, -0.26, -0.271, -0.277, -0.278, -0.27, -0.252, -0.175, 0.025, 0.2, 0.305, 0.38, 0.418, 0.44, 0.449, 0.451, 0.44, 0.4, 0.347, 0.284, 0.2,
    0.07, 0.0, -0.053, -0.12, -0.17, -0.193, -0.212, -0.227, -0.239, -0.241, -0.22, -0.16, 0.05, 0.35, 0.43, 0.452, 0.454, 0.438, 0.4, 0.325, 0.215, 0.038, -0.056, -0.115, -0.15,
    0.13, 0.062, -0.003, -0.064, -0.118, -0.14, -0.16, -0.178, -0.19, -0.193, -0.186, -0.175, 0.04, 0.2, 0.27, 0.312, 0.328, 0.322, 0.275, 0.17, 0.0, -0.2, -0.303, -0.376, -0.43,
    0.07, 0.058, 0.04, 0.012, -0.018, -0.035, -0.052, -0.067, -0.08, -0.088, -0.089, -0.067, -0.014, 0.045, 0.095, 0.124, 0.132, 0.116, 0.073, 0.0, -0.093, -0.18, -0.21, -0.23, -0.25,
    0.13, 0.08, 0.06, 0.045, 0.03, 0.0275, 0.025, 0.0225, 0.02, 0.0175, 0.015, 0.0125, 0.01, 0.0075, 0.005, 0.0025, 0.0, -0.01, -0.02, -0.03, -0.04, -0.05, -0.06, -0.07, -0.08,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
])
_X112B2 = np.array([
    0.6, 0.8, 1.0,
])
_X212B2 = np.array([
    4.0, 3.5, 3.0, 2.75, 2.5, 2.25, 2.0, 1.5, 1.0, 0.5, 0.25,
    0.0, -0.25, -0.5, -1.0, -1.5, -2.0, -2.5, -3.0, -3.5, -4.0,
])
_Y12B2 = np.array([
    -0.12, -0.2, -0.27, -0.285, -0.29, -0.287, -0.268, -0.178, -0.072, 0.0, 0.025, 0.04, 0.048, 0.05, 0.037, 0.011, -0.01, -0.02, -0.03, -0.039, -0.04,
    -0.085, -0.11, -0.13, -0.135, -0.135, -0.13, -0.12, -0.09, -0.052, -0.01, 0.017, 0.04, 0.058, 0.068, 0.065, 0.05, 0.03, 0.02, 0.029, 0.03, 0.02,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
])

# Figure 4.1.4.3-24A: the aspect-ratio index; Y13A is (15,9).
_X213A = np.array([
    0.0, 0.1, 0.15, 0.2, 0.25, 0.3, 0.4, 0.45, 0.5, 0.55, 0.6, 0.7, 0.8, 0.9, 1.0,
])
_X113A = np.array([
    0.25, 0.5, 1.0, 1.5, 2.0, 3.0, 4.0, 5.0, 6.0,
])
_Y13A = np.array([
    0.0, 0.17, 0.235, 0.28, 0.297, 0.282, 0.242, 0.215, 0.21, 0.227, 0.262, 0.378, 0.56, 0.8, 1.04,
    0.14, 0.24, 0.276, 0.3, 0.31, 0.298, 0.242, 0.215, 0.193, 0.193, 0.218, 0.287, 0.383, 0.5, 0.63,
    0.14, 0.24, 0.264, 0.28, 0.29, 0.277, 0.197, 0.132, 0.082, 0.045, 0.03, 0.02, 0.03, 0.063, 0.1,
    0.1, 0.2, 0.228, 0.24, 0.228, 0.2, 0.097, 0.022, -0.05, -0.095, -0.12, -0.13, -0.128, -0.125, -0.12,
    0.0, 0.14, 0.155, 0.154, 0.138, 0.105, 0.0, -0.085, -0.162, -0.195, -0.214, -0.235, -0.24, -0.24, -0.23,
    -0.1, -0.022, -0.002, 0.0, -0.02, -0.058, -0.195, -0.257, -0.3, -0.336, -0.362, -0.388, -0.4, -0.4, -0.39,
    -0.24, -0.2, -0.19, -0.198, -0.223, -0.26, -0.355, -0.392, -0.418, -0.442, -0.46, -0.487, -0.5, -0.507, -0.51,
    -0.38, -0.37, -0.382, -0.397, -0.417, -0.436, -0.483, -0.505, -0.526, -0.54, -0.552, -0.57, -0.58, -0.582, -0.58,
    -0.5, -0.513, -0.522, -0.53, -0.54, -0.55, -0.568, -0.577, -0.585, -0.59, -0.6, -0.612, -0.62, -0.625, -0.63,
])

# Figure 4.1.4.3-24B-1/-2: its centre-of-pressure increment; Y13B1 is
# (9,6), Y13B2 (9,3).
_X213B1 = np.array([
    1.0, 0.8, 0.6, 0.4, 0.2, 0.0, -0.2, -0.4, -0.6,
])
_X113B1 = np.array([
    0.2, 0.3, 0.4, 0.6, 0.8, 1.0,
])
_Y13B1 = np.array([
    1.0, 0.8, 0.6, 0.4, 0.2, 0.0, -0.2, -0.4, -0.6,
    0.64, 0.482, 0.35, 0.23, 0.12, 0.0, -0.16, -0.36, -0.59,
    0.38, 0.26, 0.15, 0.07, 0.02, 0.0, -0.1, -0.23, -0.42,
    0.12, 0.048, 0.0, -0.032, -0.035, -0.03, -0.12, -0.31, -0.56,
    0.05, -0.02, -0.05, -0.06, -0.058, -0.05, -0.068, -0.11, -0.17,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
])
_X213B2 = np.array([
    1.0, 0.8, 0.6, 0.4, 0.2, 0.0, -0.2, -0.4, -0.6,
])
_X113B2 = np.array([
    0.6, 0.8, 1.0,
])
_Y13B2 = np.array([
    0.3, 0.29, 0.29, 0.28, 0.26, 0.2, 0.23, 0.25, 0.28,
    0.14, 0.14, 0.14, 0.14, 0.14, 0.1, 0.11, 0.13, 0.14,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
])



def _two(flat, fastest: int, count: int) -> np.ndarray:
    """``count`` runs of ``fastest`` values, as TLINEX's (x2, x1) layout."""
    return np.asarray(flat, dtype=float).reshape(count, fastest).T


def _three(flat, n2: int, n1: int, n3: int) -> np.ndarray:
    """Source storage as TLIN3X's (x2, x1, x3) layout, x2 fastest."""
    return np.asarray(flat, dtype=float).reshape(n3, n1, n2).transpose(2, 1, 0)


_FIG_41415 = _three(_Y41412, 22, 5, 3)
_FIG_4226A = _three(_Y22A, 6, 6, 6)
_FIG_4226B = _three(_Y22B, 6, 6, 6)
_FIG_4322A = _two(_Y11C, 15, 4)
_FIG_4322B = _two(_Y12A, 6, 11)
_FIG_4323A = _two(_Y12B1, 25, 6)
_FIG_4323B = _two(_Y12B2, 21, 3)
_FIG_4324A = _two(_Y13A, 15, 9)
_FIG_4324B1 = _two(_Y13B1, 9, 6)
_FIG_4324B2 = _two(_Y13B2, 9, 3)


def _div(a: float, b: float) -> float:
    """IEEE division: the source forms reciprocals it may never use."""
    with np.errstate(divide='ignore', invalid='ignore'):
        return float(np.float64(a) / np.float64(b))


def _fig26(curve: float, ratio: float, taper: float) -> float:
    """Figure 4.1.4.2-26: part A for tan/beta below one, else part B at
    beta/tan.  ``ratio`` is tan(LE sweep)/beta."""
    if ratio >= 1.0:
        return float(tlin3x(_X122B, _X222B, _X322B, _FIG_4226B, curve,
                            _div(1.0, ratio), taper, 2, 0, 0, 2, 0, 0))
    return float(tlin3x(_X122A, _X222A, _X322A, _FIG_4226A, curve, ratio,
                        taper, 2, 0, 0, 2, 0, 0))


def calculate_cmalph(planform_type: float,
                     alpha_deg: Sequence[float],
                     geometry: Dict[str, float],
                     section: Dict[str, float],
                     lift: Dict[str, object],
                     flight: Dict[str, float],
                     sref: float, cbarr: float) -> Dict[str, object]:
    """Translate CMALPH: zero-lift moment, aerodynamic centre and CM curve.

    Args:
        planform_type: ``WINGIN(15)``; one of the ``cdrag`` WTYPE constants.
        alpha_deg: ``B(23)`` onward, the surface's local angles.
        geometry: ``A`` entries ``a1``, ``a2``, ``a3``, ``a5``, ``a7``,
            ``a10``, ``a16``, ``a21``, ``a23``, ``a26``, ``a27``, ``a30``,
            ``a34``, ``a37``, ``a38``, ``a40``, ``a43``, ``a62``, ``a86``,
            ``a124``, ``a125``, ``a133``, ``a164``, ``a166``, ``a167``,
            ``a168``, ``a169``, ``a171``, ``a172``, ``a173``.  Only those a
            path reads are needed.
        section: ``WINGIN`` entries ``cmo`` (61), ``cmot`` (67),
            ``twista`` (11), ``deltay`` (17), and ``xac`` (``WINGIN(I+71)``,
            the input section aerodynamic centre at this Mach, or UNUSED).
        lift: ``cla`` (``WING(101)``), ``cn`` (``WING(61)`` onward) and
            ``alpha_clmax`` (``B(43)``).
        flight: ``mach`` (``B(1)``) and ``beta`` (``B(2)``).
        sref, cbarr: ``SREF``, ``CBARR``.

    Returns:
        Dictionary with ``cm`` (``WING(41)`` onward), ``cma``
        (``WING(121)``), ``cm0`` (``B(47)``), ``xac`` (``C(6)``), the
        ``C`` work array as ``c`` (one-based keys), ``a170``, the ``a38``
        and ``a62`` the routine leaves in COMMON, and ``nonlinear``.

    Notes:
        A zero leading-edge tangent is replaced in COMMON, ``A(38)`` by
        ``1e-4`` on the straight path and ``A(62)`` by ``1e-15`` on the
        others, and stays replaced for every later routine.  Both are
        returned.

        On the non-straight path the outboard panel is sent to FWDXAC when
        swept forward only if the inboard panel was too: the ``A(86) < 0``
        test sits on the inboard FWDXAC branch, and the two aft-swept
        inboard branches jump past it.  An aft-swept inboard panel with a
        forward-swept outboard panel therefore reads Figure 4.1.4.2-26A at
        a negative abscissa, clamped to zero.  Kept, and flagged.
    """
    kind = float(planform_type)
    g = {k: float(v) for k, v in geometry.items()}
    alpha = np.asarray(alpha_deg, dtype=float)
    cn = np.asarray(lift['cn'], dtype=float)
    beta, mach = float(flight['beta']), float(flight['mach'])
    c = {}
    result = {'method': 'legacy_cmalph'}

    a170 = g['a173'] / g['a10']
    c[48] = g['a10'] / cbarr
    c[1], c[2] = float(section['cmo']), float(section['cmot'])
    if not (abs(c[1]) < 1.0e-10 or abs(c[2]) < 1.0e-10):
        c[1] = 0.50 * (c[1] + c[2])
    calm, _ = tbfunx(_XCMOM, _YCMOM, mach, 0, 0)
    calm = g['a3'] * g['a16'] * calm / (sref * cbarr)
    c[5] = (g['a7'] * g['a43']**2 / (g['a7'] + 2.0 * g['a43'])) * c[1] * calm
    twist = float(section['twista'])
    if abs(twist) >= 1.0e-20:
        c[4] = float(tlin3x(_X11412, _X21412, _X31412, _FIG_41415, g['a7'],
                            g['a40'], g['a27'], 2, 2, 0, 2, 2, 0))
        c[5] = c[5] + c[4] * twist * calm
    c[3] = calm

    c[9] = g['a7'] * g['a38']
    if kind == STRAIGHT_TAPERED:
        if g['a38'] == 0.0:
            g['a38'] = 0.0001
        c[10] = g['a38'] / beta
        c[11] = beta / g['a38']
        if g['a7'] > g['a125']:
            c[6] = (g['a30'] - 0.25 * g['a16']) / g['a10']
            xac = float(section.get('xac', UNUSED))
            if xac != UNUSED:
                c[6] = c[6] + g['a16'] * (xac - 0.25) / g['a10']
        elif g['a38'] < 0.0:
            fwd = calculate_fwdxac(c[9], g['a27'], c[10], mach)
            c[6] = fwd['xac']
            result['fwdxac'] = fwd
        else:
            c[6] = _fig26(c[9], c[10], g['a27'])
    else:
        if g['a62'] == 0.0:
            g['a62'] = 1.0e-15
        c[12] = g['a5'] * g['a62']
        c[13] = g['a62'] / beta
        c[14] = _div(1.0, c[13])
        c[15] = g['a168'] * g['a86']
        c[16] = g['a86'] / beta
        c[17] = _div(1.0, c[16])
        if g['a62'] < 0.0:
            c[18] = calculate_fwdxac(c[12], g['a26'], c[13], mach)['xac']
            outboard_forward = g['a86'] < 0.0
        else:
            c[18] = _fig26(c[12], c[13], g['a26'])
            outboard_forward = False
            if g['a86'] < 0.0:
                result['source_defect'] = \
                    'forward_swept_outboard_panel_reads_figure_26a'
        if outboard_forward:
            c[19] = calculate_fwdxac(c[15], g['a169'], c[16], mach)['xac']
        else:
            c[19] = _fig26(c[15], c[16], g['a169'])
        c[20] = (c[19] * g['a166'] / g['a10'] - g['a164'] * g['a86'] / g['a10']
                 + g['a23'] * g['a62'] / g['a10'])
        c[6] = ((g['a171'] * g['a1'] * c[18] + g['a172'] * g['a167'] * c[20]) /
                (g['a171'] * g['a1'] + g['a172'] * g['a167']))

    c[7] = (a170 - c[6]) * (g['a10'] / cbarr)
    c[8] = c[7] * float(lift['cla'])
    result.update({'cm0': c[5], 'xac': c[6], 'cma': c[8], 'a170': a170,
                   'a38': g['a38'], 'a62': g.get('a62'), 'c': c})

    if g['a7'] > 6.0 / g['a124'] or g['a34'] < 0.0:
        result.update({'cm': cn * c[7] + c[5], 'nonlinear': False})
        return result

    # Section 4.1.4.3, the low-aspect-ratio nonlinear method.
    c[21] = g['a7'] * (1.0 + g['a27']) * g['a38'] / 4.0
    c[22] = (g['a27'] + c[21] + ((1.0 + c[21] * g['a27']) /
                                 (1.0 + g['a27']))) / 3.0
    if not (g['a7'] >= 6.0 / g['a124'] or kind == STRAIGHT_TAPERED):
        ybsi = g['a23'] * (1.0 + 2.0 * g['a26']) / (3.0 * (1.0 + g['a26']))
        xble = (ybsi * g['a1'] * g['a62'] +
                (g['a23'] * g['a62'] + (g['a133'] - g['a21']) * g['a86']) *
                g['a2']) / g['a3']
        c[22] = (xble + g['a16'] / 2.0) / g['a10']
    c[23], _ = tbfunx(_X211B, _Y11B, g['a27'], 0, 0)
    c[24] = (c[23] + 1.0) * c[9]
    c[25] = float(tlinex(_X111C, _X211C, _FIG_4322A, g['a7'], c[24],
                         2, 0, 2, 1))
    c[26], _ = tbfunx(_X211A, _Y11A, float(section['deltay']), 2, 2)
    c[27] = c[25] + c[26]
    alpha_clmax = float(lift['alpha_clmax'])
    c[28] = math.sin(alpha_clmax / RAD)
    c[29] = math.tan(alpha_clmax / RAD)
    c[30] = c[27] / c[28] - c[6] / c[29]

    cm = []
    for m, a in enumerate(alpha):
        c[31] = math.sin(a / RAD)
        c[32] = math.cos(a / RAD)
        c[33] = c[31] / c[32]
        c[34] = g['a7'] * g['a37']
        c[35] = abs(c[33] / c[29])
        tanrat = c[33] > c[29]
        if m == 0:
            # The reference angle and its data, on the first pass only.
            c[36] = math.atan(math.tan(alpha_clmax / RAD) / 0.6)
            temp2 = c[29] / math.tan(c[36])
            c[47] = temp2
            c[40] = float(tlinex(_X112A, _X212A, _FIG_4322B, c[9], g['a27'],
                                 1, 0, 1, 0))
            c[37] = float(tlinex(_X113A, _X213A, _FIG_4324A, c[34], g['a27'],
                                 0, 0, 2, 0))
            c[51] = float(tlinex(_X113B2, _X213B2, _FIG_4324B2, temp2, c[37],
                                 2, 0, 2, 0))
            c[50] = float(tlinex(_X112B2, _X212B2, _FIG_4323B, temp2, c[40],
                                 2, 0, 2, 0))
            c[49] = (c[6] * math.cos(c[36]) +
                     math.sin(c[36]) * (c[30] + c[51] + c[50]))
        if a / RAD > c[36]:
            c[41] = c[22] - c[49]
            c[42] = PI / 2.0 - c[36]
            c[43] = c[41] / c[42]
            c[44] = c[49] + c[43] * (a / RAD - c[36])
            cm.append(cn[m] * (a170 - c[44]) * (g['a10'] / cbarr) + c[5])
            continue
        if tanrat:
            c[46] = c[29] / c[33]
            c[38] = float(tlinex(_X113B2, _X213B2, _FIG_4324B2, c[46], c[37],
                                 2, 0, 2, 0))
            c[39] = float(tlinex(_X112B2, _X212B2, _FIG_4323B, c[46], c[40],
                                 2, 0, 2, 0))
        else:
            c[38] = float(tlinex(_X113B1, _X213B1, _FIG_4324B1, c[35], c[37],
                                 2, 0, 2, 0))
            c[39] = float(tlinex(_X112B1, _X212B1, _FIG_4323A, c[35], c[40],
                                 2, 0, 2, 0))
        c[44] = c[6] * c[32] + (c[30] + c[39] + c[38]) * c[31]
        cm.append(cn[m] * (g['a173'] / g['a10'] - c[44]) * c[48] + c[5])
    result.update({'cm': np.array(cm), 'nonlinear': True})
    return result


def calculate_cacalc(alpha_deg: Sequence[float], cd: Sequence[float],
                     cl: Sequence[float]) -> Dict[str, np.ndarray]:
    """Translate CACALC: normal and axial force at the surface's angles."""
    alpha = np.asarray(alpha_deg, dtype=float) / RAD
    cd, cl = np.asarray(cd, dtype=float), np.asarray(cl, dtype=float)
    return {'cn': cl * np.cos(alpha) + cd * np.sin(alpha),
            'ca': cd * np.cos(alpha) - cl * np.sin(alpha)}


# The deviation from the linear lift slope, percent, past which M31O37 and
# M33O41 mark the moment unavailable.
_WING_DEVIATION = 90.0
_TAIL_DEVIATION = 7.5


def calculate_moment_overlay(surface: str,
                             free_alpha_deg: Sequence[float],
                             cd: Sequence[float], cl: Sequence[float],
                             cmalph_inputs: Dict[str, object],
                             experimental: bool = False) -> Dict[str, object]:
    """Translate overlays M31O37 (wing) and M33O41 (horizontal tail).

    CMALPH, then CACALC, then the surface's normal and axial force and its
    lift- and moment-curve slopes over the free-stream angles, with the
    moment marked unavailable (``2*UNUSED``) from the first angle at which
    the lift slope departs from linear by more than the overlay's limit.

    Args:
        surface: ``'wing'`` or ``'tail'``.
        free_alpha_deg: ``FLC(23)`` onward.
        cd, cl: The surface's ``(1)`` and ``(21)`` onward drag and lift.
        cmalph_inputs: :func:`calculate_cmalph`'s arguments.
        experimental: ``KWING`` or ``KHT``: experimental data supplied, in
            which case the slopes are retaken at zero angle and the
            linearity test is skipped.

    Returns:
        Dictionary with ``cm``, ``cn``, ``ca``, ``cla``, ``cma`` (``(41)``,
        ``(61)``, ``(81)``, ``(101)``, ``(121)`` onward) and ``cmalph``.

    Notes:
        The two overlays differ in three ways, all kept.  The wing's test
        allows 90 percent deviation where its comment says 15; the tail's
        allows 7.5.  The wing measures deviation against CMALPH's analytic
        ``CLa``, the tail against the TBFUNX slope at the first angle,
        which by then has overwritten ``HT(101)``.  And the tail's restore
        of the analytic ``CLa`` and ``CMa`` sits inside the branch that is
        skipped for a low-aspect-ratio straight tail or experimental data,
        so such a tail keeps the TBFUNX slopes at its first angle.
    """
    if surface not in ('wing', 'tail'):
        raise ValueError("surface must be 'wing' or 'tail'")
    cm_result = calculate_cmalph(**cmalph_inputs)
    alpha = np.asarray(free_alpha_deg, dtype=float)
    cd, cl = np.asarray(cd, dtype=float), np.asarray(cl, dtype=float)
    cm = np.array(cm_result['cm'], dtype=float)
    cla0 = float(cmalph_inputs['lift']['cla'])
    cma0 = cm_result['cma']

    forces = calculate_cacalc(alpha, cd, cl)
    cla = np.array([tbfunx(alpha, cl, a, 0, 0)[1] for a in alpha])
    cma = np.array([tbfunx(alpha, cm, a, 0, 0)[1] for a in alpha])
    if experimental:
        cla0 = tbfunx(alpha, cl, 0.0, 0, 0)[1]
        cma0 = tbfunx(alpha, cm, 0.0, 0, 0)[1]

    g = cmalph_inputs['geometry']
    skip = ((float(g['a7']) <= 6.0 / float(g['a124']) and
             float(cmalph_inputs['planform_type']) == STRAIGHT_TAPERED)
            or experimental)
    limit = _WING_DEVIATION if surface == 'wing' else _TAIL_DEVIATION
    reference = cla0 if surface == 'wing' else cla[0]
    if not skip:
        flagged = False
        for j in range(1, len(alpha)):
            if 100.0 * abs(cla[j] / reference - 1.0) > limit:
                flagged = True
            if flagged:
                cm[j] = NOT_AVAILABLE
                cma[j] = NOT_AVAILABLE
    if surface == 'wing' or not skip:
        cla[0], cma[0] = cla0, cma0
    return {'cm': cm, 'cn': forces['cn'], 'ca': forces['ca'], 'cla': cla,
            'cma': cma, 'cmalph': cm_result,
            'method': f'legacy_m{"31o37" if surface == "wing" else "33o41"}'}
