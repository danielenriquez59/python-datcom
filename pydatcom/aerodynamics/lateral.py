"""
Subsonic lateral-directional stability: overlays M29O35 and M17O21, and
SUBLAT.

- ``M29O35`` sets up the geometry in ``/SBETA/``: the wing's height on the
  body, the body's side projection at its quarter and three-quarter
  lengths, the vertical panels' arms, and the effective dihedral.
- ``SUBLAT`` forms the sideslip derivatives CY_beta, Cn_beta and Cl_beta
  (``(141)``, ``(161)`` and ``(181)`` onward of each block) for a lifting
  surface alone, with the body, and for the vertical panels.
- ``M17O21`` runs SUBLAT for the wing with the vertical tail and again for
  the horizontal tail with the ventral fin, then sums the configurations.

SUBLAT reads about sixty COMMON slots from nine blocks.  Its inputs are
therefore taken as mappings keyed by the source's own one-based indices,
``a[120]`` for ``A(120)``, so each line can be checked against the FORTRAN
directly; the argument docs say what the slots are.

The tables were extracted from the source DATA statements by parsing, with
``tools/fortran_data.py``, rather than by hand.

Reference: datcom-legacy/datcom_2000/m29o35.f, sublat.f, m17o21.f
"""

import math
import numpy as np
from typing import Dict, Mapping, Optional, Sequence
import logging

from pydatcom.aerodynamics.cdrag import STRAIGHT_TAPERED
from pydatcom.interactions.body_vortex import getmax
from pydatcom.utils.constants import PI, RAD, UNUSED
from pydatcom.utils.legacy_numeric import tbfunx, trapz
from pydatcom.utils.legacy_tables import tlinex, tlin3x, tlin4x

logger = logging.getLogger(__name__)

# Figure 5.1.2.1-27: wing CL_beta/CL from sweep.  Y27 is (8,5,3): half-chord
# sweep fastest, then aspect ratio, then taper.
_X127 = np.array([
    1.0, 2.0, 4.0, 6.0, 8.0,
])
_X227 = np.array([
    -20.0, 0.0, 20.0, 30.0, 40.0, 50.0, 55.0, 60.0,
])
_X327 = np.array([
    0.0, 0.5, 1.0,
])
_Y27 = np.array([
    0.0014, 0.0, -0.00125, -0.002, -0.0027, -0.0036, -0.004, -0.0044,
    0.0015, 0.0, -0.00145, -0.0022, -0.003, -0.0041, -0.005, -0.00595,
    0.0016, 0.0, -0.0016, -0.0024, -0.0033, -0.0047, -0.0057, -0.0071,
    0.0016, 0.0, -0.0016, -0.0024, -0.0035, -0.0049, -0.006, -0.0074,
    0.0016, 0.0, -0.0016, -0.0027, -0.0035, -0.0049, -0.006, -0.0074,
    0.0012, 0.0, -0.0012, -0.0019, -0.0026, -0.0034, -0.0039, -0.0044,
    0.0013, 0.0, -0.0013, -0.0021, -0.003, -0.0043, -0.00515, -0.0064,
    0.0015, 0.0, -0.0014, -0.0024, -0.0036, -0.005, -0.00605, -0.0075,
    0.00165, 0.0, -0.0016, -0.0025, -0.0038, -0.0054, -0.0066, -0.0082,
    0.0018, 0.0, -0.00175, -0.0027, -0.004, -0.0058, -0.007, -0.0089,
    0.00105, 0.0, -0.001, -0.0016, -0.0023, -0.003, -0.0035, -0.0038,
    0.0012, 0.0, -0.0013, -0.0021, -0.0031, -0.00435, -0.00505, -0.0062,
    0.0014, 0.0, -0.00165, -0.00245, -0.0036, -0.0052, -0.0061, -0.0078,
    0.00167, 0.0, -0.0017, -0.0028, -0.004, -0.00595, -0.00715, -0.009,
    0.0018, 0.0, -0.0018, -0.00295, -0.0042, -0.0062, -0.0078, -0.01,
])

# Figure 5.1.2.1-28A/B: the compressibility and aspect-ratio terms.
_X128A = np.array([
    2.0, 3.0, 4.0, 5.0, 6.0, 8.0, 10.0,
])
_X228A = np.array([
    0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95,
])
_Y28A = np.array([
    1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.995, 0.99,
    1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.01, 1.03, 1.03, 1.02,
    1.0, 1.0, 1.01, 1.015, 1.025, 1.05, 1.08, 1.09, 1.1, 1.1,
    1.0, 1.01, 1.015, 1.02, 1.05, 1.09, 1.115, 1.16, 1.2, 1.21,
    1.0, 1.01, 1.02, 1.04, 1.07, 1.12, 1.17, 1.24, 1.32, 1.36,
    1.0, 1.01, 1.05, 1.07, 1.12, 1.18, 1.27, 1.4, 1.58, 1.7,
    1.0, 1.02, 1.05, 1.1, 1.15, 1.23, 1.37, 1.54, 1.84, 2.08,
])
_X128B = np.array([
    0.0, 0.5, 1.0,
])
_X228B = np.array([
    1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 6.0, 8.0,
])
_Y28B = np.array([
    -0.0058, -0.00345, -0.00235, -0.00145, -0.001, -0.00045, -0.00025, 5e-05, 0.0004,
    -0.008, -0.00555, -0.004, -0.003, -0.00235, -0.0014, -0.001, -0.00065, -0.0002,
    -0.0113, -0.008, -0.00595, -0.00465, -0.0037, -0.00255, -0.00182, -0.00147, -0.00097,
])

# Figure 5.1.2.1-29: the dihedral effect, (9,3,3).
_X129 = np.array([
    0.0, 40.0, 60.0,
])
_X229 = np.array([
    0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 8.0, 10.0,
])
_X329 = np.array([
    0.0, 0.5, 1.0,
])
_Y29 = np.array([
    0.0, -5.2e-05, -8.8e-05, -0.00011, -0.000134, -0.000153, -0.000168, -0.00019, -0.0002,
    0.0, -4.8e-05, -8.5e-05, -0.000108, -0.000128, -0.000141, -0.000153, -0.000173, -0.000178,
    0.0, -4e-05, -7.3e-05, -9.5e-05, -0.000108, -0.000119, -0.000127, -0.000135, -0.000138,
    0.0, -5.2e-05, -9.8e-05, -0.000132, -0.000162, -0.000186, -0.000208, -0.00024, -0.00026,
    0.0, -5e-05, -9.6e-05, -0.000124, -0.00015, -0.00017, -0.000188, -0.000217, -0.00023,
    0.0, -5e-05, -8.7e-05, -0.000111, -0.000129, -0.000142, -0.000153, -0.000166, -0.00017,
    0.0, -5e-05, -9.6e-05, -0.000133, -0.000167, -0.000193, -0.000216, -0.000252, -0.00028,
    0.0, -5e-05, -9.5e-05, -0.000129, -0.000155, -0.000178, -0.000197, -0.000225, -0.000245,
    0.0, -5e-05, -8.8e-05, -0.000113, -0.000132, -0.000147, -0.000159, -0.000172, -0.00018,
])

# Figure 5.1.2.1-30A/B: the dihedral compressibility and twist terms.
_X130A = np.array([
    2.0, 4.0, 6.0, 8.0, 10.0,
])
_X230A = np.array([
    0.0, 0.2, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95,
])
_Y30A = np.array([
    1.0, 1.01, 1.018, 1.02, 1.023, 1.03, 1.04, 1.05, 1.057,
    1.0, 1.012, 1.03, 1.045, 1.06, 1.085, 1.118, 1.16, 1.19,
    1.0, 1.015, 1.045, 1.07, 1.1, 1.14, 1.197, 1.27, 1.33,
    1.0, 1.018, 1.05, 1.085, 1.125, 1.19, 1.26, 1.39, 1.485,
    1.0, 1.02, 1.058, 1.097, 1.148, 1.215, 1.325, 1.495, 1.635,
])
_X130B = np.array([
    0.0, 0.4, 0.6, 1.0,
])
_X230B = np.array([
    3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0,
])
_Y30B = np.array([
    -1.92e-05, -2.22e-05, -2.38e-05, -2.31e-05, -2.3e-05, -2.41e-05, -2.6e-05, -2.84e-05, -3.28e-05,
    -2.2e-05, -2.87e-05, -3.23e-05, -3.35e-05, -3.39e-05, -3.42e-05, -3.5e-05, -3.7e-05, -4.2e-05,
    -2.33e-05, -3e-05, -3.35e-05, -3.5e-05, -3.66e-05, -3.7e-05, -3.75e-05, -4e-05, -4.7e-05,
    -2.33e-05, -3e-05, -3.35e-05, -3.5e-05, -3.66e-05, -3.7e-05, -3.75e-05, -4e-05, -4.7e-05,
])

# Figure 5.1.2.1-31: non-uniform dihedral, Y31 = T31A-C end to end,
# (11,4,3,3): span station fastest, then sweep, aspect ratio, taper.
_X131 = np.array([
    -40.0, 0.0, 40.0, 60.0,
])
_X231 = np.array([
    0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0,
])
_X331 = np.array([
    2.0, 4.0, 8.0,
])
_X431 = np.array([
    0.0, 0.5, 1.0,
])
_Y31 = np.array([
    0.0, -0.0001, -0.00026, -0.00056, -0.001, -0.0015, -0.00217, -0.0028, -0.00336, -0.0039, -0.00435,
    0.0, -0.0001, -0.00032, -0.0006, -0.0011, -0.00176, -0.0025, -0.0032, -0.00385, -0.00438, -0.0048,
    0.0, -0.0001, -0.0003, -0.00063, -0.00117, -0.0018, -0.00256, -0.0033, -0.00395, -0.00448, -0.00496,
    0.0, -0.0001, -0.0003, -0.0006, -0.00115, -0.0017, -0.00236, -0.0031, -0.0037, -0.00417, -0.0045,
    0.0, -0.0001, -0.00025, -0.00078, -0.00153, -0.0024, -0.00335, -0.0043, -0.0051, -0.00586, -0.0066,
    0.0, -0.0001, -0.00038, -0.001, -0.00183, -0.0029, -0.00405, -0.0052, -0.00633, -0.00726, -0.0079,
    0.0, -0.0001, -0.0004, -0.00102, -0.00186, -0.00303, -0.00425, -0.0054, -0.00655, -0.00736, -0.00775,
    0.0, -0.0001, -0.0004, -0.001, -0.00177, -0.00275, -0.0038, -0.00476, -0.00565, -0.00635, -0.00653,
    0.0, -8e-05, -0.00057, -0.00145, -0.00255, -0.00374, -0.00484, -0.00585, -0.00682, -0.0077, -0.00845,
    0.0, -0.0001, -0.00063, -0.00165, -0.003, -0.00445, -0.00625, -0.00786, -0.0093, -0.0104, -0.0111,
    0.0, -0.00013, -0.0009, -0.00195, -0.00325, -0.00483, -0.00655, -0.00815, -0.0095, -0.01035, -0.0106,
    0.0, -0.00013, -0.0009, -0.00195, -0.00325, -0.0046, -0.00585, -0.00697, -0.00787, -0.0084, -0.00842,
    0.0, -0.00015, -0.0003, -0.00065, -0.0011, -0.0017, -0.0025, -0.00325, -0.004, -0.00453, -0.005,
    0.0, -0.00015, -0.00037, -0.00076, -0.00128, -0.0019, -0.00276, -0.00355, -0.00426, -0.00493, -0.00536,
    0.0, -0.00015, -0.00037, -0.00076, -0.00128, -0.0019, -0.00276, -0.00355, -0.00426, -0.0048, -0.0052,
    0.0, -0.00014, -0.00036, -0.0007, -0.00125, -0.0019, -0.0027, -0.0034, -0.004, -0.00452, -0.00508,
    0.0, -0.0001, -0.0004, -0.00096, -0.00163, -0.00255, -0.0038, -0.0051, -0.00637, -0.00743, -0.00825,
    0.0, -0.0002, -0.00065, -0.00125, -0.00215, -0.0032, -0.0045, -0.00595, -0.0072, -0.00845, -0.00935,
    0.0, -0.0002, -0.00065, -0.00125, -0.00215, -0.0032, -0.0045, -0.00595, -0.0072, -0.00823, -0.009,
    0.0, -0.0001, -0.0005, -0.00113, -0.00185, -0.00285, -0.004, -0.0051, -0.00605, -0.0069, -0.00745,
    0.0, -0.00015, -0.0006, -0.00125, -0.00215, -0.00327, -0.0047, -0.00645, -0.0084, -0.01, -0.0114,
    0.0, -0.0003, -0.00094, -0.00192, -0.00315, -0.0046, -0.0068, -0.0089, -0.0109, -0.01255, -0.0141,
    0.0, -0.0003, -0.0009, -0.00195, -0.00325, -0.00482, -0.00675, -0.00855, -0.01015, -0.01155, -0.01265,
    0.0, -0.00015, -0.0009, -0.00176, -0.00285, -0.00405, -0.0053, -0.00665, -0.0078, -0.00885, -0.0097,
    0.0, -0.00015, -0.00038, -0.0007, -0.00126, -0.0019, -0.00268, -0.00344, -0.00418, -0.00484, -0.0053,
    0.0, -0.00016, -0.0004, -0.00076, -0.0013, -0.002, -0.00288, -0.0037, -0.00445, -0.0051, -0.0056,
    0.0, -0.00015, -0.00035, -0.0007, -0.00125, -0.0019, -0.00278, -0.00356, -0.0043, -0.005, -0.0054,
    0.0, -0.00012, -0.0003, -0.00067, -0.0012, -0.0018, -0.00257, -0.0033, -0.00396, -0.00458, -0.00505,
    0.0, -0.00015, -0.0005, -0.00115, -0.0018, -0.00273, -0.00385, -0.00515, -0.0066, -0.0079, -0.0088,
    0.0, -0.00015, -0.00052, -0.0012, -0.0021, -0.0033, -0.00486, -0.00627, -0.00756, -0.0087, -0.00972,
    0.0, -0.00015, -0.00052, -0.0012, -0.0021, -0.0033, -0.00465, -0.006, -0.00728, -0.00836, -0.00925,
    0.0, -0.00015, -0.0005, -0.00116, -0.0019, -0.00285, -0.004, -0.0051, -0.0061, -0.00705, -0.0079,
    -0.0, -0.00018, -0.0006, -0.00125, -0.0021, -0.00335, -0.00505, -0.0071, -0.00915, -0.011, -0.01255,
    0.0, -0.00025, -0.0009, -0.0018, -0.0031, -0.0048, -0.007, -0.0094, -0.01145, -0.0133, -0.01475,
    0.0, -0.0003, -0.0009, -0.00195, -0.00327, -0.0048, -0.00697, -0.0089, -0.01057, -0.01205, -0.0132,
    0.0, -0.00012, -0.00073, -0.00155, -0.00265, -0.004, -0.00536, -0.0067, -0.00795, -0.009, -0.01,
])

# Figure 5.2.2.1-26: the body effect on dihedral; Y526 is (9,7).
_X1526 = np.array([
    4.0, 4.5, 5.0, 5.5, 6.0, 7.0, 8.0,
])
_X2526 = np.array([
    0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6,
])
_Y526 = np.array([
    1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.99, 0.97,
    1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.98, 0.948, 0.911,
    1.0, 1.0, 1.0, 1.0, 0.997, 0.971, 0.933, 0.883, 0.827,
    1.0, 1.0, 1.0, 0.991, 0.963, 0.922, 0.87, 0.811, 0.746,
    1.0, 1.0, 0.995, 0.97, 0.932, 0.884, 0.829, 0.764, 0.695,
    1.0, 1.0, 0.977, 0.944, 0.899, 0.845, 0.78, 0.715, 0.641,
    1.0, 0.985, 0.96, 0.921, 0.87, 0.812, 0.745, 0.67, 0.592,
])

# Figure 5.2.3.1-8A-C: the body yawing moment, chained.
_X158A = np.array([
    20.0, 14.0, 10.0, 8.0, 7.0, 6.0, 5.0, 4.0, 3.0, 2.5,
])
_X258A = np.array([
    0.2, 0.8,
])
_Y58A = np.array([
    0.1, 1.88,
    0.4, 2.21,
    0.74, 2.6,
    0.98, 2.8,
    1.3, 3.13,
    1.61, 3.5,
    2.0, 3.88,
    2.5, 4.4,
    2.99, 5.0,
    3.45, 5.4,
])
_X158B = np.array([
    0.8, 1.0, 1.2, 1.4, 1.6,
])
_X258B = np.array([
    0.0, 3.0, 6.0,
])
_Y58B = np.array([
    0.0, 2.35, 4.68,
    0.0, 3.0, 6.0,
    0.0, 3.6, 7.25,
    0.0, 4.18, 8.5,
    0.0, 4.79, 9.5,
])
_X158C = np.array([
    0.5, 0.6, 0.8, 1.0, 2.0,
])
_X258C = np.array([
    0.0, 6.0,
])
_Y58C = np.array([
    -0.00048, 0.00251,
    -0.00048, 0.0035,
    -0.00048, 0.00477,
    -0.00048, 0.00559,
    -0.00048, 0.00641,
])

# Figure 5.3.1.1-22A-D: vertical-tail effective aspect ratio.
_X122A = np.array([
    1.0, 0.6,
])
_X222A = np.array([
    0.0, 0.125, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 1.75, 2.0,
    2.25, 2.5, 3.0, 3.25, 3.5, 3.75, 4.0, 5.0, 7.0,
])
_Y22A = np.array([
    0.0, 0.4, 0.72, 0.99, 1.19, 1.32, 1.4, 1.46, 1.5, 1.51, 1.48, 1.42, 1.27, 1.21, 1.17, 1.13, 1.1, 1.04, 1.02,
    0.0, 0.7, 0.94, 1.18, 1.35, 1.46, 1.54, 1.6, 1.63, 1.64, 1.6, 1.53, 1.36, 1.28, 1.21, 1.16, 1.13, 1.06, 1.02,
])
_X122B = np.array([
    0.5, 0.6, 0.7, 0.8,
])
_X222B = np.array([
    0.0, -0.2, -0.3, -0.4, -0.5, -0.6, -0.65, -0.7, -0.75, -0.8, -0.85, -0.9, -1.0,
])
_Y22B = np.array([
    1.05, 0.94, 0.9, 0.87, 0.86, 0.87, 0.9, 0.93, 0.98, 1.06, 1.16, 1.29, 1.7,
    1.15, 1.0, 0.95, 0.9, 0.89, 0.9, 0.92, 0.96, 1.01, 1.08, 1.18, 1.31, 1.7,
    1.22, 1.05, 0.99, 0.94, 0.92, 0.92, 0.95, 0.98, 1.03, 1.1, 1.2, 1.33, 1.7,
    1.29, 1.09, 1.02, 0.97, 0.94, 0.94, 0.96, 1.0, 1.06, 1.12, 1.22, 1.36, 1.7,
])
_X5322C = np.array([
    0.0, 0.2, 0.4, 0.6, 0.7, 0.8, 0.9, 1.2, 1.4, 1.6, 2.0,
])
_Y5322C = np.array([
    0.0, 0.29, 0.52, 0.7, 0.77, 0.83, 0.87, 0.98, 1.04, 1.07, 1.13,
])
_X5322D = np.array([
    0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0,
])
_Y5322D = np.array([
    0.75, 0.75, 0.75, 0.75, 0.75, 0.835, 0.92, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
])

# Figure 5.3.1.1-24A-C: twin vertical panels.
_X5324A = np.array([
    0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0,
])
_Y5324A = np.array([
    1.5, 1.33, 1.197, 1.098, 1.03, 1.0, 1.03, 1.098, 1.197, 1.33, 1.5,
])
_X124B = np.array([
    20.0, 0.0,
])
_X224B = np.array([
    0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0,
])
_Y24B = np.array([
    0.9, 1.64, 2.19, 2.63, 2.95, 3.22, 3.44, 3.61, 3.75, 3.86,
    1.0, 1.72, 2.3, 2.75, 3.11, 3.39, 3.62, 3.82, 3.98, 4.11,
])
_X124C = np.array([
    0.2, 0.4, 0.6, 0.8, 1.0,
])
_X224C = np.array([
    0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0,
])
_Y24C = np.array([
    1.0, 0.9, 0.818, 0.741, 0.68, 0.631, 0.593, 0.563, 0.542, 0.527, 0.51,
    1.0, 0.93, 0.87, 0.825, 0.784, 0.752, 0.73, 0.71, 0.696, 0.688, 0.68,
    1.0, 0.95, 0.911, 0.882, 0.859, 0.844, 0.83, 0.82, 0.812, 0.805, 0.8,
    1.0, 0.968, 0.941, 0.923, 0.908, 0.898, 0.89, 0.883, 0.881, 0.88, 0.879,
    1.0, 0.979, 0.962, 0.954, 0.949, 0.943, 0.938, 0.935, 0.935, 0.931, 0.931,
])



def _two(flat, fastest: int, count: int) -> np.ndarray:
    """``count`` runs of ``fastest`` values, as TLINEX's (x2, x1) layout."""
    return np.asarray(flat, dtype=float).reshape(count, fastest).T


def _three(flat, n2: int, n1: int, n3: int) -> np.ndarray:
    return np.asarray(flat, dtype=float).reshape(n3, n1, n2).transpose(2, 1, 0)


_FIG27 = _three(_Y27, 8, 5, 3)
_FIG28A = _two(_Y28A, 10, 7)
_FIG28B = _two(_Y28B, 9, 3)
_FIG29 = _three(_Y29, 9, 3, 3)
_FIG30A = _two(_Y30A, 9, 5)
_FIG30B = _two(_Y30B, 9, 4)
_FIG31 = np.asarray(_Y31, dtype=float).reshape(3, 3, 4, 11).transpose(3, 2, 1, 0)
_FIG526 = _two(_Y526, 9, 7)
_FIG58A = _two(_Y58A, 2, 10)
_FIG58B = _two(_Y58B, 3, 5)
_FIG58C = _two(_Y58C, 2, 5)
_FIG22A = _two(_Y22A, 19, 2)
_FIG22B = _two(_Y22B, 13, 4)
_FIG24B = _two(_Y24B, 10, 2)
_FIG24C = _two(_Y24C, 11, 5)


def _div(dividend, divisor):
    """IEEE division, matching the source's unset-handling propagation."""
    with np.errstate(divide='ignore', invalid='ignore'):
        return float(np.float64(dividend) / np.float64(divisor))


def calculate_sublat(stb: Mapping[int, float],
                     surface: Dict[str, Sequence[float]],
                     win: Mapping[int, float], a: Mapping[int, float],
                     b: Mapping[int, float],
                     combination_cl: Sequence[float], body_cla: float,
                     aht: Mapping[int, float], avt: Mapping[int, float],
                     vtin: Mapping[int, float], vt_cla: float,
                     htin: Mapping[int, float], tvtin: Sequence[float],
                     syna: Mapping[int, float], flight: Dict[str, object],
                     bd1: float, sref: float, cbarr: float, blref: float,
                     flags: Dict[str, bool], ity: int) -> Dict[str, object]:
    """Translate SUBLAT: sideslip derivatives of one lifting surface and its
    vertical panels.

    Args:
        stb: ``/SBETA/`` as M29O35 and any earlier pass left it, one-based.
            SUBLAT reads 1-4, 9-12, 59-63, 72, 122 and 123, and on the
            straight low-aspect-ratio path a stale 96.
        surface: The surface-alone curves ``cl``, ``cm``, ``cn``
            (``(21)``, ``(41)``, ``(61)`` onward).
        win: The surface's input block (``WINGIN``): 2, 3, 4, 11, 12, 13,
            14 and ``'type'`` for 15.
        a: Its ``A`` block: 1, 3, 4, 25, 42, 43, 44, 49, 68, 70, 73, 94,
            97, 118, 120, 131, 163, 167, 168, 169, 171, 172.
        b: Its ``B`` block: 2 and the Mach-zero lift 3 onward (CLMCH0).
        combination_cl: The surface-body lift, ``BW(21)`` onward.
        body_cla: ``BODY(101)``.
        aht: 2, 16, 30, 62, 119 of the horizontal tail's ``A`` block.
        avt: 4, 50, 118, 120 of the vertical panel's ``A`` block.
        vtin: 3 and 4 of the vertical panel's input block.
        vt_cla: ``VTIN(I+20)``, its section lift slope at this Mach.
        htin: 3 and 4 of the horizontal tail's input block.
        tvtin: ``TVTIN(1)`` to ``(8)``, the twin-panel inputs.
        syna: ``/SYNTSS/``: 1, 7, 8, 18, 19 (``XCG``, ``ZH``, ``ALIH``,
            ``PHIV``, ``PHIF``).
        flight: ``mach`` (``FLC(I+2)``), ``reynolds_per_length``
            (``FLC(I+42)``) and ``alpha`` (``FLC(23)`` onward).
        bd1: ``BD(1)``, the body length.
        sref, cbarr, blref: ``/OPTION/``.
        flags: ``wgpl``, ``bo``, ``vtpl``, ``tvtpan``, ``htpl``,
            ``transn``.  For M17O21's second call these are the tail's.
        ity: 0 for the wing pass, 1 for the horizontal tail.

    Returns:
        Dictionary with ``surface`` (``cyb``, ``cnb``, ``clb`` arrays: the
        ``(141)``, ``(161)``, ``(181)`` onward of the surface block),
        ``combination`` (``cyb``, ``cnb`` scalars and ``clb`` array: the
        ``BW`` ones), ``vertical`` (the vertical panel's, scalars and an
        array), and the updated ``stb``.  Blocks the source does not reach
        are absent.

    Notes:
        With no body the source jumps from the wing-alone results straight
        to its return, so the vertical panel is not computed; and on a
        transonic pass it returns after the combination's CY_beta and
        Cn_beta.

        On the double-delta and cranked path, the two Figure 5.1.2.1-27
        results go to ``STB(132)`` and ``STB(128)``, but the absolute value
        taken for a negative half-chord sweep is applied to ``YA27``, the
        straight path's variable.  A forward-swept panel keeps its sign.
        Kept.

        The twin-panel correction runs when either ``PHIV`` or ``PHIF`` is
        set and takes the one for this pass, so a set ``PHIF`` alone scales
        the vertical tail by ``2*cos(UNUSED)**2 = 2``.  Kept.
    """
    stb = {int(k): float(v) for k, v in stb.items()}
    alpha = np.asarray(flight['alpha'], dtype=float)
    alpha_count = len(alpha)
    mach = float(flight['mach'])
    straight = float(win['type']) == STRAIGHT_TAPERED
    result: Dict[str, object] = {'method': 'legacy_sublat'}
    ya = {'27': UNUSED, '28a': UNUSED, '28b': UNUSED, '29': UNUSED,
          '30a': UNUSED, '30b': UNUSED, '31i': UNUSED, '31o': UNUSED}

    if flags['wgpl']:
        double_span = 2.0 * win[4]
        if not flags['transn']:
            cl = np.asarray(surface['cl'], dtype=float)
            cm = np.asarray(surface['cm'], dtype=float)
            cn = np.asarray(surface['cn'], dtype=float)
            stb16_num = 6.0 * a[44] * a[42]
            stb16_den = PI * a[120] * (a[120] + 4.0 * a[43])
            mach_beta = math.sqrt(1.0 - (mach * a[43])**2)
            cyb36_num = a[120] + 4.0 * a[43]
            cyb36_den = a[120] * mach_beta + 4.0 * a[43]
            cnb_a = 1.0 / (4.0 * PI * a[120])
            cnb_b = abs(a[44]) / (PI * a[120] * (a[120] + 4.0 * a[43]))
            cnb_c = a[43] - 0.5 * a[120] - a[120]**2 / (8.0 * a[43])
            cnb_e = (a[120] + 4.0 * a[43]) / (a[120] * mach_beta + 4.0 * a[43])
            cnb_num = ((a[120] * mach_beta)**2 + 4.0 * a[120] * mach_beta * a[43] -
                       8.0 * a[43] * a[43])
            cnb_den = (a[120] * a[120] + 4.0 * a[120] * a[43] -
                       8.0 * a[43] * a[43])
            cyb, cnb = np.zeros(alpha_count), np.zeros(alpha_count)
            for j in range(alpha_count):
                cl_at_mach0 = b[j + 3]
                stb[j + 16] = ((1.0 / RAD) * cl_at_mach0**2 * (stb16_num / stb16_den) -
                               .0001 * abs(stb[122]))
                if cl_at_mach0 == 0.0:
                    stb[j + 36] = 0.0
                else:
                    stb[j + 36] = (cyb36_num / cyb36_den) * (stb[j + 16] / cl_at_mach0)
                    cyb[j] = stb[j + 36] * cl[j]
                cm_arm_term = (0.0 if cn[j] == 0.0 else
                               6.0 * (cm[j] / cn[j] / cbarr) * (abs(a[42]) / a[120]))
                stb[j + 76] = (cnb_a - cnb_b * (cnb_c + cm_arm_term)) / RAD
                cnb[j] = (cnb_e * (cnb_num / cnb_den) * stb[j + 76] * cl[j]**2 *
                          double_span * sref / (blref * a[3]))

            if straight and a[120] < 1.0:
                clb = ((-2.0 / (RAD * 3.0 * a[120]) * cl -
                        stb[122] * a[120] / (6.0 * RAD**2) * (a[4] / sref))
                       * double_span / blref)
            elif straight:
                mach_a73 = mach * a[73]
                x1arg = abs(a[70])
                ar_over_a73 = a[120] / a[73]
                ya['27'] = tlin3x(_X127, _X227, _X327, _FIG27, a[120], x1arg,
                                  a[118], 0, 1, 0, 2, 2, 0)
                if a[70] < 0.0:
                    ya['27'] = abs(ya['27'])
                ya['28a'] = tlinex(_X128A, _X228A, _FIG28A, ar_over_a73, mach_a73,
                                   0, 2, 1, 2)
                ya['28b'] = tlinex(_X128B, _X228B, _FIG28B, a[118], a[120],
                                   0, 0, 0, 0)
                ya['30b'] = tlinex(_X130B, _X230B, _FIG30B, a[118], a[120],
                                   0, 2, 0, 2)
                for j in range(alpha_count):
                    stb[j + 96] = ((ya['27'] * ya['28a'] + ya['28b']) * cl[j] *
                                   sref / a[4] + win[11] * a[68] * ya['30b'])
                base = np.array([stb[j + 96] for j in range(alpha_count)])
                if win[12] != UNUSED and win[12] != 0.0:
                    stb[7] = a[131] / (2.0 * PI) * RAD
                    x1arg = math.atan(a[68] / b[2]) * RAD
                    x3arg = b[2] * a[120] / stb[7]
                    k_over_b2 = stb[7] / b[2]
                    ya['31i'] = tlin4x(_X131, _X231, _X331, _X431, _FIG31,
                                       x1arg, stb[2], x3arg, a[118],
                                       0, 0, 1, 0, 2, 2, 1, 0)
                    ya['31o'] = tlin4x(_X131, _X231, _X331, _X431, _FIG31,
                                       x1arg, stb[3], x3arg, a[118],
                                       0, 0, 1, 0, 2, 2, 1, 0)
                    clb = ((base + (ya['31o'] - ya['31i']) * k_over_b2 * win[14] /
                            RAD + ya['31i'] / RAD * k_over_b2 * win[13]) *
                           a[4] / sref * double_span / blref)
                else:
                    ya['30a'] = tlinex(_X130A, _X230A, _FIG30A, ar_over_a73, mach_a73,
                                       2, 2, 2, 2)
                    ya['29'] = tlin3x(_X129, _X229, _X329, _FIG29, x1arg,
                                      a[120], a[118], 0, 0, 0, 2, 2, 0)
                    clb = ((base + win[13] * ya['29'] * ya['30a']) * a[4] /
                           sref * double_span / blref)
            else:
                # Double delta or cranked: the two panels separately.
                stb[135] = -2.0 / (3.0 * RAD * a[163])
                stb[132] = tlin3x(_X127, _X227, _X327, _FIG27, a[163],
                                  abs(a[70]), a[25], 0, 1, 0, 2, 2, 0)
                if a[70] < 0.0:
                    ya['27'] = abs(ya['27'])
                stb[134] = tlinex(_X128A, _X228A, _FIG28A, a[163] / a[73],
                                  a[73] * mach, 0, 2, 1, 2)
                stb[133] = tlinex(_X128B, _X228B, _FIG28B, a[25], a[163],
                                  0, 0, 0, 0)
                stb[131] = -2.0 / (3.0 * RAD * a[168])
                stb[128] = tlin3x(_X127, _X227, _X327, _FIG27, a[168],
                                  abs(a[94]), a[169], 0, 1, 0, 2, 2, 0)
                if a[94] < 0.0:
                    ya['27'] = abs(ya['27'])
                stb[130] = tlinex(_X128A, _X228A, _FIG28A, a[168] / a[97],
                                  a[97] * mach, 0, 2, 1, 2)
                stb[129] = tlinex(_X128B, _X228B, _FIG28B, a[169], a[168],
                                  0, 0, 0, 0)
                stb[58] = a[171] * (a[1] / sref) + a[172] * (a[167] / sref)
                inner = (stb[135] if a[163] < 1.0 else
                         stb[132] * stb[134] + stb[133])
                stb[127] = (a[171] * (a[1] / sref) * inner *
                            (win[3] - win[2]) / win[4])
                outer = (stb[131] if a[168] < 1.0 else
                         stb[128] * stb[130] + stb[129])
                stb[126] = a[172] * (a[167] / sref) * outer
                clb = (cl / stb[58]) * (stb[127] + stb[126]) * double_span / blref

            result['surface'] = {'cyb': cyb, 'cnb': cnb,
                                 'clb': np.asarray(clb, dtype=float)}
            for slot, key in [(68, '27'), (71, '28a'), (70, '28b'),
                              (67, '29'), (66, '30a'), (69, '30b'),
                              (64, '31i'), (65, '31o')]:
                stb[slot] = float(ya[key])
            if not flags['bo']:
                result['stb'] = stb
                return result

        # The surface with the body.
        arg = 2.0 * stb[1] / stb[72]
        stb[57] = 1.0 + .49 * arg if arg >= 0.0 else 1.0 - .85 * arg
        cyb_bw = -stb[57] * body_cla - 0.0001 * abs(stb[122])
        stb[56] = 1.0 + math.log(1.0e-6 * float(flight['reynolds_per_length'])
                                 * bd1) / 4.86
        fig58a = tlinex(_X158A, _X258A, _FIG58A, bd1**2 / stb[62],
                        syna[1] / bd1, 2, 1, 2, 1)
        fig58b = tlinex(_X158B, _X258B, _FIG58B,
                        math.sqrt(stb[61] / stb[60]), fig58a, 2, 0, 2, 1)
        stb[15] = tlinex(_X158C, _X258C, _FIG58C, stb[59] / stb[123],
                         fig58b, 2, 0, 2, 1)
        cnb_bw = -stb[15] * stb[56] * (stb[62] * bd1) / (sref * blref)
        combination = {'cyb': float(cyb_bw), 'cnb': float(cnb_bw)}
        result['combination'] = combination
        if flags['transn']:
            result['stb'] = stb
            return result

        stb[14] = (1.2 * stb[4] * stb[72] * math.sqrt(a[120]) /
                   (RAD * 2.0 * win[4]**2))
        stb[13] = (-.0005 * math.sqrt(a[120]) * stb[122] *
                   (stb[72] / (2.0 * win[4]))**2)
        stb[8] = tlinex(_X1526, _X2526, _FIG526, a[120] / a[49],
                        stb[63] / (2.0 * win[4]), 0, 0, 0, 1)
        if not straight:
            clb_scale = (stb[8] * (stb[132] * stb[134] + stb[128] * stb[130]) +
                         stb[133] + stb[129])
            clb_fpart = 0.0
            clb_epart = stb[13]
        else:
            clb_scale = stb[68] * stb[71] * stb[8] + stb[70]
            clb_fpart = win[11] * a[68] * stb[69]
            clb_epart = (result['surface']['clb'][0] * blref / double_span -
                         (stb[96] - stb[13]) * a[4] / sref)
        combination['clb'] = ((np.asarray(combination_cl, dtype=float) *
                               clb_scale + clb_epart + (stb[14] + clb_fpart) * a[4] /
                               sref) * double_span / blref)

    # The vertical panels.
    if not (flags['vtpl'] or flags['tvtpan']):
        result['stb'] = stb
        return result
    cyb_v, cnb_v, clb_v = 0.0, 0.0, np.zeros(alpha_count)
    cos_alpha, sin_alpha = np.cos(alpha / RAD), np.sin(alpha / RAD)
    if flags['vtpl']:
        vt_span_ratio = vtin[4] / (2.0 * (vtin[4] - vtin[3]))
        stb[120] = tlinex(_X122A, _X222A, _FIG22A, avt[118], vt_span_ratio,
                          0, 0, 0, 0)
        stb[118] = tbfunx(_X5322D, _Y5322D, vt_span_ratio, 0, 0)[0]
        tail_height = syna[7] - ((htin[4] - htin[3]) * aht[62] + aht[30] -
                                 aht[16] / 4.0) * math.sin(syna[8] / RAD)
        if not flags['htpl'] or tail_height < 0.0:
            stb[121] = 0.0
            stb[119] = 0.0
        else:
            stb[121] = tlinex(_X122B, _X222B, _FIG22B, stb[9] / stb[10],
                              -tail_height / vtin[4], 0, 0, 2, 0)
            stb[119] = tbfunx(_X5322C, _Y5322C,
                              (aht[2] + aht[119]) / avt[4], 0, 0)[0]
        stb[116] = stb[120] * avt[120] * (1.0 + stb[119] * (stb[121] - 1.0))
        stb[117] = (0.724 + 3.06 * ((avt[4] / a[4]) / (1.0 + a[43])) +
                    .4 * stb[4] / stb[59] + .009 * a[120])
        cla_ratio_sq = stb[116]**2 / (vt_cla * RAD / (2.0 * PI))**2
        sweep_factor = 1.0 + avt[50]**2 / b[2]**2
        stb[5] = 2.0 * PI * stb[116] / (2.0 + math.sqrt(cla_ratio_sq * sweep_factor + 4.0))
        cyb_v = -stb[118] * stb[5] * stb[117] * avt[4] / (RAD * sref)
        cnb_v = -cyb_v * stb[11] / blref
        clb_v = cyb_v * (stb[12] * cos_alpha - stb[11] * sin_alpha) / blref
        if not (syna[18] == UNUSED and syna[19] == UNUSED):
            angle = syna[18] if ity == 0 else syna[19]
            twin = 2.0 * math.cos(angle / RAD)**2
            cyb_v, cnb_v, clb_v = twin * cyb_v, twin * cnb_v, twin * clb_v
            result['twin_correction'] = twin
    if flags['tvtpan']:
        stb[75] = tbfunx(_X5324A, _Y5324A, tvtin[0] / tvtin[1], 0, 0)[0]
        stb[6] = tvtin[1]**2 / tvtin[4]
        stb[116] = stb[75] * stb[6]
        stb[73] = tlinex(_X124B, _X224B, _FIG24B, tvtin[5], stb[116],
                         0, 2, 0, 2)
        stb[74] = tlinex(_X124C, _X224C, _FIG24C, tvtin[3] / bd1,
                         tvtin[2] / tvtin[1], 2, 0, 0, 0)
        tvt_cyb = -stb[74] * stb[73] * 2.0 * tvtin[4] / (sref * RAD)
        cyb_v = cyb_v + tvt_cyb
        cnb_v = cnb_v - tvt_cyb * tvtin[6] / blref
        clb_v = clb_v + tvt_cyb * (tvtin[7] * cos_alpha - tvtin[6] * sin_alpha) / blref
    result['vertical'] = {'cyb': float(cyb_v), 'cnb': float(cnb_v),
                          'clb': np.asarray(clb_v, dtype=float)}
    result['stb'] = stb
    return result


def calculate_m29o35_block(win: Mapping[int, float], a: Mapping[int, float],
                           avt: Mapping[int, float], vtin: Mapping[int, float],
                           aht: Mapping[int, float], syna: Mapping[int, float],
                           body: Dict[str, Sequence[float]], bd1: float,
                           bd66: float, a10: float, xv: float,
                           vertical_up: bool, surface_on: bool,
                           body_on: bool, panel_on: bool, htpl: bool,
                           tail_pass: bool,
                           stb123_wing: float = 0.0) -> Dict[str, object]:
    """One of M29O35's two geometry passes, into ``STB`` or ``STBH``.

    The wing pass takes the wing's blocks, ``XV`` (``SYNA(9)``) and the
    vertical tail; the tail pass takes the horizontal tail's, ``XVF``
    (``SYNA(12)``) and the ventral fin, with the panel's up/down sense
    reversed.

    Args:
        win: The surface's input block: 1, 2, 3, 4, 6, 12, 13, 14.
        a: The surface's ``A`` block: 23, 62, 86.
        avt: The panel's ``A`` block: 62, 86, 122, 136, 195.
        vtin: The panel's input block: 1, 2, 4, 5, 6, ``'type'``.
        aht: The horizontal tail's ``A``: 161.
        syna: ``/SYNTSS/`` 1-8.
        body: ``x``, ``r``, and ``zu``/``zl`` when given (else None).
        bd1, bd66: ``BD(1)`` and ``BD(66)``.
        a10: ``A(10)``, which both passes read from the *wing's* block.
        xv: ``SYNA(9)`` or ``SYNA(12)``.
        vertical_up: ``VERTUP``.
        surface_on, body_on, panel_on, htpl: The configuration flags.
        tail_pass: True for the ``STBH`` pass.
        stb123_wing: ``STB(123)`` after the wing pass, which the tail pass
            copies over its own.

    Returns:
        Dictionary with ``stb`` (one-based) and the dihedral inputs as
        M29O35 leaves them (``12``, ``13``, ``14``).

    Notes:
        The tail pass forms ``STBH(1)`` with the wing's ``A(10)``, the
        exposed root chord, where every other term is the tail's own.  The
        parallel slot is ``AHT(10)``.  Kept.
    """
    stb: Dict[int, float] = {}
    win = dict(win)
    if tail_pass and not surface_on:
        return {'stb': stb, 'dihedral': {}}
    stb[2] = (win[4] - win[12]) / win[4]
    stb[3] = 1.0
    if not surface_on:
        # Without a wing the wing pass goes straight to the tail pass.
        return {'stb': stb, 'dihedral': {}}
    height = syna[7] if tail_pass else syna[3]
    incidence = syna[8] if tail_pass else syna[4]
    x_surface = syna[6] if tail_pass else syna[2]
    if surface_on and body_on:
        sin_incidence = math.sin(incidence / RAD)
        stb[1] = (-height + (.25 * a10 + bd66) * sin_incidence -
                  (win[4] - win[3]) * math.tan(win[13] / RAD))
        stb[4] = -height + win[6] / 4.0 * sin_incidence
        stb[72] = (win[4] - win[3]) * 2.0
        stb[63] = (x_surface + bd66 + win[2] * a[86] + a[23] * a[62] +
                   win[1] / 2.0)
        x = np.asarray(body['x'], dtype=float)
        r = np.asarray(body['r'], dtype=float)
        if body.get('zu') is not None:
            zu = np.asarray(body['zu'], dtype=float)
            zl = np.asarray(body['zl'], dtype=float)
            stb[62] = trapz(zu, x, 1)[0] - trapz(zl, x, 1)[0]
            depth = zu - zl
            stb[61] = tbfunx(x, depth, bd1 * .25, 0, 0)[0]
            stb[60] = tbfunx(x, depth, bd1 * .75, 0, 0)[0]
            stb[59] = getmax(x, depth)[1]
        else:
            stb[62] = 2.0 * trapz(r, x, 1)[0]
            stb[61] = 2.0 * tbfunx(x, r, bd1 * .25, 0, 0)[0]
            stb[60] = 2.0 * tbfunx(x, r, bd1 * .75, 0, 0)[0]
            stb[59] = 2.0 * getmax(x, r)[1]
        stb[123] = 2.0 * getmax(x, r)[1]
        if panel_on:
            stb[11] = (xv - syna[1]) + avt[195] + avt[122] / 4.0
            panel_points_up = vertical_up if not tail_pass else not vertical_up
            if panel_points_up:
                stb[12] = -syna[5] + avt[136]
            else:
                stb[12] = -syna[5] - avt[136]
            if htpl and syna[7] >= 0.0:
                tail_arm_x = syna[6] + aht[161] - xv
                if float(vtin['type']) == STRAIGHT_TAPERED:
                    stb[9] = tail_arm_x - syna[7] * avt[62]
                    stb[10] = vtin[6] - (vtin[6] - vtin[1]) * syna[7] / vtin[4]
                elif syna[7] <= vtin[4] - vtin[2]:
                    stb[9] = tail_arm_x - syna[7] * avt[62]
                    stb[10] = (vtin[6] - (vtin[6] - vtin[5]) * syna[7] /
                               (vtin[4] - vtin[2]))
                else:
                    stb[9] = (tail_arm_x - (vtin[4] - vtin[2]) * avt[62] -
                              (syna[7] + vtin[2] - vtin[4]) * avt[86])
                    stb[10] = (vtin[1] + (vtin[5] - vtin[1]) *
                               (vtin[4] - syna[7]) / vtin[2])
    if win[13] == UNUSED:
        win[13] = 0.0
    if win[12] == UNUSED:
        win[12] = 0.0
    if win[14] == UNUSED:
        win[14] = win[13]
    stb[122] = (win[13] * (win[3] - win[12]) + win[14] * win[12]) / win[3]
    if tail_pass:
        # STBH(123) = STB(123): the wing pass's value, or whatever the
        # block held when the wing pass did not reach the body.
        stb[123] = stb123_wing
    return {'stb': stb, 'dihedral': {12: win[12], 13: win[13], 14: win[14]}}


def calculate_m17o21(alpha_count: int,
                     body: Dict[str, object],
                     wing: Optional[Dict[str, object]],
                     tail: Optional[Dict[str, object]],
                     transonic: bool = False) -> Dict[str, Dict[str, object]]:
    """Translate M17O21's combination pass.

    Args:
        alpha_count: ``NALPHA``.
        body: ``cyb``, ``cnb`` (``BODY(141)``, ``(161)``) and ``clb``
            (``BODY(181)`` onward).
        wing: The first SUBLAT result (wing, body, vertical tail), or None.
        tail: The second (horizontal tail, body, ventral fin), or None.
        transonic: ``TRANSN``.

    Returns:
        ``{'bv', 'bwh', 'bwv', 'bwhv'}``, each ``{'cyb', 'cnb', 'clb'}``
        with CY_beta and Cn_beta at the first angle only (the source marks
        the rest ``-UNUSED``) and Cl_beta at every angle.

    Notes:
        On a transonic pass the source overwrites these sets, and the
        wing's and tail's first values, with ``UNUSED``.
    """
    zero_block = {'cyb': 0.0, 'cnb': 0.0, 'clb': np.zeros(alpha_count)}

    def block_or_zero(result, key):
        """Missing SUBLAT output stays at the source's zero block."""
        if result is None or key not in result:
            return zero_block
        return dict(zero_block, **result[key])

    vertical_tail = block_or_zero(wing, 'vertical')
    ventral_fin = block_or_zero(tail, 'vertical')
    wing_body = block_or_zero(wing, 'combination')
    tail_surface = block_or_zero(tail, 'surface')

    def first_angle(value):
        return float(np.asarray(value, dtype=float).ravel()[0])

    def combine(*blocks):
        return {
            'cyb': sum(first_angle(b['cyb']) for b in blocks),
            'cnb': sum(first_angle(b['cnb']) for b in blocks),
            'clb': sum(np.asarray(b['clb'], dtype=float) for b in blocks),
        }

    body_block = {'cyb': body['cyb'], 'cnb': body['cnb'], 'clb': body['clb']}
    bwh = combine(wing_body, tail_surface)
    out = {
        'bv': combine(body_block, vertical_tail, ventral_fin),
        'bwh': bwh,
        'bwv': combine(wing_body, vertical_tail, ventral_fin),
        'bwhv': combine(bwh, vertical_tail, ventral_fin),
    }
    if transonic:
        # Labels 1060-1070: everything but BW-V's first CY_beta and
        # Cn_beta, which were summed before the overwrite.
        missing = np.full(alpha_count, UNUSED)
        for key in ('bv', 'bwh', 'bwhv'):
            out[key] = {'cyb': UNUSED, 'cnb': UNUSED, 'clb': missing.copy()}
        out['bwv']['clb'] = missing.copy()
    return out
