"""
CDRAG: zero-lift and lift-dependent drag of a lifting surface.

Total drag is the sum of a zero-lift part and a lift-dependent part.  The
zero-lift part is the same everywhere: a turbulent skin-friction coefficient
from Figure 4.1.5.1-26, capped by the surface-roughness cutoff Reynolds
number of Figure 4.1.5.1-27, times a thickness form factor, times the
lifting-surface correlation factor of Figure 4.1.5.1-28B, times the wetted
area ratio.

The lift-dependent part takes one of three paths:

===================  =====================================================
planform              method
===================  =====================================================
straight tapered      Oswald efficiency from the leading-edge suction of
                      Figure 4.1.5.2-53A/B, plus the twist terms of
                      Figures 4.1.5.2-42 and -48
cranked               span-weighted suction across the two panels, plus the
                      increment of Figure 4.1.5.2-54
double delta, curved  a single vortex-lift term in the local angle
===================  =====================================================

A non-straight planform computes its zero-lift drag as the sum of separate
inboard and outboard contributions, each with its own Reynolds number,
thickness and maximum-thickness sweep.

When a body is present the routine scales the supplied lift and lift-curve
slope by the wing-body carryover factor of Figure 4.3.1.2-10A, uses the
scaled values throughout, and restores the originals before returning.  This
translation takes them as inputs and never mutates the caller's arrays.

Reference: datcom-legacy/datcom_2000/cdrag.f
"""

import numpy as np
from typing import Dict, Optional, Sequence
import logging

from pydatcom.utils.constants import PI, RAD
from pydatcom.utils.legacy_numeric import tbfunx
from pydatcom.utils.legacy_tables import tlinex, tlin3x
from pydatcom.utils.table_lookup import fig26, fig53a

logger = logging.getLogger(__name__)

# The source's WTYPE entries.
STRAIGHT_TAPERED = 1.0
DOUBLE_DELTA = 2.0
CRANKED = 3.0
CURVED = 4.0

# The source replaces a zero leading-edge sweep tangent with this before
# dividing, and puts it back afterwards.
_ZERO_TANGENT = 0.00001

# Above this vortex Reynolds number the suction parameter comes from the
# Figure 4.1.5.2-53B curve instead of the 53A surface.
_SUCTION_REYNOLDS_SPLIT = 1.30e5

# The thickness form factor's leading coefficient switches here.
_THICKNESS_STATION_SPLIT = 0.30
_CAPL_FORWARD = 2.00
_CAPL_AFT = 1.20


# Figure 4.1.5.1-28B: X228B, cos of the maximum-thickness sweep.
_F28B_COS_SWEEP = np.array([
    0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0,
])

# X128B, a descending Mach axis.
_F28B_MACH = np.array([
    0.9, 0.8, 0.6, 0.25,
])

# Y28B, the lifting-surface correlation factor.
_F28B_RLS = np.array([
    1.1, 1.13, 1.17, 1.2, 1.24, 1.27, 1.3, 1.33, 1.34, 1.35, 1.36, 1.0,
    1.04, 1.08, 1.11, 1.15, 1.18, 1.21, 1.23, 1.25, 1.25, 1.26, 0.88,
    0.92, 0.96, 1.0, 1.04, 1.08, 1.11, 1.13, 1.14, 1.14, 1.15, 0.81,
    0.85, 0.89, 0.925, 0.96, 1.0, 1.03, 1.05, 1.06, 1.06, 1.07,
]).reshape((11, 4), order='F')

# Figure 4.1.5.1-28B dashed: X228BD.
_F28BD_COS_SWEEP = np.array([
    0.45, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0,
])

# X128BD.
_F28BD_MACH = np.array([
    0.9, 0.8, 0.6, 0.25,
])

# Y28BD.
_F28BD_RLS = np.array([
    1.2, 1.2, 1.24, 1.27, 1.3, 1.33, 1.34, 1.35, 1.36, 1.11, 1.11, 1.15,
    1.18, 1.21, 1.23, 1.25, 1.25, 1.26, 1.0, 1.0, 1.04, 1.08, 1.11,
    1.13, 1.14, 1.14, 1.15, 0.925, 0.925, 0.96, 1.0, 1.03, 1.05, 1.06,
    1.06, 1.07,
]).reshape((9, 4), order='F')

# Figure 4.1.5.2-42: X142, atan(tan(sweep c/4)/beta) in degrees.
_F42_ANGLE = np.array([
    60.0, 45.0, 30.0, 0.0,
])

# X242, the lift-dependent drag factor.
_F42_DRAG_FACTOR = np.array([
    2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 6.0, 8.0, 10.0, 12.0, 14.0,
])

# X342, taper ratio.
_F42_TAPER = np.array([
    0.1, 0.2, 0.25, 0.3, 0.4, 0.5, 0.6, 0.75, 1.0,
])

# Y42, assembled from the four EQUIVALENCEd quarters.
_F42 = np.array([
    0.00025, 0.0002, 0.0001, 0.0, -0.0002, -0.00028, -0.00038, -0.00055,
    -0.0008, -0.00093, -0.00097, -0.0009, -0.0004, -0.00042, -0.00046,
    -0.00052, -0.00062, -0.00069, -0.00075, -0.00088, -0.00111,
    -0.00127, -0.00127, -0.00125, -0.0009, -0.001, -0.0011, -0.00118,
    -0.00127, -0.00136, -0.00144, -0.00157, -0.00177, -0.0018, -0.00188,
    -0.00181, -0.00193, -0.00225, -0.00253, -0.00273, -0.00292,
    -0.00308, -0.00318, -0.00336, -0.0035, -0.0035, -0.00345, -0.0034,
    0.0012, 0.00135, 0.00144, 0.00149, 0.00148, 0.00145, 0.00136,
    0.00113, 0.00077, 0.00052, 0.00032, 0.000175, 0.0002, 0.0004,
    0.00053, 0.00062, 0.00067, 0.0007, 0.00068, 0.00063, 0.00042,
    0.00023, 0.00014, 0.0001, -0.00039, -0.00032, -0.00028, -0.00024,
    -0.0002, -0.00019, -0.00018, -0.0002, -0.00028, -0.00043, -0.00052,
    -0.00053, -0.0013, -0.0015, -0.00166, -0.0018, -0.00192, -0.00203,
    -0.00211, -0.0022, -0.0023, -0.00229, -0.0022, -0.00214, 0.00155,
    0.0018, 0.00197, 0.00207, 0.00209, 0.00203, 0.00192, 0.00173,
    0.00142, 0.0011, 0.00085, 0.00072, 0.00055, 0.00075, 0.00092,
    0.00106, 0.00115, 0.00123, 0.00127, 0.00125, 0.00105, 0.00085,
    0.0007, 0.00065, -0.00018, -8e-05, 0.0, 0.0001, 0.00017, 0.00023,
    0.00027, 0.00032, 0.0003, 0.00025, 0.00018, 5e-05, -0.00106,
    -0.00122, -0.00133, -0.00145, -0.00154, -0.00162, -0.00166,
    -0.00175, -0.00178, -0.00174, -0.00172, -0.0017, 0.00188, 0.00225,
    0.00248, 0.00262, 0.00267, 0.00263, 0.00257, 0.00243, 0.00212,
    0.00176, 0.00151, 0.0014, 0.00065, 0.001, 0.0013, 0.00152, 0.00165,
    0.00173, 0.00175, 0.00178, 0.0017, 0.0015, 0.00132, 0.00123, 0.0,
    0.00018, 0.00033, 0.00048, 0.00058, 0.00069, 0.00075, 0.00083,
    0.00085, 0.0008, 0.00069, 0.00064, -0.00088, -0.001, -0.00108,
    -0.00114, -0.00119, -0.00124, -0.00126, -0.0013, -0.0013, -0.00125,
    -0.00117, -0.00108, 0.00242, 0.003, 0.00342, 0.00365, 0.00369,
    0.00365, 0.00362, 0.0035, 0.00318, 0.0028, 0.00245, 0.00221,
    0.00107, 0.0015, 0.0019, 0.00218, 0.00239, 0.00255, 0.00266,
    0.00275, 0.00272, 0.0026, 0.0024, 0.00214, 0.00033, 0.00072, 0.001,
    0.00117, 0.0013, 0.00142, 0.00153, 0.00167, 0.00178, 0.00175,
    0.00165, 0.00163, -0.00057, -0.00062, -0.00063, -0.00064, -0.00064,
    -0.00064, -0.00061, -0.00055, -0.00042, -0.00031, -0.00025,
    -0.00026, 0.0029, 0.00395, 0.00427, 0.00443, 0.0045, 0.00452,
    0.0045, 0.00443, 0.00413, 0.0037, 0.00334, 0.003, 0.00135, 0.00207,
    0.0025, 0.0028, 0.00305, 0.00325, 0.0034, 0.0036, 0.00365, 0.00345,
    0.00325, 0.003, 0.0006, 0.00096, 0.0013, 0.00156, 0.0018, 0.00198,
    0.00216, 0.00238, 0.00262, 0.00265, 0.0026, 0.0025, -0.0003,
    -0.00027, -0.00022, -0.00018, -0.000125, -7e-05, 0.0, 0.00013,
    0.0003, 0.00046, 0.00057, 0.00063, 0.0033, 0.00395, 0.0045, 0.00495,
    0.0052, 0.0053, 0.00534, 0.00527, 0.00488, 0.00445, 0.00402,
    0.00378, 0.0017, 0.0024, 0.0029, 0.0033, 0.0036, 0.00385, 0.004,
    0.00428, 0.00434, 0.0042, 0.00402, 0.00375, 0.00085, 0.00127,
    0.00165, 0.002, 0.00225, 0.0025, 0.00267, 0.00298, 0.00328, 0.00335,
    0.00333, 0.00323, -0.0001, 0.0, 0.00012, 0.0002, 0.00032, 0.00046,
    0.0005, 0.00065, 0.00092, 0.0011, 0.00124, 0.0013, 0.004, 0.0053,
    0.00576, 0.006, 0.00613, 0.0062, 0.00625, 0.00626, 0.00595, 0.0054,
    0.00495, 0.0046, 0.002, 0.00285, 0.00345, 0.00393, 0.00433, 0.00465,
    0.00487, 0.00513, 0.00532, 0.00517, 0.00495, 0.00468, 0.00115,
    0.00172, 0.00218, 0.00257, 0.0029, 0.00318, 0.0034, 0.00373, 0.0041,
    0.00423, 0.00425, 0.0042, 0.0002, 0.00035, 0.0005, 0.00067, 0.00084,
    0.001, 0.00113, 0.00137, 0.00175, 0.002, 0.00212, 0.00216, 0.0047,
    0.006, 0.0068, 0.00712, 0.0073, 0.00742, 0.00748, 0.00746, 0.00712,
    0.00658, 0.00615, 0.00575, 0.0026, 0.00365, 0.00442, 0.00495,
    0.00535, 0.00567, 0.00582, 0.00625, 0.00643, 0.00636, 0.00618,
    0.00592, 0.00155, 0.00225, 0.00282, 0.0033, 0.00372, 0.00405,
    0.00435, 0.00478, 0.0053, 0.00545, 0.00543, 0.00535, 0.00055,
    0.00085, 0.00112, 0.00138, 0.0016, 0.00182, 0.00202, 0.00237,
    0.00288, 0.00317, 0.00335, 0.00345,
]).reshape((12, 4, 9), order='F')

# Figure 4.1.5.2-48: X148.
_F48_ANGLE = np.array([
    0.0, 30.0, 45.0, 60.0,
])

# X248.
_F48_DRAG_FACTOR = np.array([
    2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 10.0, 12.0, 14.0,
])

# X348.
_F48_TAPER = np.array([
    0.1, 0.2, 0.25, 0.3, 0.4, 0.5, 0.6, 0.75, 1.0,
])

# Y48, assembled from the three EQUIVALENCEd thirds.
_F48 = np.array([
    0.00084, 0.00113, 0.00136, 0.00152, 0.00162, 0.00169, 0.00175,
    0.00178, 0.00177, 0.00172, 0.00084, 0.00113, 0.00133, 0.00146,
    0.00154, 0.00158, 0.00159, 0.00158, 0.00154, 0.00149, 0.00083,
    0.00109, 0.00124, 0.00133, 0.00136, 0.00137, 0.00134, 0.00127,
    0.00121, 0.00116, 0.0008, 0.00095, 0.00101, 0.00102, 0.001, 0.00098,
    0.00093, 0.00083, 0.00075, 0.0007, 0.00086, 0.00118, 0.00145,
    0.00163, 0.00176, 0.00185, 0.00192, 0.00198, 0.002, 0.00198,
    0.00086, 0.00118, 0.00145, 0.0016, 0.0017, 0.00176, 0.0018, 0.00179,
    0.00173, 0.00169, 0.00086, 0.00118, 0.00136, 0.00147, 0.00152,
    0.00154, 0.00153, 0.00147, 0.0014, 0.00135, 0.00084, 0.00105,
    0.00113, 0.00114, 0.00112, 0.00109, 0.00104, 0.00095, 0.00086,
    0.0008, 0.00087, 0.00122, 0.00147, 0.00166, 0.00181, 0.00192,
    0.00199, 0.00206, 0.00207, 0.00203, 0.00087, 0.00122, 0.00145,
    0.00163, 0.00174, 0.0018, 0.00184, 0.00187, 0.00184, 0.00177,
    0.00087, 0.00119, 0.00137, 0.00148, 0.00155, 0.00157, 0.00157,
    0.00153, 0.001445, 0.0014, 0.00085, 0.00105, 0.00114, 0.00117,
    0.00116, 0.00113, 0.00109, 0.001, 0.00091, 0.00085, 0.00088,
    0.00121, 0.00149, 0.0017, 0.00185, 0.00195, 0.00203, 0.00211,
    0.00212, 0.00207, 0.00088, 0.00121, 0.00149, 0.00166, 0.00177,
    0.00185, 0.00189, 0.00191, 0.00188, 0.00181, 0.00088, 0.00121,
    0.00141, 0.00153, 0.00159, 0.00161, 0.00162, 0.00159, 0.00152,
    0.00146, 0.00085, 0.00107, 0.00117, 0.0012, 0.001195, 0.00116,
    0.00113, 0.00104, 0.00097, 0.00091, 0.00088, 0.00124, 0.00153,
    0.00173, 0.00189, 0.002, 0.00209, 0.00218, 0.0022, 0.00217, 0.00088,
    0.00124, 0.00149, 0.00169, 0.00182, 0.0019, 0.00195, 0.00199,
    0.00195, 0.0019, 0.00088, 0.00123, 0.00143, 0.00156, 0.00163,
    0.00167, 0.00168, 0.00165, 0.00158, 0.00149, 0.00088, 0.0011,
    0.0012, 0.00124, 0.001235, 0.00122, 0.00119, 0.00109, 0.001,
    0.00093, 0.00089, 0.00125, 0.00152, 0.00175, 0.00191, 0.00203,
    0.00211, 0.00221, 0.00225, 0.00223, 0.00089, 0.00125, 0.00152,
    0.0017, 0.00184, 0.00192, 0.00198, 0.00202, 0.002, 0.00195, 0.00089,
    0.00124, 0.00143, 0.00158, 0.00166, 0.00169, 0.00171, 0.001685,
    0.00162, 0.00154, 0.00089, 0.00113, 0.00123, 0.00126, 0.00127,
    0.00125, 0.00121, 0.00114, 0.00105, 0.00097, 0.0009, 0.00125,
    0.00154, 0.00175, 0.00192, 0.00205, 0.00213, 0.00223, 0.00226,
    0.00225, 0.0009, 0.00125, 0.00151, 0.0017, 0.00183, 0.00192,
    0.00199, 0.00203, 0.00202, 0.00197, 0.0009, 0.00125, 0.00146,
    0.00159, 0.00166, 0.0017, 0.00171, 0.00169, 0.00164, 0.00156,
    0.0009, 0.0011, 0.00123, 0.00129, 0.00128, 0.00125, 0.00121,
    0.00114, 0.00106, 0.00102, 0.00087, 0.00124, 0.00153, 0.00176,
    0.00193, 0.00205, 0.00213, 0.00223, 0.00226, 0.00224, 0.00087,
    0.00124, 0.00151, 0.0017, 0.00182, 0.00191, 0.00197, 0.002025,
    0.00203, 0.00199, 0.00087, 0.00124, 0.00143, 0.00157, 0.00167,
    0.00172, 0.00173, 0.00171, 0.00165, 0.00157, 0.00087, 0.00111,
    0.00122, 0.00126, 0.00127, 0.00126, 0.00124, 0.00115, 0.00107,
    0.00102, 0.00087, 0.00123, 0.00151, 0.00173, 0.0019, 0.00203,
    0.00212, 0.00221, 0.00225, 0.00224, 0.00087, 0.00123, 0.0015,
    0.00167, 0.0018, 0.00191, 0.00196, 0.002, 0.002, 0.00198, 0.00087,
    0.00123, 0.00142, 0.00155, 0.00164, 0.00168, 0.00169, 0.00169,
    0.00165, 0.00159, 0.00087, 0.00111, 0.0012, 0.00125, 0.00126,
    0.00125, 0.00121, 0.00115, 0.00109, 0.00103,
]).reshape((10, 4, 9), order='F')

# Figure 4.1.5.2-53B, the high-Reynolds leading-edge suction curve.
_F53B_X = np.array([
    0.0, 1.0, 2.0, 3.0, 4.0, 6.0, 8.0, 10.0,
])

# Y53B.
_F53B_Y = np.array([
    0.72956, 0.91078, 0.93209, 0.9453, 0.9585, 0.96226, 0.96601, 0.96977,
])

# Figure 4.1.5.2-54: X254, aspect ratio.
_F54_ASPECT_RATIO = np.array([
    3.0, 3.5, 4.0, 4.5, 5.25, 5.5, 5.75, 6.0, 6.5, 7.0, 7.5, 8.0, 8.5,
    9.0,
])

# X154, a descending |CL|/AR axis.
_F54_CL_OVER_AR = np.array([
    0.24, 0.22, 0.2, 0.18, 0.16, 0.14, 0.12, 0.11, 0.1, 0.09, 0.08, 0.0,
])

# Y54, the cranked-wing lift-dependent drag increment.
_F54 = np.array([
    0.07, 0.18, 0.18, 0.18, 0.18, 0.18, 0.18, 0.18, 0.18, 0.18, 0.18,
    0.18, 0.18, 0.18, 0.047, 0.12, 0.177, 0.177, 0.177, 0.177, 0.177,
    0.177, 0.177, 0.177, 0.177, 0.177, 0.177, 0.177, 0.03, 0.08, 0.124,
    0.162, 0.197, 0.197, 0.197, 0.197, 0.197, 0.197, 0.197, 0.197,
    0.197, 0.197, 0.017, 0.06, 0.091, 0.116, 0.147, 0.147, 0.147, 0.147,
    0.147, 0.147, 0.147, 0.147, 0.147, 0.147, 0.0075, 0.0425, 0.063,
    0.079, 0.0975, 0.0975, 0.0975, 0.0975, 0.0975, 0.0975, 0.0975,
    0.0975, 0.0975, 0.0975, 0.00225, 0.023, 0.039, 0.05, 0.0625, 0.0625,
    0.0625, 0.0625, 0.0625, 0.0625, 0.0625, 0.0625, 0.0625, 0.0625, 0.0,
    0.01, 0.019, 0.028, 0.0385, 0.0375, 0.039, 0.044, 0.08, 0.133,
    0.133, 0.133, 0.133, 0.133, 0.0, 0.005, 0.0125, 0.019, 0.0265,
    0.024, 0.022, 0.0215, 0.0385, 0.068, 0.11, 0.16, 0.16, 0.16, 0.0,
    0.0025, 0.008, 0.0125, 0.0175, 0.0135, 0.009, 0.0075, 0.012, 0.025,
    0.048, 0.0825, 0.13, 0.13, 0.0, 0.0, 0.0035, 0.0075, 0.01, 0.0065,
    0.003, 0.001, 0.0, 0.005, 0.012, 0.026, 0.05, 0.096, 0.0, 0.0,
    0.0015, 0.003, 0.0045, 0.001, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.01,
    0.0295, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0,
]).reshape((14, 12), order='F')

# Figure 4.1.5.1-27: the cutoff-Reynolds intercept versus Mach.
_F27_MACH = np.array([
    0.0, 1.0, 2.0, 3.0,
])

# X27I, the intercept itself.
_F27_INTERCEPT = np.array([
    1.5778, 1.67221, 1.98509, 2.28874,
])

# Figure 4.3.1.2-10A: KW(B) versus body diameter over span.
_F4312_10A_DIAMETER_RATIO = np.array([
    0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0,
])

# Y1OA.
_F4312_10A_KWB = np.array([
    1.0, 1.08, 1.16, 1.26, 1.36, 1.46, 1.56, 1.67, 1.78, 1.89, 2.0,
])


def _cutoff_reynolds(mach: float, mac: float, roughness: float) -> float:
    """The Figure 4.1.5.1-27 surface-roughness cutoff Reynolds number.

    The source forms ``D(2)=12.0*A(16)/ROUGFC`` once, from the *exposed*
    MAC, and reuses it for every panel, so an inboard or outboard panel is
    capped by a cutoff built on a length it does not have.
    """
    intercept, _ = tbfunx(_F27_MACH, _F27_INTERCEPT, mach, 0, 0)
    base = (12.0 * mac / roughness) ** 1.0482
    return float(base * 10.0 ** intercept)


def _zero_lift_drag(mach: float, reynolds: float, cutoff: float,
                    cos_max_thickness_sweep: float, thickness_ratio: float,
                    thickness_station: float, area: float, sref: float,
                    dashed: bool = False) -> Dict[str, float]:
    """One panel's zero-lift drag: friction, form factor and correlation.

    ``dashed`` selects the dashed variant of Figure 4.1.5.1-28B, which the
    source uses for the inboard panel of a non-straight planform.
    """
    reynolds_used = min(reynolds, cutoff)
    friction = float(fig26(reynolds_used, mach))

    if dashed:
        # The source passes LX1L=-1 here.  A negative mode differs from 0
        # only in whether a diagnostic is printed, which is omitted.
        rls = float(tlinex(
            _F28BD_MACH, _F28BD_COS_SWEEP, _F28BD_RLS,
            mach, cos_max_thickness_sweep, -1, 0, 0, 0,
        ))
    else:
        rls = float(tlinex(
            _F28B_MACH, _F28B_COS_SWEEP, _F28B_RLS,
            mach, cos_max_thickness_sweep, 0, 2, -1, 2,
        ))

    if thickness_station >= _THICKNESS_STATION_SPLIT:
        capl = _CAPL_AFT
    else:
        capl = _CAPL_FORWARD
    form = 1.0 + capl * thickness_ratio + 100.0 * thickness_ratio**4
    wetted_ratio = 2.0 * area / sref
    cdo = friction * form * rls * wetted_ratio

    return {
        'cdo': float(cdo),
        'friction': friction,
        'rls': rls,
        'form_factor': float(form),
        'capl': float(capl),
        'reynolds_used': float(reynolds_used),
        'reynolds_was_capped': bool(cutoff < reynolds),
    }


def _suction_parameter(vortex_reynolds: float, argument: float) -> float:
    """Leading-edge suction from Figure 4.1.5.2-53A or its 53B curve."""
    if vortex_reynolds > _SUCTION_REYNOLDS_SPLIT:
        value, _ = tbfunx(_F53B_X, _F53B_Y, argument, 0, 0)
        return float(value)
    return float(fig53a(vortex_reynolds, argument))


def calculate_cdrag(mach: float,
                    alpha_deg: Sequence[float],
                    cl: Sequence[float],
                    cla_first: float,
                    geometry: Dict[str, float],
                    section: Dict[str, float],
                    sweep: Dict[str, float],
                    sref: float,
                    roughness: float,
                    reynolds_per_length: float,
                    beta: float,
                    planform_type: float = STRAIGHT_TAPERED,
                    body_present: bool = False) -> Dict[str, object]:
    """Translate CDRAG: lifting-surface drag over an angle schedule.

    Args:
        mach: Free-stream Mach number, the source's ``FLC(I+2)``.
        alpha_deg: Angle schedule, the source's ``B(23)`` onward.  Used only
            by the double-delta and curved path.
        cl: Lift coefficient at each angle, ``AINN(J+20)``.
        cla_first: Lift-curve slope at the first angle, ``AINN(101)``.
        geometry: Planform quantities ``{'area', 'area_inboard',
            'area_outboard', 'mac', 'mac_inboard', 'mac_outboard',
            'aspect_ratio', 'aspect_ratio_outboard', 'taper_ratio',
            'taper_ratio_outboard', 'span_inboard', 'sspne', 'sspn'}``.
            These are ``A(3)``, ``A(1)``, ``A(2)``, ``A(16)``, ``A(15)``,
            ``A(17)``, ``A(7)``, ``A(6)``, ``A(27)``, ``A(28)``, ``A(23)``,
            ``AIN(3)`` and ``AIN(4)``.
        section: ``{'tovc', 'xovc', 'tovco', 'xovco', 'leri', 'lero',
            'twista'}``, the source's ``AIN(16)``, ``AIN(18)``, ``AIN(65)``,
            ``AIN(66)``, ``AIN(62)``, ``AIN(63)`` and ``AIN(11)``.
        sweep: ``{'cos_le', 'tan_le', 'cos_le_inboard', 'tan_le_inboard',
            'cos_le_outboard', 'tan_le_outboard', 'tan_c4',
            'cos_max_thickness', 'cos_max_thickness_inboard',
            'cos_max_thickness_outboard'}``, the source's ``A(37)``,
            ``A(38)``, ``A(61)``, ``A(62)``, ``A(85)``, ``A(86)``,
            ``A(44)``, ``A(178)``, ``A(190)`` and ``A(184)``.
        sref: Reference area.
        roughness: Surface roughness, the source's ``ROUGFC``.
        reynolds_per_length: The source's ``A(129)``.
        beta: Compressibility parameter, the source's ``B(2)``.
        planform_type: ``AIN(15)``.
        body_present: The source's ``BO``.  Scales the supplied lift and
            lift-curve slope by the Figure 4.3.1.2-10A carryover factor.

    Returns:
        Dictionary with ``cd`` and ``cdl`` over the schedule, the scalar
        ``cdo``, the per-panel zero-lift breakdown, and the flags described
        below.

    Raises:
        ValueError: If the reference area, roughness or exposed area is
            nonpositive, or if the arrays mismatch.

    Notes:
        Four things the source does that a reading can miss:

        ``TEMPI`` is set to the literal ``0.0`` on the cranked path where
        the outboard branch computes ``A(6)*A(28)/A(85)``.  The inboard
        analogue is never formed, so the inboard suction parameter is
        always read at a zero abscissa.  Preserved and flagged by
        ``inboard_suction_argument_is_zero``.

        The Oswald efficiency carries a factor of 1.1 on the straight
        tapered path -- ``D(30)=1.1*TEMP/(...)`` -- that the cranked path's
        otherwise identical expression does not.  Preserved, and reported
        through ``oswald_carries_source_1_1_factor``.

        The cutoff Reynolds number is built once from the exposed MAC and
        reused for the inboard and outboard panels, which have their own
        MACs.  Preserved; ``reynolds_was_capped`` reports per panel.

        ``TYPE=A(161)`` is assigned and never read, like MAXCL's ``DEL4``.
        It is not represented here.

        The comments on the inboard and outboard form factors read "IF XTI
        GREATER TAN 0.30, THEN CAPL= 2.0", which is the opposite of the
        code beneath them.  All three branches agree with each other, so
        the comments are wrong rather than the code.
    """
    alpha = np.atleast_1d(np.asarray(alpha_deg, dtype=float))
    lift = np.atleast_1d(np.asarray(cl, dtype=float))
    if lift.shape != alpha.shape:
        raise ValueError("CDRAG needs a lift array matching the schedule")

    area = float(geometry['area'])
    if min(sref, roughness, area) <= 0.0:
        raise ValueError(
            "CDRAG requires a positive reference area, roughness and "
            "exposed area",
        )

    # Figure 4.3.1.2-10A: the wing-body carryover factor.
    carryover = 1.0
    if body_present:
        diameter_ratio = (
            (float(geometry['sspn']) - float(geometry['sspne']))
            / float(geometry['sspn'])
        )
        value, _ = tbfunx(
            _F4312_10A_DIAMETER_RATIO, _F4312_10A_KWB, diameter_ratio, 0, 0,
        )
        carryover = float(value)

    lift = lift * carryover
    cla = float(cla_first) * carryover

    # Zero-tangent guards on leading-edge sweep (source replaces 0 with 1e-5).
    tan_le = float(sweep['tan_le']) or _ZERO_TANGENT
    tan_le_inboard = float(sweep['tan_le_inboard']) or _ZERO_TANGENT
    tan_le_outboard = float(sweep['tan_le_outboard']) or _ZERO_TANGENT

    mac = float(geometry['mac'])
    cutoff = _cutoff_reynolds(mach, mac, roughness)
    kind = float(planform_type)

    if kind == STRAIGHT_TAPERED:
        panel = _zero_lift_drag(
            mach, mac * reynolds_per_length, cutoff,
            float(sweep['cos_max_thickness']), float(section['tovc']),
            float(section['xovc']), area, sref)
        cdo = panel['cdo']
        panels = {'whole': panel}
        result = _straight_lift_drag(
            mach, lift, cla, geometry, section, sweep, sref, beta,
            reynolds_per_length, tan_le, area)
    else:
        inboard = _zero_lift_drag(
            mach, float(geometry['mac_inboard']) * reynolds_per_length,
            cutoff, float(sweep['cos_max_thickness_inboard']),
            float(section['tovc']), float(section['xovc']),
            float(geometry['area_inboard']), sref, dashed=True)
        outboard = _zero_lift_drag(
            mach, float(geometry['mac_outboard']) * reynolds_per_length,
            cutoff, float(sweep['cos_max_thickness_outboard']),
            float(section['tovco']), float(section['xovco']),
            float(geometry['area_outboard']), sref)
        cdo = inboard['cdo'] + outboard['cdo']
        panels = {'inboard': inboard, 'outboard': outboard}
        if kind == CRANKED:
            result = _cranked_lift_drag(
                mach, lift, cla, geometry, section, sweep, sref,
                reynolds_per_length, tan_le_inboard, tan_le_outboard, area)
        else:
            # Double delta and curved: a single vortex-lift term.  D(3) is
            # set to 1.0 at the top of the routine and never changed.
            result = {'cdl': 0.95 * lift * np.tan(alpha / RAD) * 1.0,
                      'path': 'double_delta_or_curved'}

    cdl = result['cdl']
    return {
        'cd': cdo + cdl,
        'cdl': cdl,
        'cdo': float(cdo),
        'panels': panels,
        'carryover': float(carryover),
        'cutoff_reynolds': float(cutoff),
        'path': result['path'],
        'oswald_efficiency': result.get('oswald_efficiency'),
        'suction': result.get('suction'),
        'detail': result,
        'inboard_suction_argument_is_zero': result.get(
            'inboard_suction_argument_is_zero', False,
        ),
        'oswald_carries_source_1_1_factor': kind == STRAIGHT_TAPERED,
        'method': 'legacy_cdrag',
    }


def _straight_lift_drag(mach, lift, cla, geometry, section, sweep, sref,
                        beta, reynolds_per_length, tan_le, area):
    """Labels 1010-1040: the straight tapered lift-dependent drag."""
    aspect_ratio = float(geometry['aspect_ratio'])
    taper = float(geometry['taper_ratio'])
    cos_le = float(sweep['cos_le'])

    # Vortex Reynolds number uses leading-edge radius times exposed MAC.
    rler = reynolds_per_length * float(section['leri']) * float(geometry['mac'])
    vortex_reynolds = (
        rler / abs(tan_le) * np.sqrt(1.0 - mach**2 * cos_le**2)
    )
    suction_argument = aspect_ratio * taper / cos_le
    suction = _suction_parameter(float(vortex_reynolds), float(suction_argument))

    lift_slope_term = (cla * sref / area) * RAD / aspect_ratio
    # The 1.1 appears only on this path; see calculate_cdrag's notes.
    oswald = (
        1.1 * lift_slope_term
        / (suction * lift_slope_term + (1.0 - suction) * PI)
    )

    drag_factor = aspect_ratio * beta
    angle = float(np.arctan(float(sweep['tan_c4']) / beta) * RAD)
    v42 = float(tlin3x(
        _F42_ANGLE, _F42_DRAG_FACTOR, _F42_TAPER, _F42,
        angle, drag_factor, taper, 0, 2, 0, 2, 2, 0,
    ))
    v48 = float(tlin3x(
        _F48_ANGLE, _F48_DRAG_FACTOR, _F48_TAPER, _F48,
        angle, drag_factor, taper, 0, 2, 0, 2, 2, 0,
    ))

    fig48_over_beta = v48 / beta
    twist = PI * float(section['twista']) / RAD
    twist_doubled = 2.0 * twist
    sref_over_area = sref / area

    cdl = (
        lift**2 / (PI * aspect_ratio * oswald) * sref_over_area
        + twist_doubled * lift * v42
        + twist_doubled**2 * fig48_over_beta * area / sref
    )

    return {
        'cdl': cdl,
        'path': 'straight_tapered',
        'oswald_efficiency': float(oswald),
        'suction': float(suction),
        'vortex_reynolds': float(vortex_reynolds),
        'fig42': v42,
        'fig48': v48,
        'angle': angle,
    }


def _cranked_lift_drag(mach, lift, cla, geometry, section, sweep, sref,
                       reynolds_per_length, tan_le_inboard, tan_le_outboard,
                       area):
    """Labels 1050-1100: the cranked-wing lift-dependent drag."""
    aspect_ratio = float(geometry['aspect_ratio'])
    cos_in = float(sweep['cos_le_inboard'])
    cos_out = float(sweep['cos_le_outboard'])
    sref_over_area = sref / area

    rler_inboard = (
        reynolds_per_length * float(section['leri'])
        * float(geometry['mac_inboard'])
    )
    rler_outboard = (
        reynolds_per_length * float(section['lero'])
        * float(geometry['mac_outboard'])
    )

    reynolds_in = (
        rler_inboard / abs(tan_le_inboard)
        * np.sqrt(1.0 - (mach * cos_in) ** 2)
    )
    # The source sets TEMPI to a literal 0.0 here.  See calculate_cdrag.
    suction_in = _suction_parameter(float(reynolds_in), 0.0)

    reynolds_out = (
        rler_outboard / abs(tan_le_outboard)
        * np.sqrt(1.0 - (mach * cos_out) ** 2)
    )
    argument_out = (
        float(geometry['aspect_ratio_outboard'])
        * float(geometry['taper_ratio_outboard'])
        / cos_out
    )
    suction_out = _suction_parameter(float(reynolds_out), float(argument_out))

    span_fraction = float(geometry['span_inboard']) / float(geometry['sspne'])
    suction = (
        suction_in * span_fraction + suction_out * (1.0 - span_fraction)
    )

    lift_slope_term = (cla * sref / area) * RAD / aspect_ratio
    oswald = lift_slope_term / (
        suction * lift_slope_term + (1.0 - suction) * PI
    )

    increments = np.array([
        float(tlinex(
            _F54_CL_OVER_AR, _F54_ASPECT_RATIO, _F54,
            float(abs(value) / aspect_ratio * sref_over_area),
            aspect_ratio, 0, 0, 0, 0,
        ))
        for value in lift
    ])
    cdl = (
        lift**2 / (PI * oswald * aspect_ratio) * sref_over_area
        + increments * area / sref
    )

    return {
        'cdl': cdl,
        'path': 'cranked',
        'oswald_efficiency': float(oswald),
        'suction': float(suction),
        'suction_inboard': float(suction_in),
        'suction_outboard': float(suction_out),
        'fig54': increments,
        'inboard_suction_argument_is_zero': True,
    }
