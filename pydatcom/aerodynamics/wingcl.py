"""
WINGCL: transonic wing lift and rolling moment due to sideslip.

Two of the source routine's three sections are translated here:

- **Transonic wing CL**, DATCOM Section 4.1.3.3: linear up to the stall
  onset angle, then a polynomial through the maximum-lift point.
- **Transonic wing CLB**, DATCOM Equation 5.1.2.1-C: rolling moment due to
  sideslip, interpolated on Mach between its subsonic and supersonic
  anchors.

The third section, transonic wing ``CDL`` from Section 4.1.5.2, is **not**
translated.  Its two figure tables are defective in the source: ``DEP55A``
and ``DEP55B`` are each declared and filled with 164 elements, but the
``INTERX`` call that reads them passes ``LENG = (7, 6, 4)``, which needs
7*6*4 = 168.  Both arrays are therefore read four elements past their end.
Counting the ``DATA`` statements shows why: each of the four taper groups
holds 41 values rather than 42, because its last aspect-ratio row carries
six entries where the other rows carry seven.  In FORTRAN's static storage
the overrun silently picks up whatever follows in memory, so the result at
the highest taper ratio depends on array layout rather than on the chart.

That cannot be reproduced meaningfully, and guessing the four missing values
would mean inventing chart data.  ``calculate_wingcl_cdl`` therefore raises
unless the caller supplies a completed 168-element table explicitly.

The damage in the original is wider than the last taper ratio.  The lookup
expects 42 values per taper group but the data packs 41, so the groups for
taper ratios 0.2, 0.5 and 1.0 are read shifted by one, two and three
positions, and the zero-taper group's last point takes the next group's
first value.  Against the published charts (datcom-legacy/figure
4.1.5.2-55A.png and -55B.png) the short row in every group is the
``A(t/c)**(1/3) = 2.0`` curve, missing one point at ``(M**2-1)/(t/c)**(2/3)``
of -3 or -2; TRANSLATION_STATUS.md lists the chart readings.

Every WINGCL entry point issues :class:`WingclTableWarning`, on every call.

Reference: datcom-legacy/datcom_2000/wingcl.f
"""

import numpy as np
from typing import Dict, Optional, Sequence
import logging
import warnings

from pydatcom.utils.constants import UNUSED
from pydatcom.utils.legacy_interp import interx
from pydatcom.utils.legacy_numeric import tbfunx

logger = logging.getLogger(__name__)


class WingclTableWarning(UserWarning):
    """WINGCL's transonic drag-due-to-lift tables are incomplete in the
    DATCOM source (see the module docstring)."""


# Shown on every call, not once per call site.
warnings.simplefilter('always', WingclTableWarning)

WINGCL_TABLE_WARNING = (
    "WINGCL: DATCOM's transonic drag-due-to-lift tables (Figures "
    "4.1.5.2-55A/B, DEP55A/DEP55B) hold 164 of the 168 values their lookup "
    "reads; the original program misreads them for any taper ratio above 0 "
    "and reads past the arrays at taper 1.0. The CDL section is not "
    "computed without a completed table. See TRANSLATION_STATUS.md, "
    "'A source defect left unresolved: WINGCL's CDL tables'.")


def _warn_incomplete_tables() -> None:
    warnings.warn(WINGCL_TABLE_WARNING, WingclTableWarning, stacklevel=3)
    logger.warning(WINGCL_TABLE_WARNING)

# Figure 4.1.5.2-55A/B independent grids.  LIND=7, so the three grids start
# at offsets 0, 7 and 14.
_PARM = [-4., -3., -2., -1., 0., 1., 2.,
         .5, .75, 1., 1.5, 1.75, 2.,
         0., 0., .2, .5, 1., 0., 0., 0.]
_CDL_SHAPE = (7, 6, 4)
_CDL_REQUIRED = _CDL_SHAPE[0] * _CDL_SHAPE[1] * _CDL_SHAPE[2]
_CDL_SOURCE_LENGTH = 164

# Figures 4.1.5.2-55A/B as the source's DATA statements hold them (164
# values each, 41 per taper group; see the module docstring), extracted by
# parsing with tools/fortran_data.py.  They are incomplete, so the CDL
# section still needs a completed 168-value table from the caller.
DEP55A = [
    1.2, 1.14, 1.13, 1.17, 1.2, 1.18, 1.13, 0.94, 0.88, 0.88,
    0.9, 0.95, 0.98, 0.97, 0.75, 0.7, 0.68, 0.71, 0.75, 0.78,
    0.8, 0.4, 0.4, 0.4, 0.46, 0.54, 0.55, 0.53, 0.27, 0.28,
    0.3, 0.37, 0.49, 0.48, 0.46, 0.16, 0.17, 0.27, 0.43, 0.42,
    0.4, 1.05, 1.04, 1.04, 1.05, 1.06, 1.05, 1.02, 0.8, 0.77,
    0.76, 0.78, 0.8, 0.82, 0.82, 0.6, 0.58, 0.58, 0.58, 0.59,
    0.62, 0.65, 0.36, 0.36, 0.37, 0.42, 0.47, 0.49, 0.5, 0.28,
    0.27, 0.3, 0.36, 0.42, 0.44, 0.45, 0.2, 0.24, 0.32, 0.39,
    0.41, 0.42, 1.02, 1.03, 1.02, 1.0, 0.97, 0.96, 0.97, 0.78,
    0.76, 0.73, 0.7, 0.68, 0.68, 0.7, 0.58, 0.55, 0.54, 0.52,
    0.52, 0.53, 0.56, 0.33, 0.34, 0.36, 0.39, 0.43, 0.45, 0.47,
    0.28, 0.29, 0.31, 0.34, 0.39, 0.43, 0.45, 0.23, 0.27, 0.32,
    0.38, 0.41, 0.43, 0.82, 0.84, 0.88, 0.92, 0.98, 1.02, 1.02,
    0.63, 0.63, 0.67, 0.75, 0.79, 0.82, 0.82, 0.51, 0.52, 0.56,
    0.62, 0.68, 0.7, 0.71, 0.4, 0.4, 0.41, 0.52, 0.55, 0.62,
    0.63, 0.37, 0.37, 0.37, 0.45, 0.54, 0.6, 0.62, 0.35, 0.35,
    0.42, 0.5, 0.56, 0.58,
]
DEP55B = [
    1.22, 1.16, 1.14, 1.18, 1.22, 1.19, 1.14, 1.03, 1.01, 0.98,
    0.96, 0.95, 0.96, 0.98, 0.87, 0.84, 0.81, 0.79, 0.77, 0.78,
    0.8, 0.57, 0.51, 0.49, 0.52, 0.55, 0.56, 0.54, 0.31, 0.31,
    0.31, 0.37, 0.49, 0.49, 0.45, 0.17, 0.21, 0.29, 0.43, 0.43,
    0.4, 1.11, 1.1, 1.09, 1.08, 1.07, 1.05, 1.03, 0.92, 0.88,
    0.86, 0.84, 0.84, 0.82, 0.83, 0.75, 0.73, 0.71, 0.68, 0.66,
    0.66, 0.67, 0.46, 0.42, 0.41, 0.45, 0.51, 0.52, 0.51, 0.28,
    0.29, 0.3, 0.32, 0.44, 0.46, 0.47, 0.2, 0.24, 0.32, 0.4,
    0.43, 0.44, 1.11, 1.09, 1.07, 1.03, 0.99, 0.98, 0.99, 0.9,
    0.88, 0.86, 0.83, 0.79, 0.78, 0.81, 0.7, 0.7, 0.68, 0.64,
    0.6, 0.62, 0.66, 0.46, 0.46, 0.47, 0.49, 0.5, 0.52, 0.52,
    0.4, 0.39, 0.4, 0.42, 0.45, 0.47, 0.48, 0.37, 0.35, 0.38,
    0.41, 0.44, 0.44, 0.83, 0.87, 0.9, 0.94, 0.99, 0.92, 0.93,
    0.7, 0.72, 0.74, 0.77, 0.8, 0.84, 0.87, 0.6, 0.61, 0.63,
    0.67, 0.73, 0.76, 0.78, 0.46, 0.46, 0.48, 0.52, 0.6, 0.64,
    0.66, 0.4, 0.39, 0.4, 0.47, 0.56, 0.61, 0.63, 0.37, 0.37,
    0.43, 0.53, 0.58, 0.6,
]
# Correction checked against the published chart (datcom-legacy/figure
# 4.1.5.2-55B.png): at taper 1.0 the A(t/c)**(1/3) = 0.5 curve keeps rising
# past 0.99 at (M**2-1)/(t/c)**(2/3) = 0, but the source has 0.92 and 0.93
# at +1 and +2, a transcription slip.  The chart reads 1.01 and 1.02.
DEP55B_CORRECTIONS = {128: 1.01, 129: 1.02}
for _index, _value in DEP55B_CORRECTIONS.items():
    DEP55B[_index] = _value

# The two figures are anchored at these values of AR*tan(leading-edge sweep).
_CDL_ANCHORS = np.array([0.0, 3.0])


def calculate_wingcl(alpha_schedule: Sequence[float],
                     alpha_zero_lift: float, alpha_stall_onset: float,
                     alpha_clmax: float, cla: float, clmax: float,
                     cl_supplied: Optional[Sequence[float]] = None
                     ) -> Dict[str, object]:
    """Translate WINGCL's transonic wing CL, DATCOM Section 4.1.3.3.

    Lift is linear at ``CLA`` up to the stall onset angle ``ALPHAS``, then
    follows a polynomial that reaches ``CLMAX`` at ``ACLMAX``.  Angles past
    ``ACLMAX`` are left unset: the source breaks out of its loop there.

    Args:
        alpha_schedule: Angle-of-attack schedule, degrees.
        alpha_zero_lift: ``ALPHAO``.
        alpha_stall_onset: ``ALPHAS``, where the curve leaves linear.
        alpha_clmax: ``ACLMAX``, the angle of maximum lift.
        cla: Linear lift-curve slope, per degree.
        clmax: Maximum lift coefficient.
        cl_supplied: Per-angle values already known; entries that are not
            None are left untouched, matching the source's UNUSED test.

    Returns:
        Dictionary with ``cl`` and the adjusted stall parameters.  Angles
        beyond ``ACLMAX`` come back as NaN.

    Raises:
        ValueError: If the schedule is empty or the stall angles are
            inconsistent.

    Notes:
        The source raises ``ACLMAX`` and ``CLMAX`` to the linear-extrapolated
        values whenever the supplied ones fall below the linear curve at the
        stall onset angle, so a configuration whose stated maximum lift is
        weaker than its own linear curve is corrected rather than rejected.
    """
    _warn_incomplete_tables()
    alpha = np.atleast_1d(np.asarray(alpha_schedule, dtype=float))
    if alpha.size == 0:
        raise ValueError("WINGCL needs a nonempty angle schedule")
    if alpha_stall_onset <= alpha_zero_lift:
        raise ValueError("stall onset must lie above the zero-lift angle")

    # The source's consistency correction.
    cl_reference = (alpha_stall_onset - alpha_zero_lift) * cla
    corrected = False
    if alpha_clmax < alpha_stall_onset or clmax < cl_reference:
        alpha_clmax = alpha_stall_onset
        clmax = cl_reference
        corrected = True

    nonlinear = alpha_stall_onset != alpha_clmax
    poly_a0 = clmax
    poly_a1 = 0.0
    poly_a2 = cla * (alpha_stall_onset - alpha_zero_lift) - clmax
    exponent = None
    if nonlinear:
        denominator = cla * (alpha_stall_onset - alpha_zero_lift) - clmax
        if denominator == 0.0:
            raise ValueError(
                "WINGCL's nonlinear exponent divides by "
                "CLA*(ALPHAS-ALPHAO) - CLMAX, which vanished")
        exponent = cla * (alpha_stall_onset - alpha_clmax) / denominator

    supplied = ([None] * alpha.size if cl_supplied is None
                else list(cl_supplied))
    if len(supplied) != alpha.size:
        raise ValueError("supplied CL must match the angle schedule")

    cl = np.full(alpha.size, np.nan)
    for angle_slot, angle_deg in enumerate(alpha):
        if angle_deg > alpha_clmax:
            # The source leaves this and every later angle untouched.
            break
        if supplied[angle_slot] is not None:
            cl[angle_slot] = supplied[angle_slot]
            continue
        if angle_deg <= alpha_stall_onset:
            cl[angle_slot] = cla * (angle_deg - alpha_zero_lift)
        elif nonlinear:
            stall_fraction = ((angle_deg - alpha_clmax) /
                              (alpha_stall_onset - alpha_clmax))
            cl[angle_slot] = (poly_a0 + poly_a1 * (angle_deg - alpha_clmax) +
                              poly_a2 * stall_fraction ** exponent)

    return {
        'cl': cl,
        'alpha_clmax': float(alpha_clmax),
        'clmax': float(clmax),
        'stall_corrected': corrected,
        'nonlinear': nonlinear,
        'exponent': exponent,
        'method': 'legacy_wingcl_lift',
    }


def calculate_wingcl_clb(cl: Sequence[float], mach: float, cla: float,
                         clb_subsonic: float, clb_supersonic: float,
                         cla_mach06: float, cla_mach14: float,
                         clb_supplied: Optional[Sequence[float]] = None
                         ) -> Dict[str, object]:
    """Translate WINGCL's transonic CLB, DATCOM Equation 5.1.2.1-C.

    ``CLB/CL`` is interpolated linearly on Mach between the subsonic anchor
    at 0.6 and the supersonic anchor at 1.4, each normalised by the lift
    slope at its own anchor Mach, then rescaled by the local slope squared.

    Args:
        cl: Lift coefficient at each angle.
        mach: Free-stream Mach number.
        cla: Local lift-curve slope.
        clb_subsonic: ``CLBSB``, the Mach 0.6 anchor.
        clb_supersonic: ``CLBSS``, the Mach 1.4 anchor.
        cla_mach06: ``CLA6``, the lift slope at Mach 0.6.
        cla_mach14: ``CLA14``, the lift slope at Mach 1.4.
        clb_supplied: Per-angle values already known.

    Returns:
        Dictionary with ``clb`` per angle and the ``clb_per_cl`` factor.

    Raises:
        ValueError: If either anchor lift slope is zero.
    """
    _warn_incomplete_tables()
    lift = np.atleast_1d(np.asarray(cl, dtype=float))
    if cla_mach06 == 0.0 or cla_mach14 == 0.0:
        raise ValueError("WINGCL CLB divides by the anchor lift slopes")

    clb_per_cl_mach06 = clb_subsonic / cla_mach06 ** 2
    clb_per_cl_mach14 = clb_supersonic / cla_mach14 ** 2
    anchor_separation = 0.8  # Mach 1.4 minus Mach 0.6
    clb_per_cl = ((clb_per_cl_mach14 - clb_per_cl_mach06) *
                  (mach - 0.6) / anchor_separation + clb_per_cl_mach06)
    clb_per_cl *= cla**2

    supplied = ([None] * lift.size if clb_supplied is None
                else list(clb_supplied))
    if len(supplied) != lift.size:
        raise ValueError("supplied CLB must match the CL schedule")

    clb = np.full(lift.size, np.nan)
    for angle_slot, cl_value in enumerate(lift):
        if not np.isfinite(cl_value) or abs(abs(cl_value) - UNUSED) <= 1.0e-40:
            continue
        clb[angle_slot] = (supplied[angle_slot]
                             if supplied[angle_slot] is not None
                             else clb_per_cl * cl_value)

    return {
        'clb': clb,
        'clb_per_cl': float(clb_per_cl),
        'method': 'legacy_wingcl_clb',
    }


def calculate_wingcl_cdl(mach: float, thickness_ratio: float,
                         aspect_ratio: float, taper_ratio: float,
                         sweep_le_rad: float,
                         table_55a: Sequence[float],
                         table_55b: Sequence[float]) -> Dict[str, float]:
    """Transonic wing CDL/CL^2, Section 4.1.5.2.

    The caller must supply both figure tables at their full
    ``7*6*4 = 168`` elements.  See the module docstring: the source's own
    ``DEP55A`` and ``DEP55B`` hold only 164, so the original reads four
    elements past the end of each and its result at the highest taper ratio
    depends on memory layout.  Completing the tables is a chart-digitisation
    question, not a translation one, so it is left to the caller.

    Args:
        mach: Free-stream Mach number.
        thickness_ratio: ``TOC``.
        aspect_ratio: ``AR``.
        taper_ratio: ``TAPR``.
        sweep_le_rad: Leading-edge sweep, radians, the source's ``OMEGA``.
        table_55a: Figure 4.1.5.2-55A at AR*tan(sweep) = 0, 168 elements.
        table_55b: Figure 4.1.5.2-55B at AR*tan(sweep) = 3, 168 elements.

    Returns:
        Dictionary with ``cdl_per_cl2`` and the two anchor values.

    Raises:
        ValueError: If either table is not exactly 168 elements, with the
            source's own length called out.
    """
    _warn_incomplete_tables()
    for name, table in (('55A', table_55a), ('55B', table_55b)):
        if len(table) != _CDL_REQUIRED:
            raise ValueError(
                f"Figure 4.1.5.2-{name} needs {_CDL_REQUIRED} elements for a "
                f"{_CDL_SHAPE[0]}x{_CDL_SHAPE[1]}x{_CDL_SHAPE[2]} table; got "
                f"{len(table)}. The source array holds only "
                f"{_CDL_SOURCE_LENGTH}, four short, and reads past its end.")
    if thickness_ratio <= 0.0:
        raise ValueError("CDL requires a positive thickness ratio")

    tc_one_third = thickness_ratio ** (1.0 / 3.0)
    tc_two_thirds = thickness_ratio ** (2.0 / 3.0)
    lookup = [
        (mach ** 2 - 1.0) / tc_two_thirds,
        aspect_ratio * tc_one_third,
        taper_ratio,
    ]

    anchor_at_sweep0 = interx(
        3, _PARM, lookup, list(_CDL_SHAPE), table_55a,
        lind=7, lx1l=1, lx2l=1, lx3l=1, lx1u=1, lx2u=1, lx3u=1,
    )
    anchor_at_sweep3 = interx(
        3, _PARM, lookup, list(_CDL_SHAPE), table_55b,
        lind=7, lx1l=1, lx2l=1, lx3l=1, lx1u=1, lx2u=1, lx3u=1,
    )
    anchor_values = np.array([anchor_at_sweep0, anchor_at_sweep3])
    blended = tbfunx(
        _CDL_ANCHORS, anchor_values,
        aspect_ratio * np.tan(sweep_le_rad),
        lower=1, upper=1,
    )[0]
    return {
        'cdl_per_cl2': float(blended * tc_one_third),
        'anchor_sweep0': float(anchor_at_sweep0),
        'anchor_sweep3': float(anchor_at_sweep3),
        'method': 'legacy_wingcl_cdl',
    }
