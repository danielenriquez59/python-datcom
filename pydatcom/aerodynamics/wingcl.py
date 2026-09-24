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

Reference: datcom-legacy/datcom_2000/wingcl.f
"""

import numpy as np
from typing import Dict, Optional, Sequence
import logging

from pydatcom.utils.constants import UNUSED
from pydatcom.utils.legacy_interp import interx
from pydatcom.utils.legacy_numeric import tbfunx

logger = logging.getLogger(__name__)

# Figure 4.1.5.2-55A/B independent grids.  LIND=7, so the three grids start
# at offsets 0, 7 and 14.
_PARM = [-4., -3., -2., -1., 0., 1., 2.,
         .5, .75, 1., 1.5, 1.75, 2.,
         0., 0., .2, .5, 1., 0., 0., 0.]
_CDL_SHAPE = (7, 6, 4)
_CDL_REQUIRED = _CDL_SHAPE[0] * _CDL_SHAPE[1] * _CDL_SHAPE[2]
_CDL_SOURCE_LENGTH = 164

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
    for index, angle in enumerate(alpha):
        if angle > alpha_clmax:
            # The source leaves this and every later angle untouched.
            break
        if supplied[index] is not None:
            cl[index] = supplied[index]
            continue
        if angle <= alpha_stall_onset:
            cl[index] = cla * (angle - alpha_zero_lift)
        elif nonlinear:
            stall_fraction = ((angle - alpha_clmax) /
                              (alpha_stall_onset - alpha_clmax))
            cl[index] = (poly_a0 + poly_a1 * (angle - alpha_clmax) +
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
    for index, value in enumerate(lift):
        if not np.isfinite(value) or abs(abs(value) - UNUSED) <= 1.0e-40:
            continue
        clb[index] = (supplied[index] if supplied[index] is not None
                      else clb_per_cl * value)

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
