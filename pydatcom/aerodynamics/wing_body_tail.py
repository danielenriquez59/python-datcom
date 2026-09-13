"""
Wing-body-tail coefficient buildup for PyDATCOM.

Translates the routines that combine the wing-body result with the
horizontal-tail load after downwash and dynamic-pressure loss:

- ``CLWBT``: wing-body-tail lift and lift-curve slope.
- ``CDWBT``: wing-body-tail drag, DATCOM Section 4.5.3.2.

The lift carryover factors ``KHB``/``KBH`` (tail-on-body and body-on-tail,
DATCOM Section 4.3.1.2) and their incidence counterparts ``KKHB``/``KKBH``
are produced by the not-yet-translated ``HBTRAN``/``SUPHB`` path, so they
are inputs here.  Their no-carryover values (1 and 0) reduce the buildup to
an isolated tail, which is a marked approximation rather than a translation
of the interference routines.

Reference: datcom-legacy/datcom_2000/clwbt.f, cdwbt.f
"""

import numpy as np
from typing import Dict, Optional
import logging

logger = logging.getLogger(__name__)


def calculate_clwbt(cl_wing_body: float, cl_tail_at_alpt: float,
                    cla_tail: float, alih_deg: float,
                    eps_deg: float, qoqi: float,
                    khb: float = 1.0, kbh: float = 0.0,
                    kkhb: float = 1.0, kkbh: float = 0.0,
                    cla_wing_body: Optional[float] = None,
                    deda: Optional[float] = None) -> Dict[str, float]:
    """Translate CLWBT's wing-body-tail lift.

    ``CLBWH = ((CLH-CLI)*(KHB+KBH) + CLI*(KKHB+KKBH))*QOQI + CLBW`` with
    ``CLI = CLAH*ALIH``, and
    ``CLABWH = CLABW + (KHB+KBH)*CLAH*(1-DEODA)*QOQI``.

    The source additionally carries a ``FACT(141)*FACT(J+141)*FACT(J+161)``
    vortex/flap term onto ``CLAH*ALPT``.  That term comes from the
    untranslated ``WHWB`` path and is omitted here; it is zero for a clean
    configuration with no deflected surfaces.

    Args:
        cl_wing_body: ``CLBW``, wing-body lift coefficient on SREF.
        cl_tail_at_alpt: ``CLH`` evaluated at the tail local angle ALPT.
        cla_tail: ``CLAH``, isolated tail lift-curve slope, per degree.
        alih_deg: ``ALIH``, tail incidence, degrees.
        eps_deg: Downwash angle at the tail, degrees.
        qoqi: ``QOQI``, tail dynamic-pressure ratio.
        khb, kbh: Lift carryover factors for the angle-of-attack load.
        kkhb, kkbh: Lift carryover factors for the incidence load.
        cla_wing_body: ``CLABW``, wing-body slope per degree; enables the
            ``CLABWH`` output.
        deda: ``DEODA``, downwash gradient; required with ``cla_wing_body``.

    Returns:
        Dictionary with ``cl_total``, the tail increment, and ``cla_total``
        when the slope inputs are supplied.
    """
    if qoqi < 0.0:
        raise ValueError("dynamic-pressure ratio cannot be negative")
    cli = cla_tail * alih_deg
    cl_tail_term = ((cl_tail_at_alpt - cli) * (khb + kbh) +
                    cli * (kkhb + kkbh))
    cl_tail_increment = cl_tail_term * qoqi
    cl_total = cl_tail_increment + cl_wing_body

    result = {
        'cl_total': float(cl_total),
        'cl_wing_body': float(cl_wing_body),
        'cl_tail_increment': float(cl_tail_increment),
        'cl_tail_isolated': float(cl_tail_at_alpt),
        'cli': float(cli),
        'eps_deg': float(eps_deg),
        'qoqi': float(qoqi),
        'method': 'legacy_clwbt_no_vortex_term',
    }
    if cla_wing_body is not None and deda is not None:
        result['cla_total'] = float(
            cla_wing_body + (khb + kbh) * cla_tail * (1.0 - deda) * qoqi)
        result['cla_wing_body'] = float(cla_wing_body)
    return result


def calculate_cdwbt(cd_wing_body: float, cd_tail: float, cl_tail: float,
                    eps_deg: float, qoqi: float,
                    cdo_vertical: float = 0.0) -> Dict[str, float]:
    """Translate CDWBT, DATCOM Section 4.5.3.2.

    ``WBTCD = CDOV + CDWB + QOQI*(CDH*COS(EPS) + CLH*SIN(EPS))``

    The downwash rotates the tail force into the freestream axes, so the
    tail lift contributes to drag through ``sin(EPS)``.  This is a complete
    translation of the source expression.

    Args:
        cd_wing_body: ``CDWB``, wing-body drag coefficient on SREF.
        cd_tail: ``CDH``, tail drag coefficient.
        cl_tail: ``CLH``, tail lift coefficient.
        eps_deg: Downwash angle at the tail, degrees.
        qoqi: ``QOQI``, tail dynamic-pressure ratio.
        cdo_vertical: ``CDOV``, vertical-tail zero-lift drag.

    Returns:
        Dictionary with ``cd_total`` and its components.
    """
    if qoqi < 0.0:
        raise ValueError("dynamic-pressure ratio cannot be negative")
    eps_rad = np.deg2rad(eps_deg)
    cd_tail_increment = qoqi * (cd_tail * np.cos(eps_rad) +
                                cl_tail * np.sin(eps_rad))
    cd_total = cdo_vertical + cd_wing_body + cd_tail_increment
    return {
        'cd_total': float(cd_total),
        'cd_wing_body': float(cd_wing_body),
        'cd_tail_increment': float(cd_tail_increment),
        'cdo_vertical': float(cdo_vertical),
        'eps_deg': float(eps_deg),
        'method': 'legacy_cdwbt',
    }


def calculate_tail_load(state: Dict, alpha_deg: float,
                        downwash: Dict[str, float],
                        cla_tail: Optional[float] = None) -> Dict[str, float]:
    """Horizontal-tail lift coefficient at the local flow angle.

    ``CLWBT`` forms the tail angle as ``ALPT = ALPHA - EPS`` and looks the
    tail load up in its own lift curve.  With a linear tail curve this is
    ``CLH = CLAH*(ALPT + ALIH)``.  The tail slope is normalized to SREF so
    the result adds directly to the aircraft coefficient.

    Args:
        state: State dictionary with htail geometry and SYNTHS entries.
        alpha_deg: Aircraft angle of attack, degrees.
        downwash: Result of ``calculate_downwash``.
        cla_tail: Tail lift-curve slope per degree on the tail's own area;
            computed from the exposed tail planform when omitted.

    Returns:
        Dictionary with ``cl_tail``, ``alpt_deg`` and ``cla_tail_sref``.
    """
    from pydatcom.geometry.wing import calculate_straight_exposed_geometry

    tail = calculate_straight_exposed_geometry(state, component='htail')
    sref = float(state.get('options_sref', tail['area']) or tail['area'])
    if sref <= 0.0:
        raise ValueError("SREF must be positive")

    if cla_tail is None:
        from pydatcom.aerodynamics.lift import calculate_lift_curve_slope_compressible
        mach = float(state.get('flight_mach', 0.0) or 0.0)
        cla_tail = calculate_lift_curve_slope_compressible(
            tail['aspect_ratio'], tail['taper_ratio'], mach)

    alih = float(state.get('synths_alih', 0.0) or 0.0)
    alpt = float(alpha_deg) - downwash['eps_deg']
    cla_tail_sref = cla_tail * tail['area'] / sref
    cl_tail = cla_tail_sref * (alpt + alih)
    return {
        'cl_tail': float(cl_tail),
        'alpt_deg': float(alpt),
        'cla_tail': float(cla_tail),
        'cla_tail_sref': float(cla_tail_sref),
        'tail_area': tail['area'],
        'sref': sref,
    }
