"""
Wing-body and wing-body-tail dynamic derivative buildup.

``DNPAWB`` combines the isolated wing and body dynamic derivatives through
the Section 4.3.1.2 carryover factors.  ``DNPWBT`` then adds the horizontal
tail, using the dynamic-pressure ratio and downwash gradient that ``DWASH``
supplies, and splits on the wing-to-tail span ratio.

Both routines are pure combination: every quantity they consume is produced
upstream, so these translations take them as arguments rather than emulating
the COMMON blocks.

Reference: datcom-legacy/datcom_2000/dnpawb.f, dnpwbt.f
"""

import numpy as np
from typing import Dict, Optional, Sequence
import logging

logger = logging.getLogger(__name__)

# DNPAWB marks an unavailable transonic CMAD with this sentinel rather than
# the usual UNUSED, and returns before combining when it is present.
_CMAD_UNAVAILABLE = 1000.0

# DNPWBT switches formulation at this wing-to-tail span ratio.
_SPAN_RATIO_SPLIT = 1.5


def calculate_dnpawb(clq_wing: float, cmq_wing: float,
                     clad_wing: float, cmad_wing: float,
                     clq_body: float, cmq_body: float,
                     clad_body: float, cmad_body: float,
                     khb: float, kbh: float,
                     beta_aspect_ratio: Optional[float] = None,
                     taper_ratio: float = 0.0,
                     supersonic: bool = False,
                     hypersonic: bool = False,
                     transonic: bool = False) -> Dict[str, float]:
    """Translate DNPAWB: wing-body dynamic derivatives.

    ``CLQWB = (KWB+KBW)*CLQW + CLQB`` and likewise for the other three, with
    the carryover applied to the wing part only.

    Args:
        clq_wing, cmq_wing, clad_wing, cmad_wing: Isolated wing derivatives.
        clq_body, cmq_body, clad_body, cmad_body: Isolated body derivatives.
        khb, kbh: Section 4.3.1.2 carryover factors, ``KWB`` and ``KBW``.
        beta_aspect_ratio: ``BETA*ASTRW``, gating the acceleration
            derivatives outside supersonic flow.
        taper_ratio: Exposed taper ratio ``LAMDA``; a nonzero value also
            gates the acceleration derivatives outside supersonic flow.
        supersonic, hypersonic, transonic: Regime selectors.

    Returns:
        Dictionary with ``clq`` and ``cmq`` always, and ``clad`` and ``cmad``
        when the source's gates allow them.  ``acceleration_available`` says
        which happened.

    Notes:
        The source zeroes the body acceleration derivatives in hypersonic
        flow and marks the transonic wing ``CMAD`` unavailable with a literal
        1000.  Both are preserved.
    """
    carryover = khb + kbh
    if hypersonic:
        clad_body = 0.0
        cmad_body = 0.0

    result = {
        'clq': float(carryover * clq_wing + clq_body),
        'cmq': float(carryover * cmq_wing + cmq_body),
        'carryover': float(carryover),
        'method': 'legacy_dnpawb',
    }

    # Outside supersonic flow: acceleration terms only for BETA*AR in [0, 4]
    # on an untapered wing.
    if not supersonic:
        beta_ar = beta_aspect_ratio
        if beta_ar is not None and not 0.0 <= beta_ar <= 4.0:
            result['acceleration_available'] = False
            result['acceleration_gate'] = 'beta_aspect_ratio_out_of_range'
            return result
        if taper_ratio != 0.0:
            result['acceleration_available'] = False
            result['acceleration_gate'] = 'tapered_wing'
            return result

    result['clad'] = float(carryover * clad_wing + clad_body)
    if transonic and cmad_wing == _CMAD_UNAVAILABLE:
        result['acceleration_available'] = False
        result['acceleration_gate'] = 'transonic_cmad_unavailable'
        return result

    result['cmad'] = float(carryover * cmad_wing + cmad_body)
    result['acceleration_available'] = True
    return result


def calculate_dnpwbt(clq_wing_body: float, cmq_wing_body: float,
                     clad_wing_body: float, cmad_wing_body: float,
                     qoqi: Sequence[float], deda: Sequence[float],
                     cla_tail: float, khb: float, kbh: float,
                     dxac: float, cbar: float,
                     wing_span: float, tail_span: float,
                     jet_term: Optional[Sequence[float]] = None,
                     taper_ratio: float = 0.0,
                     transonic: bool = False) -> Dict[str, object]:
    """Translate DNPWBT: wing-body-tail dynamic derivatives.

    The tail increment is ``2*(KBH+KHB)*(q/q_inf)*CLAH`` scaled by the tail
    arm in reference chords.  The source uses two formulations depending on
    the wing-to-tail span ratio:

    - ``BW/BH >= 1.5``: the arm is folded into the increment before it is
      applied, and the acceleration derivatives carry the downwash gradient.
    - ``BW/BH < 1.5``: the arm is applied outside, the jet term is added to
      the pitching increment, and the acceleration derivatives come from the
      jet term alone with the arm doubled.

    Args:
        clq_wing_body, cmq_wing_body, clad_wing_body, cmad_wing_body:
            Wing-body derivatives from DNPAWB.
        qoqi: Dynamic-pressure ratio at the tail, per angle.
        deda: Downwash gradient at the tail, per angle.
        cla_tail: Tail lift-curve slope, ``CLAH``.
        khb, kbh: Tail-body carryover factors.
        dxac: Tail aerodynamic-centre offset from the CG, ``DXAC``.
        cbar: Reference chord.
        wing_span: ``SBW``.
        tail_span: ``BH``.
        jet_term: ``DTJ`` per angle; zero when absent, as in transonic flow.
        taper_ratio: ``LAMDA``; a tapered wing suppresses the small-span-
            ratio acceleration derivatives.
        transonic: Suppresses the small-span-ratio acceleration branch.

    Returns:
        Dictionary with per-angle ``clq``, ``cmq`` and, where the source
        computes them, ``clad`` and ``cmad``; plus the span ratio and which
        branch ran.

    Raises:
        ValueError: If the spans or reference chord are nonpositive, or the
            per-angle arrays disagree in length.
    """
    qoqi = np.asarray(qoqi, dtype=float)
    deda = np.asarray(deda, dtype=float)
    if qoqi.shape != deda.shape or qoqi.ndim != 1 or qoqi.size == 0:
        raise ValueError(
            "DNPWBT needs matching nonempty q/q and de/da arrays",
        )
    if min(wing_span, tail_span, cbar) <= 0.0:
        raise ValueError("DNPWBT requires positive spans and reference chord")

    if jet_term is None:
        jet = np.zeros_like(qoqi)
    else:
        jet = np.asarray(jet_term, dtype=float)
    if jet.shape != qoqi.shape:
        raise ValueError("DNPWBT jet term must match the angle schedule")

    span_ratio = wing_span / tail_span
    tail_arm_chords = -dxac / cbar
    carryover = kbh + khb

    result = {
        'span_ratio': float(span_ratio),
        'tail_arm_chords': float(tail_arm_chords),
        'method': 'legacy_dnpwbt',
    }

    if span_ratio >= _SPAN_RATIO_SPLIT:
        # BW/BH >= 1.5: arm folded into the tail increment (SAVTIM / SAVTIX).
        tail_clq_increment = (
            2.0 * carryover * qoqi * tail_arm_chords * cla_tail
        )
        tail_cmq_increment = tail_clq_increment * tail_arm_chords
        result['branch'] = 'wide_wing'
        result['clq'] = clq_wing_body + tail_clq_increment
        result['cmq'] = cmq_wing_body - tail_cmq_increment
        result['clad'] = clad_wing_body + tail_clq_increment * deda
        result['cmad'] = cmad_wing_body - tail_cmq_increment * deda
        result['acceleration_available'] = True
        return result

    # BW/BH < 1.5: arm applied outside; jet term enters CMQ.
    tail_force_increment = 2.0 * carryover * qoqi * cla_tail
    tail_plus_jet = tail_force_increment + jet
    result['branch'] = 'narrow_wing'
    result['clq'] = clq_wing_body + tail_plus_jet * tail_arm_chords
    result['cmq'] = (
        cmq_wing_body - tail_plus_jet * tail_arm_chords ** 2
    )

    if transonic or taper_ratio != 0.0:
        result['acceleration_available'] = False
        result['acceleration_gate'] = (
            'transonic' if transonic else 'tapered_wing'
        )
        return result

    doubled_arm = tail_arm_chords * 2.0
    result['clad'] = clad_wing_body - doubled_arm * jet
    result['cmad'] = cmad_wing_body + (doubled_arm ** 2 / 2.0) * jet
    result['acceleration_available'] = True
    return result
