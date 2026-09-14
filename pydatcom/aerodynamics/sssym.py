"""
SSSYM: supersonic pitching moment and lift increments for trailing-edge flaps.

The first control-surface routine translated.  Pitching moment
effectiveness is the sum of three contributions, each scaled by the flap
span and normalised on the reference chord and area:

1. a hinge-line term carrying the flap chord
2. a sweep term through the hinge-line tangent
3. a moment-arm term from the flap's distance ahead of the CG

Lift effectiveness is a single term.  Both are then multiplied by each
deflection in the schedule.

The source computes nothing unless the flap type is trailing edge, and
reads its reference chord, leading-edge sweep and longitudinal position from
the horizontal tail when one is present and from the wing otherwise.

Reference: datcom-legacy/datcom_2000/sssym.f
"""

import numpy as np
from typing import Dict, Sequence
import logging

logger = logging.getLogger(__name__)

# The source's FTYPE selector for a trailing-edge flap.
TRAILING_EDGE = 1.0


def calculate_sssym(deflections: Sequence[float],
                    root_chord: float, sref: float,
                    span_inboard: float, span_outboard: float,
                    flap_chord: float, flap_area: float,
                    k3: float, taper_ratio_flap: float,
                    tan_hinge_line: float, tan_le: float,
                    x_surface: float, chord_inboard: float, xcg: float,
                    bcmd1: float, bcld1: float, bcld2: float,
                    flap_type: float = TRAILING_EDGE) -> Dict[str, object]:
    """Translate SSSYM: supersonic trailing-edge flap effectiveness.

    Args:
        deflections: Deflection schedule, the source's ``DELSYM``.
        root_chord: Reference chord ``CR``, from the tail when one is
            present and the wing otherwise.
        sref: Reference area.
        span_inboard: Inboard flap station, ``ALOCI``.
        span_outboard: Outboard flap station, ``ALOCO``.
        flap_chord: Flap chord ``CFI``.
        flap_area: Flap area ``SF``.
        k3: The source's ``K3`` effectiveness factor.
        taper_ratio_flap: ``TRTOFL``, entering the first moment term.
        tan_hinge_line: ``TANHL``, the hinge-line sweep tangent.
        tan_le: Leading-edge sweep tangent of the carrying surface.
        x_surface: Longitudinal position of that surface, ``XW`` or ``XH``.
        chord_inboard: ``CI``, the chord at the inboard flap station.
        xcg: Moment reference station.
        bcmd1: Moment effectiveness parameter.
        bcld1: First lift effectiveness parameter.
        bcld2: Second lift effectiveness parameter.
        flap_type: ``FTYPE``; anything but trailing edge produces nothing.

    Returns:
        Dictionary with ``cmd_total`` and ``cld`` effectiveness, the three
        moment contributions, and the per-deflection increments.  When the
        flap type is not trailing edge, ``applicable`` is False and the
        increments are absent.

    Raises:
        ValueError: If the reference chord or area is nonpositive.
    """
    if flap_type != TRAILING_EDGE:
        return {
            'applicable': False,
            'reason': 'the source computes flap increments for trailing '
                      'edge flaps only',
            'method': 'legacy_sssym',
        }
    if min(root_chord, sref) <= 0.0:
        raise ValueError("SSSYM requires a positive reference chord and area")

    deflections = np.atleast_1d(np.asarray(deflections, dtype=float))

    k1 = k3 * (1.0 + taper_ratio_flap + taper_ratio_flap**2)
    k2 = k3 * tan_hinge_line

    flap_span = 2.0 * (span_outboard - span_inboard)      # BEF
    normaliser = root_chord * sref                        # SAVE

    cmd1 = k1 * (1.0 / 3.0) * flap_span * flap_chord * bcmd1 / normaliser
    cmd2 = -k2 / 2.0 * flap_span * flap_area * bcld2 / normaliser
    # The flap's moment arm: its hinge line forward of the CG.
    arm = x_surface + span_inboard * tan_le + chord_inboard - flap_chord - xcg
    cmd3 = -k3 * arm * flap_area * bcld1 / normaliser
    cmd_total = cmd1 + cmd2 + cmd3

    cld = k3 * bcld1 * flap_area / sref

    return {
        'applicable': True,
        'cmd_total': float(cmd_total),
        'cmd1': float(cmd1),
        'cmd2': float(cmd2),
        'cmd3': float(cmd3),
        'cld': float(cld),
        'delta_cm': cmd_total * deflections,
        'delta_cl': cld * deflections,
        'flap_span': float(flap_span),
        'moment_arm': float(arm),
        'k1': float(k1),
        'k2': float(k2),
        'method': 'legacy_sssym',
    }
