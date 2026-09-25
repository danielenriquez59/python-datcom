"""
WTGEOM and SETUP1: complete lifting-surface planform geometry.

WTGEOM fills a surface's ``A`` block from its input block for every
planform, not only the straight tapered one: exposed, theoretical, inboard
and outboard areas, aspect ratios, tapers and mean aerodynamic chords, and
six-word ANGLES records (degrees, radians, sine, cosine, tangent, test) of
the sweep at the leading edge, quarter chord, half chord, trailing edge and
maximum thickness for the inboard, outboard and exposed panels.

The record positions, one-based as in the source:

=========  ========  ========  ========
station    exposed   inboard   outboard
=========  ========  ========  ========
LE         34        58        82
c/4        40        64        88
c/2        46        70        94
TE         52        76        100
max t/c    175       187       181
=========  ========  ========  ========

The block is returned as a mapping keyed by those indices, ``a[120]`` for
``A(120)``, since the rest of the program addresses it that way.
``pydatcom.geometry.wing.calculate_straight_exposed_geometry`` remains the
named-key view of the straight subset, and agrees with this.

SETUP1 completes the input sweep records ``A(106)`` (inboard) and ``A(112)``
(outboard) that WTGEOM starts from, and M02O02 halves a vertical panel's
areas and aspect ratios and offsets its vertical stations, which
:func:`vertical_panel_adjustments` translates.

Reference: datcom-legacy/datcom_2000/wtgeom.f, setup1.f, m02o02.f
"""

import numpy as np
from typing import Dict, Mapping
import logging

from pydatcom.utils.constants import UNUSED
from pydatcom.utils.legacy_numeric import angles, zerang

logger = logging.getLogger(__name__)

# Chord stations of the first four sweep records: LE, c/4, c/2, TE.
_FN = (0.0, 0.25, 0.50, 1.0)

# M02O02's INDX: the A slots halved for a vertical panel.
_HALVED = (1, 2, 3, 4, 5, 6, 7, 119, 120, 163, 167, 168)


def _eq1(a1: float, a2: float) -> float:
    """``4*a1**2/a2``: aspect ratio from a semispan and an area."""
    return 4.0 * a1**2 / a2


def _eq2(a1: float, a2: float) -> float:
    """``4*(1-a1)/(a2*(1+a1))``: the sweep-station factor."""
    return 4.0 * (1.0 - a1) / (a2 * (1.0 + a1))


def sweep_records(savsi_deg: float, savso_deg: float,
                  previous: Mapping[int, float] = None) -> Dict[int, float]:
    """Translate SETUP1's ``ANGLES(1,A(106))`` and ``ANGLES(1,A(112))``.

    Args:
        savsi_deg, savso_deg: The input inboard and outboard sweeps at the
            reference chord station, ``A(106)`` and ``A(112)``.
        previous: ``A(106)`` to ``A(117)`` as they stood, whose test words
            ANGLES consults; zero by default.

    Returns:
        ``{106: ..., 117: ...}``, the two completed records.
    """
    previous = previous or {}
    completed = {}
    for start, degrees in ((106, savsi_deg), (112, savso_deg)):
        record = [float(previous.get(start + word_offset, 0.0))
                  for word_offset in range(6)]
        record[0] = float(degrees)
        for word_offset, value in enumerate(angles(1, record)):
            completed[start + word_offset] = value
    return completed


def calculate_wtgeom(ain: Mapping[int, float],
                     a_in: Mapping[int, float]) -> Dict[str, object]:
    """Translate WTGEOM for any planform.

    Args:
        ain: The surface's input block, one-based: ``1`` CHRDTP, ``2``
            SSPNOP, ``3`` SSPNE, ``4`` SSPN, ``5`` CHRDBP, ``6`` CHRDR,
            ``9`` CHSTAT, ``66`` XOVCO.
        a_in: The ``A`` block as WTGEOM finds it: the sweep records
            ``106``-``117`` (see :func:`sweep_records`), ``174`` (the
            inboard maximum-thickness station, SETUP1's ``XOVC``), and any
            earlier contents of the record slots, whose ANGLES test words
            decide whether a record is recomputed.

    Returns:
        ``{'a': ..., 'ain': ...}``: the full ``A`` block after WTGEOM, and
        the input block, whose ``CHRDBP`` WTGEOM sets to ``CHRDTP`` for a
        single-panel surface.

    Notes:
        For a two-panel surface the exposed sweep at each station is formed
        from the area-weighted cosine of the two panels and resolved with
        ``ANGLES(4)``, which takes the positive sine, so it is never
        negative: a forward-swept two-panel surface loses the sign of its
        combined sweep.  A single-panel surface copies the inboard record
        and keeps it.

        When ``CHSTAT`` coincides with a station the input records are
        copied rather than recomputed.  At the maximum-thickness station
        the test looks at the inboard station ``A(174)`` only, so the
        outboard record is copied from the input even when the outboard
        station ``XOVCO`` differs.
    """
    a_block: Dict[int, float] = {int(k): float(v) for k, v in a_in.items()}
    surface_in = {int(k): float(v) for k, v in ain.items()}
    for a_word in range(1, 196):
        a_block.setdefault(a_word, 0.0)

    single = surface_in[2] < 10.0 * UNUSED
    a_block[21] = surface_in[4] - surface_in[2]
    a_block[23] = surface_in[3] - surface_in[2]
    a_block[19] = a_block[23] / a_block[21]
    if single:
        surface_in[5] = surface_in[1]
    a_block[25] = surface_in[5] / surface_in[6]
    a_block[10] = surface_in[6] * (a_block[25] + (1.0 - a_block[25]) * a_block[19])
    a_block[26] = surface_in[5] / a_block[10]
    a_block[28] = 1.0
    if surface_in[5] != 0.0:
        a_block[28] = surface_in[1] / surface_in[5]
    a_block[27] = a_block[26] * a_block[28]
    a_block[118] = surface_in[1] / surface_in[6]
    a_block[1] = (a_block[10] + surface_in[5]) * a_block[23]
    a_block[2] = (surface_in[5] + surface_in[1]) * surface_in[2]
    a_block[3] = a_block[1] + a_block[2]
    a_block[119] = (surface_in[6] + surface_in[5]) * a_block[21]
    a_block[4] = a_block[119] + a_block[2]
    a_block[5] = _eq1(a_block[23], a_block[1])
    a_block[6] = 1.0
    if not single:
        a_block[6] = _eq1(surface_in[2], a_block[2])
    a_block[7] = _eq1(surface_in[3], a_block[3])
    a_block[120] = _eq1(surface_in[4], a_block[4])
    ci = _eq2(a_block[26], a_block[5])
    co = 0.0 if single else _eq2(a_block[28], a_block[6])

    def record(start):
        return [a_block[start + word_offset] for word_offset in range(6)]

    def store(start, values):
        for word_offset, value in enumerate(values):
            a_block[start + word_offset] = float(value)

    xovco = surface_in.get(66, 0.0)
    for pass_index in range(5):
        if pass_index < 4:
            exposed = 34 + 6 * pass_index
            inboard, outboard = exposed + 24, exposed + 48
            outboard_delta = inboard_delta = surface_in[9] - _FN[pass_index]
        else:
            exposed, inboard, outboard = 175, 187, 181
            inboard_delta = surface_in[9] - a_block[174]
            outboard_delta = surface_in[9] - xovco
        if abs(inboard_delta) < 1.0e-6:
            # Label 1060: this station is the input station.
            store(inboard, record(106))
            store(outboard, record(112))
            combine = not single
        else:
            a_block[inboard + 4] = ci * inboard_delta + a_block[110]
            store(inboard, angles(5, record(inboard)))
            if single:
                store(outboard, zerang())
                combine = False
            else:
                a_block[outboard + 4] = co * outboard_delta + a_block[116]
                store(outboard, angles(5, record(outboard)))
                combine = True
        if combine:
            a_block[exposed + 3] = ((a_block[1] * a_block[inboard + 3] + a_block[2] * a_block[outboard + 3])
                              / a_block[3])
            store(exposed, angles(4, record(exposed)))
        else:
            store(exposed, record(inboard))

    a_block[15] = 2.0 * a_block[10] * (1.0 + a_block[26] + a_block[26]**2) / (3.0 * (1.0 + a_block[26]))
    a_block[121] = 2.0 * surface_in[6] * (1.0 + a_block[25] + a_block[25]**2) / (3.0 * (1.0 + a_block[25]))
    a_block[17] = 2.0 * surface_in[5] * (1.0 + a_block[28] + a_block[28]**2) / (3.0 * (1.0 + a_block[28]))
    a_block[16] = (a_block[1] * a_block[15] + a_block[2] * a_block[17]) / a_block[3]
    a_block[122] = (a_block[119] * a_block[121] + a_block[2] * a_block[17]) / a_block[4]
    a_block[32] = a_block[23] * (1.0 + 2.0 * a_block[26]) / (3.0 * (1.0 + a_block[26]))
    a_block[33] = surface_in[2] * (1.0 + 2.0 * a_block[28]) / (3.0 * (1.0 + a_block[28])) + a_block[23]
    a_block[31] = (a_block[1] * a_block[32] + a_block[2] * a_block[33]) / a_block[3]
    a_block[30] = a_block[16] / 2.0 + (a_block[1] * a_block[32] * a_block[62] + a_block[2] * (
        a_block[23] * a_block[62] + (a_block[33] - a_block[23]) * a_block[86])) / a_block[3]
    a_block[18] = (a_block[23] * a_block[62] + surface_in[2] * a_block[86]) / a_block[10]
    le_inboard = a_block[23] * a_block[62]
    le_tip = surface_in[2] * a_block[86] + le_inboard
    te_inboard = a_block[23] * a_block[80]
    te_tip = surface_in[2] * a_block[104] + te_inboard
    a_block[29] = (a_block[10] + max(0.0, te_tip, te_inboard) -
             min(0.0, le_tip, le_inboard))
    a_block[162] = a_block[23] / surface_in[3]
    a_block[163] = 4.0 * (a_block[21]**2) / a_block[119]
    a_block[164] = a_block[23] / 2.0
    a_block[165] = a_block[164] + surface_in[2]
    if not single:
        a_block[166] = surface_in[1] + a_block[165] * ((surface_in[5] - surface_in[1]) / surface_in[2])
        a_block[167] = (a_block[166] + surface_in[1]) * a_block[165]
        a_block[168] = 4.0 * (a_block[165]**2) / a_block[167]
        a_block[169] = surface_in[1] / a_block[166]
    a_block[130] = a_block[21] * (1.0 + 2.0 * a_block[25]) / (3.0 * (1.0 + a_block[25]))
    a_block[133] = 0.0
    if not single:
        a_block[133] = surface_in[2] * (1.0 + 2.0 * a_block[28]) / (3.0 * (1.0 + a_block[28])) + a_block[21]
    a_block[136] = (a_block[119] * a_block[130] + a_block[2] * a_block[133]) / a_block[4]
    a_block[195] = (a_block[119] * a_block[130] * a_block[62] + a_block[2] * (
        a_block[21] * a_block[62] + (a_block[133] - a_block[21]) * a_block[86])) / a_block[4]
    a_block[161] = a_block[122] / 4.0 + a_block[195]
    return {'a': a_block, 'ain': surface_in, 'single_panel': single,
            'method': 'legacy_wtgeom'}


def vertical_panel_adjustments(a: Mapping[int, float],
                               z_offset: float) -> Dict[int, float]:
    """M02O02's treatment of a vertical panel's ``A`` block after WTGEOM.

    WTGEOM treats the panel as a two-sided surface; M02O02 halves its
    areas and aspect ratios (every listed slot not ``UNUSED``) and adds the
    panel's vertical offset ``ZV`` or ``ZVF`` to the MAC stations ``A(130)``,
    ``A(133)`` and ``A(136)``.
    """
    adjusted = {int(word): float(value) for word, value in a.items()}
    for word in _HALVED:
        if adjusted.get(word, 0.0) != UNUSED:
            adjusted[word] = adjusted.get(word, 0.0) / 2.0
    for word in (130, 133, 136):
        adjusted[word] = adjusted.get(word, 0.0) + z_offset
    return adjusted


def calculate_setup1(alpha_deg, wing_incidence: float, tail_incidence: float,
                     body_alpha_zero: float, xovc: float,
                     savsi_deg: float, savso_deg: float,
                     previous: Mapping[int, float] = None) -> Dict[str, object]:
    """Translate SETUP1 for one lifting surface's block and the schedules.

    SETUP1 sets ``A(174) = XOVC`` and completes the input sweep records
    (:func:`sweep_records`) for each surface, and forms the local angle
    schedules: ``B(23)`` onward ``= FLC + ALIW``, ``BHT(23)`` onward
    ``= FLC + ALIH``, and the body's ``BD(255)`` onward ``= FLC + BD(81)``.
    It also completes an ``A(138)`` record that no routine reads, which is
    omitted.
    """
    alpha = np.asarray(alpha_deg, dtype=float)
    records = sweep_records(savsi_deg, savso_deg, previous)
    records[174] = float(xovc)
    return {
        'a': records,
        'wing_local_alpha': alpha + float(wing_incidence),
        'tail_local_alpha': alpha + float(tail_incidence),
        'body_alpha': alpha + float(body_alpha_zero),
        'method': 'legacy_setup1',
    }
