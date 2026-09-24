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
    out = {}
    for start, degrees in ((106, savsi_deg), (112, savso_deg)):
        record = [float(previous.get(start + k, 0.0)) for k in range(6)]
        record[0] = float(degrees)
        for k, value in enumerate(angles(1, record)):
            out[start + k] = value
    return out


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
    a: Dict[int, float] = {int(k): float(v) for k, v in a_in.items()}
    ain = {int(k): float(v) for k, v in ain.items()}
    for key in range(1, 196):
        a.setdefault(key, 0.0)

    single = ain[2] < 10.0 * UNUSED
    a[21] = ain[4] - ain[2]
    a[23] = ain[3] - ain[2]
    a[19] = a[23] / a[21]
    if single:
        ain[5] = ain[1]
    a[25] = ain[5] / ain[6]
    a[10] = ain[6] * (a[25] + (1.0 - a[25]) * a[19])
    a[26] = ain[5] / a[10]
    a[28] = 1.0
    if ain[5] != 0.0:
        a[28] = ain[1] / ain[5]
    a[27] = a[26] * a[28]
    a[118] = ain[1] / ain[6]
    a[1] = (a[10] + ain[5]) * a[23]
    a[2] = (ain[5] + ain[1]) * ain[2]
    a[3] = a[1] + a[2]
    a[119] = (ain[6] + ain[5]) * a[21]
    a[4] = a[119] + a[2]
    a[5] = _eq1(a[23], a[1])
    a[6] = 1.0
    if not single:
        a[6] = _eq1(ain[2], a[2])
    a[7] = _eq1(ain[3], a[3])
    a[120] = _eq1(ain[4], a[4])
    ci = _eq2(a[26], a[5])
    co = 0.0 if single else _eq2(a[28], a[6])

    def record(start):
        return [a[start + k] for k in range(6)]

    def store(start, values):
        for k, value in enumerate(values):
            a[start + k] = float(value)

    xovco = ain.get(66, 0.0)
    for i in range(5):
        if i < 4:
            exposed = 34 + 6 * i
            inboard, outboard = exposed + 24, exposed + 48
            tmpo = tmpi = ain[9] - _FN[i]
        else:
            exposed, inboard, outboard = 175, 187, 181
            tmpi = ain[9] - a[174]
            tmpo = ain[9] - xovco
        if abs(tmpi) < 1.0e-6:
            # Label 1060: this station is the input station.
            store(inboard, record(106))
            store(outboard, record(112))
            combine = not single
        else:
            a[inboard + 4] = ci * tmpi + a[110]
            store(inboard, angles(5, record(inboard)))
            if single:
                store(outboard, zerang())
                combine = False
            else:
                a[outboard + 4] = co * tmpo + a[116]
                store(outboard, angles(5, record(outboard)))
                combine = True
        if combine:
            a[exposed + 3] = ((a[1] * a[inboard + 3] + a[2] * a[outboard + 3])
                              / a[3])
            store(exposed, angles(4, record(exposed)))
        else:
            store(exposed, record(inboard))

    a[15] = 2.0 * a[10] * (1.0 + a[26] + a[26]**2) / (3.0 * (1.0 + a[26]))
    a[121] = 2.0 * ain[6] * (1.0 + a[25] + a[25]**2) / (3.0 * (1.0 + a[25]))
    a[17] = 2.0 * ain[5] * (1.0 + a[28] + a[28]**2) / (3.0 * (1.0 + a[28]))
    a[16] = (a[1] * a[15] + a[2] * a[17]) / a[3]
    a[122] = (a[119] * a[121] + a[2] * a[17]) / a[4]
    a[32] = a[23] * (1.0 + 2.0 * a[26]) / (3.0 * (1.0 + a[26]))
    a[33] = ain[2] * (1.0 + 2.0 * a[28]) / (3.0 * (1.0 + a[28])) + a[23]
    a[31] = (a[1] * a[32] + a[2] * a[33]) / a[3]
    a[30] = a[16] / 2.0 + (a[1] * a[32] * a[62] + a[2] * (
        a[23] * a[62] + (a[33] - a[23]) * a[86])) / a[3]
    a[18] = (a[23] * a[62] + ain[2] * a[86]) / a[10]
    le_inboard = a[23] * a[62]
    le_tip = ain[2] * a[86] + le_inboard
    te_inboard = a[23] * a[80]
    te_tip = ain[2] * a[104] + te_inboard
    a[29] = (a[10] + max(0.0, te_tip, te_inboard) -
             min(0.0, le_tip, le_inboard))
    a[162] = a[23] / ain[3]
    a[163] = 4.0 * (a[21]**2) / a[119]
    a[164] = a[23] / 2.0
    a[165] = a[164] + ain[2]
    if not single:
        a[166] = ain[1] + a[165] * ((ain[5] - ain[1]) / ain[2])
        a[167] = (a[166] + ain[1]) * a[165]
        a[168] = 4.0 * (a[165]**2) / a[167]
        a[169] = ain[1] / a[166]
    a[130] = a[21] * (1.0 + 2.0 * a[25]) / (3.0 * (1.0 + a[25]))
    a[133] = 0.0
    if not single:
        a[133] = ain[2] * (1.0 + 2.0 * a[28]) / (3.0 * (1.0 + a[28])) + a[21]
    a[136] = (a[119] * a[130] + a[2] * a[133]) / a[4]
    a[195] = (a[119] * a[130] * a[62] + a[2] * (
        a[21] * a[62] + (a[133] - a[21]) * a[86])) / a[4]
    a[161] = a[122] / 4.0 + a[195]
    return {'a': a, 'ain': ain, 'single_panel': single,
            'method': 'legacy_wtgeom'}


def vertical_panel_adjustments(a: Mapping[int, float],
                               z_offset: float) -> Dict[int, float]:
    """M02O02's treatment of a vertical panel's ``A`` block after WTGEOM.

    WTGEOM treats the panel as a two-sided surface; M02O02 halves its
    areas and aspect ratios (every listed slot not ``UNUSED``) and adds the
    panel's vertical offset ``ZV`` or ``ZVF`` to the MAC stations ``A(130)``,
    ``A(133)`` and ``A(136)``.
    """
    out = {int(k): float(v) for k, v in a.items()}
    for index in _HALVED:
        if out.get(index, 0.0) != UNUSED:
            out[index] = out.get(index, 0.0) / 2.0
    for index in (130, 133, 136):
        out[index] = out.get(index, 0.0) + z_offset
    return out


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
