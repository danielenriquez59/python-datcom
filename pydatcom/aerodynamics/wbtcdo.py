"""
WBTCDO: transonic wing-body-tail zero-lift drag, DATCOM Section 4.5.3.1.

Figure 4.5.3.1-19 comes in four parts that chain into one another, each
stage's output becoming the next stage's second independent variable:

- **A** sweep against thickness ratio
- **B** aspect ratio against the part A result
- **C** taper ratio against the part B result
- **D** configuration type against the part C result, giving the drag
  divergence Mach number ``MD``

A polynomial is then faired through the transonic rise, anchored on the
subsonic and supersonic drag levels the caller supplies.

Reference: datcom-legacy/datcom_2000/wbtcdo.f
"""

import numpy as np
from typing import Dict, Optional
import logging

from pydatcom.utils.legacy_interp import interx

logger = logging.getLogger(__name__)

# Part A: sweep (7) against thickness ratio (8).  LIND=8.
_PARMA = [0., 10., 20., 30., 40., 50., 60., 0.,
          .03, .04, .05, .06, .07, .08, .10, .12]
_XUNIT = [6.8, 7.25, 7.75, 8.40, 9.40, 10.8, 12.8,
          5., 5.45, 6.05, 6.8, 7.8, 9.25, 11.15,
          4.25, 4.75, 5.3, 6.1, 7.1, 8.7, 10.6,
          3.55, 3.9, 4.45, 5.3, 6.5, 8.1, 10.05,
          3.05, 3.4, 3.9, 4.75, 5.85, 7.6, 9.7,
          2.5, 2.8, 3.2, 4., 5.2, 6.95, 9.25,
          1.5, 1.8, 2.45, 3.25, 4.5, 6.3, 8.75,
          .75, 1.15, 1.75, 2.6, 3.85, 5.75, 8.20]

# Part B: aspect ratio (5) against the part A result (10).  LIND=10.
_PARMB = [2., 3., 4., 6., 8., 0., 0., 0., 0., 0.,
          0., 1., 2., 3., 4., 5., 6., 7., 8., 9.]
_YUNIT = [3.2, 2.5, 1.95, 1., .4,
          4.25, 3.75, 3., 2.1, 1.5,
          5.4, 4.8, 4., 3.25, 2.7,
          6.5, 5.9, 5.1, 4.4, 3.8,
          7.6, 7., 6.25, 5.55, 5.,
          8.6, 8., 7.4, 6.9, 6.3,
          9.55, 9.1, 8.7, 8.25, 7.7,
          10.45, 10.15, 9.9, 9.65, 9.25,
          11.15, 11.15, 11.15, 11.15, 11.15,
          12., 12.2, 12.7, 13., 14.25]

# Part C: taper ratio (3) against the part B result (3).  LIND=3.
_PARMC = [0., .5, 1., 0., 8., 16.]
_ZUNIT = [.13, .13, .13, 7.78, 8.18, 8.53, 15.43, 16.23, 16.93]

# Part D: configuration type (3) against the part C result (13).  LIND=13.
_PARMD = [1., 2., 3., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
          0., 1., 2., 3., 4., 5., 6., 7., 8., 9., 10., 11., 12.]
_TMD = [.755, .761, .77,
        .77, .78, .79,
        .788, .798, .81,
        .804, .815, .830,
        .82, .83, .85,
        .836, .848, .87,
        .85, .864, .89,
        .865, .88, .905,
        .879, .898, .925,
        .892, .912, .944,
        .905, .925, .96,
        .92, .942, .978,
        .932, .958, .992]

# The fairing anchors: the subsonic curve is referenced to Mach 0.7 and the
# supersonic one to Mach 1.1.
_SUBSONIC_ANCHOR = 0.7
_SUPERSONIC_ANCHOR = 1.1
# The source's substitution when the supersonic level is unknown.
_CD11_FROM_CD14 = 2.1381


def calculate_drag_divergence_mach(sweep_c4_deg: float,
                                   thickness_ratio: float,
                                   aspect_ratio: float, taper_ratio: float,
                                   configuration_type: float
                                   ) -> Dict[str, float]:
    """Figure 4.5.3.1-19: the drag divergence Mach number.

    Args:
        sweep_c4_deg: Quarter-chord sweep, degrees.
        thickness_ratio: ``TOC``.
        aspect_ratio: ``AR``.
        taper_ratio: ``TAPR``.
        configuration_type: ``ITYPE``, the body geometry selector.

    Returns:
        Dictionary with ``md`` and each intermediate stage.
    """
    stage_a = interx(2, _PARMA, [sweep_c4_deg, thickness_ratio], [7, 8],
                     _XUNIT, lind=8, lx1l=1, lx2l=2, lx1u=1, lx2u=2)
    stage_b = interx(2, _PARMB, [aspect_ratio, stage_a], [5, 10], _YUNIT,
                     lind=10, lx1l=1, lx2l=1, lx1u=1, lx2u=1)
    stage_c = interx(2, _PARMC, [taper_ratio, stage_b], [3, 3], _ZUNIT,
                     lind=3, lx1l=1, lx2l=1, lx1u=1, lx2u=1)
    md = interx(2, _PARMD, [configuration_type, stage_c], [3, 13], _TMD,
                lind=13, lx1l=1, lx2l=1, lx1u=1, lx2u=1)
    return {
        'md': float(md),
        'stage_a': float(stage_a),
        'stage_b': float(stage_b),
        'stage_c': float(stage_c),
        'method': 'legacy_wbtcdo_figure',
    }


def calculate_wbtcdo(mach: float, sweep_c4_deg: float,
                     thickness_ratio: float, aspect_ratio: float,
                     taper_ratio: float, configuration_type: float,
                     cd_mach06: float, cd_mach07: float,
                     cd_mach14: Optional[float] = None,
                     cd_mach11: Optional[float] = None
                     ) -> Dict[str, object]:
    """Translate WBTCDO: transonic wing-body-tail zero-lift drag.

    Args:
        mach: Free-stream Mach number.
        sweep_c4_deg: Quarter-chord sweep, degrees.
        thickness_ratio: ``TOC``.
        aspect_ratio: ``AR``.
        taper_ratio: ``TAPR``.
        configuration_type: ``ITYPE``.
        cd_mach06: ``CD6``, the drag level at Mach 0.6.
        cd_mach07: ``CD7``, at Mach 0.7.
        cd_mach14: ``CD14``, at Mach 1.4.
        cd_mach11: ``CD11``, at Mach 1.1.  When absent it is taken as
            ``2.1381 * CD14``.

    Returns:
        Dictionary with ``cdo``, ``md`` and the fairing coefficients, or
        ``cdo`` of None in the case the source returns without computing.

    Raises:
        ValueError: If a fairing denominator vanishes.

    Notes:
        The source returns early, leaving ``CDO`` at whatever it held, when
        the Mach 1.1 level is unknown, the Mach 1.4 level is also unknown,
        and the query Mach is above ``MD``.  That is reported as ``cdo`` of
        None rather than a stale value.

        The fairing reproduces its anchor exactly at each end: ``CDO``
        equals ``CD7`` at Mach 0.7 and ``CD11`` at Mach 1.1.

        The supersonic branch diverges as the query Mach approaches ``MD``
        from above.  Drag peaks near Mach 1.1 and falls by Mach 1.4, so
        ``DCD11`` is negative and the exponent comes out negative, which
        sends ``((M-MD)/(1.1-MD))**EXPN`` to infinity at the lower end of
        its own range.  That is a property of the source's formula, not of
        this translation: it is built to be read at Mach numbers
        meaningfully above the divergence Mach, and the subsonic branch has
        no such behaviour because its exponent is near unity.
    """
    figure = calculate_drag_divergence_mach(
        sweep_c4_deg, thickness_ratio, aspect_ratio, taper_ratio,
        configuration_type)
    md = figure['md']

    if cd_mach11 is None:
        if cd_mach14 is None and mach > md:
            return {
                'cdo': None,
                'md': md,
                'figure': figure,
                'reason': 'no supersonic drag level above the divergence Mach',
                'method': 'legacy_wbtcdo',
            }
        if cd_mach14 is None:
            raise ValueError(
                "WBTCDO needs CD14 to derive CD11 below the divergence Mach")
        cd_mach11 = _CD11_FROM_CD14 * cd_mach14

    supersonic = mach > md
    anchor_mach = _SUPERSONIC_ANCHOR if supersonic else _SUBSONIC_ANCHOR

    fairing_a0 = cd_mach07 + 0.002
    fairing_a1 = 0.10
    if supersonic:
        fairing_a2 = (cd_mach11 - fairing_a0 -
                      0.1 * (_SUPERSONIC_ANCHOR - md))
        drag_slope = ((cd_mach14 - cd_mach11) / 0.3
                      if cd_mach14 is not None else None)
        if drag_slope is None:
            raise ValueError("WBTCDO needs CD14 for the supersonic exponent")
        exponent_denominator = (cd_mach11 - fairing_a0 -
                                0.1 * (_SUPERSONIC_ANCHOR - md))
    else:
        fairing_a2 = -0.002 - 0.1 * (_SUBSONIC_ANCHOR - md)
        drag_slope = (cd_mach07 - cd_mach06) / 0.1
        exponent_denominator = (cd_mach07 - fairing_a0 -
                                0.1 * (_SUBSONIC_ANCHOR - md))

    if exponent_denominator == 0.0:
        raise ValueError("WBTCDO fairing exponent divides by zero")
    exponent = ((drag_slope - 0.1) * (anchor_mach - md) /
                exponent_denominator)

    if anchor_mach == md:
        raise ValueError(
            "WBTCDO fairing divides by (anchor - MD), which vanished")
    mach_offset = mach - md
    cdo = (fairing_a0 + fairing_a1 * mach_offset +
           fairing_a2 * (mach_offset / (anchor_mach - md)) ** exponent)

    return {
        'cdo': float(cdo),
        'md': md,
        'figure': figure,
        'a0': float(fairing_a0),
        'a1': fairing_a1,
        'a2': float(fairing_a2),
        'exponent': float(exponent),
        'anchor': anchor_mach,
        'supersonic': bool(supersonic),
        'cd_mach11': float(cd_mach11),
        'method': 'legacy_wbtcdo',
    }
