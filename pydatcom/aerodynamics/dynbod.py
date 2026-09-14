"""
DYNBOD: body-alone dynamic derivatives.

Two independent branches:

- Subsonic, transonic and supersonic: ``CLQ``, ``CMQ``, ``CLAD`` and
  ``CMAD`` from the body volume and centroid, driven by the static slopes
  the body routine already produced.
- Hypersonic: ``CNQ`` and ``CMQ`` from the closed-form segment expressions
  of Figures 7.2.1.1-9A and 7.2.1.2-12, summed over nose, afterbody and
  flare.  Those figures are algebraic in the source, so no table data is
  involved.

The source warns that the method is valid for a nose-cylinder only: boat
tails and flares are ignored in the first branch, and a boat tail disables
the flare term in the second.

Reference: datcom-legacy/datcom_2000/dynbod.f
"""

import numpy as np
from typing import Dict, Optional, Sequence
import logging

from pydatcom.utils.constants import PI
from pydatcom.utils.legacy_numeric import trapz

logger = logging.getLogger(__name__)


def calculate_dynbod_subsonic(x: Sequence[float], s: Sequence[float],
                              cla_body: float, cma_body: float,
                              body_length: float, body_area: float,
                              base_area: float, xcg: float,
                              sref: float, cbar: float) -> Dict[str, float]:
    """Translate DYNBOD's subsonic, transonic and supersonic branch.

    ``CLQ = 2*CLA*(1-h)`` and ``CLAD = 2*CLA*v`` with ``h = XCG/LB`` and
    ``v = VB/(SB*LB)``; the moment derivatives carry the volume centroid.

    Args:
        x: Body station coordinates.
        s: Cross-sectional area at each station.
        cla_body: Body lift-curve slope, the source's ``CLAB``.
        cma_body: Body moment slope, ``CMAB``.
        body_length: ``LB``, nose plus afterbody length.
        body_area: ``SB``, the body reference area.
        base_area: ``SBASE``, used for the reference-area correction.
        xcg: Moment reference station.
        sref: Aircraft reference area.
        cbar: Aircraft reference chord.

    Returns:
        Dictionary with ``clq``, ``cmq``, ``clad``, ``cmad``, the body
        volume and centroid, and the intermediate ratios.

    Raises:
        ValueError: If the geometry or references are unusable, or the
            configuration makes the source denominator vanish.

    Notes:
        ``CLAD`` is scaled by ``LB/CBARR`` to the first power while the
        other three take it squared.  That asymmetry is in the source and is
        dimensionally correct: the acceleration derivatives are normalised
        on a different power of the reference length.
    """
    x = np.asarray(x, dtype=float)
    s = np.asarray(s, dtype=float)
    if x.ndim != 1 or x.shape != s.shape or len(x) < 2:
        raise ValueError("DYNBOD requires at least two matching body stations")
    if min(body_length, body_area, sref, cbar) <= 0.0:
        raise ValueError("DYNBOD requires positive length, area and references")

    volume = float(trapz(s, x)[0])                       # VB
    if volume <= 0.0:
        raise ValueError("DYNBOD requires a positive body volume")
    first_moment = float(trapz(s * x, x)[0])
    centroid = first_moment / volume                     # XC

    cg_fraction = xcg / body_length                      # SAVE
    volume_ratio = volume / (body_area * body_length)    # TEMP
    denominator = 1.0 - cg_fraction - volume_ratio
    if abs(denominator) < 1e-12:
        raise ValueError(
            "DYNBOD's moment denominator (1 - XCG/LB - VB/(SB*LB)) vanishes "
            "for this configuration")

    centroid_arm = centroid / body_length - cg_fraction
    clq = 2.0 * cla_body * (1.0 - cg_fraction)
    cmq = (2.0 * cma_body *
           ((1.0 - cg_fraction)**2 - volume_ratio * centroid_arm) /
           denominator)
    clad = 2.0 * cla_body * volume_ratio
    cmad = 2.0 * cma_body * (volume_ratio * centroid_arm) / denominator

    # Reference correction.  CLAD takes the first power of LB/CBARR.
    area_ratio = base_area / sref
    length_ratio = body_length / cbar
    return {
        'clq': float(clq * area_ratio * length_ratio**2),
        'cmq': float(cmq * area_ratio * length_ratio**2),
        'clad': float(clad * area_ratio * length_ratio),
        'cmad': float(cmad * area_ratio * length_ratio**2),
        'volume': volume,
        'centroid': centroid,
        'cg_fraction': float(cg_fraction),
        'volume_ratio': float(volume_ratio),
        'method': 'legacy_dynbod_subsonic',
    }


def _segment_cnq(half_angle: float, taper: float, diameter: float,
                 scale: float) -> float:
    """Figure 7.2.1.1-9A for one conical segment.

    A segment of unit taper has no projected area change and contributes
    nothing, which is the source's ``LAM == 1`` branch.
    """
    if taper == 1.0:
        return 0.0
    tangent = np.tan(half_angle)
    if tangent == 0.0:
        raise ValueError(
            "Figure 7.2.1.1-9A divides by tan(theta); a cylindrical segment "
            "needs a unit taper ratio instead of a zero half angle")
    cosine = np.cos(half_angle)
    term = (0.66667 / tangent *
            (2.0 * (1.0 - taper**3) -
             3.0 * taper * cosine**2 * (1.0 - taper**2)))
    return float(term * scale * diameter**3)


def _segment_cmq(half_angle: float, taper: float, diameter: float,
                 scale: float) -> float:
    """Figure 7.2.1.2-12 for one conical segment."""
    if taper == 1.0:
        return 0.0
    sine = np.sin(half_angle)
    if sine == 0.0:
        raise ValueError(
            "Figure 7.2.1.2-12 divides by sin(theta)^2; a cylindrical "
            "segment needs a unit taper ratio instead of a zero half angle")
    cosine = np.cos(half_angle)
    a = 6.0 * taper**2 * (1.0 - taper**2)
    b = -8.0 * taper * (1.0 - taper**3)
    c = 3.0 * (1.0 - taper**4)
    term = (-(a * cosine**4 + b * cosine**2 + c) / (6.0 * sine**2) -
            taper**4 / 2.0)
    return float(term * scale * diameter**4)


def calculate_dynbod_hypersonic(theta_nose: float, theta_afterbody: float,
                                theta_flare: float,
                                d_nose: float, d1: float, d2: float,
                                nose_length: float, afterbody_length: float,
                                cna: Sequence[float], cma: Sequence[float],
                                xcg: float, sref: float,
                                cbar: float) -> Dict[str, float]:
    """Translate DYNBOD's hypersonic branch: CNQ and CMQ.

    Args:
        theta_nose: Nose half angle, radians, the source's ``THETAN``.
        theta_afterbody: Afterbody half angle, ``THETAA``.
        theta_flare: Flare half angle, ``THETAF``.  Zero or negative
            disables the flare term, as a boat tail would.
        d_nose: Nose base diameter, ``DN``.
        d1: Afterbody base diameter.
        d2: Flare base diameter.
        nose_length: ``LN``.
        afterbody_length: ``LA``.
        cna: Static normal-force slopes per segment: nose, afterbody, flare.
        cma: Static moment slopes per segment, same order.
        xcg: Moment reference station.
        sref: Reference area.
        cbar: Reference chord.

    Returns:
        Dictionary with ``cnq`` and ``cmq`` and their per-segment parts.

    Raises:
        ValueError: If the references are nonpositive or a segment's angle
            and taper are inconsistent.
    """
    if min(sref, cbar) <= 0.0:
        raise ValueError("DYNBOD requires positive SREF and CBARR")
    if len(cna) != 3 or len(cma) != 3:
        raise ValueError("DYNBOD needs three segment slopes: nose, "
                         "afterbody and flare")

    boat_tailed = theta_flare <= 0.0
    if boat_tailed:
        logger.warning(
            "DYNBOD hypersonic method is not valid for bodies with "
            "boattails; boattail effects ignored")
        theta_flare = 0.0
    if d2 < 1.0e-10:
        d2 = 0.01

    # LAMN is fixed at zero by the source: the nose closes to a point.
    taper_nose = 0.0
    taper_afterbody = d_nose / d1 if d1 else 1.0
    taper_flare = 1.0 if boat_tailed else (d1 / d2 if d2 else 1.0)

    scale = PI / (4.0 * sref * cbar)
    cnq_nose = _segment_cnq(theta_nose, taper_nose, d_nose, scale)
    cnq_after = _segment_cnq(theta_afterbody, taper_afterbody, d1, scale)
    cnq_flare = _segment_cnq(theta_flare, taper_flare, d2, scale)

    arm_nose = xcg / cbar                                 # NN
    arm_after = (xcg - nose_length) / cbar                # NA
    arm_flare = (xcg - afterbody_length) / cbar           # NF

    cnq_segments = (
        cnq_nose - 2.0 * arm_nose * cna[0],
        cnq_after - 2.0 * arm_after * cna[1],
        cnq_flare - 2.0 * arm_flare * cna[2],
    )

    moment_scale = scale / cbar
    cmq_nose = _segment_cmq(theta_nose, taper_nose, d_nose, moment_scale)
    cmq_after = _segment_cmq(theta_afterbody, taper_afterbody, d1,
                             moment_scale)
    cmq_flare = _segment_cmq(theta_flare, taper_flare, d2, moment_scale)

    cmq_segments = (
        cmq_nose - 2.0 * arm_nose * cma[0] + arm_nose * cnq_nose -
        2.0 * arm_nose**2 * cna[0],
        cmq_after - 2.0 * arm_after * cma[1] + arm_after * cnq_after -
        2.0 * arm_after**2 * cna[1],
        cmq_flare - 2.0 * arm_flare * cma[2] + arm_flare * cnq_flare -
        2.0 * arm_flare**2 * cna[2],
    )

    return {
        'cnq': float(sum(cnq_segments)),
        'cmq': float(sum(cmq_segments)),
        'cnq_segments': tuple(float(v) for v in cnq_segments),
        'cmq_segments': tuple(float(v) for v in cmq_segments),
        'taper_ratios': (taper_nose, float(taper_afterbody),
                         float(taper_flare)),
        'boat_tailed': boat_tailed,
        'method': 'legacy_dynbod_hypersonic',
    }
