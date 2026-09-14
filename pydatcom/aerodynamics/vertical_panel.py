"""
Supersonic zero-lift drag of a vertical tail or ventral fin.

Translates ``VRTCDO`` and ``VFCDO``.  The two source routines are identical
line for line apart from the COMMON block offsets that select which
surface's geometry they read -- the ventral fin reads a later slice of
``/VTI/`` and ``/VTDATA/`` -- so one translation serves both, with the
surface passed in rather than selected by offset.

The buildup is skin friction plus wave drag:

- Skin friction uses the panel MAC Reynolds number capped by the Figure
  4.1.5.1-27 roughness cutoff, then ``FIG26``.  A straight tapered panel
  uses one MAC; any other planform uses inboard and outboard panels
  separately and area-weights the result.
- Wave drag takes the sharp leading-edge form when ``KSHARP`` is supplied
  and the round leading-edge form otherwise, each selecting a ``beta`` or
  ``tan(sweep)`` denominator on the sonic leading-edge condition.

Reference: datcom-legacy/datcom_2000/vrtcdo.f, vfcdo.f
"""

import numpy as np
from typing import Dict, Optional
import logging

from pydatcom.utils.legacy_numeric import tbfunx
from pydatcom.utils.table_lookup import fig26

logger = logging.getLogger(__name__)

# Figure 4.1.5.1-27: roughness-limited Reynolds number intercept vs Mach.
_FIG_41510_27_MACH = np.array([0.0, 1.0, 2.0, 3.0])
_FIG_41510_27_CEPT = np.array([1.57780, 1.67221, 1.98509, 2.28874])

# The source caps the friction Mach lookup at 3.
_MACH_CAP = 3.0


def _friction_coefficient(mac: float, reynolds_per_length: float,
                          mach_lookup: float, roughness: float) -> Dict:
    """One panel's skin-friction coefficient with the roughness cutoff."""
    reynolds = mac * reynolds_per_length
    cutoff = None
    if roughness > 0.0:
        cutoff = (12.0 * mac / roughness)**1.0482 * 10.0**_cept(mach_lookup)
        reynolds = min(reynolds, cutoff)
    return {
        'cf': float(fig26(reynolds, mach_lookup)),
        'reynolds_used': float(reynolds),
        'roughness_cutoff_reynolds': cutoff,
    }


def _cept(mach_lookup: float) -> float:
    value, _ = tbfunx(_FIG_41510_27_MACH, _FIG_41510_27_CEPT,
                      mach_lookup, lower=0, upper=0)
    return float(value)


def calculate_vertical_panel_cdo(
        mach: float, reynolds_per_length: float, sref: float,
        mac_inboard: float, area_inboard: float,
        root_chord: float, tip_chord: float, semispan: float,
        tan_le: float, cos_le: float,
        thickness_ratio: float, leading_edge_radius: float,
        roughness: float = 1.6e-4,
        straight: bool = True,
        mac_outboard: Optional[float] = None,
        area_outboard: Optional[float] = None,
        break_chord: Optional[float] = None,
        span_outboard: Optional[float] = None,
        tan_le_outboard: Optional[float] = None,
        cos_le_outboard: Optional[float] = None,
        tan_te_outboard: Optional[float] = None,
        leading_edge_radius_outboard: Optional[float] = None,
        ksharp: Optional[float] = None) -> Dict[str, float]:
    """Translate VRTCDO/VFCDO: supersonic zero-lift drag of a vertical panel.

    Args:
        mach: Free-stream Mach number, above one.
        reynolds_per_length: Reynolds number per unit length, ``RNFS``.
        sref: Aircraft reference area.
        mac_inboard: Inboard panel MAC, ``CBARI``.
        area_inboard: Inboard exposed area, ``SISTAR``; for a straight panel
            this is the whole panel area ``SRSTAR``.
        root_chord: Panel root chord, ``CR``.
        tip_chord: Panel tip chord, ``CT``.
        semispan: Panel semispan, ``SSPN`` for a straight panel and
            ``SPANS`` for a cranked one.
        tan_le: Tangent of the inboard leading-edge sweep.
        cos_le: Cosine of the inboard leading-edge sweep.
        thickness_ratio: Effective thickness ratio, ``TCEFF``.
        leading_edge_radius: Inboard leading-edge radius parameter, ``LERI``.
        roughness: Surface roughness height, ``RUFF``.
        straight: Whether the planform is straight tapered.
        mac_outboard, area_outboard, break_chord, span_outboard,
        tan_le_outboard, cos_le_outboard, tan_te_outboard,
        leading_edge_radius_outboard: Outboard panel geometry, required when
            ``straight`` is false.
        ksharp: Sharp leading-edge wave-drag factor.  When omitted the round
            leading-edge form is used, as the source's UNUSED test selects.

    Returns:
        Dictionary with ``cdo`` and its friction and wave components.

    Raises:
        ValueError: For subsonic Mach, nonpositive references, or missing
            outboard geometry on a non-straight planform.

    Notes:
        The source computes its Mach lookup ``RACH`` only inside the
        roughness branch, then calls ``FIG26`` with it unconditionally.  With
        ``RUFF`` exactly zero that call reads an uninitialised variable.  The
        translation always forms the capped lookup, which is what every
        nonzero-roughness run does and the only defensible reading.
    """
    if mach <= 1.0:
        raise ValueError("VRTCDO/VFCDO require supersonic Mach")
    if min(sref, mac_inboard, area_inboard) <= 0.0:
        raise ValueError("VRTCDO/VFCDO require positive SREF, MAC and area")
    if roughness < 0.0:
        raise ValueError("roughness height cannot be negative")

    beta = np.sqrt(mach**2 - 1.0)
    mach_lookup = min(float(mach), _MACH_CAP)

    inboard = _friction_coefficient(mac_inboard, reynolds_per_length,
                                    mach_lookup, roughness)
    if straight:
        cd_friction = inboard['cf'] * area_inboard / sref * 2.0
        outboard = None
    else:
        required = (mac_outboard, area_outboard, break_chord, span_outboard,
                    tan_le_outboard, cos_le_outboard, tan_te_outboard,
                    leading_edge_radius_outboard)
        if any(value is None for value in required):
            raise ValueError(
                "a non-straight vertical panel needs the full outboard "
                "geometry: MAC, area, break chord, outboard span, leading "
                "and trailing edge sweeps, and leading-edge radius")
        outboard = _friction_coefficient(mac_outboard, reynolds_per_length,
                                         mach_lookup, roughness)
        cd_friction = ((inboard['cf'] * area_inboard +
                        outboard['cf'] * area_outboard) / sref * 2.0)

    # Wave-drag geometry.  A straight panel uses its own trapezoid; a
    # cranked one rebuilds the outboard trapezoid from the break chord.
    if straight:
        tangent = tan_le
        cosine = cos_le
        area = (root_chord + tip_chord) / 2.0 * semispan
        le_radius = leading_edge_radius * ((root_chord + tip_chord) / 2.0)
    else:
        tangent = tan_le_outboard
        cosine = cos_le_outboard
        span_inboard = semispan - span_outboard
        root_equivalent = break_chord + span_inboard * (tan_le_outboard -
                                                        tan_te_outboard)
        area = (root_equivalent + tip_chord) * semispan * 0.5
        le_radius = leading_edge_radius_outboard * (
            (root_chord + tip_chord + 2.0 * break_chord) / 4.0)

    if tangent == 0.0:
        raise ValueError(
            "wave drag divides by tan(leading-edge sweep); an unswept panel "
            "needs the source's UNUSED floor rather than exactly zero")
    sonic_leading_edge = beta / tangent >= 1.0
    denominator = beta if sonic_leading_edge else tangent

    if ksharp is not None:
        cd_leading_edge = 0.0
        cd_wave = ksharp * thickness_ratio**2 * area / sref / denominator
        edge = 'sharp'
    else:
        blunt = (1.28 * mach**3 * cosine**6 /
                 (1.0 + mach**3 * cosine**3))
        cd_leading_edge = blunt * 2.0 * le_radius * semispan / (sref * cosine)
        cd_wave = (cd_leading_edge +
                   16.0 * thickness_ratio**2 * area / (3.0 * sref) /
                   denominator)
        edge = 'round'

    return {
        'cdo': float(cd_friction + cd_wave),
        'cd_friction': float(cd_friction),
        'cd_wave': float(cd_wave),
        'cd_leading_edge': float(cd_leading_edge),
        'cf_inboard': inboard['cf'],
        'cf_outboard': outboard['cf'] if outboard else None,
        'reynolds_used': inboard['reynolds_used'],
        'roughness_cutoff_reynolds': inboard['roughness_cutoff_reynolds'],
        'sonic_leading_edge': bool(sonic_leading_edge),
        'leading_edge': edge,
        'wave_drag_area': float(area),
        'method': 'legacy_vrtcdo',
    }
