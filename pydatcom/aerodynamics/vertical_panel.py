"""
Zero-lift drag of a vertical tail or ventral fin.

Translates four source routines with two functions:

- ``calculate_vertical_panel_cdo`` covers the supersonic ``VRTCDO`` and
  ``VFCDO``: skin friction plus wave drag.
- ``calculate_vertical_panel_drag`` covers the subsonic ``VTDRAG`` and
  ``VFDRAG``: skin friction times a form factor and the Figure 4.1.5.1-28B
  lifting-surface correction.

Within each pair the two source routines are identical apart from the COMMON
block offsets that select which surface's geometry they read -- the ventral
fin reads a later slice of ``/VTI/`` and ``/VTDATA/`` -- so one translation
serves both, with the surface passed in rather than selected by offset.
That equivalence was established by diffing the files, not assumed from
their names.

Both supersonic and subsonic paths cap the panel MAC Reynolds number with
the Figure 4.1.5.1-27 roughness cutoff before calling ``FIG26``, and both
handle straight tapered and cranked planforms, the latter computing inboard
and outboard panels separately.

Reference: datcom-legacy/datcom_2000/vrtcdo.f, vfcdo.f, vtdrag.f, vfdrag.f
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
        average_chord = (root_chord + tip_chord) / 2.0
        area = average_chord * semispan
        le_radius = leading_edge_radius * average_chord
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
    wave_drag_denominator = beta if sonic_leading_edge else tangent

    if ksharp is not None:
        cd_leading_edge = 0.0
        cd_wave = (ksharp * thickness_ratio ** 2 * area / sref /
                   wave_drag_denominator)
        edge = 'sharp'
    else:
        blunt = (1.28 * mach ** 3 * cosine ** 6 /
                 (1.0 + mach ** 3 * cosine ** 3))
        cd_leading_edge = blunt * 2.0 * le_radius * semispan / (sref * cosine)
        cd_wave = (cd_leading_edge +
                   16.0 * thickness_ratio ** 2 * area / (3.0 * sref) /
                   wave_drag_denominator)
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


# Figure 4.1.5.1-28B: lifting-surface correction (R)LS.  vtdrag.f DATA
# X128B / X228B / Y28B.  The Mach axis is descending in the source and is
# kept that way; TLINEX handles either direction.
_FIG_41510_28B_MACH = [0.9, .80, .60, .25]
_FIG_41510_28B_COS = [0.5, .55, .60, .65, .70, .75, .80, .85, .90, .95, 1.0]
_FIG_41510_28B = np.array([
    [1.1, 1.13, 1.17, 1.20, 1.24, 1.27, 1.3, 1.33, 1.34, 1.35, 1.36],
    [1.0, 1.04, 1.08, 1.11, 1.15, 1.18, 1.21, 1.23, 1.25, 1.25, 1.26],
    [0.88, .92, .96, 1.0, 1.04, 1.08, 1.11, 1.13, 1.14, 1.14, 1.15],
    [0.81, .85, .89, .925, .96, 1.0, 1.03, 1.05, 1.06, 1.06, 1.07],
]).T

# Figure 4.1.5.1-28B dashed: the inboard-panel variant.  The source writes
# two repeated entries of the last column as "2*1.06".
_FIG_41510_28BD_MACH = [0.9, .80, .60, .25]
_FIG_41510_28BD_COS = [0.45, .65, .70, .75, .80, .85, .90, .95, 1.0]
_FIG_41510_28BD = np.array([
    [1.20, 1.20, 1.24, 1.27, 1.30, 1.33, 1.34, 1.35, 1.36],
    [1.11, 1.11, 1.15, 1.18, 1.21, 1.23, 1.25, 1.25, 1.26],
    [1.0, 1.0, 1.04, 1.08, 1.11, 1.13, 1.14, 1.14, 1.15],
    [0.925, 0.925, .96, 1.0, 1.03, 1.05, 1.06, 1.06, 1.07],
]).T

# The form-factor coefficient switches on maximum-thickness chord station.
_THICKNESS_STATION_SPLIT = 0.30


def _lifting_surface_factor(mach: float, cos_sweep_tmax: float,
                            dashed: bool) -> float:
    """Figure 4.1.5.1-28B, solid or dashed variant."""
    from pydatcom.utils.legacy_tables import tlinex
    if dashed:
        return float(tlinex(_FIG_41510_28BD_MACH, _FIG_41510_28BD_COS,
                            _FIG_41510_28BD, mach, cos_sweep_tmax,
                            0, 0, 0, 0))
    return float(tlinex(_FIG_41510_28B_MACH, _FIG_41510_28B_COS,
                        _FIG_41510_28B, mach, cos_sweep_tmax, 0, 2, 0, 2))


def _form_factor(thickness_ratio: float, thickness_station: float) -> float:
    """``1 + L*(t/c) + 100*(t/c)^4`` with the source's L switch."""
    coefficient = (1.20 if thickness_station >= _THICKNESS_STATION_SPLIT
                   else 2.00)
    return 1.0 + coefficient * thickness_ratio + 100.0 * thickness_ratio**4


def calculate_vertical_panel_drag(
        mach: float, reynolds_per_length: float, sref: float,
        mac: float, area: float, cos_sweep_tmax: float,
        thickness_ratio: float, thickness_station: float,
        roughness: float = 1.6e-4,
        straight: bool = True,
        mac_inboard: Optional[float] = None,
        area_inboard: Optional[float] = None,
        cos_sweep_tmax_inboard: Optional[float] = None,
        mac_outboard: Optional[float] = None,
        area_outboard: Optional[float] = None,
        cos_sweep_tmax_outboard: Optional[float] = None,
        thickness_ratio_outboard: Optional[float] = None,
        thickness_station_outboard: Optional[float] = None
) -> Dict[str, float]:
    """Translate VTDRAG/VFDRAG: subsonic vertical panel zero-lift drag.

    ``CDO = CF * (1 + L*(t/c) + 100*(t/c)^4) * (R)LS * 2*S/SREF``

    As with VRTCDO and VFCDO, vtdrag.f and vfdrag.f are identical apart from
    COMMON offsets and whitespace, so this covers both surfaces.

    A cranked planform computes inboard and outboard panels separately and
    sums them.  The inboard panel uses the dashed variant of Figure
    4.1.5.1-28B; the outboard panel and any straight panel use the solid one.

    Args:
        mach: Free-stream Mach number, below one.
        reynolds_per_length: Reynolds number per unit length.
        sref: Aircraft reference area.
        mac: Panel MAC, used when ``straight``.
        area: Panel exposed area, used when ``straight``.
        cos_sweep_tmax: Cosine of the sweep at maximum thickness.
        thickness_ratio: ``t/c``.
        thickness_station: Chord station of maximum thickness, ``XOVC``.
        roughness: Surface roughness height.
        straight: Whether the planform is straight tapered.
        mac_inboard, area_inboard, cos_sweep_tmax_inboard: Inboard geometry.
        mac_outboard, area_outboard, cos_sweep_tmax_outboard,
        thickness_ratio_outboard, thickness_station_outboard: Outboard
            geometry.  All are required when ``straight`` is false.

    Returns:
        Dictionary with ``cdo`` and, for a cranked panel, its two parts.

    Raises:
        ValueError: For supersonic Mach, nonpositive references, or missing
            panel geometry.

    Notes:
        The source writes the L switch two different ways -- the straight
        branch defaults to 2.0 and drops to 1.2 at or above a 0.30 station,
        the inboard branch defaults to 1.2 and rises to 2.0 below it -- which
        are logically identical.  One form is used here.
    """
    if mach >= 1.0:
        raise ValueError("VTDRAG/VFDRAG are subsonic; Mach must be < 1")
    if sref <= 0.0:
        raise ValueError("VTDRAG/VFDRAG require a positive SREF")
    if roughness <= 0.0:
        raise ValueError("VTDRAG/VFDRAG require a positive roughness height")

    def panel(panel_mac, panel_area, cosine, tc, station, dashed):
        if min(panel_mac, panel_area) <= 0.0:
            raise ValueError("panel MAC and area must be positive")
        friction = _friction_coefficient(panel_mac, reynolds_per_length,
                                         mach, roughness)
        rls = _lifting_surface_factor(mach, cosine, dashed)
        cdo = (friction["cf"] * _form_factor(tc, station) * rls *
               2.0 * panel_area / sref)
        return cdo, friction["cf"], rls

    if straight:
        cdo, cf, rls = panel(mac, area, cos_sweep_tmax, thickness_ratio,
                             thickness_station, dashed=False)
        return {
            "cdo": float(cdo), "cf": cf, "lifting_surface_factor": rls,
            "form_factor": _form_factor(thickness_ratio, thickness_station),
            "method": "legacy_vtdrag",
        }

    required = (mac_inboard, area_inboard, cos_sweep_tmax_inboard,
                mac_outboard, area_outboard, cos_sweep_tmax_outboard,
                thickness_ratio_outboard, thickness_station_outboard)
    if any(value is None for value in required):
        raise ValueError(
            "a cranked vertical panel needs inboard and outboard MAC, area, "
            "max-thickness sweep cosine, thickness ratio and station")

    inboard_cdo, inboard_cf, inboard_rls = panel(
        mac_inboard, area_inboard, cos_sweep_tmax_inboard,
        thickness_ratio, thickness_station, dashed=True)
    outboard_cdo, outboard_cf, outboard_rls = panel(
        mac_outboard, area_outboard, cos_sweep_tmax_outboard,
        thickness_ratio_outboard, thickness_station_outboard, dashed=False)
    return {
        "cdo": float(inboard_cdo + outboard_cdo),
        "cdo_inboard": float(inboard_cdo),
        "cdo_outboard": float(outboard_cdo),
        "cf_inboard": inboard_cf, "cf_outboard": outboard_cf,
        "lifting_surface_factor_inboard": inboard_rls,
        "lifting_surface_factor_outboard": outboard_rls,
        "method": "legacy_vtdrag",
    }


def m08o10_panel_block(cd0: float, alpha_count: int) -> Dict[str, np.ndarray]:
    """Translate M08O10's output setup for a vertical panel (VT or VF).

    After VTDRAG or VFDRAG, the panel's only first-angle data are its
    zero-lift drag ``DVT(20)`` or ``DVF(20)``, with zero lift, moment,
    normal force and slopes; every later angle is marked ``-UNUSED``.
    """
    unused = -1.0e-30
    block = {key: np.full(alpha_count, unused)
             for key in ('cd', 'cl', 'cm', 'cn', 'ca', 'cla', 'cma')}
    block['cd'][0] = cd0
    for key in ('cl', 'cm', 'cn', 'cla', 'cma'):
        block[key][0] = 0.0
    return block
