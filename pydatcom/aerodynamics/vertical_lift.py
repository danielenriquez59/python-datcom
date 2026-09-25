"""
VTLIFT and VFLIFT: supersonic normal-force slope of a vertical panel.

The two routines are identical apart from whitespace and the COMMON offsets
that select the vertical tail or the ventral fin, so one translation serves
both.  The figure layer lives in :mod:`vertical_lift_figures`; this module
is the branch logic that selects among those figures.

There are two paths.  A straight tapered panel goes through a single
theoretical slope and one correction factor.  Anything else is decomposed
into three pieces that are summed and then scaled:

1. the *basic wing*, the outboard panel extended inboard to the centreline
2. the *glove*, the triangle the inboard leading edge adds ahead of it
3. the *extension*, present only when the inboard and outboard trailing
   edges differ by at least four degrees

The basic-wing and glove evaluations share a single block of source code
that is entered twice, so the same figure sequence runs with different
geometry each time.

Three source defects are reproduced rather than corrected; see
:func:`calculate_vtlift`.

Reference: datcom-legacy/datcom_2000/vtlift.f, vflift.f
"""

import math
import numpy as np
from typing import Dict, Optional
import logging

from pydatcom.utils.constants import RAD
from pydatcom.aerodynamics.vertical_lift_figures import (
    fig4132_56a, fig4132_56g, fig4132_60a, fig4132_60b, fig4132_61,
    fig4132_62, fig4132_63)

logger = logging.getLogger(__name__)

# The source's WTYPE(1), the straight tapered planform selector.
STRAIGHT_TAPERED = 1.0

# The source replaces a zero leading-edge tangent with this before dividing.
_ZERO_TANGENT = 0.00001

# Trailing edges closer than this contribute no extension component.
_EXTENSION_SWEEP_THRESHOLD = 4.0

# LGB(1) as it stands at the second figure 4.1.3.2-56A call; see the note in
# calculate_vtlift.  The grid itself holds 23 entries.
_LEAKED_FIRST_LENGTH = 12
_FULL_FIRST_LENGTH = 23


def calculate_vtlift(mach: float,
                     panel: Dict[str, float],
                     sweep: Dict[str, float],
                     sref: float,
                     spans: Optional[float] = None,
                     planform_type: float = STRAIGHT_TAPERED,
                     ksharp: Optional[float] = None) -> Dict[str, object]:
    """Translate VTLIFT/VFLIFT: supersonic vertical-panel normal-force slope.

    Args:
        mach: Free-stream Mach number, the source's ``FLC(I+2)``.  Must be
            supersonic.
        panel: Panel geometry as ``{'area', 'aspect_ratio', 'taper_ratio',
            'taper_ratio_exposed', 'chrdtp', 'chrdbp', 'sspnop',
            'delta_y'}``.  These are ``AVT(3)``, ``AVT(7)``, ``AVT(25)``,
            ``AVT(27)``, ``VTIN(1)``, ``VTIN(5)``, ``VTIN(2)`` and
            ``VTIN(17)``.
        sweep: Sweep data as ``{'sweple_deg', 'cosle', 'tanle', 'tanleo',
            'tantei', 'tanteo', 'swtei_deg', 'swteo_deg'}``, the source's
            ``AVT(58)``, ``AVT(61)``, ``AVT(62)``, ``AVT(86)``, ``AVT(80)``,
            ``AVT(104)``, ``AVT(76)`` and ``AVT(100)``.  ``tanle`` is the
            inboard leading-edge tangent; the source's ``TANLE`` and
            ``TANLEI`` are EQUIVALENCEd to the same slot.
        sref: Reference area, the source's ``SW``.
        spans: The source's ``SPANS``.  Required on the non-straight path;
            see the note below about which block it reads.
        planform_type: ``VTIN(15)``.  Anything but straight tapered takes
            the three-component path.
        ksharp: ``VTIN(71)``.  ``None`` is the source's UNUSED sentinel and
            selects the round leading-edge curve of figure 4.1.3.2-62.

    Returns:
        Dictionary with ``cna``, the intermediate quantities of whichever
        path ran, and the flags described below.

    Raises:
        ValueError: If the Mach number is not supersonic, if the reference
            area is nonpositive, or if ``spans`` is missing on the
            non-straight path.

    Notes:
        Three source defects are preserved, because correcting any of them
        would change results:

        ``SPANS`` is EQUIVALENCEd to ``WINGIN(3)`` -- the *wing's* exposed
        semispan -- inside a routine whose every other planform quantity
        comes from ``VTIN``, the vertical panel's own block.  It is used as
        the panel semispan throughout the non-straight path, setting the
        basic-wing area, aspect ratio and taper.  The parallel slot is
        ``VTIN(3)``.  Nothing in DATCOM stages vertical-panel geometry into
        ``/WINGI/`` before the call, so the read is literal.  The caller
        supplies the value and ``spans_reads_wing_block`` records it.

        ``LGB(1)`` is set to 23 for the first figure 4.1.3.2-56A call, then
        overwritten with 12 by the figure 4.1.3.2-62 lookup that follows.
        The second 56A call, in the ``TANLEO >= TANLEI`` branch, reuses
        ``LGB`` without restoring it, so that lookup sees only the first 12
        of the 23 first-variable grid entries -- capping its axis at 1.0
        instead of 30.  The glove pass re-enters at the label that resets
        ``LGB(1)``, so only this one call is affected.
        ``figure_56a_second_call_uses_truncated_grid`` flags it.

        ``CNT2`` divides by ``BETA`` for a subsonic leading edge and by
        ``TANLEO`` for a supersonic one.  Every other use of figure
        4.1.3.2-56A in the routine does the opposite, ``CNTHRY=BCNA/TA``
        with ``IF(SUPLE)CNTHRY=BCNA/BETA``.  The two conventions are exactly
        swapped.  ``cnt2_divisor_inverted`` flags it and ``cnt2_consistent``
        reports the value the prevailing convention would have given.
    """
    if not np.isfinite(mach) or mach <= 1.0:
        raise ValueError(
            f"VTLIFT is the supersonic panel routine; Mach {mach} is not "
            "supersonic")
    if sref <= 0.0:
        raise ValueError("VTLIFT requires a positive reference area")

    beta = math.sqrt(mach**2 - 1.0)
    # The source doubles the panel aspect ratio on entry and halves it again
    # before returning: AVT(7) holds one side, the figures want the
    # image-plane equivalent.
    ar = 2.0 * float(panel['aspect_ratio'])

    tanle = float(sweep['tanle']) or _ZERO_TANGENT
    tanleo = float(sweep['tanleo']) or _ZERO_TANGENT
    tanlei = tanle          # TANLE and TANLEI are the same COMMON slot.

    if float(planform_type) == STRAIGHT_TAPERED:
        return _straight_tapered(beta, ar, tanle, panel, sweep, sref)
    if spans is None:
        raise ValueError(
            "VTLIFT's non-straight path needs SPANS; the source reads it "
            "from WINGIN(3), the wing's exposed semispan")
    return _three_component(beta, ar, tanlei, tanleo, float(spans), panel,
                            sweep, sref, ksharp)


def _straight_tapered(beta: float, ar: float, tanle: float,
                      panel: Dict[str, float], sweep: Dict[str, float],
                      sref: float) -> Dict[str, object]:
    """Labels 1000-1090: the straight tapered panel."""
    bovert = beta / tanle
    cosle = float(sweep['cosle'])
    supersonic_le = bovert > 1.0

    # Figure 4.1.3.2-60A or -60B, the correction to theoretical slope.
    if not supersonic_le:
        deltyt = float(panel['delta_y']) / cosle
        cncnt = fig4132_60a(bovert, deltyt)
        correction_figure = '4.1.3.2-60A'
        second_argument = deltyt
    else:
        deltdt = float(math.atan(float(panel['delta_y']) /
                                 (5.85 * cosle)) * RAD)
        cncnt = fig4132_60b(1.0 / bovert, deltdt)
        correction_figure = '4.1.3.2-60B'
        second_argument = deltdt

    rectangular = (float(panel['taper_ratio']) == 1.0 and
                   float(sweep['sweple_deg']) == 0.0)
    if rectangular:
        if ar * beta <= 1.0:
            # Figure 4.1.3.2-56G.
            cnaa = fig4132_56g(ar * beta)
            cnthry = cnaa * ar
            bcna = None
            slope_figure = '4.1.3.2-56G'
        else:
            bcna = 4.0 - 2.0 * (1.0 / (ar * beta))
            cnthry = bcna / beta
            slope_figure = 'closed form'
    else:
        bcna = fig4132_56a(bovert, ar * tanle,
                           float(panel['taper_ratio_exposed']))
        cnthry = bcna / beta if supersonic_le else bcna / tanle
        slope_figure = '4.1.3.2-56A'

    cna = cnthry * cncnt * float(panel['area']) / (sref * RAD)
    return {
        'cna': float(cna),
        'path': 'straight_tapered',
        'beta': float(beta),
        'bovert': float(bovert),
        'supersonic_leading_edge': bool(supersonic_le),
        'rectangular': bool(rectangular),
        'cncnt': float(cncnt),
        'cnthry': float(cnthry),
        'bcna': None if bcna is None else float(bcna),
        'correction_figure': correction_figure,
        'slope_figure': slope_figure,
        'correction_second_argument': float(second_argument),
        'spans_reads_wing_block': False,
        'figure_56a_second_call_uses_truncated_grid': False,
        'cnt2_divisor_inverted': False,
        'method': 'legacy_vtlift',
    }


def _component(bovert: float, beta: float, tangent: float,
               ar_tangent: float, taper: float, round_le: bool,
               first_length: int = _FULL_FIRST_LENGTH):
    """Label 1110-1130: one pass of the shared figure sequence.

    The source enters this block twice, once for the basic wing and once
    for the glove, with ``TA``, ``BOVERT`` and ``SUPLE`` reset between.
    """
    supersonic_le = bovert > 1.0
    bcna = fig4132_56a(bovert, ar_tangent, taper, first_length=first_length)
    cnthry = bcna / beta if supersonic_le else bcna / tangent
    cle = fig4132_62(bovert, sharp=not round_le)
    return bcna, cnthry, cle, supersonic_le


def _three_component(beta: float, ar: float, tanlei: float, tanleo: float,
                     spans: float, panel: Dict[str, float],
                     sweep: Dict[str, float], sref: float,
                     ksharp: Optional[float]) -> Dict[str, object]:
    """Labels 1100-1180: basic wing, glove and extension."""
    chrdtp = float(panel['chrdtp'])
    chrdbp = float(panel['chrdbp'])
    spanin = spans - float(panel['sspnop'])
    tanteo = float(sweep['tanteo'])

    crbw = chrdbp + spanin * (tanleo - tanteo)
    sbw = (crbw + chrdtp) * spans
    if crbw == 0.0 or sbw == 0.0:
        raise ValueError(
            "VTLIFT's basic-wing root chord and area must be nonzero")
    arbw = 4.0 * spans**2 / sbw
    tapbw = chrdtp / crbw
    round_le = ksharp is None

    # ---- basic wing ----------------------------------------------------
    bovert = beta / tanleo
    bcna, cnthry, cle, supersonic_le = _component(
        bovert, beta, tanleo, arbw * tanleo, tapbw, round_le)
    clebw = cle

    truncated = False
    inverted = False
    cnt2 = None
    cnt2_consistent = None
    if tanleo < tanlei:
        cnabw = cnthry * sbw / sref * clebw
    else:
        # The subtracted delta-wing component.  S2 is a triangle of semispan
        # SPANIN whose leading edge is TANLEO, so A2*TANLEO is identically
        # 4.0 -- it looks variable but never is.
        delta_wing_area = spanin ** 2 * tanleo
        delta_wing_aspect = 4.0 * spanin ** 2 / delta_wing_area
        # LGB(1) is still 12 here, left over from the figure 4.1.3.2-62
        # lookup above.  See the note in calculate_vtlift.
        bcna2 = fig4132_56a(
            bovert, delta_wing_aspect * tanleo, 0.0,
            first_length=_LEAKED_FIRST_LENGTH,
        )
        truncated = True
        # The source's divisor convention here is the inverse of CNTHRY's.
        cnt2 = bcna2 / tanleo if bovert > 1.0 else bcna2 / beta
        cnt2_consistent = bcna2 / beta if bovert > 1.0 else bcna2 / tanleo
        inverted = True
        cnabw = ((cnthry * sbw / sref - cnt2 * delta_wing_area / sref) *
                 clebw)

    # ---- glove ----------------------------------------------------------
    crglv = tanlei * spanin
    sglv = crglv * spanin
    if sglv == 0.0:
        raise ValueError("VTLIFT's glove area must be nonzero")
    arglv = 4.0 * spanin**2 / sglv
    bovert_glove = beta / tanlei
    _, cnthry_glove, cleglv, _ = _component(
        bovert_glove, beta, tanlei, arglv * tanlei, 0.0, round_le)
    cnaglv = cnthry_glove * cleglv * sglv / sref

    # ---- extension ------------------------------------------------------
    swtei = float(sweep['swtei_deg'])
    swteo = float(sweep['swteo_deg'])
    if abs(swtei - swteo) < _EXTENSION_SWEEP_THRESHOLD:
        cnae = 0.0
        cn1 = cn2 = None
    else:
        be = 2.0 * spanin
        cn1 = fig4132_63(beta / tanlei, float(sweep['tantei']) / tanlei)
        cn2 = fig4132_63(beta / tanlei, tanteo / tanlei)
        cnae = (cn1 - cn2) * be**2 / sref

    # ---- total ----------------------------------------------------------
    rkl = fig4132_61(cnaglv / (cleglv * beta))
    cna = rkl * (cnabw + cnaglv + cnae) / RAD
    return {
        'cna': float(cna),
        'path': 'three_component',
        'beta': float(beta),
        'bovert': float(bovert),
        'supersonic_leading_edge': bool(supersonic_le),
        'spanin': float(spanin),
        'crbw': float(crbw), 'sbw': float(sbw), 'arbw': float(arbw),
        'tapbw': float(tapbw),
        'cnabw': float(cnabw),
        'cnaglv': float(cnaglv),
        'cnae': float(cnae),
        'rkl': float(rkl),
        'clebw': float(clebw), 'cleglv': float(cleglv),
        'sglv': float(sglv), 'arglv': float(arglv),
        'cn1': None if cn1 is None else float(cn1),
        'cn2': None if cn2 is None else float(cn2),
        'cnt2': None if cnt2 is None else float(cnt2),
        'cnt2_consistent': (None if cnt2_consistent is None
                            else float(cnt2_consistent)),
        'round_leading_edge': bool(round_le),
        'spans_reads_wing_block': True,
        'figure_56a_second_call_uses_truncated_grid': truncated,
        'cnt2_divisor_inverted': inverted,
        'method': 'legacy_vtlift',
    }
