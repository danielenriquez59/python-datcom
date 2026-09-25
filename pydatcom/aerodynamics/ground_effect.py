"""
GRDEFF: ground effects on the wing-body-tail buildup.

The routine works by shifting the angle-of-attack schedule rather than the
coefficients.  It computes the wing's average height above the ground, reads
an induced-incidence increment out of Section 4.7.1's figures, and then
*re-reads* the free-air lift curve at the shifted angles.  Pitching moment,
drag and the normal/axial components follow from that new lift.

The five figures translated here are

=====================  =================================================
figure                  what it supplies
=====================  =================================================
4.7.1-14                the ``X`` factor, versus DX/(b/2) and h/(b/2)
4.7.1-15                ``LOLOM1``, versus CL/cos and h/CR
4.7.1-17                the flap increment ``DDCLF``, versus h(MAC)/CR
4.7.1-18A               the effective span ratio ``BWOB``
4.7.1-21                the low-aspect-ratio term ``BW``
=====================  =================================================

Figure 4.7.1-14 and 4.7.1-17 are read only on the high-aspect-ratio path
(theoretical AR >= 3); figure 4.7.1-21 only on the low-aspect-ratio path.

The source's ``CIOM`` save/restore block is a multiple-ground-height
bookkeeping device: on the first height it copies the free-air drag, lift,
moment and the two slope arrays aside, and on every height it copies them
back so each height starts from free air.  The saved set is exactly the set
the routine accumulates into.  This translation takes free-air arrays as
inputs and never mutates them, which has the same effect.

Reference: datcom-legacy/datcom_2000/grdeff.f
"""

import numpy as np
from typing import Dict, Optional, Sequence
import logging

from pydatcom.utils.constants import PI, RAD, UNUSED
from pydatcom.utils.legacy_tables import tlinex
from pydatcom.utils.legacy_interp import interx
from pydatcom.utils.legacy_numeric import tbfunx

logger = logging.getLogger(__name__)

# The source's TYPE selector for a straight tapered planform, DATA STRA.
STRAIGHT_TAPERED = 1.0

# Above this theoretical aspect ratio the high-AR path runs, below it the
# low-AR path.  The source writes IF(A(120).LT.3.0).
_LOW_ASPECT_RATIO = 3.0

# FTYPE values that read figure 4.7.1-17 at all, and the two curves.
_FLAP_TYPES_CURVE_A = (3, 4)
_FLAP_TYPE_CURVE_B = 5

# ---------------------------------------------------------------------------
# Figure 4.7.1-14:  X218 = HWOB2, X118 = DX, Y18 = X
# ---------------------------------------------------------------------------
_F14_HWOB2 = np.array([0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0])
_F14_DXOB2 = np.array([1.0, 0.5, 0.2, 0.0, -0.2, -0.5, -1.0])
_F14_X = np.array([
    1.4, 1.01, 0.79, 0.62, 0.50, 0.40, 0.33, 0.27, 0.22, 0.18, 0.15,
    1.27, 0.90, 0.69, 0.54, 0.42, 0.33, 0.28, 0.22, 0.18, 0.15, 0.12,
    1.11, 0.78, 0.58, 0.45, 0.35, 0.28, 0.22, 0.18, 0.15, 0.12, 0.099,
    1.00, 0.66, 0.50, 0.38, 0.30, 0.23, 0.19, 0.16, 0.12, 0.10, 0.085,
    0.86, 0.55, 0.40, 0.31, 0.24, 0.19, 0.16, 0.13, 0.10, 0.085, 0.080,
    0.72, 0.41, 0.29, 0.21, 0.16, 0.13, 0.10, 0.080, 0.070, 0.060, 0.050,
    0.60, 0.30, 0.19, 0.12, 0.085, 0.065, 0.045, 0.035, 0.025, 0.020,
    0.0200,
]).reshape((11, 7), order='F')

# ---------------------------------------------------------------------------
# Figure 4.7.1-15:  X219 = HWCOCR, X119 = CLOCOS, Y19 = LOLOM1
# ---------------------------------------------------------------------------
_F15_HWCOCR = np.array([.3, .4, .6, .8, 1., 1.2, 1.4, 1.6, 1.8, 2., 2.2, 2.4])
_F15_CLOCOS = np.array([0.0, 5.0, 10.0, 15.0, 18.0, 20.0, 22.0, 24.0, 36.0])
_F15_LOLOM1 = np.array([
    .40, .27, .145, .090, .060, .04, .030, .022, .015, .013, .010, .010,
    .25, .16, .065, .035, .017, .007, 0.0, -.005, -.007, -.010, -.012, -.016,
    .11, .060, 0.0, -.020, -.030, -.030, -.030, -.030, -.030, -.030, -.030,
    -.030,
    -.020, -.040, -.065, -.070, -.070, -.065, -.060, -.055, -.050, -.048,
    -.045,
    -.040, -.080, -.10, -.115, -.11, -.098, -.085, -.080, -.070, -.065, -.060,
    -.060, -.055, -.125, -.135, -.140, -.125, -.115, -.10, -.090, -.085, -.075,
    -.070, -.065, -.063, -.160, -.165, -.165, -.15, -.13, -.12, -.105, -.095,
    -.085, -.080, -.070, -.067, -.20, -.20, -.190, -.170, -.148, -.130, -.115,
    -.10, -.090, -.085, -.077, -.073, -.20, -.20, -.20, -.20, -.20, -.19, -.17,
    -.155, -0.138, -.127, -.118, -.11,
]).reshape((12, 9), order='F')

# ---------------------------------------------------------------------------
# Figure 4.7.1-18A:  X222A = 1/taper, X122A = aspect ratio, Y22A = BWOB
# ---------------------------------------------------------------------------
_F18A_INVERSE_TAPER = np.array([1.0, 1.5, 2.0, 3.0, 4.0, 5.0])
_F18A_ASPECT_RATIO = np.array([4.0, 6.0, 8.0, 10.0])
_F18A_BWOB = np.array([
    .825, .79, .77, .745, .725, .715,
    .853, .805, .77, .740, .715, .70,
    .88, .82, .77, .73, .703, .69,
    .895, .825, .77, .725, .695, .68,
]).reshape((6, 4), order='F')

# ---------------------------------------------------------------------------
# Figure 4.7.1-21:  X225 = HWOCBR, X125 = wing CL, Y25 = BW
# ---------------------------------------------------------------------------
_F21_HWOCBR = np.array([.2, .3, .4, .5, .6, .7, .8, .9, 1.0, 1.1, 1.2])
_F21_CL = np.array([0.0, .2, .4, .6, .8, 1.0, 1.2, 1.4, 1.6])
_F21_BW = np.array([
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    .92, .59, .41, .31, .23, .19, .15, .11, .09, .08, .07,
    1.92, 1.13, .80, .60, .45, .35, .25, .20, .17, .15, .12,
    2.45, 1.65, 1.15, .81, .60, .43, .35, .23, .20, .16, .12,
    2.6, 2.15, 1.42, 1.0, .70, .50, .37, .25, .19, .15, .14,
    2.6, 2.5, 1.7, 1.15, .78, .52, .37, .22, .14, .09, .04,
    2.6, 2.6, 1.85, 1.20, .76, .46, .29, .18, .08, 0.0, -.02,
    2.6, 2.6, 1.95, 1.20, .72, .39, .19, .03, -.02, -.09, -.12,
    2.6, 2.6, 2.10, 1.16, .52, .20, .03, -.07, -.16, -.22, -.26,
]).reshape((11, 9), order='F')

# ---------------------------------------------------------------------------
# Figure 4.7.1-17:  the flap DEL(DEL-CL) increment versus h(MAC)/CR
# ---------------------------------------------------------------------------
_F17_HMAC_OVER_CR = np.array([0., .1, .2, .3, .4, .5, .6, .7, .8, .9, 1., 1.1,
                              1.2])
_F17_CURVE_A = np.array([-.145, -.125, -.1, -.08, -.062, -.05, -.038, -.029,
                         -.02, -.014, -.005, 0., 0.])
_F17_CURVE_B = np.array([.098, .08, .062, .051, .04, .032, .024, .016, .01,
                         .006, 0., 0., 0.])


def _unused(value: Optional[float]) -> bool:
    """Whether a dihedral entry is the source's UNUSED sentinel."""
    return value is None or float(value) == UNUSED


def ground_effect_geometry(ground_height: float,
                           wing: Dict[str, float],
                           synthesis: Dict[str, float],
                           sweep: Dict[str, float],
                           wing_theoretical: Dict[str, float],
                           tail: Optional[Dict[str, float]] = None,
                           tail_theoretical: Optional[Dict[str, float]] = None
                           ) -> Dict[str, object]:
    """Translate GRDEFF's geometric parameter block.

    Args:
        ground_height: ``GRDHT``, the landing-gear height above the ground.
        wing: ``/WINGI/`` entries as ``{'sspnop', 'sspn', 'chrdr', 'sspndd',
            'dhdadi', 'dhdado', 'planform_type'}``.  A dihedral of ``None``
            stands for the source's UNUSED sentinel.
        synthesis: ``/SYNTSS/`` entries ``{'zw', 'aliw'}``, plus ``{'zh',
            'alih'}`` when a horizontal tail is present.
        sweep: Quarter-chord sweep data, ``{'tan_c4', 'tan_c4_inboard',
            'tan_c4_outboard', 'zero_sweep_station'}``.  These are the
            source's ``A(44)``, ``A(68)``, ``A(92)`` and ``A(8)``.
        wing_theoretical: ``{'mac', 'y_mac', 'mac_c4'}``, the source's
            ``A(122)``, ``A(136)`` and ``A(161)``.
        tail: ``/HTI/`` entries ``{'sspn', 'sspndd', 'dhdadi', 'dhdado'}``.
        tail_theoretical: ``{'y_mac', 'mac_c4'}``, ``AHT(136)``/``AHT(161)``.

    Returns:
        Dictionary with ``dx``, ``dxob2``, ``h75cr``, ``hw``, ``hwob2``,
        ``hwcr4``, ``hwcocr``, ``hwmacx``, ``hwmac4``, ``htmacx``,
        ``htmac4``, ``r``, ``sigma``, ``hwocbr`` and ``t``.  The tail
        entries are absent when no tail is supplied.

    Raises:
        ValueError: If the root chord, semispan or theoretical MAC is
            nonpositive.

    Notes:
        The ``DX`` branch taken when the planform is not straight tapered
        and ``SSPNOP <= 0.25*SSPN`` multiplies a *span* by ``A(8)``.  But
        ``A(8)`` is set by ``INFTGM`` under the banner "DETERMINE LOCATION
        OF ZERO SWEEP ANGLE" and takes values ``1.``, ``.25*CP`` and
        ``A(110)/A(13)+CHSTAT`` -- a chord fraction, not a tangent.  Its two
        sibling branches both multiply the same span by a sweep tangent
        (``A(44)`` for the straight planform, ``A(68)`` for the inboard
        panel), and when ``SSPNOP <= 0.25*SSPN`` the three-quarter-semispan
        station lies wholly inside the inboard panel, so the parallel
        construction calls for ``A(68)``.  ``A(8)`` against ``A(68)`` is a
        single-character slip.

        The source form is preserved, because changing it changes results.
        ``dx_uses_zero_sweep_station`` flags when that branch ran, and
        ``dx_with_inboard_tangent`` reports what the parallel branches would
        have produced.
    """
    sspn = float(wing['sspn'])
    chrdr = float(wing['chrdr'])
    sspnop = float(wing.get('sspnop', 0.0) or 0.0)
    sspndd = float(wing.get('sspndd', 0.0) or 0.0)
    if sspn <= 0.0 or chrdr <= 0.0:
        raise ValueError("GRDEFF requires a positive root chord and semispan")
    mac = float(wing_theoretical['mac'])
    if mac <= 0.0:
        raise ValueError("GRDEFF requires a positive theoretical MAC")

    dhdadi = wing.get('dhdadi')
    dhdado = wing.get('dhdado')
    tan_di = 0.0 if _unused(dhdadi) else float(np.tan(float(dhdadi) / RAD))
    tan_do = 0.0 if _unused(dhdado) else float(np.tan(float(dhdado) / RAD))

    straight = float(wing.get('planform_type', STRAIGHT_TAPERED) or
                     STRAIGHT_TAPERED) == STRAIGHT_TAPERED

    # ---- DELTAX --------------------------------------------------------
    zero_sweep_used = False
    if straight:
        dx = 0.5 * chrdr - float(sweep['tan_c4']) * 0.75 * sspn
    elif sspnop > 0.25 * sspn:
        dx = (0.5 * chrdr
              - float(sweep['tan_c4_outboard']) * (sspnop - 0.25 * sspn)
              - float(sweep['tan_c4_inboard']) * (sspn - sspnop))
    else:
        # See the note above: the source multiplies a span by A(8).
        dx = 0.5 * chrdr - 0.75 * sspn * float(sweep['zero_sweep_station'])
        zero_sweep_used = True
    dxob2 = dx / sspn

    # ---- average elevation of the wing above the ground ----------------
    aliw = float(synthesis['aliw'])
    tan_incidence = float(np.tan(aliw / RAD))
    h75cr = ground_height + float(synthesis['zw']) - 0.75 * chrdr * tan_incidence
    incidence_term = dx * tan_incidence * 0.5

    if _unused(dhdadi) and _unused(dhdado):
        hw = h75cr + incidence_term
    elif not _unused(dhdadi) and _unused(dhdado):
        hw = h75cr + 0.375 * sspn * tan_di + incidence_term
    elif _unused(dhdadi) and not _unused(dhdado):
        if sspndd <= 0.25 * sspn:
            hw = h75cr + incidence_term
        else:
            hw = (h75cr + 0.5 * tan_do * (sspndd - 0.25 * sspn) +
                  incidence_term)
    else:
        if sspndd <= 0.25 * sspn:
            hw = h75cr + 0.375 * sspn * tan_di + incidence_term
        else:
            hw = (h75cr + 0.5 * ((sspn - sspndd) * tan_di +
                                 (sspndd - 0.25 * sspn) * tan_do) +
                  incidence_term)
    hwob2 = hw / sspn

    hwcr4 = h75cr + 0.5 * chrdr * tan_incidence
    hwcocr = hwcr4 / chrdr

    result = {
        'dx': float(dx),
        'dxob2': float(dxob2),
        'h75cr': float(h75cr),
        'hw': float(hw),
        'hwob2': float(hwob2),
        'hwcr4': float(hwcr4),
        'hwcocr': float(hwcocr),
        'dx_uses_zero_sweep_station': zero_sweep_used,
        'dx_with_inboard_tangent': float(
            0.5 * chrdr - 0.75 * sspn * float(sweep.get('tan_c4_inboard', 0.0))
        ) if zero_sweep_used else None,
    }

    # ---- MAC elevations -------------------------------------------------
    # The source computes HWMACX/HWMAC4 only inside its IF(HTPL) block, but
    # the figure 4.7.1-17 flap lookup reads HWMAC4 whenever the aspect ratio
    # is at least 3 and FTYPE is 3, 4 or 5, with no tail test.  HWMAC4 is
    # GR(9) in a COMMON block, so with no tail that read returns zero on the
    # first case and the previous case's value on every later one.  The
    # formula needs only wing quantities, so it is evaluated unconditionally
    # here; that agrees with the source wherever the source defines it.
    hwmacx = (ground_height + float(synthesis['zw']) -
              float(wing_theoretical['mac_c4']) * tan_incidence)
    result['hwmacx'] = float(hwmacx)
    result['hwmac4'] = float(_mac_elevation(
        hwmacx, float(wing_theoretical['y_mac']), sspn, sspndd,
        dhdadi, dhdado, tan_di, tan_do))
    result['hwmac4_computed_without_tail'] = tail is None

    if tail is not None and tail_theoretical is not None:
        alih = float(synthesis['alih'])
        tan_incidence_h = float(np.tan(alih / RAD))
        htmacx = (ground_height + float(synthesis['zh']) -
                  float(tail_theoretical['mac_c4']) * tan_incidence_h)
        tail_dhdadi = tail.get('dhdadi')
        tail_dhdado = tail.get('dhdado')
        result['htmacx'] = float(htmacx)
        result['htmac4'] = float(_mac_elevation(
            htmacx, float(tail_theoretical['y_mac']), float(tail['sspn']),
            float(tail.get('sspndd', 0.0) or 0.0), tail_dhdadi, tail_dhdado,
            0.0 if _unused(tail_dhdadi) else float(np.tan(float(tail_dhdadi) / RAD)),
            0.0 if _unused(tail_dhdado) else float(np.tan(float(tail_dhdado) / RAD))))

    # ---- lift-related height parameters --------------------------------
    result['r'] = float((1.0 + hwob2**2)**0.5 - hwob2)
    result['sigma'] = float(np.exp(-2.48 * hwob2**0.768))
    hwocbr = hw / mac
    result['hwocbr'] = float(hwocbr)
    result['t'] = float((RAD / (8.0 * PI)) *
                        (hwocbr / (hwocbr**2 + 1.0 / 64.0)))
    return result


def _mac_elevation(mac_x: float, y_mac: float, sspn: float, sspndd: float,
                   dhdadi, dhdado, tan_di: float, tan_do: float) -> float:
    """The MAC quarter-chord elevation, GRDEFF's labels 1100-1180.

    The wing and the horizontal tail share this dihedral branch structure
    exactly, differing only in which COMMON block supplies the arguments.
    """
    if _unused(dhdadi) and _unused(dhdado):
        return mac_x
    if not _unused(dhdadi) and _unused(dhdado):
        return mac_x + y_mac * tan_di
    if _unused(dhdadi) and not _unused(dhdado):
        if y_mac <= sspn - sspndd:
            return mac_x
        return mac_x + (y_mac + sspndd - sspn) * tan_do
    if y_mac <= sspn - sspndd:
        return mac_x + y_mac * tan_di
    return (mac_x + (sspn - sspndd) * tan_di +
            (y_mac + sspndd - sspn) * tan_do)


def figure_4711_14(dxob2: float, hwob2: float) -> float:
    """Figure 4.7.1-14: the ``X`` induced-incidence factor."""
    return float(tlinex(_F14_DXOB2, _F14_HWOB2, _F14_X, dxob2, hwob2,
                        0, 0, 0, 0))


def figure_4711_15(clocos: float, hwcocr: float) -> float:
    """Figure 4.7.1-15: ``LOLOM1`` versus CL/cos^2 and the c/4 height."""
    return float(tlinex(_F15_CLOCOS, _F15_HWCOCR, _F15_LOLOM1, clocos,
                        hwcocr, 0, 0, 2, 2))


def figure_4711_17(flap_type: int, hwmac4_over_cr: float) -> float:
    """Figure 4.7.1-17: ``DDCLF``, the flap lift-increment correction.

    FTYPE 3 and 4 read one curve, FTYPE 5 the other; every other flap type
    leaves the increment at zero, as the source's range test does.
    """
    kind = int(flap_type)
    if kind in _FLAP_TYPES_CURVE_A:
        curve = _F17_CURVE_A
    elif kind == _FLAP_TYPE_CURVE_B:
        curve = _F17_CURVE_B
    else:
        return 0.0
    return float(interx(1, _F17_HMAC_OVER_CR, [hwmac4_over_cr], [13], curve,
                        lind=13, lx1l=1, lx1u=1))


def figure_4711_18a(aspect_ratio: float, inverse_taper: float) -> float:
    """Figure 4.7.1-18A: the effective span ratio ``BWOB``."""
    return float(tlinex(_F18A_ASPECT_RATIO, _F18A_INVERSE_TAPER, _F18A_BWOB,
                        aspect_ratio, inverse_taper, 2, 0, 2, 1))


def figure_4711_21(wing_cl: float, hwocbr: float) -> float:
    """Figure 4.7.1-21: the low-aspect-ratio ``BW`` term."""
    return float(tlinex(_F21_CL, _F21_HWOCBR, _F21_BW, wing_cl, hwocbr,
                        0, 0, 2, 2))


def ground_effect_incidence(geometry: Dict[str, object],
                            alpha_deg: Sequence[float],
                            wing: Dict[str, float],
                            sweep: Dict[str, float],
                            wing_theoretical: Dict[str, float],
                            wing_alone: Dict[str, Sequence[float]],
                            wing_body: Dict[str, Sequence[float]],
                            flap: Optional[Dict[str, float]] = None,
                            has_horizontal_tail: bool = False
                            ) -> Dict[str, object]:
    """Translate GRDEFF's angle-of-attack increment, both aspect-ratio paths.

    Args:
        geometry: The result of :func:`ground_effect_geometry`.
        alpha_deg: The free-air angle schedule, ``FLC(23)`` onward.
        wing: ``/WINGI/`` entries; ``tovc`` is read on the low-AR path.
        sweep: Quarter-chord sweep data; ``cos_c4``, ``cos_c4_inboard`` and
            ``cos_c4_outboard`` are the source's ``A(43)``, ``A(67)`` and
            ``A(91)``.
        wing_theoretical: ``{'aspect_ratio'}``, the source's ``A(120)``.
        wing_alone: ``{'cl'}``, the ``/IWING/`` lift array, plus
            ``'flap_dcl'`` for the source's scalar ``WING(L+200)``.
        wing_body: ``{'cl', 'cla'}``, free-air ``BWI(21)`` and ``BWI(101)``.
        flap: ``{'type', 'deflection'}``, the source's ``F(17)`` and
            ``DELTA(L)``.  Omit when no flap is deflected.
        has_horizontal_tail: The source's ``HTPL``.

    Returns:
        Dictionary with ``dalpha`` and ``alphwg``, the figure values used,
        and ``path`` naming which aspect-ratio branch ran.

    Notes:
        The flap lift increment ``WING(L+200)`` is added to the wing lift
        only when there is *no* horizontal tail.  That asymmetry is the
        source's and is preserved.
    """
    alpha = np.atleast_1d(np.asarray(alpha_deg, dtype=float))
    aspect_ratio = float(wing_theoretical['aspect_ratio'])
    sspn = float(wing['sspn'])
    chrdr = float(wing['chrdr'])
    cl_wing = np.atleast_1d(np.asarray(wing_alone['cl'], dtype=float))
    cl_body_wing = np.atleast_1d(np.asarray(wing_body['cl'], dtype=float))
    cla_body_wing = np.atleast_1d(np.asarray(wing_body['cla'], dtype=float))
    for name, array in (('wing CL', cl_wing), ('body-wing CL', cl_body_wing),
                        ('body-wing CLA', cla_body_wing)):
        if array.shape != alpha.shape:
            raise ValueError(f"GRDEFF needs a {name} matching the schedule")

    r = float(geometry['r'])
    dalpha = np.zeros_like(alpha)
    alphwg = np.zeros_like(alpha)

    if aspect_ratio >= _LOW_ASPECT_RATIO:
        factor_x = figure_4711_14(float(geometry['dxob2']),
                                  float(geometry['hwob2']))
        straight = float(wing.get('planform_type', STRAIGHT_TAPERED) or
                         STRAIGHT_TAPERED) == STRAIGHT_TAPERED
        if straight:
            cosl4 = float(sweep['cos_c4'])
        else:
            sspnop = float(wing.get('sspnop', 0.0) or 0.0)
            cosl4 = ((float(sweep['cos_c4_outboard']) * sspnop +
                      float(sweep['cos_c4_inboard']) * (sspn - sspnop)) / sspn)

        ddclf = 0.0
        deflection = 0.0
        if flap is not None:
            # The source takes IFTYPE=F(17)+0.5, a rounded integer.
            flap_type = int(float(flap['type']) + 0.5)
            ddclf = figure_4711_17(flap_type,
                                   float(geometry['hwmac4']) / chrdr)
            deflection = float(flap.get('deflection', 0.0) or 0.0)

        clocos = np.zeros_like(alpha)
        lolom1 = np.zeros_like(alpha)
        flap_dcl = float(wing_alone.get('flap_dcl', 0.0) or 0.0)
        for angle_slot in range(alpha.size):
            body_wing_cl = cl_body_wing[angle_slot]
            if not has_horizontal_tail:
                body_wing_cl = body_wing_cl + flap_dcl
            clocos[angle_slot] = (RAD * cl_wing[angle_slot] /
                                  (2.0 * PI * cosl4**2))
            lolom1[angle_slot] = figure_4711_15(
                float(clocos[angle_slot]), float(geometry['hwcocr']))
            dalpha[angle_slot] = (
                -(9.12 / aspect_ratio + 7.16 * chrdr / (2.0 * sspn)) *
                body_wing_cl * factor_x
                - (aspect_ratio * chrdr /
                   (4.0 * cla_body_wing[0] * sspn)) *
                lolom1[angle_slot] * body_wing_cl * r)
            alphwg[angle_slot] = (
                alpha[angle_slot] + dalpha[angle_slot] -
                ddclf * deflection**2 /
                (2500.0 * cla_body_wing[angle_slot]))
        return {
            'path': 'high_aspect_ratio',
            'dalpha': dalpha,
            'alphwg': alphwg,
            'x': float(factor_x),
            'cosl4': float(cosl4),
            'ddclf': float(ddclf),
            'clocos': clocos,
            'lolom1': lolom1,
            'method': 'legacy_grdeff',
        }

    # ---- low-aspect-ratio path -----------------------------------------
    hwocbr = float(geometry['hwocbr'])
    sigma = float(geometry['sigma'])
    t = float(geometry['t'])
    k = (RAD * 0.0030 * hwocbr *
         (1.0 / (hwocbr**2 + 1.0 / 64.0)**2 +
          1.0 / (hwocbr**2 + 9.0 / 64.0)**2))
    # WINGIN(16) is TOVC by the /WINGI/ declaration in inputc.f.
    tovc = float(wing.get('tovc', 0.0) or 0.0)
    bw = np.zeros_like(alpha)
    for angle_slot in range(alpha.size):
        bw[angle_slot] = figure_4711_21(
            float(cl_wing[angle_slot]), hwocbr)
        dalpha[angle_slot] = (
            -18.24 * cl_body_wing[angle_slot] * sigma / aspect_ratio +
            r * t * cl_body_wing[angle_slot]**2 /
            (RAD * cla_body_wing[0]) -
            r * bw[angle_slot] + k * tovc)
        alphwg[angle_slot] = alpha[angle_slot] + dalpha[angle_slot]
    return {
        'path': 'low_aspect_ratio',
        'dalpha': dalpha,
        'alphwg': alphwg,
        'k': float(k),
        'bw': bw,
        'method': 'legacy_grdeff',
    }


def ground_effect_tail(geometry: Dict[str, object],
                       alpha_deg: Sequence[float],
                       wing: Dict[str, float],
                       wing_theoretical: Dict[str, float],
                       wing_body: Dict[str, Sequence[float]],
                       wing_body_tail: Dict[str, Sequence[float]],
                       downwash: Dict[str, Sequence[float]]
                       ) -> Dict[str, object]:
    """Translate GRDEFF's tail block, labels 1260-1280.

    The downwash angle is reduced by the classic image-vortex ratio
    ``(b^2 + 4(h_t - h_w)^2) / (b^2 + 4(h_t + h_w)^2)`` over the effective
    span from figure 4.7.1-18A.

    Args:
        geometry: The result of :func:`ground_effect_geometry`, which must
            carry the tail entries.
        alpha_deg: The free-air angle schedule.
        wing: ``/WINGI/`` entries; ``sspn`` sets the effective span.
        wing_theoretical: ``{'aspect_ratio', 'taper_ratio'}``, the source's
            ``A(120)`` and ``A(118)``.  Figure 4.7.1-18A is entered on
            ``1/A(118)``.
        wing_body: ``{'cl'}``, free-air ``BWI(21)``.
        wing_body_tail: ``{'cl'}``, free-air ``BWH(21)``.
        downwash: ``{'angle'}``, the source's ``DWASH(21)`` block.

    Returns:
        Dictionary with ``bwob``, ``beff``, ``ddwash``, ``clht`` and
        ``alphat``.

    Raises:
        ValueError: If the theoretical taper ratio is zero, which the
            source would divide by.
    """
    alpha = np.atleast_1d(np.asarray(alpha_deg, dtype=float))
    taper = float(wing_theoretical['taper_ratio'])
    if taper == 0.0:
        raise ValueError(
            "GRDEFF enters figure 4.7.1-18A on 1/taper; the theoretical "
            "taper ratio A(118) cannot be zero")
    bwob = figure_4711_18a(float(wing_theoretical['aspect_ratio']),
                           1.0 / taper)
    beff = bwob * 2.0 * float(wing['sspn'])

    angle = np.atleast_1d(np.asarray(downwash['angle'], dtype=float))
    cl_body_wing = np.atleast_1d(np.asarray(wing_body['cl'], dtype=float))
    cl_body_wing_tail = np.atleast_1d(
        np.asarray(wing_body_tail['cl'], dtype=float))
    htmac4 = float(geometry['htmac4'])
    hwmac4 = float(geometry['hwmac4'])
    ratio = ((beff**2 + 4.0 * (htmac4 - hwmac4)**2) /
             (beff**2 + 4.0 * (htmac4 + hwmac4)**2))
    ddwash = angle * ratio
    clht = cl_body_wing_tail - cl_body_wing
    alphat = alpha - ddwash
    return {
        'bwob': float(bwob),
        'beff': float(beff),
        'ratio': float(ratio),
        'ddwash': ddwash,
        'clht': clht,
        'alphat': alphat,
        'method': 'legacy_grdeff',
    }


_BUILDUPS = ('wing_body', 'wing_body_vertical', 'wing_body_tail',
             'wing_body_tail_vertical')


def _free_air(buildup: Dict[str, Dict[str, Sequence[float]]], name: str,
              count: int) -> Dict[str, np.ndarray]:
    """One buildup's free-air arrays, copied so the caller's are untouched."""
    block = buildup.get(name)
    if block is None:
        return {key: np.zeros(count)
                for key in ('cd', 'cl', 'cm', 'cla', 'cma')}
    out = {}
    for key in ('cd', 'cl', 'cm', 'cla', 'cma'):
        array = np.atleast_1d(np.asarray(block.get(key, np.zeros(count)),
                                         dtype=float)).copy()
        if array.size != count:
            raise ValueError(
                f"GRDEFF needs {name}.{key} to match the angle schedule")
        out[key] = array
    return out


def calculate_grdeff(ground_height: float,
                     alpha_deg: Sequence[float],
                     wing: Dict[str, float],
                     synthesis: Dict[str, float],
                     sweep: Dict[str, float],
                     wing_theoretical: Dict[str, float],
                     wing_alone: Dict[str, Sequence[float]],
                     buildup: Dict[str, Dict[str, Sequence[float]]],
                     tail: Optional[Dict[str, float]] = None,
                     tail_theoretical: Optional[Dict[str, float]] = None,
                     downwash: Optional[Dict[str, Sequence[float]]] = None,
                     flap: Optional[Dict[str, float]] = None,
                     has_vertical_panel: bool = False) -> Dict[str, object]:
    """Translate GRDEFF: ground effects on the full buildup.

    A horizontal tail is taken to be present exactly when ``tail``,
    ``tail_theoretical`` and ``downwash`` are all supplied, which is the
    source's ``HTPL``.

    Args:
        ground_height: ``GRDHT``.
        alpha_deg: The free-air angle schedule, ``FLC(23)`` onward.
        wing: ``/WINGI/`` entries; see :func:`ground_effect_geometry` and
            :func:`ground_effect_incidence`.
        synthesis: ``/SYNTSS/`` entries ``{'xcg', 'zw', 'aliw'}`` plus
            ``{'xh', 'zh', 'alih'}`` with a tail.
        sweep: Quarter-chord sweep data.
        wing_theoretical: ``{'aspect_ratio', 'taper_ratio', 'mac', 'y_mac',
            'mac_c4'}``.
        wing_alone: ``{'cd', 'cl'}`` from ``/IWING/``, plus ``'flap_dcl'``.
        buildup: Free-air ``{'wing_body', 'wing_body_vertical',
            'wing_body_tail', 'wing_body_tail_vertical'}``, each a dict of
            ``cd``, ``cl``, ``cm``, ``cla`` and ``cma``.  These correspond
            to ``BWI``, ``BWV``, ``BWH`` and ``BWHV``, and are never
            mutated.
        tail: ``/HTI/`` entries for the horizontal tail.
        tail_theoretical: ``{'y_mac', 'mac_c4'}`` for the tail.
        downwash: ``{'slope', 'angle'}``, the ``DWASH(1)`` and ``DWASH(21)``
            blocks.
        flap: ``{'type', 'deflection'}``.
        has_vertical_panel: The source's ``VTPL .OR. VFPL .OR. TVTPAN``.

    Returns:
        Dictionary with a ``buildup`` of the four in-ground-effect blocks
        (each carrying ``cd``, ``cl``, ``cm``, ``cn``, ``ca``, ``cla`` and
        ``cma``), the ``geometry`` and ``incidence`` sub-results, the tail
        sub-result when a tail is present, and the increments ``dclwbg``,
        ``dcmwbg``, ``dclhtg``, ``dcmhtg`` and ``dcdlwg``.

    Raises:
        ValueError: If the schedules mismatch, or if the shifted angle
            schedule is not strictly increasing.  The source hands
            ``ALPHWG`` and ``ALPHAT`` to ``TBFUNX`` as abscissas; a large
            enough incidence increment can reverse them, which the source
            would silently interpolate through.

    Notes:
        The source aliases three of its working arrays onto the buildup
        blocks through EQUIVALENCE: ``CLG`` onto ``BWH(21)``, ``CMWBG``
        onto ``BWI(41)`` and ``CMG`` onto ``BWH(41)``.  So the lines that
        read like fresh assignments are in-place accumulations, and the
        body-wing-tail lift is written as ``BWV(J+20) + CLHTG(J)`` -- built
        on the *vertical-panel* buildup's lift, not the body-wing one.
        With a vertical panel present those agree, because ``BWV(J+20)`` is
        assigned ``CLWBG(J)`` earlier in the same loop.  Without one the
        source reads a slot its own restore block never refreshes; here the
        caller's ``wing_body_vertical`` lift is used and
        ``tail_lift_read_stale_vertical`` flags it.
    """
    alpha = np.atleast_1d(np.asarray(alpha_deg, dtype=float))
    count = alpha.size
    has_tail = (tail is not None and tail_theoretical is not None and
                downwash is not None)

    geometry = ground_effect_geometry(
        ground_height, wing, synthesis, sweep, wing_theoretical,
        tail=tail if has_tail else None,
        tail_theoretical=tail_theoretical if has_tail else None)

    blocks = {name: _free_air(buildup, name, count) for name in _BUILDUPS}
    bwi, bwv = blocks['wing_body'], blocks['wing_body_vertical']
    bwh, bwhv = blocks['wing_body_tail'], blocks['wing_body_tail_vertical']

    incidence = ground_effect_incidence(
        geometry, alpha, wing, sweep, wing_theoretical, wing_alone,
        {'cl': bwi['cl'], 'cla': bwi['cla']}, flap=flap,
        has_horizontal_tail=has_tail)
    alphwg = incidence['alphwg']

    tail_result = None
    if has_tail:
        tail_result = ground_effect_tail(
            geometry, alpha, wing, wing_theoretical, {'cl': bwi['cl']},
            {'cl': bwh['cl']}, downwash)

    # ---- loop 1300: re-read the lift curves at the shifted angles ------
    free_air_cl = bwi['cl'].copy()
    clwbg = np.array([tbfunx(alphwg, free_air_cl, angle_deg, 1, 2)[0]
                      for angle_deg in alpha])
    dclwbg = clwbg - free_air_cl
    if has_vertical_panel:
        bwv['cl'] = clwbg.copy()
    stale_vertical = False
    clhtg = None
    if has_tail:
        clhtg = np.array([
            tbfunx(tail_result['alphat'], tail_result['clht'],
                   angle_deg, 1, 2)[0]
            for angle_deg in alpha])
        # CLG is EQUIVALENCEd onto BWH(21); the sum is built on BWV.
        stale_vertical = not has_vertical_panel
        bwh['cl'] = bwv['cl'] + clhtg
        if has_vertical_panel:
            bwhv['cl'] = bwh['cl'].copy()
    # Loop 1310.
    bwi['cl'] = clwbg

    # ---- loop 1350: pitching moment, then drag -------------------------
    dxcp = float(bwi['cma'][0] / bwi['cla'][0])      # BWI(121)/BWI(101)
    dcmwbg = dxcp * dclwbg
    bwi['cm'] = bwi['cm'] + dcmwbg                   # CMWBG is BWI(41)
    if has_vertical_panel:
        bwv['cm'] = bwi['cm'].copy()

    dclhtg = dcmhtg = None
    if has_tail:
        lh = (float(synthesis['xh']) + float(tail_theoretical['mac_c4']) -
              float(synthesis['xcg']))
        lhocbr = lh / float(wing_theoretical['mac'])
        slope = np.atleast_1d(np.asarray(downwash['slope'], dtype=float))
        dclhtg = clhtg - tail_result['clht']
        dcmhtg = -dclhtg * lhocbr * slope
        bwh['cm'] = bwh['cm'] + dcmhtg               # CMG is BWH(41)
        if has_vertical_panel:
            bwhv['cm'] = bwh['cm'].copy()

    aspect_ratio = float(wing_theoretical['aspect_ratio'])
    sigma = float(geometry['sigma'])
    mirror_r = float(geometry['r'])
    mirror_t = float(geometry['t'])
    cd_wing = np.atleast_1d(np.asarray(wing_alone['cd'], dtype=float))
    cl_wing = np.atleast_1d(np.asarray(wing_alone['cl'], dtype=float))
    induced = sigma * cl_wing**2 / (PI * aspect_ratio)
    dcdlwg = (-induced -
              (cd_wing - induced) * mirror_r * mirror_t * cl_wing / RAD)
    bwi['cd'] = bwi['cd'] + dcdlwg
    if has_vertical_panel:
        bwv['cd'] = bwv['cd'] + dcdlwg
    if has_tail:
        bwh['cd'] = bwh['cd'] + dcdlwg
        if has_vertical_panel:
            bwhv['cd'] = bwhv['cd'] + dcdlwg

    # ---- loops 1360-1390: normal and axial force, then the slopes ------
    sin_a = np.sin(alpha / RAD)
    cos_a = np.cos(alpha / RAD)
    for block in (bwi, bwv, bwh, bwhv):
        block['cn'] = block['cl'] * cos_a + block['cd'] * sin_a
        block['ca'] = block['cd'] * cos_a - block['cl'] * sin_a
        block['cla'] = np.array([
            tbfunx(alpha, block['cl'], angle_deg, 0, 0)[1]
            for angle_deg in alpha])
        block['cma'] = np.array([
            tbfunx(alpha, block['cm'], angle_deg, 0, 0)[1]
            for angle_deg in alpha])

    result = {
        'buildup': blocks,
        'geometry': geometry,
        'incidence': incidence,
        'dclwbg': dclwbg,
        'dcmwbg': dcmwbg,
        'dcdlwg': dcdlwg,
        'dxcp': dxcp,
        'clwbg': clwbg,
        'has_horizontal_tail': has_tail,
        'tail_lift_read_stale_vertical': stale_vertical,
        'method': 'legacy_grdeff',
    }
    if has_tail:
        result['tail'] = tail_result
        result['clhtg'] = clhtg
        result['dclhtg'] = dclhtg
        result['dcmhtg'] = dcmhtg
    return result
