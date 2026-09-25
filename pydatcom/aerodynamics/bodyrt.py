"""
BODYRT: axisymmetric body-alone lift, drag and pitching moment.

A complete translation of the subsonic body routine, now possible because
every one of its five dependencies is translated: ``EQSPC1``, ``FIG26``,
``GETMAX``, ``TBFUNX`` and ``TRAPZ``.

The routine has two parts.  The setup, at source labels 1000 to 1090, builds
the potential-flow lift-curve slope from the apparent-mass factor of Figure
4.2.1.1-20, the moment slope from the resampled area distribution, and the
zero-lift drag from skin friction and base drag.  The angle loop at label
1100 adds Allen-Perkins crossflow through Figures 4.2.1.2-35A and -35B and
rotates the result into wind axes.

Reference: datcom-legacy/datcom_2000/bodyrt.f
"""

import numpy as np
from typing import Dict, Sequence
import logging

from pydatcom.utils.constants import PI, RAD, UNUSED
from pydatcom.utils.legacy_numeric import tbfunx, trapz
from pydatcom.utils.legacy_interp import eqspc1
from pydatcom.utils.table_lookup import fig26
from pydatcom.interactions.body_vortex import getmax

logger = logging.getLogger(__name__)

# Figure 4.2.1.1-20: apparent mass factor k2-k1 against fineness ratio.
_FIG_42110_20_X = np.array([4., 5., 6., 8., 10., 12., 14., 16., 18., 20.])
_FIG_42110_20_Y = np.array([.77, .825, .865, .91, .94, .955, .965, .97,
                            .973, .975])

# Figure 4.2.1.2-35A: steady-state crossflow drag proportionality factor
# against body fineness ratio.
_FIG_42120_35A_X = np.array([2., 4., 8., 12., 16., 20., 24., 28.])
_FIG_42120_35A_Y = np.array([.56, .6, .66, .71, .74, .76, .775, .79])

# Figure 4.2.1.2-35B: crossflow drag of a circular cylinder against
# crossflow Mach number.  The source writes three repeated entries as
# "3*1.8".
_FIG_42120_35B_X = np.array([0., .2, .3, .36, .4, .5, .6, .7, .77, .8, .86,
                             .9, .98, 1.])
_FIG_42120_35B_Y = np.array([1.2, 1.2, 1.21, 1.23, 1.27, 1.36, 1.5, 1.67,
                             1.75, 1.77, 1.8, 1.8, 1.8, 1.79])

# Figure 4.1.5.1-27: roughness-limited Reynolds number intercept vs Mach.
_FIG_41510_27_MACH = np.array([0.0, 1.0, 2.0, 3.0])
_FIG_41510_27_CEPT = np.array([1.57780, 1.67221, 1.98509, 2.28874])

# The source resamples onto a fixed 20-station grid.
_STATIONS = 20


def _body_reference_station(x, s):
    """Source labels 1000 to 1030: the reference station X1.

    ``X1`` is where the area distribution closes down fastest.  If the body
    never contracts, the source uses the last station instead.
    """
    x = np.asarray(x, dtype=float)
    s = np.asarray(s, dtype=float)
    contracts = any(
        s[station] < s[station - 1]
        for station in range(1, len(s))
    )
    if not contracts:
        return float(x[-1]), False
    # BD(K+174) = -dS/dx at every station, with TBFUNX end modes 2 and 1.
    slopes = np.array([
        -tbfunx(x, s, station, lower=2, upper=1)[1] for station in x
    ])
    x1, _, _ = getmax(x, slopes)
    return float(x1), True


def calculate_bodyrt(x: Sequence[float], s: Sequence[float],
                     p: Sequence[float], r: Sequence[float],
                     alpha_deg: Sequence[float], mach: float,
                     reynolds_per_length: float, sref: float, cbar: float,
                     xcg: float, roughness: float = 1.6e-4,
                     transonic: bool = False) -> Dict[str, object]:
    """Translate BODYRT for an axisymmetric body.

    Args:
        x: Station coordinates from the nose, length ``NX``.
        s: Cross-sectional area at each station.
        p: Perimeter at each station.
        r: Body half-width (radius) at each station.
        alpha_deg: Angle-of-attack schedule, degrees.
        mach: Free-stream Mach number, the source's ``B(1)``.
        reynolds_per_length: Reynolds number per unit length, ``FLC(M+42)``.
        sref: Reference area.
        cbar: Reference chord.
        xcg: Moment reference station, the source's ``BD(33)``.
        roughness: Surface roughness height, ``ROUGFC``.
        transonic: When set, the source returns after the drag buildup and
            computes no angle-dependent coefficients.

    Returns:
        Dictionary with ``cla`` and ``cma`` (per degree), the drag buildup,
        and per-angle arrays ``cn``, ``cm``, ``cd``, ``cl`` and ``ca`` unless
        ``transonic`` is set.

    Raises:
        ValueError: If the geometry is too short, inconsistent, or has a
            nonpositive reference quantity.
    """
    x = np.asarray(x, dtype=float)
    s = np.asarray(s, dtype=float)
    p = np.asarray(p, dtype=float)
    r = np.asarray(r, dtype=float)
    if x.ndim != 1 or len(x) < 2:
        raise ValueError("BODYRT needs at least two body stations")
    if not (x.shape == s.shape == p.shape == r.shape):
        raise ValueError("BODYRT station arrays must have matching lengths")
    if min(sref, cbar) <= 0.0:
        raise ValueError("BODYRT requires positive SREF and CBARR")
    if roughness <= 0.0:
        raise ValueError("BODYRT requires a positive roughness height")

    # --- Setup (source labels 1000–1090) ---
    length = float(x[-1])              # BD(1)
    base_area = float(s[-1])           # BD(57)
    roughness_length = 12.0 * length / roughness   # BD(55)
    nose_length = length               # BD(5): body alone sets LNOSE = L

    # GETMAX: maximum cross-sectional area.
    _, max_area, _ = getmax(x, s)      # BD(3) -> BD(56)
    if max_area <= 0.0:
        raise ValueError("BODYRT requires a positive maximum cross-section")
    # Boat-tailed bodies have their base area floored at 30% of maximum.
    base_area = max(base_area, 0.30 * max_area)

    x1, contracts = _body_reference_station(x, s)
    x0 = 0.378 * length + 0.527 * x1   # BD(7)

    if x0 <= nose_length:
        area_ref = float(tbfunx(x, s, x0, lower=0, upper=0)[0])  # BD(6)
        tmp1, tmp2, tmp3, tmp5 = area_ref, max_area, length, x0
    else:
        nose_area = float(tbfunx(x, s, nose_length, lower=0, upper=0)[0])
        tmp1 = tmp2 = nose_area
        tmp3 = tmp5 = nose_length

    # Effective fineness ratio feeding the apparent-mass factor.
    fineness_effective = tmp3 / np.sqrt(tmp2 * 4.0 / PI)          # TMP4
    apparent_mass, _ = tbfunx(_FIG_42110_20_X, _FIG_42110_20_Y,
                              fineness_effective, lower=2, upper=1)  # BD(9)
    cla = 2.0 * apparent_mass * tmp1 / (RAD * sref)               # BODY(101)

    # Wetted-perimeter integral over the equally spaced grid.
    # eqspc1 names its resampled array 'se' whatever was passed in.
    perimeter = eqspc1(x, p, _STATIONS)
    perimeter_integral = float(trapz(perimeter['se'], perimeter['xe'])[0])

    # The source temporarily substitutes (TMP5, TMP1) at the first station
    # at or beyond TMP5, resamples, then restores the originals.
    splice_station = (int(np.argmax(x >= tmp5)) if np.any(x >= tmp5)
                      else len(x) - 1)
    x_work, s_work = x.copy(), s.copy()
    x_work[splice_station], s_work[splice_station] = tmp5, tmp1
    # IL is the last station up to that point where the area still changes.
    last_changing = 1
    for station in range(1, splice_station + 1):
        if s_work[station] - s_work[station - 1] != 0.0:
            last_changing = station
    count = max(last_changing + 1, 2)
    area = eqspc1(x_work[:count], s_work[:count], _STATIONS)

    # CMA: the first moment of the area slope about the reference station.
    const = 2.0 * apparent_mass / (RAD * sref * cbar)
    moment_integral = float(trapz(area['dsedx'] * area['xe'], area['xe'])[0])
    cma = (xcg * cla / cbar) - const * moment_integral            # BODY(121)

    # Planform integrals aft of the reference station, for crossflow.  The
    # source calls EQSPC1(X(L),R(L),...) while X(L) still holds the
    # substituted TMP5 (it is restored only after the drag buildup), so the
    # integral starts at TMP5, with the station's own radius.
    planform = eqspc1(x_work[splice_station:], r[splice_station:], _STATIONS)
    planform_area = float(trapz(planform['se'], planform['xe'])[0])   # BD(88)
    planform_moment = float(
        trapz(planform['se'] * planform['xe'], planform['xe'])[0])    # RXDFI

    # Skin friction with the source roughness cutoff.
    reynolds = reynolds_per_length * length                        # BD(90)
    cept, _ = tbfunx(_FIG_41510_27_MACH, _FIG_41510_27_CEPT,
                     mach, lower=0, upper=0)
    cutoff = roughness_length**1.0482 * 10.0**cept                 # BD(91)
    reynolds_used = min(reynolds, cutoff)
    friction_mach = 0.60 if transonic else mach
    cf = fig26(reynolds_used, friction_mach)                        # BD(92)

    base_diameter = np.sqrt(base_area * 4.0 / PI)                   # BD(86)
    max_diameter = np.sqrt(max_area * 4.0 / PI)                     # BD(85)
    fineness = length / max_diameter                                # BD(75)

    # Friction drag on the wetted area, with the form factor, then base drag.
    cd_friction = (cf * (1.0 + 60.0 / fineness**3 + 0.0025 * fineness) *
                   perimeter_integral / max_area)                   # BD(59)
    cd_base = 0.029 * ((base_diameter / max_diameter)**3 /
                       np.sqrt(cd_friction) * max_area / sref)      # BD(60)
    cd_friction = cd_friction * max_area / sref
    cd_zero_lift = cd_friction + cd_base                            # BD(61)

    result = {
        'cla': float(cla),
        'cma': float(cma),
        'base_area': float(base_area),                              # BD(57)
        'cd_zero_lift': float(cd_zero_lift),
        'cd_friction': float(cd_friction),
        'cd_base': float(cd_base),
        'cf': float(cf),
        'reynolds_used': float(reynolds_used),
        'roughness_cutoff_reynolds': float(cutoff),
        'apparent_mass_factor': float(apparent_mass),
        'fineness_ratio': float(fineness),
        'fineness_effective': float(fineness_effective),
        'planform_area': planform_area,
        'planform_moment': planform_moment,
        'perimeter_integral': perimeter_integral,
        'reference_station': float(x1),
        'body_contracts': contracts,
        'method': 'legacy_bodyrt',
    }
    if transonic:
        # The source returns here for the transonic pass.
        return result

    # --- Angle loop (source label 1100) ---
    angles = np.atleast_1d(np.asarray(alpha_deg, dtype=float))
    cn, cm, cd, cl, ca = (np.empty(len(angles)) for _ in range(5))
    cdc = np.empty(len(angles))

    # Figure 4.2.1.2-35A depends only on geometry, so it is read once.
    crossflow_drag, _ = tbfunx(
        _FIG_42120_35A_X, _FIG_42120_35A_Y, fineness, lower=2, upper=2,
    )

    for angle_index, angle in enumerate(angles):
        sin_a = np.sin(angle / RAD)
        sin_a_squared = sin_a ** 2
        cos_a = np.cos(angle / RAD)

        cn_potential = cla * angle                                  # BD(J+154)
        eta, _ = tbfunx(
            _FIG_42120_35B_X, _FIG_42120_35B_Y,
            mach * abs(sin_a), lower=0, upper=0,
        )
        cdc[angle_index] = eta
        sign = 1.0 if angle >= 0.0 else -1.0

        cn_viscous = (
            2.0 * sin_a_squared * crossflow_drag * eta
            * planform_area / sref * sign
        )
        cn[angle_index] = cn_potential + cn_viscous
        cm[angle_index] = (
            cma * angle
            - 2.0 * sin_a_squared * eta * crossflow_drag
            * (planform_moment - xcg * planform_area)
            / (cbar * sref) * sign
        )
        cd[angle_index] = cd_zero_lift + (cn_potential + cn_viscous) * sin_a

        # Source labels 1100: the rotation as written.
        cl[angle_index] = cn[angle_index] * cos_a + cd[angle_index] * sin_a
        ca[angle_index] = cd[angle_index] * cos_a - cn[angle_index] * sin_a

    result.update({
        'alpha_deg': angles,
        'cn': cn,
        'cm': cm,
        'cd': cd,
        'cl': cl,
        'ca': ca,
        'crossflow_drag_factor': float(crossflow_drag),
        # Under the source's own names, which BODYJM uses: ETA is BD(76),
        # read above as crossflow_drag, and CDC is BD(J+134), read above
        # per angle as eta.
        'bd76_eta': float(crossflow_drag),
        'bd135_cdc': cdc,
    })
    return result


def calculate_bodyjm(x: Sequence[float], s: Sequence[float],
                     alpha_deg: Sequence[float], bodyrt: Dict[str, object],
                     body_length: float, xcg: float, sref: float,
                     cbar: float, ellipticity: float = UNUSED
                     ) -> Dict[str, object]:
    """Translate BODYJM: body forces by Jorgensen's method.

    ``CN = (A3*F(a)*CN/CNS + SP*eta*cdc*sin|sin|*CN/CNN)/SREF`` with
    ``F(a) = sin(2a)*cos(a/2)``, the axial force ``CD0*cos(a)**2``, and CL
    and CD by rotation.  A circular section (``ELLIP`` unset or one) takes
    its potential terms from BODYRT's slopes; an elliptic one from the
    volume and base area, scaled by the Jorgensen cross-section ratios.

    Args:
        x, s: Body stations and areas; ``REQ = sqrt(S/pi)``.
        alpha_deg: ``FLC(23)`` onward.
        bodyrt: BODYRT's result, for ``cla`` and ``cma`` (``BODY(101)``,
            ``BODY(121)``), ``base_area`` ``BD(57)``, ``cd_zero_lift``
            ``BD(61)``, and ``bd76_eta`` and ``bd135_cdc``.
        body_length: ``BD(1)``.
        xcg: ``XCG`` from ``/SYNTSS/``; BODYRT uses ``BD(33)``.
        sref, cbar: ``SREF``, ``CBARR``.
        ellipticity: ``ELLIP``; ``UNUSED`` or less is taken as one, and
            written back as such.

    Returns:
        Dictionary with ``cd``, ``cl``, ``cm``, ``cn``, ``ca`` (the
        ``BODY`` slots it overwrites), their potential and viscous parts,
        ``volume`` ``BD(34)``, ``centroid`` ``BD(35)``, ``planform_area``
        ``BD(275)`` and the ``ellipticity`` left in COMMON.
    """
    x = np.asarray(x, dtype=float)
    req = np.sqrt(np.asarray(s, dtype=float) / PI)

    volume = trapz(req, x, -1)[0]
    moment = trapz(req * x, x, 1)[0]
    planform = 2.0 * trapz(req, x, 1)[0]
    centroid = 2.0 * moment / planform

    base = float(bodyrt['base_area'])
    a1 = (volume - base * (body_length - xcg)) / (sref * cbar)
    a2 = (planform * (xcg - centroid)) / (sref * cbar)
    a3 = base

    ellip = float(ellipticity)
    if ellip <= UNUSED:
        ellip = 1.0
    if ellip == 1.0:
        a1 = float(bodyrt['cma']) * RAD / 2.0
        a3 = float(bodyrt['cla']) * sref * RAD / 2.0

    aob = 1.0 / ellip if ellip < 1.0 else ellip
    cnocns = aob if ellip < 1.0 else 1.0 / aob
    cnocnn = 1.0
    if ellip < 1.0:
        e = 1.0 - 1.0 / aob ** 2
        cnocnn = 1.5 * np.sqrt(aob) * (
            -1.0 / aob ** 2 / e ** 1.5 * np.log(aob * (1.0 + np.sqrt(e)))
            + 1.0 / e
        )
    elif ellip > 1.0:
        cnocnn = 1.5 * np.sqrt(1.0 / aob) * (
            aob ** 2 / (aob ** 2 - 1.0) ** 1.5
            * np.arctan(np.sqrt(aob ** 2 - 1.0))
            - 1.0 / (aob ** 2 - 1.0)
        )

    eta = float(bodyrt['bd76_eta'])
    cdc = np.asarray(bodyrt['bd135_cdc'], dtype=float)
    cd0 = float(bodyrt['cd_zero_lift'])
    alp = np.asarray(alpha_deg, dtype=float) / RAD
    sin_alpha, cos_alpha = np.sin(alp), np.cos(alp)
    jorgensen_f = np.sin(2.0 * alp) * np.cos(alp / 2.0)

    cn_pot = a3 * jorgensen_f * cnocns / sref
    cn_vis = (
        planform * eta * cdc * sin_alpha * np.abs(sin_alpha) * cnocnn / sref
    )
    cm_pot = a1 * jorgensen_f * cnocns
    cm_vis = a2 * eta * cdc * sin_alpha * np.abs(sin_alpha) * cnocnn
    cn = cn_pot + cn_vis
    cm = cm_pot + cm_vis
    ca = cd0 * cos_alpha ** 2

    return {
        'cn': cn,
        'cm': cm,
        'ca': ca,
        'cl': cn * cos_alpha - ca * sin_alpha,
        'cd': ca * cos_alpha + cn * sin_alpha,
        'cn_potential': cn_pot,
        'cn_viscous': cn_vis,
        'cm_potential': cm_pot,
        'cm_viscous': cm_vis,
        'volume': float(volume),
        'centroid': float(centroid),
        'planform_area': float(planform),
        'ellipticity': ellip,
        'method': 'legacy_bodyjm',
    }


def calculate_m06o06(x: Sequence[float], s: Sequence[float],
                     p: Sequence[float], r: Sequence[float],
                     alpha_deg: Sequence[float], alpha_zero: float,
                     mach: float, reynolds_per_length: float, sref: float,
                     cbar: float, blref: float, xcg_bd33: float,
                     xcg: float, roughness: float = 1.6e-4,
                     method: float = 1.0, ellipticity: float = UNUSED,
                     experimental: bool = False) -> Dict[str, object]:
    """Translate overlay M06O06: the subsonic axisymmetric body.

    BODYRT, then BODYJM for ``METHOD`` above 1.5, then the slope pass:
    CLa and CMa by TBFUNX over the schedule from the second angle (every
    angle with experimental data), ``CY_b = -CLa``,
    ``Cn_b = -(CBARR/BLREF)*CMa``, ``Cl_b = 0``, and CN and CA by rotation.

    Args:
        x, s, p, r: ``/BODYI/`` stations.
        alpha_deg: ``FLC(23)`` onward.
        alpha_zero: ``BD(81)``; BODYRT works at ``BD(J+254) = FLC+BD(81)``.
        mach: ``B(1)``.
        reynolds_per_length: ``FLC(M+42)``.
        sref, cbar, blref: ``/OPTION/``.
        xcg_bd33: ``BD(33)``, BODYRT's moment reference.
        xcg: ``XCG``, BODYJM's.
        roughness: ``ROUGFC``.
        method: ``METHOD``, ``BODYIN(128)``.
        ellipticity: ``ELLIP``, ``BODYIN(129)``.
        experimental: ``KBODY``, which also retakes CM0 at zero lift.

    Returns:
        Dictionary with ``cd``, ``cl``, ``cm``, ``cn``, ``ca``, ``cla``,
        ``cma``, ``cyb``, ``cnb``, ``clb`` (the ``BODY`` block), and the
        ``bodyrt`` and ``bodyjm`` results.

    Notes:
        BODYRT stores its normal force in ``BODY(21)`` onward, the slot the
        rest of the program reads as lift, and this overlay then rotates it
        as lift.  With BODYJM the slot does hold lift.  Both are as
        executed.
    """
    alpha = np.asarray(alpha_deg, dtype=float)
    rt = calculate_bodyrt(x, s, p, r, alpha + alpha_zero, mach,
                          reynolds_per_length, sref, cbar, xcg_bd33,
                          roughness)
    cd, cl, cm = (np.array(rt[k], dtype=float) for k in ('cd', 'cn', 'cm'))
    jm = None
    if method > 1.5:
        jm = calculate_bodyjm(x, s, alpha, rt, float(np.asarray(x)[-1]), xcg,
                              sref, cbar, ellipticity)
        cd, cl, cm = jm['cd'].copy(), jm['cl'].copy(), jm['cm'].copy()
    result = body_slope_pass(alpha, cd, cl, cm, rt['cla'], rt['cma'], cbar,
                             blref, experimental)
    result.update({'bodyrt': rt, 'bodyjm': jm, 'method': 'legacy_m06o06'})
    return result


def body_slope_pass(alpha_deg: Sequence[float], cd: Sequence[float],
                    cl: Sequence[float], cm: Sequence[float], cla0: float,
                    cma0: float, cbar: float, blref: float,
                    experimental: bool = False) -> Dict[str, object]:
    """The closing pass M04O04 and M06O06 share.

    CLa and CMa by TBFUNX over the schedule from the second angle (every
    angle with experimental data, which also retakes CM0 at zero lift),
    the first keeping the method's analytic slopes; ``CY_b = -CLa``,
    ``Cn_b = -(CBARR/BLREF)*CMa``, ``Cl_b = 0``; CN and CA by rotation.
    """
    alpha = np.asarray(alpha_deg, dtype=float)
    cd, cl, cm = (np.asarray(v, dtype=float) for v in (cd, cl, cm))
    cla = np.zeros(len(alpha))
    cma = np.zeros(len(alpha))
    cla[0], cma[0] = cla0, cma0
    cm0 = None
    if experimental:
        cm0 = tbfunx(cl, cm, 0.0, 1, 1, ordered=False)[0]
    for angle_index in range(len(alpha)):
        if angle_index == 0 and not experimental:
            continue
        cla[angle_index] = tbfunx(alpha, cl, alpha[angle_index], 0, 0)[1]
        cma[angle_index] = tbfunx(alpha, cm, alpha[angle_index], 0, 0)[1]

    cos_alpha = np.cos(alpha / RAD)
    sin_alpha = np.sin(alpha / RAD)
    result = {
        'cd': cd,
        'cl': cl,
        'cm': cm,
        'cn': cl * cos_alpha + cd * sin_alpha,
        'ca': cd * cos_alpha - cl * sin_alpha,
        'cla': cla,
        'cma': cma,
        'cyb': -cla,
        'cnb': -(cbar / blref) * cma,
        'clb': np.zeros(len(alpha)),
    }
    if cm0 is not None:
        result['cm0'] = float(cm0)
    return result
