"""
Subsonic downwash at the horizontal tail for PyDATCOM.

Translates the geometry setup and the angle-of-attack downwash routine that
feed the wing-body-tail buildup:

- ``INFTGM``: downwash synthesizing dimensions A(193), A(194), A(24), A(12),
  A(11) and the tail-height/tail-arm pair they define.
- ``DWASH``: Figure 4.4.1-68A leading/trailing-edge separation criterion,
  Figure 4.4.1-66 effective aspect ratio, the Figure 4.4.1-67 vortex
  downwash gradient, Figure 4.4.1-68B tail-position factor, and the
  ``TWASH=2`` DATCOM Section 4.4.1 low-speed gradient.

Reference: datcom-legacy/datcom_2000/inftgm.f, dwash.f
"""

import numpy as np
from typing import Dict, Optional, Sequence
import logging

from pydatcom.utils.constants import PI, RAD, UNUSED
from pydatcom.utils.legacy_numeric import tbfunx, trapz
from pydatcom.utils.legacy_tables import tlinex
from pydatcom.geometry.wing import calculate_straight_exposed_geometry

logger = logging.getLogger(__name__)

# Figure 4.4.1-68A: type of flow separation as a function of airfoil
# leading-edge sharpness parameter and wing sweep, subsonic speeds.
# dwash.f DATA X4157A / Y4157A.
_FIG_4417_68A_SWEEP = np.array([0., 19., 20., 25., 30., 35., 40., 45.,
                                50., 54.])
_FIG_4417_68A_DELTAY = np.array([1.5, 1.5, 1.5, 1.57, 1.68, 1.81, 2.0,
                                 2.25, 2.75, 3.0])

# Figure 4.4.1-68B: average downwash acting on an aft lifting surface at
# low speeds.  dwash.f DATA X157B1, X157B2, Y4157B.  The source call is
# TLINEX(X157B1,X157B2,Y4157B,4,14,AOBVO2,BHOBV,...), so X157B1 is the
# a/(b_v/2) axis and X157B2 is the b_h/b_v axis.  Y4157B is declared
# (14,4): 14 rows of b_h/b_v by 4 columns of a/(b_v/2), filled column by
# column by the DATA statement.
_FIG_4417_68B_AOBV = np.array([0.0, 0.2, 0.6, 1.0])
_FIG_4417_68B_BHBV = np.array([.2, .3, .4, .5, .6, .7, .8, .833, 1.0,
                               1.1, 1.2, 1.3, 1.4, 1.5])
_FIG_4417_68B_TABLE = np.array([
    [1., 1.03, 1.068, 1.11, 1.176, 1.26, 1.36, 1.34, 1.20, 1.11, 1.0,
     .88, .74, .6],
    [.96, .96, .965, .976, .98, 1.016, 1.036, 1.04, .946, .88, .8,
     .69, .58, .46],
    [.74, .73, .72, .71, .69, .67, .644, .638, .582, .56, .52, .48,
     .44, .4],
    [.5, .5, .5, .49, .476, .46, .448, .44, .4, .372, .352, .32, .3, .27],
]).T


def fig4417_68a(sweep_c4_deg: float) -> float:
    """Figure 4.4.1-68A: limiting DELTAY for leading-edge separation.

    Returns the tabulated sharpness parameter at the given quarter-chord
    sweep.  ``DWASH`` sets ``LESEPA`` when this value exceeds the section
    DELTAY, i.e. leading-edge separation is predicted.

    Args:
        sweep_c4_deg: Quarter-chord sweep of the exposed panel, degrees.

    Returns:
        Limiting DELTAY.
    """
    value, _ = tbfunx(_FIG_4417_68A_SWEEP, _FIG_4417_68A_DELTAY,
                      float(sweep_c4_deg), lower=0, upper=2)
    return float(value)


def fig4417_68b(a_over_bv_half: float, bh_over_bv: float) -> float:
    """Figure 4.4.1-68B: average downwash factor on the aft surface.

    Args:
        a_over_bv_half: ``2*a/b_v``, the tail height above the vortex sheet
            divided by the vortex semispan.
        bh_over_bv: ``b_h/b_v``, tail span over vortex span.

    Returns:
        DEBODE, clamped at zero as ``DWASH`` does.
    """
    value = tlinex(_FIG_4417_68B_AOBV, _FIG_4417_68B_BHBV,
                   _FIG_4417_68B_TABLE,
                   float(a_over_bv_half), float(bh_over_bv),
                   1, 2, 1, 2)
    return float(max(value, 0.0))


def calculate_downwash_geometry(state: Dict) -> Dict[str, float]:
    """Translate INFTGM's downwash synthesizing dimensions.

    Computes the streamwise tail arm ``A(24)`` and the tail height above the
    extended wing chord plane ``A(12)``, including the hinge-axis rotation of
    the tail reference point and both surface incidences.

    Args:
        state: State dictionary with wing, htail and SYNTHS entries.

    Returns:
        Dictionary with ``tail_arm`` (A(24)), ``tail_height`` (A(12)),
        ``tail_angle`` (A(11), radians), and the A(193)/A(194) intermediates.

    Raises:
        ValueError: If either surface lacks complete straight-taper geometry.
    """
    wing = calculate_straight_exposed_geometry(state, component='wing')
    tail = calculate_straight_exposed_geometry(state, component='htail')

    chrdr = wing['theoretical_root_chord']
    aliw = np.deg2rad(float(state.get('synths_aliw', 0.0) or 0.0))
    alih = np.deg2rad(float(state.get('synths_alih', 0.0) or 0.0))
    xw = float(state.get('synths_xw', 0.0) or 0.0)
    zw = float(state.get('synths_zw', 0.0) or 0.0)
    xh = float(state.get('synths_xh', 0.0) or 0.0)
    zh = float(state.get('synths_zh', 0.0) or 0.0)

    # INFTGM rotates the tail reference point about the hinge axis first.
    hinax = state.get('synths_hinax')
    if hinax is not None:
        hinax = float(hinax)
        zh = zh + np.sin(alih) * (hinax - xh)
        xh = hinax * (1.0 - np.cos(alih)) + xh * np.cos(alih)

    # XBRSTH = ATH(30)-ATH(16)/4 is the exposed MAC quarter chord aft of the
    # exposed root leading edge; DXSTAR carries it to the centerline.
    xbrsth = tail['mac_c4_location']
    dxstar = (tail['theoretical_semispan'] - tail['semispan']) * tail['tan_le']
    dxbh = xbrsth + dxstar

    a193 = xh - xw - chrdr * np.cos(aliw)
    a194 = a193 + dxbh * np.cos(alih)
    zph = zh - dxbh * np.sin(alih) - zw + chrdr * np.sin(aliw)
    dlh = zph * np.tan(aliw)
    tail_arm = (a194 - dlh) * np.cos(aliw)
    tail_height = zph / np.cos(aliw) + (a194 - dlh) * np.sin(aliw)
    tail_angle = np.arctan2(tail_height, tail_arm) if tail_arm != 0.0 else 0.0

    return {
        'tail_arm': float(tail_arm),
        'tail_height': float(tail_height),
        'tail_angle': float(tail_angle),
        'a193': float(a193),
        'a194': float(a194),
        'zph': float(zph),
        'dxbh': float(dxbh),
        'xh_rotated': float(xh),
        'zh_rotated': float(zh),
    }


def calculate_downwash_gradient_441(state: Dict) -> Dict[str, float]:
    """Translate DWASH's TWASH=2 branch: DATCOM Section 4.4.1 gradient.

    ``de/da = 4.44 * (K_A * K_lambda * K_H * sqrt(cos(sweep_c4)))**1.19``

    with ``K_A = 1/A - 1/(1+A**1.7)``, ``K_lambda = (10-3*taper)/7`` and
    ``K_H = (1-|h_H/b|)/(2*l_H/b)**(1/3)``, all on the exposed panel.  This
    gradient is independent of angle of attack.

    The tail arm follows the source expression at dwash.f:
    ``XLH=(XH-AHT(161)*COS(ALIH))-(XW-A(161)*COS(ALIW))``.  Both bundled
    sources agree on it.  Note that it *subtracts* each surface's
    theoretical MAC quarter-chord offset where the published Section 4.4.1
    method adds them; the source form is preserved and the offsets are
    returned so a caller can inspect the difference.

    Args:
        state: State dictionary with wing, htail and SYNTHS entries.

    Returns:
        Dictionary with ``deda`` and the K-factors and lengths behind it.

    Raises:
        ValueError: If the geometry is incomplete or gives a nonpositive arm.
    """
    wing = calculate_straight_exposed_geometry(state, component='wing')
    tail = calculate_straight_exposed_geometry(state, component='htail')
    geometry = calculate_downwash_geometry(state)

    aspect_ratio = wing['aspect_ratio']
    taper = wing['taper_ratio']
    sspne = wing['semispan']
    if aspect_ratio <= 0.0 or sspne <= 0.0:
        raise ValueError("Section 4.4.1 downwash requires positive exposed "
                         "aspect ratio and semispan")

    aliw = np.deg2rad(float(state.get('synths_aliw', 0.0) or 0.0))
    alih = np.deg2rad(float(state.get('synths_alih', 0.0) or 0.0))
    xw = float(state.get('synths_xw', 0.0) or 0.0)
    xh = float(state.get('synths_xh', 0.0) or 0.0)

    xlh = ((xh - tail['mac_c4_theoretical'] * np.cos(alih)) -
           (xw - wing['mac_c4_theoretical'] * np.cos(aliw)))
    if xlh <= 0.0:
        raise ValueError("Section 4.4.1 downwash requires a positive tail "
                         f"arm; source XLH expression gave {xlh:g}")

    cos_c4 = 1.0 / np.sqrt(1.0 + wing['tan_c4']**2)
    xka = 1.0 / aspect_ratio - 1.0 / (1.0 + aspect_ratio**1.7)
    xkl = (10.0 - 3.0 * taper) / 7.0
    xkh = ((1.0 - abs(0.5 * geometry['tail_height'] / sspne)) /
           (xlh / sspne)**(1.0 / 3.0))
    product = xka * xkl * xkh * np.sqrt(cos_c4)
    deda = 4.44 * product**1.19 if product > 0.0 else 0.0

    return {
        'deda': float(deda),
        'k_a': float(xka),
        'k_lambda': float(xkl),
        'k_h': float(xkh),
        'cos_sweep_c4': float(cos_c4),
        'tail_arm_xlh': float(xlh),
        'tail_height': geometry['tail_height'],
        'wing_mac_c4_offset': wing['mac_c4_theoretical'],
        'tail_mac_c4_offset': tail['mac_c4_theoretical'],
        'method': 'legacy_dwash_section_441',
    }


def calculate_dwash(alpha_deg: Sequence[float],
                    wing_alone: Dict[str, Sequence[float]],
                    wing: Dict[str, float],
                    wing_geometry: Dict[str, float],
                    synthesis: Dict[str, float],
                    tail: Dict[str, float],
                    tail_geometry: Dict[str, float],
                    sref: float,
                    kwb: float) -> Dict[str, object]:
    """Translate DWASH: subsonic downwash at the horizontal tail.

    For each angle in the schedule the routine builds a local grid of ``NA``
    angles from the wing zero-lift angle up to the local wing angle, forms
    the downwash gradient at each, and integrates it with TRAPZ into the
    downwash angle.  The gradient comes from one of two methods:

    - ``TWASH != 2``: the Figure 4.4.1-66 effective aspect ratio and the
      Figure 4.4.1-67 vortex gradient, with the integrated angle scaled by
      the Figure 4.4.1-68B average-downwash factor at the local angle.
    - ``TWASH = 2``: the Section 4.4.1 low-speed gradient of
      :func:`calculate_downwash_gradient_441`, constant in angle.

    Both branches compute the vortex height and span (``FACT(J+61)`` and
    ``FACT(J+81)``), since ``WBTAIL`` reads the span whichever method is
    selected.  Each is the value at the last grid point, the local angle.

    A final pass replaces each stored gradient with the TBFUNX slope of the
    downwash angle over the schedule, end modes 1 and 1.  The integrand
    value at the local angle is returned separately as ``gradient_local``.

    Args:
        alpha_deg: ``FLC(23)`` onward, the angle-of-attack schedule.
        wing_alone: ``{'alpha', 'cl'}``: the wing lift curve ``WING(21)``
            onward, on the SREF basis, over its local angles ``B(23)``
            onward (the schedule plus ALIW).
        wing: ``/WINGI/`` entries ``sspn``, ``sspnop``, ``sspndd``,
            ``chrdtp``, ``chrdr``, ``chrdbp``, ``dhdadi``, ``dhdado``,
            ``deltay``, ``twash``, plus ``sspne`` for the TWASH=2 branch.
        wing_geometry: ``A`` block entries: ``area`` A(3), ``aspect_ratio``
            A(7), ``taper_ratio`` A(27), ``taper_ratio_theoretical``
            A(118), ``sweep_c4_deg`` A(40), ``cos_c4`` A(43), ``tan_c4``
            A(44), ``tan_le`` A(62), ``mac_c4_theoretical`` A(161);
            ``alpha_zero_lift`` B(49), the start of every integration grid;
            ``alpha_zero_lift_reference`` A(126) and ``alpha_clmax_reference``
            A(127), the Mach-zero pair CLMCH0 stores, which scale the
            effective aspect ratio.
        synthesis: ``/SYNTSS/`` entries ``aliw``, ``xw``, ``xh``, ``alih``.
        tail: ``/HTI/`` entries; only ``sspn`` (HTIN(4)) is read.
        tail_geometry: ``tail_arm`` A(24), ``tail_height`` A(12) and
            ``a22`` A(22) from INFTGM, plus the tail's
            ``mac_c4_theoretical`` AHT(161).
        sref: ``SREF``.
        kwb: ``WB(2)``, the wing-in-presence-of-body factor K_W(B).

    Returns:
        Dictionary with ``angle`` (``DWASHI(21)`` onward, degrees),
        ``gradient`` (``DWASHI(41)`` onward after the TBFUNX pass),
        ``gradient_local``, ``vortex_height`` (``FACT(62)`` onward),
        ``vortex_span`` (``FACT(82)`` onward), ``debode``, ``a20`` (the
        effective tail arm DWASH writes back to A(20)) and
        ``leading_edge_separation``.

    Raises:
        ValueError: If the reference angles coincide, or the TWASH=2 tail
            arm is not positive.
    """
    alpha = np.asarray(alpha_deg, dtype=float)
    ang = np.asarray(wing_alone['alpha'], dtype=float)
    cl_table = np.asarray(wing_alone['cl'], dtype=float)
    nalpha = len(alpha)

    idwash = int(float(wing.get('twash', 0.0)) + 0.5)
    sspn = float(wing['sspn'])
    sspnop = float(wing.get('sspnop', 0.0))
    sspndd = float(wing.get('sspndd', 0.0))
    chrdtp = float(wing['chrdtp'])
    chrdr = float(wing['chrdr'])
    chrdbp = float(wing.get('chrdbp', chrdtp))
    dhdadi = float(wing.get('dhdadi', 0.0))
    dhdado = float(wing.get('dhdado', 0.0))

    area = float(wing_geometry['area'])
    aspect_ratio = float(wing_geometry['aspect_ratio'])
    taper = float(wing_geometry['taper_ratio'])
    cos_c4 = float(wing_geometry['cos_c4'])
    tan_c4 = float(wing_geometry['tan_c4'])
    tan_le = float(wing_geometry['tan_le'])
    alpha_zero = float(wing_geometry['alpha_zero_lift'])
    alpha_zero_ref = float(wing_geometry.get('alpha_zero_lift_reference',
                                             alpha_zero))
    alpha_clmax_ref = float(wing_geometry['alpha_clmax_reference'])
    if alpha_clmax_ref == alpha_zero_ref:
        raise ValueError("DWASH divides by A(127)-A(126), the Mach-zero "
                         "stall and zero-lift angles, which coincide")

    tail_arm = float(tail_geometry['tail_arm'])
    tail_height = float(tail_geometry['tail_height'])
    a22 = float(tail_geometry['a22'])
    corr = float(synthesis.get('aliw', 0.0))

    # Figure 4.4.1-68A at the quarter-chord sweep.
    lesepa = fig4417_68a(wing_geometry['sweep_c4_deg']) > float(wing['deltay'])
    bvrutp = (0.78 + 0.10 * (float(wing_geometry['taper_ratio_theoretical'])
                             - 0.40) + 0.003 * float(wing_geometry['sweep_c4_deg']))
    sratio = sref / area
    tl2ob = tail_arm / sspn
    lex = -1 if alpha_zero_ref < ang[0] else 1

    # The TWASH=2 gradient is independent of angle, so form it once.
    if idwash == 2:
        aliw_rad = corr / RAD
        alih_rad = float(synthesis.get('alih', 0.0)) / RAD
        sspne = float(wing['sspne'])
        xlh = ((float(synthesis['xh']) -
                float(tail_geometry['mac_c4_theoretical']) * np.cos(alih_rad)) -
               (float(synthesis['xw']) -
                float(wing_geometry['mac_c4_theoretical']) * np.cos(aliw_rad)))
        if xlh <= 0.0:
            raise ValueError("DWASH's Section 4.4.1 branch requires a "
                             f"positive tail arm; XLH is {xlh:g}")
        xka = 1.0 / aspect_ratio - 1.0 / (1.0 + aspect_ratio**1.7)
        xkl = (10.0 - 3.0 * taper) / 7.0
        xkh = ((1.0 - abs(0.5 * tail_height / sspne)) /
               (xlh / sspne)**(1.0 / 3.0))
        deda_441 = 4.44 * (xka * xkl * xkh * np.sqrt(cos_c4))**1.19

    angle = np.zeros(nalpha)
    gradient_local = np.zeros(nalpha)
    vortex_height = np.zeros(nalpha)
    vortex_span = np.zeros(nalpha)
    debode_out = np.zeros(nalpha)
    a20 = 0.0
    for j in range(nalpha):
        alp = alpha[j] + corr
        na = min(max(int(abs(alp) + 1.5), 2), 21)
        xna = float(na - 1)
        grid = np.array([alpha_zero + k * (alp - alpha_zero) / xna
                         for k in range(na)])
        deda = np.zeros(na)
        for k in range(na):
            bj22 = alpha_zero_ref + k * (alp - alpha_zero_ref) / xna
            adoad = abs((bj22 - alpha_zero_ref) /
                        (alpha_clmax_ref - alpha_zero_ref))
            clwj, claw = tbfunx(ang, cl_table, bj22, lex, 1)
            if lex < 0 and bj22 < ang[0]:
                clwj = clwj * (bj22 - alpha_zero_ref) / (ang[0] - alpha_zero_ref)
            clwj = sratio * clwj * kwb

            # Figure 4.4.1-66: effective aspect ratio and span.
            aeefoa = (taper / 3.0 + (1.0 - taper / 3.0) *
                      (1.4 - (1.4 - cos_c4 + .04 * tan_c4) * adoad))
            if adoad < .4 / (1.4 - cos_c4 + .04 * tan_c4):
                aeefoa = 1.0
            aeff = aspect_ratio * aeefoa
            beffob = 2.0 * aeefoa / (1.0 + taper + aeefoa * (1.0 - taper))

            # Figure 4.4.1-67 with compressibility.
            dedai = 1.62 * claw / (PI * aeff) * RAD * sratio * kwb
            tzob = tl2ob + np.sqrt(0.5 * (-1.0 + np.sqrt(1.0 + 4.0 / aeff**2)))
            dedav = dedai + (1.0 - dedai) / (aeff * tzob * np.sqrt(1.0 + tzob**2))
            beff = 2.0 * sspn * beffob

            # Effective tail length.
            bdff = sspn - sspndd
            beffo2 = beff / 2.0
            if sspnop <= UNUSED:
                cteff = (chrdtp - chrdr) / sspn * beffo2 + chrdr
            else:
                cteff = (((chrdtp - chrdbp) / sspnop) *
                         (sspnop - (sspn - beffo2)) + chrdbp)
            a20 = tail_arm - (beffo2 * tan_le + cteff / 4.0) + chrdr
            wfact = tail_arm + a22 if lesepa else a20
            bvru = bvrutp * beff
            drop = wfact * (bj22 / RAD - 0.41 * clwj / (PI * aeff))
            if beffo2 > bdff:
                height = tail_height - drop - (
                    bdff * np.tan(dhdadi / RAD) +
                    (beffo2 - bdff) * np.tan(dhdado / RAD))
            else:
                height = tail_height - drop - .5 * beff * np.tan(dhdadi / RAD)
            if clwj == 0.0:
                span = beff
            else:
                eru = 0.56 * aspect_ratio / clwj
                span = beff - (beff - bvru) * np.sqrt(abs(a20 / (sspn * eru)))

            # Figure 4.4.1-68B.
            debode = fig4417_68b(abs(2.0 * height / span),
                                 2.0 * float(tail['sspn']) / span)
            deda[k] = deda_441 if idwash == 2 else dedav

        vortex_height[j] = height
        vortex_span[j] = span
        debode_out[j] = debode
        gradient_local[j] = deda[-1]
        angle[j] = trapz(deda, grid, 1)[0]
        if idwash != 2:
            angle[j] *= debode

    # Label 1060: the stored gradient becomes the slope of the angle.
    gradient = np.array([tbfunx(alpha, angle, a, 1, 1)[1] for a in alpha])

    return {
        'angle': angle,
        'gradient': gradient,
        'gradient_local': gradient_local,
        'vortex_height': vortex_height,
        'vortex_span': vortex_span,
        'debode': debode_out,
        'a20': float(a20),
        'leading_edge_separation': bool(lesepa),
        'method': ('legacy_dwash_section_441' if idwash == 2
                   else 'legacy_dwash_vortex'),
    }


def calculate_dyprls(state: Dict, alpha_deg: float, cd0_wing: float,
                     cl_wing: Optional[float] = None,
                     eps_rad: Optional[float] = None) -> Dict[str, float]:
    """Translate DYPRLS: dynamic-pressure loss in the wing wake.

    The wake has a cosine-squared velocity profile of half-width ``ZWOCB``
    and centerline loss ``DQOQ0``, both in MAC units:

    ``q/q_inf = 1 - DQOQ0 * cos(pi/2 * ZOCB/ZWOCB)**2``

    and no loss at all once the surface lies outside the wake.  The source
    takes the wake deflection ``EJ`` either from the wing load, as
    ``1.62*CL/(pi*A)``, or from the translated downwash when ``KEPSLN`` is
    set; supplying ``eps_rad`` selects the latter.

    Args:
        state: State dictionary with wing and htail geometry.
        alpha_deg: Angle of attack, degrees.
        cd0_wing: ``CDOW``, wing zero-lift drag coefficient.
        cl_wing: Wing lift coefficient on the aircraft SREF basis; required
            unless ``eps_rad`` is given.
        eps_rad: Downwash angle in radians, selecting the ``KEPSLN`` branch.

    Returns:
        Dictionary with ``qoqi``, the wake half-width, the centerline loss
        and the surface offset, all in MAC units.

    Raises:
        ValueError: If the geometry or the drag input is unusable.
    """
    if cd0_wing < 0.0:
        raise ValueError("DYPRLS requires a nonnegative wing zero-lift drag")
    if cl_wing is None and eps_rad is None:
        raise ValueError("DYPRLS needs either a wing CL or a downwash angle")

    wing = calculate_straight_exposed_geometry(state, component='wing')
    geometry = calculate_downwash_geometry(state)
    area = wing['area']
    mac = wing['mac']
    sref = float(state.get('options_sref', area) or area)
    if min(area, mac, sref) <= 0.0:
        raise ValueError("DYPRLS requires positive exposed area, MAC and SREF")

    # GAMMA is A(11), the inclination of the line from the wing to the
    # tail, and I2 is A(24), INFTGM's streamwise tail arm, exactly as
    # M09O11 passes them.
    curve = calculate_dyprls_curve(
        cd0_wing, geometry['tail_arm'], mac, [0.0 if cl_wing is None
                                             else float(cl_wing)],
        wing['aspect_ratio'], geometry['tail_angle'], [float(alpha_deg)],
        sref, area,
        None if eps_rad is None else [float(eps_rad) * RAD])
    return {key: (value[0] if isinstance(value, np.ndarray) else value)
            for key, value in curve.items()}


def calculate_dyprls_curve(cdow: float, i2: float, cbar: float,
                           cl_wing: Sequence[float], aspect_ratio: float,
                           gamma: float, alpha_deg: Sequence[float],
                           sref: float, area: float,
                           downwash_deg: Optional[Sequence[float]] = None
                           ) -> Dict[str, object]:
    """Translate DYPRLS with its source arguments, at every angle.

    ``q/q_inf = 1 - DQOQ0 * cos(pi/2 * ZOCB/ZWOCB)**2`` inside the wake and
    1 outside it, with the wake deflection from the wing load,
    ``1.62*CL/(pi*A)`` on the exposed-area basis, or from DWASH when
    ``KEPSLN`` is set.

    Args:
        cdow: ``CDOW``, the wing zero-lift drag (M09O11 passes ``B(46)``).
        i2: ``I2``, M09O11's ``A(24)``: INFTGM's streamwise tail arm.  The
            routine forms the distance along the line to the tail itself,
            through ``cos(GAMMA)``.
        cbar: ``CBAR``, ``A(16)``, the exposed MAC.
        cl_wing: ``CLJW``, ``WING(21)`` onward, on the SREF basis.
        aspect_ratio: ``AW``, ``A(7)``.
        gamma: ``GAMMA``, ``A(11)``, radians.
        alpha_deg: ``ALPHA``, ``B(23)`` onward: the wing's local angles.
        sref: ``SREF``.
        area: ``A(3)``, the exposed area.
        downwash_deg: ``DWASH(21)`` onward, which selects the ``KEPSLN``
            branch.

    Returns:
        Dictionary of arrays: ``qoqi`` (``DWASH(1)`` onward) and the wake
        quantities behind it.

    Notes:
        With zero drag the wake has no width and the source divides by
        zero; outside or at its edge the result is 1, as here.
    """
    alpha = np.asarray(alpha_deg, dtype=float)
    cl = np.asarray(cl_wing, dtype=float)
    fact = 1.62 / (PI * aspect_ratio)
    if downwash_deg is not None:
        ej = np.asarray(downwash_deg, dtype=float) / RAD
    else:
        ej = fact * cl * sref / area
    angle = alpha / RAD
    i2ocb = i2 * np.cos(gamma - angle + ej) / (np.cos(gamma) * cbar)
    zwocb = 0.68 * np.sqrt(cdow * (i2ocb + 0.15) * sref / area)
    dqoq0 = 2.42 * np.sqrt(cdow * sref / area) / (i2ocb + 0.3)
    zocb = i2ocb * np.tan(ej + gamma - angle)
    with np.errstate(divide='ignore', invalid='ignore'):
        ratio = zocb / zwocb
        inside = np.abs(ratio) < 1.0
        qoqi = np.where(inside, 1.0 - dqoq0 * np.cos(0.5 * PI * ratio)**2,
                        1.0)
    return {
        'qoqi': qoqi, 'wake_half_width': zwocb, 'centerline_loss': dqoq0,
        'surface_offset': zocb, 'streamwise_distance': i2ocb,
        'wake_deflection': ej, 'in_wake': inside,
        'method': 'legacy_dyprls',
    }


def calculate_downwash(state: Dict, alpha_deg: float,
                       cl_wing: Optional[float] = None) -> Dict[str, float]:
    """Downwash angle and gradient at the horizontal tail.

    Uses the translated Section 4.4.1 gradient.  ``DWASH`` integrates the
    gradient from the wing zero-lift angle to the local angle with TRAPZ, so
    for the constant 4.4.1 gradient the downwash angle is
    ``eps = deda * (alpha - alpha_0L)``.

    Args:
        state: State dictionary.
        alpha_deg: Angle of attack, degrees.
        cl_wing: Unused by the 4.4.1 branch; accepted so callers can pass the
            wing load that the vortex branch requires.

    Returns:
        Dictionary with ``eps_deg``, ``deda`` and the dynamic-pressure ratio.
    """
    gradient = calculate_downwash_gradient_441(state)
    alpha_zero = float(state.get('wing_alpha_zero_lift', 0.0) or 0.0)
    eps_deg = gradient['deda'] * (float(alpha_deg) - alpha_zero)

    # QOQI: an explicit input wins; otherwise use the translated DYPRLS wake
    # model when a wing zero-lift drag is available, and fall back to no
    # loss when it is not.
    supplied = state.get('htail_qoqi')
    wake = None
    if supplied is not None:
        qoqi = float(supplied)
        qoqi_method = 'supplied'
    else:
        cd0_wing = state.get('wing_cdo')
        if cd0_wing is None:
            qoqi = 1.0
            qoqi_method = 'no_loss_default'
        else:
            wake = calculate_dyprls(state, alpha_deg, float(cd0_wing),
                                    cl_wing=cl_wing,
                                    eps_rad=np.deg2rad(eps_deg))
            qoqi = wake['qoqi']
            qoqi_method = wake['method']

    result = dict(gradient)
    result.update({
        'eps_deg': float(eps_deg),
        'alpha_deg': float(alpha_deg),
        'alpha_zero_lift': alpha_zero,
        'qoqi': float(qoqi),
        'qoqi_method': qoqi_method,
    })
    if wake is not None:
        result['wake'] = wake
    return result
