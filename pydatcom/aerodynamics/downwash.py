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
from typing import Dict, Optional
import logging

from pydatcom.utils.legacy_numeric import tbfunx
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
    del cl_wing
    gradient = calculate_downwash_gradient_441(state)
    alpha_zero = float(state.get('wing_alpha_zero_lift', 0.0) or 0.0)
    eps_deg = gradient['deda'] * (float(alpha_deg) - alpha_zero)
    result = dict(gradient)
    result.update({
        'eps_deg': float(eps_deg),
        'alpha_deg': float(alpha_deg),
        'alpha_zero_lift': alpha_zero,
        'qoqi': float(state.get('htail_qoqi', 1.0) or 1.0),
    })
    return result
