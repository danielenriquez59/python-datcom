"""
Body vortex effects on a lifting surface, DATCOM Section 4.3.1.3.

Translates the routines that supply the vortex term ``WBCLB`` and ``CLWBT``
apply to the surface load at high angle of attack:

- ``GETMAX``: maximum cross-sectional area and its station.
- ``ALI``:    vortex interference factor by the circle-theorem image system.
- ``BODOWG``: Figures 4.3.1.3-13A, -13B, -14 and -15, giving the vortex
  separation station, its lateral and vertical position, and the vortex
  strength ``Gamma/(2*pi*a*V)``.

The source gates the whole effect on ``|alpha| >= 6`` degrees and on a
positive downstream distance, so the term vanishes for cruise angles.

Reference: datcom-legacy/datcom_2000/bodowg.f, ali.f, getmax.f, wbclb.f
"""

import numpy as np
from typing import Dict, Sequence, Tuple
import logging

from pydatcom.utils.legacy_numeric import tbfunx
from pydatcom.utils.constants import RAD

logger = logging.getLogger(__name__)

# Figure 4.3.1.3-13A: bodowg.f DATA X1313A / Y1313A.
_FIG_13A_ALPHA = np.array([6.8, 7.2, 7.5, 8., 8.4, 9., 9.6, 10.4, 11.2,
                           12.1, 13.4, 15., 16., 17.1, 18., 20.])
_FIG_13A_XOR = np.array([20., 19., 18., 17., 16., 15., 14., 13., 12., 11.,
                         10., 9., 8.5, 8., 7.7, 7.])

# Figure 4.3.1.3-13B: bodowg.f DATA X1313B / Y1313B.
_FIG_13B_XD = np.array([0., .5, 1., 2., 2.5, 3.3, 4., 4.7, 5.5, 6., 7.])
_FIG_13B_ZOR = np.array([.86, 1.02, 1.2, 1.52, 1.65, 1.8, 1.92, 2., 2.09,
                         2.14, 2.23])

# Figure 4.3.1.3-14: bodowg.f DATA X31314 / Y31314.  The source writes the
# final four entries with the FORTRAN repeat count "4*.765".
_FIG_14_XD = np.array([0., .5, 1., 1.5, 2., 2.5, 3., 3.5, 4., 4.5, 5., 5.5,
                       6., 6.5, 7.])
_FIG_14_YOR = np.array([.5, .57, .62, .66, .69, .71, .72, .735, .75, .755,
                        .76, .765, .765, .765, .765])

# Figure 4.3.1.3-15: bodowg.f DATA X31315 / Y31315.
_FIG_15_XD = np.array([0., .5, 1., 1.5, 2., 2.5, 3., 3.5, 4., 4.5, 5., 5.5,
                       6., 6.5, 7.])
_FIG_15_GAMMA = np.array([.4, .5, .6, .7, .81, .91, 1.02, 1.13, 1.24, 1.35,
                          1.47, 1.6, 1.74, 1.88, 2.02])

# BODOWG skips the vortex system below this angle of attack.
_ALPHA_CUTOFF_DEG = 6.0


def getmax(x: Sequence[float], s: Sequence[float]) -> Tuple[float, float, int]:
    """Translate GETMAX: the maximum of S and the X and index at it.

    The source keeps the first occurrence of the maximum: a later station is
    adopted only when it is strictly greater.

    Args:
        x: Station coordinates.
        s: Values to maximize, normally cross-sectional area.

    Returns:
        ``(xmax, smax, index)`` with a zero-based index.

    Raises:
        ValueError: If the inputs are empty or mismatched.
    """
    x = np.asarray(x, dtype=float)
    s = np.asarray(s, dtype=float)
    if x.ndim != 1 or x.shape != s.shape or len(x) == 0:
        raise ValueError("GETMAX requires nonempty matching one-dimensional arrays")
    index = int(np.argmax(s))
    return float(x[index]), float(s[index]), index


def ali(q: float, r: float, xst: float, xrt: float, xslt: float) -> float:
    """Translate the ALI vortex interference factor.

    Sums four contributions: the vortex at ``(r, q)``, its mirror across the
    plane of symmetry, and the two circle-theorem images inside the body of
    radius ``xrt``.

    Args:
        q: Vortex height above the surface plane, length units.
        r: Vortex lateral position, length units.
        xst: Surface semispan.
        xrt: Body radius.
        xslt: Surface taper ratio.

    Returns:
        The interference factor ``ALI``.

    Notes:
        The source forms ``sign(FST)*H*atan(|FST|)``, which is undefined when
        ``FST`` is exactly zero.  That expression equals ``H*atan(FST)`` for
        every nonzero argument, so the identity is used instead; it removes
        the singular case without changing any valid source result.
    """
    if xst == xrt:
        raise ValueError("ALI requires a nonzero exposed semispan (XST != XRT)")
    if xslt == -1.0:
        raise ValueError("ALI is singular at a taper ratio of -1")

    ff = float(r)
    f = float(r)
    h = float(q)
    al = np.zeros(4)
    for mk in range(4):
        denominator = h**2 + (f - xrt)**2
        if denominator == 0.0:
            # The vortex sits exactly on the body surface.  The source
            # divides by zero here too; reject it rather than propagate an
            # infinity.  BODOWG never produces it: Figure 4.3.1.3-13B keeps
            # z/r at or above 0.86 and Figure 4.3.1.3-14 keeps y/r below 1.
            raise ValueError("ALI is singular for a vortex on the body surface")
        ab = np.log((h**2 + (f - xst)**2) / denominator)
        if h == 0.0:
            bb = xst - xrt
        else:
            bb = (h * np.arctan((f - xst) / h) -
                  h * np.arctan((f - xrt) / h) + (xst - xrt))
        term = ((xst - xrt * xslt) - f * (1.0 - xslt)) / (2.0 * (xst - xrt))
        al[mk] = term * ab - ((1.0 - xslt) / (xst - xrt)) * bb
        f = -f
        if mk == 1:
            # Circle-theorem image: reflect the pair inside the body radius.
            rt2 = xrt * xrt
            h2 = h * h
            f = (f * rt2) / (f**2 + h2)
            h = (h * rt2) / (ff**2 + h2)
    return float(2.0 / (1.0 + xslt) * (al[0] - al[1] - al[2] + al[3]))


def fig4313_13a(alpha_deg: float) -> float:
    """Figure 4.3.1.3-13A: vortex separation station ``x/r``."""
    value, _ = tbfunx(_FIG_13A_ALPHA, _FIG_13A_XOR, float(alpha_deg),
                      lower=0, upper=0)
    return float(value)


def fig4313_13b(xd: float) -> float:
    """Figure 4.3.1.3-13B: vortex height ``z/r``."""
    value, _ = tbfunx(_FIG_13B_XD, _FIG_13B_ZOR, float(xd), lower=0, upper=0)
    return float(value)


def fig4313_14(xd: float) -> float:
    """Figure 4.3.1.3-14: vortex lateral position ``y/r``."""
    value, _ = tbfunx(_FIG_14_XD, _FIG_14_YOR, float(xd), lower=0, upper=1)
    return float(value)


def fig4313_15(xd: float) -> float:
    """Figure 4.3.1.3-15: vortex strength ``Gamma/(2*pi*a*V)``."""
    value, _ = tbfunx(_FIG_15_XD, _FIG_15_GAMMA, float(xd), lower=0, upper=1)
    return float(value)


def calculate_bodowg(alpha_deg: float, x_quarter_chord: float,
                     body_radius: float, semispan: float,
                     taper_ratio: float) -> Dict[str, float]:
    """Translate BODOWG's per-angle body vortex effect on a surface.

    Args:
        alpha_deg: Angle of attack, degrees.  The source uses ``|alpha|``.
        x_quarter_chord: ``XCBO4``, the surface quarter-chord station
            measured aft of the body nose.
        body_radius: ``RCREO2``, from the maximum cross-sectional area as
            ``sqrt(Smax/pi)``.
        semispan: ``BWO2``, the surface semispan.
        taper_ratio: ``TRAT``, the surface taper ratio.

    Returns:
        Dictionary with ``ivbw`` (interference factor), ``go2pav`` (vortex
        strength), the intermediate ``xd``, and ``active``.

    Raises:
        ValueError: If the body radius is not positive.
    """
    if body_radius <= 0.0:
        raise ValueError("BODOWG requires a positive body radius")

    inactive = {
        'ivbw': 0.0, 'go2pav': 0.0, 'xd': 0.0, 'active': False,
        'method': 'legacy_bodowg',
    }
    magnitude = abs(float(alpha_deg))
    if magnitude < _ALPHA_CUTOFF_DEG:
        return inactive

    xrt = x_quarter_chord / body_radius
    xor = fig4313_13a(magnitude)
    xd = magnitude * (xrt - xor) / RAD
    if xd <= 0.0:
        return dict(inactive, xd=float(xd))

    zor = fig4313_13b(xd)
    yor = fig4313_14(xd)
    gamma = fig4313_15(xd)
    z0 = zor * body_radius
    y0 = yor * body_radius
    return {
        'ivbw': ali(z0, y0, semispan, body_radius, taper_ratio),
        'go2pav': float(gamma),
        'xd': float(xd),
        'z0': float(z0),
        'y0': float(y0),
        'active': True,
        'method': 'legacy_bodowg',
    }


def body_vortex_lift_increment(state: Dict, alpha_deg: float,
                               cla_surface: float,
                               component: str = 'htail') -> Dict[str, float]:
    """The ``FACT`` vortex term that WBCLB and CLWBT add to the surface load.

    ``WBCLB`` forms ``FACT(1)*FACT(J+1)*FACT(J+21)*CLAW*(ALP-ALIW)`` and
    applies it only when ``FACT(1) >= 1/3``; ``CLWBT`` applies the identical
    product to ``CLAH*ALPT``.  ``FACT(1)`` is the body-radius to semispan
    ratio and the other two factors are BODOWG's ``IVBW`` and ``GO2PAV``.

    Args:
        state: State dictionary with the surface planform.
        alpha_deg: Local surface flow angle, degrees.
        cla_surface: Surface lift-curve slope per degree, on the basis the
            caller wants the increment in.
        component: Surface key prefix.

    Returns:
        Dictionary with ``increment`` and the factors behind it.  The
        increment is zero whenever the source gate is closed.
    """
    from pydatcom.geometry.wing import calculate_straight_exposed_geometry
    from pydatcom.interactions.carryover import body_semispan_ratio

    ratio = body_semispan_ratio(state, component)
    geometry = calculate_straight_exposed_geometry(state, component=component)
    closed = {
        'increment': 0.0, 'ratio': ratio, 'ivbw': 0.0, 'go2pav': 0.0,
        'gated': True, 'method': 'legacy_wbclb_vortex_term',
    }
    # WBCLB skips the term entirely below a third; the body is too small to
    # shed a vortex system that reaches the surface.
    if ratio < 1.0 / 3.0:
        return closed

    body_radius = ratio * geometry['theoretical_semispan']
    origin = {'wing': 'synths_xw', 'htail': 'synths_xh'}.get(
        component, 'synths_xw')
    x_quarter_chord = (float(state.get(origin, 0.0) or 0.0) +
                       geometry['mac_c4_theoretical'])
    # WBCLB calls BODOWG with WINGIN(4), the theoretical semispan measured
    # from the centerline, because ALI integrates from the body radius out
    # to the tip.  The taper over that interval is the exposed taper ratio.
    vortex = calculate_bodowg(alpha_deg, x_quarter_chord, body_radius,
                              geometry['theoretical_semispan'],
                              geometry['taper_ratio'])
    if not vortex['active']:
        return dict(closed, gated=False, xd=vortex['xd'])
    increment = (ratio * vortex['ivbw'] * vortex['go2pav'] *
                 cla_surface * float(alpha_deg))
    return {
        'increment': float(increment),
        'ratio': ratio,
        'ivbw': vortex['ivbw'],
        'go2pav': vortex['go2pav'],
        'xd': vortex['xd'],
        'gated': False,
        'method': 'legacy_wbclb_vortex_term',
    }
