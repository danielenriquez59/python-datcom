"""
Hypersonic aerodynamic calculations for PyDATCOM.

Implements hypersonic flow methods (Mach > 5):
- Newtonian impact theory
- Modified Newtonian method
- Hypersonic shock relations

Reference: datcom.f HYPBOD, HYPFLP, HYPROP subroutines
"""

import numpy as np
from typing import Dict, Optional
import logging

from pydatcom.utils.legacy_numeric import tbfunx

logger = logging.getLogger(__name__)


def _body_arrays(state: Dict) -> Optional[tuple[np.ndarray, np.ndarray]]:
    """Return validated station arrays required by legacy HYPBOD."""
    x = np.asarray(state.get('body_x', []), dtype=float)
    r = np.asarray(state.get('body_r', []), dtype=float)
    nx = int(state.get('body_nx', min(x.size, r.size)) or 0)
    if nx < 2 or x.size < nx or r.size < nx:
        return None
    x, r = x[:nx], r[:nx]
    if (not np.all(np.isfinite(x)) or not np.all(np.isfinite(r)) or
            np.any(np.diff(x) <= 0.0) or x[-1] <= 0.0 or np.any(r < 0.0)):
        return None
    return x, r


def _hypbod_coefficients(state: Dict, alpha_deg: float, mach: float,
                         x: np.ndarray, r: np.ndarray) -> Dict[str, float]:
    """Translate HYPBOD lines 246--310 for one angle of attack."""
    sref = float(state.get('options_sref', 1.0) or 1.0)
    cbar = float(state.get('options_cbarr', 1.0) or 1.0)
    xcg = float(state.get('synths_xcg', 0.0) or 0.0)
    if sref <= 0.0 or cbar <= 0.0:
        raise ValueError('options_sref and options_cbarr must be positive')

    # HYPBOD calls TBFUNX at every station and uses its local parabolic
    # derivative, even though the radius value interpolation is linear.
    drdx = np.array([tbfunx(x, r, station)[1] for station in x])
    theta = np.arctan(drdx)
    theta[-1] = theta[-2]  # HYPBOD: THETA(NX)=THETA(NX-1)

    angle = abs(np.deg2rad(alpha_deg))
    sa, ca = np.sin(angle), np.cos(angle)
    ta = np.tan(angle)
    ktheta = np.empty_like(r)
    kaf = np.empty_like(r)

    for n, local_theta in enumerate(theta):
        tn = np.tan(local_theta)
        cnn, sn = np.cos(local_theta), np.sin(local_theta)
        if angle > abs(local_theta):
            phe = np.arccos(np.clip(tn / ta, -1.0, 1.0))
        else:
            phe = 0.0 if local_theta > 0.0 else np.pi

        sp, cp = np.sin(phe), np.cos(phe)
        arg1 = (2.0 / 3.0) * (cnn * sa) ** 2 * sp * (cp ** 2 + 2.0)
        arg2 = 4.0 * sn * cnn * ca * sa * (
            np.pi / 2.0 - 0.5 * sp * cp - phe / 2.0)
        arg3 = 2.0 * (sn * ca) ** 2 * sp
        ktheta[n] = arg1 + arg2 + arg3

        arg6 = 2.0 * (ca * sn) ** 2 * tn * (np.pi - phe)
        arg7 = 4.0 * ca * sa * sp * sn ** 2
        arg8 = cnn * sn * sa ** 2 * (np.pi - phe - sp * cp)
        kaf[n] = arg6 + arg7 + arg8

    # HYPBOD integrates against X/RLB and multiplies by K*RLB/SR.  The
    # cancellation below keeps every moment arm dimensional until /CBAR.
    impact_factor = 1.833 * (1.0 - 0.4545 / mach ** 2)
    normal_integral = np.trapezoid(ktheta * r, x)
    axial_integral = np.trapezoid(kaf * r, x)
    moment_integral = np.trapezoid(ktheta * r * (xcg - x), x)

    sign = -1.0 if alpha_deg < 0.0 else 1.0
    cn = sign * impact_factor * normal_integral / sref
    cm = sign * impact_factor * moment_integral / (sref * cbar)
    ca_force = impact_factor * axial_integral / sref

    signed_alpha = np.deg2rad(alpha_deg)
    cl = cn * np.cos(signed_alpha) - ca_force * np.sin(signed_alpha)
    cd = ca_force * np.cos(signed_alpha) + cn * np.sin(signed_alpha)
    xcp = xcg - cm * cbar / cn if abs(cn) > 1e-14 else np.nan
    return {
        'cl': cl, 'cd': cd, 'cm': cm, 'cn': cn, 'ca': ca_force,
        'cp_max': impact_factor, 'xcp': xcp,
        'method': 'legacy_hypbod_nasa_tn_d176',
        'regime': 'hypersonic', 'mach': mach, 'alpha': alpha_deg,
    }


def calculate_hypersonic_coefficients(state: Dict, alpha_deg: float,
                                      mach: float) -> Dict[str, float]:
    """Calculate hypersonic coefficients at one angle of attack.

    Body stations select the legacy HYPBOD force integration.  The legacy
    wing/body assembly needs intermediate geometry absent from the current
    Python state, so a marked modified-Newtonian fallback remains for such
    configurations.
    """
    if mach <= 0.0:
        raise ValueError('mach must be positive')
    body = _body_arrays(state)
    if body is not None:
        return _hypbod_coefficients(state, alpha_deg, mach, *body)

    alpha = np.deg2rad(alpha_deg)
    sign = -1.0 if alpha_deg < 0.0 else 1.0
    cp_max = min(2.0, 1.84 + 0.032 * max(mach - 5.0, 0.0))
    cn = sign * cp_max * np.sin(alpha) ** 2
    ca_force = 0.2
    cl = cn * np.cos(alpha) - ca_force * np.sin(alpha)
    cd = ca_force * np.cos(alpha) + cn * np.sin(alpha)

    chord = float(state.get('wing_chrdtp', 0.0) or
                  state.get('wing_chrdr', 0.0) or
                  state.get('options_cbarr', 1.0) or 1.0)
    xle = float(state.get('synths_xw', 0.0) or 0.0)
    xcp = xle + 0.5 * chord
    xcg = float(state.get('synths_xcg', 0.0) or 0.0)
    cbar = float(state.get('options_cbarr', chord) or chord)
    cm = cn * (xcg - xcp) / cbar if cbar > 0.0 else 0.0
    return {
        'cl': cl, 'cd': cd, 'cm': cm, 'cn': cn, 'ca': ca_force,
        'cp_max': cp_max, 'xcp': xcp,
        'method': 'modified_newtonian_fallback',
        'regime': 'hypersonic', 'mach': mach, 'alpha': alpha_deg,
    }

