"""
SUBPAW / SUBPAH: subsonic lifting-surface dynamic derivatives.

Gives ``CLQ``, ``CMQ``, ``CLAD`` and ``CMAD`` for a wing or a horizontal
tail.  The two source routines are identical apart from whitespace and the
COMMON offsets that select which surface they read, so one translation
serves both.

Pitching derivatives are closed-form, with a compressibility correction
applied above Mach 0.2.  The acceleration derivatives add two small
one-variable figure lookups and are gated: the source computes them only for
an untapered planform below Mach 1 with ``beta * aspect_ratio`` at or below
4.

Reference: datcom-legacy/datcom_2000/subpaw.f, subpah.f
"""

import numpy as np
from typing import Dict, Optional
import logging

from pydatcom.utils.constants import PI
from pydatcom.utils.legacy_interp import interx

logger = logging.getLogger(__name__)

# Figure 7.1.4.1-6, for the lift acceleration derivative.
_FIG_71410_6_X = [0., .25, .5, 1.5, 2.25, 2.75, 3.25, 4.]
_FIG_71410_6_Y = [0., .005, .014, .065, .095, .112, .125, .14]

# Figure 7.1.4.2-8, for the moment acceleration derivative.
_FIG_71420_8_X = [0., .25, .5, 1.25, 1.75, 2.25, 2.75, 3.25, 4.]
_FIG_71420_8_Y = [0., .003, .008, .033, .047, .057, .066, .072, .077]

# Above this Mach the source applies its compressibility correction to CMQ.
_COMPRESSIBILITY_MACH = 0.2
# The acceleration derivatives are limited to this beta*AR.
_BETA_ASPECT_LIMIT = 4.0


def calculate_subpaw(mach: float, cla: float, section_cla_compressible: float,
                     cos_sweep_c4: float, tan_sweep_c4: float,
                     aspect_ratio: float, area: float, sref: float,
                     mac: float, cbarr: float,
                     dcm_dcl: float, dcm_dclq: float,
                     xac_root_fraction: float = 0.0,
                     dxcg: float = 0.0,
                     taper_ratio: float = 0.0) -> Dict[str, object]:
    """Translate SUBPAW/SUBPAH: subsonic surface dynamic derivatives.

    Args:
        mach: Free-stream Mach number.
        cla: Surface lift-curve slope, the source's ``CCLAL``.
        section_cla_compressible: ``CLASM``, the section slope already
            multiplied by ``sqrt(|1 - M^2|)``.
        cos_sweep_c4: Cosine of the quarter-chord sweep.
        tan_sweep_c4: Tangent of the quarter-chord sweep.
        aspect_ratio: Exposed aspect ratio, ``ASTRW``.
        area: Exposed surface area, ``SW``.
        sref: Aircraft reference area.
        mac: Exposed mean aerodynamic chord, ``MACOE``.
        cbarr: Reference chord.
        dcm_dcl: ``DCMDCL``.
        dcm_dclq: ``DCMCLQ``.
        xac_root_fraction: ``XACCR``, used by the acceleration derivatives.
        dxcg: CG offset used to shift ``CMAD``.
        taper_ratio: ``LAMDA``; anything nonzero blocks the acceleration
            derivatives.

    Returns:
        Dictionary with ``clq`` and ``cmq`` always, and ``clad`` and ``cmad``
        when the source's gates allow them.  ``acceleration_available`` says
        which happened and ``acceleration_gate`` why.

    Raises:
        ValueError: If the references are nonpositive.
    """
    if min(sref, cbarr, mac, area) <= 0.0:
        raise ValueError("SUBPAW requires positive areas and chords")
    # The CMQ correction takes sqrt(1 - M^2 cos^2), which the source
    # evaluates before its own supersonic test, so a fast swept case
    # would silently produce a NaN there.
    if mach > _COMPRESSIBILITY_MACH and abs(mach * cos_sweep_c4) >= 1.0:
        raise ValueError(
            "SUBPAW's CMQ correction needs M*cos(sweep) < 1; got "
            f"{mach * cos_sweep_c4:.4f}. The source reaches its sqrt "
            "before the supersonic branch test and would return NaN.")

    # Pitching derivatives.
    clq = (0.5 + 2.0 * dcm_dclq) * cla * mac / cbarr

    cmq_incompressible = (
        -0.7 * section_cla_compressible * cos_sweep_c4 *
        (aspect_ratio * (0.5 * dcm_dcl + 2.0 * dcm_dcl**2) /
         (aspect_ratio + 2.0 * cos_sweep_c4) +
         aspect_ratio**3 * tan_sweep_c4**2 /
         ((aspect_ratio + 6.0 * cos_sweep_c4) * 24.0) + 0.125) *
        area / sref * (mac / cbarr)**2)

    if mach > _COMPRESSIBILITY_MACH:
        beta_sweep = np.sqrt(1.0 - mach**2 * cos_sweep_c4**2)
        swept_cubic = aspect_ratio**3 * tan_sweep_c4**2
        cmq = ((swept_cubic / (aspect_ratio * beta_sweep + 6.0 * cos_sweep_c4)
                + 3.0 / beta_sweep) /
               (swept_cubic / (aspect_ratio + 6.0 * cos_sweep_c4) + 3.0)) * \
            cmq_incompressible
    else:
        cmq = cmq_incompressible

    result = {
        'clq': float(clq),
        'cmq': float(cmq),
        'cmq_incompressible': float(cmq_incompressible),
        'method': 'legacy_subpaw',
    }

    # Acceleration derivative gates, in the source's order.
    if mach > 1.0:
        result['acceleration_available'] = False
        result['acceleration_gate'] = 'supersonic'
        return result
    if taper_ratio != 0.0:
        result['acceleration_available'] = False
        result['acceleration_gate'] = 'tapered_surface'
        return result

    beta = np.sqrt(1.0 - mach**2)
    beta_aspect = beta * aspect_ratio
    if beta_aspect > _BETA_ASPECT_LIMIT:
        result['acceleration_available'] = False
        result['acceleration_gate'] = 'beta_aspect_ratio_above_four'
        return result

    lift_factor = interx(1, _FIG_71410_6_X, [beta_aspect], [8],
                         _FIG_71410_6_Y)
    # The source guards both figure results against a zero beta.
    clg = (lift_factor * PI * aspect_ratio / (-2.0 * beta**2)
           if beta != 0.0 else 0.0)
    surface_to_ref = area * mac / (sref * cbarr)
    clad = (1.5 * xac_root_fraction * cla * sref / area +
            clg / 19.1) * surface_to_ref

    moment_factor = interx(1, _FIG_71420_8_X, [beta_aspect], [9],
                           _FIG_71420_8_Y)
    cmog = (moment_factor * aspect_ratio * 0.5 * PI / beta**2
            if beta != 0.0 else 0.0)
    cmad_plain = ((-2.53125 * xac_root_fraction**2 * cla * sref / area +
                   0.078532 * cmog) * surface_to_ref * mac / cbarr)
    cmad = cmad_plain + (dxcg / cbarr) * clad

    result.update({
        'clad': float(clad),
        'cmad': float(cmad),
        'cmad_before_cg_shift': float(cmad_plain),
        'clg': float(clg),
        'cmog': float(cmog),
        'beta_aspect_ratio': float(beta_aspect),
        'acceleration_available': True,
    })
    return result
