"""
SUBWBT: subsonic wing-body-tail lateral dynamic derivatives.

Adds the vertical tail and ventral fin increments to ``CYP``, ``CNP`` and
``CNR``.  Both panels use the same geometry: an effective vertical arm that
rotates with angle of attack, and a moment arm about the CG.

The horizontal tail contributes nothing to these three derivatives, so the
source zeroes its entries and makes the body-wing-tail result identical to
body-wing.

Two source defects are preserved here rather than corrected; see
:func:`calculate_subwbt`.

Reference: datcom-legacy/datcom_2000/subwbt.f
"""

import numpy as np
from typing import Dict, Optional, Sequence
import logging

from pydatcom.utils.constants import RAD

logger = logging.getLogger(__name__)


def _panel(alpha_deg, arm_x: float, arm_z: float, cyb: float,
           blref: float, ventral: bool, cnr_cyb: float):
    """One panel's CYP, CNP and CNR increments."""
    alpha = np.atleast_1d(np.asarray(alpha_deg, dtype=float))
    cos_alpha = np.cos(alpha / RAD)
    sin_alpha = np.sin(alpha / RAD)

    effective = arm_z * cos_alpha - arm_x * sin_alpha          # ZEE
    moment_arm = arm_x * cos_alpha + arm_z * sin_alpha         # ABCDE

    if ventral:
        # The source writes (2.*ZEE-ZPF) here where the vertical tail branch
        # writes 2.*(ZEE-ZP).  See the note in calculate_subwbt.
        cyp = (2.0 * effective - arm_z) * cyb / blref
    else:
        cyp = 2.0 * (effective - arm_z) * cyb / blref

    cnp = -(2.0 * moment_arm * (effective - arm_z)) * cyb / blref**2
    # cnr_cyb is the sideslip derivative the source actually multiplies by,
    # which for the ventral fin is the vertical tail's.
    cnr = 2.0 * moment_arm**2 * cnr_cyb / blref**2
    return cyp, cnp, cnr


def calculate_subwbt(alpha_deg: Sequence[float],
                     cyp_wing_body: Sequence[float],
                     cnp_wing_body: Sequence[float],
                     cnr_wing_body: Sequence[float],
                     blref: float,
                     vertical_tail: Optional[Dict[str, float]] = None,
                     ventral_fin: Optional[Dict[str, float]] = None
                     ) -> Dict[str, object]:
    """Translate SUBWBT: subsonic lateral dynamic derivative buildup.

    Args:
        alpha_deg: Angle-of-attack schedule, degrees.
        cyp_wing_body: Wing-body ``CYP`` at each angle.
        cnp_wing_body: Wing-body ``CNP``.
        cnr_wing_body: Wing-body ``CNR``.
        blref: Lateral reference length.
        vertical_tail: Panel data as ``{'arm_x', 'arm_z', 'cyb'}``, the
            source's ``LP``, ``ZP`` and ``DYBV``.  Omit when absent.
        ventral_fin: The same for the ventral fin: ``LPF``, ``ZPF``,
            ``DYBF``.

    Returns:
        Dictionary with the total ``cyp``, ``cnp`` and ``cnr``, the per-panel
        increments, and the body-wing-tail and body-wing-vertical-tail
        results the source also fills.

    Raises:
        ValueError: If the arrays mismatch or the reference length is
            nonpositive.

    Notes:
        Two source defects are reproduced rather than corrected, so that the
        translation matches the original's numbers:

        The ventral fin ``CYP`` is written ``(2.*ZEE-ZPF)`` where the
        vertical tail writes ``2.*(ZEE-ZP)``.  A missing pair of parentheses
        is the obvious reading, but changing it would change results, so the
        source form stands and ``ventral_cyp_uses_source_grouping`` records
        it.

        The ventral fin ``CNR`` is multiplied by ``DYBV``, the *vertical
        tail's* sideslip derivative, where every other ventral term uses
        ``DYBF``.  The same substitution appears in ``CLRDER``, which writes
        ``VF(J+280)=2.*DCYBV*...`` with ventral geometry, so it is a
        recurring transcription slip in the original rather than a one-off.
        When no vertical tail is supplied the translation has no ``DYBV`` to
        use and falls back to ``DYBF``, which is flagged in the result.
    """
    alpha = np.atleast_1d(np.asarray(alpha_deg, dtype=float))
    arrays = [np.atleast_1d(np.asarray(value, dtype=float))
              for value in (cyp_wing_body, cnp_wing_body, cnr_wing_body)]
    if any(item.shape != alpha.shape for item in arrays):
        raise ValueError("SUBWBT needs wing-body arrays matching the schedule")
    if blref <= 0.0:
        raise ValueError("SUBWBT requires a positive lateral reference length")
    base_cyp, base_cnp, base_cnr = arrays

    zero = np.zeros_like(alpha)
    vt = {'cyp': zero, 'cnp': zero, 'cnr': zero}
    vf = {'cyp': zero, 'cnp': zero, 'cnr': zero}
    cnr_substitution = False

    if vertical_tail is not None:
        cyp, cnp, cnr = _panel(alpha, vertical_tail['arm_x'],
                               vertical_tail['arm_z'], vertical_tail['cyb'],
                               blref, ventral=False,
                               cnr_cyb=vertical_tail['cyb'])
        vt = {'cyp': cyp, 'cnp': cnp, 'cnr': cnr}

    if ventral_fin is not None:
        # The source's DYBV substitution in the ventral CNR term.
        if vertical_tail is not None:
            cnr_cyb = vertical_tail['cyb']
        else:
            cnr_cyb = ventral_fin['cyb']
            cnr_substitution = True
        cyp, cnp, cnr = _panel(alpha, ventral_fin['arm_x'],
                               ventral_fin['arm_z'], ventral_fin['cyb'],
                               blref, ventral=True, cnr_cyb=cnr_cyb)
        vf = {'cyp': cyp, 'cnp': cnp, 'cnr': cnr}

    total_cyp = base_cyp + vt['cyp'] + vf['cyp']
    total_cnp = base_cnp + vt['cnp'] + vf['cnp']
    total_cnr = base_cnr + vt['cnr'] + vf['cnr']

    return {
        'cyp': total_cyp,
        'cnp': total_cnp,
        'cnr': total_cnr,
        'vertical_tail': vt,
        'ventral_fin': vf,
        # The horizontal tail contributes nothing to these three.
        'horizontal_tail': {'cyp': zero, 'cnp': zero, 'cnr': zero},
        # Body-wing-tail equals body-wing; body-wing-vertical equals the total.
        'wing_body_tail': {'cyp': base_cyp, 'cnp': base_cnp,
                           'cnr': base_cnr},
        'wing_body_vertical': {'cyp': total_cyp, 'cnp': total_cnp,
                               'cnr': total_cnr},
        'ventral_cyp_uses_source_grouping': ventral_fin is not None,
        'ventral_cnr_fell_back_to_own_cyb': cnr_substitution,
        'method': 'legacy_subwbt',
    }
