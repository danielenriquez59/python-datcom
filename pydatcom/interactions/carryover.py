"""
Wing/tail-body lift carryover factors, DATCOM Section 4.3.1.2.

Translates the interference ratios that the wing-body and wing-body-tail
buildups apply to the isolated-surface load:

- Figure 4.3.1.2-10  ``KWB``/``KBW``: carryover of the angle-of-attack load.
- Figure 4.3.1.2-12A ``KKWB``/``KKBW``: carryover of the incidence load.

``K_W(B)`` is the ratio of the lift on the surface in the presence of the
body to the isolated-surface lift; ``K_B(W)`` is the lift carried over onto
the body by the surface.  Both are tabulated against the body-radius to
semispan ratio ``r/s = (SSPN-SSPNE)/SSPN``.  The same figures serve the
wing-body and horizontal-tail-body cases, which is why ``CLWBT`` reads them
from ``SHB`` where ``WBTRAN`` reads them from ``SWB``.

Reference: datcom-legacy/datcom_2000/hbtran.f (Figure 4.3.1.2-10, labels
1030-1060), wbclb.f (Figure 4.3.1.2-12A), clwbt.f, wbtran.f
"""

import math

import numpy as np
from typing import Dict
import logging

from pydatcom.utils.constants import PI, RAD
from pydatcom.utils.legacy_numeric import tbfunx

logger = logging.getLogger(__name__)

# Figure 4.3.1.2-10: hbtran.f DATA TFIG10 / DKWB10 / DKBW10.
_FIG_431210_RATIO = np.array([0.0, .1, .2, .3, .4, .5, .6, .7, .8, .9, 1.0])
_FIG_431210_KWB = np.array([1.0, 1.08, 1.16, 1.26, 1.36, 1.46, 1.56, 1.67,
                            1.78, 1.89, 2.0])
_FIG_431210_KBW = np.array([0.0, .13, .29, .45, .62, .80, 1.0, 1.22, 1.45,
                            1.70, 2.0])

# Figure 4.3.1.2-12A: wbclb.f DATA X12A / Y12A1 / Y12A2.  The source writes
# the four repeated entries of Y12A1 as the FORTRAN repeat count "4*.94".
_FIG_431212A_RATIO = np.array([0., .1, .2, .3, .4, .5, .6, .7, .8, .9, 1.])
_FIG_431212A_KKWB = np.array([1., .97, .95, .94, .94, .94, .94, .95, .96,
                              .98, .99])
_FIG_431212A_KKBW = np.array([0., .11, .21, .31, .41, .51, .60, .70, .80,
                              .90, 1.0])


def body_semispan_ratio(state: Dict, component: str = 'wing') -> float:
    """Body-radius to semispan ratio ``r/s`` for a lifting surface.

    Both source call sites reduce to the same quantity: ``WBCLB`` forms
    ``FACT(1)=YB/WINGIN(4)`` with ``YB=SSPN-SSPNE``, and ``HBTRAN`` forms
    ``(SPAN-SPANS)/SPAN`` for KWB and ``DD/(2*SPAN)`` with ``DD=2*(SPAN-
    SPANS)`` for KBW.

    Args:
        state: State dictionary.
        component: Surface key prefix, ``'wing'`` or ``'htail'``.

    Returns:
        ``r/s`` in [0, 1].

    Raises:
        ValueError: If the surface spans are missing or inconsistent.
    """
    sspn = float(state.get(f'{component}_sspn', 0.0) or 0.0)
    sspne = float(state.get(f'{component}_sspne', sspn) or 0.0)
    if sspn <= 0.0:
        raise ValueError(f"{component} SSPN must be positive")
    if not 0.0 <= sspne <= sspn:
        raise ValueError(f"{component} SSPNE must lie in (0, SSPN]")
    return (sspn - sspne) / sspn


def fig4312_10(ratio: float) -> Dict[str, float]:
    """Figure 4.3.1.2-10: angle-of-attack lift carryover.

    Args:
        ratio: Body-radius to semispan ratio ``r/s``.

    Returns:
        Dictionary with ``kwb`` and ``kbw``.
    """
    kwb, _ = tbfunx(_FIG_431210_RATIO, _FIG_431210_KWB, float(ratio),
                    lower=0, upper=0)
    kbw, _ = tbfunx(_FIG_431210_RATIO, _FIG_431210_KBW, float(ratio),
                    lower=0, upper=0)
    return {'kwb': float(kwb), 'kbw': float(kbw)}


def fig4312_12a(ratio: float) -> Dict[str, float]:
    """Figure 4.3.1.2-12A: incidence lift carryover.

    Args:
        ratio: Body-radius to semispan ratio ``r/s``.

    Returns:
        Dictionary with ``kkwb`` and ``kkbw``.
    """
    kkwb, _ = tbfunx(_FIG_431212A_RATIO, _FIG_431212A_KKWB, float(ratio),
                     lower=0, upper=0)
    kkbw, _ = tbfunx(_FIG_431212A_RATIO, _FIG_431212A_KKBW, float(ratio),
                     lower=0, upper=0)
    return {'kkwb': float(kkwb), 'kkbw': float(kkbw)}


def calculate_carryover_factors(state: Dict,
                                component: str = 'wing') -> Dict[str, float]:
    """All four Section 4.3.1.2 carryover factors for a lifting surface.

    Args:
        state: State dictionary with the surface's SSPN and SSPNE.
        component: ``'wing'`` for KWB/KBW/KKWB/KKBW, ``'htail'`` for the
            KHB/KBH/KKHB/KKBH that ``CLWBT`` consumes.

    Returns:
        Dictionary with ``kwb``, ``kbw``, ``kkwb``, ``kkbw``, the ``ratio``
        they were read at, and ``method``.
    """
    ratio = body_semispan_ratio(state, component)
    factors = fig4312_10(ratio)
    factors.update(fig4312_12a(ratio))
    factors['ratio'] = ratio
    factors['method'] = 'legacy_fig_4312_10_and_12a'
    return factors


def intkbw(mach: float, sweep_le_deg: float, root_chord: float,
           body_diameter: float, afterbody_length: float):
    """Translate INTKBW: supersonic K_B(W) and its centre by integration.

    The carryover of wing lift onto the body behind a supersonic wing
    (Section 4.3.1.2), integrated over the body's lifting region with
    QUADIN on a 50 by 50 grid: the supersonic-leading-edge form (``ACOS``)
    when ``beta*cot(LE) > 1`` and the subsonic one (``SQRT``) otherwise.
    The region ends at the afterbody length ``DX``, or the Mach line if
    that comes first.

    Args:
        mach: ``MACH``; at or below 1 the source returns without setting
            its outputs.
        sweep_le_deg: ``OLE``, the leading-edge sweep in degrees.
        root_chord: ``CR``.
        body_diameter: ``D``.
        afterbody_length: ``DX``.

    Returns:
        ``(kbw, xac)``, ``xac`` in root chords; ``None`` for Mach 1 or
        below.

    Reference: datcom-legacy/datcom_2000/intkbw.f
    """
    from pydatcom.utils.legacy_numeric import quadin
    if mach <= 1.0:
        return None
    nn = 50
    xx = float(nn - 1)
    m = abs(1.0 / math.tan(sweep_le_deg / RAD))
    dn = body_diameter / xx
    beta = math.sqrt(mach**2 - 1.0)
    cr, dx = root_chord, afterbody_length
    if dx < body_diameter * beta - cr:
        dn = (cr + dx) / (xx * beta)
    supersonic_edge = beta * m > 1.0
    save, savm = [], []
    for row in range(nn):
        n = dn * row
        ll = beta * n
        ul = min(cr + ll, cr + dx)
        de = (ul - ll) / xx
        data, datm = [], []
        for col in range(nn):
            e = ll + de * col
            if supersonic_edge:
                val = 1. / (beta * m) if n == 0.0 else \
                    (e / beta + beta * m * n) / (n + m * e)
                val = min(val, 1.0)
                value = math.acos(val)
            else:
                val = 1. / (beta * m) if n == 0.0 else \
                    (e / beta - n) / (n + m * e)
                val = max(val, 0.0)
                value = math.sqrt(val)
            data.append(value)
            datm.append(e * value)
        save.append(quadin(data, de))
        savm.append(quadin(datm, de))
    kbw = quadin(save, dn)
    xac = quadin(savm, dn) / (kbw * cr)
    if supersonic_edge:
        kbw = (8. * beta * m / (PI * (body_diameter * cr / 2.) *
                                math.sqrt(beta**2 * m**2 - 1.)) * kbw)
    else:
        kbw = (16. * (beta * m)**1.5 / (PI * (body_diameter * cr / 2.) *
                                         (beta * m + 1.)) * kbw)
    return kbw, xac
