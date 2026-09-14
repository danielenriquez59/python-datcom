"""
Section and planform parameters: DELY and ARCLSS.

``DELY`` produces ``DELTAY``, the leading-edge sharpness parameter that both
``MAXCL`` and ``CLRDER`` take as their primary independent variable.  It is
the difference in upper-surface ordinate between the 6 percent and 0.15
percent chord stations, in percent chord.

``ARCLSS`` classifies the aspect ratio through Figure 4.1.3.4-24B and forms
the two derived ratios that the low-aspect-ratio lift path uses.

Reference: datcom-legacy/datcom_2000/dely.f, arclss.f
"""

import numpy as np
from typing import Dict, Sequence
import logging

from pydatcom.utils.legacy_numeric import tbfunx

logger = logging.getLogger(__name__)

# The two chord stations DELY samples, as fractions of chord.
_DELY_STATIONS = (0.0015, 0.0600)

# Figure 4.1.3.4-24B: aspect-ratio classification against taper ratio.
_FIG_41340_24B_TAPER = np.array([0., .1, .2, .3, .4, .5, .6, .7, .8, .9, 1.])
_FIG_41340_24B_VALUE = np.array([0., .225, .47, .496, .43, .32, .21, .125,
                                 .075, .0475, .0])


def calculate_dely(x: Sequence[float],
                   thickness: Sequence[float]) -> Dict[str, float]:
    """Translate DELY: the leading-edge sharpness parameter DELTAY.

    ``DELTAY = (y at 6% chord - y at 0.15% chord) * 100``

    Args:
        x: Chord stations as fractions of chord.
        thickness: Upper-surface ordinate at each station, the source's
            ``THN`` array.

    Returns:
        Dictionary with ``deltay`` and the two sampled ordinates.

    Raises:
        ValueError: If the arrays are mismatched or too short.

    Notes:
        The two stations are sampled with TBFUNX end mode 0, so a section
        whose coordinates start beyond 0.15 percent chord is clamped to its
        first point rather than extrapolated.
    """
    x = np.asarray(x, dtype=float)
    thickness = np.asarray(thickness, dtype=float)
    if x.ndim != 1 or x.shape != thickness.shape or x.size < 2:
        raise ValueError("DELY needs matching chord and ordinate arrays")

    near, far = (float(tbfunx(x, thickness, station, lower=0, upper=0)[0])
                 for station in _DELY_STATIONS)
    return {
        'deltay': (far - near) * 100.0,
        'y_at_0015': near,
        'y_at_06': far,
        'stations': _DELY_STATIONS,
        'method': 'legacy_dely',
    }


def calculate_arclss(taper_ratio: float, cos_sweep_le: float,
                     reference: float) -> Dict[str, float]:
    """Translate ARCLSS: aspect-ratio classification.

    Figure 4.1.3.4-24B gives a taper-dependent factor; the source then forms
    ``A(124) = (factor + 1) * cos(LE sweep)`` and ``A(125) = reference /
    A(124)``.

    Args:
        taper_ratio: Exposed taper ratio, the source's ``A(27)``.
        cos_sweep_le: Cosine of the leading-edge sweep, ``A(37)``.
        reference: The quantity the source divides by ``A(124)``, its
            ``A(128)``.

    Returns:
        Dictionary with ``factor`` (A(123)), ``classified`` (A(124)) and
        ``ratio`` (A(125)).

    Raises:
        ValueError: If the classified aspect ratio vanishes.
    """
    factor, _ = tbfunx(_FIG_41340_24B_TAPER, _FIG_41340_24B_VALUE,
                       float(taper_ratio), lower=0, upper=0)
    classified = (factor + 1.0) * cos_sweep_le
    if classified == 0.0:
        raise ValueError("ARCLSS ratio divides by a vanishing (A123+1)*cos")
    return {
        'factor': float(factor),
        'classified': float(classified),
        'ratio': float(reference / classified),
        'method': 'legacy_arclss',
    }
