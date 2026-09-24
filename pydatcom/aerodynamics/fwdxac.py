"""
FWDXAC: aerodynamic centre of forward-swept wings, Figure 4.1.4.2-26.

The figure gives the aerodynamic centre in root chords aft of the root
chord leading edge as a function of three things: the curve parameter
``A*tan(sweep)``, the abscissa, and the exposed taper ratio.  The abscissa
is ``tan(sweep)/beta`` while that is at most one and ``beta/tan(sweep)``
beyond, so it runs 0 to 1 and back to 0 across each speed regime, reaching
zero at Mach one from both sides.  That gives four 216-value tables, one per
half of the abscissa in each regime:

====================  ======  ==================  ====================
table                  regime  abscissa            source arrays
====================  ======  ==================  ====================
``TYSUBL``             sub     tan/beta <= 1       SUBT1 to SUBT6
``TYSUBR``             sub     beta/tan  < 1       SUBT7 to SUBT12
``TYSUPL``             super   beta/tan  < 1       SUPT1 to SUPT6
``TYSUPR``             super   tan/beta <= 1       SUPT7 to SUPT12
====================  ======  ==================  ====================

Each table is six taper blocks of six curves by six abscissa points, looked
up with TLIN3X under end modes (2,0,0) and (2,0,0): quadratic extrapolation
along the curve parameter, clamping along the other two.

The table values were extracted from the source DATA statements by parsing
rather than by hand.

Reference: datcom-legacy/datcom_2000/fwdxac.f
"""

import numpy as np
from typing import Dict
import logging

from pydatcom.utils.legacy_tables import tlin3x

logger = logging.getLogger(__name__)

# TCURVE, TX and TGRAPH: the curve parameter A*tan(sweep), the abscissa and
# the taper ratio.
_CURVE = np.array([1., 2., 3., 4., 5., 6.])
_ABSCISSA = np.array([0., .2, .4, .6, .8, 1.])
_TAPER = np.array([0., .2, .25, .33, .5, 1.])


def _table(flat):
    """Reshape a source table into TLIN3X's ``(x2, x1, x3)`` layout.

    Within each taper block the abscissa varies fastest, each run of six
    being one curve, which is the ``Y(NX2,NX1,NX3)`` order the translated
    TLIN3X takes.
    """
    return np.asarray(flat, dtype=float).reshape(6, 6, 6).transpose(2, 1, 0)


_SUBSONIC_TAN_BETA = _table([
    # taper 0.0: SUBT1
    0.08, 0.04, 0.02, 0.01, 0.00, -0.01,
    -0.01, -0.03, -0.05, -0.06, -0.07, -0.08,
    -0.09, -0.11, -0.12, -0.13, -0.14, -0.15,
    -0.18, -0.19, -0.20, -0.21, -0.22, -0.22,
    -0.27, -0.28, -0.29, -0.29, -0.29, -0.30,
    -0.36, -0.37, -0.37, -0.38, -0.38, -0.38,
    # taper 0.2: SUBT2
    0.05, 0.03, 0.02, 0.01, -0.01, -0.01,
    -0.06, -0.07, -0.08, -0.09, -0.10, -0.10,
    -0.18, -0.19, -0.19, -0.20, -0.20, -0.20,
    -0.30, -0.31, -0.31, -0.31, -0.31, -0.31,
    -0.42, -0.42, -0.42, -0.42, -0.42, -0.41,
    -0.54, -0.54, -0.54, -0.53, -0.53, -0.52,
    # taper 0.25: SUBT3
    0.04, 0.03, 0.01, 0.00, -0.01, -0.02,
    -0.08, -0.09, -0.10, -0.11, -0.11, -0.11,
    -0.21, -0.21, -0.21, -0.21, -0.22, -0.22,
    -0.34, -0.33, -0.33, -0.33, -0.33, -0.33,
    -0.46, -0.46, -0.46, -0.45, -0.44, -0.44,
    -0.59, -0.58, -0.58, -0.57, -0.56, -0.55,
    # taper 0.33: SUBT4
    0.04, 0.03, 0.01, 0.00, -0.01, -0.02,
    -0.10, -0.10, -0.11, -0.11, -0.12, -0.12,
    -0.24, -0.24, -0.23, -0.23, -0.23, -0.23,
    -0.38, -0.37, -0.37, -0.36, -0.36, -0.36,
    -0.52, -0.51, -0.50, -0.49, -0.49, -0.48,
    -0.66, -0.65, -0.64, -0.63, -0.62, -0.61,
    # taper 0.5: SUBT5
    0.03, 0.02, 0.01, 0.00, -0.01, -0.02,
    -0.14, -0.13, -0.13, -0.13, -0.14, -0.14,
    -0.31, -0.30, -0.29, -0.28, -0.28, -0.27,
    -0.47, -0.45, -0.44, -0.43, -0.43, -0.42,
    -0.64, -0.62, -0.60, -0.59, -0.57, -0.56,
    -0.82, -0.78, -0.76, -0.75, -0.73, -0.72,
    # taper 1.0: SUBT6
    0.01, 0.02, 0.01, 0.00, -0.01, -0.02,
    -0.24, -0.21, -0.20, -0.20, -0.19, -0.19,
    -0.50, -0.45, -0.43, -0.41, -0.40, -0.39,
    -0.75, -0.69, -0.65, -0.63, -0.62, -0.60,
    -1.00, -0.94, -0.90, -0.86, -0.84, -0.82,
    -1.25, -1.19, -1.13, -1.09, -1.06, -1.04,
])

_SUBSONIC_BETA_TAN = _table([
    # taper 0.0: SUBT7
    -0.04, -0.03, -0.03, -0.02, -0.02, -0.01,
    -0.10, -0.09, -0.09, -0.09, -0.08, -0.08,
    -0.16, -0.16, -0.16, -0.15, -0.15, -0.15,
    -0.23, -0.23, -0.23, -0.23, -0.22, -0.22,
    -0.29, -0.30, -0.30, -0.30, -0.30, -0.30,
    -0.37, -0.37, -0.37, -0.37, -0.38, -0.38,
    # taper 0.2: SUBT8
    -0.05, -0.05, -0.04, -0.03, -0.02, -0.01,
    -0.12, -0.12, -0.12, -0.11, -0.11, -0.10,
    -0.21, -0.21, -0.21, -0.21, -0.21, -0.20,
    -0.30, -0.30, -0.30, -0.31, -0.31, -0.31,
    -0.40, -0.40, -0.40, -0.41, -0.41, -0.41,
    -0.49, -0.50, -0.51, -0.51, -0.52, -0.52,
    # taper 0.25: SUBT9
    -0.05, -0.04, -0.04, -0.03, -0.03, -0.02,
    -0.12, -0.12, -0.12, -0.12, -0.12, -0.11,
    -0.22, -0.22, -0.22, -0.22, -0.22, -0.22,
    -0.32, -0.32, -0.32, -0.32, -0.33, -0.33,
    -0.42, -0.42, -0.43, -0.43, -0.43, -0.44,
    -0.52, -0.53, -0.53, -0.54, -0.55, -0.55,
    # taper 0.33: SUBT10
    -0.04, -0.04, -0.03, -0.03, -0.02, -0.02,
    -0.14, -0.14, -0.13, -0.13, -0.12, -0.12,
    -0.23, -0.23, -0.24, -0.23, -0.23, -0.23,
    -0.34, -0.35, -0.35, -0.35, -0.35, -0.36,
    -0.46, -0.46, -0.46, -0.47, -0.48, -0.48,
    -0.57, -0.57, -0.58, -0.59, -0.60, -0.61,
    # taper 0.5: SUBT11
    -0.05, -0.04, -0.04, -0.03, -0.02, -0.02,
    -0.15, -0.15, -0.15, -0.14, -0.14, -0.14,
    -0.26, -0.27, -0.27, -0.27, -0.27, -0.27,
    -0.39, -0.40, -0.40, -0.41, -0.41, -0.42,
    -0.52, -0.53, -0.54, -0.55, -0.55, -0.56,
    -0.66, -0.67, -0.68, -0.70, -0.71, -0.72,
    # taper 1.0: SUBT12
    -0.07, -0.06, -0.05, -0.04, -0.03, -0.02,
    -0.18, -0.19, -0.19, -0.19, -0.19, -0.19,
    -0.35, -0.36, -0.37, -0.37, -0.38, -0.39,
    -0.53, -0.54, -0.56, -0.57, -0.59, -0.60,
    -0.73, -0.74, -0.76, -0.78, -0.80, -0.82,
    -0.94, -0.95, -0.97, -0.99, -1.02, -1.04,
])

_SUPERSONIC_BETA_TAN = _table([
    # taper 0.0: SUPT1
    -0.04, -0.04, -0.04, -0.03, -0.03, -0.02,
    -0.10, -0.10, -0.09, -0.08, -0.07, -0.05,
    -0.16, -0.16, -0.15, -0.14, -0.13, -0.11,
    -0.23, -0.23, -0.22, -0.21, -0.20, -0.18,
    -0.29, -0.29, -0.28, -0.27, -0.26, -0.24,
    -0.37, -0.37, -0.36, -0.35, -0.33, -0.31,
    # taper 0.2: SUPT2
    -0.04, -0.04, -0.04, -0.04, -0.03, -0.01,
    -0.12, -0.12, -0.11, -0.10, -0.09, -0.07,
    -0.20, -0.20, -0.20, -0.19, -0.18, -0.15,
    -0.30, -0.30, -0.29, -0.28, -0.26, -0.23,
    -0.40, -0.40, -0.39, -0.37, -0.34, -0.30,
    -0.49, -0.49, -0.48, -0.46, -0.45, -0.42,
    # taper 0.25: SUPT3
    -0.04, -0.04, -0.04, -0.03, -0.02, -0.01,
    -0.12, -0.12, -0.11, -0.10, -0.08, -0.05,
    -0.22, -0.22, -0.21, -0.20, -0.18, -0.15,
    -0.31, -0.31, -0.30, -0.29, -0.27, -0.23,
    -0.42, -0.42, -0.40, -0.38, -0.34, -0.31,
    -0.51, -0.51, -0.50, -0.47, -0.43, -0.38,
    # taper 0.33: SUPT4
    -0.04, -0.04, -0.04, -0.03, -0.01, 0.02,
    -0.13, -0.13, -0.12, -0.10, -0.08, -0.02,
    -0.23, -0.23, -0.22, -0.20, -0.18, -0.14,
    -0.34, -0.34, -0.33, -0.31, -0.28, -0.24,
    -0.45, -0.45, -0.43, -0.41, -0.38, -0.35,
    -0.56, -0.56, -0.54, -0.52, -0.49, -0.45,
    # taper 0.5: SUPT5
    -0.05, -0.05, -0.05, -0.04, -0.02, 0.02,
    -0.14, -0.14, -0.12, -0.10, -0.07, -0.03,
    -0.26, -0.26, -0.24, -0.23, -0.21, -0.16,
    -0.39, -0.39, -0.37, -0.35, -0.32, -0.28,
    -0.51, -0.51, -0.49, -0.46, -0.43, -0.39,
    -0.65, -0.65, -0.63, -0.62, -0.60, -0.58,
    # taper 1.0: SUPT6
    -0.08, -0.08, -0.08, -0.04, 0.02, 0.08,
    -0.17, -0.17, -0.15, -0.12, -0.09, -0.06,
    -0.34, -0.34, -0.32, -0.29, -0.27, -0.24,
    -0.51, -0.51, -0.48, -0.46, -0.43, -0.41,
    -0.71, -0.71, -0.68, -0.56, -0.64, -0.63,
    -0.92, -0.92, -0.90, -0.87, -0.85, -0.84,
])

_SUPERSONIC_TAN_BETA = _table([
    # taper 0.0: SUPT7
    0.22, 0.20, 0.10, 0.03, 0.00, -0.02,
    0.13, 0.13, 0.12, 0.04, -0.02, -0.05,
    0.04, 0.04, 0.05, -0.01, -0.07, -0.11,
    -0.04, -0.04, -0.03, -0.07, -0.14, -0.18,
    -0.13, -0.13, -0.12, -0.11, -0.20, -0.24,
    -0.22, -0.22, -0.21, -0.17, -0.26, -0.31,
    # taper 0.2: SUPT8
    0.20, 0.19, 0.13, 0.06, 0.02, -0.01,
    0.08, 0.09, 0.03, 0.03, -0.03, -0.07,
    -0.04, -0.04, -0.02, -0.04, -0.10, -0.15,
    -0.17, -0.16, -0.14, -0.12, -0.19, -0.23,
    -0.28, -0.28, -0.26, -0.21, -0.25, -0.30,
    -0.40, -0.39, -0.38, -0.36, -0.37, -0.42,
    # taper 0.25: SUPT9
    0.19, 0.20, 0.15, 0.08, 0.02, -0.01,
    0.06, 0.07, 0.08, 0.04, -0.01, -0.05,
    -0.07, -0.06, -0.03, -0.05, -0.10, -0.15,
    -0.19, -0.18, -0.16, -0.13, -0.19, -0.23,
    -0.32, -0.31, -0.29, -0.29, -0.27, -0.31,
    -0.44, -0.44, -0.43, -0.40, -0.33, -0.38,
    # taper 0.33: SUPT10
    0.18, 0.19, 0.16, 0.11, 0.06, 0.02,
    0.00, 0.05, 0.06, 0.07, 0.04, -0.02,
    -0.09, -0.08, -0.06, -0.05, -0.10, -0.14,
    -0.19, -0.18, -0.16, -0.16, -0.19, -0.24,
    -0.37, -0.37, -0.34, -0.29, -0.30, -0.35,
    -0.51, -0.51, -0.49, -0.45, -0.41, -0.45,
    # taper 0.5: SUPT11
    0.18, 0.19, 0.16, 0.11, 0.06, 0.02,
    0.02, 0.03, 0.05, 0.04, 0.00, -0.03,
    -0.15, -0.14, -0.11, -0.09, -0.12, -0.16,
    -0.32, -0.31, -0.28, -0.23, -0.24, -0.28,
    -0.48, -0.48, -0.45, -0.40, -0.35, -0.39,
    -0.65, -0.65, -0.62, -0.57, -0.55, -0.58,
    # taper 1.0: SUPT12
    0.19, 0.20, 0.19, 0.16, 0.12, 0.08,
    -0.05, -0.03, -0.01, -0.01, -0.03, -0.06,
    -0.30, -0.29, -0.25, -0.20, -0.22, -0.24,
    -0.55, -0.54, -0.50, -0.42, -0.39, -0.41,
    -0.80, -0.79, -0.75, -0.67, -0.67, -0.63,
    -1.05, -1.04, -1.00, -0.92, -0.84, -0.84,
])

_TABLES = {
    'TYSUBL': _SUBSONIC_TAN_BETA,
    'TYSUBR': _SUBSONIC_BETA_TAN,
    'TYSUPL': _SUPERSONIC_BETA_TAN,
    'TYSUPR': _SUPERSONIC_TAN_BETA,
}


def _lookup(table: str, curve: float, abscissa: float, taper: float) -> float:
    return float(tlin3x(
        _CURVE, _ABSCISSA, _TAPER, _TABLES[table],
        curve, abscissa, taper, 2, 0, 0, 2, 0, 0,
    ))


def calculate_fwdxac(atnswp: float, taper: float, factor: float,
                     mach: float) -> Dict[str, object]:
    """Translate FWDXAC: aerodynamic centre from Figure 4.1.4.2-26.

    Args:
        atnswp: ``ATNSWP``, aspect ratio times the tangent of the sweep.
        taper: ``TAPER``, the exposed taper ratio.
        factor: ``FACTOR``, ``tan(sweep)/beta``.
        mach: ``MACH``.  Only its side of one matters.

    Returns:
        Dictionary with ``xac``, the source's result in root chords from the
        root chord leading edge, the ``table`` it came from and the
        ``abscissa`` used.  A subsonic lookup with ``tan/beta > 1`` also
        reports ``source_defect`` and ``xac_tysubr``; see Notes.

    Notes:
        The subsonic branch for ``tan(sweep)/beta > 1`` looks up ``TYSUPR``,
        the supersonic tan/beta table, at ``beta/tan``.  ``TYSUBR`` is filled
        with the subsonic beta/tan data and read by no branch at all.  The
        data settles which is meant: ``TYSUBR`` equals ``TYSUBL`` exactly at
        the abscissa of one and ``TYSUPL`` to within 0.02 at zero, so it
        alone joins its neighbours on both sides, while ``TYSUPR`` in its
        place leaves a step of up to 0.20 root chords at ``tan/beta = 1``.
        This is a one-letter slip, but the executable reads ``TYSUPR``, so
        the source result is returned and the ``TYSUBR`` value is reported
        alongside for comparison.

        SUPT6 row five reads ``-.71, -.71, -.68, -.56, -.64, -.63``.  The
        ``-.56`` breaks a row that otherwise falls smoothly and looks like
        a digitisation slip for about ``-.66``.  Chart data is not
        invented here, so it stands.

        The source sets ``FACTI=1.E10`` for a zero ``FACTOR``, but that
        value only reaches a lookup when ``|FACTOR| > 1``, so it is never
        used.
    """
    abs_factor = abs(float(factor))
    abs_atnswp = abs(float(atnswp))
    subsonic = float(mach) < 1.0
    beta_over_tan = abs_factor > 1.0
    abscissa = 1.0 / abs_factor if beta_over_tan else abs_factor

    if subsonic and not beta_over_tan:
        table = 'TYSUBL'
    elif subsonic and beta_over_tan:
        table = 'TYSUPR'
    elif beta_over_tan:
        table = 'TYSUPL'
    else:
        table = 'TYSUPR'

    taper = float(taper)
    result = {
        'xac': _lookup(table, abs_atnswp, abscissa, taper),
        'table': table,
        'abscissa': abscissa,
        'method': 'legacy_fwdxac',
    }

    if subsonic and beta_over_tan:
        result['source_defect'] = 'subsonic_beta_over_tan_reads_TYSUPR'
        result['xac_tysubr'] = _lookup('TYSUBR', abs_atnswp, abscissa, taper)

    return result
