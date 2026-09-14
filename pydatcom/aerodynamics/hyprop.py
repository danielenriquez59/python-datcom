"""
HYPROP: equilibrium real-gas flow properties.

Table lookups per DATCOM Figures 6.3.1-31, -37, -43 and -49, giving free
stream pressure, temperature, Mach number and density behind the shock.
Each figure comes in six altitude variants A to F, covering sea level
through 300,000 feet in 50,000 foot steps.  The routine reads the two
variants bracketing the requested altitude and interpolates between them.

Two source defects were found while translating this; both are described on
:func:`calculate_hyprop` and neither is silently reproduced.

Reference: datcom-legacy/datcom_2000/hyprop.f
"""

import numpy as np
from typing import Dict, Optional, Sequence
import logging

from pydatcom.utils.legacy_tables import tlinex

logger = logging.getLogger(__name__)

# The six altitude variants.  The source's comments place block B at
# 100,000 feet and block F at 300,000, so A is 50,000 and the letters run in
# 50,000 foot steps.  Its band index IT then brackets an altitude between
# (IT-1)*50000 and IT*50000, which are the altitudes of blocks IT-1 and IT.
_LEVELS = 'ABCDEF'
_LEVEL_ALTITUDE = [50000.0, 100000.0, 150000.0, 200000.0, 250000.0, 300000.0]
_BAND = 50000.0


# Figure 6.3.1-31: pressure.
_F131A_X1 = [
    7, 5,
]
_F131A_X2 = [
    0.0239, 2.04, 4.15, 6.08, 8.05, 10.2, 12.1, 14.1, 16.2, 18.1, 20.1,
]
_F131A_Y = [
    1.06, 1.47, 2.03, 2.66, 3.47, 4.49, 5.63, 6.99, 8.47, 10.1, 11.8,
    1.06, 1.33, 1.69, 2.14, 2.59, 3.14, 3.75, 4.48, 5.28, 6.15, 7.04,
]
_F131B_X1 = [
    10, 7, 5,
]
_F131B_X2 = [
    0.0273, 2, 3.95, 6.04, 8, 9.92, 12.1, 14, 16.1, 18.1, 20.1,
]
_F131B_Y = [
    1.16, 1.53, 2.55, 3.58, 5.25, 7.12, 9.35, 11.9, 14.8, 17.8, 21.2,
    1.16, 1.38, 1.95, 2.68, 3.46, 4.49, 5.67, 7, 8.56, 10.1, 11.8,
    1.16, 1.29, 1.63, 1.99, 2.54, 3.15, 3.86, 4.57, 5.35, 6.11, 6.93,
]
_F131C_X1 = [
    15, 10,
]
_F131C_X2 = [
    -0.0309, 2.02, 3.99, 5.97, 7.96, 9.96, 12, 13.9, 16, 17.9, 19.9,
]
_F131C_Y = [
    1.4, 1.84, 3.56, 5.75, 8.98, 13.1, 17.8, 23, 28.7, 34.6, 41.4,
    1.39, 1.6, 2.17, 3.4, 5.14, 6.89, 9.05, 11.8, 14.4, 17.4, 20.7,
]
_F131D_X1 = [
    30, 20, 15,
]
_F131D_X2 = [
    0.0972, 2.11, 4.08, 6.02, 8.05, 10.1, 12.1, 14.2, 16.1, 18, 20,
]
_F131D_Y = [
    0.864, 3.28, 7.83, 16.7, 27.3, 41.5, 55.9, 75, 94.2, 114, 136,
    0.561, 2.36, 4.64, 8.78, 13.8, 21, 28.6, 37.9, 46.9, 57.3, 68.7,
    0.594, 1.44, 3.02, 5.85, 9.6, 13.6, 18, 23.4, 28.5, 34.9, 40.9,
]
_F131E_X1 = [
    30, 20,
]
_F131E_X2 = [
    0.00668, 2.02, 4.1, 6.04, 8.04, 10, 12.1, 14, 16, 18.1, 20.1,
]
_F131E_Y = [
    1.58, 3.97, 9.37, 17.1, 28.5, 43, 60.7, 79.8, 99.1, 120, 138,
    1.58, 2.35, 4.99, 8.94, 13.5, 20.9, 29.9, 39.2, 49.7, 62, 74.8,
]
_F131F_X1 = [
    30, 20,
]
_F131F_X2 = [
    0, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20,
]
_F131F_Y = [
    1, 4, 9.8, 18, 30, 44, 61.2, 79, 98.2, 118, 138,
    1, 2, 5, 9.1, 15, 22.3, 30.05, 40, 50.2, 61, 73,
]

# Figure 6.3.1-37: temperature.
_F137A_X1 = [
    7, 5,
]
_F137A_X2 = [
    -0.0116, 2, 4.07, 6.04, 7.97, 10, 12, 14, 16, 18, 20,
]
_F137A_Y = [
    1.01, 1.12, 1.24, 1.36, 1.52, 1.69, 1.89, 2.11, 2.35, 2.61, 2.89,
    1.01, 1.08, 1.16, 1.25, 1.34, 1.44, 1.56, 1.7, 1.84, 1.98, 2.13,
]
_F137B_X1 = [
    10, 7, 5,
]
_F137B_X2 = [
    -0.028, 2.03, 4.08, 6.08, 8, 10, 12.1, 14.1, 16.1, 18.1, 20,
]
_F137B_Y = [
    0.998, 1.13, 1.31, 1.54, 1.83, 2.12, 2.48, 2.88, 3.34, 3.83, 4.31,
    0.998, 1.09, 1.19, 1.32, 1.49, 1.68, 1.88, 2.1, 2.35, 2.58, 2.82,
    0.996, 1.07, 1.13, 1.23, 1.32, 1.44, 1.55, 1.69, 1.8, 1.94, 2.1,
]
_F137C_X1 = [
    15, 10,
]
_F137C_X2 = [
    -0.00997, 2.04, 4.03, 6.03, 8.07, 10.1, 12.1, 14.1, 16.1, 18.1, 20.1,
]
_F137C_Y = [
    1.13, 1.23, 1.58, 1.96, 2.5, 3.06, 3.75, 4.48, 5.32, 6.32, 7.29,
    1.11, 1.21, 1.44, 1.62, 1.87, 2.15, 2.46, 2.78, 3.23, 3.76, 4.24,
]
_F137D_X1 = [
    30, 20, 15,
]
_F137D_X2 = [
    0, 2.09, 4.09, 6.14, 8.09, 10, 12, 14, 16.1, 18.1, 20,
]
_F137D_Y = [
    1, 1.57, 2.49, 3.8, 5.5, 7.32, 8.95, 10.42, 11.8, 12.98, 14.1,
    1, 1.3, 1.81, 2.5, 3.34, 4.31, 5.41, 6.67, 7.89, 9.06, 10.16,
    1, 1.2, 1.52, 2, 2.47, 3.1, 3.82, 4.61, 5.5, 6.42, 7.35,
]
_F137E_X1 = [
    30, 20,
]
_F137E_X2 = [
    0.128, 2.04, 4.13, 6.08, 8.11, 10.1, 12.1, 14.1, 16.1, 18.1, 20,
]
_F137E_Y = [
    1.04, 1.63, 2.52, 3.79, 5.71, 7.91, 10.2, 12.2, 13.8, 15, 15.8,
    1.02, 1.37, 1.84, 2.54, 3.41, 4.48, 5.72, 7.09, 8.51, 10, 11.6,
]
_F137F_X1 = [
    30, 20,
]
_F137F_X2 = [
    0.0366, 1.99, 3.99, 6.06, 8.01, 10.1, 12.1, 14.2, 16.2, 18.1, 20.2,
]
_F137F_Y = [
    1.05, 1.68, 2.5, 3.88, 5.68, 7.96, 10.2, 11.7, 13.1, 14, 14.6,
    1.09, 1.36, 1.91, 2.58, 3.46, 4.59, 5.78, 7.16, 8.6, 9.93, 11.4,
]

# Figure 6.3.1-43: mach.
_F143A_X1 = [
    7, 5,
]
_F143A_X2 = [
    0.0463, 9.91, 20,
]
_F143A_Y = [
    7.04, 5.29, 3.56, 5.08, 4.04, 2.98,
]
_F143B_X1 = [
    10, 7, 5,
]
_F143B_X2 = [
    5.7e-07, 2, 4.04, 6.04, 8.03, 10.1, 12.1, 14, 16.1, 18.2, 20.1,
]
_F143B_Y = [
    10, 9.28, 8.58, 7.97, 7.31, 6.64, 5.98, 5.49, 5.06, 4.57, 4.17,
    7.02, 6.65, 6.27, 5.9, 5.53, 5.19, 4.83, 4.49, 4.21, 3.94, 3.69,
    5.05, 4.83, 4.61, 4.42, 4.15, 3.96, 3.77, 3.62, 3.38, 3.19, 2.95,
]
_F143C_X1 = [
    15, 10,
]
_F143C_X2 = [
    0.0501, 2.03, 4.01, 6.08, 8.02, 10.1, 12, 14.1, 16, 18, 20,
]
_F143C_Y = [
    15, 13.5, 12.1, 10.7, 9.52, 8.5, 7.67, 6.93, 6.28, 5.8, 5.38,
    10, 9.34, 8.63, 7.96, 7.38, 6.7, 6.23, 5.77, 5.32, 4.93, 4.64,
]
_F143D_X1 = [
    30, 20, 15,
]
_F143D_X2 = [
    0.0397, 2.11, 4, 6.03, 8.02, 9.91, 12, 14, 15.9, 17.9, 20,
]
_F143D_Y = [
    30, 23.7, 18.6, 15.3, 13.1, 11.5, 10.5, 9.62, 8.96, 8.3, 7.96,
    20, 17, 14.7, 12.7, 11, 9.73, 8.66, 7.67, 7.01, 6.4, 6.02,
    14.9, 13.4, 11.9, 10.5, 9.37, 8.43, 7.54, 6.79, 6.18, 5.66, 5.23,
]
_F143E_X1 = [
    30, 20,
]
_F143E_X2 = [
    0.000137, 1.93, 3.95, 6, 8.01, 10, 12, 14, 15.9, 18, 20,
]
_F143E_Y = [
    32, 24, 18.9, 15.1, 12.6, 10.8, 9.74, 8.83, 8.37, 7.82, 7.6,
    20.2, 17.4, 14.9, 12.5, 10.7, 9.28, 8.41, 7.45, 6.81, 6.21, 5.76,
]
_F143F_X1 = [
    30, 20,
]
_F143F_X2 = [
    -0.0155, 2.02, 4.12, 5.99, 7.96, 9.99, 12, 13.9, 15.9, 17.9, 19.8,
]
_F143F_Y = [
    29.9, 23.7, 18.7, 15.1, 12.4, 10.8, 9.71, 9.05, 8.42, 8.08, 7.84,
    19.9, 17.1, 14.5, 12.4, 10.7, 9.24, 8.03, 7.27, 6.56, 6.13, 5.79,
]

# Figure 6.3.1-49: density.
_F149A_X1 = [
    7, 5,
]
_F149A_X2 = [
    0.0479, 2.04, 4.13, 6.12, 8.06, 10.2, 12.1, 14.1, 16.1, 18.1, 20.1,
]
_F149A_Y = [
    1.02, 1.21, 1.4, 1.56, 1.69, 1.8, 1.89, 1.91, 1.92, 1.89, 1.83,
    1.01, 1.15, 1.29, 1.41, 1.52, 1.61, 1.69, 1.74, 1.76, 1.77, 1.75,
]
_F149B_X1 = [
    10, 7, 5,
]
_F149B_X2 = [
    0.0221, 1.93, 3.93, 6.06, 8, 9.93, 12.1, 14, 16, 18.1, 20.1,
]
_F149B_Y = [
    0.998, 1.29, 1.54, 1.75, 1.88, 1.96, 1.97, 1.92, 1.84, 1.77, 1.71,
    0.998, 1.21, 1.39, 1.56, 1.7, 1.78, 1.86, 1.91, 1.93, 1.91, 1.82,
    0.993, 1.14, 1.27, 1.41, 1.51, 1.6, 1.66, 1.72, 1.76, 1.77, 1.76,
]
_F149C_X1 = [
    15, 10,
]
_F149C_X2 = [
    0.0835, 2.07, 4.08, 6.11, 8.12, 10.1, 12.2, 14.1, 16.1, 18.2, 20.1,
]
_F149C_Y = [
    1.02, 1.46, 1.77, 1.97, 2.03, 2.04, 1.96, 1.83, 1.7, 1.59, 1.51,
    1.02, 1.3, 1.55, 1.77, 1.92, 2, 2.02, 1.96, 1.89, 1.82, 1.76,
]
_F149D_X1 = [
    30, 20, 15,
]
_F149D_X2 = [
    0.0763, 2.01, 4.01, 6.04, 7.99, 10.1, 12.1, 14, 16, 18.1, 19.9,
]
_F149D_Y = [
    1.01, 1.8, 2.03, 1.95, 1.68, 1.58, 1.55, 1.56, 1.57, 1.6, 1.64,
    1.01, 1.58, 1.91, 2.03, 1.99, 1.87, 1.69, 1.56, 1.49, 1.46, 1.46,
    1.01, 1.44, 1.78, 1.96, 2, 1.97, 1.91, 1.82, 1.71, 1.6, 1.5,
]
_F149E_X1 = [
    20, 30,
]
_F149E_X2 = [
    0.0232, 2.01, 3.95, 6, 8.05, 10.1, 12, 14.1, 16.1, 18, 20.1,
]
_F149E_Y = [
    1, 1.77, 2, 1.86, 1.64, 1.42, 1.27, 1.2, 1.25, 1.37, 1.41,
    1, 1.53, 1.93, 2, 1.93, 1.78, 1.58, 1.44, 1.33, 1.25, 1.18,
]
_F149F_X1 = [
    20, 30,
]
_F149F_X2 = [
    0.0213, 2, 4.05, 5.97, 8.02, 9.98, 12, 14, 16, 18, 20,
]
_F149F_Y = [
    1, 1.76, 1.99, 1.86, 1.65, 1.44, 1.29, 1.27, 1.38, 1.5, 1.56,
    1, 1.54, 1.92, 1.99, 1.89, 1.73, 1.58, 1.43, 1.33, 1.25, 1.23,
]


_FIGURES = {

    'pressure': {
        'A': (_F131A_X1, _F131A_X2, _F131A_Y),
        'B': (_F131B_X1, _F131B_X2, _F131B_Y),
        'C': (_F131C_X1, _F131C_X2, _F131C_Y),
        'D': (_F131D_X1, _F131D_X2, _F131D_Y),
        'E': (_F131E_X1, _F131E_X2, _F131E_Y),
        'F': (_F131F_X1, _F131F_X2, _F131F_Y),
    },
    'temperature': {
        'A': (_F137A_X1, _F137A_X2, _F137A_Y),
        'B': (_F137B_X1, _F137B_X2, _F137B_Y),
        'C': (_F137C_X1, _F137C_X2, _F137C_Y),
        'D': (_F137D_X1, _F137D_X2, _F137D_Y),
        'E': (_F137E_X1, _F137E_X2, _F137E_Y),
        'F': (_F137F_X1, _F137F_X2, _F137F_Y),
    },
    'mach': {
        'A': (_F143A_X1, _F143A_X2, _F143A_Y),
        'B': (_F143B_X1, _F143B_X2, _F143B_Y),
        'C': (_F143C_X1, _F143C_X2, _F143C_Y),
        'D': (_F143D_X1, _F143D_X2, _F143D_Y),
        'E': (_F143E_X1, _F143E_X2, _F143E_Y),
        'F': (_F143F_X1, _F143F_X2, _F143F_Y),
    },
    'density': {
        'A': (_F149A_X1, _F149A_X2, _F149A_Y),
        'B': (_F149B_X1, _F149B_X2, _F149B_Y),
        'C': (_F149C_X1, _F149C_X2, _F149C_Y),
        'D': (_F149D_X1, _F149D_X2, _F149D_Y),
        'E': (_F149E_X1, _F149E_X2, _F149E_Y),
        'F': (_F149F_X1, _F149F_X2, _F149F_Y),
    },
}


# Figure 6.3.1-43A is declared with a 3-entry second grid and 6 dependent
# values, but its call passes NX2=11 and so needs 22.  See below.
_SHORT_FIGURE = ('mach', 'A')
_SHORT_DECLARED = 6
_SHORT_REQUIRED = 22


def _lookup(prop: str, level: str, mach: float, alpha: float,
            replacement: Optional[Sequence[float]] = None) -> float:
    """One figure read, with the source's TLINEX end modes of zero."""
    x1, x2, y = _FIGURES[prop][level]
    if (prop, level) == _SHORT_FIGURE:
        if replacement is None:
            raise ValueError(
                f"Figure 6.3.1-43A is declared with {_SHORT_DECLARED} "
                f"dependent values but its source call requires "
                f"{_SHORT_REQUIRED}; the original reads past the end of the "
                "array. Supply a completed table to use this path.")
        y = list(replacement)
        x2 = list(_FIGURES['pressure']['A'][1])
    grid = np.asarray(y, dtype=float).reshape(
        (len(x2), len(x1)), order='F')
    return float(tlinex(x1, x2, grid, mach, alpha, 0, 0, 0, 0))


def calculate_hyprop(altitude: float, mach: float, alpha_deg: float,
                     mach_figure_43a: Optional[Sequence[float]] = None
                     ) -> Dict[str, object]:
    """Translate HYPROP: equilibrium real-gas properties behind the shock.

    Args:
        altitude: Geometric altitude, feet.
        mach: Free-stream Mach number.
        alpha_deg: Angle of attack, degrees.
        mach_figure_43a: A completed 22-value Figure 6.3.1-43A, needed only
            when the altitude band reaches the sea-level variant.

    Returns:
        Dictionary with ``pressure``, ``temperature``, ``mach`` and
        ``density``, plus the two bracketing levels and the interpolation
        fraction.

    Raises:
        ValueError: If the sea-level Mach figure is needed and no completed
            table is supplied.

    Notes:
        Two defects in the source, neither reproduced silently.

        Figure 6.3.1-43A is declared ``X2143A(3), Y3143A(6)`` but its call
        passes ``NX2=11``, so the original reads 16 elements past the end of
        a 6-element array.  Every sibling figure in that position holds 22
        values.  Rather than fabricate the missing sixteen, this path
        requires a completed table from the caller, as the WINGCL drag
        tables do.

        In the highest-altitude branch the source assigns the temperature
        figure's result to ``P(K)`` rather than ``T(K)``, so pressure is
        overwritten with a temperature reading and the temperature is never
        computed there.  That is plainly a typo; reproducing it would
        propagate a wrong pressure.  The translation assigns each figure to
        its own property and records the divergence here.
    """
    # The source's IT: band 1 up to 50,000 feet, then one band per 50,000.
    height = max(float(altitude), 0.0)
    band = 1 if height <= _BAND else min(int(np.ceil(height / _BAND)), 7)
    upper = min(band - 1, len(_LEVELS) - 1)
    lower = max(band - 2, 0)

    # The bracket the source interpolates over is (IT-1)*50000 to IT*50000.
    low_altitude = (band - 1) * _BAND
    high_altitude = band * _BAND
    fraction = (height - low_altitude) / (high_altitude - low_altitude)

    values = {}
    for prop in ('pressure', 'temperature', 'mach', 'density'):
        low = _lookup(prop, _LEVELS[lower], mach, alpha_deg, mach_figure_43a)
        high = _lookup(prop, _LEVELS[upper], mach, alpha_deg, mach_figure_43a)
        values[prop] = low + (high - low) * fraction

    return {
        **values,
        'band': band,
        'level_lower': _LEVELS[lower],
        'level_upper': _LEVELS[upper],
        'fraction': float(fraction),
        'method': 'legacy_hyprop',
    }
