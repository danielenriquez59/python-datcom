"""
Supersonic downwash, Figure 4.4.1-76.

The figure is split across five source routines by wing taper ratio, each a
chain of table lookups in which every stage feeds the next:

- ``SDWA``  Figure 4.4.1-76A, taper 0.00
- ``SDWB``  Figure 4.4.1-76B, taper 1.00
- ``SDWC``  Figure 4.4.1-76C, taper 0.25 and 0.50
- ``SDWD``  Figure 4.4.1-76D, taper 0.00, 0.25 and 0.50
- ``SDWE``  Figure 4.4.1-76E, general taper

Each stage is an ``INTERX`` two-variable lookup whose dependent tables use
the flat ``TABLE(LIND,4)`` packing, so the grid for variable ``k`` starts at
offset ``k*LIND`` regardless of how many entries the previous grid used.
The source end modes are carried through per call.

Reference: datcom-legacy/datcom_2000/sdwa.f, sdwb.f, sdwc.f, sdwd.f, sdwe.f
"""

import numpy as np
from typing import Dict
import logging

from pydatcom.utils.legacy_interp import interx
from pydatcom.utils.legacy_numeric import tbfunx

logger = logging.getLogger(__name__)

# --------------------------------------------------------------------------
# Figure 4.4.1-76A, taper 0.00.  sdwa.f
# --------------------------------------------------------------------------
_XA1 = [.90, 1.0, 1.2, 1.4, 1.6, 2.4, 0.0, 0.1, 0.3, 0.5]
_YA1 = [.90, 1.18, 1.65, 1.95, 2.10, 2.34,
        .55, 0.84, 1.30, 1.60, 1.74, 2.00,
        0.0, 0.00, 0.50, 0.87, 1.10, 1.50,
        0.0, 0.00, 0.00, 0.30, 0.61, 1.00]
_XA2 = [0.0, 3.0, 0.0, 4.0, 8.0, 12.]
_YA2 = [0.0, 4.45, 0.0, 3.00, 0.0, 2.50]
_XA3 = [0.0, 3.0, 0.0, 0.0, 0.0, 0.0, 0.3, 0.4, 0.6, 0.8]
_YA3 = [0.00, 0.60, -.04, 0.54, -.09, 0.44, -.17, 0.35, -.28, 0.26]

# --------------------------------------------------------------------------
# Figure 4.4.1-76B, taper 1.00.  sdwb.f
# --------------------------------------------------------------------------
_XB1 = [1.0, 1.2, 1.4, 1.5, 1.6, 2.0, 2.2, 2.4, 0.0, 0.1, 0.3, 0.5]
_YB1 = [0.32, 0.74, 1.21, 1.39, 1.44, 1.56, 1.62, 1.77,
        0.11, 0.60, 1.10, 1.25, 1.30, 1.42, 1.50, 1.62,
        0.00, 0.20, 0.76, 1.00, 1.10, 1.34, 1.42, 1.47,
        0.00, 0.00, 0.00, 0.50, 0.75, 1.06, 1.12, 1.20]
_XB2 = [0.00, .8, 1.70, 3.0, 2.00, 4., 8.00, 12.]
_YB2 = [0.00, .60, 4.43, 12.7,
        0.00, .80, 1.70, 3.00,
        0.00, .40, 0.80, 1.40,
        0.00, .27, 0.50, 0.94]
_XB3 = [0.0, 5.0, 0.0, 0.0, .15, .30]
_YB3 = [0.0, 1.00, 0.0, 1.10, 0.0, 1.20]


def calculate_sdwa(x: float, y: float, z: float,
                   beta_aspect: float) -> Dict[str, float]:
    """Translate SDWA: Figure 4.4.1-76A downwash, taper 0.00.

    Three chained lookups: 76A1 on ``(x, z)``, 76A2 on that result and the
    compressible aspect-ratio parameter, then 76A3 on that result and ``y``.

    Args:
        x: Streamwise station.
        y: Spanwise station.
        z: Vertical station.
        beta_aspect: ``AB``, the compressible aspect-ratio parameter.

    Returns:
        Dictionary with ``sdw`` and the two intermediate stage values.
    """
    stage1 = interx(2, _XA1, [x, z], [6, 4], _YA1, lind=6,
                    lx1l=2, lx2l=0, lx1u=1, lx2u=2)
    stage2 = interx(2, _XA2, [stage1, beta_aspect], [2, 3], _YA2, lind=3,
                    lx1l=0, lx2l=1, lx1u=1, lx2u=1)
    sdw = interx(2, _XA3, [stage2, y], [2, 5], _YA3, lind=5,
                 lx1l=0, lx2l=0, lx1u=1, lx2u=2)
    return {
        'sdw': float(sdw),
        'stage1': float(stage1),
        'stage2': float(stage2),
        'taper': 0.00,
        'method': 'legacy_sdwa',
    }


def calculate_sdwb(x: float, y: float, z: float,
                   beta_aspect: float) -> Dict[str, float]:
    """Translate SDWB: Figure 4.4.1-76B downwash, taper 1.00.

    Args:
        x: Streamwise station.
        y: Spanwise station.
        z: Vertical station.
        beta_aspect: ``AB``, the compressible aspect-ratio parameter.

    Returns:
        Dictionary with ``sdw`` and the two intermediate stage values.
    """
    stage1 = interx(2, _XB1, [x, z], [8, 4], _YB1, lind=8,
                    lx1l=2, lx2l=0, lx1u=1, lx2u=2)
    stage2 = interx(2, _XB2, [stage1, beta_aspect], [4, 4], _YB2, lind=4,
                    lx1l=0, lx2l=2, lx1u=1, lx2u=2)
    sdw = interx(2, _XB3, [stage2, y], [2, 3], _YB3, lind=3,
                 lx1l=0, lx2l=0, lx1u=1, lx2u=2)
    return {
        'sdw': float(sdw),
        'stage1': float(stage1),
        'stage2': float(stage2),
        'taper': 1.00,
        'method': 'legacy_sdwb',
    }


# --------------------------------------------------------------------------
# Figure 4.4.1-76C, taper 0.25 and 0.50.  sdwc.f
# --------------------------------------------------------------------------
_XC1 = [.900, .950, 1.12, 1.40, 1.80, 2.20, 2.40, 0.00, 0.10, 0.30, 0.50]
_YC1 = [0.80, 1.00, 1.30, 1.63, 1.96, 2.08, 2.10,
        0.00, 0.70, 1.00, 1.39, 1.70, 1.81, 1.80,
        0.00, 0.00, 0.33, 0.85, 1.20, 1.40, 1.40,
        0.00, 0.00, 0.00, 0.38, 0.79, 0.99, 1.00]
# The source writes the padding between the two packed grids as "4*0.".
_XC2 = [0.00, 2.50, 0., 0., 0., 0.,
        2.70, 3.20, 5.40, 6.40, 10.7, 12.8]
_YC2 = [0.0, 4.60, 0.0, 3.55, 0.0, 2.90, 0.0, 2.50, 0.0, 1.84, 0.0, 1.59]
_XC3 = np.array([0.0, 3.00, 3.50, 3.75, 4.00])
_YC3A = np.array([0.0, 3.00, 3.50, 3.75, 4.00])
_YC3B = np.array([0.0, 2.35, 3.00, 3.50, 3.50])
_XC4 = [0.0, 4.00, 0.0, .30]
_YC4 = [0.0, .780, 0.0, .820]

# --------------------------------------------------------------------------
# Figure 4.4.1-76D, taper 0.00, 0.25 and 0.50.  sdwd.f
# --------------------------------------------------------------------------
_XD1 = [1.0, 1.4, 1.8, 2.2, 2.4, 0.0, 0.1, 0.3, 0.5]
_YD1 = [1.68, 2.77, 3.38, 3.60, 3.62,
        1.05, 2.11, 2.78, 3.00, 3.00,
        0.00, 1.20, 1.91, 2.23, 2.30,
        0.00, 0.20, 1.10, 1.45, 1.51]
_XD2 = [0.0, .2, .3, 1.0, 4.0, 2.0, 4., 8., 12.]
_YD2 = [0.00, .100, .200, 1.00, 7.00,
        0.00, .200, .300, 1.00, 4.00,
        0.00, .600, .800, 1.20, 2.40,
        0.00, .400, .600, 1.00, 1.60]
_XD3 = np.array([0.0, 5.5])
_YD3A = np.array([0.0, 6.3])
_YD3B = np.array([0.0, 5.5])
_YD3C = np.array([0.0, 5.0])
_XD4 = [0.0, 0.6, 1.0, 1.5, 6.0, 0.0, .15, .30]
_YD4 = [0.00, .120, .200, .300, 1.20,
        .040, .040, .180, .280, 1.14,
        .040, .040, .120, .260, .960]

# --------------------------------------------------------------------------
# Figure 4.4.1-76E, general taper.  sdwe.f
#
# XE1/YE1 are identical to XD1/YD1, and XE4/YE4 to XD4/YD4; the source
# repeats both DATA blocks rather than sharing them.  They are kept separate
# here so each routine reads its own figure.
# --------------------------------------------------------------------------
_XE1 = [1.0, 1.4, 1.8, 2.2, 2.4, 0.0, 0.1, 0.3, 0.5]
_YE1 = [1.68, 2.77, 3.38, 3.60, 3.62,
        1.05, 2.11, 2.78, 3.00, 3.00,
        0.00, 1.20, 1.91, 2.23, 2.30,
        0.00, 0.20, 1.10, 1.45, 1.51]
_XE2 = [0.0, 1.0, 2.0, 3.0, 4.0, 4.0, 5.0]
_YE2 = [0.00, 0.92, 1.97, 3.00, 4.10,
        0.00, 1.63, 3.36, 5.20, 7.30]
_XE3 = [0.00, 5.50, 0.25, 0.50]
_YE3 = [0.0, 5.5, 0.0, 5.0]
_XE4 = [0.0, 0.6, 1.0, 1.5, 6.0, 0.0, .15, .30]
_YE4 = [0.00, .120, .200, .300, 1.20,
        .040, .040, .180, .280, 1.14,
        .040, .040, .120, .260, .960]


def calculate_sdwc(x: float, y: float, z: float,
                   beta_aspect: float) -> Dict[str, object]:
    """Translate SDWC: Figure 4.4.1-76C downwash, taper 0.25 and 0.50.

    The first two stages are shared; the third splits on taper through two
    TBFUNX curves of Figure 4.4.1-76C3, and the fourth is common again.

    Args:
        x: Streamwise station.
        y: Spanwise station.
        z: Vertical station.
        beta_aspect: ``AB``, the compressible aspect-ratio parameter.

    Returns:
        Dictionary with ``sdw``, a two-element array for taper 0.25 and
        0.50, and the shared intermediate stages.
    """
    stage1 = interx(2, _XC1, [x, z], [7, 4], _YC1, lind=7,
                    lx1l=2, lx2l=0, lx1u=1, lx2u=2)
    stage2 = interx(2, _XC2, [stage1, beta_aspect], [2, 6], _YC2, lind=6,
                    lx1l=0, lx2l=2, lx1u=1, lx2u=2)
    results = []
    for curve in (_YC3A, _YC3B):
        stage3 = tbfunx(_XC3, curve, stage2, lower=0, upper=1)[0]
        results.append(interx(2, _XC4, [stage3, y], [2, 2], _YC4, lind=2,
                              lx1l=0, lx2l=1, lx1u=1, lx2u=1))
    return {
        'sdw': np.array(results),
        'stage1': float(stage1),
        'stage2': float(stage2),
        'taper': (0.25, 0.50),
        'method': 'legacy_sdwc',
    }


def calculate_sdwd(x: float, y: float, z: float,
                   beta_aspect: float) -> Dict[str, object]:
    """Translate SDWD: Figure 4.4.1-76D downwash, taper 0.00, 0.25 and 0.50.

    Args:
        x: Streamwise station.
        y: Spanwise station.
        z: Vertical station.
        beta_aspect: ``AB``, the compressible aspect-ratio parameter.

    Returns:
        Dictionary with ``sdw``, a three-element array for taper 0.00, 0.25
        and 0.50, and the shared intermediate stages.
    """
    stage1 = interx(2, _XD1, [x, z], [5, 4], _YD1, lind=5,
                    lx1l=2, lx2l=0, lx1u=1, lx2u=2)
    stage2 = interx(2, _XD2, [stage1, beta_aspect], [5, 4], _YD2, lind=5,
                    lx1l=0, lx2l=2, lx1u=1, lx2u=2)
    results = []
    for curve in (_YD3A, _YD3B, _YD3C):
        stage3 = tbfunx(_XD3, curve, stage2, lower=0, upper=1)[0]
        results.append(interx(2, _XD4, [stage3, y], [5, 3], _YD4, lind=5,
                              lx1l=0, lx2l=1, lx1u=1, lx2u=1))
    return {
        'sdw': np.array(results),
        'stage1': float(stage1),
        'stage2': float(stage2),
        'taper': (0.00, 0.25, 0.50),
        'method': 'legacy_sdwd',
    }


def calculate_sdwe(x: float, y: float, z: float, beta_aspect: float,
                   taper_ratio: float) -> Dict[str, float]:
    """Translate SDWE: Figure 4.4.1-76E downwash at a general taper.

    Four chained stages, with taper entering at the third rather than
    selecting between separate curves.

    Args:
        x: Streamwise station.
        y: Spanwise station.
        z: Vertical station.
        beta_aspect: ``AB``, the compressible aspect-ratio parameter.
        taper_ratio: Wing taper ratio.

    Returns:
        Dictionary with ``sdw`` and the three intermediate stage values.
    """
    stage1 = interx(2, _XE1, [x, z], [5, 4], _YE1, lind=5,
                    lx1l=2, lx2l=0, lx1u=1, lx2u=2)
    stage2 = interx(2, _XE2, [stage1, beta_aspect], [5, 2], _YE2, lind=5,
                    lx1l=0, lx2l=2, lx1u=1, lx2u=2)
    stage3 = interx(2, _XE3, [stage2, taper_ratio], [2, 2], _YE3, lind=2,
                    lx1l=0, lx2l=1, lx1u=1, lx2u=1)
    sdw = interx(2, _XE4, [stage3, y], [5, 3], _YE4, lind=5,
                 lx1l=0, lx2l=1, lx1u=1, lx2u=1)
    return {
        'sdw': float(sdw),
        'stage1': float(stage1),
        'stage2': float(stage2),
        'stage3': float(stage3),
        'method': 'legacy_sdwe',
    }
