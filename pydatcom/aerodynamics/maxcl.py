"""
MAXCL: section maximum lift coefficient at flight conditions.

A base value from the leading-edge sharpness parameter and the
maximum-thickness station, plus three increments:

- ``DEL1`` for camber, a three-variable lookup
- ``DEL2`` for a maximum-thickness station away from 0.30, camber only
- ``DEL3`` for Reynolds number

This is the first translated routine to exercise INTERX's three-variable
branch, whose memory layout was settled against a compiled probe of the
legacy TLIN3X.

Reference: datcom-legacy/datcom_2000/maxcl.f
"""

import numpy as np
from typing import Dict, Optional, Sequence
import logging

from pydatcom.utils.legacy_interp import interx
from pydatcom.utils.legacy_numeric import tbfunx

logger = logging.getLogger(__name__)

# Base CLMAX against (DELTAY, XOVC).  LIND=13, LGH=(13,4).
_PARM58 = [0., 1., 1.1, 1.25, 1.50, 2.00, 2.25, 2.50, 3.00, 3.50, 4.00,
           4.50, 5.00, .3, .35, .4, .45, 0., 0., 0., 0., 0., 0., 0., 0., 0.]
_EVAL58 = [.8, .8, .81, .85, .98, 1.2, 1.31, 1.43, 1.58, 1.59, 1.55, 1.47, 1.42,
           .8, .8, .81, .85, .98, 1.2, 1.31, 1.43, 1.51, 1.52, 1.48, 1.41, 1.41,
           .8, .8, .81, .85, .98, 1.2, 1.31, 1.39, 1.44, 1.45, 1.43, 1.35, 1.35,
           .8, .8, .81, .85, .98, 1.2, 1.29, 1.33, 1.35, 1.35, 1.35, 1.35, 1.35]

# Camber increment against (DELTAY, COVC, CXVC).  LIND=10, LGH=(10,4,4).
# The third grid starts at offset 2*LIND = 20.
_PARM59 = [0., 1., 1.5, 2., 2.5, 3., 3.5, 4., 4.5, 5.,
           0., .02, .04, .06, 0., 0., 0., 0., 0., 0.,
           .15, .30, .40, .50, 0., 0., 0., 0., 0., 0.]
_EVAL59 = [
    0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
    0., .12, .23, .32, .23, .12, .07, .07, .06, .03,
    0., .27, .36, .40, .30, .16, .1, .1, .1, .06,
    0., .27, .36, .40, .30, .16, .1, .1, .1, .06,
    0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
    0., 0., .17, .22, .25, .15, .08, .06, .08, .02,
    0., .15, .27, .40, .33, .18, .08, .06, .07, .01,
    0., .53, .64, .56, .38, .22, .10, .06, .06, 0.,
    0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
    0., 0., .05, .19, .23, .14, .07, .06, .06, .02,
    0., .16, .27, .37, .34, .19, .10, .07, .09, .05,
    0., .39, .50, .50, .40, .25, .13, .10, .13, .09,
    0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
    0., 0., .04, .10, .15, .14, .09, .07, .10, .08,
    0., .10, .19, .27, .30, .22, .15, .13, .18, .12,
    0., .26, .36, .43, .45, .30, .21, .20, .23, .20,
]

# Thickness-station increment against (DELTAY, XOVC).  LIND=10, LGH=(10,3).
_PARM60 = [0., .25, 1., 1.75, 2., 2.75, 3., 3.5, 4., 4.5,
           .35, .4, .45, 0., 0., 0., 0., 0., 0., 0.]
_EVAL60 = [0., .16, .19, .20, .18, .08, .05, -.01, -.03, -.04,
           0., .16, .19, .20, .18, .08, .05, .03, .04, .02,
           0., .16, .19, .20, .18, .08, .08, .10, .12, .08]

# Reynolds increment against (DELTAY, RN).  LIND=8, LGH=(8,4).
_PARM61 = [0., 1., 1.5, 2., 2.5, 3., 3.5, 4.,
           3.e06, 6.e06, 9.e06, 2.5e07, 0., 0., 0., 0.]
_EVAL61 = [0., -.09, -.13, -.13, -.09, -.11, -.20, -.24,
           0., 0., -.08, -.07, -.03, -.03, -.07, -.10,
           0., 0., 0., 0., 0., 0., 0., 0.,
           0., .12, .14, .07, .01, .11, .18, .20]

# A fourth increment the source computes and never uses; see below.
_PARM62 = np.array([0., 1., 1.5, 2., 2.5, 3., 3.5, 4.])
_EVAL62 = np.array([0., 0., -.19, -.37, -.46, -.5, -.5, -.5])

# The source snaps a near-0.30 thickness station onto exactly 0.30 before
# testing it, so DEL2 is skipped for a station that only nearly matches.
_XOVC_SNAP = 0.30
_XOVC_SNAP_TOLERANCE = 1.0e-5


def calculate_maxcl(delta_y: float, thickness_station: float,
                    reynolds_per_length: Sequence[float], chord: float,
                    camber: bool = False,
                    camber_ratio: float = 0.0,
                    camber_station: float = 0.0) -> Dict[str, object]:
    """Translate MAXCL: section CLMAX over a Mach schedule.

    ``CLMAX = CLBASE + DEL1 + DEL2 + DEL3``

    Args:
        delta_y: Leading-edge sharpness parameter ``DELTAY``.
        thickness_station: ``XOVC``, the maximum-thickness chord station.
        reynolds_per_length: Reynolds number per unit length at each Mach
            condition; the source reads ``FLC(I+42)``.
        chord: Reference chord used to form the section Reynolds number.
        camber: Whether the section is cambered, enabling DEL1 and DEL2.
        camber_ratio: ``COVC``, used by DEL1.
        camber_station: ``CXVC``, used by DEL1.

    Returns:
        Dictionary with ``clmax`` per Mach condition and the base and
        increment values behind it.

    Raises:
        ValueError: If the chord is nonpositive or the schedule is empty.

    Notes:
        The source also evaluates a fourth increment ``DEL4`` from Figure
        PARM62/EVAL62 inside the Mach loop and then never adds it to
        ``CLMAX0``.  It is computed and reported here as ``del4_unused`` so
        the dead path stays visible, but like the source it does not enter
        the result.

        ``EVAL59`` is declared with 160 elements while the call passes
        ``LDEP=200``.  Only 10*4*4 = 160 are read, so the mismatch is
        harmless; the array is carried at its true length.
    """
    reynolds_per_length = np.atleast_1d(
        np.asarray(reynolds_per_length, dtype=float))
    if reynolds_per_length.size == 0:
        raise ValueError("MAXCL needs at least one flight condition")
    if chord <= 0.0:
        raise ValueError("MAXCL requires a positive chord")

    station = float(thickness_station)
    if abs(station - _XOVC_SNAP) <= _XOVC_SNAP_TOLERANCE:
        station = _XOVC_SNAP

    base = interx(2, _PARM58, [delta_y, station], [13, 4], _EVAL58, lind=13)

    del1 = 0.0
    del2 = 0.0
    if camber:
        del1 = interx(3, _PARM59, [delta_y, camber_ratio, camber_station],
                      [10, 4, 4], _EVAL59, lind=10)
        if station != _XOVC_SNAP:
            del2 = interx(2, _PARM60, [delta_y, station], [10, 3],
                          _EVAL60, lind=10)

    clmax = np.empty(reynolds_per_length.size)
    del3_values = np.empty(reynolds_per_length.size)
    for index, per_length in enumerate(reynolds_per_length):
        reynolds = per_length * chord
        # The source substitutes a nominal Reynolds number below unity.
        lookup = 9.0e6 if reynolds < 1.0 else reynolds
        del3 = interx(2, _PARM61, [delta_y, lookup], [8, 4], _EVAL61, lind=8)
        del3_values[index] = del3
        clmax[index] = base + del1 + del2 + del3

    del4_unused = float(tbfunx(_PARM62, _EVAL62, float(delta_y),
                               lower=0, upper=0)[0])

    return {
        'clmax': clmax,
        'clbase': float(base),
        'del1': float(del1),
        'del2': float(del2),
        'del3': del3_values,
        'del4_unused': del4_unused,
        'thickness_station_used': station,
        'method': 'legacy_maxcl',
    }
