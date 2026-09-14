"""
PTCP: control-surface pressure ratio and centre of pressure.

For one Mach-cone region of a deflected control surface the routine builds a
pressure distribution ``PP`` along a normalised station, integrates it twice
-- once for area and once for its first moment -- and reports the running
pressure ratio ``PSTR`` and the centre-of-pressure parameter ``TCP`` at a
requested station.

There are four paths, chosen by the region number and by whether the region
sits at the tip or the root:

=========================  ==========  =======  ==================
region (``CASENO``)         location    points   what varies
=========================  ==========  =======  ==================
1, 2, 7, 8                  tip             19   leading-edge cone
1, 2, 7, 8                  root            19   leading-edge cone
3, 4, 5, 6                  tip             10   a second generator
3, 4, 5, 6                  root            10   a second generator
=========================  ==========  =======  ==================

Regions 3 and 4 take the trailing-edge sweep as that second generator;
regions 5 and 6 take the hinge-line sweep.  The ten-point paths stop at a
station of 1.0, while the nineteen-point paths continue to 10.0 in unit
steps.

The four paths do not share their quadrature weights.  Each has its own
``AREA`` set, its own ``AMON`` moment set, and its own expression for
``TCP``, and the source reaches them through a dense sequence of label
jumps in which the same weight blocks are entered from several directions.

Reference: datcom-legacy/datcom_2000/ptcp.f
"""

import numpy as np
from typing import Dict
import logging

from pydatcom.utils.legacy_numeric import arccos, tbfunx

logger = logging.getLogger(__name__)

# The source's TOR selector.
TIP = 1.0
ROOT = 2.0

# Regions strictly between 2 and 7 use the ten-point path, the source's
# ONLY10 = CASENO.GT.2 .AND. CASENO.LT.7.
_TEN_POINT_REGIONS = (3, 4, 5, 6)

# Regions 3 and 4 take the trailing-edge generator, 5 and 6 the hinge line.
_TRAILING_EDGE_REGIONS = (3, 4)

# The station grid: tenths to 1.0, then unit steps to 10.0.
_STATIONS = np.array([.1, .2, .3, .4, .5, .6, .7, .8, .9, 1.,
                      2., 3., 4., 5., 6., 7., 8., 9., 10.])


def _area_first(pp: np.ndarray, root: bool) -> float:
    """AREA(1), whose weights differ between the tip and the root."""
    if root:
        return .130786 * pp[0] - .045340 * pp[1]
    return .113301 * pp[0] - .032975 * pp[1]


def _area_two_to_four(pp: np.ndarray) -> list:
    """AREA(2) to AREA(4), shared by every path at label 1140."""
    return [
        .037500 * pp[0] + .079167 * pp[1] - .020833 * pp[2] + .004167 * pp[3],
        -.004167 * pp[0] + .054167 * pp[1] + .054167 * pp[2] - .004167 * pp[3],
        .004167 * pp[0] - .020833 * pp[1] + .079167 * pp[2] + .037500 * pp[3],
    ]


def _area_five_to_ten_staggered(pp: np.ndarray) -> list:
    """AREA(5) to AREA(10) on the ten-point tip path.

    This is the block the source falls into when neither of label 1140's
    two tests fires, and its weights are not those of label 1200.
    """
    return [
        .041667 * pp[3] + .066667 * pp[4] - .008333 * pp[5],
        -.008333 * pp[3] + .066667 * pp[4] + .041667 * pp[5],
        .037500 * pp[5] + .079167 * pp[6] - .020833 * pp[7] + .004167 * pp[8],
        -.004167 * pp[5] + .054167 * pp[6] + .054167 * pp[7] - .004167 * pp[8],
        .004167 * pp[5] - .020833 * pp[6] + .079167 * pp[7] + .037500 * pp[8],
        -.038494 * pp[7] + .121105 * pp[8] + .017389 * pp[9],
    ]


def _area_five_to_ten_regular(pp: np.ndarray) -> list:
    """AREA(5) to AREA(10) at label 1200, two clean four-point panels."""
    return [
        .037500 * pp[3] + .079167 * pp[4] - .020833 * pp[5] + .004167 * pp[6],
        -.004167 * pp[3] + .054167 * pp[4] + .054167 * pp[5] - .004167 * pp[6],
        .004167 * pp[3] - .020833 * pp[4] + .079167 * pp[5] + .037500 * pp[6],
        .037500 * pp[6] + .079167 * pp[7] - .020833 * pp[8] + .004167 * pp[9],
        -.004167 * pp[6] + .054167 * pp[7] + .054167 * pp[8] - .004167 * pp[9],
        .004167 * pp[6] - .020833 * pp[7] + .079167 * pp[8] + .037500 * pp[9],
    ]


def _area_eleven_to_nineteen(pp: np.ndarray) -> list:
    """AREA(11) to AREA(19), the unit-step tail at label 1020."""
    out = []
    for base in (9, 12, 15):
        out.extend([
            .375000 * pp[base] + .791667 * pp[base + 1] -
            .208333 * pp[base + 2] + .041667 * pp[base + 3],
            -.041667 * pp[base] + .541667 * pp[base + 1] +
            .541667 * pp[base + 2] - .041667 * pp[base + 3],
            .041667 * pp[base] - .208333 * pp[base + 1] +
            .791667 * pp[base + 2] + .375000 * pp[base + 3],
        ])
    return out


def _moment_ten_point_a(pp: np.ndarray) -> list:
    """AMON(1) to AMON(10) reached by falling through label 1140."""
    return [
        .007917 * pp[0] - .004167 * pp[1] + .001250 * pp[2],
        .005417 * pp[0] + .010833 * pp[1] - .001250 * pp[2],
        -.002083 * pp[0] + .015833 * pp[1] + .011250 * pp[2],
        .011250 * pp[2] + .031667 * pp[3] - .010416 * pp[4] + .002500 * pp[5],
        -.001250 * pp[2] + .021667 * pp[3] + .027083 * pp[4] - .002500 * pp[5],
        .001250 * pp[2] - .008333 * pp[3] + .039583 * pp[4] + .022500 * pp[5],
        .022500 * pp[5] + .055417 * pp[6] - .016667 * pp[7] + .003750 * pp[8],
        -.002500 * pp[5] + .037917 * pp[6] + .043334 * pp[7] - .003750 * pp[8],
        .002500 * pp[5] - .014583 * pp[6] + .063334 * pp[7] + .033750 * pp[8],
        -.030795 * pp[7] + .108995 * pp[8] + .017389 * pp[9],
    ]


def _moment_ten_point_b(pp: np.ndarray) -> list:
    """AMON(1) to AMON(10) at label 1220."""
    return [
        .006667 * pp[0] - .001667 * pp[1],
        .006667 * pp[0] + .008333 * pp[1],
        .007500 * pp[1] + .023750 * pp[2] - .008333 * pp[3] + .002083 * pp[4],
        -.000833 * pp[1] + .016250 * pp[2] + .021667 * pp[3] - .002083 * pp[4],
        .000833 * pp[1] - .006250 * pp[2] + .031667 * pp[3] + .018750 * pp[4],
        .018750 * pp[4] + .047500 * pp[5] - .014583 * pp[6] + .003333 * pp[7],
        -.002083 * pp[4] + .032500 * pp[5] + .037917 * pp[6] - .003333 * pp[7],
        .002083 * pp[4] - .012500 * pp[5] + .055417 * pp[6] + .030000 * pp[7],
        .033333 * pp[7] + .060000 * pp[8] - .008333 * pp[9],
        -.006667 * pp[7] + .060000 * pp[8] + .041667 * pp[9],
    ]


def _moment_eleven_to_nineteen(pp: np.ndarray) -> list:
    """AMON(11) to AMON(19) at label 1040."""
    return [
        .37500 * pp[9] + 1.583333 * pp[10] - .625000 * pp[11] +
        .166667 * pp[12],
        -.041667 * pp[9] + 1.083333 * pp[10] + 1.625000 * pp[11] -
        .166667 * pp[12],
        .041667 * pp[9] - .416667 * pp[10] + 2.375000 * pp[11] +
        1.500000 * pp[12],
        1.5000 * pp[12] + 3.958333 * pp[13] - 1.250000 * pp[14] +
        .291667 * pp[15],
        -.166667 * pp[12] + 2.708333 * pp[13] + 3.250000 * pp[14] -
        .291667 * pp[15],
        .166667 * pp[12] - 1.041667 * pp[13] + 4.750000 * pp[14] +
        2.625 * pp[15],
        2.625 * pp[15] + 6.333333 * pp[16] - 1.875000 * pp[17] +
        .416667 * pp[18],
        -.291667 * pp[15] + 4.333333 * pp[16] + 4.875000 * pp[17] -
        .416667 * pp[18],
        .291667 * pp[15] - 1.666667 * pp[16] + 7.12500 * pp[17] +
        3.75 * pp[18],
    ]


def calculate_ptcp(station: float, region: int, location: float,
                   tan_le: float, beta: float,
                   tan_te: float = 0.0,
                   tan_hinge_line: float = 0.0) -> Dict[str, object]:
    """Translate PTCP: pressure ratio and centre of pressure for one region.

    Args:
        station: ``RSP``, the station at which both outputs are read.
        region: ``CASENO``, the Mach-cone region number 1 to 8.
        location: ``TOR``; :data:`TIP` or :data:`ROOT`.
        tan_le: Leading-edge sweep tangent, the source's ``A(62)``.
        beta: Compressibility parameter, the source's ``SPR(1)``.
        tan_te: Trailing-edge sweep tangent, ``A(80)``.  Used by regions 3
            and 4.
        tan_hinge_line: Hinge-line sweep tangent, ``SPR(14)``.  Used by
            regions 5 and 6.

    Returns:
        Dictionary with ``pressure_ratio`` and ``centre_of_pressure``, the
        source's ``PSTROT`` and ``TCPOUT``, together with the distributions
        behind them and which path ran.

    Raises:
        ValueError: If beta is zero, or the region is outside 1 to 8.

    Notes:
        The ten-point tip and root paths differ in the sign of the
        generator term in the centre-of-pressure denominator: the tip forms
        ``(SUMA-SUMM)/(SUMA + TGEN*SUMM)`` and the root
        ``(SUMA-SUMM)/(SUMA - TGEN*SUMM)``.  The nineteen-point paths use
        ``SUMA/(SUMA+SUMM)`` and involve no generator at all.

        The source computes ``SUMA(11)=SUMA(10)+AREA(11)`` and then
        immediately recomputes the same element as the first pass of its
        ``DO 1030 K=11,19`` loop.  Harmless, and not reproduced.

        Every quadrature element in the routine integrates exactly what it
        should -- each ``AREA`` panel summing to its station width and each
        ``AMON`` panel to that interval's integral of ``r dr`` -- with two
        exceptions.  ``AREA(1)`` is a deliberate endpoint rule for the
        leading-edge singularity and sums to neither width.  And
        ``AMON(10)`` on the ten-point tip path sums to 0.095589 where its
        interval's first moment is 0.095, a 0.6 percent discrepancy that
        ``AREA(10)`` of the same block and ``AMON(10)`` of the other set
        both avoid.  Since its last coefficient is shared with
        ``AREA(10)``, the pair looks derived together for the endpoint, so
        this is recorded rather than labelled a defect.  Both are pinned by
        test.
    """
    if beta == 0.0:
        raise ValueError("PTCP divides the sweep tangents by beta")
    region = int(region)
    if not 1 <= region <= 8:
        raise ValueError(f"PTCP region {region} is outside 1 to 8")
    tan_le_beta = tan_le / beta
    root = float(location) == ROOT
    ten_point = region in _TEN_POINT_REGIONS

    if ten_point:
        generator = (tan_te if region in _TRAILING_EDGE_REGIONS
                     else tan_hinge_line) / beta
        if root:
            # Label 1180.
            squared = tan_le_beta**2
            numerator = 2.0 * (1.0 - squared)
            pp = np.array([
                arccos(numerator /
                       (1.0 - squared *
                        ((1.0 - n) / (1.0 - generator * n))**2) - 1.0) / np.pi
                for n in np.arange(1, 11) * 0.1])
        else:
            # The ten-point tip block after label 1110.
            one_plus = 1.0 + tan_le_beta
            two_plus = 2.0 + tan_le_beta + generator
            difference = tan_le_beta - generator
            pp = np.array([
                arccos((one_plus - two_plus * n) /
                       (one_plus - difference * n)) / np.pi
                for n in np.arange(1, 11) * 0.1])
    else:
        generator = 0.0
        if root:
            # Label 1070.
            squared = tan_le_beta**2
            one_minus_two = 1.0 - 2.0 * squared
            pp = np.array([
                arccos((squared + one_minus_two * (1.0 + r)**2) /
                       ((1.0 + r)**2 - squared)) / np.pi
                for r in _STATIONS])
        else:
            # The nineteen-point tip block after label 1000.
            one_plus = 1.0 + tan_le_beta
            pp = np.array([arccos((one_plus - r) / (one_plus + r)) / np.pi
                           for r in _STATIONS])

    area = [_area_first(pp, root)] + _area_two_to_four(pp)
    if ten_point and not root:
        # Falls through label 1140 into the staggered block.
        area += _area_five_to_ten_staggered(pp)
        moment = _moment_ten_point_a(pp)
        path = 'ten_point_tip'
    else:
        # Label 1200.
        area += _area_five_to_ten_regular(pp)
        moment = _moment_ten_point_b(pp)
        path = 'ten_point_root' if ten_point else (
            'nineteen_point_root' if root else 'nineteen_point_tip')
    if not ten_point:
        area += _area_eleven_to_nineteen(pp)
        moment += _moment_eleven_to_nineteen(pp)

    cumulative_area = np.cumsum(np.array(area))
    cumulative_moment = np.cumsum(np.array(moment))
    count = 10 if ten_point else 19
    stations = _STATIONS[:count]

    pressure = cumulative_area / stations
    if not ten_point:
        centre = cumulative_area / (cumulative_area + cumulative_moment)
    elif root:
        centre = ((cumulative_area - cumulative_moment) /
                  (cumulative_area - generator * cumulative_moment))
    else:
        centre = ((cumulative_area - cumulative_moment) /
                  (cumulative_area + generator * cumulative_moment))

    pressure_at, _ = tbfunx(stations, pressure, float(station), 0, 0)
    centre_at, _ = tbfunx(stations, centre, float(station), 0, 0)
    return {
        'pressure_ratio': float(pressure_at),
        'centre_of_pressure': float(centre_at),
        'path': path,
        'points': count,
        'stations': stations,
        'pp': pp,
        'area': np.array(area),
        'moment': np.array(moment),
        'cumulative_area': cumulative_area,
        'cumulative_moment': cumulative_moment,
        'pressure_distribution': pressure,
        'centre_distribution': centre,
        'generator': float(generator),
        'method': 'legacy_ptcp',
    }
