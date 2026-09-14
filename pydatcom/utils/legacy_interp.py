"""
INTERX and the resampling routines built directly on it.

``INTERX`` is a thin dispatcher over the dimension-specific routines:
``TLIN1X`` for one independent variable, ``TLINEX`` for two, ``TLIN3X`` for
three.  The four-variable branch calling ``TLIN4X`` is commented out in the
source with the note "TLIN4X CALL DELETED TO SAVE CORE", so a four-variable
request is an error in the original as well.

52 routines call ``INTERX``, at 256 one-variable, 294 two-variable and 68
three-variable call sites, so this is the single highest-leverage utility in
the legacy codebase.

``EQSPC1`` and ``EQSPCE`` resample body station data onto an equally spaced
X grid using ``INTERX``, then take the slope of the resampled curve with
``TBFUNX``.  ``BODYRT`` and ``BODOPT`` both call them as setup.

Reference: datcom-legacy/datcom_2000/interx.f, tlin1x.f, tlinex.f,
tlin3x.f, eqspc1.f, eqspce.f
"""

import numpy as np
from typing import Dict, Sequence
import logging

from pydatcom.utils.legacy_tables import tlin1x, tlinex, tlin3x
from pydatcom.utils.legacy_numeric import tbfunx

logger = logging.getLogger(__name__)


def interx(n_independent: int, table, var: Sequence[float],
           length: Sequence[int], dep, lind: int = None,
           lx1l: int = 0, lx2l: int = 0, lx3l: int = 0, lx4l: int = 0,
           lx1u: int = 0, lx2u: int = 0, lx3u: int = 0,
           lx4u: int = 0) -> float:
    """Translate INTERX: linear interpolation in up to three variables.

    The dependent table is stored with the *first* independent variable
    varying fastest, as the source header spells out::

        DEP(1) = F(X1,Y1,Z1)   DEP(2) = F(X2,Y1,Z1)
        DEP(5) = F(X1,Y2,Z1)   DEP(13)= F(X1,Y1,Z2)

    which is FORTRAN column-major order over ``DEP(NX1,NX2,NX3)``.

    Note that ``INTERX`` swaps the first two variables when it calls
    ``TLINEX``: it passes ``TABLE(1,2)`` as that routine's ``X1`` and
    ``TABLE(1,1)`` as its ``X2``, with the query values and end modes
    swapped to match.  That swap is preserved here.

    The three-variable table uses the same first-variable-fastest
    ordering.  ``TLIN3X`` declares ``Y(NX1,NX2,NX3)`` but only ever
    touches ``Y`` by handing slices to ``TLINEX``, which declares
    ``Y(NX2,NX1)``, so that declaration is dead and the effective layout
    cannot be settled by reading the source.  It was established by
    compiling and running the legacy routine: with ``Y(i)=i`` over an
    asymmetric 2x3x2 table the returned flat indices fit
    ``j + (i-1)*NX2 + (k-1)*NX1*NX2``, so X2 varies fastest.  Composed
    with INTERX's own swap this reduces to the ordering documented above.

    Args:
        n_independent: ``NIND``, the number of independent variables.
        table: Independent variable tables.  A flat array laid out as the
            source ``TABLE(LIND,4)``, in which case ``lind`` gives the
            column stride; a 2-D array with one column per variable; or a
            sequence of per-variable grids.
        var: Query value for each independent variable.
        length: Used length of each independent variable table.
        dep: Dependent variable values in the source ordering above.
        lind: ``LIND``, the declared first dimension of ``TABLE``.  Source
            call sites pack every variable's grid into one flat array, so
            variable ``k`` starts at offset ``k*lind``.  Required for a flat
            table with more than one independent variable.
        lx1l, lx2l, lx3l, lx4l: Lower end modes per variable.
        lx1u, lx2u, lx3u, lx4u: Upper end modes per variable.

    Returns:
        The interpolated value, the source's ``ANS``.

    Raises:
        ValueError: For an unsupported ``NIND`` or an inconsistent table.
    """
    if n_independent not in (1, 2, 3, 4):
        raise ValueError("INTERX supports 1 to 4 independent variables")
    if n_independent == 4:
        # The source comments out the TLIN4X call entirely.
        raise ValueError(
            "INTERX has no four-variable branch; the source deleted the "
            "TLIN4X call to save core")
    grids = _columns(table, n_independent, length, lind)
    values = np.asarray(dep, dtype=float).ravel()
    query = [float(v) for v in var]

    if n_independent == 1:
        size = int(length[0])
        if values.size < size:
            raise ValueError("DEP is shorter than the independent table")
        return float(tlin1x(grids[0], values[:size], query[0], lx1l, lx1u))

    if n_independent == 2:
        n1, n2 = int(length[0]), int(length[1])
        if values.size < n1 * n2:
            raise ValueError(
                f"DEP needs {n1 * n2} values for a {n1}x{n2} table")
        # DEP(NX1,NX2) in column-major order; grid[i, j] = F(x1_i, x2_j).
        grid = values[:n1 * n2].reshape((n1, n2), order='F')
        # The source swap: TLINEX's X1 is the second variable.
        return float(tlinex(grids[1], grids[0], grid,
                            query[1], query[0],
                            lx2l, lx1l, lx2u, lx1u))

    n1, n2, n3 = (int(length[0]), int(length[1]), int(length[2]))
    if values.size < n1 * n2 * n3:
        raise ValueError(
            f"DEP needs {n1 * n2 * n3} values for a {n1}x{n2}x{n3} table")
    # DEP(NX1,NX2,NX3) in column-major order, first variable fastest.  That
    # is exactly TLIN3X's (len(x2), len(x1), len(x3)) once INTERX's own
    # first/second variable swap is applied.
    cube = values[:n1 * n2 * n3].reshape((n1, n2, n3), order='F')
    return float(tlin3x(grids[1], grids[0], grids[2], cube,
                        query[1], query[0], query[2],
                        lx2l, lx1l, lx3l, lx2u, lx1u, lx3u))


def _columns(table, n_independent: int, length: Sequence[int],
             lind: int = None):
    """Independent variable grids, trimmed to their used lengths.

    Source call sites declare one flat ``DATA`` array holding every
    variable's grid end to end, and let INTERX's ``TABLE(LIND,4)`` dummy
    argument slice it into columns of stride ``LIND``.  A trailing sentinel
    such as 999999. often pads the last column.
    """
    if len(length) < n_independent:
        raise ValueError("LENGTH must cover every independent variable")
    try:
        array = np.asarray(table, dtype=float)
    except (ValueError, TypeError):
        # A ragged sequence of per-variable grids, which numpy cannot make
        # rectangular.  Each entry is taken as one variable's grid.
        array = None
    if array is None or array.dtype == object:
        return [np.asarray(table[index], dtype=float)[:int(length[index])]
                for index in range(n_independent)]
    flat = array.ndim == 1

    if flat and n_independent > 1 and lind is None:
        raise ValueError(
            "a flat TABLE with more than one independent variable needs "
            "LIND, the declared column stride")

    grids = []
    for index in range(n_independent):
        size = int(length[index])
        if size < 1:
            raise ValueError(f"variable {index + 1} has a nonpositive length")
        if flat:
            start = index * int(lind) if lind else 0
            column = array[start:start + size]
        elif array.ndim == 2:
            # Source layout TABLE(LIND,4): one column per variable.
            if array.shape[1] <= index:
                raise ValueError(
                    f"TABLE has no column for variable {index + 1}")
            column = array[:, index]
        else:
            column = np.asarray(table[index], dtype=float)
        if column.size < size:
            raise ValueError(
                f"variable {index + 1} table is shorter than LENGTH "
                f"({column.size} < {size})")
        grids.append(np.asarray(column[:size], dtype=float))
    return grids


def _equal_spaced_stations(x, n_equal: int):
    """The equally spaced X grid EQSPC1 and EQSPCE build.

    The source accumulates ``XE(I)=XE(I-1)+XIN`` rather than forming
    ``X(1)+(I-1)*XIN``, and pins both endpoints to the original ones.  The
    accumulation is preserved: it is what the legacy routines actually do,
    and the two differ in the last bits.
    """
    x = np.asarray(x, dtype=float)
    count = int(n_equal)
    if count < 2:
        raise ValueError("EQSPC1/EQSPCE need at least two output stations")
    stations = np.empty(count)
    step = (x[-1] - x[0]) / float(count - 1)
    stations[0] = x[0]
    for index in range(1, count - 1):
        stations[index] = stations[index - 1] + step
    stations[count - 1] = x[-1]
    return stations


def _resample(x, values, stations):
    """Interior values via INTERX, endpoints taken from the source arrays."""
    x = np.asarray(x, dtype=float)
    values = np.asarray(values, dtype=float)
    count = len(stations)
    out = np.empty(count)
    out[0] = values[0]
    out[count - 1] = values[-1]
    for index in range(1, count - 1):
        out[index] = interx(1, x, [stations[index]], [len(x)], values)
    return out


def eqspc1(x, s, n_equal: int) -> Dict[str, np.ndarray]:
    """Translate EQSPC1: resample one curve onto equally spaced stations.

    Args:
        x: Original station coordinates, length ``NP``.
        s: Values at those stations.
        n_equal: ``NE``, the number of equally spaced output stations.

    Returns:
        Dictionary with ``xe``, ``se`` and ``dsedx``.

    Raises:
        ValueError: If the inputs are empty, mismatched or fewer than two
            output stations are requested.

    Notes:
        ``DSEDX`` comes from ``TBFUNX`` evaluated on the *resampled* table,
        so it is a local quadratic slope of ``(XE, SE)``, not the derivative
        of the piecewise-linear interpolant that produced ``SE``.  A single
        input station is the source's degenerate branch: every output
        station takes that value and a zero slope.
    """
    x = np.asarray(x, dtype=float)
    s = np.asarray(s, dtype=float)
    if x.ndim != 1 or x.shape != s.shape or len(x) == 0:
        raise ValueError("EQSPC1 requires nonempty matching one-dimensional arrays")
    count = int(n_equal)
    if count < 2:
        raise ValueError("EQSPC1 needs at least two output stations")

    if len(x) == 1:
        return {
            'xe': np.full(count, x[0]),
            'se': np.full(count, s[0]),
            'dsedx': np.zeros(count),
        }

    stations = _equal_spaced_stations(x, count)
    resampled = _resample(x, s, stations)
    slopes = np.array([tbfunx(stations, resampled, station, lower=0, upper=0)[1]
                       for station in stations])
    return {'xe': stations, 'se': resampled, 'dsedx': slopes}


def eqspce(x, r, p, s, n_equal: int) -> Dict[str, np.ndarray]:
    """Translate EQSPCE: resample three curves onto equally spaced stations.

    The four-array body form of ``EQSPC1``: ``R``, ``P`` and ``S`` are each
    resampled over the same new X grid, and only ``S`` gets a slope.

    Args:
        x: Original station coordinates, length ``NP``.
        r: First dependent array, normally body radius.
        p: Second dependent array, normally perimeter.
        s: Third dependent array, normally cross-sectional area.
        n_equal: ``NE``, the number of equally spaced output stations.

    Returns:
        Dictionary with ``xe``, ``re``, ``pe``, ``se`` and ``dsedx``.

    Raises:
        ValueError: If the inputs are empty, mismatched or fewer than two
            output stations are requested.
    """
    arrays = [np.asarray(value, dtype=float) for value in (x, r, p, s)]
    x, r, p, s = arrays
    if x.ndim != 1 or len(x) == 0 or any(a.shape != x.shape for a in arrays):
        raise ValueError("EQSPCE requires nonempty matching one-dimensional arrays")
    count = int(n_equal)
    if count < 2:
        raise ValueError("EQSPCE needs at least two output stations")

    if len(x) == 1:
        return {
            'xe': np.full(count, x[0]), 're': np.full(count, r[0]),
            'pe': np.full(count, p[0]), 'se': np.full(count, s[0]),
            'dsedx': np.zeros(count),
        }

    stations = _equal_spaced_stations(x, count)
    resampled = {name: _resample(x, values, stations)
                 for name, values in (('re', r), ('pe', p), ('se', s))}
    resampled['dsedx'] = np.array([
        tbfunx(stations, resampled['se'], station, lower=0, upper=0)[1]
        for station in stations])
    resampled['xe'] = stations
    return resampled
