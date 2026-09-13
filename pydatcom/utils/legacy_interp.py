"""
INTERX: the general table interpolator used throughout Digital DATCOM.

``INTERX`` is a thin dispatcher over the dimension-specific routines:
``TLIN1X`` for one independent variable, ``TLINEX`` for two, ``TLIN3X`` for
three.  The four-variable branch calling ``TLIN4X`` is commented out in the
source with the note "TLIN4X CALL DELETED TO SAVE CORE", so a four-variable
request is an error in the original as well.

52 routines call ``INTERX``, at 256 one-variable, 294 two-variable and 68
three-variable call sites, so this is the single highest-leverage utility in
the legacy codebase.

Reference: datcom-legacy/datcom_2000/interx.f, tlin1x.f, tlinex.f, tlin3x.f
"""

import numpy as np
from typing import Sequence
import logging

from pydatcom.utils.legacy_tables import tlin1x, tlinex

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
        NotImplementedError: For ``NIND`` of 3, pending TLIN3X; see below.
    """
    if n_independent not in (1, 2, 3, 4):
        raise ValueError("INTERX supports 1 to 4 independent variables")
    if n_independent == 4:
        # The source comments out the TLIN4X call entirely.
        raise ValueError(
            "INTERX has no four-variable branch; the source deleted the "
            "TLIN4X call to save core")
    if n_independent == 3:
        # Reported before the table is validated: an unsupported branch is
        # more useful to hear about than a shape complaint about its inputs.
        raise NotImplementedError(
            "INTERX's three-variable branch needs TLIN3X, which is not yet "
            "translated. Its source declares Y(NX1,NX2,NX3) but hands slices "
            "to TLINEX, which declares Y(NX2,NX1); those layouts disagree "
            "unless NX1 == NX2, so the ordering must be settled against a "
            "compiled run before the branch can be trusted.")

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

    raise AssertionError("unreachable: NIND is validated above")


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
    array = np.asarray(table, dtype=float)
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
