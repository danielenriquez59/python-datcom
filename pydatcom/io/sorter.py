"""
Integer table sort: SORTER.

SORTER orders the rows of an ``N`` by ``M`` INTEGER array ``A`` on column
``ICOL`` and then its columns on row ``IROW``, each ascending, by repeated
minimum search: each chosen entry is overwritten with ``BMAX+BMAX`` so it
is not chosen again.  That sentinel is formed in 32-bit INTEGER arithmetic
and is reproduced here with its wraparound.

Reference: datcom-legacy/datcom_2000/sorter.f
"""

from typing import List, Optional, Sequence, Tuple


def _i32(x: int) -> int:
    return (x + 2 ** 31) % 2 ** 32 - 2 ** 31


def _order(values: Sequence[int]) -> List[int]:
    """The source's selection: the 1-based positions in ascending order,
    the lowest position first among equals."""
    b = list(values)
    n = len(b)
    bmax = b[0]
    for scan_index in range(n, 0, -1):
        if bmax <= b[scan_index - 1]:
            bmax = b[scan_index - 1]
    sentinel = _i32(bmax + bmax)
    order = []
    # The first search starts from the first entry, later ones from the
    # sentinel; KK carries over when no entry equals the minimum.
    bmin, kk = b[0], 1
    for _ in range(n):
        for scan_index in range(n, 0, -1):
            if bmin >= b[scan_index - 1]:
                bmin = b[scan_index - 1]
            if bmin == b[scan_index - 1]:
                kk = scan_index
        b[kk - 1] = sentinel
        order.append(kk)
        bmin = sentinel
    return order


def sorter(a: Sequence[Sequence[int]], irow: int, icol: int,
           iflag: Optional[int] = None
           ) -> Tuple[List[List[int]], List[List[int]], List[int], List[int],
                      Optional[int]]:
    """Translate SORTER.

    Args:
        a: The array, ``N`` rows of ``M`` INTEGER values.  irow, icol: The
            row and column to sort on; zero or less skips that sort.
            iflag: ``IFLAG`` as it stood, returned unchanged for a 1 by 1
            array, which the source leaves alone.

    Returns:
        ``(a, b, ir, ic, iflag)``: ``A`` (rewritten with the row-sorted
        table when both sorts run), ``B`` (the sorted table; ``None`` when
        neither sort runs, since the source leaves it unset), the 1-based
        row and column orders, and ``IFLAG``: 1 when the order is
        unchanged, 0 otherwise.
    """
    a = [[_i32(int(v)) for v in row] for row in a]
    n = len(a)
    m = len(a[0]) if n else 0
    if n <= 0 or m <= 0:
        return a, None, [], [], 1
    ir, ic = list(range(1, n + 1)), list(range(1, m + 1))
    if n == m and n < 2:
        return a, None, ir, ic, iflag
    b = None
    if icol > 0:
        ir = _order([row[icol - 1] for row in a])
        b = [list(a[ll - 1]) for ll in ir]
        if irow > 0:
            a = [list(row) for row in b]
    if irow > 0:
        ic = _order(a[irow - 1])
        b = [[row[ll - 1] for ll in ic] for row in a]
    unchanged = ir == list(range(1, n + 1)) and ic == list(range(1, m + 1))
    return a, b, ir, ic, 1 if unchanged else 0
