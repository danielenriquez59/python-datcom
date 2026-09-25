"""
Small output routines: DMPARY, PRCSID, SWRITE and EXIT.

Each returns the records it writes to unit 6, rendered through
:func:`pydatcom.io.fortran_format.fortran_write` from the same FORMAT the
source uses (assembled at run time in DMPARY and SWRITE).  Column 1 is the
carriage-control character, written as is.

Reference: datcom-legacy/datcom_2000/dmpary.f, prcsid.f, swrite.f, exit.f
"""

import struct
from typing import IO, Iterable, List, Sequence, Tuple, Union

from pydatcom.io.fortran_format import fortran_write
from pydatcom.utils.constants import UNUSED

Value = Union[float, str]

_DMPARY_A = ['(5(5', '(5(4', '(5(3', '(5(2']
_DMPARY_B = ['X,A1', 'X,A2', 'X,A3', 'X,A4']
_DMPARY_TAIL = [',1H(', ',I3,', '2H)=', '1P  ', ',E12', '.5))']


def dmpary(array: Sequence[float], name: str, nlet: int) -> List[str]:
    """Translate DMPARY: dump ``NAME(I)=value``, five to a line.

    ``name`` is the Hollerith word (up to four characters); ``nlet`` is
    how many of its characters to print.  A blank line comes first.
    """
    fmt = ''.join([_DMPARY_A[nlet - 1], _DMPARY_B[nlet - 1]] + _DMPARY_TAIL)
    items = []
    for word_index, value in enumerate(array, 1):
        items += [name.ljust(4), word_index, float(value)]
    return fortran_write('(1H )') + fortran_write(fmt, items)


def prcsid(head: bool, idcse: Sequence[str]) -> List[str]:
    """Translate PRCSID: print the case title centred on a 132 line.

    ``idcse`` holds the 74 ``/CASEID/`` words, one character each.
    Nothing is printed without ``HEAD`` or for a blank title.
    """
    if not head:
        return []
    chars = [c if c else ' ' for c in idcse]
    marked = [column for column, ch in enumerate(chars, 1) if ch != ' ']
    if not marked:
        return []
    ifb, ilb = marked[0], marked[-1]
    first = 65 - (ilb - ifb + 1) // 2
    ipr = [' '] * (first - 1) + chars[ifb - 1:ilb]
    return fortran_write('(1X,132A1)', ipr)


_FORSUB = [',2X ', ',A4,', '1X  ', ',2X ', ',A4,', '2X  ', ',3X ', ',A4,',
           '2X  ', ',3X ', ',A4,', '3X  ', ',4X ', ',A4,', '3X  ', ',4X ',
           ',A4,', '4X  ', ',5X ', ',A4,', '4X  ']


def swrite(last: int, form: Sequence[str], icont: Sequence[int],
           columns: Sequence[Sequence[float]], nrpeat: int,
           ndmf: bool = False, naf: bool = False
           ) -> Tuple[List[str], bool, bool]:
    """Translate SWRITE: print rows of up to 14 columns, showing
    ``UNUSED`` as ``NDM``, ``-UNUSED`` as blank and ``2*UNUSED`` as ``NA``.

    Args:
        last: The number of columns.  form: The caller's FORMAT as its
            CHARACTER*4 words; column ``I`` occupies words ``4I-1`` to
            ``4I+1`` (1-based), which are swapped for an ``A4`` field of
            the width ``icont[I]`` while that column prints a marker.
        columns: The column arrays.  nrpeat: The number of rows.
        ndmf, naf: The flags as they stood; set when ``NDM`` or ``NA`` is
            printed.

    Returns:
        The records, and ``NDMF`` and ``NAF``.
    """
    words = list(form)
    ippat = [0] * 15
    lines: List[str] = []
    for row in range(nrpeat):
        z: List[Value] = [0.0] * 15
        for column_index in range(1, last + 1):
            z[column_index] = float(columns[column_index - 1][row])
        for column_index in range(1, last + 1):
            if z[column_index] == UNUSED:
                icpat, z[column_index], ndmf = 1, 'NDM ', True
            elif z[column_index] == -UNUSED:
                icpat, z[column_index] = 1, '    '
            elif z[column_index] == 2 * UNUSED:
                icpat, z[column_index], naf = 1, 'NA  ', True
            else:
                icpat = 0
            if icpat == ippat[column_index]:
                continue
            k = 4 * column_index - 1
            if ippat[column_index] == 0:
                j2 = 3 * (icont[column_index - 1] - 6) - 2
                words[k - 1:k + 2] = _FORSUB[j2 - 1:j2 + 2]
            else:
                words[k - 1:k + 2] = form[k - 1:k + 2]
            ippat[column_index] = icpat
        lines += fortran_write(''.join(words), z[1:last + 1])
    return lines, ndmf, naf


def exit_units(units: Iterable[IO]) -> None:
    """Translate EXIT: close the run's files (units 5, 6, 8-14)."""
    for handle in units:
        handle.close()


def _int_as_real(v: int) -> float:
    """An INTEGER word read as a single-precision REAL (MESSGE's RL/ML
    EQUIVALENCE, in the source's word size)."""
    return struct.unpack('<f', struct.pack('<i', int(v)))[0]


def messge(rout: Sequence[str], mess: Sequence[str],
           x: Sequence[Sequence[float]], msscl: Sequence, nf: int,
           iovly: int) -> List[str]:
    """Translate MESSGE: the extrapolation diagnostic written to unit 12.

    Args:
        rout: The calling routine's name, two Hollerith words.
        mess: The message words (Hollerith).  x: The (up to four)
            independent-variable tables ``X1``-``X4``.
        msscl: The 21 message-control words (1-based list): 1-2 Hollerith,
            3 ``NSTQ`` (message words), 4, 5 ``NSTP`` (variables), and per
            variable ``L`` words ``4L+2`` (a code), ``4L+3`` (the index of
            the last table entry), ``4L+4`` and ``4L+5``.
        nf: ``NF`` from ``/OVERLY/``; negative suppresses the message.
        iovly: The overlay number.

    Returns:
        The records.  Kept as executed: the codes in ``MSSCL(4)`` and
        ``MSSCL(4L+2)`` share storage with REAL words through an
        EQUIVALENCE and are printed as REALs, so their bits show as tiny
        denormal numbers.
    """
    if nf < 0:
        return []
    m = msscl
    nstq, nstp = int(m[3]), int(m[5])
    rl = [0.0] * 14
    rl[1] = _int_as_real(m[4])
    for step in range(1, min(nstp, 4) + 1):
        rl[step + 1] = _int_as_real(m[4 * step + 2])
        table = x[step - 1]
        rl[step + 5] = float(table[0])
        rl[step + 9] = float(table[int(m[4 * step + 3]) - 1])
    lines = fortran_write('(3I3)', [iovly, nstq, nstp])
    lines += fortran_write('(24A4)', [m[1], m[2], rout[0], rout[1]] +
                           list(mess[:nstq]))
    for step in range(1, nstp + 1):
        lines += fortran_write('(3E12.5,2I2)',
                               [rl[step + 1], rl[step + 5], rl[step + 9],
                                int(m[4 * step + 4]), int(m[4 * step + 5])])
    lines += fortran_write('(E12.5)', [rl[1]])
    return lines
