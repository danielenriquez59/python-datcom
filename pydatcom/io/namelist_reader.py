"""
DATCOM's namelist reader: NAMER.

NAMER scans cards for ``$NAME``, then reads ``VAR=value,...`` and
``VAR(n)=value,...`` assignments (with ``k*value`` repeats and value
continuation across cards) into a namelist's storage block, up to the
closing ``$``.  A syntax error prints the card with a caret under the
offending column and stops the run; here that raises
:class:`NamelistInputError` carrying the printed records.

Reference: datcom-legacy/datcom_2000/namer.f
"""

from typing import IO, List, MutableSequence, Sequence

from pydatcom.io.fortran_format import fortran_write
from pydatcom.io.namelist_io import (extrst, findch, findvn, readcd, reptct,
                                     skipbl, todec, toint, tolog)

# The error marker.  datcom_2000/namer.f has lost it (``4H   /`` does not
# compile); the earlier datcom.f has ``4H^``.
CARET = '^'
_HEAD = '26H0*** NAMELIST INPUT ERROR.,'
_ERRORS = {
    1310: _HEAD + '51H  ILLEGAL/INCORRECT SPECIFICATION OF NAMELIST NAME.',
    1320: _HEAD + '36H  NO EQUALS FOLLOWING VARIABLE NAME.',
    1330: _HEAD + '51H  NO CLOSING RIGHT PARENTHESIS IN ARRAY DEFINITION.',
    1350: _HEAD + '32H  VARIABLE NAME NOT IN NAMELIST.',
    1360: _HEAD + ('57H  ARRAY SUBSCRIPT OR NUMBER OF CONSTANTS EXCEEDS '
                   'VARIABLE,32H DIMENSION OR COMMON BLOCK SIZE.'),
    1370: _HEAD + '34H  ILLEGAL/INVALID ARRAY SUBSCRIPT.',
    1380: _HEAD + '48H  REPEAT COUNT EXCEEDS VARIABLE ARRAY DIMENSION.',
    1390: _HEAD + '49H  CONSTANT DOES NOT MATCH TYPE OF INPUT REQUIRED.',
}
_F1420 = '(49H0*** EXECUTION TERMINATING DUE TO NAMELIST ERROR.)'
_DIGITS = '0123456789+-.'


class NamelistInputError(Exception):
    """NAMER's error stop; ``lines`` holds everything it printed."""

    def __init__(self, lines: List[str]):
        super().__init__(lines[-3] if len(lines) >= 3 else 'namelist error')
        self.lines = lines


def _pad(chars: Sequence[str]) -> List[str]:
    return list(chars) + [' '] * (80 - len(chars))


def _findch(kol, char, kcol):
    """FINDCH from a column that may be past the card: a column beyond 80
    matches nothing and ends the search one column on."""
    return findch(kol, char, kcol) if kcol <= 80 else kcol + 1


def _skipbl(kol, kcol):
    return skipbl(kol, kcol) if kcol <= 80 else kcol


def namer(kand: str, unit: IO[str], nlname: Sequence[str],
          vname: Sequence[str], lenvn: Sequence[int],
          vdime: Sequence[int], comblk: MutableSequence, loc: Sequence[int],
          state: dict = None) -> dict:
    """Translate NAMER: read one namelist into ``comblk``.

    Args:
        kand: The delimiter, ``$``.  unit: The card file.
        nlname: The namelist name's characters.  vname, lenvn: The
            variable names, packed, and their lengths.  vdime: Each
            variable's dimension, negative for LOGICAL.  comblk: The
            storage block, 1-based, written in place.  loc: Each
            variable's 1-based position in ``comblk`` (below 1: read and
            checked, not stored).
        state: The saved locals ``vtype``, ``nvn`` and ``found``.

    Returns:
        ``ieof`` (the file ended first) and ``lines`` (records printed:
        none unless an error stops the read).

    Raises:
        NamelistInputError: on a syntax error, as the source's STOP.

    Notes:
        Kept as executed: only the first ``len(nlname)`` characters of a
        namelist name are compared; skipping another namelist resumes at
        the card after the first ``$`` found, which may be the start of the
        wanted namelist; and the INTEGER branches are unreachable (the
        type is LOGICAL or REAL from the sign of the dimension).
    """
    st = state if state is not None else {}
    st.setdefault('vtype', 2)
    st.setdefault('nvn', 0)
    st.setdefault('found', False)
    maxcom = len(comblk) - 1
    lines: List[str] = []
    kol: List[str] = [' '] * 80
    kerr: List[str] = [' '] * 80
    at = lambda k: kol[k - 1] if 1 <= k <= 80 else '\0'  # noqa: E731

    def fail(label, column):
        kerr[column - 1] = CARET
        lines.extend(fortran_write('(' + _ERRORS[label] + ',/,1X,80A1)',
                                   kol))
        lines.extend(fortran_write('(1X,80A1)', kerr))
        lines.extend(fortran_write(_F1420))
        raise NamelistInputError(lines)

    def read():
        nonlocal kol
        kol, eof = readcd(unit, kol)
        return eof

    label = 1000
    iend = search = False
    icol = lcol = iecol = ioff = kname = kcons = irept = 0
    iname = idim = iconst = None
    while True:
        if label == 1000:
            iend = False
            kerr = [' '] * 80
            label = 1020
        elif label == 1020:
            if read():
                return {'ieof': True, 'lines': lines}
            icol = skipbl(kol, 1)
            if icol > 80 or at(icol) != kand:
                continue
            icol += 1
            lcol = iecol = icol
            lcol = _findch(kol, ' ', lcol)
            if lcol == icol:
                fail(1310, iecol)
            inln = _pad(extrst(kol, icol, lcol - 1))
            if any(inln[column] != ch
                   for column, ch in enumerate(nlname)):
                label = 1050
                continue
            icol = lcol
            label = 1060
        elif label == 1050:
            icol = findch(kol, kand, icol)
            if icol <= 80:
                label = 1000
                continue
            if read():
                return {'ieof': True, 'lines': lines}
            icol = 1
        elif label == 1060:
            search, ioff = True, 0
            iname, idim, kerr = [' '] * 80, [' '] * 80, [' '] * 80
            icol = _skipbl(kol, icol)
            if icol >= 81:
                if read():
                    return {'ieof': True, 'lines': lines}
                icol = skipbl(kol, 1)
            if at(icol) == kand or iend:
                return {'ieof': False, 'lines': lines}
            lcol = iecol = icol
            lcol = _findch(kol, '=', lcol)
            if lcol >= 81:
                fail(1320, iecol)
            kname = iecol = icol
            iname = _pad(extrst(kol, icol, lcol - 1))
            kk = findch(iname, '(', 1)
            if kk <= 80:
                ll = findch(iname, ')', kk)
                if ll >= 81:
                    fail(1330, iecol)
                idim = _pad(extrst(iname, kk + 1, ll - 1))
                iecol = kname + kk
                for column in range(kk, ll + 1):
                    iname[column - 1] = ' '
                ioff, ierr = toint(idim)
                if ierr != 0 or ioff < 1:
                    fail(1370, iecol)
                ioff -= 1
            icol = lcol + 1
            label = 1140
        elif label == 1140:
            if at(icol) == kand:
                iend = True
            lcol = kcons = iecol = icol
            lcol = _findch(kol, ',', lcol)
            if lcol >= 81:
                lcol = _findch(kol, kand, icol)
                if lcol < 81:
                    iend = True
            iconst = _pad(extrst(kol, icol, lcol - 1))
            iecol = kcons
            irept, ierr = reptct(iconst)
            if not (ierr == 0 and irept > 0):
                fail(1390, iecol)
            icol = lcol + 1
            iecol = kname
            if search:
                st['nvn'], st['found'] = findvn(lenvn, iname, vname)
                if not st['found']:
                    fail(1350, iecol)
            nvn = st['nvn']
            iecol = kcons
            dim = abs(vdime[nvn - 1])
            if ioff + 1 > dim:
                fail(1370, iecol)
            if vdime[nvn - 1] < 0:
                st['vtype'] = 0
            if vdime[nvn - 1] > 0:
                st['vtype'] = 2
            vtype = st['vtype']
            value, ierr = (tolog(iconst) if vtype == 0 else
                           toint(iconst) if vtype == 1 else todec(iconst))
            if ierr != 0:
                fail(1390, iecol)
            for repeat_index in range(1, irept + 1):
                if repeat_index + ioff > dim:
                    fail(1380, kcons)
                if loc[nvn - 1] < 1:
                    continue
                ii = loc[nvn - 1] + repeat_index + ioff - 1
                if ii > maxcom:
                    fail(1360, iecol)
                comblk[ii] = value
            ioff += irept
            icol = _skipbl(kol, icol)
            if iend:
                return {'ieof': False, 'lines': lines}
            if icol >= 81:
                if read():
                    return {'ieof': True, 'lines': lines}
                icol = skipbl(kol, 1)
            if at(icol) in _DIGITS:
                search = False
                label = 1140
            else:
                label = 1060
