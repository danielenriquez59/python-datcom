"""
Namelist reader internals: FINDCH, SKIPBL, EXTRST, TODEC, FINDVN, the
value printers FORINT, FORLOG and FORREA, the assignment helpers SUBINT,
SUBLOG and SUBREA, and the card reader READCD.

Cards are sequences of single characters, column ``k`` at index ``k-1``.

Reference: datcom-legacy/datcom_2000/findch.f, skipbl.f, extrst.f,
todec.f, findvn.f, forint.f, forlog.f, forrea.f, subint.f, sublog.f,
subrea.f, readcd.f
"""

import struct
from typing import IO, List, Sequence, Tuple

from pydatcom.io.fortran_format import fortran_write


def _i32(x: int) -> int:
    return (x + 2 ** 31) % 2 ** 32 - 2 ** 31


def _powi(x: float, m: int) -> float:
    """libgcc's ``__powidf2``, which gfortran uses for ``X**N`` with an
    INTEGER ``N``: square-and-multiply, then a reciprocal for ``N < 0``."""
    n = -m if m < 0 else m
    y = x if n % 2 else 1.0
    while n > 1:
        n >>= 1
        x = x * x
        if n % 2:
            y = y * x
    return 1.0 / y if m < 0 else y


def findch(kol: Sequence[str], char: str, kcol: int) -> int:
    """Translate FINDCH: the first column from ``kcol`` holding ``char``,
    or 81."""
    while kol[kcol - 1] != char:
        kcol += 1
        if kcol >= 81:
            break
    return kcol


def skipbl(kol: Sequence[str], kcol: int) -> int:
    """Translate SKIPBL: the first non-blank column from ``kcol``, or 81."""
    while kol[kcol - 1] == ' ':
        kcol += 1
        if kcol >= 81:
            break
    return kcol


def extrst(kol: Sequence[str], first: int, last: int) -> List[str]:
    """Translate EXTRST: columns ``first`` to ``last``."""
    return list(kol[first - 1:last])


def todec(alpha: Sequence[str]) -> Tuple[float, int]:
    """Translate TODEC: decode a number from 80 characters.

    Returns ``(ANS, IERR)``; ``IERR`` is 1 for a stray character in the
    mantissa.  Blanks are skipped in the mantissa; the exponent after ``E``
    stops at its first character that is not a sign or digit.  The digit
    accumulators are 32-bit INTEGERs and the powers of ten are formed as
    gfortran forms ``10.**N``.
    """
    ierr, frac = 0, False
    lvalue = rvalue = power = pten = signp = 0
    signv = 0.0
    icol = 0
    exponent = False
    for column in range(1, 81):
        icol = column
        c = alpha[column - 1]
        if c == ' ':
            continue
        if c == 'E':
            exponent = True
            break
        if c in '+-':
            if signv != 0.:
                break
            signv = 1. if c == '+' else -1.
            continue
        if c == '.':
            frac = True
            continue
        if not c.isdigit() or len(c) != 1:
            ierr = 1
            break
        if frac:
            rvalue = _i32(10 * rvalue + int(c))
            pten += 1
        else:
            lvalue = _i32(10 * lvalue + int(c))
    if exponent and icol + 1 <= 80:
        for column in range(icol + 1, 81):
            c = alpha[column - 1]
            if c in '+-':
                if signp != 0:
                    break
                signp = 1 if c == '+' else -1
                continue
            if len(c) == 1 and c.isdigit():
                power = _i32(10 * power + int(c))
                continue
            break
    if signv == 0.:
        signv = 1.
    if signp == 0:
        signp = 1
    ans = signv * (float(lvalue) + float(rvalue) * _powi(10., -pten)) * \
        _powi(10., signp * power)
    return ans, ierr


def findvn(lenvn: Sequence[int], iname: Sequence[str],
           vname: Sequence[str]) -> Tuple[int, bool]:
    """Translate FINDVN: find the name in columns 1-80 of ``iname`` among
    the names packed in ``vname`` with lengths ``lenvn``; the rest of the
    80 columns must be blank.  Returns ``(NVN, FOUND)``, ``NVN`` 0 if not
    found."""
    indx, maxc = 1, 0
    for nvn, length in enumerate(lenvn, 1):
        if nvn > 1:
            indx += maxc
        maxc = length
        if any(iname[column - 1] != vname[column + indx - 2]
               for column in range(1, maxc + 1)):
            continue
        if any(iname[column - 1] != ' '
               for column in range(maxc + 1, 81)):
            continue
        return nvn, True
    return 0, False


def forint(name: Sequence[str], values: Sequence[int]) -> List[str]:
    """Translate FORINT: print an INTEGER namelist variable."""
    return fortran_write('(1H0,8A1,8(I14,1H,)/(9X,8(I14,1H,)))',
                         list(name[:8]) + [int(v) for v in values])


def forlog(name: Sequence[str], values: Sequence[bool]) -> List[str]:
    """Translate FORLOG: print a LOGICAL namelist variable."""
    return fortran_write('(1H0,8A1,8(13X,L1,1H,)/(9X,8(13X,L1,1H,)))',
                         list(name[:8]) + [bool(v) for v in values])


def _as_real(v) -> float:
    """A word printed as REAL.  A LOGICAL word (NAMEW passes the whole
    block to FORREA) is its bits read as a single-precision REAL."""
    if isinstance(v, bool):
        return struct.unpack('<f', struct.pack('<i', int(v)))[0]
    return float(v)


def forrea(name: Sequence[str], values: Sequence[float]) -> List[str]:
    """Translate FORREA: print a REAL namelist variable."""
    return fortran_write('(1H0,8A1,8(G14.7,1H,)/(9X,8(G14.7,1H,)))',
                         list(name[:8]) + [_as_real(v) for v in values])


def subint(lans: int) -> int:
    """Translate SUBINT: ``LLOC=LANS`` (the caller stores the result)."""
    return lans


sublog = subint
subrea = subint


def readcd(unit: IO[str], kol: Sequence[str]) -> Tuple[List[str], bool]:
    """Translate READCD: read one 80-column card.

    Returns the card (80 characters, blank-padded) and ``IEOF``.  At end
    of file the source retries once, then returns the card unchanged with
    ``IEOF`` set.
    """
    for _ in range(2):
        line = unit.readline()
        if line:
            text = line.rstrip('\r\n')
            return list(text[:80].ljust(80)), False
    return list(kol), True


def toint(alpha: Sequence[str]) -> Tuple[int, int]:
    """Translate TOINT: decode an integer.

    Returns ``(IANS, IERR)``: the value rounded half away from zero, and
    1 when the field has a decimal point, 2 when TODEC finds a stray
    character (the later test wins).
    """
    ierr = 0
    ier1 = 1 if findch(alpha, '.', 1) < 81 else 0
    ans, ier2 = todec(alpha)
    if ier1 > 0:
        ierr = 1
    if ier2 > 0:
        ierr = 2
    sign = -1.0 if ans < 0. else 1.0
    return int(sign * (abs(ans) + 0.5)), ierr


def namew(kand: str, nlname: Sequence[str], vname: Sequence[str],
          lenvn: Sequence[int], vdime: Sequence[int], comblk: Sequence,
          loc: Sequence[int], vtype: int = 0) -> Tuple[List[str], int]:
    """Translate NAMEW: print a namelist, one variable at a time.

    ``vdime`` gives each variable's length, negative for LOGICAL; ``loc``
    its 1-based position in ``comblk`` (below 1: not printed).  Returns
    the records and ``VTYPE``: a variable with ``vdime`` 0 is printed with
    the type left by the previous one (the saved local), and the INTEGER
    branch is never selected.
    """
    lines = fortran_write('(1H1,8A1)', [kand] + list(nlname))
    indx, maxc = 1, 0
    for nvn, length in enumerate(lenvn, 1):
        if nvn > 1:
            indx += maxc
        maxc = length
        kol = [' '] * 8
        kol[:maxc] = vname[indx - 1:indx - 1 + maxc]
        kol[7] = '='
        kk = loc[nvn - 1]
        if kk < 1:
            continue
        jj = abs(vdime[nvn - 1])
        if vdime[nvn - 1] < 0:
            vtype = 0
        if vdime[nvn - 1] > 0:
            vtype = 2
        values = comblk[kk - 1:kk - 1 + jj]
        if vtype == 0:
            lines += forlog(kol, values)
        if vtype == 1:
            lines += forint(kol, values)
        if vtype == 2:
            lines += forrea(kol, values)
    lines += fortran_write('(1H0,A1,3HEND)', [kand])
    return lines, vtype


def _at(alpha: Sequence[str], k: int) -> str:
    """Column ``k``; past the card reads as a character matching nothing."""
    return alpha[k - 1] if 1 <= k <= len(alpha) else '\0'


def tolog(alpha: Sequence[str]) -> Tuple[bool, int]:
    """Translate TOLOG: decode ``.TRUE.``/``.T.`` or ``.FALSE.``/``.F.``
    after leading blanks.  Returns ``(LANS, IERR)``, ``IERR`` 1 when
    neither is found."""
    icol = skipbl(alpha, 1)
    if icol >= 81 or _at(alpha, icol) != '.':
        return False, 1
    if _at(alpha, icol + 1) == 'T' and _at(alpha, icol + 2) == '.':
        return True, 0
    if all(_at(alpha, icol + 1 + k) == c for k, c in enumerate('TRUE')):
        return (True, 0) if _at(alpha, icol + 5) == '.' else (False, 1)
    if _at(alpha, icol + 1) == 'F' and _at(alpha, icol + 2) == '.':
        return False, 0
    if all(_at(alpha, icol + 1 + k) == c for k, c in enumerate('FALSE')):
        return (False, 0) if _at(alpha, icol + 6) == '.' else (False, 1)
    return False, 1


def reptct(iconst: List[str]) -> Tuple[int, int]:
    """Translate REPTCT: split a repeat count ``n*`` off a constant.

    ``iconst`` (80 columns) is edited in place: with a ``*`` the text after
    it is shifted to column 1, the last columns keeping their old
    characters.  Returns ``(IREPT, IERR)``: the count (1 without one, or
    when it decodes as 0) and TOINT's error code.
    """
    ierr, irept = 0, 1
    icol = findch(iconst, '*', 1)
    if icol >= 81:
        return irept, ierr
    icol -= 1
    itcon = extrst(iconst, 1, icol) + [' '] * (80 - icol)
    irr, ierr = toint(itcon)
    if irr != 0:
        irept = irr
    icol += 2
    for j, i in enumerate(range(icol, 81), 1):
        iconst[j - 1] = iconst[i - 1]
    return irept, ierr
