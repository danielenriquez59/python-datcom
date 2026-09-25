"""
Namelist card syntax checks: VNAME, LVALUE and RVALUE.

TESTOR checks each namelist card before it is read.  VNAME matches a
variable name (and an array subscript) at column ``L``; LVALUE and RVALUE
scan the logical and real values after it, counting the values (``NDML``)
and the syntax faults (``NF``) they find.

A card is a string, one character per source word; column ``L`` is
``card[L-1]``.  The source reads past column 80 in several places, into
whatever follows the card in the caller's buffer; here a column beyond the
string reads as ``'\\0'``, a character that matches nothing, so a scan
that runs off the card stops as it would at any unexpected character.

Reference: datcom-legacy/datcom_2000/vname.f, lvalue.f, rvalue.f
"""

from typing import Sequence, Tuple

from pydatcom.utils.constants import KAND

# TESTOR's NUMBER table: the digits, then + - . * E.
NUMBER = '0123456789+-.*E'
BLANK, EQUAL, COMMA = ' ', '=', ','


def _col(card: str, k: int) -> str:
    return card[k - 1] if 1 <= k <= len(card) else '\0'


def _nmtest(card: str, l: int, key: str) -> bool:
    return all(_col(card, l + i) == key[i] for i in range(len(key)))


def vname(card: str, l: int, names: Sequence[str], nf: int,
          array: bool = False, ndms: int = 0
          ) -> Tuple[int, int, bool, bool, int, int]:
    """Translate VNAME: match a variable name, and its subscript.

    Args:
        card: The card.  l: The column to start at.  names: The candidate
            names, in the source's order (``NAMES``/``LEN``).  nf: The fault
            count so far.  array, ndms: ``ARRAY`` and ``NDMS`` as they stood;
            an unmatched name leaves them.

    Returns:
        ``(l, i, found, array, ndms, nf)``: the column after the name (and
        subscript), the 1-based index of the name (0 if none), whether it
        was found, whether it was subscripted, the subscript less one
        (``NDMS``), and the fault count.

    Raises:
        ValueError: when a malformed subscript is followed by no ``=`` on
            the card; the source scans on through memory.
    """
    for i, name in enumerate(names, 1):
        if _nmtest(card, l, name) and \
                _col(card, l + len(name)) in (BLANK, EQUAL, '('):
            break
    else:
        return l, 0, False, array, ndms, nf
    l += len(name)
    ndms = 0
    if _col(card, l) != '(':
        return l, i, True, False, ndms, nf
    l += 1
    while True:
        if l >= 79:
            if _col(card, l) != ')':
                nf += 1
            break
        c = _col(card, l)
        if c in NUMBER[:10]:
            ndms = 10 * ndms + NUMBER.index(c)
            l += 1
            continue
        if c == ')':
            break
        nf += 1
        while True:
            l += 1
            if l > len(card):
                raise ValueError("VNAME: no '=' after a bad subscript")
            if _col(card, l) == EQUAL:
                return l, i, True, True, ndms, nf
    l += 1
    ndms -= 1
    if ndms < 0:
        nf += 1
        ndms = 0
    return l, i, True, True, ndms, nf


def lvalue(card: str, l: int, ndml: int, nf: int) -> Tuple[int, int, int]:
    """Translate LVALUE: scan logical values from column ``l``.

    Returns ``(l, ndml, nf)``: where the scan stopped, the value count and
    the fault count.  Kept as executed: at columns 80 and on only the digit
    0 is recognised as a repeat count, and each blank after a value counts
    as a fault.
    """
    while True:
        mult, star = 0, False
        if l >= 81:
            break
        if _col(card, l) == BLANK:
            l += 1
            continue
        while True:
            j = 0
            for jj in range(1, 15):
                if 11 <= jj <= 13:
                    continue
                if _col(card, l) == NUMBER[jj - 1]:
                    j = jj
                    break
                if l >= 80:
                    break
            if not j:
                break
            if star:
                nf += 1
            if j == 14:
                star = True
            if j <= 10:
                mult = 10 * mult + j - 1
            l += 1
        if not star and mult > 0:
            nf += 1
        if star and mult == 0:
            nf += 1
        if mult == 0:
            mult = 1
        if _nmtest(card, l, '.TRUE.'):
            l += 6
        elif _nmtest(card, l, '.FALSE.'):
            l += 7
        else:
            break
        while l < 80 and _col(card, l) == BLANK:
            nf += 1
            l += 1
        if _col(card, l) not in (COMMA, KAND):
            nf += 1
        if _col(card, l) == COMMA:
            l += 1
        ndml += mult
    if ndml == 0:
        nf += 1
    return l, ndml, nf


def rvalue(card: str, l: int, ndml: int, nf: int) -> Tuple[int, int, int]:
    """Translate RVALUE: scan real values from column ``l``.

    Returns ``(l, ndml, nf)``.  Kept as executed: a value without a decimal
    point or exponent counts as a fault, as does each blank inside the
    value list, and the scan does not stop at column 80.
    """
    while True:
        end = sign = exp = star = dec = False
        mult = powr = 0
        while True:
            first = l
            if l >= 81:
                return l, ndml, nf + (ndml == 0)
            if _col(card, l) != BLANK:
                break
            l += 1
        while True:
            c = _col(card, l)
            j = NUMBER.index(c) + 1 if c in NUMBER else 16
            if j <= 10:
                if not star:
                    mult = 10 * mult + j - 1
                if exp:
                    powr = 10 * powr + j - 1
            elif j <= 12:
                if sign:
                    nf += 1
                if not (l == first or _col(card, l - 1) in 'E*'):
                    nf += 1
                sign = True
            elif j <= 13:
                if dec:
                    nf += 1
                if exp:
                    nf += 1
                dec = True
            elif j <= 14:
                nf += star + sign + exp + dec
                if mult == 0:
                    nf += 1
                    mult = 1
                star = True
            elif j <= 15 and l != first:
                if not dec:
                    nf += 1
                if exp:
                    nf += 1
                exp, sign = True, False
            elif c == BLANK:
                nf += 1
            elif c == KAND and l == first:
                return l, ndml, nf + (ndml == 0)
            elif c not in (COMMA, KAND):
                return l, ndml, nf + (ndml == 0)
            else:
                end = True
            if l == first and exp:
                return l, ndml, nf + (ndml == 0)
            if end:
                break
            l += 1
        ndml += 1
        if star:
            ndml = ndml - 1 + mult
        if not dec and not exp:
            nf += 1
        if exp and powr == 0:
            nf += 1
        if _col(card, l) == KAND:
            return l, ndml, nf + (ndml == 0)
        l += 1
