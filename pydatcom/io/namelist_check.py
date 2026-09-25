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

from typing import List, Sequence, Tuple

from pydatcom.io.fortran_format import fortran_write

from pydatcom.utils.constants import KAND

# TESTOR's NUMBER table: the digits, then + - . * E.
NUMBER = '0123456789+-.*E'
BLANK, EQUAL, COMMA = ' ', '=', ','


def _col(card: str, k: int) -> str:
    return card[k - 1] if 1 <= k <= len(card) else '\0'


def _nmtest(card: str, l: int, key: str) -> bool:
    return all(_col(card, l + offset) == key[offset]
               for offset in range(len(key)))


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
    for name_index, name in enumerate(names, 1):
        if _nmtest(card, l, name) and \
                _col(card, l + len(name)) in (BLANK, EQUAL, '('):
            break
    else:
        return l, 0, False, array, ndms, nf
    l += len(name)
    ndms = 0
    if _col(card, l) != '(':
        return l, name_index, True, False, ndms, nf
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
                return l, name_index, True, True, ndms, nf
    l += 1
    ndms -= 1
    if ndms < 0:
        nf += 1
        ndms = 0
    return l, name_index, True, True, ndms, nf


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


_F1170 = '(1X,80A1,3X,33H** ERROR ** UNKNOWN NAMELIST NAME)'
_F1180 = '(1X,80A1)'
_F1190 = '(2X,A1,3HEND,78X,32H** MISSING NAMELIST TERMINATION ,8HADDED **)'
_F1200 = '(1X,A1,3HEND,75A1)'
_F1210 = '(1X,80A1,3X,39H** ERROR ** NO NAMELIST NAME FOLLOWING ,A2)'
_F1220 = ('(1X,80A1,3X,11H** ERROR **,I3,2H*A,I3,2H*B,I3,2H*C,'
          'I3,2H*D,I3,2H*E,I3,2H*F)')


def testor(card: str, l: int, nam: int, k: int, ier: int,
           names: Sequence[str], ldm: Sequence[int], state: dict) -> dict:
    """Translate TESTOR: check one namelist card and echo it.

    Args:
        card: The card (80 columns; any columns past 80 are what the
            source would read beyond it).  l: The column to start at.
        nam: The namelist-name count (below 2: no namelist).  k: 1 for
            the namelist's first card, 2 for a continuation.  ier: nonzero
            for an unknown namelist name.
        names, ldm: The namelist's variable names and dimensions
            (negative for LOGICAL).
        state: TESTOR's saved locals ``end`` (initially true), ``ndml``,
            ``ndm`` and ``is``; updated in place.

    Returns:
        ``l``, ``unit6`` and ``unit11`` (the records written) and the six
        fault counts ``na`` (unknown names) to ``nf`` (value syntax).

    Notes:
        Kept as executed: the ``.TRUE.``/``.FALSE.`` skips compare with the
        undeclared ``TRUE`` and ``FALSE``, zero words, and never match; a
        continuation card that starts with a number adds its values to the
        previous variable's count.
    """
    st = state
    st.setdefault('end', True)
    st.setdefault('ndml', 0)
    st.setdefault('ndm', 0)
    st.setdefault('is', 0)
    kol = [_col(card, column) for column in range(1, 81)]
    out6: List[str] = []
    out11: List[str] = []
    n = {'na': 0, 'nb': 0, 'nc': 0, 'nd': 0, 'ne': 0, 'nf': 0}
    ndms = 0
    i = 0
    array = False

    def result():
        out11.extend(fortran_write('(80A1)', kol))
        return dict(n, l=l, unit6=out6, unit11=out11)

    if ier != 0:
        out6 += fortran_write(_F1170, kol)
        return result()
    if nam < 2:
        out6 += fortran_write(_F1180, kol)
        st['end'] = True
        return result()
    if k == 1 and not st['end']:
        out6 += fortran_write(_F1190, [KAND])
        out11 += fortran_write(_F1200, [KAND] + [' '] * 75)
        st['end'] = True
    if k == 2 and st['end']:
        # KAND is a 4-character Hollerith word; A2 prints '$ '.
        out6 += fortran_write(_F1210, kol + [KAND.ljust(4)])
        return result()
    label = 1040
    while True:
        if label == 1040:
            if k == 1:
                st['end'] = False
            if _col(card, l) == BLANK:
                l += 1
                label = 1150 if l >= 81 else 1040
                continue
            label = 1050
        elif label == 1050:
            if _col(card, l) == KAND:
                st['end'] = True
            if st['end']:
                label = 1150
                continue
            l, i, found, array, ndms, n['nf'] = vname(
                card, l, names, n['nf'], array, ndms)
            if found:
                label = 1090
                continue
            if k != 1:
                ndms = st['ndml']
                st['ndml'] = 0
                if _col(card, l) in NUMBER[:14]:
                    label = 1120
                    continue
            st['is'] = i
            n['na'] += 1
            label = 1075
        elif label == 1075:
            if l >= 81:
                label = 1150
                continue
            if _col(card, l) == KAND:
                label = 1050
                continue
            if _col(card, l) in (BLANK, EQUAL, COMMA):
                label = 1077
                continue
            l += 1
        elif label == 1077:
            if l >= 81:
                label = 1150
                continue
            c = _col(card, l)
            if c in NUMBER or c in (BLANK, EQUAL, COMMA):
                l += 1
                continue
            if _col(card, l - 1) == NUMBER[14]:
                l -= 1
            label = 1040
        elif label == 1090:
            st['is'] = i
            st['ndm'] = abs(ldm[i - 1])
            st['ndml'] = 0
            if array and st['ndm'] == 1:
                array = False
                n['nc'] += 1
            while l < 80 and _col(card, l) == BLANK:
                l += 1
            if _col(card, l) != EQUAL:
                n['nb'] += 1
            if _col(card, l) == EQUAL:
                l += 1
            label = 1120
        elif label == 1120:
            i = st['is']
            if i == 0:
                label = 1075
                continue
            if l >= 81:
                label = 1150
                continue
            if _col(card, l) == BLANK:
                l += 1
                continue
            if ldm[i - 1] < 0:
                l, st['ndml'], n['nf'] = lvalue(card, l, st['ndml'],
                                                n['nf'])
            if ldm[i - 1] > 0:
                l, st['ndml'], n['nf'] = rvalue(card, l, st['ndml'],
                                                n['nf'])
            if st['ndml'] == 0:
                label = 1077
                continue
            if st['ndm'] == 1 and st['ndml'] > 1:
                n['nd'] += 1
            st['ndml'] = st['ndml'] + ndms
            if st['ndm'] > 1 and st['ndml'] > st['ndm']:
                n['ne'] += 1
            label = 1040 if l <= 80 else 1150
        elif label == 1150:
            nerr = sum(n.values())
            if nerr == 0:
                out6 += fortran_write(_F1180, kol)
            else:
                out6 += fortran_write(_F1220, kol + [n['na'], n['nb'],
                                                     n['nc'], n['nd'],
                                                     n['ne'], n['nf']])
            return result()


# NMLIST's routing: namelist number -> the check table it uses.
NMLIST_TABLES = {1: 'FLTCON', 2: 'FLTCON', 3: 'OPTINS', 4: 'BODY',
                 5: 'PLNF', 6: 'SCHR', 7: 'SYNTHS', 8: 'PLNF', 9: 'SCHR',
                 10: 'PLNF', 11: 'SCHR', 12: 'PRPOWR', 13: 'JETPWR',
                 14: 'LARWB', 15: 'GRNDEF', 16: 'TVTPAN', 17: 'EXPR',
                 18: 'SYMFLP', 19: 'ASYFLP', 20: 'HYPEFF', 21: 'TRNJET',
                 22: 'CONTAB', 23: 'SCHR', 24: 'PLNF'}


def nmlist(card: str, l: int, name: int, k: int, ier: int,
           tables: dict, state: dict) -> dict:
    """Translate NMLIST: check a card against namelist ``name``'s table.

    ``tables`` maps the table names of :data:`NMLIST_TABLES` to
    ``(names, ldm)``.  An unknown namelist (``ier`` 1) is routed as
    namelist 1.  Returns TESTOR's result, with ``name`` as used.
    """
    if ier == 1:
        name = 1
    names, ldm = tables[NMLIST_TABLES[name]]
    result = testor(card, l, name, k, ier, names, ldm, state)
    return dict(result, name=name)
