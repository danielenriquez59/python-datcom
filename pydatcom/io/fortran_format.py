"""
A FORTRAN FORMAT interpreter for formatted output, matching gfortran.

DATCOM's printers write through FORMAT statements, some of them assembled
at run time (DMPARY, SWRITE).  :func:`fortran_write` renders a format and
an item list to the records gfortran would write, so the printers can be
translated line for line and checked against the compiled source.

Supported: ``nX``, ``Tn``, ``TLn``, ``TRn``, ``nHtext`` and quoted text,
``/``, ``kP``, ``Aw``, ``Lw``, ``Iw[.m]``, ``Fw.d``, ``Ew.d[Ee]``,
``Gw.d``, repeat counts and nested groups, and format reversion.  Carriage
control is not interpreted: column 1 is written like any other, as gfortran
does.  Character items are ``str``; a Hollerith word is its characters.

Checked against gfortran by ``tools/probes/fortran_format.py``.
"""

import math
from decimal import ROUND_HALF_EVEN, Decimal, localcontext
from typing import Iterable, List, Sequence, Tuple, Union

Item = Union[int, float, str, bool]


class _Node:
    """A parsed descriptor: ``kind`` plus its parameters."""

    def __init__(self, kind, **kw):
        self.kind = kind
        self.__dict__.update(kw)


def _parse(fmt: str) -> List[_Node]:
    s = fmt.strip()
    if not (s.startswith('(') and s.endswith(')')):
        raise ValueError(f'format must be parenthesised: {fmt!r}')
    nodes, pos = _parse_list(s, 1)
    if pos != len(s) - 1:
        raise ValueError(f'unbalanced format: {fmt!r}')
    return nodes


def _number(s, i):
    j = i
    while j < len(s) and (s[j].isdigit() or s[j] == ' '):
        j += 1
    digits = s[i:j].replace(' ', '')
    return (int(digits) if digits else None), j


def _parse_list(s, i):
    nodes = []
    while i < len(s):
        c = s[i]
        if c in ' ,':
            i += 1
            continue
        if c == ')':
            return nodes, i
        if c == '/':
            nodes.append(_Node('/'))
            i += 1
            continue
        if c == ':':
            nodes.append(_Node(':'))
            i += 1
            continue
        if c in '\'"':
            j = i + 1
            text = []
            while True:
                if s[j] == c:
                    if j + 1 < len(s) and s[j + 1] == c:
                        text.append(c)
                        j += 2
                        continue
                    break
                text.append(s[j])
                j += 1
            nodes.append(_Node('lit', text=''.join(text)))
            i = j + 1
            continue
        sign = 1
        if c in '+-':
            sign = -1 if c == '-' else 1
            i += 1
        n, i = _number(s, i)
        c = s[i].upper()
        if c == '(':
            inner, i = _parse_list(s, i + 1)
            nodes.append(_Node('group', count=n or 1, nodes=inner))
            i += 1
            continue
        if c == 'H':
            nodes.append(_Node('lit', text=s[i + 1:i + 1 + n]))
            i += 1 + n
            continue
        if c == 'P':
            nodes.append(_Node('P', k=sign * (n or 0)))
            i += 1
            continue
        if c == 'X':
            nodes.append(_Node('X', n=n or 1))
            i += 1
            continue
        if c == 'T':
            kind = 'T'
            if s[i + 1].upper() in 'LR':
                kind = 'T' + s[i + 1].upper()
                i += 1
            w, i = _number(s, i + 1)
            nodes.append(_Node(kind, n=w))
            continue
        if c in 'AILFEGD':
            w, i = _number(s, i + 1)
            d = m = e = None
            if i < len(s) and s[i] == '.':
                d, i = _number(s, i + 1)
                if c == 'I':
                    m, d = d, None
            if i < len(s) and s[i].upper() == 'E' and c in 'EG':
                e, i = _number(s, i + 1)
            node = _Node(c if c != 'D' else 'E', w=w, d=d, m=m, e=e)
            nodes.extend([node] * (n or 1) if n else [node])
            continue
        raise ValueError(f'unsupported format item at {s[i:]!r}')
    raise ValueError('unterminated format')


def _stars(w):
    return '*' * w


def _fmt_i(v, w, m=None):
    v = int(v)
    digits = str(abs(v))
    if m is not None:
        digits = digits.rjust(m, '0') if not (m == 0 and v == 0) else ''
    text = ('-' if v < 0 else '') + digits
    return _stars(w) if len(text) > w else text.rjust(w)


def _special(v, w):
    if math.isnan(v):
        text = 'NaN'
    else:
        text = ('-' if v < 0 else '+' if w > 8 else '') + \
            ('Infinity' if w >= 9 else 'Inf')
        if v > 0 and w <= 8:
            text = 'Inf'
    return _stars(w) if len(text) > w else text.rjust(w)


def _fmt_f(v, w, d, k=0):
    v = float(v)
    if not math.isfinite(v):
        return _special(v, w)
    # The scale factor shifts the exact decimal value, not the binary one.
    with localcontext() as ctx:
        ctx.prec = 800
        scaled = Decimal(abs(v)).scaleb(k)
        body = format(scaled.quantize(Decimal(1).scaleb(-d),
                                      rounding=ROUND_HALF_EVEN), 'f')
    if d == 0:
        body += '.'
    neg = math.copysign(1.0, v) < 0
    if body.startswith('0') and len(body) > 1 and body[1] == '.':
        if len(body) + neg > w:
            body = body[1:]
    text = ('-' if neg else '') + body
    if len(text) > w:
        return _stars(w)
    return text.rjust(w)


def _exp_part(expo, e):
    if e is None:
        if abs(expo) <= 99:
            return 'E' + ('-' if expo < 0 else '+') + f'{abs(expo):02d}'
        if abs(expo) <= 999:
            return ('-' if expo < 0 else '+') + f'{abs(expo):03d}'
        return None
    if abs(expo) >= 10 ** e:
        return None
    return 'E' + ('-' if expo < 0 else '+') + f'{abs(expo):0{e}d}'


def _fmt_e(v, w, d, k=0, e=None):
    v = float(v)
    if not math.isfinite(v):
        return _special(v, w)
    neg = math.copysign(1.0, v) < 0
    a = abs(v)
    if k <= 0:
        sig = d + k
        if sig <= 0:
            return _stars(w)
    else:
        sig = d + 1
    if a == 0.0:
        digits, expo = '0' * sig, 0
    else:
        mant, ex = f'{a:.{sig - 1}e}'.split('e')
        digits = mant.replace('.', '')
        expo = int(ex) + 1          # value = 0.digits * 10**expo
    # The printed mantissa is the value times 10**k.
    exp_val = expo - k if a != 0.0 else 0
    if k <= 0:
        body = '0.' + '0' * (-k) + digits
    else:
        body = digits[:k] + '.' + digits[k:]
    tail = _exp_part(exp_val, e)
    if tail is None:
        return _stars(w)
    text = ('-' if neg else '') + body + tail
    if len(text) > w and body.startswith('0.'):
        text = ('-' if neg else '') + body[1:] + tail
    return _stars(w) if len(text) > w else text.rjust(w)


def _fmt_g(v, w, d, k=0, e=None):
    v = float(v)
    if not math.isfinite(v):
        return _special(v, w)
    a = abs(v)
    n = 4 if e is None else e + 2
    if a == 0.0:
        return _fmt_f(v, w - n, d - 1) + ' ' * n
    # The F form applies when 0.1 - 0.5*10**-(d+1) <= |v| < 10**d - 0.5.
    if 0.1 - 0.5 * 10.0 ** (-d - 1) <= a < 10.0 ** d - 0.5:
        for integer_digits in range(0, d + 1):
            if a < (10.0 ** integer_digits
                    - 0.5 * 10.0 ** (integer_digits - d)):
                return (_fmt_f(v, w - n, d - integer_digits) + ' ' * n)
    return _fmt_e(v, w, d, k, e)


def _fmt_a(v, w):
    text = v if isinstance(v, str) else str(v)
    if w is None:
        return text
    return text[:w] if len(text) >= w else text.rjust(w)


def _fmt_l(v, w):
    return ('T' if v else 'F').rjust(w or 2)


class _Record:
    def __init__(self):
        self.chars: List[str] = []
        self.pos = 0
        self.high = 0
        self.skip_from = None

    def skip(self, n):
        """X or TR: move right; the gap is blanked if text follows."""
        if self.skip_from is None:
            self.skip_from = self.pos
        self.pos += n

    def tab(self, pos):
        self.skip_from = None
        self.pos = max(0, pos)

    def put(self, text):
        if self.skip_from is not None:
            gap = self.pos - self.skip_from
            self.skip_from, self.pos = None, self.skip_from
            self.put(' ' * gap)
        end = self.pos + len(text)
        if len(self.chars) < end:
            self.chars.extend(' ' * (end - len(self.chars)))
        self.chars[self.pos:end] = list(text)
        self.pos = end
        self.high = max(self.high, end)

    def line(self):
        return ''.join(self.chars[:self.high])


def fortran_write(fmt: str, items: Iterable[Item] = ()) -> List[str]:
    """Render ``items`` through the FORMAT ``fmt``; return the records."""
    nodes = _parse(fmt)
    items = list(items)
    state = {'i': 0, 'k': 0}
    records: List[str] = []
    rec = _Record()
    # The reversion point: the last top-level group, or the whole format.
    groups = [j for j, nd in enumerate(nodes) if nd.kind == 'group']
    revert = groups[-1] if groups else 0

    def run(seq) -> bool:
        """Process ``seq``; False when output must stop."""
        nonlocal rec
        for nd in seq:
            if nd.kind == 'group':
                for _ in range(nd.count):
                    if not run(nd.nodes):
                        return False
                continue
            if nd.kind == 'lit':
                rec.put(nd.text)
            elif nd.kind in ('X', 'TR'):
                rec.skip(nd.n)
            elif nd.kind == 'T':
                rec.tab(nd.n - 1)
            elif nd.kind == 'TL':
                rec.tab(rec.pos - nd.n)
            elif nd.kind == 'P':
                state['k'] = nd.k
            elif nd.kind == '/':
                records.append(rec.line())
                rec = _Record()
            elif nd.kind == ':':
                if state['i'] >= len(items):
                    return False
            else:
                if state['i'] >= len(items):
                    return False
                v = items[state['i']]
                state['i'] += 1
                k = state['k']
                if nd.kind == 'A':
                    rec.put(_fmt_a(v, nd.w))
                elif nd.kind == 'L':
                    rec.put(_fmt_l(v, nd.w))
                elif nd.kind == 'I':
                    rec.put(_fmt_i(v, nd.w, nd.m))
                elif nd.kind == 'F':
                    rec.put(_fmt_f(v, nd.w, nd.d, k))
                elif nd.kind == 'E':
                    rec.put(_fmt_e(v, nd.w, nd.d, k, nd.e))
                elif nd.kind == 'G':
                    rec.put(_fmt_g(v, nd.w, nd.d, k, nd.e))
        return True

    if run(nodes):
        while state['i'] < len(items):
            records.append(rec.line())
            rec = _Record()
            before = state['i']
            if not run(nodes[revert:]) or state['i'] == before:
                break
    records.append(rec.line())
    return records
