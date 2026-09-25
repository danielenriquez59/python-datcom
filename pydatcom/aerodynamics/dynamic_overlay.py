"""
The dynamic-derivative overlay M46O56.

M46O56 runs the dynamic-derivative chain (DYNBOD, DNPAWB, DNPWBT, SUBWBT,
CLRDER), marks every per-angle result word beyond the first angle
``-UNUSED`` (the derivatives are printed at one angle), copies the
wing-body words into the combinations, and dumps the blocks requested.

Blocks are 1-based lists (index 0 unused), mutated in place, as the
callables that stand for the routines also mutate them.

Reference: datcom-legacy/datcom_2000/m46o56.f
"""

from typing import Callable, Dict, List, MutableSequence, Sequence

from pydatcom.io.fortran_format import fortran_write
from pydatcom.io.printers import dmpary
from pydatcom.utils.constants import UNUSED

Block = MutableSequence[float]

# (flag, block, dump name, letters, words) in the source's order.
_IDEAL_DUMPS = (('dpibdy', 'body', 'BODY', 4, 400),
                ('dpiwg', 'wing', 'WING', 4, 400),
                ('dpiht', 'ht', 'HT', 2, 380), ('dpivt', 'vt', 'VT', 2, 380),
                ('dpivf', 'vf', 'VF', 2, 380), ('dpibw', 'bw', 'BW', 2, 380),
                ('dpibh', 'bh', 'BH', 2, 380), ('dpibv', 'bv', 'BV', 2, 380),
                ('dpibwh', 'bwh', 'BWH', 3, 380),
                ('dpibwv', 'bwv', 'BWV', 3, 380),
                ('dpitot', 'bwhv', 'BWHV', 4, 380),
                ('dpipwr', 'power', 'POWR', 4, 200),
                ('dpidwh', 'dwsh', 'DWSH', 4, 60))
# The header is printed for every ideal dump flag except DPIVF.
_HEADER_FLAGS = ('dpibdy', 'dpiwg', 'dpiht', 'dpivt', 'dpibw', 'dpibh',
                 'dpibv', 'dpibwh', 'dpibwv', 'dpitot', 'dpipwr', 'dpidwh')


def m46o56(flags: Dict[str, bool], nalpha: int, mach: float,
           blocks: Dict[str, Block],
           run: Dict[str, Callable[[], object]]) -> Dict[str, object]:
    """M46O56: the dynamic derivatives, their marking and dumps.

    Args:
        flags: The configuration flags ``bo``, ``wgpl``, ``htpl``,
            ``vtpl``, ``vfpl``, ``subson``, ``transn``, ``hypers`` and the
            dump flags ``dmpcse``, ``dpdyn``, ``dpdynh`` and those named in
            ``_IDEAL_DUMPS``; missing flags are false.
        nalpha: The number of angles of attack.  mach: ``FLC(I+2)``.
        blocks: ``body``, ``wing``, ``ht``, ``vt``, ``vf``, ``bw``, ``bh``,
            ``bv``, ``bwh``, ``bwv``, ``bwhv``, ``power``, ``dwsh``, ``dyn``
            and ``dynh``, 1-based.
        run: Callables for ``DYNBOD``, ``DNPAWB``, ``DNPWBT``, ``SUBWBT``
            and ``CLRDER``.

    Returns:
        ``ran`` (the routines called, in order) and ``lines`` (the dump
        records).

    Notes:
        Kept as executed: the ventral fin dump (``DPIVF``) does not bring
        the header line with it.
    """
    f = {k: bool(v) for k, v in flags.items()}
    g = lambda k: f.get(k, False)  # noqa: E731
    b = blocks
    ran: List[str] = []

    def call(name):
        ran.append(name)
        run[name]()

    if g('bo'):
        call('DYNBOD')
        if not (g('hypers') and mach > 1.4):
            if g('wgpl'):
                call('DNPAWB')
            if g('wgpl') and g('htpl'):
                call('DNPWBT')
            if (g('vtpl') or g('vfpl')) and g('subson'):
                call('SUBWBT')
    if not g('hypers'):
        call('CLRDER')
    m = -UNUSED
    for angle_slot in range(2, nalpha + 1):
        for off in range(200, 360, 20):
            b['body'][angle_slot + off] = m
        b['ht'][angle_slot + 200] = b['ht'][angle_slot + 220] = m
        if g('subson'):
            b['ht'][angle_slot + 300] = b['ht'][angle_slot + 320] = b['ht'][angle_slot + 340] = m
        if g('transn'):
            b['ht'][angle_slot + 240] = b['ht'][angle_slot + 260] = m
        for name in ('vt', 'vf', 'bw'):
            for off in (200, 220, 240, 260):
                b[name][angle_slot + off] = m
        b['bwh'][angle_slot + 200] = b['bwh'][angle_slot + 220] = m
        if g('transn'):
            b['bwh'][angle_slot + 240] = b['bwh'][angle_slot + 260] = m
        if not g('hypers'):
            for off in (200, 220, 240, 260, 300, 320, 340):
                b['bh'][angle_slot + off] = m
            for off in (200, 220, 240, 260, 280, 300, 320, 340):
                b['bv'][angle_slot + off] = m
        b['bwhv'][angle_slot + 200] = b['bwhv'][angle_slot + 220] = m
        if g('transn'):
            b['bwhv'][angle_slot + 240] = b['bwhv'][angle_slot + 260] = m
    for angle_slot in range(1, nalpha + 1):
        for off in (300, 320, 340):
            b['bwh'][angle_slot + off] = b['bw'][angle_slot + off]
        if (g('vtpl') or g('vfpl')) and g('subson'):
            b['bwv'][angle_slot + 280] = b['bw'][angle_slot + 280] + b['vt'][angle_slot + 280] + \
                b['vf'][angle_slot + 280]
            b['bwhv'][angle_slot + 280] = b['bwh'][angle_slot + 280] + b['vt'][angle_slot + 280] + \
                b['vf'][angle_slot + 280]
        for off in (200, 220, 240, 260):
            b['bwv'][angle_slot + off] = b['bw'][angle_slot + off]
    lines: List[str] = []
    if g('dmpcse') or g('dpdyn'):
        lines += dmpary(b['dyn'][1:214], 'DYN', 3)
    if g('dmpcse') or g('dpdynh'):
        lines += dmpary(b['dynh'][1:214], 'DYNH', 4)
    if any(g(k) for k in _HEADER_FLAGS):
        lines += fortran_write('(55H0**** THE FOLLOWING ARE IDEAL OUTPUT '
                               'MATRIX ARRAYS ****)')
    for flag, name, label, nlet, words in _IDEAL_DUMPS:
        if g(flag):
            lines += dmpary(b[name][1:words + 1], label, nlet)
    return {'ran': ran, 'lines': lines}
