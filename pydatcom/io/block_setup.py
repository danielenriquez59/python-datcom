"""
Block initialisation and set-up: INIZ, INITZ1, INITZ2, CLEARA, SYNDIM,
SECI and SECO.

Blocks are 1-based lists (index 0 unused), edited in place, named as the
source's COMMON arrays.

Reference: datcom-legacy/datcom_2000/iniz.f, initz1.f, initz2.f,
cleara.f, syndim.f, seci.f, seco.f
"""

import math
from typing import Dict, List, MutableSequence, Sequence

from pydatcom.utils.constants import RAD, UNUSED

Block = MutableSequence[float]
BLANK = '    '

# INITZ1's blocks and their lengths, as it sets them.
INITZ1_BLOCKS = {'c': 51, 'd': 55, 'hb': 39, 'wb': 39, 'fact': 182,
                 'wbt': 156, 'pw': 315, 'cht': 51, 'dht': 55, 'dvf': 55,
                 'dvt': 55, 'stbh': 135, 'stb': 135, 'sbd': 229, 'gr': 303,
                 'shb': 61, 'swb': 61, 'dwa': 237, 'trah': 108, 'tra': 108,
                 'body': 400, 'wing': 400, 'ht': 380, 'vt': 380, 'vf': 380,
                 'bw': 380, 'bh': 380, 'bv': 380, 'bwh': 380, 'bwv': 380,
                 'bwhv': 380, 'power': 200, 'dwash': 60}


def _fill(block: Block, first: int, last: int, value) -> None:
    for word in range(first, last + 1):
        block[word] = value


def iniz(bwh: Block, bwhv: Block, bh: Block) -> None:
    """INIZ: ``BWH`` 1-47, ``BWHV`` 1-7 and ``BH`` 1-144 to ``UNUSED``."""
    _fill(bh, 1, 144, UNUSED)
    _fill(bwh, 1, 47, UNUSED)
    _fill(bwhv, 1, 7, UNUSED)


def initz1() -> Dict[str, List[float]]:
    """INITZ1: every block it clears, set to ``UNUSED``."""
    return {name: [0.0] + [UNUSED] * size
            for name, size in INITZ1_BLOCKS.items()}


def initz2(powr: Block, fcm: Block, supdw: Block, body: Block, wing: Block,
           ht: Block, vt: Block, vf: Block) -> None:
    """INITZ2: ``/POWR/`` 1-315, ``/SUPWH/`` 1-282, ``/SUPDW/`` 1-93 and
    the result words from 201 to ``UNUSED``.

    Kept as executed: ``VF(I+200)`` is set twice and ``VT`` never, so the
    vertical tail's result words keep their values.
    """
    _fill(powr, 1, 315, UNUSED)
    _fill(fcm, 1, 282, UNUSED)
    _fill(supdw, 1, 93, UNUSED)
    _fill(body, 201, 400, UNUSED)
    _fill(wing, 201, 400, UNUSED)
    _fill(ht, 201, 380, UNUSED)
    _fill(vf, 201, 380, UNUSED)


def cleara() -> Dict[str, object]:
    """CLEARA: the figure-tracking state, cleared.

    Returns the cleared words by name: the 34 figure arrays (121 words,
    zero), ``LFIGN``, ``LFIGO``, ``LFIGS`` (zero), the extrapolation
    records, the names set to blank, ``IOVLY=999`` and the counters, and
    ``IFIG`` (20 by 121 blanks).  Kept as executed: of ``LDUM`` only word
    121 is cleared.
    """
    zeros = [0.0] * 121
    out: Dict[str, object] = {f'afig{fig_index:02d}': list(zeros)
                              for fig_index in range(1, 15)}
    out.update({f'jfig{fig_index:02d}': [0] * 121
                for fig_index in range(1, 21)})
    out.update(lfign=[0] * 121, lfigo=[0] * 121, lfigs=[0] * 121,
               ldum121=0, iexcd=[0] * 4, xll=[0.0] * 4, xul=[0.0] * 4,
               iextrl=[[0, 0] for _ in range(4)],
               iextru=[[0, 0] for _ in range(4)], xval=[0.0] * 4,
               rout=[BLANK] * 2, msscl=[BLANK] * 2, ifigst=[BLANK] * 20,
               iovly=999, iovl=0, nstq=0, nstp=0, finalr=0.0, nfig=0,
               ifign=0, routl=BLANK, msscll=BLANK,
               ifig=[[BLANK] * 121 for _ in range(20)])
    return out


def syndim(bd: Block, a: Block, aht: Block, sspn: float, sspne: float,
           syna: Sequence[float], htin: Sequence[float]) -> None:
    """SYNDIM: the wing's body-axis dimensions in ``/BDATA/`` and the
    moment arms ``A(173)`` and ``AHT(173)``.

    ``syna`` and ``htin`` are 1-based; ``sspn``/``sspne`` are the wing
    semispans (``/WINGI/`` words 4 and 3).
    """
    bd[87] = 2. * (sspn - sspne)
    bd[78] = math.sin(bd[77] / RAD)
    bd[79] = math.cos(bd[77] / RAD)
    bd[80] = bd[78] / bd[79]
    bd[66] = 0.5 * bd[87] * a[62] * bd[79]
    bd[67] = bd[33] - (bd[65] + bd[66])
    a[173] = bd[67]
    bd[68] = bd[74] - bd[65] * bd[80]
    xcg, xh, alih = syna[1], syna[6], syna[8]
    aht[173] = xcg - (xh + (htin[4] - htin[3]) * aht[62] *
                      math.cos(alih / RAD))


def seci(a: Sequence[float], typein: Sequence, straight: bool
         ) -> Dict[str, object]:
    """SECI: unpack an airfoil-section input block.

    Args:
        a: The 380-word section block, 1-based.  typein: The planform's
            input words, 1-based (``TYPEIN(15)`` tested by ``straight``).
        straight: Whether ``TYPEIN(15)`` is ``STRA``.

    Returns:
        ``atype`` (negated for a planform that is not straight-tapered),
        ``l`` (the station count), ``naca`` (80 words), ``x``, and either
        ``yu``/``yl`` (type 1) or ``cam``/``thn`` (type 2), each ``l``
        long, ``cla`` (20), and ``cbar``, the planform's reference chord.
    """
    out: Dict[str, object] = {}
    atype = a[1]
    l = 60
    out['naca'] = list(a[2:82])
    if atype > UNUSED:
        l = int(a[82] + 0.5)
        out['x'] = list(a[83:83 + l])
        if atype == 1.:
            out['yu'] = list(a[133:133 + l])
            out['yl'] = list(a[183:183 + l])
        elif atype == 2.:
            out['cam'] = list(a[133:133 + l])
            out['thn'] = [v / 2. for v in a[183:183 + l]]
    if not straight:
        atype = -atype
    out['cla'] = [typein[20 + mach_slot] for mach_slot in range(1, 21)]
    chrdtp, sspnop, sspne, sspn, chrdbp, chrdr = (typein[word]
                                                  for word in range(1, 7))
    if sspnop <= 10. * UNUSED:
        chrdbp = chrdtp
        sspnop = 0.
    tapri = chrdbp / chrdr
    crei = chrdr * (tapri + (1. - tapri) * (sspne - sspnop) /
                    (sspn - sspnop))
    tapri = chrdbp / crei
    tapro = 0.
    if chrdbp != 0.:
        tapro = chrdtp / chrdbp
    creo = chrdbp
    cbari = 2. * crei * (1. + tapri + tapri ** 2) / (3. * (1. + tapri))
    cbaro = 2. * creo * (1. + tapro + tapro ** 2) / (3. * (1. + tapro))
    areai = (sspn - sspnop) * crei * (1. + tapri)
    areao = sspnop * creo * (1. + tapro)
    out['cbar'] = (cbari * areai + cbaro * areao) / (areai + areao)
    out.update(atype=atype, l=l)
    return out


def seco(a: Block, camber: float, atype: float, nmach: int,
         s: Dict[str, object]) -> None:
    """SECO: fill the section block's unset words from the computed
    section.

    ``a`` is the 162-word block, 1-based, edited in place.  ``s`` holds
    the section values: ``tovc``, ``deltay``, ``xovc``, ``cli``, ``ai``,
    ``cla``, ``clmax``, ``xac`` (per Mach), ``cmco4``, ``rho``,
    ``clmax0``, ``cla0``, ``alo``, ``covc``.  ``nmach`` is
    ``FLC(1)+0.5`` truncated.
    """
    def fill(word, value):
        if a[word] == UNUSED:
            a[word] = value

    fill(16, 2. * s['tovc'])
    fill(17, s['deltay'])
    fill(18, s['xovc'])
    fill(19, s['cli'])
    fill(20, s['ai'])
    for mach_slot in range(1, nmach + 1):
        fill(mach_slot + 20, s['cla'][mach_slot - 1])
        fill(mach_slot + 40, s['clmax'][mach_slot - 1])
        fill(mach_slot + 71, s['xac'][mach_slot - 1])
    fill(61, s['cmco4'])
    fill(62, s['rho'])
    a[64] = camber
    fill(68, s['clmax0'])
    fill(69, s['cla0'])
    fill(10, s['alo'])
    fill(93, s['covc'])
    if atype < 0.:
        fill(63, s['rho'])
        fill(65, 2. * s['tovc'])
        fill(66, s['xovc'])
        fill(67, s['cmco4'])


def m51o63(initze: int, run_initz1, run_initz2, powr: Block) -> List[str]:
    """M51O63: INITZ1 (``INITZE`` 1), INITZ2 (2), or ``/POWR/`` 1-315 to
    ``UNUSED`` (3).  Returns the routines called."""
    ran = []
    if initze == 1:
        ran.append('INITZ1')
        run_initz1()
    if initze == 2:
        ran.append('INITZ2')
        run_initz2()
    if initze == 3:
        _fill(powr, 1, 315, UNUSED)
    return ran
