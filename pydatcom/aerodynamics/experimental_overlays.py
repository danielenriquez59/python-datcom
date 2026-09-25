"""
Overlays around the experimental-data substitution: M19O23, M21O25,
M28O34, and EXPDAT, which stages a Mach number's experimental namelists.

Blocks are 1-based lists (index 0 unused), edited in place; the routines
the overlays call are passed as callables.

Reference: datcom-legacy/datcom_2000/m19o23.f, m21o25.f, m28o34.f,
expdat.f
"""

import math
from typing import Callable, Dict, IO, Iterable, List, MutableSequence

from pydatcom.io.fortran_format import fortran_write
from pydatcom.utils.constants import RAD, UNUSED
from pydatcom.utils.legacy_numeric import tbfunx

Block = MutableSequence[float]


def m19o23(option: Block, a: Block, wingin: Block, body: Block,
           alpha: List[float], nalpha: int, mach: float, tsmach: float,
           run_sypbod: Callable[[], object],
           run_exsubt: Callable[[], object]) -> None:
    """M19O23: the supersonic body.

    Fills unset reference quantities (``SREF`` from ``A(4)``, ``CBARR``
    from ``A(122)``, ``ROUGFC`` 1.6E-4, ``BLREF`` twice the semispan),
    runs SYPBOD and EXSUBT, then above ``TSMACH`` forms the body's
    lift- and moment-curve slopes by TBFUNX at each angle of attack, the
    normal and axial force, and the sideslip slopes ``-CLA`` and
    ``-CMA*CBARR/BLREF``.

    ``option`` holds ``SREF``, ``CBARR``, ``ROUGFC``, ``BLREF``; ``alpha``
    the angles of attack (degrees).
    """
    if option[1] == UNUSED:
        option[1] = a[4]
    if option[2] == UNUSED:
        option[2] = a[122]
    if option[3] == UNUSED:
        option[3] = 1.6e-4
    if option[4] == UNUSED:
        option[4] = 2.0 * wingin[4]
    run_sypbod()
    run_exsubt()
    if mach < tsmach:
        return
    bfact = option[2] / option[4]
    al = [float(v) for v in alpha[:nalpha]]
    cl = [body[20 + angle_slot] for angle_slot in range(1, nalpha + 1)]
    cm = [body[40 + angle_slot] for angle_slot in range(1, nalpha + 1)]
    for angle_slot in range(1, nalpha + 1):
        body[angle_slot + 100] = tbfunx(al, cl, al[angle_slot - 1])[1]
        body[angle_slot + 120] = tbfunx(al, cm, al[angle_slot - 1])[1]
        ca, sa = (math.cos(al[angle_slot - 1] / RAD),
                  math.sin(al[angle_slot - 1] / RAD))
        body[angle_slot + 60] = (body[angle_slot + 20] * ca
                                 + body[angle_slot] * sa)
        body[angle_slot + 80] = (body[angle_slot] * ca
                                 - body[angle_slot + 20] * sa)
        body[angle_slot + 140] = -body[angle_slot + 100]
        body[angle_slot + 160] = -bfact * body[angle_slot + 120]
        body[angle_slot + 180] = 0.0


def m21o25(kepsln: bool, run_exsubt: Callable[[], object],
           run_sdwash: Callable[[int], object]) -> List[str]:
    """M21O25: SDWASH (with ``KEY`` 1 when experimental downwash is
    given, after a first EXSUBT), then EXSUBT.  Returns the calls."""
    ran = []
    if kepsln:
        ran.append('EXSUBT')
        run_exsubt()
    ran.append('SDWASH')
    run_sdwash(1 if kepsln else 0)
    ran.append('EXSUBT')
    run_exsubt()
    return ran


def m28o34(kdwash: List[bool], nalpha: int, bwh: Block, bwhv: Block,
           run_exsubt: Callable[[], object],
           run_supwbt: Callable[[int, int, int], object]) -> None:
    """M28O34: EXSUBT, SUPWBT with the experimental-downwash keys
    (``KDEODA``, ``KQOQIN``, ``KEPSLN``), then the total-configuration
    words ``BWHV(J+20)``, ``(J+100)`` and ``(J+120)`` from ``BWH``."""
    run_exsubt()
    key = [1 if k else 0 for k in kdwash[:3]]
    run_supwbt(key[2], key[0], key[1])
    for angle_slot in range(1, nalpha + 1):
        for off in (20, 100, 120):
            bwhv[angle_slot + off] = bwh[angle_slot + off]


def expdat(nlist: List[int], mach_index: int, unit8: IO[str],
           unit10: IO[str], blocks: Dict[str, Block], list_flag: bool,
           run_exsubt: Callable[[], object],
           echo_xnam23: Callable[[], List[str]]) -> Dict[str, object]:
    """Translate EXPDAT: stage this Mach number's experimental namelists.

    Args:
        nlist: The saved namelist index: each entry ``1000*cards + mach``.
        mach_index: ``I``, the current Mach number.  unit8, unit10: The
            saved experimental cards and the scratch file EXSUBT rereads.
        blocks: ``body``, ``wing``, ``ht``, ``vt``, ``bw``, ``dwash`` as
            filled by the reads (1-based, from their first data word).
        list_flag: ``LIST``: echo the namelists.
        run_exsubt: EXSUBT.  echo_xnam23: XNAM23's echo; returns its
            records.

    Returns:
        ``nnames``, the flags ``mdata``, ``kbody``, ``kwing``, ``kht``,
        ``kvt``, ``kwb``, ``kdwash`` (three), and ``lines``.
    """
    unit8.seek(0)
    unit10.seek(0)
    unit10.truncate()
    flags = dict(mdata=False, kbody=False, kwing=False, kht=False,
                 kvt=False, kwb=False, kdwash=[False, False, False])
    nnames = 0
    for entry in nlist:
        ncard = entry // 1000
        imach = entry - 1000 * ncard
        for _ in range(ncard):
            card = unit8.readline()
            if imach == mach_index:
                text = card.rstrip('\r\n')[:80]
                unit10.write(text.ljust(80).rstrip() + '\n')
        if imach == mach_index:
            nnames += 1
    lines: List[str] = []
    if nnames:
        flags['mdata'] = True
        run_exsubt()
        for checklist_slot in range(1, 6):
            probe_word = (1 + 20 * (checklist_slot - 1)
                          + (40 if checklist_slot >= 4 else 0))
            for key, name in (('kbody', 'body'), ('kwing', 'wing'),
                              ('kht', 'ht'), ('kvt', 'vt'), ('kwb', 'bw')):
                if blocks[name][probe_word] != UNUSED:
                    flags[key] = True
            if checklist_slot <= 3 and blocks['dwash'][probe_word] != UNUSED:
                flags['kdwash'][checklist_slot - 1] = True
    if flags['mdata'] and list_flag:
        lines += fortran_write('(1H1)')
        lines += echo_xnam23()
        lines += fortran_write('(1H1)')
    return dict(flags, nnames=nnames, lines=lines)


def m48o60(run_expdat: Callable[[], object]) -> object:
    """M48O60: the overlay that runs EXPDAT."""
    return run_expdat()
