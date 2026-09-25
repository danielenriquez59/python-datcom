"""
Overlays that run one routine and dump its blocks: M11O13 and M42O52.

Reference: datcom-legacy/datcom_2000/m11o13.f, m42o52.f
"""

from typing import Callable, Dict, List, Sequence, Tuple

from pydatcom.io.printers import dmpary
from pydatcom.utils.constants import UNUSED


def m11o13(run_grdeff: Callable[[List[float]], Sequence[float]],
           dump: bool) -> Tuple[List[float], List[str]]:
    """M11O13: GRDEFF on a ``GR`` block (``/SUPWH/`` 1-303) set to
    ``UNUSED``, dumped when ``DPGR`` or ``DMPCSE`` is set.

    ``run_grdeff`` receives the block and returns it as GRDEFF leaves it.
    Returns the block and the dump's records.
    """
    gr = list(run_grdeff([UNUSED] * 303))
    return gr, (dmpary(gr, 'GR', 2) if dump else [])


_M42O52_DUMPS = (('f', '   F'), ('body', 'BODY'), ('ht', '  HT'),
                 ('hyp', ' HYP'), ('vt', '  VT'), ('wing', 'WING'))


def m42o52(run_hypflp: Callable[[], Dict[str, Sequence[float]]],
           dump: bool) -> Tuple[Dict[str, Sequence[float]], List[str]]:
    """M42O52: HYPFLP, then with ``DMPCSE`` or ``DPHYP`` a dump of
    ``F`` (16 words), ``BODY``, ``HT``, ``HYP`` (``/BDATA/`` 1-80), ``VT``
    and ``WING`` (200 words each).

    ``run_hypflp`` returns those blocks by the names ``f``, ``body``,
    ``ht``, ``hyp``, ``vt`` and ``wing``.
    """
    blocks = run_hypflp()
    lines: List[str] = []
    if dump:
        for key, name in _M42O52_DUMPS:
            lines += dmpary(blocks[key], name, 4)
    return blocks, lines
