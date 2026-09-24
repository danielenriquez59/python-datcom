"""
Experimental-data namelist sizes: XPERNM, its overlay M34O42, and the
TEST name match.

XPERNM rereads the case's input deck (tape 8) and records, in the
thousands of each ``NLIST`` word, how many cards belong to each ``EXPR``
namelist: the header card and every card after it up to the next ``EXPR``
header, whatever those cards are.

Reference: datcom-legacy/datcom_2000/xpernm.f, m34o42.f, tbtrn.f (TEST)
"""

from typing import Dict, List, Sequence

_KEXP = 'EXPR'


def test_name(kol: str, key: str) -> bool:
    """Translate TEST: ``key`` matches the start of ``kol`` and the next
    column is blank."""
    n = len(key)
    kol = kol.ljust(n + 1)
    return kol[:n] == key and kol[n] == ' '


def xpernm(cards: Sequence[str], nlist: Sequence[int],
           kand: str = '$') -> Dict[str, object]:
    """Translate XPERNM: count the cards of each experimental namelist.

    Args:
        cards: The input deck's lines (tape 8).
        nlist: ``NLIST(1..100)`` as the namelist reader left it.
        kand: The namelist delimiter, ``KAND`` (``$``).

    Returns:
        ``nlist`` (each word below 1000 gains 1000 times its card count),
        ``klist`` (the number of ``EXPR`` namelists, capped at 100),
        ``warning`` (more than 100 were found).  With no ``EXPR`` header
        the source adds the count to ``NLIST(0)``, which is ``KLIST`` in
        ``/EXPER/``, and then overwrites ``KLIST`` with 0.

    Notes:
        Kept as executed: the cards before the first header are dropped,
        and a namelist's count runs to the next ``EXPR`` header, so it
        includes any other namelists between.
    """
    words: List[int] = [0] + [int(v) for v in nlist]
    k = ier = ncards = 0
    first = True
    for card in cards:
        card = card[:80].ljust(80)
        if card[1] == kand and test_name(card[2:], _KEXP):
            if not (first or ier > 0) and words[k] < 1000:
                words[k] += 1000 * ncards
            k += 1
            ncards = 0
            first = False
            if k > 100:
                ier = 1
        ncards += 1
    if ier <= 0 and words[k] < 1000:
        words[k] += 1000 * ncards
    return {'nlist': words[1:101], 'klist': 100 if ier else k,
            'warning': bool(ier), 'method': 'legacy_xpernm'}


def m34o42(klist: int, cards: Sequence[str], nlist: Sequence[int],
           kand: str = '$'):
    """M34O42: XPERNM runs only when the reader found an ``EXPR`` list."""
    if klist >= 1:
        return xpernm(cards, nlist, kand)
    return None
