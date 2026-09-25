"""
Supersonic control and Mach-shadow overlays: M41O51, M53O65 and M56O70.

Reference: datcom-legacy/datcom_2000/m41o51.f, m53o65.f, m56o70.f
"""

from typing import Callable, Dict, Mapping, Optional

from pydatcom.aerodynamics.body_shadow import calculate_bdarea
from pydatcom.aerodynamics.mach_shadow import calculate_vtarea
from pydatcom.utils.constants import UNUSED


def m41o51(flap_type: float, run_sshing: Callable[[], object]):
    """M41O51: SSHING unless the flap type is 5 or 6, then ``WING(252..
    260)`` set to -UNUSED.  Returns SSHING's result (or ``None``) and the
    ``WING`` words."""
    result = None
    if flap_type != 5.0 and flap_type != 6.0:
        result = run_sshing()
    return result, {250 + word_index: -UNUSED for word_index in range(2, 11)}


def m53o65(control_type: float, deltal, deltar,
           run_spryaw: Callable[[], object]):
    """M53O65: SPRYAW unless the control type is 3, then ``HT(201..210)``
    set to the differential deflections ``F(19..28)-F(29..38)``."""
    result = None
    if control_type != 3.0:
        result = run_spryaw()
    return result, {200 + word_index: deltal[word_index - 1] -
                    deltar[word_index - 1]
                    for word_index in range(1, 11)}


def m56o70(data: Mapping[str, object]) -> Dict[str, object]:
    """M56O70: the Mach-shadow areas for the vertical tail, ventral fin and
    body, where not supplied.

    Args:
        data: By name: ``i`` (the Mach index), ``mach``, ``vtpl``,
            ``vfpl``, ``htpl``, ``vertup``, ``syna`` (``/SYNTSS/`` words
            1-15), ``vtin``, ``vfin`` (the panels' ``VTIN``/``VFIN`` words,
            including 94+I, 114+I, 134+I), ``avt``, ``avf``, ``vt_common``
            (``AVT`` 59, 77, 83, 101), ``wing``, ``tail`` (as
            :func:`mach_shadow.calculate_vtarea`), ``htin`` (3, 4 and
            94+I, 114+I, 134+I), ``aht`` (10, 16, 30, 62), ``body_x``,
            ``body_r``.

    Returns:
        ``vtin``, ``vfin``, ``htin`` (the words as left), ``ig`` (1 when
        BDAREA aborts: the lateral derivatives are deleted and a message
        printed), and ``syna``, restored.

    Notes:
        Kept as executed: each panel's condition reads
        ``PL .AND. A .OR. B .OR. C``, so VTAREA runs for an absent panel
        whenever its second or third area word is unset.  VTAREA's turn of
        the tail to the wing's incidence stays in ``XH``, ``ZH`` for
        BDAREA, and the overlay restores both afterwards.
    """
    mach_index = int(data['i'])
    syna = {int(k): float(v) for k, v in data['syna'].items()}
    saved = (syna[6], syna[7])
    vertup = bool(data['vertup'])
    out = {'ig': 0}
    for key, a_key, up, xv, zv, present in (
            ('vtin', 'avt', vertup, 9, 14, data['vtpl']),
            ('vfin', 'avf', not vertup, 12, 15, data['vfpl'])):
        words = {int(k): float(v) for k, v in data[key].items()}
        flag = (bool(present) and words[mach_index + 94] == UNUSED) or \
            words[mach_index + 114] == UNUSED or \
            words[mach_index + 134] == UNUSED
        if flag:
            r = calculate_vtarea(
                words, data[a_key], up, syna[xv], syna[zv], data['mach'],
                data['wing'], data['tail'],
                {word: syna[word] for word in (2, 3, 4, 6, 7, 8)},
                data['htpl'],
                {int(k): v for k, v in data['vt_common'].items()}, mach_index,
                {134 + mach_index: words[134 + mach_index]})
            words[94 + mach_index] = r['svwb']
            words[134 + mach_index] = r['svhb']
            words[114 + mach_index] = r['svb']
            syna[6], syna[7] = r['syna6'], r['syna7']
        out[key] = words
    htin = {int(k): float(v) for k, v in data['htin'].items()}
    if data['htpl'] and not (htin[94 + mach_index] != UNUSED and
                             htin[114 + mach_index] != UNUSED and
                             htin[134 + mach_index] != UNUSED):
        b = calculate_bdarea(data['body_x'], data['body_r'], data['mach'],
                             {1: syna[1], 6: syna[6], 7: syna[7],
                              8: syna[8]}, htin, data['aht'])
        if b['abort']:
            out['ig'] = 1
        else:
            htin[114 + mach_index], htin[134 + mach_index] = b['sb'], b['s']
            htin[94 + mach_index] = b['xbar']
    out['htin'] = htin
    syna[6], syna[7] = saved
    out['syna'] = syna
    return out
