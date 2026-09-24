"""
Input unit conversion and scaling: CONV.

CONV converts the case's dimensional inputs from the units named on the
``DIM`` card (feet, inches, metres or centimetres) to DATCOM's internal
feet, applies the ``SCALE`` factor to lengths and areas, then shifts every
longitudinal station so the body starts at x = 0.

Arrays are 1-based lists (index 0 unused) named as in the source.

Reference: datcom-legacy/datcom_2000/conv.f
"""

from typing import Dict, Mapping

from pydatcom.utils.constants import UNUSED

# DIM codes 1-4: FT, IN, M, CM, printed as two A1 characters.
_UNIT_NAMES = {1: 'FT', 2: 'IN', 3: 'M ', 4: 'CM'}

# Arrays CONV reads, with their lengths.
_ARRAYS = {'rn': 20, 'alt': 20, 'pinf': 20, 'tinf': 20, 'vinf': 20,
           'grdh': 10, 'syna': 19, 'x': 20, 's': 20, 'p': 20, 'r': 20,
           'zu': 20, 'zl': 20, 'bl': 3, 'wgin': 101, 'htin': 154,
           'vtin': 154, 'tvtin': 8, 'vfin': 154, 'pwin': 29, 'lbin': 21,
           'f': 138}
_SCALARS = ('wt', 'rougfc', 'sref', 'cbarr', 'blref', 'xnx')


def calculate_conv(data: Mapping[str, object]) -> Dict[str, object]:
    """Translate CONV: convert and scale the dimensional inputs.

    Args:
        data: ``idim`` (1 ft, 2 in, 3 m, 4 cm), ``scale``, the flags
            ``symfp``, ``asyfp``, ``trajet``, ``hypef``; the arrays named
            in ``_ARRAYS`` as 0-based sequences of their source lengths;
            the scalars ``wt``, ``rougfc``, ``sref``, ``cbarr``, ``blref``,
            ``xnx``; and ``bd``, the ``/BDATA/`` words 33, 65, 74 and 82.

    Returns:
        The arrays as converted, 1-based with a leading pad; the scalars;
        ``bd`` with word 11 set to ``X0``, the original first body
        station; and ``message``, the line CONV prints.
    """
    v = {n: [0.0] + [float(x) for x in data[n]] for n in _ARRAYS}
    s = {n: float(data[n]) for n in _SCALARS}
    bd = {int(k): float(x) for k, x in data['bd'].items()}
    idim, scale = int(data['idim']), float(data['scale'])
    xl = xa = xr = xp = xt = xw = xf = 1.0
    if idim == 2:
        xl, xa, xp = 12.0, 144.0, 144.0
    if idim >= 3:
        xl, xa, xr, xp = 0.3048, 0.09290304, 0.3048, 0.0208854
        xt, xw, xf = 1.8, 0.2248089, 2.54
    if idim == 4:
        xl, xa, xp = 30.48, 929.0304, 208.854
    ascale = scale ** 2
    length = lambda y: y * scale / xl  # noqa: E731
    area = lambda y: y * ascale / xa  # noqa: E731

    def conv(name, index, fn):
        if v[name][index] != UNUSED:
            v[name][index] = fn(v[name][index])

    def conv_s(name, fn):
        if s[name] != UNUSED:
            s[name] = fn(s[name])

    def conv_bd(word, fn):
        if bd[word] != UNUSED:
            bd[word] = fn(bd[word])

    if not (idim == 1 and scale == 1.0):
        for i in range(1, 21):
            conv('rn', i, lambda y: y * xr)
            conv('alt', i, lambda y: y / xl)
            conv('pinf', i, lambda y: y * xp)
            conv('tinf', i, lambda y: y * xt)
            conv('vinf', i, lambda y: y / xl)
            if i <= 10:
                conv('grdh', i, lambda y: y / xl)
        conv_s('wt', lambda y: y * xw)
        conv_s('rougfc', lambda y: y / xf)
        conv_s('sref', area)
        conv_s('cbarr', length)
        conv_s('blref', length)
        for i in range(1, 20):
            if i not in (4, 8, 10, 13, 18, 19):
                conv('syna', i, length)
        for word in (33, 65, 74, 82):
            conv_bd(word, length)
        for i in range(1, 21):
            for name in ('x', 'p', 'r', 'zu', 'zl'):
                conv(name, i, length)
            conv('s', i, area)
            if i <= 3:
                conv('bl', i, length)
        for i in list(range(1, 7)) + [12]:
            for name in ('wgin', 'htin', 'vtin', 'vfin'):
                conv(name, i, length)
        for i in range(95, 155):
            conv('vtin', i, area)
            conv('vfin', i, area)
            conv('htin', i, area if i >= 115 else length)
        for i in (4, 5, 6, 8, 9, 10, 16, 17, 18, 24, 27, 28):
            conv('pwin', i, length)
        conv('pwin', 19, area)
        conv('pwin', 21, lambda y: y / xl)
        conv('pwin', 22, lambda y: y * xt)
        conv('pwin', 23, lambda y: y * xt)
        conv('pwin', 25, lambda y: y * xp)
        conv('pwin', 26, lambda y: y * xp)
        for i in range(1, 5):
            conv('tvtin', i, length)
        conv('tvtin', 5, area)
        conv('tvtin', 7, length)
        conv('tvtin', 8, length)
        if data['symfp']:
            for i in range(12, 126):
                if (16 <= i <= 38 or 61 <= i <= 84 or 105 <= i <= 114 or
                        i == 117):
                    continue
                conv('f', i, length)
            conv('f', 133, lambda y: y * xl)
            conv('f', 134, lambda y: y * xl)
            conv('f', 135, lambda y: y * xp)
        if data['asyfp']:
            for i in range(12, 16):
                conv('f', i, length)
        if data['trajet']:
            for i in range(12, 22):
                conv('f', i, lambda y: y * xw)
            conv('f', 34, lambda y: y / xl)
            conv('f', 38, lambda y: y / xl)
        if data['hypef']:
            for i in (1, 2, 4):
                conv('f', i, lambda y: y / xl)
        for i in (2, 4, 9, 11, 18, 19):
            conv('lbin', i, lambda y: y / xa)
        for i in (1, 8, 10, 12, 13, 15, 20, 21):
            conv('lbin', i, lambda y: y / xl)

    x0 = v['x'][1]
    bd[11] = x0
    if abs(x0) > UNUSED:
        for i in range(1, int(s['xnx'] + 0.5) + 1):
            v['x'][i] = v['x'][i] - x0
        for i in (1, 2, 6, 9):
            v['syna'][i] = v['syna'][i] - x0
        for i in (11, 12):
            if v['syna'][i] != UNUSED:
                v['syna'][i] = v['syna'][i] - x0
        for i in (4, 16, 18):
            v['pwin'][i] = v['pwin'][i] - x0
        bd[33] = bd[33] - x0
        bd[65] = bd[65] - x0
    message = (f'INPUT DIMENSIONS ARE IN {_UNIT_NAMES[idim]}, '
               f'SCALE FACTOR IS {scale:6.4f}')
    return dict(v, **s, bd=bd, message=message)
