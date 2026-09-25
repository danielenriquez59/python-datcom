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
    v = {name: [0.0] + [float(x) for x in data[name]] for name in _ARRAYS}
    s = {name: float(data[name]) for name in _SCALARS}
    bd = {int(word): float(value) for word, value in data['bd'].items()}
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

    def conv(name, word, fn):
        if v[name][word] != UNUSED:
            v[name][word] = fn(v[name][word])

    def conv_s(name, fn):
        if s[name] != UNUSED:
            s[name] = fn(s[name])

    def conv_bd(word, fn):
        if bd[word] != UNUSED:
            bd[word] = fn(bd[word])

    if not (idim == 1 and scale == 1.0):
        for word in range(1, 21):
            conv('rn', word, lambda y: y * xr)
            conv('alt', word, lambda y: y / xl)
            conv('pinf', word, lambda y: y * xp)
            conv('tinf', word, lambda y: y * xt)
            conv('vinf', word, lambda y: y / xl)
            if word <= 10:
                conv('grdh', word, lambda y: y / xl)
        conv_s('wt', lambda y: y * xw)
        conv_s('rougfc', lambda y: y / xf)
        conv_s('sref', area)
        conv_s('cbarr', length)
        conv_s('blref', length)
        for word in range(1, 20):
            if word not in (4, 8, 10, 13, 18, 19):
                conv('syna', word, length)
        for word in (33, 65, 74, 82):
            conv_bd(word, length)
        for station_index in range(1, 21):
            for name in ('x', 'p', 'r', 'zu', 'zl'):
                conv(name, station_index, length)
            conv('s', station_index, area)
            if station_index <= 3:
                conv('bl', station_index, length)
        for word in list(range(1, 7)) + [12]:
            for name in ('wgin', 'htin', 'vtin', 'vfin'):
                conv(name, word, length)
        for word in range(95, 155):
            conv('vtin', word, area)
            conv('vfin', word, area)
            conv('htin', word, area if word >= 115 else length)
        for word in (4, 5, 6, 8, 9, 10, 16, 17, 18, 24, 27, 28):
            conv('pwin', word, length)
        conv('pwin', 19, area)
        conv('pwin', 21, lambda y: y / xl)
        conv('pwin', 22, lambda y: y * xt)
        conv('pwin', 23, lambda y: y * xt)
        conv('pwin', 25, lambda y: y * xp)
        conv('pwin', 26, lambda y: y * xp)
        for word in range(1, 5):
            conv('tvtin', word, length)
        conv('tvtin', 5, area)
        conv('tvtin', 7, length)
        conv('tvtin', 8, length)
        if data['symfp']:
            for word in range(12, 126):
                if (16 <= word <= 38 or 61 <= word <= 84 or 105 <= word <= 114
                        or word == 117):
                    continue
                conv('f', word, length)
            conv('f', 133, lambda y: y * xl)
            conv('f', 134, lambda y: y * xl)
            conv('f', 135, lambda y: y * xp)
        if data['asyfp']:
            for word in range(12, 16):
                conv('f', word, length)
        if data['trajet']:
            for word in range(12, 22):
                conv('f', word, lambda y: y * xw)
            conv('f', 34, lambda y: y / xl)
            conv('f', 38, lambda y: y / xl)
        if data['hypef']:
            for word in (1, 2, 4):
                conv('f', word, lambda y: y / xl)
        for word in (2, 4, 9, 11, 18, 19):
            conv('lbin', word, lambda y: y / xa)
        for word in (1, 8, 10, 12, 13, 15, 20, 21):
            conv('lbin', word, lambda y: y / xl)

    x0 = v['x'][1]
    bd[11] = x0
    if abs(x0) > UNUSED:
        for station_index in range(1, int(s['xnx'] + 0.5) + 1):
            v['x'][station_index] = v['x'][station_index] - x0
        for word in (1, 2, 6, 9):
            v['syna'][word] = v['syna'][word] - x0
        for word in (11, 12):
            if v['syna'][word] != UNUSED:
                v['syna'][word] = v['syna'][word] - x0
        for word in (4, 16, 18):
            v['pwin'][word] = v['pwin'][word] - x0
        bd[33] = bd[33] - x0
        bd[65] = bd[65] - x0
    message = (f'INPUT DIMENSIONS ARE IN {_UNIT_NAMES[idim]}, '
               f'SCALE FACTOR IS {scale:6.4f}')
    return dict(v, **s, bd=bd, message=message)
