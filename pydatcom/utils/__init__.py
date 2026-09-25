"""Utility functions for PyDATCOM - math, atmosphere, interpolation, etc."""

from pydatcom.utils.constants import PI, DEG, RAD, UNUSED, get_constants_dict
from pydatcom.utils.math_utils import arcsin, area1, area2, det4
from pydatcom.utils.atmosphere import Atmosphere
from pydatcom.utils.interpolation import asmint, bilinear_interp
from pydatcom.utils.table_lookup import (
    fig26, fig53a, fig60b, fig68, angdet, get_table_manager
)
from pydatcom.utils.legacy_numeric import (arccos, quad, sign, trapz, tbfunx,
                                           tranf)
from pydatcom.utils.legacy_tables import (
    glook, switch, tlin1x, tlinex, tlin3x, tlin4x,
)
from pydatcom.utils.legacy_interp import interx, eqspc1, eqspce
from pydatcom.utils.tranac import tranac
from pydatcom.utils.packed_tables import (yup, unpack_table, tlip1x, tlip2x,
                                         tlip3x, intep3)

__all__ = [
    'PI', 'DEG', 'RAD', 'UNUSED', 'get_constants_dict',
    'arcsin', 'arccos', 'area1', 'area2', 'det4', 'sign',
    'Atmosphere',
    'asmint', 'bilinear_interp',
    'fig26', 'fig53a', 'fig60b', 'fig68', 'angdet', 'get_table_manager',
    'quad', 'trapz', 'tbfunx', 'tranf', 'glook', 'switch',
    'tlin1x', 'tlinex', 'tlin3x', 'tlin4x',
    'interx', 'eqspc1', 'eqspce', 'tranac',
    'yup', 'unpack_table', 'tlip1x', 'tlip2x', 'tlip3x', 'intep3',
]

