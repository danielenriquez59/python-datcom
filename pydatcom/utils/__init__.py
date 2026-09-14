"""Utility functions for PyDATCOM - math, atmosphere, interpolation, etc."""

from pydatcom.utils.constants import PI, DEG, RAD, UNUSED, get_constants_dict
from pydatcom.utils.math_utils import (
    arcsin, arccos, area1, area2, det4, sign, linear_interp
)
from pydatcom.utils.atmosphere import Atmosphere
from pydatcom.utils.interpolation import (
    asmint, bilinear_interp, TableInterpolator
)
from pydatcom.utils.table_lookup import (
    fig26, fig53a, fig60b, fig68, angdet, get_table_manager
)
from pydatcom.utils.legacy_numeric import quad, trapz, tbfunx, tranf
from pydatcom.utils.legacy_tables import glook, switch, tlin1x, tlinex
from pydatcom.utils.legacy_interp import interx, eqspc1, eqspce
from pydatcom.utils.tranac import tranac

__all__ = [
    'PI', 'DEG', 'RAD', 'UNUSED', 'get_constants_dict',
    'arcsin', 'arccos', 'area1', 'area2', 'det4', 'sign', 'linear_interp',
    'Atmosphere',
    'asmint', 'bilinear_interp', 'TableInterpolator',
    'fig26', 'fig53a', 'fig60b', 'fig68', 'angdet', 'get_table_manager',
    'quad', 'trapz', 'tbfunx', 'tranf', 'glook', 'switch', 'tlin1x', 'tlinex',
    'interx', 'eqspc1', 'eqspce', 'tranac',
]

