"""
Regression tests for the VTLIFT / VFLIFT figure layer.

Source of truth: datcom-legacy/datcom_2000/vtlift.f, vflift.f.

Only the figures are translated so far; the branch logic that selects among
them is not. These tests re-parse the FORTRAN DATA statements independently
and compare every value, and check that each figure's grid split matches the
dependent-table size its source INTERX call declares.
"""

import pathlib
import re
import sys

import numpy as np
import pytest

from pydatcom.aerodynamics import vertical_lift_figures as figures
from pydatcom.aerodynamics.vertical_lift_figures import (
    fig4132_56a, fig4132_56g, fig4132_60a, fig4132_60b,
    fig4132_61, fig4132_62, fig4132_63,
)

_ROOT = pathlib.Path(__file__).resolve().parent.parent
_SOURCE = _ROOT / 'datcom-legacy' / 'datcom_2000' / 'vtlift.f'

pytestmark = pytest.mark.skipif(
    not _SOURCE.exists(), reason="legacy vtlift.f not present")


def _parse(names):
    """Parse DATA blocks, including the name-then-slash continuation form."""
    lines = _SOURCE.read_text(errors='replace').splitlines()
    out = {}
    for name in names:
        collected, active, seen_slash = [], False, False
        for line in lines:
            if not active:
                if re.match(r'\s*DATA\s+' + name + r'\s*(/|$)', line):
                    active = True
                    if '/' in line:
                        line = line.split('/', 1)[1]
                        seen_slash = True
                    else:
                        continue
                else:
                    continue
            else:
                line = line[6:] if len(line) > 6 else ''
                if not seen_slash:
                    if '/' not in line:
                        continue
                    line = line.split('/', 1)[1]
                    seen_slash = True
            if '/' in line:
                collected.append(line.split('/')[0])
                break
            collected.append(line)
        raw = re.sub(r'E\s+(\d)', r'E+\1', ''.join(collected))
        values = []
        for token in raw.split(','):
            token = token.strip()
            if not token:
                continue
            if re.match(r'^\d+\*', token):
                count, value = token.split('*', 1)
                values.extend([float(value)] * int(count))
            else:
                values.append(float(token))
        out[name] = values
    return out


# --------------------------------------------------------------------------
# Every table against the source file
# --------------------------------------------------------------------------

_PAIRS = [
    ('_FIG_56A_GRID', 'T13246'), ('_FIG_56G_GRID', 'G13246'),
    ('_FIG_56G', 'DG3246'), ('_FIG_60A_GRID', 'A1350'),
    ('_FIG_60A', 'DA50'), ('_FIG_60B_GRID', 'B1350'),
    ('_FIG_60B', 'DB50'), ('_FIG_61_GRID', 'T13251'),
    ('_FIG_61', 'D13251'), ('_FIG_62_SHARP_GRID', 'S13252'),
    ('_FIG_62_SHARP', 'DSHP52'), ('_FIG_62_ROUND_GRID', 'R13252'),
    ('_FIG_62_ROUND', 'DRND52'), ('_FIG_63_GRID', 'T13253'),
    ('_FIG_63', 'D13253'),
]


@pytest.mark.parametrize("attribute,block", _PAIRS,
                         ids=[name for name, _ in _PAIRS])
def test_table_matches_the_source(attribute, block):
    expected = _parse([block])[block]
    assert getattr(figures, attribute) == pytest.approx(expected, rel=1e-12)


def test_the_large_three_variable_table_matches_the_source():
    """D13246 is assembled from six DUMY blocks totalling 1104 values."""
    parsed = _parse([f'DUMY{index}' for index in range(1, 7)])
    expected = [value for index in range(1, 7)
                for value in parsed[f'DUMY{index}']]
    assert len(expected) == 1104
    assert figures._FIG_56A == pytest.approx(expected, rel=1e-12)


# --------------------------------------------------------------------------
# Shapes self-validate against the source INTERX declarations
# --------------------------------------------------------------------------

def test_each_grid_split_matches_its_dependent_table_size():
    """The product of the grid lengths must equal the declared LDEP."""
    assert np.prod(figures._FIG_56A_SHAPE) == len(figures._FIG_56A) == 1104
    assert np.prod(figures._FIG_60A_SHAPE) == len(figures._FIG_60A) == 72
    assert np.prod(figures._FIG_60B_SHAPE) == len(figures._FIG_60B) == 88
    assert np.prod(figures._FIG_63_SHAPE) == len(figures._FIG_63) == 140


def test_grid_arrays_are_long_enough_for_their_lind_strides():
    """Each packed grid must reach its last column at stride LIND."""
    assert len(figures._FIG_56A_GRID) >= 2 * 23 + figures._FIG_56A_SHAPE[2]
    assert len(figures._FIG_60A_GRID) >= 9 + figures._FIG_60A_SHAPE[1]
    assert len(figures._FIG_60B_GRID) >= 11 + figures._FIG_60B_SHAPE[1]
    assert len(figures._FIG_63_GRID) >= 14 + figures._FIG_63_SHAPE[1]


# --------------------------------------------------------------------------
# Lookups
# --------------------------------------------------------------------------

def test_three_variable_figure_reaches_its_corners():
    """The first and last stored values sit at opposite grid corners."""
    assert fig4132_56a(0.0, 0.25, 0.0) == pytest.approx(figures._FIG_56A[0])
    assert fig4132_56a(30.0, 6.0, 1.0) == pytest.approx(figures._FIG_56A[-1])


def test_one_variable_figures_return_their_endpoints():
    for function, grid, table in (
            (fig4132_56g, figures._FIG_56G_GRID, figures._FIG_56G),
            (fig4132_61, figures._FIG_61_GRID, figures._FIG_61)):
        assert function(grid[0]) == pytest.approx(table[0])
        assert function(grid[len(table) - 1]) == pytest.approx(table[-1])


def test_leading_edge_suction_has_two_distinct_curves():
    """Figure 4.1.3.2-62 keeps separate sharp and round variants."""
    assert figures._FIG_62_SHARP != figures._FIG_62_ROUND
    assert fig4132_62(0.5, sharp=True) != fig4132_62(0.5, sharp=False)


def test_every_lookup_returns_a_finite_number():
    assert np.isfinite(fig4132_56a(5.0, 2.0, 0.3))
    assert np.isfinite(fig4132_56g(0.5))
    assert np.isfinite(fig4132_60a(0.5, 1.0))
    assert np.isfinite(fig4132_60b(0.5, 1.0))
    assert np.isfinite(fig4132_61(0.5))
    assert np.isfinite(fig4132_62(0.5))
    assert np.isfinite(fig4132_63(0.5, 1.0))
