"""
Regression tests for the INTERX translation.

Source of truth: datcom-legacy/datcom_2000/interx.f.  The fixtures are real
call sites from hbtran.f rather than invented tables, so they exercise the
flat TABLE(LIND,4) unpacking and the first/second variable swap exactly as
the legacy callers do.

Kept deliberately small: INTERX dispatches to TLIN1X and TLINEX, which have
their own translation tests in test_legacy_tables.py.  What needs covering
here is the dispatch, the table unpacking and the argument swap.
"""

import numpy as np
import pytest

from pydatcom.utils.legacy_interp import interx
from pydatcom.interactions.carryover import fig4312_10

# hbtran.f Figure 4.3.1.2-10: CALL INTERX(1,TFIG10,VAR,LGH,DKWB10,KWB,11,11,...)
_TFIG10 = [0.0, .1, .2, .3, .4, .5, .6, .7, .8, .9, 1.0]
_DKWB10 = [1.0, 1.08, 1.16, 1.26, 1.36, 1.46, 1.56, 1.67, 1.78, 1.89, 2.0]

# hbtran.f Figure 4.3.2.2-37A: LGH=(8,3), LIND=8, LDEP=24, LX1U=1.
# T4337A packs both grids end to end; the 999999. is the source sentinel.
_T4337A = [0.0, .4, .8, 1.2, 1.6, 2., 2.4, 2.8, .1, 1.0, 999999.]
_D4337A = [.5, .72, .900, 1.08, 1.24, 1.39, 1.53, 1.68,
           .5, .72, .910, 1.09, 1.25, 1.41, 1.57, 1.72,
           .5, .73, .920, 1.11, 1.27, 1.43, 1.59, 1.74]


@pytest.mark.parametrize("ratio,expected", list(zip(_TFIG10, _DKWB10)))
def test_one_variable_returns_source_table(ratio, expected):
    """NIND=1 dispatches to TLIN1X and reproduces the source grid."""
    assert interx(1, _TFIG10, [ratio], [11], _DKWB10) == pytest.approx(
        expected, abs=1e-12)


def test_one_variable_agrees_with_independent_translation():
    """INTERX and the hand-written Figure 4.3.1.2-10 lookup must agree."""
    for ratio in np.linspace(0.0, 1.0, 51):
        assert interx(1, _TFIG10, [ratio], [11], _DKWB10) == pytest.approx(
            fig4312_10(ratio)['kwb'], abs=1e-12)


@pytest.mark.parametrize("v1,v2,expected", [
    (0.0, 0.1, .5), (2.8, 0.1, 1.68), (0.0, 1.0, .5), (2.8, 1.0, 1.72),
    (1.2, 0.1, 1.08), (1.6, 1.0, 1.25), (2.0, 0.1, 1.39),
])
def test_two_variables_return_source_table(v1, v2, expected):
    """NIND=2 unpacks TABLE(LIND,4), applies the swap, and hits every node."""
    assert interx(2, _T4337A, [v1, v2], [8, 3], _D4337A,
                  lind=8, lx1u=1) == pytest.approx(expected, abs=1e-12)


def test_two_variables_interpolate_between_nodes():
    """A midpoint in the first variable falls between its bracketing nodes."""
    value = interx(2, _T4337A, [0.2, 0.1], [8, 3], _D4337A, lind=8, lx1u=1)
    assert .5 < value < .72


def test_dependent_table_ordering_is_first_variable_fastest():
    """DEP(1)=F(X1,Y1), DEP(2)=F(X2,Y1) as the source header specifies.

    Transposing the dependent table must change the answer, which pins the
    column-major unpacking rather than letting a symmetric fixture hide it.
    """
    table = [0.0, 1.0, 0.0, 10.0]
    dep = [1.0, 2.0, 3.0, 4.0]
    assert interx(2, table, [1.0, 0.0], [2, 2], dep, lind=2) == pytest.approx(2.0)
    assert interx(2, table, [0.0, 10.0], [2, 2], dep, lind=2) == pytest.approx(3.0)


def test_flat_table_requires_lind_for_multiple_variables():
    with pytest.raises(ValueError, match="LIND"):
        interx(2, _T4337A, [1.0, 0.5], [8, 3], _D4337A)


def test_four_variables_rejected_as_in_source():
    """The source comments out the TLIN4X call to save core."""
    with pytest.raises(ValueError, match="TLIN4X"):
        interx(4, _TFIG10, [0.5] * 4, [11] * 4, _DKWB10)


def test_three_variables_reports_the_unresolved_layout():
    """TLIN3X is pending: its Y layout conflicts with TLINEX's."""
    with pytest.raises(NotImplementedError, match="TLIN3X"):
        interx(3, _TFIG10, [0.5] * 3, [11] * 3, _DKWB10, lind=11)


def test_short_dependent_table_is_rejected():
    with pytest.raises(ValueError):
        interx(2, _T4337A, [1.0, 0.5], [8, 3], _D4337A[:10], lind=8)
