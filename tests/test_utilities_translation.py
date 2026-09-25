"""
Regression tests for the small numerical utilities: ARCSIN, AREA1, DET4,
SLEQ, QUADIN, MACH2, SIMUL2, TLINVS and INTER3.

Source of truth: the routines of the same names in
datcom-legacy/datcom_2000/.  Checked against a compiled probe
(tools/probes/utilities.py).
"""

import json
import pathlib

import numpy as np
import pytest

from pydatcom.interactions.carryover import intkbw
from pydatcom.utils.legacy_numeric import (
    inter3, mach2, quadin, simul2, simul4, sleq, tlinvs,
)
from pydatcom.utils.math_utils import arcsin, area1, det4

_ROOT = pathlib.Path(__file__).resolve().parent.parent
_PROBE = json.loads((_ROOT / 'tests' / 'fixtures' / 'probes' /
                     'utilities.json').read_text())
_IN = _PROBE['inputs']


def _outputs():
    it = iter(_PROBE['outputs'])
    return {name: [next(it) for _ in range(count)] for name, count in [
        ('arcsin', len(_IN['arcsin'])), ('area1', len(_IN['area1'])),
        ('det4', len(_IN['det4'])), ('sleq', len(_IN['sleq'])),
        ('quadin', len(_IN['quadin'])), ('mach2', len(_IN['mach2'])),
        ('simul2', len(_IN['simul2'])),
        ('tlinvs', len(_IN['tlinvs']['queries'])),
        ('inter3', len(_IN['inter3']['queries'])),
        ('simul4', len(_IN['simul4'])), ('intkbw', len(_IN['intkbw']))]}


_OUT = _outputs()


def test_arcsin():
    for v, (r,) in zip(_IN['arcsin'], _OUT['arcsin']):
        assert arcsin(v) == pytest.approx(r, rel=1e-12, abs=1e-15)


def test_area1():
    for v, (r,) in zip(_IN['area1'], _OUT['area1']):
        assert area1(np.array(v['x'], float), np.array(v['y'], float),
                     v['nsum']) == pytest.approx(r, rel=1e-12, abs=1e-14)


def test_det4():
    for v, (r,) in zip(_IN['det4'], _OUT['det4']):
        assert det4(np.array(v)) == pytest.approx(r, rel=1e-12, abs=1e-12)


def test_det4_is_exactly_zero_for_an_exactly_singular_matrix():
    # SIMUL4 takes the source's divide-by-zero path only on an exact zero.
    # Row 4 = row 1 + 2 row 2: np.linalg.det gives about 3e-14 here.
    rows = np.array([[5., -4., -1., -1.], [4., -3., 0., -3.],
                     [-5., 3., -5., -2.], [13., -10., -1., -7.]])
    assert np.linalg.det(rows) != 0.0
    assert det4(rows) == 0.0


def test_sleq():
    for v, r in zip(_IN['sleq'], _OUT['sleq']):
        x, ok = sleq(v['a'], v['b'])
        if r[0] == -7.0:
            assert not ok and x is None
        else:
            np.testing.assert_allclose(x, r, rtol=1e-11, atol=1e-13)


def test_sleq_gives_up_on_a_regular_system_after_a_zero_pivot():
    """The recovery rotates the partly reduced rows; for this regular
    system that never clears the zero pivot, and SLEQ gives up."""
    v = _IN['sleq'][3]
    assert abs(np.linalg.det(np.array(v['a']))) > 1.0
    assert sleq(v['a'], v['b']) == (None, False)


def test_quadin():
    for v, (r,) in zip(_IN['quadin'], _OUT['quadin']):
        assert quadin(v['y'], v['h']) == pytest.approx(r, rel=1e-12,
                                                       abs=1e-15)


def test_mach2():
    for v, (r, ier) in zip(_IN['mach2'], _OUT['mach2']):
        mach, code = mach2(v)
        assert code == ier
        assert mach == pytest.approx(r, rel=1e-10)


def test_simul2():
    for v, r in zip(_IN['simul2'], _OUT['simul2']):
        assert simul2(v['x'], v['c1'], v['c2']) == pytest.approx(
            tuple(r), rel=1e-12)


def test_tlinvs():
    t = _IN['tlinvs']
    for q, (r,) in zip(t['queries'], _OUT['tlinvs']):
        assert tlinvs(t['x1'], t['x2'], t['y'], q[0], q[1]) == \
            pytest.approx(r, rel=1e-12)


def test_inter3():
    tabs = [(t['x1'], t['x2'], t['y']) for t in _IN['inter3']['tables']]
    for q, (r,) in zip(_IN['inter3']['queries'], _OUT['inter3']):
        assert inter3(q[0], q[1], q[2], tabs) == pytest.approx(r, rel=1e-12)


def test_quadin_is_exact_for_quartics_and_its_remainders():
    x = np.linspace(0.0, 1.2, 13)
    assert quadin(x**4, 0.1) == pytest.approx(1.2**5 / 5, rel=1e-12)
    for n in (6, 7, 8):          # trapezoid, Simpson, 3/8 remainders
        y = np.ones(n)
        assert quadin(y, 0.5) == pytest.approx(0.5 * (n - 1))


def test_simul4():
    for v, r in zip(_IN['simul4'], _OUT['simul4']):
        np.testing.assert_allclose(simul4(v['coff'], v['eq']), r, rtol=1e-11)
        a = np.array(v['coff']).reshape(4, 4)   # COFF(4(i-1)+m)
        np.testing.assert_allclose(np.linalg.solve(a, v['eq']), r,
                                   rtol=1e-9)


def test_intkbw():
    for v, r in zip(_IN['intkbw'], _OUT['intkbw']):
        assert intkbw(*v) == pytest.approx(tuple(r), rel=1e-12)
    assert intkbw(0.9, 45.0, 10.0, 3.0, 5.0) is None
    # Both leading-edge forms are reached.
    edges = [np.sqrt(m**2 - 1) / np.tan(np.radians(le)) > 1.0
             for m, le, *_ in _IN['intkbw']]
    assert any(edges) and not all(edges)


def test_supcm0_is_tracm0_with_another_output_word():
    root = _ROOT / 'datcom-legacy' / 'datcom_2000'

    def body(name):
        lines = [l for l in (root / f'{name}.f').read_text().splitlines()
                 if l[:1] not in 'cC*' and 'SUBROUTINE' not in l]
        return sorted(' '.join(l.split()) for l in lines)
    wing = [l.replace('TRA(74)', 'TRA(73)').replace('TRAH(74)', 'TRAH(73)')
            for l in body('tracm0')]
    assert sorted(wing) == body('supcm0')
