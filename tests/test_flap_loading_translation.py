"""
Regression tests for AGENR and GDELTA, the flap spanwise loading.

Source of truth: datcom-legacy/datcom_2000/agenr.f, gdelta.f.

Checked against a compiled probe (tools/probes/gdelta.py) on the three
loading curves, the all-moving tail's ``TCD(43..46)`` and the ``BOCH`` it
forms.
"""

import json
import pathlib
import sys

import numpy as np
import pytest

from pydatcom.aerodynamics import flap_loading as module
from pydatcom.aerodynamics.flap_loading import agenr, calculate_gdelta
from pydatcom.utils.legacy_numeric import simul4

_ROOT = pathlib.Path(__file__).resolve().parent.parent
sys.path.insert(0, str(_ROOT / 'tools'))
from fortran_data import parse  # noqa: E402

_PROBE = json.loads(
    (_ROOT / 'tests' / 'fixtures' / 'probes' / 'gdelta.json').read_text())


def _run(c):
    return calculate_gdelta(c['efi'], c['efo'], c['boch'], c['sb'],
                            c['asyfp'], c['tail'])


@pytest.mark.parametrize("case", range(len(_PROBE)))
def test_matches_compiled_routine(case):
    c, o = _PROBE[case]['inputs'], _PROBE[case]['outputs']
    r = _run(c)
    out = o['R']
    if c['asyfp']:
        np.testing.assert_allclose(r['gdh'], out[42:46], rtol=1e-12)
        np.testing.assert_allclose(r['boch'], o['B'], rtol=1e-14)
    else:
        np.testing.assert_allclose(r['gd1'] + r['gd2'] + r['gd3'],
                                   out[:42], rtol=1e-12, atol=1e-15)
        assert out[42:46] == [0.0] * 4
        assert o['B'] == c['boch']


def test_solutions_satisfy_the_influence_equations():
    """Each SIMUL4 solution solves AGENR's system for its right-hand side."""
    a = np.array(agenr([6.0, 6.5, 7.5, 9.0], 30.0)).reshape(4, 4)
    eq = module._EQ
    for i in range(0, 16, 4):
        x = simul4(a.ravel(), eq[i:i + 4])
        np.testing.assert_allclose(a @ x, eq[i:i + 4], atol=1e-12)


def test_full_span_flap_curves_agree():
    """A flap from the root to the tip loads like the full-span curve."""
    r = _run(next(p['inputs'] for p in _PROBE
                  if p['inputs']['efi'] == 0.0 and p['inputs']['efo'] == 1.0))
    np.testing.assert_allclose(r['gd3'], r['gd1'], atol=1e-12)
    assert max(abs(v) for v in r['gd2']) == 0.0


def test_tail_path_skips_the_tip_station():
    """The all-moving tail's chords are formed at SD(4)..SD(1); SD(5), the
    tip, is never used."""
    c = next(p['inputs'] for p in _PROBE if p['inputs']['asyfp'])
    t = c['tail']
    arg1 = (t['tante'] - t['tanle']) * t['bsto2']
    assert _run(c)['boch'][0] == pytest.approx(
        2. * t['bsto2'] / (t['crh'] + 0.924 * arg1), rel=1e-15)


def test_tables_match_the_source():
    for stem in ('agenr', 'gdelta'):
        for name, values in parse(stem).items():
            if name in ('ZRX', 'GD') or (name.startswith('I') and
                                         len(values) == 1):
                continue
            assert getattr(module, '_' + name) == values, name
