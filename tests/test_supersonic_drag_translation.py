"""
Regression tests for SUPDRG and M18O22, the supersonic wing drag.

Source of truth: datcom-legacy/datcom_2000/supdrg.f, m18o22.f.

Checked against a compiled probe (tools/probes/supdrg.py) on every
``/SUPWH/`` word the routine sets.  The probe's cases run in one program,
so the replay carries SUPDRG's saved local ``RACH`` from case to case.
"""

import json
import pathlib
import sys

import numpy as np
import pytest

from pydatcom.aerodynamics import supersonic as sup
from pydatcom.aerodynamics.supersonic_drag import (calculate_supdrg,
                                                   m18o22_options)
from pydatcom.utils.constants import UNUSED

_ROOT = pathlib.Path(__file__).resolve().parent.parent
sys.path.insert(0, str(_ROOT / 'tools'))
from fortran_data import parse  # noqa: E402

_PROBE = json.loads(
    (_ROOT / 'tests' / 'fixtures' / 'probes' / 'supdrg.json').read_text())

_WORDS = {'beta': 1, 'bovert': 2, 'cdw': 79, 'cdo': 80, 'dragc': 81, 'p': 82,
          'cfo': 83, 'cfi': 84, 'rno': 85, 'rni': 86, 'cdf': 87, 'cf': 88,
          'rlcoff': 89, 'rnn': 90}


def _replay():
    rach, out = 0.0, []
    for p in _PROBE:
        c = p['inputs']
        stale = dict(c['stale'], rach=rach)
        r = calculate_supdrg(c['mach'], c['win'], c['a'], c['sref'],
                             c['roughness'], c['sbw'], stale)
        if c['roughness'] != 0.0:
            rach = min(c['mach'], 3.0)
        out.append(r)
    return out


_RESULTS = _replay()


@pytest.mark.parametrize("case", range(len(_PROBE)))
def test_matches_compiled_routine(case):
    o, r = _PROBE[case]['outputs'], _RESULTS[case]
    slg = o['SLG']
    set_words = {i for i in _WORDS.values()}
    for key, index in _WORDS.items():
        if key in r:
            assert r[key] == pytest.approx(slg[index - 1], rel=1e-12,
                                           abs=1e-300), key
    # Words the routine does not set stay as the driver left them.
    for index in range(1, 142):
        if index not in set_words and index not in (84, 89, 119):
            assert slg[index - 1] == 0.0, index
    assert (r['a62'], r['a86']) == pytest.approx(tuple(o['A']), abs=0.0)


def test_probe_reaches_every_branch():
    cases = [p['inputs'] for p in _PROBE]
    kinds = [float(c['win']['15']) for c in cases]
    assert 1.0 in kinds and 2.0 in kinds and 3.0 in kinds
    assert any(c['roughness'] == 0.0 for c in cases)
    assert any(r['bovert'] < 1. for r in _RESULTS)
    assert any(r['bovert'] >= 1. for r in _RESULTS)
    assert any(r['rnn'] == r['rlcoff'] for r in _RESULTS)
    assert any('rno' in r and 'rni' not in r for r in _RESULTS)  # equal MACs


def test_cranked_wing_leaves_the_snapped_sweep():
    """A zero outboard sweep is set to 1e-5 and restored only on the
    straight-wing path, so a cranked wing leaves 1e-5 in A(86)."""
    hits = [(p, r) for p, r in zip(_PROBE, _RESULTS)
            if float(p['inputs']['a']['86']) == 0.0]
    assert hits
    for p, r in hits:
        assert p['outputs']['A'][1] == pytest.approx(1e-5) == r['a86']


def test_no_roughness_reads_the_previous_mach():
    """Without roughness FIG26 is read at the Mach number the previous call
    left: replaying with the case's own Mach changes the friction."""
    k = next(i for i, p in enumerate(_PROBE)
             if p['inputs']['roughness'] == 0.0)
    c = _PROBE[k]['inputs']
    own = calculate_supdrg(c['mach'], c['win'], c['a'], c['sref'], 0.0,
                           c['sbw'], dict(c['stale'], rach=c['mach']))
    assert own['cf'] != _RESULTS[k]['cf']
    assert _PROBE[k]['outputs']['SLG'][87] == pytest.approx(
        _RESULTS[k]['cf'], rel=1e-12)


def test_m18o22_options():
    o = m18o22_options(UNUSED, UNUSED, 0.0, 250.0, 7.5)
    assert o == {'sref': 250.0, 'cbarr': 7.5, 'roughness': 1.6e-4}
    o = m18o22_options(300.0, 8.0, 2e-4, 250.0, 7.5)
    assert o == {'sref': 300.0, 'cbarr': 8.0, 'roughness': 2e-4}
    # UNUSED (1e-30) is below 1e-10, so it too takes the default.
    assert m18o22_options(1., 1., UNUSED, 0., 0.)['roughness'] == 1.6e-4


def test_tables_match_the_source():
    source = parse('supdrg')
    np.testing.assert_array_equal(sup._FIG_415258_X, source['T15258'])
    np.testing.assert_array_equal(sup._FIG_415258_SHARP, source['DS5258'])
    np.testing.assert_array_equal(sup._FIG_415258_ROUND, source['DR5258'])
    np.testing.assert_array_equal(sup._FIG_415127_MACH, source['X27M'])
    np.testing.assert_array_equal(sup._FIG_415127_CEPT, source['X27I'])
