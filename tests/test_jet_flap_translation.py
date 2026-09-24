"""
Regression tests for JETFP and M55O67, the jet-flap increments.

Source of truth: datcom-legacy/datcom_2000/jetfp.f, m55o67.f.

Checked against a compiled probe of the overlay (tools/probes/jetfp.py)
on every ``WING`` word it sets and on ``JEANGL``.  The probe's cases run
in one program, so the replay carries JETFP's saved ``ETAT``.
"""

import json
import pathlib
import sys

import numpy as np
import pytest

from pydatcom.aerodynamics import jet_flap as module
from pydatcom.aerodynamics.jet_flap import calculate_jetfp

_ROOT = pathlib.Path(__file__).resolve().parent.parent
sys.path.insert(0, str(_ROOT / 'tools'))
from fortran_data import parse  # noqa: E402

_PROBE = json.loads(
    (_ROOT / 'tests' / 'fixtures' / 'probes' / 'jetfp.json').read_text())


def _replay():
    etat, out = 0.0, []
    for p in _PROBE:
        c = dict(p['inputs'])
        c['stale'] = dict(c['stale'], etat=etat)
        r = calculate_jetfp(c)
        etat = r['etat']
        out.append((c, r))
    return out


_RESULTS = _replay()


@pytest.mark.parametrize("case", range(len(_PROBE)))
def test_matches_compiled_overlay(case):
    c, r = _RESULTS[case]
    ours = r['deccl'] + r['delcm'] + r['dclmax'] + r['clab'] + [r['jeangl']]
    np.testing.assert_allclose(ours, _PROBE[case]['outputs']['R'],
                               rtol=1e-12, atol=1e-15)
    assert r['terminated'] == ('_text' in _PROBE[case]['outputs'])


def test_probe_reaches_every_branch():
    kinds = {int(c['f']['74']) for c, r in _RESULTS if not r['skipped']}
    assert kinds == {1, 2, 3, 4, 5}
    assert any(r['skipped'] for _, r in _RESULTS)
    assert any(r['terminated'] for _, r in _RESULTS)
    assert any(r['dclmax'][0] != 9.0 for _, r in _RESULTS)
    assert any(c['a']['118'] == 1.0 for c, _ in _RESULTS)
    assert any(c['f']['17'] > 5 for c, _ in _RESULTS)


def test_combination_flap_stores_no_lift():
    """Type 4 forms its section increments but leaves WING(201..) as it
    was."""
    c, r = next((c, r) for c, r in _RESULTS if c['f']['74'] == 4.0)
    assert r['deccl'] == c['stale']['deccl']


def test_pure_jet_flap_reads_the_saved_turning_efficiency():
    """ETAT is set only for the EBF; a later pure jet flap's moment uses
    the value it left."""
    k = next(i for i, (c, r) in enumerate(_RESULTS)
             if c['f']['74'] == 1.0 and c['stale']['etat'] != 0.0
             and not r['skipped'])
    c, r = _RESULTS[k]
    fresh = calculate_jetfp(dict(c, stale=dict(c['stale'], etat=0.0)))
    assert fresh['delcm'] != r['delcm']


def test_tables_match_the_source():
    for name, values in parse('jetfp').items():
        np.testing.assert_array_equal(getattr(module, '_' + name), values,
                                      err_msg=name)
