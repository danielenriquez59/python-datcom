"""
Regression tests for TRANJT, the hypersonic transverse jet sizing.

Source of truth: datcom-legacy/datcom_2000/tranjt.f.

Checked against a compiled probe (tools/probes/tranjt.py) on the 150
``JET`` words and the ``JETA`` words set.
"""

import json
import pathlib
import sys

import numpy as np
import pytest

from pydatcom.aerodynamics import hypersonic_jet as module
from pydatcom.aerodynamics.hypersonic_jet import calculate_tranjt

_ROOT = pathlib.Path(__file__).resolve().parent.parent
sys.path.insert(0, str(_ROOT / 'tools'))
from fortran_data import parse  # noqa: E402

_PROBE = json.loads(
    (_ROOT / 'tests' / 'fixtures' / 'probes' / 'tranjt.json').read_text())


@pytest.mark.parametrize("case", range(len(_PROBE)))
def test_matches_compiled_routine(case):
    c, o = _PROBE[case]['inputs'], _PROBE[case]['outputs']
    r = calculate_tranjt(c)
    assert r['jet'] == pytest.approx(o['JET'][:147], rel=1e-9, abs=1e-14)
    assert [r['qinf'], r['cf0'], r['veoa'], r['fjmax'], r['pjmax'],
            r['dt']] == pytest.approx(o['JETA'][:6], rel=1e-10)
    nt = int(c['nt'])
    assert r['k'] == pytest.approx(o['JETA'][16:16 + nt], rel=1e-10)


def test_centres_of_pressure_overwrite_the_last_weights():
    """XCP starts at JET(138), inside WEIGHT (JET(131..140)): with nine
    schedule points WEIGHT(8), (9) are replaced by XCP(1), (2)."""
    c, o = next((p['inputs'], p['outputs']) for p in _PROBE
                if p['inputs']['nt'] > 7)
    r = calculate_tranjt(c)
    assert o['JET'][137] == pytest.approx(r['xcp'][0], rel=1e-10)
    assert o['JET'][137] != pytest.approx(r['weight'][7], rel=1e-6)


def test_probe_reaches_every_branch():
    rs = [(p['inputs'], calculate_tranjt(p['inputs'])) for p in _PROBE]
    assert any(any(c['laminar']) for c, _ in rs)
    turbulent_m1 = [m for c, r in rs for m, lam in zip(r['m1'],
                                                       c['laminar'])
                    if not lam]
    assert min(turbulent_m1) <= 5.0 < max(turbulent_m1)


def test_tables_match_the_source():
    for name, values in parse('tranjt').items():
        if name.startswith('Q'):
            continue
        np.testing.assert_array_equal(getattr(module, '_' + name), values,
                                      err_msg=name)
