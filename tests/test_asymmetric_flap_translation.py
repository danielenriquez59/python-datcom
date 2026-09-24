"""
Regression tests for LATFLP, the rolling and yawing moments of asymmetric
flaps, spoilers and a differentially deflected horizontal tail.

Source of truth: datcom-legacy/datcom_2000/latflp.f.

Checked against a compiled probe (tools/probes/latflp.py) on every word of
``/FLAPIN/``, ``FLA``, ``HT`` 201-230, ``WING`` 201-400 and ``BODY``
201-400.
"""

import json
import math
import pathlib
import sys

import numpy as np
import pytest

from pydatcom.aerodynamics import asymmetric_flap as module
from pydatcom.aerodynamics.asymmetric_flap import calculate_latflp

_ROOT = pathlib.Path(__file__).resolve().parent.parent
sys.path.insert(0, str(_ROOT / 'tools'))
from fortran_data import parse  # noqa: E402

_PROBE = json.loads(
    (_ROOT / 'tests' / 'fixtures' / 'probes' / 'latflp.json').read_text())


def _inputs(c):
    s = dict(c['surface'])
    if c['transn']:
        s['clasec'] = c['win69'] / .8
    return dict(c, surface=s)


_RESULTS = [calculate_latflp(_inputs(p['inputs'])) for p in _PROBE]


@pytest.mark.parametrize("case", range(len(_PROBE)))
def test_matches_compiled_routine(case):
    o, r = _PROBE[case]['outputs'], _RESULTS[case]
    for tag, ours in (('F', r['f'][1:]), ('FLA', r['fla'][1:]),
                      ('HT', r['ht201']), ('WING', r['clrol']),
                      ('BODY', r['cn'])):
        np.testing.assert_allclose(ours, o[tag], rtol=1e-9, atol=1e-14,
                                   err_msg=tag)


def test_probe_reaches_every_branch():
    c = [p['inputs'] for p in _PROBE]
    assert {p['stype'] for p in c} == {1.0, 2.0, 3.0, 4.0, 5.0}
    spans = {(p['f'][14] - p['f'][13]) / p['surface']['bo2'] for p in c
             if p['stype'] in (1.0, 2.0, 3.0) and
             abs(p['surface']['sweple'] - p['surface']['swepte']) > 4.}
    assert any(b <= .4 for b in spans) and any(b > .8 for b in spans)
    assert any(.4 < b <= .6 for b in spans)
    assert any(.6 < b <= .8 for b in spans)


def test_tip_superposition_restores_the_flap():
    """A flap short of the tip is built from two tip panels; the inner
    chord and the stored loadings come back as they were."""
    k = next(i for i, p in enumerate(_PROBE) if p['inputs']['stype'] == 4.
             and p['inputs']['f'][14] / p['inputs']['surface']['bo2'] < .98)
    assert _RESULTS[k]['f'][12] == _PROBE[k]['inputs']['f'][11]
    assert _RESULTS[k]['fla'][3] == 0.0          # BCLOKO stored as zero


def test_rolling_moment_is_antisymmetric():
    """Swapping the left and right deflections reverses the moment."""
    c = _inputs(_PROBE[1]['inputs'])
    f = list(c['f'])
    f[18:28], f[28:38] = f[28:38], f[18:28]
    swapped = calculate_latflp(dict(c, f=f))
    n = int(c['f'][15])
    for j in range(10, 10 + n):
        assert math.isclose(swapped['ht201'][j], -_RESULTS[1]['ht201'][j],
                            rel_tol=1e-12)


def test_tables_match_the_source():
    for name, values in parse('latflp').items():
        np.testing.assert_array_equal(getattr(module, '_' + name), values,
                                      err_msg=name)


def test_overlay_runs_latflp():
    c = _inputs(_PROBE[0]['inputs'])
    assert module.m52o64(c) == calculate_latflp(c)
