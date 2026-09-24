"""
Regression tests for SYPBOD, the supersonic body alone.

Source of truth: datcom-legacy/datcom_2000/sypbod.f.

Checked against a compiled probe (tools/probes/sypbod.py) on every
``SBD`` word, the ``BODY`` curves and the ``BD`` potential/viscous split.
The probe's cases run in one program, so the replay carries the saved
``RACH``.
"""

import json
import math
import pathlib
import sys

import numpy as np
import pytest

from pydatcom.aerodynamics import supersonic_body as module
from pydatcom.aerodynamics.supersonic_body import calculate_sypbod

_ROOT = pathlib.Path(__file__).resolve().parent.parent
sys.path.insert(0, str(_ROOT / 'tools'))
from fortran_data import parse  # noqa: E402

_PROBE = json.loads(
    (_ROOT / 'tests' / 'fixtures' / 'probes' / 'sypbod.json').read_text())
_STALE = {'cda': 0.0, 'thetat': 0.0, 'thetaf': 0.0}


def _replay():
    rach, out = 0.0, []
    for p in _PROBE:
        r = calculate_sypbod(dict(p['inputs'], stale=dict(_STALE,
                                                          rach=rach)))
        rach = r['rach']
        out.append(r)
    return out


_RESULTS = _replay()


@pytest.mark.parametrize("case", range(len(_PROBE)))
def test_matches_compiled_routine(case):
    o, r = _PROBE[case]['outputs'], _RESULTS[case]
    sbd = o['SBD']
    for k, v in r['sbd'].items():
        assert v == pytest.approx(sbd[k - 1], rel=1e-9, abs=1e-14), k
    for k in range(1, 230):
        if k not in r['sbd']:
            assert sbd[k - 1] == 0.0, k
    for key, off in (('cd', 0), ('cl', 20), ('cm', 40), ('cn', 60),
                     ('ca', 80)):
        if key in r:
            np.testing.assert_allclose(r[key], o['BODY'][off:off + len(
                r[key])], rtol=1e-9, atol=1e-14, err_msg=key)
    for key, off in (('req', 0), ('cnpot', 20), ('cnvis', 40),
                     ('cmpot', 60), ('cmvis', 80)):
        if key in r:
            np.testing.assert_allclose(r[key], o['BD'][off:off + len(
                r[key])], rtol=1e-9, atol=1e-14, err_msg=key)


def test_probe_reaches_every_branch():
    cs = [p['inputs'] for p in _PROBE]
    assert any(c['body']['btail'] == 0.0 for c in cs)
    assert any(c['body']['bnose'] == 1.0 for c in cs)
    assert any(c['body']['rla'] == 0.0 and c['body']['btail'] != 0.0
               for c in cs)
    assert any(c['transn'] for c in cs)
    assert any(c['rough'] == 0.0 for c in cs)
    assert {c['body']['ellip'] < 1.0 for c in cs} == {True, False}
    rs = _RESULTS
    assert any(16 in r['sbd'] for r in rs) and any(14 in r['sbd']
                                                   for r in rs)
    assert any(122 in r['sbd'] for r in rs)          # Figure 50/55
    assert any(122 not in r['sbd'] for r in rs)      # Figure 60


def test_afterbody_drag_word_is_not_reset():
    """With no afterbody wave drag CDAB is zeroed but CDA keeps what SBD(119)
    held, and CD0 includes it."""
    c = next(p['inputs'] for p in _PROBE if p['inputs']['body']['btail']
             == 0.0)
    base = calculate_sypbod(dict(c, stale=dict(_STALE, rach=2.0)))
    stale = calculate_sypbod(dict(c, stale=dict(_STALE, rach=2.0,
                                                cda=0.01)))
    assert stale['sbd'][118] == 0.0
    assert stale['sbd'][124] - base['sbd'][124] == pytest.approx(0.01)


def test_tables_match_the_source():
    source = parse('sypbod')
    for name, values in source.items():
        if name.startswith('DUMY') or name == 'I27':
            continue
        np.testing.assert_array_equal(getattr(module, '_' + name), values,
                                      err_msg=name)
    np.testing.assert_array_equal(module._D4220A,
                                  source['DUMY1'] + source['DUMY2'])
    np.testing.assert_array_equal(module._D4350,
                                  source['DUMYA'] + source['DUMYB'])
    np.testing.assert_array_equal(module._D4355,
                                  source['DUMYC'] + source['DUMYD'])
    assert math.isclose(module._X27I[0], 1.5778)
