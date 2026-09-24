"""
Regression tests for TRNYRL and M40O50, the transonic flap increments.

Source of truth: datcom-legacy/datcom_2000/trnyrl.f, m40o50.f.

Checked against a compiled probe of the overlay (tools/probes/trnyrl.py)
on every ``WING``, ``HT``, ``BODY`` and ``TRN`` word it touches.
"""

import json
import pathlib
import sys

import numpy as np
import pytest

from pydatcom.aerodynamics import transonic_flap as module
from pydatcom.aerodynamics.transonic_flap import (calculate_trnyrl,
                                                  m40o50_words)

_ROOT = pathlib.Path(__file__).resolve().parent.parent
sys.path.insert(0, str(_ROOT / 'tools'))
from fortran_data import parse  # noqa: E402

_PROBE = json.loads(
    (_ROOT / 'tests' / 'fixtures' / 'probes' / 'trnyrl.json').read_text())


def _expected(c):
    r = calculate_trnyrl(c)
    wing = {201 + k: v for k, v in enumerate(c['delcl6'])}
    wing.update({k: 0.0 for k in range(205, 261)})
    for key, start in (('delcl', 201), ('claldl', 241)):
        wing.update({start + k: v for k, v in enumerate(r.get(key, []))})
    wing.update(m40o50_words(c['asyfp']))
    ht = list(c['clrlm6']) + [0.0] * 6
    body = list(c['cnym6']) + [0.0] * 8
    for key in ('clafs', 'clrolt'):
        if key in r:
            ht[:len(r[key])] = r[key]
    if 'cnafs' in r:
        body[:len(r['cnafs'])] = r['cnafs']
    trn = [r.get(k, 0.0) for k in ('encepe', 'yh', 'etaqrs', 'cldelc',
                                   'cldalc', 'kbh', 'khb')]
    return r, [wing[k] for k in range(201, 261)], ht, body, trn


@pytest.mark.parametrize("case", range(len(_PROBE)))
def test_matches_compiled_overlay(case):
    c, o = _PROBE[case]['inputs'], _PROBE[case]['outputs']
    _, wing, ht, body, trn = _expected(c)
    for ours, tag in ((wing, 'WING'), (ht, 'HT'), (body, 'BODY'),
                      (trn, 'TRN')):
        np.testing.assert_allclose(ours, o[tag], rtol=1e-12, atol=0,
                                   err_msg=tag)


def test_probe_reaches_every_branch():
    rs = [calculate_trnyrl(p['inputs']) for p in _PROBE]
    cs = [p['inputs'] for p in _PROBE]
    assert any('delcl' in r and r['transl'] for r in rs)
    assert any('delcl' in r and not r['transl'] for r in rs)
    assert any('cnafs' in r for r in rs)
    assert any('clafs' in r and 'cnafs' not in r for r in rs)
    assert any('cldelc' in r for r in rs) and any('cldalc' in r for r in rs)
    ratios = [(c['sspn'] - c['sspne']) / c['sspn'] for c, r in zip(cs, rs)
              if 'cldelc' in r]
    assert min(ratios) < 0.165 and max(ratios) > 0.3


def test_figure_22_reads_four_of_six_points():
    """The call passes NP=4, so beyond d/b = 0.3 the factor is the parabola
    through the table's second to fourth points."""
    c = next(p['inputs'] for p in _PROBE
             if (p['inputs']['sspn'] - p['inputs']['sspne']) /
             p['inputs']['sspn'] > 0.4 and p['inputs']['mach'] < 1)
    eta = calculate_trnyrl(c)['etaqrs']
    assert eta == pytest.approx(0.554, abs=1e-12)
    from pydatcom.utils.legacy_numeric import tbfunx
    full = tbfunx(module._X12222, module._Y21222, 0.5, 1, 2)[0]
    assert full == pytest.approx(0.464, abs=1e-12)


def test_translating_flap_slope_counts_the_basic_slope_twice():
    c = next(p['inputs'] for p in _PROBE
             if not p['inputs']['asyfp'] and not p['inputs']['htpl']
             and p['inputs']['flap']['type'] == 2.0)
    r = calculate_trnyrl(c)
    claw = c['wing_cla']
    assert r['claldl'][0] == pytest.approx(claw * (2. + c['cfact'][0]))


def test_tables_match_the_source():
    source = parse('trnyrl')
    for name, values in source.items():
        if name.startswith('I') and len(values) == 1:
            continue
        assert getattr(module, '_' + name) == values, name
