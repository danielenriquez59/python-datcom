"""
Regression tests for HYPBOD and M26O32, the hypersonic body.

Source of truth: datcom-legacy/datcom_2000/hypbod.f, m26o32.f.

Checked against a compiled probe of the overlay (tools/probes/hypbod.py)
on every ``SBD`` and ``BODY`` word it sets.
"""

import json
import pathlib
import sys

import numpy as np
import pytest

from pydatcom.aerodynamics import hypersonic_body as module
from pydatcom.utils.constants import PI
from pydatcom.aerodynamics.hypersonic_body import (calculate_hypbod,
                                                   m26o32_slopes)

_ROOT = pathlib.Path(__file__).resolve().parent.parent
sys.path.insert(0, str(_ROOT / 'tools'))
from fortran_data import parse  # noqa: E402

_PROBE = json.loads(
    (_ROOT / 'tests' / 'fixtures' / 'probes' / 'hypbod.json').read_text())

_CLOSE = dict(rtol=1e-11, atol=1e-15)


def _run(c):
    return calculate_hypbod(c['alpha'], c['mach'], c['body'], c['xcg'],
                            c['sref'], c['cbar'], c['stale'])


@pytest.mark.parametrize("case", range(len(_PROBE)))
def test_matches_compiled_overlay(case):
    c, o = _PROBE[case]['inputs'], _PROBE[case]['outputs']
    r = _run(c)
    sbd, body = o['SBD'], o['BODY']
    na = len(c['alpha'])
    for index, value in r['sbd'].items():
        assert value == pytest.approx(sbd[index - 1], rel=1e-11,
                                      abs=1e-15), index
    for key, start in [('theta', 141), ('lx', 161), ('intgcn', 181),
                       ('intgcm', 201)]:
        np.testing.assert_allclose(r[key], sbd[start - 1:start - 1 +
                                               len(r[key])], **_CLOSE)
    for index in range(1, 230):
        if index not in r['sbd'] and not 141 <= index <= 220:
            assert sbd[index - 1] == 0.0, index
    for key, start in [('cd', 1), ('cl', 21), ('cm', 41), ('cn', 61),
                       ('ca', 81)]:
        np.testing.assert_allclose(r[key], body[start - 1:start - 1 + na],
                                   err_msg=key, **_CLOSE)
    s = m26o32_slopes(c['alpha'], r['cl'], r['cm'], c['cbar'], c['blref'])
    for key, start in [('cla', 101), ('cma', 121), ('b141', 141),
                       ('b161', 161), ('b181', 181)]:
        np.testing.assert_allclose(s[key], body[start - 1:start - 1 + na],
                                   rtol=1e-10, atol=1e-14, err_msg=key)


def test_probe_reaches_every_branch():
    results = [(p['inputs'], _run(p['inputs'])['sbd']) for p in _PROBE]
    assert any(c['body']['bnose'] == 1.0 for c, _ in results)
    assert any(c['body']['bnose'] != 1.0 for c, _ in results)
    assert any(s[131] != 0.0 for _, s in results)            # centre flare
    assert any(c['body']['rla'] != 0.0 and s[131] == 0.0     # cylinder
               for c, s in results)
    assert any(c['body']['rla'] == 0.0 for c, _ in results)
    assert any(16 in s for _, s in results)                   # tail flare
    assert any(14 in s and s[3] != 0.0 for _, s in results)   # boattail
    assert any(s[3] == 0.0 for _, s in results)               # no tail


def test_no_tail_stores_the_stale_tail_angle():
    """With no tail THETAT is never computed, and SBD(14) takes whatever
    SBD(135) held."""
    for p in _PROBE:
        c = p['inputs']
        s = _run(c)['sbd']
        if s[3] == 0.0:
            assert s[14] == s[135] == c['stale']['thetat']
            assert p['outputs']['SBD'][13] == c['stale']['thetat']
            return
    pytest.fail('no tail-less case')


def test_centre_section_moment_lacks_rad():
    """CMAA = CMAAF*AAA*D1/(CBAR*SR) + ...: the tail's CMAT divides the same
    term by RAD as well."""
    c = _PROBE[0]['inputs']
    s = _run(c)['sbd']
    d1 = s[5]
    aaa = PI * d1**2 / 4.
    own = s[132] * aaa * d1 / (c['cbar'] * c['sref'])
    arm = s[133] * (c['xcg'] - c['body']['rln']) / c['cbar']
    assert s[134] == pytest.approx(own + arm, rel=1e-12)
    assert own != 0.0
    assert s[134] != pytest.approx(own / 57.2957795 + arm, rel=1e-3)


def test_tables_match_the_source():
    source = parse('hypbod')
    for name in ('T425', 'D425', 'T422', 'D422'):
        np.testing.assert_array_equal(getattr(module, '_' + name),
                                      source[name], err_msg=name)
