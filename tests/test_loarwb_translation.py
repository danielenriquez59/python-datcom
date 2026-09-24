"""
Regression tests for LOARWB and M14O16, the low-aspect-ratio wing-body.

Source of truth: datcom-legacy/datcom_2000/loarwb.f, m14o16.f.

Checked against a compiled probe (tools/probes/loarwb.py), including every
``/SUPDW/`` word the routine sets; the tables are re-parsed from source.
"""

import json
import pathlib
import sys

import numpy as np
import pytest

from pydatcom.aerodynamics import loarwb as module
from pydatcom.aerodynamics.loarwb import calculate_loarwb, calculate_m14o16
from pydatcom.utils.constants import UNUSED

_ROOT = pathlib.Path(__file__).resolve().parent.parent
sys.path.insert(0, str(_ROOT / 'tools'))
from fortran_data import parse  # noqa: E402

_PROBE = json.loads(
    (_ROOT / 'tests' / 'fixtures' / 'probes' / 'loarwb.json').read_text())


def _run(c, routine=calculate_m14o16):
    return routine(c['alpha'], {int(k): v for k, v in c['lbin'].items()},
                   c['mach'], c['reynolds'], c['roughness'], c['stale_xocrb'])


@pytest.mark.parametrize("case", range(len(_PROBE)))
def test_matches_compiled_routine(case):
    c, o = _PROBE[case]['inputs'], _PROBE[case]['outputs']
    r = _run(c)
    for tag, key in [('CD', 'cd'), ('CL', 'cl'), ('CM', 'cm'), ('CN', 'cn'),
                     ('CA', 'ca'), ('CLA', 'cla'), ('CMA', 'cma'),
                     ('KYB', 'kyb'), ('KNB', 'knb'), ('KLB', 'klb'),
                     ('XCP', 'xcp')]:
        np.testing.assert_allclose(r[key], o[tag], rtol=1e-11, atol=1e-16,
                                   err_msg=tag)
    for index, value in r['lb'].items():
        assert value == pytest.approx(o['LB'][index - 1], rel=1e-11,
                                      abs=1e-16), f"LB({index})"
    assert r['xcg'] == o['XCG'][0]


def test_probe_reaches_every_branch():
    inputs = [{int(k): v for k, v in p['inputs']['lbin'].items()}
              for p in _PROBE]
    assert any(i[3] == UNUSED for i in inputs)
    assert any(i[3] != UNUSED for i in inputs)
    assert any(i[1] == 0.0 for i in inputs) and any(i[1] != 0.0
                                                    for i in inputs)
    assert any(i[14] for i in inputs) and any(not i[14] for i in inputs)
    assert any(i[17] for i in inputs) and any(not i[17] for i in inputs)
    capped = [_run(p['inputs'], calculate_loarwb)['lb'][84] <
              p['inputs']['reynolds'] * float(p['inputs']['lbin']['8'])
              for p in _PROBE]
    assert any(capped) and not all(capped)


def test_sharp_nose_keeps_the_stale_blunt_shift():
    """XOCRB is only computed for a round nose; otherwise LB(119) keeps
    whatever it held and enters the centre of pressure."""
    c = next(p['inputs'] for p in _PROBE if not p['inputs']['lbin']['17'])
    r = _run(c, calculate_loarwb)
    assert r['lb'][119] == c['stale_xocrb']
    shifted = _run(dict(c, stale_xocrb=c['stale_xocrb'] + 0.1),
                   calculate_loarwb)
    assert shifted['lb'][116] == pytest.approx(r['lb'][116] + 0.1)


def test_m14o16_recovers_loarwb_normal_and_axial_force():
    for p in _PROBE:
        a = _run(p['inputs'], calculate_loarwb)
        b = _run(p['inputs'])
        np.testing.assert_allclose(b['cn'], a['cn'], rtol=1e-12, atol=1e-18)
        np.testing.assert_allclose(b['ca'], a['ca'], rtol=1e-12, atol=1e-18)


def test_tables_match_the_source():
    source = parse('loarwb')
    for name, values in source.items():
        if name.startswith(('IN', 'N')) or name == 'I27':
            continue
        np.testing.assert_array_equal(getattr(module, '_' + name), values,
                                      err_msg=name)
