"""
Regression tests for LIFTFP, the flap lift, lift-curve slope and maximum
lift.

Source of truth: datcom-legacy/datcom_2000/liftfp.f.

Checked against a compiled probe (tools/probes/liftfp.py) on every word of
``/FLAPIN/``, ``FLP``, ``/SUPWH/`` 282-287 and ``WING`` 201-250.  The
probe's cases run in one program, so the replay carries the saved chord
factors ``CFACTR``.
"""

import json
import pathlib
import sys

import numpy as np
import pytest

from pydatcom.aerodynamics import flap_lift as module
from pydatcom.aerodynamics.flap_lift import calculate_liftfp

_ROOT = pathlib.Path(__file__).resolve().parent.parent
sys.path.insert(0, str(_ROOT / 'tools'))
from fortran_data import parse  # noqa: E402

_PROBE = json.loads(
    (_ROOT / 'tests' / 'fixtures' / 'probes' / 'liftfp.json').read_text())


def _inputs(c):
    s = dict(c['surface'])
    if c['transn']:
        s['cla'] = c['tra70']
        s['clasec'] = c['win69'] / .8
    return dict(c, surface=s)


def _replay():
    state, out = {'cfactr': [0.0] * 4}, []
    for p in _PROBE:
        r = calculate_liftfp(dict(_inputs(p['inputs']), state=state))
        state = r['state']
        out.append(r)
    return out


_RESULTS = _replay()


@pytest.mark.parametrize("case", range(len(_PROBE)))
def test_matches_compiled_routine(case):
    o, r = _PROBE[case]['outputs'], _RESULTS[case]
    for tag, ours in (('FLP', r['flp'][1:]), ('F', r['f'][1:]),
                      ('FCM', r['fcm282']), ('WING', r['wing'])):
        np.testing.assert_allclose(ours, o[tag], rtol=1e-9, atol=1e-14,
                                   err_msg=tag)


def test_probe_reaches_every_flap_type():
    types = {p['inputs']['ftype'] for p in _PROBE}
    assert types == {float(t) for t in range(1, 9)}
    assert any(p['inputs']['htpl'] for p in _PROBE)
    assert any(p['inputs']['transn'] for p in _PROBE)


def test_section_lift_path_stays_on():
    """Once a deflection gives SDCL, later deflections take the
    experimental path too: the third deflection's DCL follows the SDCL
    formula even though its SDCL is UNUSED."""
    k = next(i for i, p in enumerate(_PROBE)
             if p['inputs']['f'][18] == 1e-30 and p['inputs']['f'][19]
             != 1e-30)
    wing = _RESULTS[k]['wing']
    assert wing[2] == pytest.approx(_PROBE[k]['outputs']['WING'][2])
    assert abs(wing[2]) < 1e-20          # SDCL = UNUSED gives ~0


def test_zero_deflection_is_nudged():
    c = next(p for p in _PROBE if p['inputs']['f'][2] == 0.0)
    assert c['outputs']['F'][2] == pytest.approx(0.01)


def test_tables_match_the_source():
    for name, values in parse('liftfp').items():
        if name == 'ITRANS' or (name.startswith('F') and name[1].isdigit()):
            continue
        np.testing.assert_array_equal(getattr(module, '_' + name), values,
                                      err_msg=name)
