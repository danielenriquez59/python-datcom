"""
Regression tests for PTINT1 and VTAREA, the vertical-panel area in the
Mach shadows of the wing and horizontal tail.

Source of truth: datcom-legacy/datcom_2000/ptint1.f, vtarea.f.

Checked against a compiled probe (tools/probes/vtarea.py) on the three
shadowed areas and the tail position the routine leaves.
"""

import json
import pathlib

import pytest

import pydatcom.aerodynamics.mach_shadow as module
from pydatcom.aerodynamics.mach_shadow import calculate_vtarea

_ROOT = pathlib.Path(__file__).resolve().parent.parent
_PROBE = json.loads(
    (_ROOT / 'tests' / 'fixtures' / 'probes' / 'vtarea.json').read_text())


def _run(c, vt_common=None):
    common = vt_common or {int(k): v for k, v in c['vt_common'].items()}
    return calculate_vtarea(c['vtin'], c['avt'], c['vertup'], c['xv'],
                            c['zv'], c['mach'], c['wing'], c['tail'],
                            c['syna'], c['htpl'], common, 2,
                            {136: c['stale134']})


@pytest.mark.parametrize("case", range(len(_PROBE)))
def test_matches_compiled_routine(case):
    c, o = _PROBE[case]['inputs'], _PROBE[case]['outputs']['R']
    r = _run(c)
    ours = [r['svwb'], r['svhb'], r['svb'], r['syna6'], r['syna7']]
    assert ours == pytest.approx(o, rel=1e-13, abs=1e-13)


def test_probe_reaches_every_shadow_shape(monkeypatch):
    seen = set()
    original = module.ptint1

    def spy(*args, **kwargs):
        r = original(*args, **kwargs)
        seen.add((r['nsum'], r['k'] if r['nsum'] == 5 else None, args[6]))
        return r
    monkeypatch.setattr(module, 'ptint1', spy)
    for p in _PROBE:
        _run(p['inputs'])
    for vertup in (True, False):
        shapes = {s for s, _, v in seen if v == vertup}
        assert {0, 2, 4, 5, 6} <= shapes, vertup
        assert {(5, 3, vertup), (5, 4, vertup)} <= seen
    assert {1, 3, 7} <= {s for s, _, _ in seen}


def test_ventral_fin_reads_the_vertical_tails_sweeps():
    """PTINT1 takes the sweeps from /VTDATA/'s AVT, so the ventral fin's
    shadow depends on the vertical tail's planform: with the fin's own
    sweeps in their place the compiled result is not reproduced."""
    k = len(_PROBE) - 1
    c, o = _PROBE[k]['inputs'], _PROBE[k]['outputs']['R']
    assert c['ventral']
    own = {59: 0.2, 77: -0.3, 83: 0.25, 101: -0.2}
    assert _run(c)['svwb'] == pytest.approx(o[0], rel=1e-13)
    assert _run(c, own)['svwb'] != pytest.approx(o[0], rel=1e-6)


def test_tail_incidence_turn_moves_the_tail():
    """With the tail's incidence different from the wing's, VTAREA turns
    the tail about its quarter chord and leaves the new XH, ZH in COMMON."""
    for p in _PROBE:
        c, o = p['inputs'], p['outputs']['R']
        s = c['syna']
        if c['htpl'] and s['4'] != s['8']:
            assert (o[3], o[4]) != (s['6'], s['7'])
            return
    pytest.fail('no turned tail')
