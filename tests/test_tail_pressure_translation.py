"""
Regression tests for DPRESR, the supersonic dynamic pressure at the tail.

Source of truth: datcom-legacy/datcom_2000/dpresr.f.

Checked against a compiled probe (tools/probes/dpresr.py) on the pressure
ratio, the Mach number and the ``DWA`` words set.
"""

import json
import pathlib

import pytest

from pydatcom.aerodynamics.tail_pressure import _nu, calculate_dpresr

_ROOT = pathlib.Path(__file__).resolve().parent.parent
_PROBE = json.loads(
    (_ROOT / 'tests' / 'fixtures' / 'probes' / 'dpresr.json').read_text())


def _run(c):
    return calculate_dpresr(c['zj'], c['zwake'], c['alpha'], c['dwangl'],
                            c['mach'], c['cr'], c['a12'], c['rl2'],
                            c['slope'])


@pytest.mark.parametrize("case", range(len(_PROBE)))
def test_matches_compiled_routine(case):
    c, o = _PROBE[case]['inputs'], _PROBE[case]['outputs']['R']
    r = _run(c)
    qq = c['qq_in'] if r['qqinfy'] is None else r['qqinfy']
    ours = [qq, r['mj'], r['dle'], r['deltaz'], r['xsur'],
            r.get('theta1', 0.0), r.get('delte', 0.0), r.get('thete', 0.0)]
    for name, a, b in zip(('qqinfy', 'mj', 'dle', 'deltaz', 'xsur',
                           'theta1', 'delte', 'thete'), ours, o):
        assert a == pytest.approx(b, rel=1e-12, abs=1e-14), name


def test_probe_reaches_every_branch():
    results = [(p['inputs'], _run(p['inputs'])) for p in _PROBE]
    assert any(r['dle'] >= 0 and c['zj'] > 0 for c, r in results)
    assert any(r['dle'] >= 0 and c['zj'] < 0 for c, r in results)
    assert any(r['dle'] < 0 and c['zj'] > 0 for c, r in results)
    assert any(r['dle'] < 0 and c['zj'] < 0 for c, r in results)
    assert any(r['qqinfy'] is None for _, r in results)          # M(1) < 1
    live = [(c, r) for c, r in results if r['qqinfy'] is not None]
    assert any(c['zj'] / r['z'][5] > 1 for c, r in live)         # TBFUNX
    assert any(c['zj'] / r['z'][5] <= 1 for c, r in live)        # blend
    assert any('thete' not in r for _, r in live)                # TE fan
    assert any('thete' in r and r['delte'] != abs(
        -c['alpha'] + c['slope'][5]) for c, r in live)            # detached
    # Shocks after the leading edge: a slope that rises along the chord.
    assert any(any(s[i] < s[i + 1] for i in range(5)) and r['dle'] < 0
               for s, r in ((c['slope'], r) for c, r in live))


def test_free_stream_expansion_angle_is_truncated():
    """KNUINF is implicitly INTEGER: at Mach 2 the free-stream angle
    26.38 degrees enters the leading-edge expansion as 26."""
    assert 26.3 < _nu(2.0) < 26.4
    c = next(p['inputs'] for p in _PROBE
             if _run(p['inputs'])['dle'] >= 0 and p['inputs']['mach'] == 2.0)
    r = _run(c)
    from pydatcom.utils.legacy_numeric import mach2
    assert r['m'][0] == pytest.approx(mach2(26.0 + r['dle'])[0], rel=1e-14)
    assert r['m'][0] != pytest.approx(mach2(_nu(2.0) + r['dle'])[0],
                                      rel=1e-6)
