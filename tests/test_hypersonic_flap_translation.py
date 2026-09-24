"""
Regression tests for HYPFLP, the hypersonic flap increments.

Source of truth: datcom-legacy/datcom_2000/hypflp.f.

Checked against a compiled probe (tools/probes/hypflp.py) on the normal-
and axial-force and hinge-moment increments.  The probe's cases run in one
program, so the replay carries the saved ``PHE`` and ``CPI2``.  SIMUL2's
crossing search stops at 0.1 percent, which turns last-digit differences
into parts in 10^7 on a few separated cases.
"""

import json
import math
import pathlib

import pytest

import pydatcom.aerodynamics.hypersonic_flap as module
from pydatcom.aerodynamics.hypersonic_flap import calculate_hypflp

_ROOT = pathlib.Path(__file__).resolve().parent.parent
_PROBE = json.loads(
    (_ROOT / 'tests' / 'fixtures' / 'probes' / 'hypflp.json').read_text())


def _replay():
    stale, out = {'phe': 0.0, 'cpi2': 0.0}, []
    for p in _PROBE:
        r = calculate_hypflp(dict(p['inputs'], stale=stale))
        stale = r['stale']
        out.append(r)
    return out


_RESULTS = _replay()


@pytest.mark.parametrize("case", range(len(_PROBE)))
def test_matches_compiled_routine(case):
    o, r = _PROBE[case]['outputs'], _RESULTS[case]
    for key, tag in (('dcn', 'DCN'), ('dca', 'DCA'), ('dcmcn', 'DCMCN'),
                     ('dcmca', 'DCMCA')):
        assert r[key], key
        for k, v in r[key].items():
            assert v == pytest.approx(o[tag][k], rel=1e-6, abs=1e-12), \
                (key, k)
    # Angles outside 0-20 degrees are skipped.
    alphas = _PROBE[case]['inputs']['alpha']
    for j, a in enumerate(alphas):
        if a < 0 or a > 20:
            assert all(o['DCN'][10 * j + n] == 0.0 for n in range(10))


def test_probe_reaches_every_branch(monkeypatch):
    seen = []
    original = module.simul2

    def spy(*args):
        r = original(*args)
        seen.append(r[0])
        return r
    monkeypatch.setattr(module, 'simul2', spy)
    stale = {'phe': 0.0, 'cpi2': 0.0}
    counts = []
    for p in _PROBE:
        seen.clear()
        r = calculate_hypflp(dict(p['inputs'], stale=stale))
        stale = r['stale']
        counts.append((p['inputs']['laminar'], len(seen),
                       sum(x == -1000. for x in seen)))
    total = [sum(1 for a in p['inputs']['alpha'] if 0 <= a <= 20) *
             int(p['inputs']['f']['16']) for p in _PROBE]
    assert any(lam and n < t for (lam, n, _), t in zip(counts, total))
    assert any(not lam and 0 < n < t for (lam, n, _), t in zip(counts,
                                                                total))
    assert any(le for _, _, le in counts)          # separation at the L.E.


def test_attached_flow_does_not_depend_on_the_stale_plateau():
    """The stale CPIP an attached deflection reads cancels: its increments
    are those of the flap pressure CPI2 alone."""
    k = len(_PROBE) - 1
    c = _PROBE[k]['inputs']
    r = _RESULTS[k]
    assert c['f']['5'] == 1.0
    cf, sr = c['f']['4'], c['sref']
    cosdf = math.cos(math.radians(1.0) * 57.29577951308232 / 57.2957795)
    dcn = r['dcn'][10]              # second angle, attached, after separation
    dca = r['dca'][10]
    ratio = dca / dcn
    assert ratio == pytest.approx(math.sin(1.0 / 57.2957795) / cosdf,
                                  rel=1e-9)
    assert cf > 0 and sr > 0
