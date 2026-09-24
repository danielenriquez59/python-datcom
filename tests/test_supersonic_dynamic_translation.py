"""
Regression tests for CALCA, the supersonic wing acceleration derivatives.

Source of truth: datcom-legacy/datcom_2000/calca.f.

Checked against a compiled probe (tools/probes/calca.py) on the four
``DYN`` words it sets; the tables are re-parsed from the source.
"""

import json
import pathlib
import sys

import pytest

from pydatcom.aerodynamics import supersonic_dynamic as module
from pydatcom.aerodynamics.supersonic_dynamic import calculate_calca

_ROOT = pathlib.Path(__file__).resolve().parent.parent
sys.path.insert(0, str(_ROOT / 'tools'))
from fortran_data import parse  # noqa: E402

_PROBE = json.loads(
    (_ROOT / 'tests' / 'fixtures' / 'probes' / 'calca.json').read_text())


def _run(c):
    return calculate_calca(c['mach'], c['a'], c['win'])


@pytest.mark.parametrize("case", range(len(_PROBE)))
def test_matches_compiled_routine(case):
    c, o = _PROBE[case]['inputs'], _PROBE[case]['outputs']['R']
    r = _run(c)
    if r is None:
        assert o == [-5.0] * 4
        return
    ours = [r['dyn'][k] for k in (16, 22, 26, 27)]
    assert ours == pytest.approx(o, rel=1e-11)


def test_probe_reaches_every_branch():
    rs = [_run(p['inputs']) for p in _PROBE]
    live = [r for r in rs if r is not None]
    assert any(r is None for r in rs)
    assert any(0 < r['gg'] < 1 for r in live)
    assert any(r['gg'] >= 1 for r in live)


def test_tables_match_the_source():
    source = parse('calca')
    for name in ('XEGIN', 'YEGIN', 'XGM', 'YGM', 'XKGIN', 'YKGIN', 'XEDM',
                 'YEDM'):
        assert getattr(module, '_' + name) == source[name], name
