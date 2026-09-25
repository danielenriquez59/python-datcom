"""
Regression tests for M46O56, the dynamic-derivative overlay.

Source of truth: datcom-legacy/datcom_2000/m46o56.f.

Checked against the compiled overlay (tools/probes/m46o56.py) with its
five routines stubbed to announce themselves: the call sequence, the dump
records and every word of the fifteen blocks.
"""

import json
import pathlib

import pytest

from pydatcom.aerodynamics.dynamic_overlay import m46o56

_PROBE = json.loads((pathlib.Path(__file__).resolve().parent / 'fixtures' /
                     'probes' / 'm46o56.json').read_text())
_CALLEES = ('DYNBOD', 'DNPAWB', 'DNPWBT', 'SUBWBT', 'CLRDER')


def _run(case):
    blocks = {k: [0.0] + [step * n for n in range(1, size + 1)]
              for k, (_, size, step) in _PROBE['blocks'].items()}
    run = {name: (lambda: None) for name in _CALLEES}
    r = m46o56(case['flags'], case['nalpha'], case['mach'], blocks, run)
    return r, blocks


@pytest.mark.parametrize("case", range(len(_PROBE['cases'])))
def test_matches_compiled_overlay(case):
    r, blocks = _run(_PROBE['cases'][case])
    rec = _PROBE['records'][case]
    assert ['CALL ' + n for n in r['ran']] + r['lines'] == rec['text']
    for key, (name, _, _) in _PROBE['blocks'].items():
        assert blocks[key][1:] == pytest.approx(rec['blocks'][name],
                                                rel=1e-12), name


def test_ventral_fin_dump_has_no_header():
    """DPIVF alone dumps VF without the ideal-arrays header line."""
    r, _ = _run({'flags': {'hypers': True, 'dpivf': True}, 'nalpha': 1,
                 'mach': 6.0})
    assert not any('IDEAL OUTPUT' in line for line in r['lines'])
    assert any(line.startswith('    VF(') for line in r['lines'])
