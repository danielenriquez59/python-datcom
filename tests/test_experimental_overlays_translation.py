"""
Regression tests for M19O23, M21O25, M28O34 and EXPDAT.

Source of truth: datcom-legacy/datcom_2000/m19o23.f, m21o25.f, m28o34.f,
expdat.f.

M19O23 is checked against a compiled probe (tools/probes/m19o23.py) with
SYPBOD and EXSUBT stubbed, on the reference quantities and body words
1-200.
"""

import io
import json
import pathlib

import pytest

from pydatcom.aerodynamics.experimental_overlays import (expdat, m19o23,
                                                         m21o25, m28o34)
from pydatcom.utils.constants import UNUSED

_PROBE = json.loads((pathlib.Path(__file__).resolve().parent / 'fixtures' /
                     'probes' / 'm19o23.json').read_text())


@pytest.mark.parametrize("case", range(len(_PROBE)))
def test_m19o23_matches_compiled_overlay(case):
    c, o = _PROBE[case]['inputs'], _PROBE[case]['outputs']
    body = [0.0] + [((k % 37) - 11) * 0.0625 for k in range(1, 401)]
    option = [0.0] + list(c['option'])
    a, wingin = [0.0] * 196, [0.0] * 101
    a[4], a[122], wingin[4] = 310.0, 7.5, 16.0
    m19o23(option, a, wingin, body, c['alpha'], len(c['alpha']), c['mach'],
           c['tsmach'], lambda: None, lambda: None)
    assert option[1:] == pytest.approx(o['OPT'], rel=1e-12)
    assert body[1:201] == pytest.approx(o['BODY'], rel=1e-12, abs=1e-14)


def test_m21o25_and_m28o34_sequences():
    calls = []
    assert m21o25(True, lambda: calls.append('X'),
                  lambda key: calls.append(key)) == ['EXSUBT', 'SDWASH',
                                                     'EXSUBT']
    assert calls == ['X', 1, 'X']
    bwh = [0.0] + [float(k) for k in range(1, 381)]
    bwhv = [0.0] * 381
    seen = []
    m28o34([True, False, True], 2, bwh, bwhv, lambda: None,
           lambda *keys: seen.append(keys))
    assert seen == [(1, 1, 0)]
    assert bwhv[21] == 21.0 and bwhv[122] == 122.0 and bwhv[23] == 0.0


def test_expdat_stages_this_mach_and_sets_flags():
    unit8 = io.StringIO('CARD A1\nCARD A2\nCARD B1\nCARD C1\n')
    unit10 = io.StringIO()
    blocks = {name: [0.0] + [UNUSED] * 200 for name in
              ('body', 'wing', 'ht', 'vt', 'bw', 'dwash')}
    blocks['wing'][21] = 0.5          # CL given
    blocks['dwash'][41] = 0.1         # d(eps)/d(alpha) given
    r = expdat([2001, 1002, 1001], 1, unit8, unit10, blocks, True,
               lambda: None, lambda: ['ECHO'])
    assert unit10.getvalue().split('\n')[:3] == ['CARD A1', 'CARD A2',
                                                 'CARD C1']
    assert r['nnames'] == 2 and r['mdata'] and r['kwing']
    assert not r['kbody'] and r['kdwash'] == [False, False, True]
    assert r['lines'] == ['1', 'ECHO', '1']
