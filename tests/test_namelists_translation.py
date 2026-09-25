"""
Regression tests for NAMER (the namelist reader) and XNAM1-XNAM23 (the
input namelists).

Source of truth: datcom-legacy/datcom_2000/namer.f and xnam1.f ...
xnam23.f.

NAMER is checked against a compiled probe (tools/probes/namer.py) on six
reads and twelve error stops; the build copy restores the caret that the
shipped namer.f has lost.  Each XNAM routine is checked
(tools/probes/xnam.py) on reading and echoing a namelist that sets every
variable: the printed records and every COMMON word, bit for bit, in the
source's own word size.
"""

import io
import json
import pathlib
import struct

import pytest

from pydatcom.io.namelist_reader import NamelistInputError, namer
from pydatcom.io.namelist_check import nmlist as run_nmlist
from pydatcom.io.namelists import NAMELISTS, exsubt, xnam

_FIXTURES = pathlib.Path(__file__).resolve().parent / 'fixtures' / 'probes'
_NAMER = json.loads((_FIXTURES / 'namer.json').read_text())
_XNAM = json.loads((_FIXTURES / 'xnam.json').read_text())


def _read(cards):
    unit = io.StringIO(''.join(c.ljust(80)[:80] + '\n' for c in cards))
    block = [0.0] + [_NAMER['mark']] * 161
    try:
        r = namer('$', unit, list('FLTCON'), list(''.join(_NAMER['names'])),
                  _NAMER['len'], _NAMER['ldm'], block, _NAMER['loc'])
        return r['lines'], r['ieof'], block
    except NamelistInputError as e:
        return e.lines, None, None


@pytest.mark.parametrize("case", range(len(_NAMER['good'])))
def test_namer_reads_as_compiled(case):
    lines, eof, block = _read(_NAMER['good'][case])
    rec = _NAMER['good_records'][case]
    assert lines == rec['text'] and eof == rec['eof']
    for k in range(1, 162):
        if isinstance(block[k], bool):
            assert k == 161 and block[k] == rec['hypers']
        else:
            assert block[k] == pytest.approx(rec['block'][k - 1], rel=1e-12)


@pytest.mark.parametrize("case", range(len(_NAMER['bad'])))
def test_namer_error_stops_as_compiled(case):
    lines, _, _ = _read(_NAMER['bad'][case])
    assert lines == _NAMER['bad_records'][case]['text']


def _bits(v):
    if isinstance(v, bool):
        return f'{int(v):X}'
    return f'{struct.unpack("<I", struct.pack("<f", v))[0]:X}'


@pytest.mark.parametrize("n", range(1, 24))
def test_xnam_matches_compiled_routine(n):
    rec = _XNAM[str(n)]
    blocks = {}
    for bi, (block, arrays) in enumerate(rec['commons']):
        words = [0.0]
        for ai, (_, size) in enumerate(arrays):
            size = 324 if (block == 'VTI' and size == 316) else size
            words += [k * 0.25 + 1000 * (bi + 1) + ai
                      for k in range(1, size + 1)]
        blocks[block] = words
    text = ''.join(c.ljust(80) + '\n' for c in rec['cards'])
    state = {}
    read = xnam(n, 0, blocks, {NAMELISTS[n]['unit']: io.StringIO(text)},
                state)
    echo = xnam(n, 1, blocks, {}, state)
    assert read['lines'] == rec['read']
    assert echo['lines'] == rec['echo']
    for block, hexes in rec['words'].items():
        assert [_bits(v) for v in blocks[block][1:]] == hexes, block


def test_xnam21_copies_past_its_block():
    """XNAM21 moves VTI words 317-324, beyond the 316 it declares."""
    blocks = {'VTI': [0.0] + [float(k) for k in range(1, 325)],
              'IVF': [0.0] * 364}
    card = io.StringIO(' $VFSCHR $\n')
    xnam(21, 0, blocks, {9: card})
    assert blocks['VTI'][324] == 324.0


def test_xnam_skipped_words_are_not_written_back():
    """XNAM4 reads its tenth planform word but never stores it."""
    blocks = {'WINGI': [0.0] * 102, 'WINGD': [0.0] * 139}
    d = NAMELISTS[4]
    k = d['iequ'].index(10)
    name = d['names'][k]
    card = io.StringIO(f' $WGPLNF {name}=9.$\n')
    xnam(4, 0, blocks, {9: card})
    assert blocks['WINGI'][10] == 0.0


def test_exsubt_rereads_each_name():
    calls = []
    assert exsubt(0, True, 3, lambda: calls.append('rewind'),
                  lambda: calls.append('xnam23')) == 3
    assert calls == ['rewind'] + ['xnam23'] * 3
    assert exsubt(-1, True, 3, None, None) == 0
    assert exsubt(0, False, 3, None, None) == 0


def test_nmlist_routes_to_the_check_table():
    names = {'FLTCON': (['NMACH', 'MACH'], [1, 20]),
             'PLNF': (['CHRDR'], [1])}
    r = run_nmlist(' $WGPLNF CHRDR=2.$'.ljust(80), 9, 5, 1, 0, names, {})
    assert r['name'] == 5 and r['unit6'][0].rstrip().endswith('CHRDR=2.$')
    r = run_nmlist(' $XX NMACH=1.$'.ljust(80), 5, 7, 1, 1, names, {})
    assert r['name'] == 1 and 'UNKNOWN NAMELIST NAME' in r['unit6'][0]
