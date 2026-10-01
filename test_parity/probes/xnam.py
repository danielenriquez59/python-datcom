"""
Probe XNAM1-XNAM23, the input namelists, and save the fixture.

Run from the repository root: ``python test_parity/probes/xnam.py``.

Each namelist runs in its own program: every COMMON word it declares is
filled with a marker, a namelist file assigns every variable an integer
value (``.TRUE.`` for LOGICALs), then the routine reads it (IOP 0) and
echoes it (IOP 1), and every block word is printed.  The builds use the
source's own word size (no REAL*8 promotion), so LOGICAL words sit in REAL
words as the source intends; the integer values are exact at that size.
``/VTI/`` is declared with 324 words so that XNAM21's copy of eight words
past its 316 lands somewhere definite.
"""

import json
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import probe  # noqa: E402

LEGACY = Path(__file__).resolve().parent.parent.parent / 'datcom-legacy' / \
    'datcom_2000'
FIXTURE = (Path(__file__).resolve().parent.parent.parent / 'tests' /
           'fixtures' / 'probes' / 'xnam.json')
ROUTINES = ['namer', 'namew', 'readcd', 'skipbl', 'findch', 'extrst',
            'toint', 'todec', 'reptct', 'findvn', 'tolog', 'sublog',
            'subint', 'subrea', 'forint', 'forlog', 'forrea']
PATCH = {'namer': [('DATA CARET / 4H\n   /', 'DATA CARET / 4H^   /')],
         'readcd': [('DO 1000 J=1,2', 'DO 1000 J=1,1')]}
sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent))
from pydatcom.io.namelists import NAMELISTS  # noqa: E402


def commons(n):
    """``[(block, [(array, size), ...]), ...]`` as XNAMn declares them."""
    text = (LEGACY / f'xnam{n}.f').read_text()
    stmts, cur = [], None
    for ln in text.splitlines():
        if not ln or ln[0] in 'cC':
            continue
        ln = ln[:72].ljust(72)
        if ln[5] not in ' 0' and cur is not None:
            cur += ln[6:]
            continue
        if cur:
            stmts.append(cur)
        cur = ln[6:]
    stmts.append(cur)
    out = []
    for s in stmts:
        m = re.match(r'\s*COMMON\s*/\s*(\w+)\s*/(.*)', s)
        if not m or m.group(1) == 'CONSNT':
            continue
        parts = re.findall(r'(\w+)\s*(?:\(\s*(\d+)\s*\))?', m.group(2))
        out.append((m.group(1), [(a, int(k) if k else 1) for a, k in parts]))
    return out


def cards(n):
    d = NAMELISTS[n]
    items = []
    for k, (name, dim) in enumerate(zip(d['names'], d['dims'])):
        value = '.TRUE.' if dim < 0 else f'{k + 2}.'
        if name == 'TYPE' and n in (7, 20):
            value = '7.'            # past 4: reported and reset
        items.append(f'{name}={value},')
    items[-1] = items[-1][:-1] + '$'
    out, line = [], f' ${d["name"]} '
    for item in items:
        if len(line) + len(item) > 78:
            out.append(line)
            line = '  '
        line += item
    out.append(line)
    return out


def driver(n) -> str:
    blocks = commons(n)
    lines = ['      PROGRAM PROBE',
             '      COMMON /CONSNT/ PI,DEG,UNUSED,RAD,KAND',
             '      CHARACTER*80 C', '      INTEGER KD',
             "      KD=TRANSFER('$   ',KD)", '      KAND=KD',
             '      UNUSED=1.E-30']
    decl = []
    for block, arrays in blocks:
        arrays = [(a, 324 if (block == 'VTI' and s == 316) else s)
                  for a, s in arrays]
        decl.append(f'      COMMON /{block}/ ' +
                    ','.join(f'{a}({s})' for a, s in arrays))
    lines[1:1] = decl
    label = 100
    for bi, (block, arrays) in enumerate(blocks):
        for ai, (a, s) in enumerate(arrays):
            s = 324 if (block == 'VTI' and s == 316) else s
            label += 1
            lines.append(f'      DO {label} K=1,{s}')
            lines.append(f'{label:5d} {a}(K)=K*0.25+{1000 * (bi + 1) + ai}.')
    unit = NAMELISTS[n]['unit']
    lines.append(f"      OPEN({unit},FILE='cards.txt',STATUS='UNKNOWN')")
    for c in cards(n):
        c = c.ljust(80)
        lines.append(f"      C='{c[:50]}'")
        lines.append(f"      C(51:80)='{c[50:]}'")
        lines.append(f"      WRITE({unit},'(A)') C")
    lines.append(f'      REWIND {unit}')
    lines.append(f"      WRITE(6,'(A)') '@@READ'")
    lines.append(f'      CALL XNAM{n}(0)')
    lines.append(f"      WRITE(6,'(A)') '@@ECHO'")
    lines.append(f'      CALL XNAM{n}(1)')
    lines.append(f"      WRITE(6,'(A)') '@@BLOCKS'")
    for block, arrays in blocks:
        for a, s in arrays:
            s = 324 if (block == 'VTI' and s == 316) else s
            lines.append(f"      WRITE(6,'(A,{s}Z9)') '{block}',({a}(K),"
                         f'K=1,{s})')
    lines += ['      STOP', '      END', '']
    return '\n'.join(lines)


def main():
    saved = probe.FLAGS
    probe.FLAGS = saved.replace(' -fdefault-real-8', '').replace(
        ' -fdefault-double-8', '')
    out = {}
    try:
        for n in range(1, 24):
            text = probe.run(f'xnam{n}', driver(n), ROUTINES + [f'xnam{n}'],
                             layout_patches=PATCH)
            head, blocks = text.split('@@BLOCKS\n')
            read, echo = head.split('@@READ\n')[1].split('@@ECHO\n')
            words = {}
            for line in blocks.strip('\n').split('\n'):
                name, *hexes = line.split()
                words.setdefault(name, []).extend(hexes)
            out[n] = {'cards': cards(n), 'commons': commons(n),
                      'read': [ln for ln in read.split('\n') if ln],
                      'echo': [ln for ln in echo.split('\n') if ln],
                      'words': words}
    finally:
        probe.FLAGS = saved
    FIXTURE.write_text(json.dumps(out))
    print(FIXTURE)


if __name__ == '__main__':
    main()
