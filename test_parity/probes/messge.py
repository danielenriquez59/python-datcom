"""
Probe MESSGE, the extrapolation diagnostic written to unit 12, and save
the records.

Run from the repository root: ``python test_parity/probes/messge.py``.
MESSGE overlays INTEGER and REAL words through an EQUIVALENCE, so it is
built in the source's own word size (no REAL*8 promotion).  The usual
MESSGE stub is left out.
"""

import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import probe  # noqa: E402

FIXTURE = (Path(__file__).resolve().parent.parent.parent / 'tests' /
           'fixtures' / 'probes' / 'messge.json')
X = [[0.5, 1.0, 1.5, 2.0], [10.0, 20.0, 30.0], [-4.0, -2.0], [0.25, 0.75]]
CASES = [
    {'nstp': 1, 'nf': 0, 'iovly': 12, 'mess': ['EXTR', 'APOL', 'ATED']},
    {'nstp': 2, 'nf': 3, 'iovly': 7, 'mess': ['ONE ']},
    {'nstp': 4, 'nf': 1, 'iovly': 45, 'mess': ['A   ', 'B   ', 'C   ',
                                               'D   ', 'E   ']},
    {'nstp': 3, 'nf': -1, 'iovly': 2, 'mess': ['NONE']},
]


def msscl(c):
    m = [0] * 22
    m[3], m[4], m[5] = len(c['mess']), 5, c['nstp']
    for step in range(1, 5):
        m[4 * step + 2] = 10 * step + 1
        m[4 * step + 3] = len(X[step - 1])
        m[4 * step + 4] = step
        m[4 * step + 5] = step + 1
    return m


def driver():
    lines = ['      PROGRAM PROBE', '      COMMON /OVERLY/ NN(5),NF,NM(2),'
             'IOVLY', '      DIMENSION ROUT(2),MESS(5),MSSCL(21)',
             '      DIMENSION X1(4),X2(3),X3(2),X4(2)',
             '      INTEGER ROUT,MESS',
             "      OPEN(12,FILE='unit12.txt',STATUS='UNKNOWN')",
             "      ROUT(1)=TRANSFER('TLIN',ROUT(1))",
             "      ROUT(2)=TRANSFER('EX  ',ROUT(2))",
             "      MSSCL(1)=TRANSFER('MSG1',MSSCL(1))",
             "      MSSCL(2)=TRANSFER('MSG2',MSSCL(2))"]
    for name, table in zip(('X1', 'X2', 'X3', 'X4'), X):
        for k, v in enumerate(table):
            lines.append(f'      {name}({k + 1})={v!r}')
    for c in CASES:
        m = msscl(c)
        lines.append("      WRITE(12,'(A)') '@@CASE'")
        for k in range(3, 22):
            lines.append(f'      MSSCL({k})={m[k]}')
        for k, w in enumerate(c['mess']):
            lines.append(f"      MESS({k + 1})=TRANSFER('{w}',MESS(1))")
        lines.append(f"      NF={c['nf']}")
        lines.append(f"      IOVLY={c['iovly']}")
        lines.append('      CALL MESSGE(ROUT,MESS,X1,X2,X3,X4,MSSCL)')
    lines += ['      REWIND 12', '   20 READ(12,\'(A)\',END=30) C',
              "      WRITE(6,'(A)') C", '      GO TO 20',
              '   30 STOP', '      END', '']
    lines.insert(5, '      CHARACTER*96 C')
    return '\n'.join(lines)


def main():
    saved = probe.FLAGS
    probe.FLAGS = saved.replace(' -fdefault-real-8', '').replace(
        ' -fdefault-double-8', '')
    try:
        out = probe.run('messge', driver(), ['messge'], stubs='')
    finally:
        probe.FLAGS = saved
    records = []
    for chunk in out.split('@@CASE')[1:]:
        records.append([ln.rstrip() for ln in chunk.split('\n')[1:]
                        if ln.strip()])
    FIXTURE.write_text(json.dumps({'x': X, 'cases': CASES,
                                   'msscl': [msscl(c) for c in CASES],
                                   'records': records}, indent=1))
    print(FIXTURE)


if __name__ == '__main__':
    main()
