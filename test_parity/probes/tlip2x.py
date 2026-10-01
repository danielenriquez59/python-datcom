"""
Probe TLIP2X, TLINEX over a packed table, and save the fixture.

Run from the repository root: ``python test_parity/probes/tlip2x.py``.
The ordinates are packed two to a word in YUP's format by :func:`pack`.
"""

import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from probe import assign, parse_records, run, save  # noqa: E402

ROUTINES = ['tlip2x', 'tlip1x', 'yup', 'switch', 'glook', 'quad']


def _digits(v):
    """Three digits and a decimal exponent (0-9, either sign)."""
    if v == 0:
        return 0, 0
    for e in range(-9, 10):
        d = abs(v) / 10.0 ** e
        if d <= 999 and abs(d - round(d)) < 1e-9:
            return int(round(d)), e
    raise ValueError(v)


def pack(values):
    """YUP's packing: pairs of three-digit numbers in one INTEGER."""
    words = []
    vals = list(values) + ([0.0] if len(values) % 2 else [])
    for a, b in zip(vals[::2], vals[1::2]):
        da, ea = _digits(a)
        db, eb = _digits(b)
        code = 0
        if ea < 0 and eb < 0:
            code = 1
        elif eb < 0:
            code = 2
        elif ea < 0:
            code = 3
        w = db + 1000 * da + 10 ** 6 * abs(eb) + 10 ** 7 * abs(ea) + \
            10 ** 8 * code + 10 ** 9 * (1 if b < 0 else 0)
        words.append(-w if a < 0 else w)
    return words


X1 = [0.0, 1.0, 2.0]
X2 = [0.0, 0.5, 1.0, 2.0]
# Y(NX2, NX1), column-major.
Y = [1.0, 1.5, 2.25, 4.0,
     -3.0, -2.5, 0.125, 1.25,
     0.05, 0.4, 7.5, 12.0]
Y3 = Y + [2.0 * v for v in Y]          # a second X3 slice
CASES = [
    (0.5, 0.25, 0, 0, 0, 0), (1.5, 1.5, 0, 0, 0, 0), (1.0, 0.5, 0, 0, 0, 0),
    (-0.5, 0.3, 1, 0, 1, 0), (2.5, 0.3, 2, 0, 2, 0), (2.5, 2.5, 1, 1, 1, 1),
    (0.7, -0.4, 2, 2, 2, 2), (0.7, 3.0, 0, 2, 0, 2), (3.0, -1.0, 0, 0, 0, 0),
    (2.0, 2.0, 0, 0, 0, 0),
]


def driver():
    py, py3 = pack(Y), pack(Y3)
    lines = ['      PROGRAM PROBE',
             f'      DIMENSION X1(3),X2(4),NY({len(py)}),NY3({len(py3)}),'
             'NXX(7),NXS(7)',
             '      DIMENSION ROUT(2),MESS(2)']
    lines += [assign(f'X1({k + 1})', v) for k, v in enumerate(X1)]
    lines += [assign(f'X2({k + 1})', v) for k, v in enumerate(X2)]
    lines += [f'      NY({k + 1})={w}' for k, w in enumerate(py)]
    lines += [f'      NY3({k + 1})={w}' for k, w in enumerate(py3)]
    lines += ['      NXX(1)=3'] + [f'      NXX({k})=0' for k in range(2, 8)]
    lines += ['      NXS(1)=-4', '      NXS(2)=3', '      NXS(3)=2',
              '      NXS(4)=0', '      NXS(5)=0', '      NXS(6)=2',
              '      NXS(7)=0']
    for q1, q2, l1, l2, u1, u2 in CASES:
        lines.append("      WRITE(6,'(A)') 'CASE'")
        for tag, arr, nxx in (('PLAIN', 'NY', 'NXX'), ('SLICE', 'NY3',
                                                        'NXS')):
            lines.append(assign('A1', q1))
            lines.append(assign('A2', q2))
            lines.append(f'      CALL TLIP2X(X1,X2,{arr},{nxx},4,A1,A2,YA,'
                         f'{l1},{l2},{u1},{u2},')
            lines.append('     1            MESS,1,ROUT)')
            lines.append(f"      WRITE(6,'(A,ES25.16)') '{tag}',YA")
    lines += ['      STOP', '      END', '']
    return '\n'.join(lines)


def main():
    records = parse_records(run('tlip2x', driver(), ROUTINES))
    print(save('tlip2x', {'x1': X1, 'x2': X2, 'y': Y, 'y3': Y3,
                          'packed': pack(Y), 'packed3': pack(Y3),
                          'cases': CASES, 'records': records}))


if __name__ == '__main__':
    main()
