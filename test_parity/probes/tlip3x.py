"""
Probe TLIP3X (TLIN3X over a packed table) and INTEP3 (interpolation
between packed charts), and save the fixture.

Run from the repository root: ``python test_parity/probes/tlip3x.py``.
INTEP3 runs its cases in sequence: a LAMDA outside [0, 1] reuses the chart
pair of the call before.
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
sys.path.insert(0, str(Path(__file__).resolve().parent))
from probe import assign, parse_records, run, save  # noqa: E402
from tlip2x import X1, X2, Y, pack  # noqa: E402

ROUTINES = ['tlip3x', 'intep3', 'tlip2x', 'tlip1x', 'yup', 'switch',
            'glook', 'quad']
X3 = [0.0, 4.0]
Y3 = Y + [2.0 * v - 1.0 for v in Y]
POINTS = [(0.5, 0.25, 1.0), (1.5, 1.5, 3.0), (2.5, 0.3, 5.0),
          (-0.5, 2.5, -1.0), (1.0, 0.5, 4.0)]
MODES = [(0, 0, 0, 0, 0, 0), (1, 1, 1, 1, 1, 1), (2, 2, 2, 2, 2, 2)]
SCALE = {'A': 1.0, 'B': 2.0, 'C': 4.0, 'D': 8.0, 'E': 10.0}
LAMDAS = [0.0, 0.1, 0.25, 0.3, 0.5, 0.6, 0.75, 0.9, 1.0, 1.2, -0.1]


def chart(key):
    return [v * SCALE[key] for v in Y]


def driver():
    p3 = pack(Y3)
    lines = ['      PROGRAM PROBE',
             f'      DIMENSION X1(3),X2(4),X3(2),NY3({len(p3)}),NXX(7)',
             '      DIMENSION ROUT(2),MESS(2),NDX(7)']
    for key in SCALE:
        lines.append(f'      DIMENSION N{key}(6)')
    lines += [assign(f'X1({k + 1})', v) for k, v in enumerate(X1)]
    lines += [assign(f'X2({k + 1})', v) for k, v in enumerate(X2)]
    lines += [assign(f'X3({k + 1})', v) for k, v in enumerate(X3)]
    lines += [f'      NY3({k + 1})={w}' for k, w in enumerate(p3)]
    for key in SCALE:
        lines += [f'      N{key}({k + 1})={w}'
                  for k, w in enumerate(pack(chart(key)))]
    lines += ['      NXX(1)=3'] + [f'      NXX({k})=0' for k in range(2, 8)]
    for q1, q2, q3 in POINTS:
        for m in MODES:
            lines.append("      WRITE(6,'(A)') 'CASE'")
            lines.append(assign('A1', q1))
            lines.append(assign('A2', q2))
            lines.append(assign('A3', q3))
            lines.append('      CALL TLIP3X(X1,X2,X3,NY3,3,4,2,A1,A2,A3,YA,'
                         + ','.join(map(str, m)) + ',')
            lines.append('     1            MESS,1,ROUT)')
            lines.append("      WRITE(6,'(A,ES25.16)') 'T3',YA")
    for n2d in (4, 1):
        for lam in LAMDAS:
            lines.append("      WRITE(6,'(A)') 'CASE'")
            lines.append(assign('LAM', lam))
            lines.append('      CALL INTEP3(0.8D0,0.3D0,LAM,'
                         'X1,X2,NA,NXX,4,MESS,')
            lines.append('     1  X1,X2,NB,NXX,4,MESS,X1,X2,NC,NXX,4,MESS,')
            lines.append(f'     2  X1,X2,ND,NXX,{n2d},MESS,'
                         'X1,X2,NE,NXX,4,MESS,ANS)')
            lines.append("      WRITE(6,'(A,ES25.16)') 'I3',ANS")
    lines.insert(3, '      REAL LAM')
    lines += ['      STOP', '      END', '']
    return '\n'.join(lines)


def main():
    records = parse_records(run('tlip3x', driver(), ROUTINES))
    charts = {k.lower(): pack(chart(k)) for k in SCALE}
    print(save('tlip3x', {'x1': X1, 'x2': X2, 'x3': X3,
                          'packed3': pack(Y3), 'charts': charts,
                          'points': POINTS, 'modes': MODES,
                          'lamdas': LAMDAS, 'records': records}))


if __name__ == '__main__':
    main()
