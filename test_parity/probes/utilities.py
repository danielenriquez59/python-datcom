"""
Probe the small numerical utilities and save the fixture: ARCSIN, AREA1,
DET4, SLEQ, QUADIN, MACH2, SIMUL2, TLINVS and INTER3.

Run from the repository root: ``python test_parity/probes/utilities.py``.
"""

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from probe import assign, parse_records, run, save  # noqa: E402

ROUTINES = ['arcsin', 'area1', 'det4', 'sleq', 'quadin', 'mach2', 'simul2',
            'simul4', 'intkbw',
            'tlinvs', 'inter3', 'tbfunx', 'quad', 'tlinex', 'tlin1x',
            'glook', 'switch']

HEADER = """\
      PROGRAM PROBE
      COMMON /CONSNT/ PI,DEG,UNUSED,RAD
      DIMENSION X(4),Y(4),A16(16),A(6,7),XS(6),B(6),YQ(20)
      DIMENSION XC(12),C1(12),C2(12)
      DIMENSION T1(4),T2(3),TY(12),Q(3)
      DIMENSION RA1(3),RA2(2),RY(6),RB1(3),RB2(2),RYB(6)
      DIMENSION RC1(3),RC2(2),RYC(6),RD1(3),RD2(2),RYD(6)
      DIMENSION RE1(3),RE2(2),RYE(6)
      DIMENSION EQ(4),UNK(4)
      PI=3.141592654
      DEG=0.01745329
      UNUSED=1.E-30
      RAD=57.2957795
"""


def arrays(name, values):
    return [assign(f'{name}({k + 1})', v) for k, v in enumerate(values)]


def cases():
    rng = np.random.default_rng(99)
    c = {
        'arcsin': [0.0, 0.3, -0.7, 0.999999, 1.0, -1.0, 1.2],
        'area1': [
            {'x': [0, 4, 3, 0.5], 'y': [0, 0, 2, 2.5], 'nsum': 3},
            {'x': [0, 4, 3, 0.5], 'y': [0, 0, 2, 2.5], 'nsum': 4},
            {'x': [1, 5, 2, -1], 'y': [0, 1, 4, 3], 'nsum': 6},
            {'x': [1, 5, 2, -1], 'y': [0, 1, 4, 3], 'nsum': 5}],
        'det4': [list(rng.normal(size=16)) for _ in range(3)] +
                [list(np.arange(1.0, 17.0))],
        'sleq': [
            {'a': rng.normal(size=(4, 4)).tolist(), 'b': [1., 2., 3., 4.]},
            {'a': [[0., 2., 1.], [1., 1., 1.], [2., 0., 3.]],
             'b': [3., 3., 5.]},                        # zero first pivot
            {'a': [[1., 2.], [2., 4.]], 'b': [1., 2.]},  # singular
            {'a': [[2., 1., 0., 0.], [0., 0., 3., 1.], [1., 0., 0., 2.],
                   [0., 1., 1., 1.]], 'b': [1., 0., 2., 1.]}],
        'quadin': [{'y': list(np.sin(np.linspace(0, 2, n))), 'h': h}
                   for n, h in ((1, .1), (2, .3), (3, .2), (4, .25),
                                (5, .1), (6, .2), (7, .15), (8, .1),
                                (9, .3), (13, .05), (5, 0.0))],
        'mach2': [0.0, -1.0, 5.0, 26.38, 60.0, 102.3, 129.9, 131.0],
        'simul2': [
            {'x': [0, 1, 2, 3, 4], 'c1': [0, 1, 2, 3, 4],
             'c2': [3, 2.5, 2.2, 2.1, 2.0]},
            {'x': [0, 1, 2, 3], 'c1': [1, 1, 1, 1], 'c2': [0, 0.5, 1, 2]},
            {'x': [0, 1, 2, 3], 'c1': [1, 2, 3, 4], 'c2': [0, 0, 0, 0]},
            {'x': [0, 2, 4, 6, 8, 10], 'c1': [0.0, 0.5, 1.8, 2.2, 2.4, 2.5],
             'c2': [2.0, 1.9, 1.7, 1.2, 0.4, -0.5]}],
    }
    # TLINVS on a table falling with X1 (columns) and varying with X2.
    t1 = [1.0, 2.0, 4.0, 8.0]
    t2 = [0.0, 0.5, 1.0]
    ty = [[10. / (a * (1 + .3 * b)) for a in t1] for b in t2]   # (x2, x1)
    c['tlinvs'] = {'x1': t1, 'x2': t2, 'y': ty, 'queries': [
        (0.25, 3.0), (0.0, 20.0), (0.0, 0.5), (0.5, 4.0), (1.0, 2.0),
        (1.0, 12.0), (1.0, 1.0), (0.75, 1.7), (-0.2, 3.0), (0.3, 9.5)]}
    tabs = []
    for k in range(5):
        x1, x2 = [0.0, 1.0, 2.0], [0.0, 3.0]
        y = (np.array([[1.0, 2.0, 4.0], [2.0, 3.0, 7.0]]) * (1 + 0.1 * k)
             ).tolist()
        tabs.append({'x1': x1, 'x2': x2, 'y': y})
    c['inter3'] = {'tables': tabs, 'queries': [
        (0.5, 1.0, 5e4), (0.5, 1.0, 1e5), (1.5, 2.0, 3e5), (1.5, 2.0, 1e6),
        (0.2, 0.5, 5e6), (1.0, 3.0, 2e7), (1.7, 1.0, 5e8), (0.4, 2.0, 3e9)]}
    c['simul4'] = [{'coff': list(rng.normal(size=16)),
                    'eq': list(rng.normal(size=4))} for _ in range(3)]
    c['intkbw'] = [(1.6, 45.0, 10.0, 3.0, 5.0), (1.3, 60.0, 10.0, 3.0, 5.0),
                   (2.5, 30.0, 8.0, 4.0, 1.0), (1.8, 70.0, 6.0, 2.0, 20.0),
                   (1.1, 50.0, 12.0, 3.5, 0.5)]
    return c


def driver(c) -> str:
    lines = [HEADER.rstrip('\n')]
    n = 0

    def case():
        nonlocal n
        lines.append(f"      WRITE(6,'(A,I4)') 'CASE',{n}")
        n += 1
    for v in c['arcsin']:
        case()
        lines.append(assign('V', v))
        lines.append('      R=ARCSIN(V)')
        lines.append("      WRITE(6,'(A,ES25.16)') 'R',R")
    for v in c['area1']:
        case()
        lines += arrays('X', v['x']) + arrays('Y', v['y'])
        lines.append(f"      CALL AREA1(X,Y,R,{v['nsum']})")
        lines.append("      WRITE(6,'(A,ES25.16)') 'R',R")
    for v in c['det4']:
        case()
        lines += arrays('A16', v)
        lines.append('      CALL DET4(A16,R)')
        lines.append("      WRITE(6,'(A,ES25.16)') 'R',R")
    for v in c['sleq']:
        case()
        m = len(v['b'])
        lines.append(f'      DO {5000 + n} K=1,6')
        lines.append(f' {5000 + n} XS(K)=-7.')
        for i in range(m):
            for j in range(m):
                lines.append(assign(f'A({i + 1},{j + 1})', v['a'][i][j]))
        lines += arrays('B', v['b'])
        # SLEQ indexes A(N,M) with its own leading dimension N.
        lines.append(f'      CALL SLEQPR(A,XS,B,{m})')
        lines.append(f"      WRITE(6,'(A,6ES25.16)') 'R',(XS(K),K=1,{m})")
    for v in c['quadin']:
        case()
        lines += arrays('YQ', v['y'])
        lines.append(assign('H', v['h']))
        lines.append(f"      CALL QUADIN(YQ,{len(v['y'])},H,R)")
        lines.append("      WRITE(6,'(A,ES25.16)') 'R',R")
    for v in c['mach2']:
        case()
        lines.append(assign('V', v))
        lines.append('      CALL MACH2(V,R,IER)')
        lines.append("      WRITE(6,'(A,ES25.16,I3)') 'R',R,IER")
    for v in c['simul2']:
        case()
        npt = len(v['x'])
        lines += arrays('XC', v['x']) + arrays('C1', v['c1']) + \
            arrays('C2', v['c2'])
        lines.append(f'      CALL SIMUL2(XC,C1,C2,{npt},R,S)')
        lines.append("      WRITE(6,'(A,2ES25.16)') 'R',R,S")
    t = c['tlinvs']
    flat = [t['y'][i2][i1] for i1 in range(4) for i2 in range(3)]
    for q in t['queries']:
        case()
        lines += arrays('T1', t['x1']) + arrays('T2', t['x2']) + \
            arrays('TY', flat)
        lines.append(assign('XA2', q[0]))
        lines.append(assign('ZA', q[1]))
        lines.append('      CALL TLINVS(T1,T2,TY,4,3,R,XA2,ZA)')
        lines.append("      WRITE(6,'(A,ES25.16)') 'R',R")
    names = [('RA1', 'RA2', 'RY'), ('RB1', 'RB2', 'RYB'),
             ('RC1', 'RC2', 'RYC'), ('RD1', 'RD2', 'RYD'),
             ('RE1', 'RE2', 'RYE')]
    for q in c['inter3']['queries']:
        case()
        for (n1, n2, ny), tab in zip(names, c['inter3']['tables']):
            flat = [tab['y'][i2][i1] for i1 in range(3) for i2 in range(2)]
            lines += arrays(n1, tab['x1']) + arrays(n2, tab['x2']) + \
                arrays(ny, flat)
        lines.append(assign('G1', q[0]))
        lines.append(assign('G2', q[1]))
        lines.append(assign('RL', q[2]))
        lines += [
            '      CALL INTER3(G1,G2,RL,RA1,RA2,RY,3,2,Q,RB1,RB2,RYB,3,2,Q,',
            '     1  RC1,RC2,RYC,3,2,Q,RD1,RD2,RYD,3,2,Q,RE1,RE2,RYE,3,2,',
            '     2  Q,R)',
            "      WRITE(6,'(A,ES25.16)') 'R',R"]
    for v in c['simul4']:
        case()
        lines += arrays('A16', v['coff']) + arrays('EQ', v['eq'])
        lines.append('      CALL SIMUL4(A16,EQ,UNK)')
        lines.append("      WRITE(6,'(A,4ES25.16)') 'R',UNK")
    for v in c['intkbw']:
        case()
        for name, value in zip(('XM', 'OLE', 'CR', 'D', 'DX'), v):
            lines.append(assign(name, value))
        lines.append('      CALL INTKBW(XM,OLE,CR,D,DX,R,S)')
        lines.append("      WRITE(6,'(A,2ES25.16)') 'R',R,S")
    lines += ['      STOP', '      END',
              '      SUBROUTINE SLEQPR(A,X,B,N)',
              '      DIMENSION A(6,7),X(6),B(6),W(6,7)',
              '      DO 10 I=1,N',
              '      DO 10 J=1,N',
              '   10 W(I,J)=A(I,J)',
              '      CALL SLEQ2(W,X,B,N)',
              '      RETURN', '      END',
              '      SUBROUTINE SLEQ2(W,X,B,N)',
              '      DIMENSION W(6,7),X(6),B(6),P(42)',
              '      K=0',
              '      DO 20 J=1,N+1',
              '      DO 20 I=1,N',
              '      K=K+1',
              '   20 P(K)=W(I,J)',
              '      CALL SLEQ(P,X,B,N,N+1)',
              '      RETURN', '      END', '']
    return '\n'.join(lines)


def main():
    c = cases()
    records = parse_records(run('utilities', driver(c), ROUTINES))
    print(save('utilities', {'inputs': c,
                             'outputs': [r.get('R') for r in records],
                             'text': [r.get('_text') for r in records]}))


if __name__ == '__main__':
    main()
