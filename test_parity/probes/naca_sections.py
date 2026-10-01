"""
Probe the NACA section coordinate routines (COORD1, COORD4, COORD5,
COORD6, CORD4M, CORD5M, XYCORD) and save the fixture.

Run from the repository root: ``python test_parity/probes/naca_sections.py``.
"""

import math
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from probe import assign, parse_records, run, save  # noqa: E402

ROUTINES = ['coord1', 'coord4', 'coord5', 'coord6', 'cord4m', 'cord5m',
            'xycord', 'sleq', 'arccos', 'tbfunx', 'quad']

COMMONS = """\
      COMMON /CONSNT/ PI,DEG,UNUSED,RAD
      COMMON /IBODY/  PB,NACA(80)
      COMMON /IWING/  PW, X(60)
      COMMON / IHT /  PHT, XU(60),XL(60),YU(60),YL(60)
      COMMON / IVT /  PVT, YUN(60),YLN(60)
      COMMON / IBW /  PBW,L,I,J,K,II,JJ,KK,III,JJJ
      COMMON / IBH /  PBH, THN(60),CAM(60)
      COMMON /IBWH/   PBWH,ALPHAI,ALPHAO,AII
      COMMON /IBWHV/  PBWHV, RHO,T,DELTAY,XOVC,TOVC,ZM,ZP
      REAL NACA
"""

X = [0.0, .0025, .005, .0125, .025, .05, .075, .1, .15, .2, .25, .3, .35,
     .4, .45, .5, .55, .6, .65, .7, .75, .8, .85, .9, .95, 1.0]

DIGITS = ('i', 'j', 'k', 'ii', 'jj', 'kk', 'iii', 'jjj')


def case(routine, **digits):
    d = {k: digits.get(k, 0) for k in DIGITS}
    return {'routine': routine, 'digits': d, 'x': X}


def cases():
    yu = [0.06 * math.sqrt(x) * (1 - x) + 0.02 * x * (1 - x) for x in X]
    yl = [-0.05 * math.sqrt(x) * (1 - x) + 0.01 * x * (1 - x) for x in X]
    return [
        case('COORD4', i=2, j=4, k=1, ii=2),             # 2412
        case('COORD4', i=0, j=0, k=1, ii=2),             # 0012
        case('COORD4', i=4, j=4, k=1, ii=5),             # 4415
        case('COORD5', i=2, j=3, k=0, ii=1, jj=2),       # 23012
        case('COORD5', i=2, j=3, k=1, ii=1, jj=2),       # 23112, reflexed
        case('COORD5', i=2, j=4, k=0, ii=1, jj=5),       # 24015
        case('CORD4M', i=0, j=0, k=1, ii=2, kk=6, iii=4),     # 0012-64
        case('CORD4M', i=2, j=4, k=1, ii=2, kk=6, iii=3),     # 2412-63
        case('CORD5M', i=2, j=3, k=0, ii=1, jj=2, iii=3, jjj=3),
        case('CORD5M', i=2, j=3, k=1, ii=1, jj=2, iii=6, jjj=4),
        case('COORD1', j=6, ii=2, jj=1, kk=2),            # 16-212
        case('COORD1', j=6, ii=0, jj=0, kk=9),            # 16-009
        case('COORD6', j=5, jj=2, kk=1, iii=2),           # 65-212
        case('COORD6', j=4, ii=1, jj=0, kk=1, iii=0),     # 64A010
        case('COORD6', j=3, jj=4, kk=1, iii=5, jjj=5),    # 63-415 a=0.5
        case('COORD6', j=6, jj=2, kk=1, iii=5, jjj=8),    # 66-215 a=0.8
        case('COORD6', j=5, ii=1, jj=2, kk=1, iii=2),     # 65A212
        {'routine': 'XYCORD', 'digits': {}, 'x': X, 'yu': yu, 'yl': yl},
    ]


def driver(all_cases) -> str:
    lines = ['      PROGRAM PROBE', COMMONS.rstrip('\n'),
             '      PI=3.141592654', '      DEG=0.01745329',
             '      UNUSED=1.E-30', '      RAD=57.2957795']
    for n, c in enumerate(all_cases):
        lines.append(f"      WRITE(6,'(A,I4)') 'CASE',{n}")
        lines.append(f'      DO {3000 + n} M=1,60')
        for arr in ('XU', 'XL', 'YUN', 'YLN', 'THN', 'CAM', 'X'):
            lines.append(f'         {arr}(M)=0.')
        lines.append(f' {3000 + n} CONTINUE')
        lines.append(f"      L={len(c['x'])}")
        for m, x in enumerate(c['x']):
            lines.append(assign(f'X({m + 1})', x))
        if c['routine'] == 'XYCORD':
            for m, (a, b) in enumerate(zip(c['yu'], c['yl'])):
                lines.append(assign(f'YU({m + 1})', a))
                lines.append(assign(f'YL({m + 1})', b))
            lines.append('      CALL XYCORD(0,0)')
        else:
            for name, v in c['digits'].items():
                lines.append(f'      {name.upper()}={v}')
            lines.append(f"      CALL {c['routine']}")
        for arr in ('XU', 'XL', 'YUN', 'YLN', 'THN', 'CAM'):
            lines.append(f"      WRITE(6,'(A,60ES25.16)') '{arr}',")
            lines.append(f'     1({arr}(M),M=1,L)')
        lines.append("      WRITE(6,'(A,7ES25.16)') 'S',RHO,T,ZM,ZP,ALPHAI,"
                     "ALPHAO,")
        lines.append('     1AII')
    lines += ['      STOP', '      END', '']
    return '\n'.join(lines)


def main():
    all_cases = cases()
    records = parse_records(run('naca_sections', driver(all_cases),
                                ROUTINES))
    print(save('naca_sections', [{'inputs': c, 'outputs': r}
                                 for c, r in zip(all_cases, records)]))


if __name__ == '__main__':
    main()
