"""
Probe BDAREA (with PTINT2 and AREA2), the body area in the horizontal
tail's Mach zone, and save the fixture.

Run from the repository root: ``python test_parity/probes/bdarea.py``.
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from probe import assign, parse_records, run, save  # noqa: E402

ROUTINES = ['bdarea', 'ptint2', 'area2']

COMMONS = """\
      COMMON /OVERLY/ NLOG,NMACH,I,NALPHA,IG
      COMMON /FLGTCD/ FLC(73)
      COMMON /CONSNT/ PI,DEG,UNUSED,RAD
      COMMON /HTI/    HTIN(154)
      COMMON /HTDATA/ AHT(195)
      COMMON /BODYI/  XNX,XB(20),SB(20),PB(20),RB(20)
      COMMON /SYNTSS/ SYNA(10)
      LOGICAL ABORT
"""

MACH_INDEX = 2
XB = [0.0, 2.0, 5.0, 10.0, 15.0, 20.0, 25.0, 30.0, 35.0, 40.0]
RB = [0.0, 0.8, 1.3, 1.8, 2.0, 2.0, 2.0, 1.9, 1.7, 1.5]


def case(mach=2.0, xh=30.0, zh=0.5, alih=0.0, chord=3.0, xb=None, rb=None):
    return {'mach': mach, 'xb': xb or XB, 'rb': rb or RB,
            'syna': {1: 20.0, 6: xh, 7: zh, 8: alih},
            'htin': {3: 4.0, 4: 5.0}, 'aht': {10: chord, 62: 0.7}}


CASES = [
    case(),
    case(mach=1.3),
    case(mach=3.0, xh=35.0),
    case(xh=36.0, chord=4.0),
    case(alih=3.0),
    case(zh=-1.0, alih=-2.0),
    case(mach=1.2, xh=20.0),
    # Found by sampling: each area shape on each side, and an abort.
    case(mach=3.0, xh=17.7, zh=-0.71, alih=-3.0, chord=5.0),
    case(mach=1.3, xh=21.4, zh=-1.99, chord=2.0),
    case(mach=1.6, xh=36.8, zh=-1.25, alih=2.0, chord=2.0),
    case(mach=3.0, xh=32.5, zh=1.49, alih=2.0, chord=2.0),
    case(mach=4.0, xh=33.9, zh=-0.76, alih=-3.0, chord=3.0),
    case(mach=4.0, xh=33.8, zh=1.25, chord=5.0),
    case(mach=4.0, xh=37.7, zh=-2.35, alih=2.0, chord=2.0),
]


def cases():
    return CASES


def driver(all_cases) -> str:
    lines = ['      PROGRAM PROBE', COMMONS.rstrip('\n'),
             '      PI=3.141592654', '      DEG=0.01745329',
             '      UNUSED=1.E-30', '      RAD=57.2957795',
             f'      I={MACH_INDEX}']
    m = MACH_INDEX
    for n, c in enumerate(all_cases):
        lines.append(f"      WRITE(6,'(A,I4)') 'CASE',{n}")
        lines.append(assign(f'HTIN({94 + m})', -1.0))
        lines.append(assign(f'HTIN({114 + m})', -2.0))
        lines.append(assign(f'HTIN({134 + m})', -3.0))
        lines.append(assign(f'FLC({m + 2})', c['mach']))
        lines.append(f"      XNX={len(c['xb'])}.")
        for k, (x, r) in enumerate(zip(c['xb'], c['rb'])):
            lines.append(assign(f'XB({k + 1})', x))
            lines.append(assign(f'RB({k + 1})', r))
        for k, v in c['syna'].items():
            lines.append(assign(f'SYNA({k})', v))
        for k, v in c['htin'].items():
            lines.append(assign(f'HTIN({k})', v))
        for k, v in c['aht'].items():
            lines.append(assign(f'AHT({k})', v))
        lines.append('      CALL BDAREA(ABORT)')
        lines.append('      IA=0')
        lines.append('      IF(ABORT) IA=1')
        lines.append("      WRITE(6,'(A,I2,3ES25.16)') 'R',IA,")
        lines.append(f'     1HTIN({114 + m}),HTIN({134 + m}),HTIN({94 + m})')
    lines += ['      STOP', '      END', '']
    return '\n'.join(lines)


def main():
    all_cases = cases()
    records = parse_records(run('bdarea', driver(all_cases), ROUTINES))
    payload = []
    for c, r in zip(all_cases, records):
        c = dict(c)
        for key in ('syna', 'htin', 'aht'):
            c[key] = {str(k): v for k, v in c[key].items()}
        payload.append({'inputs': c, 'outputs': r})
    print(save('bdarea', payload))


if __name__ == '__main__':
    main()
