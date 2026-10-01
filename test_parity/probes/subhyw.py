"""
Probe SUBHYW, the horizontal tail's subsonic roll damping, and save the
fixture.

Run from the repository root: ``python test_parity/probes/subhyw.py``.
The cases run in one program (INTEP3 keeps its chart pair between calls).
"""

import math
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import probe  # noqa: E402
from probe import assign, parse_records, save  # noqa: E402

ROUTINES = ['subhyw', 'intep3', 'tlip2x', 'tlip1x', 'yup', 'switch',
            'glook', 'quad', 'tbfunx', 'interx', 'tlin1x', 'tlinex',
            'tlin3x']
COMMONS = """\
      COMMON /OVERLY/ NLOG,NMACH,IM,NALPHA
      COMMON /OPTION/ SR,CBARR,RUFF,BLREF
      COMMON /CONSNT/ PI,DR,UNUSED,RAD
      COMMON /FLGTCD/ FLC(93)
      COMMON /SYNTSS/ SYNA(19)
      COMMON /HTDATA/ A(195),B(49)
      COMMON /HTI/    WINGIN(77)
      COMMON /IHT/    PWING, WING(380)
      COMMON /BDATA/  BD(300),DYN(213)
      COMMON /WHAERO/ DUM1(106),C(51)
      COMMON /SBETA/  DUM2(351),STB(135)
      COMMON /FLOLOG/ FLTC,OPTI,BO,WGPL,WGSC,SYNT,HTPL,HTSC,VTPL,VTSC,
     1                HEAD,PRPOWR,JETPOW,LOASRT,TVTPAN,
     2                SUPERS,SUBSON,TRANSN,HYPERS,
     3                SYMFP,ASYFP,TRIMC,TRIM
      LOGICAL FLTC,OPTI,BO,WGPL,WGSC,SYNT,HTPL,HTSC,VTPL,VTSC,
     1        HEAD,PRPOWR,JETPOW,LOASRT,TVTPAN,
     2        SUPERS,SUBSON,TRANSN,HYPERS,
     3        SYMFP,ASYFP,TRIMC,TRIM
"""
ALPHA = [-4.0, 0.0, 4.0, 8.0, 12.0]
CL = {'straddle': [-0.3, 0.02, 0.35, 0.66, 0.9],
      'positive': [0.1, 0.4, 0.7, 0.95, 1.1],
      'negative': [-1.0, -0.8, -0.55, -0.3, -0.1]}
CD = [0.03, 0.02, 0.025, 0.04, 0.07]


def case(mach, taper, cl='straddle', sweep=35.0, ar=4.0, htpl=True,
         bo=False, sspne=5.0):
    return {'mach': mach, 'taper': taper, 'cl': cl, 'sweep': sweep,
            'ar': ar, 'htpl': htpl, 'bo': bo, 'sspne': sspne,
            'beta': math.sqrt(1. - mach ** 2)}


CASES = [case(0.3, 0.0), case(0.5, 0.2), case(0.6, 0.4, 'positive'),
         case(0.7, 0.6, 'negative'), case(0.4, 0.9, sweep=10.0, ar=6.0),
         case(0.8, 1.0, sweep=55.0, ar=3.0), case(0.2, 0.3, 'positive'),
         case(0.5, 0.5, htpl=True, bo=True, sspne=3.0),
         case(0.6, 0.75, bo=True, sspne=5.5)]


def driver():
    lines = ['      PROGRAM PROBE', COMMONS.rstrip('\n'),
             '      PI=3.141592654D0', '      UNUSED=1.D-30',
             '      RAD=57.2957795D0', '      IM=1', '      SR=320.0D0',
             '      BLREF=36.0D0', f'      NALPHA={len(ALPHA)}']
    lines += [assign(f'FLC({23 + k})', v) for k, v in enumerate(ALPHA)]
    for n, c in enumerate(CASES):
        lines.append(f"      WRITE(6,'(A,I4)') 'CASE',{n}")
        lines += [f'      DO {100 + n} K=1,380', '      WING(K)=K*0.001D0',
                  '      IF(K.LE.213) DYN(K)=K*0.002D0',
                  f'{100 + n:5d} CONTINUE']
        lines.append(f"      HTPL=.{'TRUE' if c['htpl'] else 'FALSE'}.")
        lines.append(f"      BO=.{'TRUE' if c['bo'] else 'FALSE'}.")
        lines.append(assign('FLC(3)', c['mach']))
        for k, v in enumerate(CL[c['cl']]):
            lines.append(assign(f'WING({21 + k})', v))
        for k, v in enumerate(CD):
            lines.append(assign(f'WING({1 + k})', v))
        for k, v in {4: 60.0, 27: c['taper'], 64: c['sweep'],
                     68: math.tan(math.radians(c['sweep'])),
                     120: c['ar']}.items():
            lines.append(assign(f'A({k})', v))
        lines += [assign('B(2)', c['beta']),
                  assign('WINGIN(3)', c['sspne']), assign('WINGIN(4)', 6.5),
                  assign('WINGIN(21)', 0.105), assign('SYNA(5)', 1.2),
                  assign('SYNA(7)', 2.0), assign('STB(122)', 3.0)]
        lines.append('      CALL SUBHYW')
        lines.append("      WRITE(6,'(A,60ES25.16)') 'DYN',(DYN(K),"
                     'K=40,99)')
        lines.append("      WRITE(6,'(A,20ES25.16)') 'HT',(WING(K),"
                     'K=281,300)')
    lines += ['      STOP', '      END', '']
    return '\n'.join(lines)


def main():
    records = parse_records(probe.run('subhyw', driver(), ROUTINES))
    print(save('subhyw', {'cases': CASES, 'alpha': ALPHA, 'cl': CL,
                          'cd': CD, 'records': records}))


if __name__ == '__main__':
    main()
