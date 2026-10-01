"""
Probe SUPCLD (wing) and SUPHLD (horizontal tail), the supersonic
CL-alpha-dot, and save the fixture.

Run from the repository root: ``python test_parity/probes/supcld.py``.
Each routine runs in its own program with its own COMMON blocks; every
word of the blocks it touches is filled with a marker first.
"""

import math
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import probe  # noqa: E402
from probe import assign, parse_records, save  # noqa: E402

ROUTINES = ['calca', 'arcsin', 'interx', 'tlin1x', 'tlinex', 'tlin3x',
            'glook', 'switch', 'quad', 'tlip3x', 'tlip2x', 'tlip1x', 'yup']


def case(mach, tanle, taper, ar=3.0):
    beta = math.sqrt(mach ** 2 - 1.)
    return {'mach': mach, 'sr': 320.0, 'cbarr': 8.5,
            'a': {3: 280.0, 7: ar, 9: 0.2, 10: 11.0, 16: 7.8, 27: taper,
                  29: 12.0, 58: math.degrees(math.atan(tanle)), 62: tanle},
            'wingin': {1: 3.0, 3: 12.5, 4: 14.0},
            'slg': {1: beta, 134: 0.42}}


CASES = [case(1.8, 1.0, 0.5), case(2.5, 0.3, 0.2), case(1.4, 2.0, 0.0),
         case(1.4, 2.0, 0.1), case(1.4, 2.0, 0.5), case(1.6, 3.0, 0.3, 2.0),
         case(1.3, 2.5, 0.2, 4.0), case(1.5, 0.0, 0.5)]

WING = """\
      COMMON /OVERLY/ NLOG,NMACH,IM,NALPHA
      COMMON /OPTION/ SR,CBARR,RUFF,BLREF
      COMMON /CONSNT/ PI,DR,UNUSED,RAD
      COMMON /FLGTCD/ FLC(93)
      COMMON /SUPWH/  SLG(141)
      COMMON /POWR/   DYN(213)
      COMMON /WINGD/  A(195)
      COMMON /IWING/  PWING, WING(400)
      COMMON /WINGI/  WINGIN(101)
"""
TAIL = """\
      COMMON /OVERLY/ NLOG,NMACH,IM,NALPHA
      COMMON /OPTION/ SR,CBARR,RUFF,BLREF
      COMMON /CONSNT/ PI,DR,UNUSED,RAD
      COMMON /FLGTCD/ FLC(93)
      COMMON /SUPWH/  XSLG(141),SLG(141)
      COMMON /BDATA/  XBD(300),DYN(213)
      COMMON /HTDATA/ A(195)
      COMMON /IHT/    PWING, WING(380)
      COMMON /HTI/    WINGIN(101)
"""


def driver(commons, call):
    lines = ['      PROGRAM PROBE', commons.rstrip('\n'),
             '      PI=3.141592654D0', '      UNUSED=1.D-30',
             '      RAD=57.2957795D0', '      IM=1']
    for n, c in enumerate(CASES):
        lines.append(f"      WRITE(6,'(A,I4)') 'CASE',{n}")
        lines += [f'      DO {100 + n} K=1,213', '      DYN(K)=K*0.01D0',
                  '      IF(K.LE.195) A(K)=K*0.03D0',
                  '      IF(K.LE.141) SLG(K)=K*0.02D0',
                  '      IF(K.LE.101) WINGIN(K)=K*0.04D0',
                  f'{100 + n:5d} CONTINUE', '      WING(241)=-7.0D0']
        lines.append(assign('FLC(3)', c['mach']))
        lines.append(assign('SR', c['sr']))
        lines.append(assign('CBARR', c['cbarr']))
        for block, name in (('a', 'A'), ('wingin', 'WINGIN'),
                            ('slg', 'SLG')):
            for k, v in c[block].items():
                lines.append(assign(f'{name}({k})', v))
        lines.append(assign('DYN(9)', 0.8))
        lines.append(assign('DYN(10)', 1.3))
        lines.append(f'      CALL {call}')
        lines.append("      WRITE(6,'(A,213ES25.16)') 'DYN',DYN")
        lines.append("      WRITE(6,'(A,195ES25.16)') 'A',A")
        lines.append("      WRITE(6,'(A,ES25.16)') 'CLAD',WING(241)")
    lines += ['      STOP', '      END', '']
    return '\n'.join(lines)


def main():
    out = {}
    for name, commons in (('supcld', WING), ('suphld', TAIL)):
        out[name] = parse_records(probe.run(name, driver(commons,
                                                         name.upper()),
                                            ROUTINES + [name]))
    print(save('supcld', {'cases': CASES, 'records': out}))


if __name__ == '__main__':
    main()
