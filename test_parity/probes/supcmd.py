"""
Probe SUPCMD (wing) and SUPHMD (horizontal tail), the supersonic
Cm-alpha-dot, and save the fixture.

Run from the repository root: ``python test_parity/probes/supcmd.py``.
The cases are SUPCLD's (test_parity/probes/supcld.py).
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
sys.path.insert(0, str(Path(__file__).resolve().parent))
import probe  # noqa: E402
from probe import assign, parse_records, save  # noqa: E402
from supcld import CASES, TAIL, WING  # noqa: E402

ROUTINES = ['interx', 'tlin1x', 'tlinex', 'tlin3x', 'glook', 'switch',
            'quad', 'tlip3x', 'tlip2x', 'tlip1x', 'yup']


def driver(commons, call):
    lines = ['      PROGRAM PROBE', commons.rstrip('\n'),
             '      PI=3.141592654D0', '      UNUSED=1.D-30',
             '      RAD=57.2957795D0', '      IM=1']
    for n, c in enumerate(CASES):
        lines.append(f"      WRITE(6,'(A,I4)') 'CASE',{n}")
        lines += [f'      DO {100 + n} K=1,213', '      DYN(K)=K*0.01D0',
                  '      IF(K.LE.195) A(K)=K*0.03D0',
                  '      IF(K.LE.141) SLG(K)=K*0.02D0',
                  f'{100 + n:5d} CONTINUE', '      WING(241)=-7.0D0',
                  '      WING(261)=3.5D0']
        lines.append(assign('FLC(3)', c['mach']))
        lines.append(assign('SR', c['sr']))
        lines.append(assign('CBARR', c['cbarr']))
        for k, v in c['a'].items():
            lines.append(assign(f'A({k})', v))
        for k, v in c['slg'].items():
            lines.append(assign(f'SLG({k})', v))
        lines.append(assign('A(173)', 1.25))
        lines.append(f'      CALL {call}')
        lines.append("      WRITE(6,'(A,213ES25.16)') 'DYN',DYN")
        lines.append("      WRITE(6,'(A,195ES25.16)') 'A',A")
        lines.append("      WRITE(6,'(A,2ES25.16)') 'RES',WING(241),"
                     'WING(261)')
    lines += ['      STOP', '      END', '']
    return '\n'.join(lines)


def main():
    out = {}
    for name, commons in (('supcmd', WING), ('suphmd', TAIL)):
        out[name] = parse_records(probe.run(name, driver(commons,
                                                         name.upper()),
                                            ROUTINES + [name]))
    print(save('supcmd', {'cases': CASES, 'records': out}))


if __name__ == '__main__':
    main()
