"""
Probe CALCA, the supersonic wing acceleration derivatives, and save the
fixture.

Run from the repository root: ``python test_parity/probes/calca.py``.
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from probe import assign, parse_records, run, save  # noqa: E402

ROUTINES = ['calca', 'arcsin', 'interx', 'tlin1x', 'tlinex', 'tlin3x',
            'glook', 'switch', 'tbfunx', 'quad']

COMMONS = """\
      COMMON /OVERLY/ NLOG,NMACH,IM,NALPHA
      COMMON /CONSNT/ PI,DR,UNUSED,RAD
      COMMON /FLGTCD/ FLC(93)
      DIMENSION DYN(213),A(195),WINGIN(77)
"""


def cases():
    out = []
    for mach, sweep, ar, taper in ((1.2, 60., 2.5, 0.3), (1.5, 60., 2.0, 0.2),
                                   (2.0, 70., 1.8, 0.1), (1.4, 55., 3.0, 0.4),
                                   (2.5, 75., 1.5, 0.0), (1.1, 45., 4.0, 0.5),
                                   (1.8, 65., 2.2, 0.25), (3.0, 72., 1.2, 0.)):
        out.append({'mach': mach,
                    'a': {3: 180.0, 7: ar, 9: 12.0, 27: taper, 58: sweep},
                    'win': {1: 12.0 * taper, 3: 9.0, 4: 10.0}})
    return out


def driver(all_cases) -> str:
    lines = ['      PROGRAM PROBE', COMMONS.rstrip('\n'),
             '      PI=3.141592654', '      DR=0.01745329',
             '      UNUSED=1.E-30', '      RAD=57.2957795', '      IM=1']
    for n, c in enumerate(all_cases):
        lines.append(f"      WRITE(6,'(A,I4)') 'CASE',{n}")
        for k in (16, 22, 26, 27):
            lines.append(f'      DYN({k})=-5.')
        lines.append(assign('FLC(3)', c['mach']))
        for k, v in c['a'].items():
            lines.append(assign(f'A({k})', v))
        for k, v in c['win'].items():
            lines.append(assign(f'WINGIN({k})', v))
        lines.append('      CALL CALCA(DYN,A,WINGIN)')
        lines.append("      WRITE(6,'(A,4ES25.16)') 'R',DYN(16),DYN(22),"
                     "DYN(26),")
        lines.append('     1DYN(27)')
    lines += ['      STOP', '      END', '']
    return '\n'.join(lines)


def main():
    all_cases = cases()
    records = parse_records(run('calca', driver(all_cases), ROUTINES))
    payload = []
    for c, r in zip(all_cases, records):
        c['a'] = {str(k): v for k, v in c['a'].items()}
        c['win'] = {str(k): v for k, v in c['win'].items()}
        payload.append({'inputs': c, 'outputs': r})
    print(save('calca', payload))


if __name__ == '__main__':
    main()
