"""
Probe DFLCON, the supersonic control derivatives, and save the fixture.

Run from the repository root: ``python test_parity/probes/dflcon.py``.
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from probe import assign, parse_records, run, save  # noqa: E402

ROUTINES = ['dflcon', 'arccos']

COMMONS = """\
      COMMON /CONSNT/ PI,DR,UNUSED,RAD
      COMMON /POWR/   SPR(59)
      LOGICAL INBORD,TAPERD
"""


def cases():
    out = []
    for taperd, fl in ((True, 0.6), (True, 0.85), (True, 1.3), (False, 2.5),
                       (False, 4.0), (False, 1.6)):
        for inbord in (True, False):
            for a, d in ((0.3, 0.1), (0.5, -0.2), (-0.2, -0.5),
                         (0.7, 0.45), (0.1, 0.05)):
                out.append({'a': a, 'd': d, 'fl': fl, 'beta': 1.3,
                            'inbord': inbord, 'taperd': taperd})
    out.append({'a': 0.3, 'd': 1.2, 'fl': 0.7, 'beta': 1.1, 'inbord': True,
                'taperd': True})                      # |D| > 1: no result
    out.append({'a': 0.3, 'd': 0.1, 'fl': 0.0, 'beta': 2.0,
                'inbord': True, 'taperd': True})       # FLD = 0
    return out


def driver(all_cases) -> str:
    lines = ['      PROGRAM PROBE', COMMONS.rstrip('\n'),
             '      PI=3.141592654', '      DR=0.01745329',
             '      UNUSED=1.E-30', '      RAD=57.2957795']
    for n, c in enumerate(all_cases):
        lines.append(f"      WRITE(6,'(A,I4)') 'CASE',{n}")
        for k in (12, 13, 17, 18):
            lines.append(f'      SPR({k})=-9.')
        lines.append(assign('AA', c['a']))
        lines.append(assign('DD', c['d']))
        lines.append(assign('FL', c['fl']))
        lines.append(assign('BETA', c['beta']))
        lines.append(f"      INBORD=.{'TRUE' if c['inbord'] else 'FALSE'}.")
        lines.append(f"      TAPERD=.{'TRUE' if c['taperd'] else 'FALSE'}.")
        lines.append('      CALL DFLCON(AA,DD,FL,BETA,INBORD,TAPERD)')
        lines.append("      WRITE(6,'(A,4ES25.16)') 'R',SPR(12),SPR(13),")
        lines.append('     1SPR(17),SPR(18)')
    lines += ['      STOP', '      END', '']
    return '\n'.join(lines)


def main():
    all_cases = cases()
    records = parse_records(run('dflcon', driver(all_cases), ROUTINES))
    print(save('dflcon', [{'inputs': c, 'outputs': r}
                          for c, r in zip(all_cases, records)]))


if __name__ == '__main__':
    main()
