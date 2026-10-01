"""
Probe SPRYAW (with DFLCON), the supersonic control roll and yaw, and save
the fixture.

Run from the repository root: ``python test_parity/probes/spryaw.py``.
"""

import math
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from probe import assign, parse_records, run, save  # noqa: E402

ROUTINES = ['spryaw', 'dflcon', 'arccos', 'tbfunx', 'quad', 'tlinex',
            'tlin1x', 'glook', 'switch']

COMMONS = """\
      COMMON /OVERLY/ NLOG,NMACH,IM,NALPHA
      COMMON /OPTION/ SR,CBARR,RUFF,BLREF
      COMMON /FLGTCD/ FLC(93)
      COMMON /WINGI/  WINGIN(77)
      COMMON /WINGD/  A(195)
      COMMON /HTI/    HTIN(131)
      COMMON /CONSNT/ PI,DR,UNUSED,RAD
      COMMON /FLAPIN/ F(69)
      COMMON /IHT/    PHT,HT(380)
      COMMON /IBODY/  PBODY,BODY(400)
      COMMON /IWING/  PWING,WING(400)
      COMMON /POWR/   SPR(59)
      COMMON /WBHCAL/ STP(155)
      COMMON /HTDATA/ AHT(195),BHT(49)
"""

ALPHA = [-4.0, 2.0, 8.0]


def case(stype=4.0, mach=2.0, aloci=3.0, aloco=7.0, cfi=2.0, cfo=1.5):
    f = {11: 0.02, 12: cfi, 13: cfo, 14: aloci, 15: aloco, 16: 3.0,
         17: 1.0, 18: stype}
    for j, (dl, dr, ds) in enumerate(((5.0, -5.0, 0.02), (12.0, -8.0, 0.05),
                                      (-10.0, 6.0, 0.1))):
        f[19 + j], f[29 + j], f[39 + j] = dl, dr, ds
    return {
        'alpha': ALPHA, 'mach': mach, 'f': f,
        'win': {2: 0.0, 3: 9.0, 4: 10.0},
        'a': {3: 150.0, 26: 0.4, 28: 0.5, 62: 0.3, 70: 15.0, 76: -5.0,
              77: -0.087},
        'htin': {3: 4.0, 4: 5.0}, 'sh': 30.0, 'cnahs': 0.045,
        'cladeg': 0.052, 'ivbh': [0.0, 0.2, 0.5],
        'gamma': [0.0, 0.3, 0.6], 'sref': 150.0, 'blref': 20.0,
        'spr': [0.011 * (k + 1) for k in range(59)], 'sae025': 0.0,
    }


def cases():
    return [
        case(),                                  # flap, tapered
        case(cfo=2.0, aloco=9.0),                # untapered, at the tip
        case(stype=1.0),                         # spoiler
        case(stype=5.0),                         # all-moving tail
        case(stype=3.0),                         # returns at once
        case(mach=1.05, aloci=1.0, aloco=2.0, cfi=4.0, cfo=1.0),
        case(mach=3.2, aloci=6.0, aloco=8.9),
        case(stype=1.0, mach=1.6, aloci=7.0, aloco=8.0),
    ]


def driver(all_cases) -> str:
    lines = ['      PROGRAM PROBE', COMMONS.rstrip('\n'),
             '      PI=3.141592654', '      DR=0.01745329',
             '      UNUSED=1.E-30', '      RAD=57.2957795', '      IM=1']
    for n, c in enumerate(all_cases):
        lines.append(f"      WRITE(6,'(A,I4)') 'CASE',{n}")
        lines.append(f'      DO {3000 + n} K=1,260')
        lines.append('         WING(K)=0.')
        lines.append('         BODY(K)=0.')
        lines.append('         HT(K)=0.')
        lines.append(f' {3000 + n} CONTINUE')
        lines.append(f"      NALPHA={len(c['alpha'])}")
        for j, al in enumerate(c['alpha']):
            lines.append(assign(f'FLC({23 + j})', al))
            lines.append(assign(f'STP({133 + j})', c['ivbh'][j]))
            lines.append(assign(f'STP({111 + j})', c['gamma'][j]))
        lines.append(assign('FLC(3)', c['mach']))
        for name, key in (('F', 'f'), ('WINGIN', 'win'), ('A', 'a'),
                          ('HTIN', 'htin')):
            for k, v in c[key].items():
                lines.append(assign(f'{name}({k})', v))
        for name, key in (('AHT(3)', 'sh'), ('HT(101)', 'cnahs'),
                          ('WING(101)', 'cladeg'), ('SR', 'sref'),
                          ('BLREF', 'blref')):
            lines.append(assign(name, c[key]))
        for k, v in enumerate(c['spr']):
            lines.append(assign(f'SPR({k + 1})', v))
        lines.append('      CALL SPRYAW')
        for tag, text in (('SPR', 'SPR'), ('HT', '(HT(K),K=211,230)'),
                          ('BODY', '(BODY(K),K=201,230)'),
                          ('WING', '(WING(K),K=201,230)')):
            lines.append(f"      WRITE(6,'(A,59ES25.16)') '{tag}',{text}")
    lines += ['      STOP', '      END', '']
    return '\n'.join(lines)


def main():
    all_cases = cases()
    records = parse_records(run('spryaw', driver(all_cases), ROUTINES))
    payload = []
    for c, r in zip(all_cases, records):
        for key in ('f', 'win', 'a', 'htin'):
            c[key] = {str(k): v for k, v in c[key].items()}
        payload.append({'inputs': c, 'outputs': r})
    print(save('spryaw', payload))


if __name__ == '__main__':
    main()
