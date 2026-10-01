"""
Probe SUPDRG, the supersonic wing drag, and save the fixture.

The cases run in one program, in order: SUPDRG's local ``RACH`` is saved
between calls (the build uses ``-fno-automatic``, as the program does), so
a case without roughness reads the Mach number the previous case left.

Run from the repository root: ``python test_parity/probes/supdrg.py``.
"""

import math
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from probe import assign, parse_records, run, save  # noqa: E402

ROUTINES = ['supdrg', 'fig26', 'interx', 'tlin1x', 'tlinex', 'tlin3x',
            'glook', 'switch', 'tbfunx', 'quad']

COMMONS = """\
      COMMON /OVERLY/ NLOG,NMACH,N,NALPHA,IG
      COMMON /CONSNT/ PI,DEG,UNUSED,RAD
      COMMON /FLGTCD/ FLC(93)
      COMMON /OPTION/ SR,CBARR,RUFF,BLREF
      COMMON /WINGD/  A(195),B(48)
      COMMON /WINGI/  WINGIN(77)
      COMMON /SUPWH/  SLG(141)
      DIMENSION WT(4)
      DATA WT /4HSTRA,4HDOUB,4HCRAN,4HCURV/
"""

UNUSED = 1.0e-30


def case(mach=2.0, kind=1, ksharp=UNUSED, ruff=1.6e-4, sweep_i=45.0,
         sweep_o=30.0, cbari=8.0, cbaro=4.0, tceff=0.05):
    ti, to = math.tan(math.radians(sweep_i)), math.tan(math.radians(sweep_o))
    return {
        'mach': mach, 'sref': 300.0, 'roughness': ruff, 'sbw': 140.0,
        'win': {1: 2.5, 3: 13.0, 5: 5.0, 6: 11.0, 15: kind, 62: 0.004,
                63: 0.003, 70: tceff, 71: ksharp},
        'a': {1: 90.0, 2: 60.0, 3: 150.0, 7: 2.6, 10: 10.5, 15: cbari,
              17: cbaro, 18: 0.9, 61: math.cos(math.atan(ti)), 62: ti,
              85: math.cos(math.atan(to)), 86: to, 129: 2.0e6},
        'stale': {'rlcoff': 1.5e7, 'cfi': 0.0031},
    }


def cases():
    return [
        case(),                                      # straight, round
        case(mach=1.3, ksharp=5.0),                  # sharp, BOVERT < 1
        case(mach=3.5, ksharp=5.33),                 # sharp, Mach cap
        case(ruff=0.0),                              # no roughness: stale RACH
        case(kind=3),                                # cranked, round
        case(kind=2, ksharp=6.0, mach=1.5),          # double delta, sharp
        case(kind=3, cbari=6.0, cbaro=6.0),          # equal MACs: stale CFI
        case(kind=3, sweep_o=0.0),                   # TANLEO = 0 left 1e-5
        case(sweep_i=0.0, mach=1.1),                 # TANLEI = 0 restored
        case(ruff=0.05, mach=4.0),                   # cutoff below RN
        case(kind=3, ruff=0.0, mach=2.2),            # cranked, no roughness
        case(mach=1.02),                             # Fig. 58 lower end
    ]


def driver(all_cases) -> str:
    lines = ['      PROGRAM PROBE', COMMONS.rstrip('\n'),
             '      PI=3.141592654', '      DEG=0.01745329',
             '      UNUSED=1.E-30', '      RAD=57.2957795', '      N=1']
    for n, c in enumerate(all_cases):
        lines.append(f"      WRITE(6,'(A,I4)') 'CASE',{n}")
        lines.append(f'      DO {3000 + n} K=1,141')
        lines.append(f' {3000 + n} SLG(K)=0.')
        lines.append(assign('FLC(3)', c['mach']))
        lines.append(assign('SR', c['sref']))
        lines.append(assign('RUFF', c['roughness']))
        lines.append(assign('SLG(119)', c['sbw']))
        lines.append(assign('SLG(89)', c['stale']['rlcoff']))
        lines.append(assign('SLG(84)', c['stale']['cfi']))
        for k, v in c['win'].items():
            if k == 15:
                lines.append(f'      WINGIN(15)=WT({int(v)})')
            else:
                lines.append(assign(f'WINGIN({k})', v))
        for k, v in c['a'].items():
            lines.append(assign(f'A({k})', v))
        lines.append('      CALL SUPDRG')
        lines.append("      WRITE(6,'(A,141ES25.16)') 'SLG',")
        lines.append('     1(SLG(J),J=1,141)')
        lines.append("      WRITE(6,'(A,2ES25.16)') 'A',A(62),A(86)")
    lines += ['      STOP', '      END', '']
    return '\n'.join(lines)


def main():
    all_cases = cases()
    records = parse_records(run('supdrg', driver(all_cases), ROUTINES))
    payload = []
    for c, r in zip(all_cases, records):
        c['win'] = {str(k): v for k, v in c['win'].items()}
        c['a'] = {str(k): v for k, v in c['a'].items()}
        payload.append({'inputs': c, 'outputs': r})
    print(save('supdrg', payload))


if __name__ == '__main__':
    main()
