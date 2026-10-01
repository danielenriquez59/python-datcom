"""
Probe TRANJT, the hypersonic transverse jet sizing, and save the fixture.

Run from the repository root: ``python test_parity/probes/tranjt.py``.
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from probe import assign, parse_records, run, save  # noqa: E402

ROUTINES = ['tranjt', 'inter3', 'simul2', 'trapz', 'tbfunx', 'quad',
            'tlinex', 'tlin1x', 'glook', 'switch']

COMMONS = """\
      COMMON /OVERLY/ NLOG,NMACH,IM,NALPHA
      COMMON /CONSNT/ PI,DR,UNUSED,RAD
      COMMON /FLGTCD/ FLC(93)
      COMMON /IBODY/  PBODY,JET(200)
      COMMON /SUPDW/  JETA(26)
      COMMON /FLAPIN/ F(48)
      LOGICAL LAMNR(10)
      REAL JET,JETA
      EQUIVALENCE (LAMNR(1),F(39))
"""


def case(mach=8.0, nt=5, laminar=True, rl=1.0e6, pinf=0.1, phe=0.0,
         me=3.0, force=500.0):
    alpha = [2.0, 4.0, 6.0, 8.0, 10.0, 12.0, 14.0, 5.0, 7.0, 9.0][:nt]
    lam = ([laminar] * nt if isinstance(laminar, bool) else
           list(laminar))
    return {'mach': mach, 'rl': rl, 'pinf': pinf, 'nt': float(nt),
            'time': [0.5 * k for k in range(nt)],
            'fc': [force * (1.0 + 0.1 * k) for k in range(nt)],
            'alpha': alpha, 'me': me, 'isp': 250.0, 'span': 1.0,
            'phe': phe, 'gp': 1.4, 'cc': 0.95, 'l': 20.0,
            'laminar': lam}


def cases():
    return [
        case(),
        case(laminar=False),
        case(mach=12.0, nt=9, rl=5.0e6, phe=20.0),
        case(mach=6.0, laminar=[True, False, True, False, True]),
        case(mach=16.0, rl=2.0e7, me=4.0, force=2000.0),
        case(mach=10.0, laminar=False, nt=8, pinf=0.02),
    ]


def driver(all_cases) -> str:
    lines = ['      PROGRAM PROBE', COMMONS.rstrip('\n'),
             '      PI=3.141592654', '      DR=0.01745329',
             '      UNUSED=1.E-30', '      RAD=57.2957795', '      IM=1']
    for n, c in enumerate(all_cases):
        lines.append(f"      WRITE(6,'(A,I4)') 'CASE',{n}")
        lines.append(f'      DO {3000 + n} K=1,200')
        lines.append('         JET(K)=0.')
        lines.append('         IF(K.LE.26) JETA(K)=0.')
        lines.append(f' {3000 + n} CONTINUE')
        lines.append(assign('FLC(3)', c['mach']))
        lines.append(assign('FLC(43)', c['rl']))
        lines.append(assign('FLC(74)', c['pinf']))
        lines.append(assign('F(11)', c['nt']))
        for k in range(int(c['nt'])):
            lines.append(assign(f'F({1 + k})', c['time'][k]))
            lines.append(assign(f'F({12 + k})', c['fc'][k]))
            lines.append(assign(f'F({22 + k})', c['alpha'][k]))
            lines.append(f"      LAMNR({k + 1})="
                         f".{'TRUE' if c['laminar'][k] else 'FALSE'}.")
        for idx, key in ((32, 'me'), (33, 'isp'), (34, 'span'), (35, 'phe'),
                         (36, 'gp'), (37, 'cc'), (38, 'l')):
            lines.append(assign(f'F({idx})', c[key]))
        lines.append('      CALL TRANJT')
        lines.append("      WRITE(6,'(A,150ES25.16)') 'JET',(JET(K),K=1,150)")
        lines.append("      WRITE(6,'(A,26ES25.16)') 'JETA',JETA")
    lines += ['      STOP', '      END', '']
    return '\n'.join(lines)


def main():
    all_cases = cases()
    records = parse_records(run('tranjt', driver(all_cases), ROUTINES))
    print(save('tranjt', [{'inputs': c, 'outputs': r}
                          for c, r in zip(all_cases, records)]))


if __name__ == '__main__':
    main()
