"""
Probe M55O67 (JETFP), the jet-flap increments, and save the fixture.

The cases run in one program, in order: JETFP's local ``ETAT`` is saved
between calls, so the replay carries it.

Run from the repository root: ``python test_parity/probes/jetfp.py``.
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from probe import assign, parse_records, run, save  # noqa: E402

ROUTINES = ['m55o67', 'jetfp', 'interx', 'tlin1x', 'tlinex', 'tlin3x',
            'glook', 'switch', 'tbfunx', 'quad']

COMMONS = """\
      COMMON /IWING/  PWING,WING(400)
      COMMON /POWER/  P(12),PWIN(26)
      COMMON /OPTION/ SREF,CBARR
      COMMON /CONSNT/ PI,DEG,UNUSED,RAD
      COMMON /SYNTSS/ XCG,XW,ZW,ALIW,ZCG,XH,ZH,ALIH,XV,VERTUP,HINAX,
     1                XVF,SCALE,ZV,ZVF,YV,YF,PHIV,PHIF
      COMMON /FLAPIN/ F(84)
      COMMON /WINGD/  A(195)
      COMMON /WINGI/  WINGIN(100)
      COMMON /OVERLY/ IJKDUM(8),NOVLY
      LOGICAL VERTUP
"""

UNUSED = 1.0e-30


def case(jetflp=1.0, ftype=2.0, cmu=0.5, jevloc=-2.0, jeangl=UNUSED,
         tapr=0.4, jelloc=9.0):
    f = {12: 2.5, 13: 1.8, 14: 4.0, 15: 16.0, 16: 3.0, 17: ftype,
         63: cmu, 74: jetflp}
    for n, (d, cpi, cpo, dj, ej) in enumerate(zip(
            (10.0, 20.0, 40.0), (10.5, 11.0, 11.5), (7.0, 7.3, 7.6),
            (15.0, 30.0, 50.0), (20.0, 35.0, 45.0))):
        f.update({1 + n: d, 39 + n: cpi, 49 + n: cpo, 64 + n: dj,
                  75 + n: ej})
    return {
        'f': f,
        'jet': {'aietlj': 2.0, 'jevloc': jevloc, 'jealoc': 8.0,
                'jeangl': jeangl, 'jelloc': jelloc, 'jerad': 1.0},
        'a': {4: 300.0, 38: 0.5, 118: tapr, 120: 8.0},
        'win': {1: 4.0, 3: 23.0, 4: 25.0, 6: 10.0, 11: -2.0, 16: 0.12},
        'claub': [0.09, 0.095, 0.1],
        'position': {'xcg': 20.0, 'xw': 15.0, 'zw': 0.0, 'zcg': 0.5},
        'sref': 300.0, 'cbarr': 8.0,
        'stale': {'deccl': [7.0] * 3, 'delcm': [8.0] * 3,
                  'dclmax': [9.0] * 3, 'clab': [6.0] * 3},
    }


def cases():
    return [
        case(),
        case(jetflp=2.0),
        case(jetflp=3.0, ftype=4.0),                # EBF, double slotted
        case(jetflp=3.0),
        case(jetflp=4.0),                           # combination
        case(jetflp=5.0),
        case(ftype=1.0),                            # plain: C' = C
        case(jetflp=2.0, ftype=6.0),                # no moment
        case(jetflp=3.0, ftype=4.0, jevloc=10.0),   # jet misses the flap
        case(cmu=UNUSED),
        case(tapr=1.0, jeangl=10.0),
        case(jetflp=3.0, ftype=4.0, jelloc=20.0),
        case(jetflp=3.0, ftype=4.0, jelloc=3.0),
    ]


def driver(all_cases) -> str:
    lines = ['      PROGRAM PROBE', COMMONS.rstrip('\n'),
             '      PI=3.141592654', '      DEG=0.01745329',
             '      UNUSED=1.E-30', '      RAD=57.2957795']
    for n, c in enumerate(all_cases):
        lines.append(f"      WRITE(6,'(A,I4)') 'CASE',{n}")
        for k, v in c['f'].items():
            lines.append(assign(f'F({k})', v))
        for name, index in [('aietlj', 1), ('jevloc', 5), ('jealoc', 6),
                            ('jeangl', 8), ('jelloc', 12), ('jerad', 15)]:
            lines.append(assign(f'PWIN({index})', c['jet'][name]))
        for k, v in c['a'].items():
            lines.append(assign(f'A({k})', v))
        for k, v in c['win'].items():
            lines.append(assign(f'WINGIN({k})', v))
        for k, v in enumerate(c['claub']):
            lines.append(assign(f'WING({101 + k})', v))
        for name, v in c['position'].items():
            lines.append(assign(name.upper(), v))
        lines.append(assign('SREF', c['sref']))
        lines.append(assign('CBARR', c['cbarr']))
        for key, start in [('deccl', 201), ('delcm', 211),
                           ('dclmax', 221), ('clab', 241)]:
            for k, v in enumerate(c['stale'][key]):
                lines.append(assign(f'WING({start + k})', v))
        lines.append('      CALL M55O67')
        lines.append("      WRITE(6,'(A,13ES25.16)') 'R',")
        lines.append('     1(WING(K),K=201,203),(WING(K),K=211,213),')
        lines.append('     2(WING(K),K=221,223),(WING(K),K=241,243),PWIN(8)')
    lines += ['      STOP', '      END', '']
    return '\n'.join(lines)


def main():
    all_cases = cases()
    records = parse_records(run('jetfp', driver(all_cases), ROUTINES))
    payload = []
    for c, r in zip(all_cases, records):
        for key in ('f', 'a', 'win'):
            c[key] = {str(k): v for k, v in c[key].items()}
        payload.append({'inputs': c, 'outputs': r})
    print(save('jetfp', payload))


if __name__ == '__main__':
    main()
