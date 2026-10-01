"""
Probe M26O32 (HYPBOD and the slope pass) and save the fixture.

Run from the repository root: ``python test_parity/probes/hypbod.py``.
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from probe import assign, parse_records, run, save  # noqa: E402

ROUTINES = ['m26o32', 'hypbod', 'interx', 'tlin1x', 'tlinex', 'tlin3x',
            'glook', 'switch', 'tbfunx', 'quad', 'trapz']

COMMONS = """\
      COMMON /OVERLY/ NLOG,NMACH,I,NALPHA,IG,IJKDUM(3),NOVLY
      COMMON /CONSNT/ PI,DEG,UNUSED,RAD
      COMMON /FLGTCD/ FLC(93)
      COMMON /OPTION/ SR,CBAR,ROUGHC,BLREF
      COMMON /BODYI/  NXX,X(20),S(20),P(20),R(20),ZU(20),ZL(20),BNOSE,
     1                BTAIL,RLN,RLA,DS
      COMMON /SYNTSS/ XCG,XW,ZW,ALIW,ZCG,XH,ZH,ALIH,XV,VERTUP,HINAX,
     1                XVF,SCALE,ZV,ZVF,YV,YF,PHIV,PHIF
      COMMON /SUPBOD/ SBD(229)
      COMMON /IBODY/  PB, BODY(400)
      REAL NXX
"""

X = [0.0, 1.0, 2.5, 4.0, 6.0, 8.0, 10.0, 12.0, 14.0, 16.0]
ALPHA = [-10.0, -2.0, 0.0, 2.0, 5.0, 10.0, 20.0, 35.0]


def case(r, rln=6.0, rla=6.0, ds=0.0, bnose=0.0, mach=8.0, x=None,
         alpha=None):
    return {'alpha': alpha or ALPHA, 'mach': mach, 'xcg': 7.5,
            'sref': 3.2, 'cbar': 2.0, 'blref': 2.4,
            'body': {'x': x or X, 'r': r, 'rln': rln, 'rla': rla, 'ds': ds,
                     'bnose': bnose},
            'stale': {'thetaa': 0.123, 'thetat': -0.0456}}


OGIVE = [0.0, 0.3, 0.55, 0.72, 0.8, 0.82, 0.85, 0.9, 1.0, 1.1]
CONE = [0.05, 0.18, 0.37, 0.56, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8]
BOAT = [0.0, 0.3, 0.55, 0.72, 0.8, 0.8, 0.8, 0.8, 0.7, 0.6]


def cases():
    return [
        case(OGIVE),                                  # both flares
        case(CONE, bnose=1.0, ds=0.1, rla=6.0),       # cone, cylinder
        case(BOAT, rla=0.0),                          # no centre, boattail
        case(OGIVE, rla=10.0),                        # no tail
        case(OGIVE, ds=0.2, mach=5.0),                # blunted ogive
        case(CONE, bnose=1.0, ds=0.1, rla=10.0, mach=12.0),
        case(OGIVE, rln=4.0, rla=5.0, alpha=[0.0, 4.0, 8.0]),
    ]


def driver(all_cases) -> str:
    lines = ['      PROGRAM PROBE', COMMONS.rstrip('\n'),
             '      PI=3.141592654', '      DEG=0.01745329',
             '      UNUSED=1.E-30', '      RAD=57.2957795', '      I=1']
    for n, c in enumerate(all_cases):
        lines.append(f"      WRITE(6,'(A,I4)') 'CASE',{n}")
        lines.append(f'      DO {3000 + n} K=1,400')
        lines.append('         BODY(K)=0.')
        lines.append(f'         IF(K.LE.229) SBD(K)=0.')
        lines.append(f' {3000 + n} CONTINUE')
        na = len(c['alpha'])
        lines.append(f'      NALPHA={na}')
        for j, a in enumerate(c['alpha']):
            lines.append(assign(f'FLC({23 + j})', a))
        lines.append(assign('FLC(3)', c['mach']))
        for name in ('xcg', 'sref', 'cbar', 'blref'):
            target = {'xcg': 'XCG', 'sref': 'SR', 'cbar': 'CBAR',
                      'blref': 'BLREF'}[name]
            lines.append(assign(target, c[name]))
        b = c['body']
        lines.append(f"      NXX={len(b['x'])}.")
        for k, (x, r) in enumerate(zip(b['x'], b['r'])):
            lines.append(assign(f'X({k + 1})', x))
            lines.append(assign(f'R({k + 1})', r))
        for name, key in [('RLN', 'rln'), ('RLA', 'rla'), ('DS', 'ds'),
                          ('BNOSE', 'bnose')]:
            lines.append(assign(name, b[key]))
        lines.append(assign('SBD(130)', c['stale']['thetaa']))
        lines.append(assign('SBD(135)', c['stale']['thetat']))
        lines.append('      CALL M26O32')
        lines.append("      WRITE(6,'(A,229ES25.16)') 'SBD',")
        lines.append('     1(SBD(J),J=1,229)')
        lines.append("      WRITE(6,'(A,200ES25.16)') 'BODY',")
        lines.append('     1(BODY(J),J=1,200)')
    lines += ['      STOP', '      END', '']
    return '\n'.join(lines)


def main():
    all_cases = cases()
    records = parse_records(run('hypbod', driver(all_cases), ROUTINES))
    print(save('hypbod', [{'inputs': c, 'outputs': r}
                          for c, r in zip(all_cases, records)]))


if __name__ == '__main__':
    main()
