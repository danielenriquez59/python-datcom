"""
Probe SYPBOD, the supersonic body alone, and save the fixture.

The cases run in one program, in order: SYPBOD's local ``RACH`` is saved
between calls.

Run from the repository root: ``python test_parity/probes/sypbod.py``.
"""

import math
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from probe import assign, parse_records, run, save  # noqa: E402

ROUTINES = ['sypbod', 'interx', 'tlin1x', 'tlinex', 'tlin3x', 'glook',
            'switch', 'tbfunx', 'quad', 'trapz', 'fig26']

COMMONS = """\
      COMMON /SUPBOD/ SBD(229)
      COMMON /OVERLY/ NLOG,NMACH,M,NALPHA
      COMMON /OPTION/ SR,CRBAR,RUFF,BLREF
      COMMON /IBODY/  PBODY,BODY(400)
      COMMON /BODYI/  NXX,XCOOR(20),S(20),PERIM(20),RADIUS(20),
     1                ZU(20),ZL(20),BNOSE,BTAIL,RLN,RLA,DS,
     2                ITYPE,METHOD,ELLIP
      COMMON /SYNTSS/ XCG,XW,ZW,ALIW,ZCG,XH,ZH,ALIH,XV,VERTUP,HINAX,
     1                XVF,SCALE,ZV,ZVF,YV,YF,PHIV,PHIF
      COMMON /CONSNT/ PI,DEG,UNUSED,RAD
      COMMON /FLGTCD/ FLC(73)
      COMMON /FLOLOG/ FLTC,OPTI,BO,WGPL,WGSC,SYNT,HTPL,HTSC,VTPL,VTSC,
     1                HEAD,PRPOWR,JETPOW,LOASRT,TVTPAN,
     2                SUPERS,SUBSON,TRANSN,HYPERS
      COMMON /BDATA/  BD(762)
      LOGICAL         FLTC,OPTI,BO,WGPL,WGSC,SYNT,HTPL,HTSC,VTPL,VTSC,
     1                HEAD,PRPOWR,JETPOW,LOASRT,TVTPAN,
     2                SUPERS,SUBSON,TRANSN,HYPERS
      REAL NXX
"""

UNUSED = 1.0e-30
ALPHA = [-4.0, 0.0, 4.0, 8.0, 16.0]
X = [0.0, 2.0, 5.0, 10.0, 15.0, 20.0, 25.0, 30.0, 35.0, 40.0]
BOAT = [0.0, 0.8, 1.3, 1.8, 2.0, 2.0, 2.0, 1.9, 1.7, 1.5]
FLARE = [0.0, 0.8, 1.3, 1.8, 2.0, 2.0, 2.0, 2.2, 2.5, 2.8]
CYL = [0.0, 0.8, 1.3, 1.8, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0]


def case(profile=None, bnose=0.0, btail=2.0, rln=10.0, rla=15.0,
         mach=2.0, ellip=UNUSED, transn=False, rough=1.6e-4):
    profile = profile or BOAT
    body = {'x': X, 'r': profile,
            's': [math.pi * r * r for r in profile],
            'perim': [2 * math.pi * r for r in profile],
            'bnose': bnose, 'btail': btail, 'rln': rln, 'rla': rla,
            'ellip': ellip}
    return {'alpha': ALPHA, 'mach': mach, 'rl': 3.0e6, 'transn': transn,
            'body': body, 'xcg': 20.0, 'sref': 12.5, 'cbarr': 4.0,
            'rough': rough}


def cases():
    return [
        case(),                                        # ogive, boattail
        case(profile=FLARE),                           # flare
        case(bnose=1.0, btail=1.0),                    # cone, conical tail
        case(profile=CYL, btail=0.0, rla=30.0),        # no tail
        case(rla=0.0),                                 # tail, no centre
        case(mach=3.5, ellip=1.5),                     # elliptic
        case(mach=1.6, ellip=0.6, rough=0.0),          # stale RACH
        case(transn=True),
        case(profile=CYL, btail=0.0, rla=0.0, rln=40.0),
        case(mach=5.0, bnose=1.0, profile=FLARE, btail=1.0),
    ]


def driver(all_cases) -> str:
    lines = ['      PROGRAM PROBE', COMMONS.rstrip('\n'),
             '      PI=3.141592654', '      DEG=0.01745329',
             '      UNUSED=1.E-30', '      RAD=57.2957795']
    for n, c in enumerate(all_cases):
        b = c['body']
        lines.append(f"      WRITE(6,'(A,I4)') 'CASE',{n}")
        lines.append(f'      DO {3000 + n} K=1,229')
        lines.append('         SBD(K)=0.')
        lines.append('         BODY(K)=0.')
        lines.append(f' {3000 + n} CONTINUE')
        lines.append(f"      NALPHA={len(c['alpha'])}")
        for j, a in enumerate(c['alpha']):
            lines.append(assign(f'FLC({23 + j})', a))
        lines.append(assign('FLC(3)', c['mach']))
        lines.append(assign('FLC(43)', c['rl']))
        lines.append(f"      TRANSN=.{'TRUE' if c['transn'] else 'FALSE'}.")
        lines.append(f"      NXX={len(b['x'])}.")
        for k, (x, r, s, p) in enumerate(zip(b['x'], b['r'], b['s'],
                                             b['perim'])):
            lines.append(assign(f'XCOOR({k + 1})', x))
            lines.append(assign(f'RADIUS({k + 1})', r))
            lines.append(assign(f'S({k + 1})', s))
            lines.append(assign(f'PERIM({k + 1})', p))
        for name in ('bnose', 'btail', 'rln', 'rla', 'ellip'):
            lines.append(assign(name.upper(), b[name]))
        lines.append(assign('XCG', c['xcg']))
        lines.append(assign('SR', c['sref']))
        lines.append(assign('CRBAR', c['cbarr']))
        lines.append(assign('RUFF', c['rough']))
        lines.append('      CALL SYPBOD(1)')
        lines.append("      WRITE(6,'(A,229ES25.16)') 'SBD',SBD")
        lines.append("      WRITE(6,'(A,100ES25.16)') 'BODY',"
                     "(BODY(K),K=1,100)")
        lines.append("      WRITE(6,'(A,100ES25.16)') 'BD',"
                     "(BD(K),K=536,635)")
        lines.append("      WRITE(6,'(A,ES25.16)') 'E',ELLIP")
    lines += ['      STOP', '      END', '']
    return '\n'.join(lines)


def main():
    all_cases = cases()
    records = parse_records(run('sypbod', driver(all_cases), ROUTINES))
    print(save('sypbod', [{'inputs': c, 'outputs': r}
                          for c, r in zip(all_cases, records)]))


if __name__ == '__main__':
    main()
