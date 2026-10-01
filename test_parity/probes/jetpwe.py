"""
Probe M30O36 (JETPWE, FG6115 and the slope pass) and save the fixture.

The cases run in one program, in order: JETPWE's local ``COSAIH`` is
saved between calls (``-fno-automatic``, as the program is built), so a
case without a tail reads the previous case's value.

Run from the repository root: ``python test_parity/probes/jetpwe.py``.
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from probe import assign, parse_records, run, save  # noqa: E402

ROUTINES = ['m30o36', 'jetpwe', 'fg6115', 'tlinvs', 'tlinex', 'tlin1x',
            'glook', 'switch', 'tbfunx', 'quad']

COMMONS = """\
      COMMON /POWER/  XX(12),AIETLJ,NENGSJ,THSTCJ,JIALOC,JEVLOC,JEALOC,
     1                JINLTA,JEANGL,JEVELO,AMBTMP,JESTMP,JELLOC,JETOTP,
     2                AMBSTP,JERAD
      COMMON /POWR/   PW(315)
      COMMON /FLOLOG/ FLTC,OPTI,BO,WGPL,WGSC,SYNT,HTPL,HTSC,VTPL,VTSC,
     1                HEAD,PRPOWR,JETPOW,LOASRT
      COMMON /CONSNT/ PI,DEG,UNUSED,RAD
      COMMON /SYNTSS/ XCG,XW,ZW,ALIW,ZCG,XH,ZH,ALIH,XV,VERTUP,HINAX,
     1                XVF,SCALE,ZV,ZVF,YV,YF,PHIV,PHIF
      COMMON /OVERLY/ NLOG,NMACH,M,NALPHA,IG,IJKDUM(3),NOVLY
      COMMON /FLGTCD/ FLC(93)
      COMMON /OPTION/ SR,CBARR
      COMMON /WINGD/  A(195),B(49)
      COMMON /HTDATA/ AHT(195)
      COMMON /WINGI/  WINGIN(77)
      COMMON /HTI/    HTIN(131)
      COMMON /IHT/    PHT, HT(380)
      COMMON /IDWASH/ PDWASH,DWASH(60)
      COMMON /IPOWER/ PP, POWER(200)
      LOGICAL VERTUP
      LOGICAL FLTC,OPTI,BO,WGPL,WGSC,SYNT,HTPL,HTSC,VTPL,VTSC,
     1        HEAD,PRPOWR,JETPOW,LOASRT
      REAL NENGSJ,JIALOC,JEVLOC,JEALOC,JINLTA,JEANGL,JEVELO,JESTMP,
     1     JELLOC,JETOTP,JERAD
"""

UNUSED = 1.0e-30
ALPHA = [-4.0, 0.0, 4.0, 8.0, 12.0]
JET = {'aietlj': 2.0, 'nengsj': 2.0, 'thstcj': 0.05, 'jialoc': 10.0,
       'jevloc': -1.0, 'jealoc': 30.0, 'jinlta': 4.0, 'jeangl': 15.0,
       'jevelo': 800.0, 'ambtmp': 520.0, 'jestmp': 1200.0, 'jelloc': 4.0,
       'jetotp': 3000.0, 'ambstp': 2116.0, 'jerad': 1.0}


def case(mach=0.6, htpl=True, outboard=UNUSED, **jet):
    return {'alpha': ALPHA, 'mach': mach, 'htpl': htpl,
            'jet': dict(JET, **jet),
            'win': {1: 3.0, 2: outboard, 4: 20.0, 5: 6.0, 6: 10.0},
            'a': {38: 0.3, 86: 0.5, 120: 6.0, 134: -1.0},
            'position': {'xcg': 20.0, 'xw': 12.0, 'aliw': 1.0, 'zcg': 0.0,
                         'xh': 40.0, 'zh': 2.0, 'alih': -1.0},
            'tail': {'span': 6.0, 'cla': 0.06, 'xbarr': 1.2,
                     'q': [0.95 - 0.005 * a for a in ALPHA]},
            'sref': 300.0, 'cbarr': 8.0, 'comp1': -0.37}


def cases():
    return [
        case(),                                    # subsonic jet, part A
        case(mach=0.2),                            # parts A and B
        case(mach=0.1),                            # parts B and C
        case(mach=0.07),                           # parts B and C
        case(mach=0.05),                           # part C
        case(jerad=0.5),                           # far downstream
        case(jevelo=1500.0),                       # supersonic jet
        case(htpl=False),                          # no tail: stale COSAIH
        case(outboard=8.0),                        # inboard jet, cranked
        case(outboard=8.0, jelloc=15.0),           # outboard jet
        case(nengsj=1.0, jelloc=0.0),              # single engine
        case(mach=0.1242),                         # part B alone
        case(jevelo=1300.0, jetotp=6000.0, mach=0.9),
    ]


def driver(all_cases) -> str:
    lines = ['      PROGRAM PROBE', COMMONS.rstrip('\n'),
             '      PI=3.141592654', '      DEG=0.01745329',
             '      UNUSED=1.E-30', '      RAD=57.2957795', '      M=1']
    for n, c in enumerate(all_cases):
        lines.append(f"      WRITE(6,'(A,I4)') 'CASE',{n}")
        lines.append(f'      DO {3000 + n} K=1,315')
        lines.append('         PW(K)=0.')
        lines.append('         IF(K.LE.200) POWER(K)=0.')
        lines.append(f' {3000 + n} CONTINUE')
        lines.append(assign('PW(63)', c['comp1']))
        lines.append(f"      HTPL=.{'TRUE' if c['htpl'] else 'FALSE'}.")
        na = len(c['alpha'])
        lines.append(f'      NALPHA={na}')
        for j, a in enumerate(c['alpha']):
            lines.append(assign(f'FLC({23 + j})', a))
            lines.append(assign(f'DWASH({1 + j})', c['tail']['q'][j]))
        lines.append(assign('FLC(3)', c['mach']))
        for k, v in c['jet'].items():
            lines.append(assign(k.upper(), v))
        for k, v in c['win'].items():
            lines.append(assign(f'WINGIN({k})', v))
        for k, v in c['a'].items():
            lines.append(assign(f'A({k})', v))
        for k, v in c['position'].items():
            lines.append(assign(k.upper(), v))
        t = c['tail']
        lines.append(assign('HTIN(4)', t['span']))
        lines.append(assign('HT(101)', t['cla']))
        lines.append(assign('AHT(161)', t['xbarr']))
        lines.append(assign('SR', c['sref']))
        lines.append(assign('CBARR', c['cbarr']))
        lines.append('      CALL M30O36')
        lines.append("      WRITE(6,'(A,155ES25.16)') 'PW',")
        lines.append('     1(PW(K),K=1,155)')
        lines.append("      WRITE(6,'(A,140ES25.16)') 'POWER',")
        lines.append('     1(POWER(K),K=1,140)')
        lines.append("      WRITE(6,'(A,ES25.16)') 'WIN1',WINGIN(1)")
    lines += ['      STOP', '      END', '']
    return '\n'.join(lines)


def main():
    all_cases = cases()
    records = parse_records(run('jetpwe', driver(all_cases), ROUTINES))
    payload = []
    for c, r in zip(all_cases, records):
        c['win'] = {str(k): v for k, v in c['win'].items()}
        c['a'] = {str(k): v for k, v in c['a'].items()}
        payload.append({'inputs': c, 'outputs': r})
    print(save('jetpwe', payload))


if __name__ == '__main__':
    main()
