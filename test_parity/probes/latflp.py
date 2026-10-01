"""
Probe LATFLP, the rolling and yawing moments of asymmetric flaps, spoilers
and a differentially deflected tail, and save the fixture.

Run from the repository root: ``python test_parity/probes/latflp.py``.
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from probe import assign, parse_records, run, save  # noqa: E402

ROUTINES = ['latflp', 'interx', 'tlin1x', 'tlinex', 'tlin3x', 'tlin4x',
            'glook', 'switch', 'tbfunx', 'quad']

COMMONS = """\
      COMMON /OVERLY/ NLOG,NMACH,M,NALPHA
      COMMON /CONSNT/ PI,DR,UNUSED,RAD
      COMMON /OPTION/ SREF,CBARR,ROUGFC,BLREF
      COMMON /FLGTCD/ FLC(93)
      COMMON /IHT/    PHT,HT(380)
      COMMON /HTI/    HTIN(131)
      COMMON /HTDATA/ AHT(195),BHT(49)
      COMMON /IBODY/  PBODY,BODY(400)
      COMMON /WINGI/  WINGIN(77)
      COMMON /POWR/   PW(59),FLA(45)
      COMMON /FLAPIN/ F(69)
      COMMON /WINGD/  A(195),B(49)
      COMMON /WBHCAL/ WBT(155)
      COMMON /IDWASH/ PDWASH,DWASH(60)
      COMMON /IWING/  PWING,WING(400)
      COMMON /SUPDW/  DW(35),TCD(58)
      COMMON /FLOLOG/ FLTC,OPTI,BO,WGPL,WGSC,SYNT,HTPL,HTSC,VTPL,VTSC,
     1                HEAD,PRPOWR,JETPOW,LOASRT,TVTPAN,
     2                SUPERS,SUBSON,TRANSN,HYPERS,
     3                SYMFP,ASYFP,TRIMC,TRIM
      LOGICAL         FLTC,OPTI,BO,WGPL,WGSC,SYNT,HTPL,HTSC,VTPL,VTSC,
     1                HEAD,PRPOWR,JETPOW,LOASRT,TVTPAN,
     2                SUPERS,SUBSON,TRANSN,HYPERS,
     3                SYMFP,ASYFP,TRIMC,TRIM
"""

UNUSED = 1.0e-30
M = 1
SURFACE = {'bo2': 13.0, 'cr': 9.0, 'tovc': 0.1, 'clasec': 0.105,
           'area': 250.0, 'cbarex': 7.0, 'sweple': 30.0, 'tanle': 0.577,
           'cosc4': 0.9, 'tanc4': 0.45, 'swepte': 10.0, 'coste': 0.985,
           'tante': 0.176, 'taprw': 0.4, 'aw': 6.0, 'cla': 0.08}
# Where each surface word lives: (array, index).
WORDS = {'bo2': ('WINGIN', 4), 'cr': ('WINGIN', 6), 'tovc': ('WINGIN', 16),
         'clasec': ('WINGIN', 20 + M), 'area': ('A', 4),
         'cbarex': ('A', 16), 'sweple': ('A', 58), 'tanle': ('A', 62),
         'cosc4': ('A', 67), 'tanc4': ('A', 68), 'swepte': ('A', 76),
         'coste': ('A', 79), 'tante': ('A', 80), 'taprw': ('A', 118),
         'aw': ('A', 120), 'cla': ('WING', 101)}
HTAIL = {'bo2h': 6.0, 'bo2hst': 5.2, 'sh': 60.0, 'clah': 0.07}


def case(stype, bif=2.0, bof=8.0, ndelta=3, nalpha=4, transn=False,
         surface=None, unused_clw=()):
    s = dict(SURFACE, **(surface or {}))
    f = [0.0] * 69
    for j in range(ndelta):
        f[j] = 0.05 + 0.01 * j                    # DDOC
        f[18 + j] = 5.0 + 7.0 * j                 # DELTAL
        f[28 + j] = -2.0 - 3.0 * j                # DELTAR
        f[38 + j] = 0.06 + 0.02 * j               # DSOC
        f[48 + j] = 0.55 + 0.05 * j               # XSOC
        f[59 + j] = 0.07 + 0.02 * j               # HSOC
    f[10], f[11], f[12], f[13] = 0.03, 2.2, 1.6, bif
    f[14], f[15], f[17], f[58] = bof, float(ndelta), stype, 0.72
    s['clw'] = [UNUSED if n in unused_clw else 0.1 + 0.07 * n
                for n in range(nalpha)]
    return {'stype': stype, 'transn': transn, 'nalpha': nalpha,
            'mach': 0.4, 'rl': 2.5e6, 'sref': 300.0, 'blref': 26.0,
            'alpha': [-2.0 + 2.0 * n for n in range(nalpha)],
            'surface': s, 'htail': HTAIL, 'win69': 0.09,
            'gd': [0.3, 0.25, 0.2, 0.1],
            'dedalp': [0.3 + 0.01 * n for n in range(nalpha)],
            'rivbh': [0.9 - 0.01 * n for n in range(nalpha)],
            'gamvr': [0.2 + 0.02 * n for n in range(nalpha)],
            'f': f,
            'fla': [2.0 * k for k in range(1, 46)],
            'ht201': [0.125 * k for k in range(1, 31)],
            'clrol': [0.5 * k for k in range(1, 201)],
            'cn': [0.25 * k for k in range(1, 201)]}


def cases():
    unswept = {'sweple': 12.0, 'swepte': 10.0}
    return [
        case(4.0, bif=2.0, bof=12.9),                   # flap to the tip
        case(4.0, bif=3.0, bof=9.0),                    # tip superposition
        case(4.0, bif=1.0, bof=6.0, transn=True, unused_clw=(1, 3)),
        case(4.0, bif=2.0, bof=12.9, unused_clw=(0,)),
        case(4.0, bif=2.5, bof=10.0, ndelta=10, nalpha=20),
        case(1.0, surface=unswept),                     # spoiler, Fig 10
        case(1.0, bif=6.0, bof=10.0),                   # Fig 11A0
        case(2.0, bif=3.0, bof=9.0),                    # Fig 11A1
        case(3.0, bif=1.0, bof=10.0),                   # Fig 11A2, slot
        case(3.0, bif=2.0, bof=7.0, surface=unswept),
        case(2.0, bif=0.5, bof=12.5),                   # Fig 11A3
        case(1.0, bif=0.5, bof=12.5, transn=True),
        case(5.0, nalpha=6),                            # differential tail
        case(5.0, ndelta=10, nalpha=20),
    ]


def driver(all_cases) -> str:
    lines = ['      PROGRAM PROBE', COMMONS.rstrip('\n'),
             '      PI=3.141592654', '      DR=0.01745329',
             '      UNUSED=1.E-30', '      RAD=57.2957795', f'      M={M}']
    for n, c in enumerate(all_cases):
        lines.append(f"      WRITE(6,'(A,I4)') 'CASE',{n}")
        lines.append(f"      TRANSN=.{'TRUE' if c['transn'] else 'FALSE'}.")
        lines.append(f"      NALPHA={c['nalpha']}")
        lines.append(assign(f'FLC({M + 2})', c['mach']))
        lines.append(assign(f'FLC({M + 42})', c['rl']))
        for k, v in enumerate(c['alpha']):
            lines.append(assign(f'FLC({23 + k})', v))
        lines.append(assign('SREF', c['sref']))
        lines.append(assign('BLREF', c['blref']))
        for key, (name, idx) in WORDS.items():
            lines.append(assign(f'{name}({idx})', c['surface'][key]))
        lines.append(assign('WINGIN(69)', c['win69']))
        for k, v in enumerate(c['surface']['clw']):
            lines.append(assign(f'WING({21 + k})', v))
        lines.append(assign('HTIN(4)', c['htail']['bo2h']))
        lines.append(assign('HTIN(3)', c['htail']['bo2hst']))
        lines.append(assign('AHT(3)', c['htail']['sh']))
        lines.append(assign('HT(101)', c['htail']['clah']))
        for k, v in enumerate(c['gd']):
            lines.append(assign(f'TCD({43 + k})', v))
        for k in range(c['nalpha']):
            lines.append(assign(f'DWASH({41 + k})', c['dedalp'][k]))
            lines.append(assign(f'WBT({68 + k})', c['rivbh'][k]))
            lines.append(assign(f'WBT({46 + k})', c['gamvr'][k]))
        for k, v in enumerate(c['f']):
            lines.append(assign(f'F({k + 1})', v) if v else
                         f'      F({k + 1})=0.')
        lines += [f'      DO {10 + n} K=1,200',
                  '      IF(K.LE.45) FLA(K)=K*2.0',
                  '      IF(K.LE.30) HT(200+K)=K*0.125',
                  '      WING(200+K)=K*0.5',
                  '      BODY(200+K)=K*0.25',
                  f'{10 + n:5d} CONTINUE']
        lines.append('      CALL LATFLP')
        lines.append("      WRITE(6,'(A,69ES25.16)') 'F',F")
        lines.append("      WRITE(6,'(A,45ES25.16)') 'FLA',FLA")
        lines.append("      WRITE(6,'(A,30ES25.16)') 'HT',(HT(K),K=201,"
                     "230)")
        lines.append("      WRITE(6,'(A,200ES25.16)') 'WING',(WING(K),"
                     "K=201,400)")
        lines.append("      WRITE(6,'(A,200ES25.16)') 'BODY',(BODY(K),"
                     "K=201,400)")
    lines += ['      STOP', '      END', '']
    return '\n'.join(lines)


def main():
    all_cases = cases()
    records = parse_records(run('latflp', driver(all_cases), ROUTINES))
    print(save('latflp', [{'inputs': c, 'outputs': r}
                          for c, r in zip(all_cases, records)]))


if __name__ == '__main__':
    main()
