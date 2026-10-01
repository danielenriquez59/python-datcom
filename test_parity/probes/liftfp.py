"""
Probe LIFTFP, the flap lift, lift-curve slope and maximum lift, and save
the fixture.

Run from the repository root: ``python test_parity/probes/liftfp.py``.
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from probe import assign, parse_records, run, save  # noqa: E402

ROUTINES = ['liftfp', 'interx', 'tlin1x', 'tlinex', 'tlin3x', 'glook',
            'switch', 'tbfunx', 'quad']

COMMONS = """\
      COMMON /OVERLY/ NLOG,NMACH,M,NALPHA
      COMMON /FLGTCD/ FLC(93)
      COMMON /FLOLOG/ FLTC,OPTI,BO,WGPL,WGSC,SYNT,HTPL,HTSC,VTPL,VTSC,
     1                HEAD,PRPOWR,JETPOW,LOASRT,TVTPAN,
     2                SUPERS,SUBSON,TRANSN,HYPERS,
     3                SYMFP,ASYFP,TRIMC,TRIM
      COMMON /CONSNT/ PI,DR,UNUSED,RAD
      COMMON /FLAPIN/ F(116)
      COMMON /POWR/   PW(104),FLP(189)
      COMMON /OPTION/ SREF,CBARR,RUFF,BLREF
      COMMON /WINGI/  WINGIN(77)
      COMMON /SUPWH/  FCM(287)
      COMMON /HTI/    HTIN(131)
      COMMON /WINGD/  A(195)
      COMMON /HTDATA/ AHT(195)
      COMMON /IWING/  PWING, WING(400)
      COMMON /IHT/    PHT, HT(380)
      COMMON /SBETA/  STB(135), TRA(108), TRAH(108)
      LOGICAL   FLTC,OPTI,BO,WGPL,WGSC,SYNT,HTPL,HTSC,VTPL,VTSC,
     1          HEAD,PRPOWR,JETPOW,LOASRT,TVTPAN,
     2          SUPERS,SUBSON,TRANSN,HYPERS,
     3          SYMFP,ASYFP,TRIMC,TRIM
"""

UNUSED = 1.0e-30
M = 1
SURFACE = {'cbarex': 7.0, 'bo2': 13.0, 'tante': -0.05, 'tanle': 0.35,
           'cr': 9.0, 'tapri': 0.4, 'bstro2': 4.0, 'tanteo': -0.02,
           'tanleo': 0.4, 'cb': 5.0, 'tapro': 0.6, 'aw': 7.0,
           'btheo': 15.0, 'tovc': 0.12, 'tovco': 0.1, 'cosc4': 0.95,
           'clasec': 0.105, 'cla': 0.08}
WORDS = {'cbarex': ('A', 16), 'bo2': ('IN', 3), 'tante': ('A', 80),
         'tanle': ('A', 62), 'cr': ('IN', 6), 'tapri': ('A', 25),
         'bstro2': ('IN', 2), 'tanteo': ('A', 104), 'tanleo': ('A', 86),
         'cb': ('IN', 5), 'tapro': ('A', 28), 'aw': ('A', 120),
         'btheo': ('IN', 4), 'tovc': ('IN', 16), 'tovco': ('IN', 66),
         'cosc4': ('A', 67), 'clasec': ('IN', 20 + M)}


def case(ftype=1.0, bif=2.0, bof=8.0, sdcl=None, htpl=False, transn=False,
         deltas=(10.0, 20.0, 0.0), tanphe=0.02):
    f = [0.0] * 116
    for i, d in enumerate(deltas):
        f[i] = d
        f[18 + i] = UNUSED if sdcl is None or sdcl[i] is None else sdcl[i]
        f[38 + i] = 9.5 + 0.3 * i
        f[48 + i] = 6.2 + 0.2 * i
        f[84 + i] = 8.0 + 0.1 * i
        f[94 + i] = 5.0 + 0.1 * i
        f[104 + i] = 15.0 + 5.0 * i
    f[10], f[11], f[12], f[13], f[14] = tanphe, 2.2, 1.6, bif, bof
    f[15], f[16] = float(len(deltas)), ftype
    f[114], f[115] = 0.9, 0.7
    return {'ftype': ftype, 'htpl': htpl, 'transn': transn, 'm': M,
            'mach': 0.4, 'rl': 2.5e6, 'sref': 300.0, 'surface': SURFACE,
            'tra70': 0.07, 'win69': 0.09, 'f': f,
            'flp': [0.001 * (k + 1) for k in range(189)],
            'fcm282': [0.3 + 0.01 * k for k in range(6)],
            'wing': [0.5 + 0.001 * k for k in range(50)]}


def cases():
    out = [case(float(t)) for t in range(1, 9)]
    out += [
        case(1.0, sdcl=(0.4, 0.7, 0.2)),               # section lift given
        case(2.0, sdcl=(0.4, 0.7, 0.2)),               # ...translating flap
        case(1.0, sdcl=(None, 0.7, None)),             # switches mid-loop
        case(3.0, htpl=True, bif=0.5, bof=4.0),
        case(1.0, transn=True, bif=10.0, bof=14.0),    # outboard panel
        case(4.0, deltas=(5.0, 40.0, 60.0)),           # double slotted
        case(5.0, deltas=(-10.0, 30.0, 45.0)),
    ]
    return out


def driver(all_cases) -> str:
    lines = ['      PROGRAM PROBE', COMMONS.rstrip('\n'),
             '      PI=3.141592654', '      DR=0.01745329',
             '      UNUSED=1.E-30', '      RAD=57.2957795', f'      M={M}']
    for n, c in enumerate(all_cases):
        a, win, blk, tra = (('AHT', 'HTIN', 'HT', 'TRAH') if c['htpl'] else
                            ('A', 'WINGIN', 'WING', 'TRA'))
        lines.append(f"      WRITE(6,'(A,I4)') 'CASE',{n}")
        lines.append(f"      HTPL=.{'TRUE' if c['htpl'] else 'FALSE'}.")
        lines.append(f"      TRANSN=.{'TRUE' if c['transn'] else 'FALSE'}.")
        lines.append(assign(f'FLC({M + 2})', c['mach']))
        lines.append(assign(f'FLC({M + 42})', c['rl']))
        lines.append(assign('SREF', c['sref']))
        for key, (kind, idx) in WORDS.items():
            name = a if kind == 'A' else win
            lines.append(assign(f'{name}({idx})', c['surface'][key]))
        lines.append(assign(f'{blk}(101)', c['surface']['cla']))
        lines.append(assign(f'{tra}(70)', c['tra70']))
        lines.append(assign(f'{win}(69)', c['win69']))
        for k, v in enumerate(c['f']):
            lines.append(assign(f'F({k + 1})', v) if v else
                         f'      F({k + 1})=0.')
        for k, v in enumerate(c['flp']):
            lines.append(assign(f'FLP({k + 1})', v))
        for k, v in enumerate(c['fcm282']):
            lines.append(assign(f'FCM({282 + k})', v))
        for k, v in enumerate(c['wing']):
            lines.append(assign(f'WING({201 + k})', v))
        lines.append('      CALL LIFTFP')
        lines.append("      WRITE(6,'(A,189ES25.16)') 'FLP',FLP")
        lines.append("      WRITE(6,'(A,116ES25.16)') 'F',F")
        lines.append("      WRITE(6,'(A,6ES25.16)') 'FCM',(FCM(K),K=282,"
                     "287)")
        lines.append("      WRITE(6,'(A,50ES25.16)') 'WING',(WING(K),K=201,"
                     "250)")
    lines += ['      STOP', '      END', '']
    return '\n'.join(lines)


def main():
    all_cases = cases()
    records = parse_records(run('liftfp', driver(all_cases), ROUTINES))
    print(save('liftfp', [{'inputs': c, 'outputs': r}
                          for c, r in zip(all_cases, records)]))


if __name__ == '__main__':
    main()
