"""
Probe FLAPCM (with GDELTA), the flap pitching moment, and save the
fixture.

The cases run in one program, in order: FLAPCM keeps its span stations
``ET`` (which it edits) and ``KINBD``/``KOUTBD`` between calls.

Run from the repository root: ``python test_parity/probes/flapcm.py``.
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from probe import assign, parse_records, run, save  # noqa: E402

ROUTINES = ['flapcm', 'gdelta', 'agenr', 'simul4', 'det4', 'tbfunx', 'quad',
            'tlinex', 'tlin1x', 'glook', 'switch', 'interx', 'tlin3x',
            'trapz']

PATCHES = {'gdelta': [('COMMON /FLOLOG/ X(20),ASYFP',
                       'COMMON /FLOLOG/ IX(20),ASYFP')]}

COMMONS = """\
      COMMON /OVERLY/ NLOG,NMACH,II,NALPHA
      COMMON /CONSNT/ PI,DR,UNUSED,RAD
      COMMON /FLOLOG/ FLTC,OPTI,BO,WGPL,WGSC,SYNT,HTPL,HTSC,VTPL,VTSC,
     1                HEAD,PRPOWR,JETPOW,LOASRT,TVTPAN,
     2                SUPERS,SUBSON,TRANSN,HYPERS,
     3                SYMFP,ASYFP,TRIMC,TRIM
      COMMON /FLGTCD/ FLC(93)
      COMMON /FLAPIN/ F(69)
      COMMON /SYNTSS/ SYNA(19)
      COMMON /POWR/   PW(104),FLP(189)
      COMMON /SUPDW/  DW(35),TCD(58)
      COMMON /SUPWH/  FCM(282)
      COMMON /WINGD/  A(195)
      COMMON /HTDATA/ AHT(195)
      COMMON /OPTION/ SREF,CBARR
      COMMON /IWING/  PWING, WING(400)
      COMMON /WINGI/  WINGIN(77)
      COMMON /HTI/    HTIN(131)
      LOGICAL FLTC,OPTI,BO,WGPL,WGSC,SYNT,HTPL,HTSC,VTPL,VTSC,
     1        HEAD,PRPOWR,JETPOW,LOASRT,TVTPAN,
     2        SUPERS,SUBSON,TRANSN,HYPERS,
     3        SYMFP,ASYFP,TRIMC,TRIM
"""

UNUSED = 1.0e-30


def case(ftype=1.0, eta1=0.2, eta5=0.7, sdcl=False, scmd=False, htpl=False,
         mach=0.4):
    f = [0.0] * 69
    for i, d in enumerate((10.0, 20.0, 30.0)):
        f[i] = d
        f[18 + i] = 0.3 + 0.1 * i if sdcl else UNUSED
        f[28 + i] = -0.05 * (i + 1) if scmd else UNUSED
        f[38 + i] = 13.0 + i
        f[48 + i] = 7.0 + 0.5 * i
    if not sdcl:
        f[18] = UNUSED
    if not scmd:
        f[28] = UNUSED
    f[11], f[12], f[15], f[16] = 2.5, 1.5, 3.0, ftype
    flp = [0.0] * 189
    flp[0], flp[4] = eta1, eta5
    flp[59] = 0.25 * (eta5 - eta1)
    for k in range(12):
        flp[109 + k] = 0.2 + 0.05 * k
        flp[149 + k] = 0.4 + 0.02 * k
    surface = {'clasec': 0.1, 'tanc4': 0.3, 'swstr': 300.0, 'bsto2': 15.0,
               'bo2': 17.0, 'tante': -0.05, 'tanle': 0.35, 'x': 10.0,
               'cr': 12.0, 'tapexp': 0.35, 'arstar': 6.0, 'alpo': -2.0,
               'cmo': -0.05}
    return {'htpl': htpl, 'asyfp': False, 'mach': mach, 'ii': 1,
            'surface': surface, 'f': f, 'flp': flp, 'xcg': 15.0,
            'sref': 300.0, 'cbarr': 8.0,
            'fcm': [UNUSED] * 282, 'tcd': [0.0] * 58}


def cases():
    return [
        case(),
        case(sdcl=True),
        case(scmd=True),
        case(ftype=2.0),
        case(ftype=5.0),
        case(ftype=6.0),
        case(eta1=0.4153, eta5=0.6549),                # on the stations
        case(htpl=True, eta1=0.1, eta5=0.5, mach=0.6),
        case(eta1=0.05, eta5=0.95),
        case(),                                         # after ET edits
    ]


def driver(all_cases) -> str:
    lines = ['      PROGRAM PROBE', COMMONS.rstrip('\n'),
             '      PI=3.141592654', '      DR=0.01745329',
             '      UNUSED=1.E-30', '      RAD=57.2957795', '      II=1',
             '      ASYFP=.FALSE.']
    for n, c in enumerate(all_cases):
        s = c['surface']
        a, blk, inn = (('AHT', 'HT', 'HTIN') if c['htpl'] else
                       ('A', 'W', 'WINGIN'))
        lines.append(f"      WRITE(6,'(A,I4)') 'CASE',{n}")
        lines.append(f"      HTPL=.{'TRUE' if c['htpl'] else 'FALSE'}.")
        lines.append(assign('FLC(3)', c['mach']))
        for k, v in enumerate(c['f']):
            if v:
                lines.append(assign(f'F({k + 1})', v))
            else:
                lines.append(f'      F({k + 1})=0.')
        for k, v in enumerate(c['flp']):
            if v:
                lines.append(assign(f'FLP({k + 1})', v))
        lines.append(f'      DO {3000 + n} K=1,282')
        lines.append('         FCM(K)=UNUSED')
        lines.append('         IF(K.LE.58) TCD(K)=0.')
        lines.append(f' {3000 + n} CONTINUE')
        for name, key in ((f'{inn}(21)', 'clasec'), (f'{a}(68)', 'tanc4'),
                          (f'{a}(3)', 'swstr'), (f'{inn}(3)', 'bsto2'),
                          (f'{inn}(4)', 'bo2'), (f'{a}(80)', 'tante'),
                          (f'{a}(62)', 'tanle'), (f'{a}(10)', 'cr'),
                          (f'{a}(27)', 'tapexp'), (f'{a}(7)', 'arstar'),
                          (f'{a}(134)', 'alpo'), (f'{inn}(61)', 'cmo')):
            lines.append(assign(name, s[key]))
        lines.append(assign('SYNA(6)' if c['htpl'] else 'SYNA(2)', s['x']))
        lines.append(assign('SYNA(1)', c['xcg']))
        lines.append(assign('SREF', c['sref']))
        lines.append(assign('CBARR', c['cbarr']))
        lines.append('      CALL FLAPCM')
        lines.append("      WRITE(6,'(A,10ES25.16)') 'D',(WING(K),K=211,"
                     "220)")
        lines.append("      WRITE(6,'(A,282ES25.16)') 'FCM',FCM")
        lines.append("      WRITE(6,'(A,58ES25.16)') 'TCD',TCD")
        lines.append("      WRITE(6,'(A,2ES25.16)') 'ETA',FLP(1),FLP(5)")
    lines += ['      STOP', '      END', '']
    return '\n'.join(lines)


def main():
    all_cases = cases()
    records = parse_records(run('flapcm', driver(all_cases), ROUTINES,
                                layout_patches=PATCHES))
    print(save('flapcm', [{'inputs': c, 'outputs': r}
                          for c, r in zip(all_cases, records)]))


if __name__ == '__main__':
    main()
