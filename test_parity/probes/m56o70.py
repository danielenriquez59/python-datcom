"""
Probe M56O70 (VTAREA for both panels and BDAREA), the Mach-shadow
executive, and save the fixture.

Run from the repository root: ``python test_parity/probes/m56o70.py``.
"""

import math
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from probe import assign, parse_records, run, save  # noqa: E402

ROUTINES = ['m56o70', 'vtarea', 'ptint1', 'area1', 'bdarea', 'ptint2',
            'area2']

COMMONS = """\
      COMMON /OVERLY/ NLOG,NMACH,I,NALPHA,IG,IJKDUM(3),NOVLY
      COMMON /VTDATA/ AVT(195), AVF(195)
      COMMON /SYNTSS/ SYNA(19)
      COMMON /HTI/    HTIN(154)
      COMMON /VTI/    VTIN(154),TVTIN(8),VFIN(154)
      COMMON /WINGD/  A(195)
      COMMON /HTDATA/ AHT(195)
      COMMON /WINGI/  WINGIN(100)
      COMMON /FLGTCD/ FLC(73)
      COMMON /BODYI/  XNX,XB(20),SB(20),PB(20),RB(20)
      COMMON /FLOLOG/ FLTC,OPTI,BO,WGPL,WGSC,SYNT,HTPL,HTSC,VTPL,VTSC,
     1                HEAD,PRPOWR,JETPOW,LOASRT,TVTPAN,SUPERS,SUBSON,
     2                TRANSN,HYPERS,SYMFP,ASYFP,TRIMC,TRIM,DAMP,
     3                HYPEF,TRAJET,BUILD,FIRST,DRCONV,PART,
     4                VFPL,VFSC,CTAB
      COMMON /CONSNT/ PI,DEG,UNUSED,RAD
      LOGICAL  FLTC,OPTI,BO,WGPL,WGSC,SYNT,HTPL,HTSC,VTPL,VTSC,
     1         HEAD,PRPOWR,JETPOW,LOASRT,TVTPAN,SUPERS,SUBSON,
     2         TRANSN,HYPERS,SYMFP,ASYFP,TRIMC,TRIM,DAMP,
     3         HYPEF,TRAJET,BUILD,FIRST,DRCONV,PART,
     4         VFPL,VFSC,CTAB
      LOGICAL VERTUP
      EQUIVALENCE (VERTUP,SYNA(10))
      DIMENSION WT(4)
      DATA WT /4HSTRA,4HDOUB,4HCRAN,4HCURV/
"""

UNUSED = 1.0e-30
I = 2
XB = [0.0, 2.0, 5.0, 10.0, 15.0, 20.0, 25.0, 30.0, 35.0, 40.0]
RB = [0.0, 0.8, 1.3, 1.8, 2.0, 2.0, 2.0, 1.9, 1.7, 1.5]


def case(vtpl=True, vfpl=True, htpl=True, alih=-3.0, vt_set=(), ht_set=False,
         xh=30.0, zh=0.5, mach=2.0):
    vtin = {1: 2.0, 3: 5.2, 4: 6.0, 5: 2.0, 15: 1.0,
            94 + I: UNUSED, 114 + I: UNUSED, 134 + I: UNUSED}
    for k in vt_set:
        vtin[k] = 3.3
    vfin = {1: 1.0, 3: 2.2, 4: 2.8, 5: 1.0, 15: 1.0,
            94 + I: UNUSED, 114 + I: UNUSED, 134 + I: UNUSED}
    htin = {3: 4.0, 4: 5.0, 94 + I: UNUSED, 114 + I: UNUSED,
            134 + I: UNUSED}
    if ht_set:
        htin.update({94 + I: 1.1, 114 + I: 2.2, 134 + I: 3.3})
    return {
        'i': I, 'mach': mach, 'vtpl': vtpl, 'vfpl': vfpl, 'htpl': htpl,
        'vertup': True,
        'syna': {1: 20.0, 2: 10.0, 3: 0.0, 4: 2.0, 6: xh, 7: zh, 8: alih,
                 9: 28.0, 12: 30.0, 14: 1.0, 15: -0.6},
        'vtin': vtin, 'vfin': vfin,
        'avt': {1: 18.2, 2: 0.0, 3: 18.2, 10: 5.0, 21: 6.0, 62: 0.8,
                86: 0.0},
        'avf': {1: 4.0, 2: 0.0, 3: 4.0, 10: 2.0, 21: 2.8, 62: 0.5,
                86: 0.0},
        'vt_common': {59: math.atan(0.8), 77: 0.2, 83: math.atan(1.0),
                      101: 0.3},
        'wing': {'span': 15.0, 'spans': 13.0, 'a62': 0.6, 'a10': 6.0},
        'tail': {'span': 5.0, 'spans': 4.0, 'a10': 3.0, 'a16': 2.5,
                 'a30': 0.5, 'a62': 0.7},
        'htin': htin, 'aht': {10: 3.0, 16: 2.5, 30: 0.5, 62: 0.7},
        'body_x': XB, 'body_r': RB,
    }


def cases():
    return [
        case(),
        case(vtpl=False, vt_set=(94 + I, 134 + I)),   # runs: 114+I unset
        case(vtpl=False, vt_set=(94 + I, 114 + I, 134 + I)),
        case(ht_set=True),
        case(xh=37.7, zh=-2.35, alih=2.0, mach=4.0),  # BDAREA aborts
        case(htpl=False, alih=0.0),
    ]


def driver(all_cases) -> str:
    lines = ['      PROGRAM PROBE', COMMONS.rstrip('\n'),
             '      PI=3.141592654', '      DEG=0.01745329',
             '      UNUSED=1.E-30', '      RAD=57.2957795', f'      I={I}']
    for n, c in enumerate(all_cases):
        lines.append(f"      WRITE(6,'(A,I4)') 'CASE',{n}")
        lines.append(f'      DO {3000 + n} K=1,154')
        lines.append('         VTIN(K)=0.')
        lines.append('         VFIN(K)=0.')
        lines.append('         HTIN(K)=0.')
        lines.append(f' {3000 + n} CONTINUE')
        lines.append(assign(f'FLC({I + 2})', c['mach']))
        for flag in ('vtpl', 'vfpl', 'htpl'):
            lines.append(f"      {flag.upper()}="
                         f".{'TRUE' if c[flag] else 'FALSE'}.")
        for k, v in c['syna'].items():
            lines.append(assign(f'SYNA({k})', v))
        lines.append('      VERTUP=.TRUE.')
        for name, key in (('VTIN', 'vtin'), ('VFIN', 'vfin'),
                          ('HTIN', 'htin'), ('AVF', 'avf'), ('AHT', 'aht')):
            for k, v in c[key].items():
                if k == 15 and name in ('VTIN', 'VFIN'):
                    lines.append(f'      {name}(15)=WT({int(v)})')
                else:
                    lines.append(assign(f'{name}({k})', v))
        for k, v in c['avt'].items():
            lines.append(assign(f'AVT({k})', v))
        for k, v in c['vt_common'].items():
            lines.append(assign(f'AVT({k})', v))
        w = c['wing']
        lines.append(assign('WINGIN(4)', w['span']))
        lines.append(assign('WINGIN(3)', w['spans']))
        lines.append(assign('A(62)', w['a62']))
        lines.append(assign('A(10)', w['a10']))
        lines.append(f"      XNX={len(c['body_x'])}.")
        for k, (x, r) in enumerate(zip(c['body_x'], c['body_r'])):
            lines.append(assign(f'XB({k + 1})', x))
            lines.append(assign(f'RB({k + 1})', r))
        lines.append('      IG=-1')
        lines.append('      CALL M56O70')
        m = I
        lines.append("      WRITE(6,'(A,12ES25.16)') 'R',")
        lines.append(f'     1VTIN({94 + m}),VTIN({114 + m}),VTIN({134 + m}),')
        lines.append(f'     2VFIN({94 + m}),VFIN({114 + m}),VFIN({134 + m}),')
        lines.append(f'     3HTIN({94 + m}),HTIN({114 + m}),HTIN({134 + m}),')
        lines.append('     4SYNA(6),SYNA(7)')
        lines.append("      WRITE(6,'(A,I4)') 'IG',IG")
    lines += ['      STOP', '      END', '']
    return '\n'.join(lines)


def main():
    all_cases = cases()
    records = parse_records(run('m56o70', driver(all_cases), ROUTINES))
    payload = []
    for c, r in zip(all_cases, records):
        for key in ('syna', 'vtin', 'vfin', 'htin', 'avt', 'avf', 'aht',
                    'vt_common'):
            c[key] = {str(k): v for k, v in c[key].items()}
        r.pop('_text', None)
        payload.append({'inputs': c, 'outputs': r})
    print(save('m56o70', payload))


if __name__ == '__main__':
    main()
